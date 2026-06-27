// DwbaInputExpander.cpp — expand a DWBA reaction description into a Ptolemy deck.
//
// Logic mirrors Cleopatra's InFileCreator.h (Ryan Tang, 2018) line for line.
// Differences from the original:
//   * masses / Z,A come from PtolemyCpp's mass table (excess(), azCode())
//     instead of Cleopatra's Isotope class;
//   * optical-model parameters come from OpticalPotentialLibrary;
//   * output is built into a std::string instead of written to a FILE*;
//   * diagnostics go to stderr instead of stdout.
// The emitted Ptolemy deck format (field widths, keyword spelling, block order)
// is reproduced exactly so existing decks remain comparable.

#include "DwbaInputExpander.h"
#include "OpticalPotentialLibrary.h"
#include "masstable.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <sstream>
#include <vector>
#include <string>

using std::string;
using std::vector;

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------
static string toUpperStr(const string& s) {
    string r = s;
    for (char& c : r) c = (char)std::toupper((unsigned char)c);
    return r;
}

// Split on runs of whitespace, dropping empty tokens (mirrors SplitStr(line," ")).
static vector<string> splitWS(const string& s) {
    vector<string> out;
    size_t i = 0, n = s.size();
    while (i < n) {
        while (i < n && std::isspace((unsigned char)s[i])) ++i;
        size_t j = i;
        while (j < n && !std::isspace((unsigned char)s[j])) ++j;
        if (j > i) out.push_back(s.substr(i, j - i));
        i = j;
    }
    return out;
}

// spdf orbital letter -> l value.
static int getLValue(const string& spdf) {
    if (spdf == "s") return 0;
    if (spdf == "p") return 1;
    if (spdf == "d") return 2;
    if (spdf == "f") return 3;
    if (spdf == "g") return 4;
    if (spdf == "h") return 5;
    if (spdf == "i") return 6;
    if (spdf == "j") return 7;
    return -1;
}

// "9/2" -> 4.5, "2" -> 2.0
static float getSpinValue(const string& spin) {
    size_t slashPos = spin.find('/');
    if (slashPos != string::npos) {
        int num = std::atoi(spin.substr(0, slashPos).c_str());
        int den = std::atoi(spin.substr(slashPos + 1).c_str());
        if (den != 0) return (float)num / (float)den;
    } else {
        return (float)std::atof(spin.c_str());
    }
    return -1.0f;
}

// Decode an isotope/particle name ("206Hg", "d", "p", "3He", ...) to Z, A using
// the mass table's azCode. Returns false on syntax error / unknown element.
static bool getIso(const string& nameIn, int& Z, int& A) {
    char buf[9];
    for (int i = 0; i < 8; ++i) buf[i] = ' ';
    buf[8] = '\0';
    string up = toUpperStr(nameIn);
    for (size_t i = 0; i < up.size() && i < 8; ++i) buf[i] = up[i];
    int iz = 0, ia = 0, rc = 0;
    azCode(buf, iz, ia, rc);
    if (rc != 0) return false;
    Z = iz; A = ia;
    return true;
}

// Atomic mass excess [MeV]; ok=false if the nuclide is not tabulated.
static bool massEx(int Z, int A, double& out) {
    int nt = 0;
    out = excess(Z, A, nt);
    return nt == 0;
}

// ---------------------------------------------------------------------------
// OM-block emitters (exact field widths from InFileCreator.h)
// ---------------------------------------------------------------------------
static void emitFullOM(std::ostringstream& os, const OMPset& p) {
    char b[128];
    std::snprintf(b, sizeof b, "v    = %7.3f    r0 = %7.3f    a = %7.3f\n", p.v, p.r0, p.a);            os << b;
    std::snprintf(b, sizeof b, "vi   = %7.3f   ri0 = %7.3f   ai = %7.3f\n", p.vi, p.ri0, p.ai);         os << b;
    std::snprintf(b, sizeof b, "vsi  = %7.3f  rsi0 = %7.3f  asi = %7.3f\n", p.vsi, p.rsi0, p.asi);      os << b;
    std::snprintf(b, sizeof b, "vso  = %7.3f  rso0 = %7.3f  aso = %7.3f\n", p.vso, p.rso0, p.aso);      os << b;
    std::snprintf(b, sizeof b, "vsoi = %7.3f rsoi0 = %7.3f asoi = %7.3f  rc0 = %7.3f\n",
                  p.vsoi, p.rsoi0, p.asoi, p.rc0);                                                       os << b;
}

// Reduced 3-line block used by the inelastic branch (rc0 trails the vsi line).
static void emitReducedOM(std::ostringstream& os, const OMPset& p) {
    char b[128];
    std::snprintf(b, sizeof b, "v    = %7.3f    r0 = %7.3f    a = %7.3f\n", p.v, p.r0, p.a);            os << b;
    std::snprintf(b, sizeof b, "vi   = %7.3f   ri0 = %7.3f   ai = %7.3f\n", p.vi, p.ri0, p.ai);         os << b;
    std::snprintf(b, sizeof b, "vsi  = %7.3f  rsi0 = %7.3f  asi = %7.3f  rc0 = %7.3f\n",
                  p.vsi, p.rsi0, p.asi, p.rc0);                                                          os << b;
}

// ---------------------------------------------------------------------------
// looksLikeDwba — auto-detection heuristic
// ---------------------------------------------------------------------------
bool DwbaExpander::looksLikeDwba(const string& input) {
    std::istringstream in(input);
    string line;
    while (std::getline(in, line)) {
        // trim leading whitespace
        size_t s = line.find_first_not_of(" \t\r\n");
        if (s == string::npos) continue;          // blank
        string t = line.substr(s);
        if (t[0] == '#' || t[0] == '!' || t[0] == '$') continue;  // comment
        // First content line decides. A DWBA line starts with
        // target(in,out)residual and has >= 7 whitespace tokens.
        vector<string> tok = splitWS(t);
        if (tok.empty()) return false;
        const string& r = tok[0];
        size_t lp = r.find('('), cm = r.find(','), rp = r.find(')');
        bool reactionShape = lp != string::npos && cm != string::npos &&
                             rp != string::npos && lp < cm && cm < rp &&
                             lp > 0 && rp + 1 < r.size();
        return reactionShape && tok.size() >= 7;
    }
    return false;
}

// ---------------------------------------------------------------------------
// expand
// ---------------------------------------------------------------------------
string DwbaExpander::expand(const string& dwbaInput,
                            double angMin, double angMax, double angStep) {
    std::ostringstream os;
    std::istringstream in(dwbaInput);

    double exN, exH;   // neutron / 1H mass excess for separation energies
    { int nt = 0; exN = excess(0, 1, nt); exH = excess(1, 1, nt); }

    int numOfReaction = 0;
    string tempLine;
    while (std::getline(in, tempLine)) {
        // strip a trailing '\r' (CRLF inputs)
        if (!tempLine.empty() && tempLine.back() == '\r') tempLine.pop_back();

        if (tempLine.substr(0, 1) == "#") continue;
        if (tempLine.size() < 5) continue;

        vector<string> str0 = splitWS(tempLine);
        if (str0.empty()) continue;

        std::fprintf(stderr, "  %s\n", tempLine.c_str());

        // ---- decode the reaction string: target(in,out)residual ----
        const string& rstr = str0[0];
        size_t lp = rstr.find('('), cm = rstr.find(',', lp == string::npos ? 0 : lp),
               rp = rstr.find(')', cm == string::npos ? 0 : cm);
        if (lp == string::npos || cm == string::npos || rp == string::npos ||
            !(lp < cm && cm < rp)) {
            std::fprintf(stderr, "\033[31m  ===> Ignored. Cannot parse reaction '%s'. \033[0m\n", rstr.c_str());
            continue;
        }
        string isoA      = rstr.substr(0, lp);                 // target
        string mass_a    = rstr.substr(lp + 1, cm - lp - 1);   // projectile in
        string mass_b    = rstr.substr(cm + 1, rp - cm - 1);   // projectile out
        string isoB      = rstr.substr(rp + 1);                // residual
        string reactionType = "(" + mass_a + "," + mass_b + ")";

        int Za, Aa, Zb, Ab, ZA, AA, ZB, AB;
        if (!getIso(mass_a, Za, Aa) || !getIso(mass_b, Zb, Ab)) {
            std::fprintf(stderr, "\033[31m  ===> Ignored. Unknown light particle. \033[0m\n");
            continue;
        }

        // ---- reaction-support checks (mirror InFileCreator) ----
        bool isReactionSupported = false;
        bool isTransferReaction  = true;

        if (Aa <= 4 && Za <= 2 && Ab <= 4 && Zb <= 2) isReactionSupported = true;
        // elastic-ish: same projectile in and out
        if (Aa == Ab && Za == Zb) isTransferReaction = false;
        // p/n-exchange not supported
        if (Aa == Ab && Za != Zb) isReactionSupported = false;
        // >= 3 nucleon transfer not supported
        int numNucleonsTransfer = Aa - Ab;
        if (std::abs(numNucleonsTransfer) >= 3) isReactionSupported = false;

        if (!isReactionSupported) {
            std::fprintf(stderr, "\033[31m  ===> Ignored. Reaction type not supported. \033[0m\n");
            continue;
        }

        if (str0.size() < 7) {
            std::fprintf(stderr, "\033[31m  ===> Ignored. Need >= 7 fields. \033[0m\n");
            continue;
        }

        // ---- remaining tokens ----
        string gsSpinparityA = str0[1];
        string gsSpinA  = gsSpinparityA.substr(0, gsSpinparityA.length() - 1);
        string gsParityA = gsSpinparityA.substr(gsSpinparityA.length() - 1);
        double gsSpinAValue = getSpinValue(gsSpinA);
        int gsParityAValue = (gsParityA == "+" ? 1 : -1);

        string orbital = str0[2];

        string spinParity = str0[3];
        string spin   = spinParity.substr(0, spinParity.length() - 1);
        string parity = spinParity.substr(spinParity.length() - 1);
        double spinValue = getSpinValue(spin);
        int parityValue = (parity == "+" ? 1 : -1);

        string Ex             = str0[4];
        string reactionEnergy = str0[5];
        string potential      = str0[6];

        // ---- target / residual masses ----
        if (!getIso(isoA, ZA, AA) || !getIso(isoB, ZB, AB)) {
            std::fprintf(stderr, "\033[31m  ===> Error! cannot decode target/residual. \033[0m\n");
            continue;
        }
        double mA, mB, ma, mb;
        if (!massEx(ZA, AA, mA) || !massEx(ZB, AB, mB) ||
            !massEx(Za, Aa, ma) || !massEx(Zb, Ab, mb)) {
            std::fprintf(stderr, "\033[31m  ===> Error! mass not found. \033[0m\n");
            continue;
        }

        // ---- A/Z balance ----
        if (ZA + Za != ZB + Zb || AA + Aa != AB + Ab) {
            std::fprintf(stderr, "\033[31m====> ERROR! A-number or Z-number not balanced. \033[0m\n");
            continue;
        }

        if (isTransferReaction && potential.length() != 2) {
            std::fprintf(stderr, "\033[31m====> ERROR! Potential input should be 2 characters! skipped. \033[0m\n");
            continue;
        }

        // ---- transferred orbital decode (transfer only) ----
        string node, jStr, lStr;
        int spdf = -1;
        if (isTransferReaction) {
            node = orbital.substr(0, 1);

            if (std::abs(Aa - Ab) == 1) {       // single-nucleon transfer
                lStr = orbital.substr(1, 1);
                jStr = orbital.substr(2);
                spdf = getLValue(lStr);

                // parity: gsParity * stateParity must equal (-1)^l
                int expectedParity = gsParityAValue * parityValue;
                int calculatedParity = ((spdf % 2) == 0) ? 1 : -1;
                if (expectedParity != calculatedParity) {
                    std::fprintf(stderr, "\033[31m ===> skipped. Parity mismatch. (expected %d but got %d) \033[0m\n",
                                 expectedParity, calculatedParity);
                    continue;
                }

                // angular momentum: spin in |gsSpin - j| .. gsSpin + j
                float jValue = getSpinValue(jStr);
                bool spinMatch = false;
                for (float s = std::fabs(gsSpinAValue - jValue); s <= (gsSpinAValue + jValue); s += 1.0f) {
                    if (std::fabs(s - spinValue) < 1e-6) { spinMatch = true; break; }
                }
                if (!spinMatch) {
                    std::fprintf(stderr, "\033[31m ===> skipped. Angular momentum mismatch. \033[0m\n");
                    continue;
                }
            }

            if (std::abs(Aa - Ab) == 2) {       // two-nucleon transfer
                size_t posEq = orbital.find('=');
                lStr = orbital.substr(posEq + 1, 1);
                spdf = std::atoi(lStr.c_str());
            }

            if (std::abs(Aa - Ab) == 0) {
                std::fprintf(stderr, "\033[31m ===> skipped. p-n exchange reaction does not support. \033[0m\n");
            }

            if (spdf == -1) {
                std::fprintf(stderr, "\033[31m ===> skipped. Not recognized orbital-label. (user input : l=%s | %s) \033[0m\n",
                             lStr.c_str(), orbital.c_str());
                continue;
            }
        }

        // ---- beam energy (MeV total or MeV/u) ----
        int pos = (int)reactionEnergy.length() - 1;
        for (int i = pos; i >= 0; --i) {
            if (std::isdigit((unsigned char)reactionEnergy[i])) { pos = i; break; }
        }
        string unit = reactionEnergy.substr(pos + 1);
        int factor = 1;
        if (unit == "MeV/u") factor = Aa;
        double totalBeamEnergy = std::atof(reactionEnergy.substr(0, pos + 1).c_str()) * factor;

        double Qvalue = ma + mA - mb - mB;

        // =====================================================================
        // write Ptolemy input
        // =====================================================================
        ++numOfReaction;
        char hdr[256];

        if (!isTransferReaction) {
            // ----------------- elastic-ish -----------------
            if (std::atof(Ex.c_str()) == 0.0) {
                std::snprintf(hdr, sizeof hdr,
                    "$============================================ ELab=%5.2f(%s+%s)%s\n",
                    totalBeamEnergy, mass_a.c_str(), isoA.c_str(), potential.c_str());
                os << hdr;
                os << "reset\n";
                std::snprintf(hdr, sizeof hdr, "CHANNEL %s + %s\n", mass_a.c_str(), isoA.c_str()); os << hdr;
                os << "r0target\n";
                std::snprintf(hdr, sizeof hdr, "ELAB = %f\n", totalBeamEnergy); os << hdr;
                std::snprintf(hdr, sizeof hdr, "JBIGA=%s\n", gsSpinparityA.c_str()); os << hdr;
                string pot1Name = potential.substr(0, 1);
                std::snprintf(hdr, sizeof hdr, "$%s\n", potentialRef(pot1Name).c_str()); os << hdr;
                emitFullOM(os, callPotential(pot1Name, AA, ZA, totalBeamEnergy, Za));
                os << "ELASTIC SCATTERING\n";
                os << ";\n";
            } else {
                // ----------------- inelastic -----------------
                std::snprintf(hdr, sizeof hdr,
                    "$============================================ Ex=%s(%s+%s|%s%s)%s,ELab=%5.2f\n",
                    Ex.c_str(), mass_a.c_str(), isoA.c_str(), spin.c_str(), parity.c_str(),
                    potential.c_str(), totalBeamEnergy);
                os << hdr;
                os << "reset\n";
                std::snprintf(hdr, sizeof hdr, "REACTION: %s%s%s(%s%s %s) ELAB=%7.3f\n",
                    isoA.c_str(), reactionType.c_str(), isoB.c_str(),
                    spin.c_str(), parity.c_str(), Ex.c_str(), totalBeamEnergy);
                os << hdr;
                os << "PARAMETERSET ineloca2 r0target\n";
                std::snprintf(hdr, sizeof hdr, "JBIGA=%s\n", gsSpinparityA.c_str()); os << hdr;
                if (str0.size() >= 8) {
                    std::snprintf(hdr, sizeof hdr, "BETA=%s\n", str0[7].c_str()); os << hdr;
                }
                string pot1Name = potential.substr(0, 1);
                string pot1Ref = potentialRef(pot1Name);
                std::snprintf(hdr, sizeof hdr, "$%s\n", pot1Ref.c_str()); os << hdr;
                os << "INCOMING\n";
                emitReducedOM(os, callPotential(pot1Name, AA, ZA, totalBeamEnergy, Za));
                os << ";\n";
                os << "OUTGOING\n";
                std::snprintf(hdr, sizeof hdr, "$%s\n", pot1Ref.c_str()); os << hdr;
                emitReducedOM(os, callPotential(pot1Name, AA, ZA, totalBeamEnergy - std::atof(Ex.c_str()), Za));
                os << ";\n";
            }
        } else {
            // ----------------- transfer -----------------
            std::snprintf(hdr, sizeof hdr,
                "$============================================ Ex=%s(%s)%s\n",
                Ex.c_str(), orbital.c_str(), potential.c_str());
            os << hdr;
            os << "reset\n";
            std::snprintf(hdr, sizeof hdr, "REACTION: %s%s%s(%s%s %s) ELAB=%7.3f\n",
                isoA.c_str(), reactionType.c_str(), isoB.c_str(),
                spin.c_str(), parity.c_str(), Ex.c_str(), totalBeamEnergy);
            os << hdr;

            // ---- projectile (light) side ----
            if (std::abs(numNucleonsTransfer) == 1) {
                if (Aa <= 2 && Za <= 1 && Ab <= 2 && Zb <= 1) {   // d or p in
                    os << "PARAMETERSET dpsb r0target \n";
                    os << "lstep=1 lmin=0 lmax=30 maxlextrap=0 asymptopia=50 \n";
                    os << "\n";
                    os << "PROJECTILE \n";
                    os << "wavefunction av18 \n";
                    os << "r0=1 a=0.5 l=0 rc0=1.2\n";
                }
                if ((3 <= Aa && Aa <= 4) || (3 <= Ab && Ab <= 4)) {   // 3He/t/a involved
                    os << "PARAMETERSET alpha3 r0target \n";
                    os << "lstep=1 lmin=0 lmax=30 maxlextrap=0 asymptopia=50 \n";
                    os << "\n";
                    os << "PROJECTILE \n";
                    os << "wavefunction phiffer \n";
                    if (Za + Zb == 2) {   // (t,d) or (d,t)
                        os << "nodes=0 l=0 jp=1/2 spfacp=1.30 v=172.88 r=0.56 a=0.69 param1=0.64 param2=1.15 rc=2.0\n";
                    }
                    if (Za + Zb == 3) {   // (3He,d) or (d,3He)
                        os << "nodes=0 l=0 jp=1/2 spfacp=1.31 v=179.94 r=0.54 a=0.68 param1=0.64 param2=1.13 rc=2.0\n";
                    }
                    if (Ab == 4) {
                        os << "nodes=0 l=0 jp=1/2 spfacp=1.61 v=202.21 r=.93 a=.66 param1=.81 param2=.87 rc=2.0 $ rc=2 is a quirk\n";
                    }
                }
            } else if (std::abs(numNucleonsTransfer) == 2) {   // two-nucleon transfer
                os << "PARAMETERSET alpha3 r0target\n";
                os << "lstep=1 lmin=0 lmax=30 maxlextrap=0 ASYMPTOPIA=40\n";
                os << "\n";
                os << "PROJECTILE\n";
                os << "wavefunction phiffer\n";
                os << "L = 0  NODES=0 R0 = 1.25  A = .65     RC0 = 1.25\n";
            }
            os << ";\n";

            // ---- target (bound state) ----
            os << "TARGET\n";
            // Separation energies of the residual B from the mass table:
            //   Sn = dM(Z,A-1) + dM(n)  - dM(Z,A)
            //   Sp = dM(Z-1,A-1) + dM(1H) - dM(Z,A)
            // If a neighbouring nuclide is absent from PtolemyCpp's table (it has
            // less exotic-isotope coverage than Cleopatra's mass20.txt), that
            // channel is treated as bound (sentinel) rather than spuriously
            // unbound, so the deck structure stays correct. The original would
            // set the binding hack (E=-.2) whenever a mass was missing.
            const double SENTINEL = 1.0e9;
            double nThreshold = SENTINEL; { double mBm1; if (massEx(ZB, AB - 1, mBm1)) nThreshold = mBm1 + exN - mB; }
            double pThreshold = SENTINEL; { double mBp1; if (massEx(ZB - 1, AB - 1, mBp1)) pThreshold = mBp1 + exH - mB; }
            double ExEnergy = std::atof(Ex.c_str());
            bool isAboveThreshold = false;
            if (ExEnergy > nThreshold || ExEnergy > pThreshold) {
                isAboveThreshold = true;
                std::fprintf(stderr, "         Ex = %.3f MeV is above thresholds; Sp = %.3f MeV, Sn = %.3f MeV\n",
                             ExEnergy, pThreshold, nThreshold);
            }

            if (std::abs(Aa - Ab) == 1) {
                std::snprintf(hdr, sizeof hdr, "JBIGA=%s\n", gsSpinparityA.c_str()); os << hdr;
                if (isAboveThreshold) {
                    std::snprintf(hdr, sizeof hdr, "nodes=%s l=%d jp=%s E=-.2 $node is n-1, set binding 200 keV \n",
                                  node.c_str(), spdf, jStr.c_str());
                } else {
                    std::snprintf(hdr, sizeof hdr, "nodes=%s l=%d jp=%s $node is n-1 \n",
                                  node.c_str(), spdf, jStr.c_str());
                }
                os << hdr;
                os << "r0=1.25 a=.65 \n";
                os << "vso=6 rso0=1.10 aso=.65 \n";
                os << "rc0=1.3 \n";
            }
            if (std::abs(Aa - Ab) == 2) {
                std::snprintf(hdr, sizeof hdr, "JBIGA=%s\n", gsSpinparityA.c_str()); os << hdr;
                if (isAboveThreshold) {
                    std::snprintf(hdr, sizeof hdr, "nodes=%s L=%d E=-.2 $node is n-1, binding by 200 keV \n",
                                  node.c_str(), spdf);
                } else {
                    std::snprintf(hdr, sizeof hdr, "nodes=%s L=%d  $node is n-1 \n", node.c_str(), spdf);
                }
                os << hdr;
            }
            os << ";\n";

            // ---- potentials ----
            string pot1Name = potential.substr(0, 1);
            std::snprintf(hdr, sizeof hdr, "INCOMING $%s\n", potentialRef(pot1Name).c_str()); os << hdr;
            emitFullOM(os, callPotential(pot1Name, AA, ZA, totalBeamEnergy, Za));
            os << ";\n";

            string pot2Name = potential.substr(1, 1);
            std::snprintf(hdr, sizeof hdr, "OUTGOING $%s\n", potentialRef(pot2Name).c_str()); os << hdr;
            double eBeam = totalBeamEnergy + Qvalue - std::atof(Ex.c_str());
            emitFullOM(os, callPotential(pot2Name, AB, ZB, eBeam, Zb));
            os << ";\n";
        }

        char ang[128];
        std::snprintf(ang, sizeof ang, "anglemin=%f anglemax=%f anglestep=%f\n", angMin, angMax, angStep);
        os << ang;
        os << ";\n";
    }

    std::fprintf(stderr, "================= end of input. Number of Reaction : %d \n", numOfReaction);
    os << "end $================================== end of input\n";

    return os.str();
}
