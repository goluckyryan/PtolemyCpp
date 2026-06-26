// set_channels.cpp — SETCHN: sets up parameters for a given channel.

#include "Reaction.h"
#include <cstdio>
#include "Constants.h"

// aMs(r) — 0-based pointer into the massesArr std::array<double, 5>
// One of three TUs with positional aMs(*this)[i] callers (others:
// parameters.cpp, probe_print.cpp); each anon-namespace copy is independent,
// so this one is 0-based (callers shift -1) while parameters.cpp stays 1-based.
namespace {
double* aMs(Reaction& r) { return r.masses.massesArr.data(); }  // 0-based
}

void Reaction::setChannel(int channelIndex, int& returnCode)
{
    // ELECTRON MASS FOR CONVERTING ATOMIC TO NUCLEAR MASSES
    static const double electronMass = 0.5110034;

    returnCode = 1;

    int jProjOld = (int)this->angMom.jProj;
    int jOld = (int)this->angMom.J;
    int lOld = this->angMom.L;
    int parityOld = this->angMom.parity;

    // Declare all variables at function scope to avoid goto-crosses-init
    double signMx = this->internalState.undefValue;
    int signZx = NOTDEF_INT;
    int projIndex = 0, tgtIndex = 0;
    double tMp = 0, tMt = 0, cmToLab = 1.0;
    double amxcProj = 0, amxcTarget = 0;

    // wasSet tracks things set to other things; CLRCHN uses it.
    // Tight range matches the actually-live slots (see InternalState.h).
    for (int i = 1; i <= 29; i++)
        this->internalState.wasSet[i] = 0;

    // boundIndex = 0 means not bound state, > 0 means bound state
    int boundIndex = 1;
    // Wrap pre-kinematic setup in do-while-false so the 5 former `goto L310;`
    // sites become `break;` to skip straight to the kinematic block below.
    do {
    if (channelIndex != 5) {
    boundIndex = 0;
    if (channelIndex == 6) break;

    // Compute mass excess = mass excess (ground state) + E*
    for (int i = 1; i <= 5; i++) {
        if (this->masses.amxcs[i] != this->internalState.undefValue) continue;
        if (this->masses.amxcgs[i] == this->internalState.undefValue) continue;
        this->masses.amxcs[i] = this->masses.amxcgs[i] + this->energies.exs[i];
    }

    // Find ZX, mX and the signed ZX, mX

    // If charge transfer is 0 or mX = 0 (inelastic), no ambiguity
    if (std::fabs(this->masses.massesArr[4]) < 0.3) signMx = 0;
    if (this->charges.zArray[5] == 0) signZx = 0;

    for (int i = 1; i <= 3; i += 2) {
        if (this->charges.zArray[i] != NOTDEF_INT && this->charges.zArray[i+1] != NOTDEF_INT) {
            signZx = (this->charges.zArray[i] - this->charges.zArray[i+1]) * (2 - i);
            if (this->charges.zArray[5] == NOTDEF_INT) this->charges.zArray[5] = std::abs(signZx);
        }
        if (aMs(*this)[i - 1] == this->internalState.undefValue || aMs(*this)[i] == this->internalState.undefValue)
            continue;
        signMx = (aMs(*this)[i - 1] - aMs(*this)[i]) * (2 - i);
        if (aMs(*this)[4] == this->internalState.undefValue) aMs(*this)[4] = std::fabs(signMx);
    }
    if (std::fabs(signMx) < 0.3) signMx = 0;

    // Supply missing mass or Z values if possible
    for (int i = 1; i <= 4; i++) {
        int mirrorIndex = i - 1 + 2 * ((i) % (2));
        int kSign = std::abs(5 - 2*i) - 2;
        if (signZx != NOTDEF_INT && this->charges.zArray[i] == NOTDEF_INT && this->charges.zArray[mirrorIndex] != NOTDEF_INT) {
            this->charges.zArray[i] = this->charges.zArray[mirrorIndex] + kSign * signZx;
        }
        if (signMx == this->internalState.undefValue) continue;
        if (aMs(*this)[i - 1] != this->internalState.undefValue) continue;
        if (aMs(*this)[mirrorIndex - 1] == this->internalState.undefValue) continue;
        aMs(*this)[i - 1] = aMs(*this)[mirrorIndex - 1] + kSign * signMx;
    }


    // Determine which is P and T
    this->internalState.stripPickup = (signZx != NOTDEF_INT) ? signZx : NOTDEF_INT;
    if (signMx != this->internalState.undefValue) this->internalState.stripPickup = (int)signMx;
    if (this->internalState.stripPickup == NOTDEF_INT)
        std::printf("\n***** WARNING:  NOT ENOUGH DATA TO DETERMIN IF"
                    " REACTION IS STRIPPING OR PICKUP - STRIPPING ASSUMED.\n");
    if (this->internalState.stripPickup == NOTDEF_INT) this->internalState.stripPickup = 1;
    if (this->internalState.stripPickup != 0)
        this->internalState.stripPickup = (this->internalState.stripPickup >= 0) ? 1 : -1;

    if (channelIndex >= 3) {
        // Scattering states - no composite
        boundIndex = 0;
        projIndex = channelIndex - 2;
        tgtIndex = projIndex + 2;
    } else {
        // Bound state, projectile is always X
        projIndex = 5;
        if (this->internalState.stripPickup < 0) {
            boundIndex = channelIndex + 1;
            tgtIndex = boundIndex + 2*channelIndex - 3;
        } else {
            tgtIndex = channelIndex + 1;
            boundIndex = tgtIndex + 2*channelIndex - 3;
        }
    }

    // Define undefined P and T stuff
    if (this->charges.zProj == NOTDEF_INT) this->charges.zProj = this->charges.zArray[projIndex];
    if (this->charges.zTarget == NOTDEF_INT) this->charges.zTarget = this->charges.zArray[tgtIndex];
    if (this->masses.massProj == this->internalState.undefValue) this->masses.massProj = aMs(*this)[projIndex - 1];
    if (this->masses.massTgt == this->internalState.undefValue) this->masses.massTgt = aMs(*this)[tgtIndex - 1];
    if (this->angMom.spinProj == this->internalState.notDefSentinel) this->angMom.spinProj = this->angMom.js[projIndex];
    if (this->angMom.spinTarget == this->internalState.notDefSentinel) this->angMom.spinTarget = this->angMom.js[tgtIndex];

    // Copy E* for P and T
    if (this->energies.exsPt[1] == this->internalState.undefValue) this->energies.exsPt[1] = this->energies.exs[projIndex];
    if (this->energies.exsPt[2] == this->internalState.undefValue) this->energies.exsPt[2] = this->energies.exs[tgtIndex];
    if (this->masses.amxgPt[1] == this->internalState.undefValue) this->masses.amxgPt[1] = this->masses.amxcgs[projIndex];
    if (this->masses.amxgPt[2] == this->internalState.undefValue) this->masses.amxgPt[2] = this->masses.amxcgs[tgtIndex];
    if (this->angMom.parityPt[1] == 0) this->angMom.parityPt[1] = this->angMom.parities[projIndex];
    if (this->angMom.parityPt[2] == 0) this->angMom.parityPt[2] = this->angMom.parities[tgtIndex];
    this->internalState.lSpcPt2 = this->internalState.lSpecs[tgtIndex];
    this->internalState.nodePt2 = this->internalState.nodesP[tgtIndex];

    // For scattering allow ranges of L, jProj, and J
    if (boundIndex == 0) break;

    // Following page or so is for bound states only
    if (this->angMom.J == this->internalState.notDefSentinel) this->angMom.J = this->angMom.js[boundIndex];
    if (this->angMom.parity == 0) this->angMom.parity = this->angMom.parities[boundIndex];
    if (this->angMom.L == NOTDEF_INT) this->angMom.L = this->internalState.lSpecs[boundIndex];
    if (this->angMom.nNodes == NOTDEF_INT) this->angMom.nNodes = this->internalState.nodesP[boundIndex];
    }

    // The following code is also used when not doing a DWBA
    if (this->angMom.L == NOTDEF_INT) this->angMom.L = this->internalState.lSpcPt2;
    if (this->angMom.nNodes == NOTDEF_INT) this->angMom.nNodes = this->internalState.nodePt2;

    // If some angular momenta are zero then we can define others
    if (this->angMom.jProj == this->internalState.notDefSentinel && this->angMom.spinTarget == 0) this->angMom.jProj = this->angMom.J;
    if (this->angMom.jProj == this->internalState.notDefSentinel && this->angMom.J == 0) this->angMom.jProj = this->angMom.spinTarget;
    if (this->angMom.jProj == this->internalState.notDefSentinel && this->angMom.L == 0) this->angMom.jProj = this->angMom.spinProj;
    if (this->angMom.jProj == this->internalState.notDefSentinel && this->angMom.spinProj == 0 &&
        this->angMom.L != NOTDEF_INT) this->angMom.jProj = 2*this->angMom.L;
    if (this->angMom.J == this->internalState.notDefSentinel && this->angMom.spinTarget == 0) this->angMom.J = this->angMom.jProj;
    if (this->angMom.L == NOTDEF_INT) {
        if (this->angMom.jProj != this->internalState.notDefSentinel && this->angMom.spinProj == 0) this->angMom.L = (int)this->angMom.jProj / 2;
        if (this->angMom.spinProj != this->internalState.notDefSentinel && this->angMom.jProj == 0) this->angMom.L = (int)this->angMom.spinProj / 2;

        // If only two possible L's, parity can tell us
        if (this->angMom.L == NOTDEF_INT) {
            int i = this->angMom.parity * this->angMom.parityPt[1] * this->angMom.parityPt[2];
            if (i == 0) break;
            i = (i + 3) / 2;
            if (this->angMom.jProj == this->internalState.notDefSentinel || this->angMom.spinProj == this->internalState.notDefSentinel) break;
            int lMin = std::abs((int)this->angMom.jProj - (int)this->angMom.spinProj) / 2;
            lMin = lMin + ((i + lMin) % (2));
            int lMax = ((int)this->angMom.jProj + (int)this->angMom.spinProj) / 2;
            if (lMax - lMin > 1) break;
            this->angMom.L = lMin;
        }
    }

    // L is defined; we may be able to patch up a parity
    { int i = +1;
    if (((this->angMom.L) % (2)) != 0) i = -1;
    if (this->angMom.parity == 0) this->angMom.parity = this->angMom.parityPt[1] * this->angMom.parityPt[2] * i;
    if (this->angMom.parityPt[1] == 0) this->angMom.parityPt[1] = i * this->angMom.parity * this->angMom.parityPt[2];
    if (this->angMom.parityPt[2] == 0) this->angMom.parityPt[2] = i * this->angMom.parity * this->angMom.parityPt[1]; }

    } while (false);  // end of pre-kinematic skip-section

    // Setup kinematic projectile and target masses
    amxcProj = 0;
    if (this->masses.amxgPt[1] != this->internalState.undefValue)
        amxcProj = (this->masses.amxgPt[1] + this->energies.exsPt[1]) / Constants::amu_MeV;
    amxcTarget = 0;
    if (this->masses.amxgPt[2] != this->internalState.undefValue)
        amxcTarget = (this->masses.amxgPt[2] + this->energies.exsPt[2]) / Constants::amu_MeV;

    // MP and MT must be separately defined
    if (this->masses.massProj == this->internalState.undefValue || this->masses.massTgt == this->internalState.undefValue) {
        std::printf("\n**** ERROR:  BOTH MP AND MT ARE REQUIRED.\n");
        returnCode = 0;
    } else {

        // Now the kinematic masses
        tMp = this->masses.massProj;
        tMt = this->masses.massTgt;
        // MASTYP == 0 always; the original `if (MASTYP == 0)` guard
        // dropped along with the dead == 2 branch above.
        { int intMass = (int)tMp;
        if (intMass == tMp) tMp = tMp + amxcProj - this->charges.zProj*(electronMass/Constants::amu_MeV); }
        { int intMass = (int)tMt;
        if (intMass == tMt) tMt = tMt + amxcTarget - this->charges.zTarget*(electronMass/Constants::amu_MeV); }

        // Define reduced mass if necessary
        if (this->masses.aM == this->internalState.undefValue)
            this->masses.aM = Constants::amu_MeV * tMp * tMt / (tMp + tMt);

        // Find eLab to eCm conversion
        this->internalState.ratMass = tMp / tMt;
        cmToLab = 1 + this->internalState.ratMass;
    }

    // Find E if necessary.
    // The "Incoming scattering state" block runs for channelIndex==6
    // unconditionally, and for channelIndex==3 (or any other channelIndex<6 not in the
    // switch) when E is undefined; cases 1/2/4/5 handle their own E setup.
    bool runIncoming = (channelIndex == 6);
    // The "E MUST BE DEFINED" bailout (printf + returnCode=0) recurs across the
    // four channelIndex cases below; share it via one helper. The trailing
    // break/err handling stays at each call site.
    auto reportEUndefined = [&] { std::printf("\n**** E MUST BE DEFINED.\n"); returnCode = 0; };
    if (channelIndex < 6 && this->energies.E == this->internalState.undefValue) {
        switch (channelIndex) {
            case 1:
            case 2:
                // Bound state; first try Q and then mass excess
                if (this->energies.Q != this->internalState.undefValue && this->internalState.eBnds[3-channelIndex] != this->internalState.undefValue) {
                    this->energies.E = (3 - 2*channelIndex) * this->internalState.stripPickup * this->energies.Q + this->internalState.eBnds[3-channelIndex];
                    break;
                }
                if (this->masses.amxcs[projIndex] == this->internalState.undefValue || this->masses.amxcs[tgtIndex] == this->internalState.undefValue
                    || this->masses.amxcs[boundIndex] == this->internalState.undefValue) { reportEUndefined(); break; }
                this->energies.E = this->masses.amxcs[boundIndex] - this->masses.amxcs[projIndex] - this->masses.amxcs[tgtIndex];
                break;
            case 4:
                // Outgoing scattering state; apply Q to eCm
                if (this->energies.eCm == this->internalState.undefValue) { reportEUndefined(); break; }
                if (this->energies.Q == this->internalState.undefValue) {
                    // Attempt to get Q from total mass excesses
                    bool err = false;
                    for (int i = 1; i <= 4; i++) {
                        if (this->masses.amxcs[i] == this->internalState.undefValue) { reportEUndefined(); err = true; break; }
                    }
                    if (err) break;
                    this->energies.Q = this->masses.amxcs[1] - this->masses.amxcs[2] + this->masses.amxcs[3] - this->masses.amxcs[4];
                }
                this->energies.E = this->energies.eCm + this->energies.Q;
                break;
            case 5: reportEUndefined(); break;
            case 3: default: runIncoming = true; break;
        }
    }

    if (runIncoming) {
        // Incoming scattering state; eCm or eLab applies
        if (this->energies.eCm != this->internalState.undefValue) this->energies.E = this->energies.eCm;
        if (this->energies.eLab != this->internalState.undefValue) {
            this->energies.E = this->energies.eLab / cmToLab;
        } else if (this->energies.E == this->internalState.undefValue) {
            std::printf("\n**** E, ECM OR ELAB MUST BE DEFINED.\n");
            returnCode = 0;
        }
    }

    // At last E is defined. For incoming wave define eCm, eLab
    if (channelIndex == 3 || channelIndex == 6) {
        this->energies.eCm = this->energies.E;
        this->internalState.wasSet[16] = (this->energies.eLab == this->internalState.undefValue && channelIndex == 6);
        this->energies.eLab = this->energies.E * cmToLab;
    }

    // For outgoing channel, define Q
    if (channelIndex == 4 && this->energies.eCm != this->internalState.undefValue)
        this->energies.Q = this->energies.E - this->energies.eCm;

    // Store P and T stuff back into 1-5 stuff
    if (channelIndex <= 4) {
        aMs(*this)[projIndex - 1] = this->masses.massProj;
        aMs(*this)[tgtIndex - 1] = this->masses.massTgt;
        this->charges.zArray[projIndex] = this->charges.zProj;
        this->charges.zArray[tgtIndex] = this->charges.zTarget;
        this->angMom.js[projIndex] = this->angMom.spinProj;
        this->angMom.js[tgtIndex] = this->angMom.spinTarget;
        this->angMom.parities[projIndex] = this->angMom.parityPt[1];
        this->angMom.parities[tgtIndex] = this->angMom.parityPt[2];
        if (boundIndex != 0) {
            this->angMom.js[boundIndex] = this->angMom.J;
            this->angMom.parities[boundIndex] = this->angMom.parity;
        }
        if (this->spec.spam != this->internalState.undefValue) {
            if (channelIndex == 1) this->spec.specAmpProj = this->spec.spam;
            if (channelIndex == 2) this->spec.specAmpTgt = this->spec.spam;
        }
    }

    this->internalState.wasSet[12] = (jProjOld != (int)this->angMom.jProj);
    this->internalState.wasSet[13] = (jOld != (int)this->angMom.J);
    this->internalState.wasSet[14] = (lOld != this->angMom.L);
    this->internalState.wasSet[29] = (this->angMom.parity != parityOld);

    // Set asymptopia if needed
    if (this->integrationGrid.asymptopia == this->internalState.undefValue) {
        this->internalState.wasSet[15] = 1;
        this->integrationGrid.asymptopia = this->integrationGrid.boundAsy;
        if (boundIndex == 0) this->integrationGrid.asymptopia = std::fabs(this->integrationGrid.scatAsy);
    }
}

// ---------------------------------------------------------------------------
// Backward-compat free function wrapper — calls the class method
// ---------------------------------------------------------------------------
