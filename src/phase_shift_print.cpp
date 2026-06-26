// phase_shift_print.cpp — PHSPRT: computes and prints elastic/reaction S-matrix
// magnitudes, phases, and unitarity.

#include "CrossSectionCalc.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include "Reaction.h"

static double modPi(double x) {
    const double subHalf = 0.4999999990;
    const double TWOPI  = 6.2831853071795864769;
    return x - TWOPI * std::trunc(x / TWOPI + (x >= 0 ? subHalf : -subHalf));
}
static double unModulo(double phasP, double phase) {
    return phasP + modPi(phase - phasP);
}

// Fortran G15.5: Gw.d uses F(w-4).k+4blanks when 0.1<=|val|<10^d, else Ew.d
static void print_G15_5(double val) {
    double absVal = std::fabs(val);
    if (absVal == 0.0 || (absVal >= 0.1 && absVal < 1.0e5)) {
        int e   = (absVal >= 1.0) ? (int)std::floor(std::log10(absVal)) + 1 : 0;
        int decimals = 5 - e; if (decimals < 0) decimals = 0;
        std::printf("%11.*f    ", decimals, val);
    } else {
        int exponent = (int)std::floor(std::log10(absVal)) + 1;
        double mantissa = val / std::pow(10.0, (double)exponent);
        long long m = (long long)(std::fabs(mantissa) * 1e5 + 0.5);
        if (m >= 100000LL) { m /= 10; exponent++; }
        std::printf("   %c0.%05lldE%+03d", (val < 0) ? '-' : ' ', m, exponent);
    }
}

// Fortran F10.5 with overflow (Fortran prints "**********" when value won't fit)
static void print_F10_5(double val) {
    if (val >= 10000.0 || val < -1000.0)
        std::printf("**********");
    else
        std::printf("%10.5f", val);
}

// ============================================================================
void CrossSectionCalc::phasePrint(bool printSwitch)
{
    int lMin   = reaction_.angMom.lMin;
    int lMax   = reaction_.angMom.lMax;
    int lSkip  = reaction_.inelastic.lSkip;
    int nSpl   = reaction_.inelastic.nSpl;
    int lInMax = reaction_.internalState.lInMax;
    double eLab = reaction_.energies.eLab;


    // Elastic unitarity accumulator (lMax+1 elements, indexed 0..lMax)
    std::vector<double> elasticUnitarityVector(lMax + 2, 0.0);
    double* elasticUnitarityPointer = elasticUnitarityVector.data();  // elasticUnitarityPointer[0..lMax] valid

    // Magnitude/phase pointer arrays (1-based: chanNumber=1..3)
    // chanNumber=2,3 (SINMAG/SOUTMAG/SINPHASE/SOUTPHAS) are local-only vectors.
    std::vector<float> elasticSMagVectors[4];  // [2] = incoming elastic, [3] = outgoing elastic
    std::vector<float> elasticSPhaseVectors[4];

    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------
    for (int chanNumber = 1; chanNumber <= 3; chanNumber++) {
        bool isElasticChannel = (chanNumber > 1);
        int  channelIndex = chanNumber - 1;

        float* smagPointer;
        float* sphasePointer;
        if (!isElasticChannel) {
            // Sized to cover the NaN scrub range in xSection (ntot_sp = nSpl*(lInMax+1));
            // phasePrint writes are at 0-based indices [0..] via smagPointer = data()
            // (the loop below uses magIndex = i-1 for this reaction branch).
            int cap = reaction_.inelastic.nSpl * (lInMax + 1) + 1;
            if (cap < reaction_.inelastic.liloSize + 1) cap = reaction_.inelastic.liloSize + 1;
            reaction_.inelastic.smagArr.assign(cap, 0.0f);
            reaction_.inelastic.sphaseArr.assign(cap, 0.0f);
            smagPointer   = reaction_.inelastic.smagArr.data();
            sphasePointer = reaction_.inelastic.sphaseArr.data();
        } else {
            int cap = reaction_.distortedWave.channel[channelIndex].nJStates * (lMax + 2);
            elasticSMagVectors[chanNumber].assign(cap, 0.0f);
            elasticSPhaseVectors[chanNumber].assign(cap, 0.0f);
            smagPointer   = elasticSMagVectors[chanNumber].data();
            sphasePointer = elasticSPhaseVectors[chanNumber].data();
        }

        // Channel parameters
        int    stateCount, lLast, lBase;
        double phaseMultiplier;

        float*  smatReactionPointer = nullptr;   // float* for reaction IRDINT
        double* smatElasticPointer = nullptr;   // double* for elastic idx_smat
        // tocPointer is 1-based base for the TOCS array — pool slot (inelastic)
        // or toceArr.data()-1 (per-channel elastic).
        const int* tocPointer = nullptr;
        if (!isElasticChannel) {
            stateCount   = nSpl;   lLast = lInMax;  lBase = lMin;  phaseMultiplier = 1.0;
            tocPointer = reaction_.inelastic.tocsPointer;  // 1-based class-owned
            smatReactionPointer = reaction_.inelastic.sMatrixArr.data();  // 0-based (accessed [2*i-2]/[2*i-1]); was PoolView<float>(IRDINT).data()
        } else {
            stateCount   = reaction_.distortedWave.channel[channelIndex].nJStates;
            lLast  = lMax;    lBase = 0;       phaseMultiplier = 0.5;
            tocPointer = reaction_.distortedWave.channel[channelIndex].toceArr.data() - 1;
            // smatElasticPointer is a 0-based ptr into smatArr (accessed [2*i-2]/[2*i-1],
            // matching the non-elastic smatReactionPointer branch above).
            smatElasticPointer = reaction_.distortedWave.channel[channelIndex].smatArr.data();
        }

        for (int kOffset = 1; kOffset <= stateCount; kOffset++) {
            if (tocPointer[4*kOffset] < 0) continue;

            double phasP = 0.0;
            // Backward sweep: lFwd=lBase..lLast → li=lLast..lBase
            for (int lFwd = lBase; lFwd <= lLast; lFwd++) {
                int li = lLast + lBase - lFwd;
                int i  = (li - lBase) * stateCount + kOffset;

                double sReal, sImag;
                if (!isElasticChannel) {
                    sReal = (double)smatReactionPointer[2*i - 2];
                    sImag = (double)smatReactionPointer[2*i - 1];
                } else {
                    sReal = smatElasticPointer[2*i - 2];
                    sImag = smatElasticPointer[2*i - 1];
                }

                // smagArr/sphaseArr (reaction branch) are 0-based: magIndex = i-1.
                // elasticSMagVectors/elasticSPhaseVectors keep their 1-based layout
                // ([0] unused) via .data(): magIndex = i.
                int magIndex = isElasticChannel ? i : i - 1;
                float sMag = (float)std::sqrt(sReal*sReal + sImag*sImag);
                smagPointer[magIndex] = sMag;
                if (sMag == 0.0f) continue;

                double phase = std::atan2(sImag, sReal);
                phase = unModulo(phasP, phase);
                phasP = phase;
                sphasePointer[magIndex] = (float)(phaseMultiplier * phase);

                if (chanNumber == 2)  // incoming elastic: accumulate unitarity
                    elasticUnitarityPointer[li] += sReal*sReal + sImag*sImag;
            }
        }
    }

    if (!printSwitch) return;

    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------
    {
        // toceArr replaces pool slots; 0-based bases (accessed [4*k-3]/[4*k-4]).
        const int* toceInPointer = reaction_.distortedWave.channel[1].toceArr.data();
        const int* toceOutPointer = reaction_.distortedWave.channel[2].toceArr.data();
        const int* tocsPointer = reaction_.inelastic.tocsPointer;  // 1-based class-owned

        // PoolView-based access; ISMG[1]/ISPH[1] stay pool, [2,3] are local vectors.
        float* smagInPointer  = elasticSMagVectors[2].data();
        float* smagOutPointer  = elasticSMagVectors[3].data();
        float* sphaseInPointer = elasticSPhaseVectors[2].data();
        float* sphaseOutPointer = elasticSPhaseVectors[3].data();

        double* sigmaInPointer = reaction_.distortedWave.channel[1].sigmaArr.data();
        double* sigmaOutPointer = reaction_.distortedWave.channel[2].sigmaArr.data();
        float* unitrPointer = reaction_.inelastic.unitrArr.data() - lMin;  // unitrPointer[li] valid for li=lMin..jMost

        int lines  = std::max(reaction_.distortedWave.channel[1].nJStates, reaction_.distortedWave.channel[2].nJStates);
        int lParityOffset  = ((reaction_.boundState.vertex[1].lBound + reaction_.boundState.vertex[2].lBound) % 2 != 0) ? 1 : 0;

        // -----------------------------------------------------------
        // Elastic partial wave S-matrix table
        // -----------------------------------------------------------
        int  lineCount = 100;
        bool isAtPageTop  = false;

        for (int li = lMin; li <= lMax; li += lSkip) {

            if (lineCount + lines > 58) {
                std::printf("1%56sP T O L E M Y\n", "");
                std::printf("%38sELASTIC PARTIAL WAVE S-MATRIX ELEMENTS AND UNITARITY\n", "");
                std::printf("0%.45s%-6s%7.2f MEV     %.65s\n",
                            &reaction_.reactStr[1], "ELAB =", eLab, &reaction_.header[1]);
                std::printf("\n");

                std::printf("%15sINCOMING ELASTIC %25sUNITARITY%26sOUTGOING ELASTIC\n",
                            "", "", "");
                std::printf("    L   L' LX  MAGNITUDE   PHASE  COULOMB  "
                            "     ELASTIC  REACTION  RESIDUAL"
                            "        L   L' LX  MAGNITUDE   PHASE  COULOMB \n");
                std::printf("\n");

                lineCount = 9;
                isAtPageTop  = true;
            }

            // Unitarity values for this li
            double reactionUnitarity = (double)unitrPointer[li];
            double elasticUnitarity = elasticUnitarityPointer[li];
            double residual  = 1.0 - elasticUnitarity - reactionUnitarity;
            double coulombPhaseIn = sigmaInPointer[li];
            int    lasO   = li + lParityOffset;
            double coulombPhaseOut = sigmaOutPointer[lasO];

            // Blank line before L block (except right after page header)
            if (!isAtPageTop && (lines > 1 || li % 10 == 0)) {
                std::printf(" \n");
                lineCount++;
            }

            // Inner loop: kOffset=1..lines, iterating incoming/outgoing together
            int kIn = 0, kOut = 0;
            for (int kOffset = 1; kOffset <= lines; kOffset++) {
                int jump = 0;
                int lxIn = 0, loIn = 0;
                double magIn = 0.0, phaseIn = 0.0;
                int lxOut = 0, loOut = 0;
                double magOut = 0.0, phaseOut = 0.0;

                // Incoming: advance kIn, skip zeros
                {
                    int kIn0 = kIn;
                    for (;;) {
                        kIn++;
                        if (kIn > reaction_.distortedWave.channel[1].nJStates) { kIn = kIn0; break; }
                        double mag = (double)smagInPointer[reaction_.distortedWave.channel[1].nJStates*li + kIn];
                        if (mag != 0.0) {
                            magIn = mag;
                            phaseIn = (double)sphaseInPointer[reaction_.distortedWave.channel[1].nJStates*li + kIn];
                            lxIn   = toceInPointer[4*kIn - 3];
                            loIn   = li + toceInPointer[4*kIn - 4];
                            jump += 2;
                            break;
                        }
                    }
                }

                // Outgoing: advance kOut, skip zeros
                {
                    int kOut0 = kOut;
                    for (;;) {
                        kOut++;
                        if (kOut > reaction_.distortedWave.channel[2].nJStates) { kOut = kOut0; break; }
                        double mag = (double)smagOutPointer[reaction_.distortedWave.channel[2].nJStates*lasO + kOut];
                        if (mag != 0.0) {
                            magOut = mag;
                            phaseOut = (double)sphaseOutPointer[reaction_.distortedWave.channel[2].nJStates*lasO + kOut];
                            lxOut   = toceOutPointer[4*kOut - 3];
                            loOut   = lasO + toceOutPointer[4*kOut - 4];
                            jump += 1;
                            break;
                        }
                    }
                }

                if (kOffset == 1) {
                    std::printf(" %4d%5d%3d%10.6f%9.3f%9.3f    ",
                                li, loIn, lxIn, magIn, phaseIn, coulombPhaseIn);
                    print_F10_5(elasticUnitarity);
                    print_F10_5(reactionUnitarity);
                    print_F10_5(residual);
                    std::printf("     %4d%5d%3d%10.6f%9.3f%9.3f\n",
                                lasO, loOut, lxOut, magOut, phaseOut, coulombPhaseOut);
                } else {
                    if (jump == 0) break;  // no more data for this li
                    if (jump == 1) {
                        std::printf("%80s%4d%5d%3d%10.6f%9.3f\n",
                                    "", lasO, loOut, lxOut, magOut, phaseOut);
                    } else if (jump == 2) {
                        // incoming only (no outgoing data)
                        // 1X + I4+I5+I3+F10.6+F9.3, then record ends (no trailing 48X)
                        std::printf(" %4d%5d%3d%10.6f%9.3f\n",
                                    li, loIn, lxIn, magIn, phaseIn);
                    } else {
                        // jump==3: both incoming and outgoing
                        // 1X + I4+I5+I3+F10.6+F9.3 + 48X + I4+I5+I3+F10.6+F9.3
                        std::printf(" %4d%5d%3d%10.6f%9.3f%48s%4d%5d%3d%10.6f%9.3f\n",
                                    li, loIn, lxIn, magIn, phaseIn,
                                    "", lasO, loOut, lxOut, magOut, phaseOut);
                    }
                }
                lineCount++;
                isAtPageTop = false;
            }
        }

        // -----------------------------------------------------------
        // Reaction S-matrix table (non-CC: NCHNDF=2, NCHN=2 only)
        // -----------------------------------------------------------

        {
            // Count valid kOffset
            int validPwCount = 0;
            for (int kOffset = 1; kOffset <= nSpl; kOffset++)
                if (tocsPointer[4*kOffset] >= 0) validPwCount++;
            int setCount = (validPwCount + 4) / 5;

            // 0-based pointers into class-owned smagArr/sphaseArr (accessed [flatIndex-1]; was ISMG[1]/ISPH[1] pool slots).
            float* smagReactionPointer   = reaction_.inelastic.smagArr.data();
            float* sphaseReactionPointer = reaction_.inelastic.sphaseArr.data();

            int kFirst = 1;
            for (int setIndex = 1; setIndex <= setCount; setIndex++) {
                // Collect up to 5 valid partial waves
                int jTts[6]={}, jTps[6]={}, lXs[6]={}, lDeltas[6]={}, kOffsets[6]={};
                int pwCount = 0;
                for (int kOffset = kFirst; kOffset <= nSpl && pwCount < 5; kOffset++) {
                    int jTt = tocsPointer[4*kOffset];
                    if (jTt < 0) continue;
                    pwCount++;
                    jTts[pwCount]   = jTt;
                    lDeltas[pwCount]  = tocsPointer[4*kOffset - 3];
                    lXs[pwCount]  = tocsPointer[4*kOffset - 2];
                    jTps[pwCount]   = tocsPointer[4*kOffset - 1];
                    kOffsets[pwCount] = kOffset;
                    if (pwCount == 5) { kFirst = kOffset + 1; break; }
                    if (kOffset == nSpl) kFirst = kOffset + 1;
                }

                // Force page break if this set won't fit on current page
                {
                    double est = 1.1 * (double)(lMax - lMin + lSkip) / lSkip;
                    if (lineCount + 8 + est > 58) lineCount = 100;
                }

                bool isAtPageTopReaction = true;

                for (int li = lMin; li <= lMax; li += lSkip) {

                    bool needHeader = (lineCount > 58) || isAtPageTopReaction;
                    if (needHeader) {
                        if (lineCount > 58) {
                            std::printf("1%47sP T O L E M Y\n", "");
                            std::printf(" %.45s%-6s%7.2f MEV     %.65s\n",
                                        &reaction_.reactStr[1], "ELAB =", eLab, &reaction_.header[1]);
                        } else {
                            std::printf(" \n");
                        }
                        std::printf("0 L IN%13sS-MATRICES FOR CHANNEL%3d"
                                    "  AND PARTIAL WAVES LABELED BY ( JP, JT, LX, L(OUT)-L(IN) )\n",
                                    "", 2);
                        std::printf("0    ");
                        for (int i = 1; i <= pwCount; i++)
                            std::printf("     (%2d/2%3d/2%3d%4d )",
                                        jTps[i], jTts[i], lXs[i], lDeltas[i]);
                        std::printf("\n");
                        std::printf("     ");
                        for (int i = 1; i <= pwCount; i++)
                            std::printf("        |S|        PHASE");
                        std::printf("\n");
                        std::printf(" \n");
                        isAtPageTopReaction = true;
                        lineCount += 8;
                        if (lineCount > 58) lineCount = 9;
                    }

                    // Blank line at multiples of 10 (not right after header)
                    if (li % 10 == 0 && !isAtPageTopReaction) {
                        std::printf(" \n");
                        lineCount++;
                    }

                    std::printf("%5d", li);
                    for (int i = 1; i <= pwCount; i++) {
                        int flatIndex = nSpl * (li - lMin) + kOffsets[i];
                        print_G15_5((double)smagReactionPointer[flatIndex - 1]);
                        std::printf("%9.3f", (double)sphaseReactionPointer[flatIndex - 1]);
                    }
                    std::printf("\n");
                    isAtPageTopReaction = false;
                    lineCount++;
                }
            }
        }
    }

}

// ---------------------------------------------------------------------------
// Backward-compat free function wrapper — calls the class method
// ---------------------------------------------------------------------------

