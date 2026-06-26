// CrossSection_anapow.cpp — ANAPOW: computes and prints analyzing powers from F amplitudes.

#include "CrossSectionCalc.h"
#include "MathTables.h"
#include "Reaction.h"
#include "Constants.h"
#include <cstdio>
#include <vector>

// ============================================================================
// ============================================================================

void CrossSectionCalc::analyzingPower(double angleMin, double angleMax, double angleStep, int jA, int jB,
            int jResidual, int jBigB, int nSpline, int lxParity,
            int debugSwitch, double eLab, const char* channelName,
            const std::vector<double>& F_in, const int* tocsBase)
{
    auto& outputInLab    = reaction_.flags.outputInLab;

    // labCm: 4-char label
    const char* labCm = (outputInLab != 0) ? "LAB." : "C.M.";

    int angleCount = (int)((angleMax - angleMin) / angleStep + 1.5);
    int numAnalyzingPowers  = ((jA + 2) * jA) / 2;

    std::vector<double> analyzingPowerVector(angleCount * numAnalyzingPowers, 0.0);

    // qCount for K=1 (first iteration of outer K loop uses qCount = jA+1-QMIN)
    int qCount = jA + (jA + 1) % 2;   // = jA + ((jA+1) % (2)); max qCount across all K
    int allocSize = nSpline * nSpline * qCount;
    std::vector<double> analyzingPowerCoefVector(allocSize, 0.0);    // was IAPCO
    allocSize = qCount + 2 * jA + 4;  // +1 extra: termK needs jA+1 elements (index 0..jA)
    std::vector<double> analyzingPowerScratchVector(allocSize, 0.0);  // was IAPSCR
    // fProductVector, crossSectionVector: local scratch — not in FRENAM or NAMLOC
    std::vector<double> fProductVector(angleCount, 0.0);
    std::vector<double> crossSectionVector(angleCount, Constants::smlNum);  // INIT8 → fill with smlNum
    int kCross = std::min(jA, 2);

    // Get allocator pointers (setLog equivalent)
    setLog(1);

    // Sub-array offsets within analyzingPowerScratchVector (0-indexed)
    // analyzingPowerScratchVector[0..qCount-1]=COEP, [qCount]=COET, [qCount+1..qCount+1+jA]=coEx, [qCount+jA+3..]=termK
    double* coEpPointer   = analyzingPowerScratchVector.data();
    double* coEtPointer   = analyzingPowerScratchVector.data() + qCount;
    double* coExPointer   = analyzingPowerScratchVector.data() + qCount + 1;
    double* termKPointer  = analyzingPowerScratchVector.data() + qCount + jA + 3;
    // tocsBase is 1-based (tocsBase[1] = first int).
    const int* tocsPointer = tocsBase;
    // F_in is the caller-owned F-amplitude buffer.
    // Same 1-based offset semantics so fPointer[2*angIndex + 2*angleCount*kOff] addresses F_in[2*angIndex + 2*angleCount*(kOff-1) - 2].
    const double* fBase = F_in.empty() ? nullptr : F_in.data() - 2 * angleCount - 2;

    // Outer K loop
    for (int k = 1; k <= jA; k++) {
        int qMin = k % 2;
        qCount = k + 1 - qMin;

        // ASSIGN 180/170 equivalent: use REAL part for even K (qMin=0), IMAG for odd K (qMin=1)
        bool useReal = (qMin == 0);

        // 0-based column offset into analyzingPowerVector for this K:
        // Fortran: LAPOW_offset = ((K+1)*(K-1) + qMin)/2 - qMin  (column-pair index)
        int apowKBase = angleCount * (((k+1)*(k-1) + qMin)/2 - qMin);

        // Call MUELCO to get angle-independent coefficients
        // ITOC: 0-based int* where ITOC[0] = first int of slot.
        // COEF → ALLOC starting at Z(IAPCO) (1-based)
        muElCoupling(k, jA, jB, jResidual, jBigB,
               qCount, nSpline, lxParity,
               const_cast<int*>(tocsBase + 1),  // ITOC
               analyzingPowerCoefVector.data(),      // COEF
               coEpPointer,        // COEP
               coEtPointer,        // COET
               coExPointer,        // coEx
               termKPointer,       // TERMK
               debugSwitch);

        // Loops over kOffP (prime) and kOff
        for (int kOffP = 1; kOffP <= nSpline; kOffP++) {
            if (tocsPointer[4*kOffP] < 0) continue;
            for (int kOff = 1; kOff <= nSpline; kOff++) {
                if (tocsPointer[4*kOff] < 0) continue;

                bool hasComputedProducts = false;

                // analyzingPowerCoefVector index: offset = -qMin + qCount*(kOff-1 + nSpline*(kOffP-1))
                int apcoBase = -qMin + qCount * ((kOff-1) + nSpline*(kOffP-1));

                // For diagonal k=kCross: accumulate cross section
                if (kOff == kOffP && k == kCross) {
                    double tCs = 2.0;
                    int mX = (tocsPointer[4*kOff - 3] + tocsPointer[4*kOff - 2] + 1) / 2;
                    if (mX == 0) tCs = 1.0;
                    {
                    const double* fPointer = fBase + 2*angleCount*kOff;  // fPointer[2*angIndex] = F_in[2*angleCount*(kOff-1) + 2*angIndex - 2]
                    for (int angIndex = 1; angIndex <= angleCount; angIndex++) {
                        double fr = fPointer[2*angIndex], fi = fPointer[2*angIndex + 1];
                        fProductVector[angIndex-1] = fr*fr + fi*fi;
                        crossSectionVector[angIndex-1] += tCs * fProductVector[angIndex-1];
                    }
                    }
                }

                // Print APCO coefficients if debugging
                if (debugSwitch) {
                    std::printf(" APCO%2d%2d", kOff, kOffP);
                    for (int q = qMin; q <= k; q++)
                        std::printf("%12.5g", analyzingPowerCoefVector[apcoBase + q]);
                    std::printf("\n");
                }

                // Loop over Q, looking for non-zero coefficient
                for (int q = qMin; q <= k; q++) {
                    if (analyzingPowerCoefVector[apcoBase + q] == 0.0) continue;

                    // Compute F*F'* products at all angles if needed
                    if (!hasComputedProducts) {
                        {
                        const double* fiPointer = fBase + 2*angleCount*kOff;
                        const double* fjPointer = fBase + 2*angleCount*kOffP;
                        for (int angIndex = 1; angIndex <= angleCount; angIndex++) {
                            double ir = fiPointer[2*angIndex], ii = fiPointer[2*angIndex+1];
                            double jr = fjPointer[2*angIndex], ji = fjPointer[2*angIndex+1];
                            if (useReal)
                                fProductVector[angIndex-1] = ir*jr + ii*ji;
                            else
                                fProductVector[angIndex-1] = jr*ii - ir*ji;
                        }
                        }
                        hasComputedProducts = true;
                    }

                    // Accumulate analyzing powers at this Q (0-based vector view)
                    int apowColOffset = apowKBase + angleCount * q;       // 0-based start of this Q column
                    double tAp = analyzingPowerCoefVector[apcoBase + q];
                    double* lapPointer = analyzingPowerVector.data() + apowColOffset;  // lapPointer[angIndex-1] = analyzingPowerVector[apowColOffset + angIndex - 1]
                    for (int angIndex = 1; angIndex <= angleCount; angIndex++)
                        lapPointer[angIndex-1] += tAp * fProductVector[angIndex-1];
                } // q loop
            } // kOff loop
        } // kOffP loop
    } // k loop

    // Divide analyzing powers by cross section
    {
        double* apow0Pointer = analyzingPowerVector.data();
        int idx = 0;
        for (int ii = 1; ii <= numAnalyzingPowers; ii++) {
            for (int angIndex = 1; angIndex <= angleCount; angIndex++) {
                apow0Pointer[idx] /= crossSectionVector[angIndex-1];
                idx++;
            }
        }
    }

    // analyzingPowerCoefVector, analyzingPowerScratchVector, fProductVector, crossSectionVector auto-destruct here (no NFREE needed)


    // Determine block size for blank-line spacing
    double angleBlock = 300.0;
    if (angleCount > 40) {
        if (angleStep <= 0.201) angleBlock = 1.0;
        else if (angleStep < 0.99) angleBlock = 5.0;
        else angleBlock = 10.0;
    }
    int numBlocks = (int)(angleBlock / angleStep + 0.5);
    const int lineMax = 57;

    // Build column label lists
    int colK[13] = {}, colQ[13] = {};
    char colIFlag[13] = {};   // ' ' or 'i'

    // 0-based column iterator into analyzingPowerVector
    int colStart = 0;
    double* apow1Pointer = analyzingPowerVector.data();  // 0-based: apow1Pointer[ii - colStart + angIndex - 1]
    int kStart  = 1;

    // Outer print loop: one page per group of K values
    while (kStart <= jA) {
        int kEnd = (kStart == 5) ? std::min(jA, 6) : std::min(jA, 4);
        if (kStart  > 5) kEnd = kStart;

        // Build column labels for K = kStart..kEnd
        int columnCount = 0;
        for (int k = kStart; k <= kEnd; k++) {
            int qMin = k % 2;
            for (int q = qMin; q <= k; q++) {
                columnCount++;
                colK[columnCount] = k;
                colQ[columnCount] = q;
                colIFlag[columnCount] = (qMin != 0) ? 'i' : ' ';
            }
        }

        int colEnd  = colStart + (columnCount - 1) * angleCount;
        int lineCount = 1000;  // Force header on first line

        // Loop through angles
        for (int angIndex = 1; angIndex <= angleCount; angIndex++) {
            double angle = angleMin + angleStep * (angIndex - 1);
            double anglePrint = angle;
            if (outputInLab != 0 && angle > reaction_.kin.aBar)
                anglePrint = 2.0 * reaction_.kin.aBar - angle;

            // Print header if needed
            if (lineCount >= lineMax) {
                std::printf("1%54sP T O L E M Y\n", "");
                std::printf("%10sANALYZING POWERS FOR THE  >>>> %-8s <<<<  CHANNEL\n",
                            "", channelName);
                // '0' = blank line
                std::printf("0%.45s%7.2f MEV     %.65s\n",
                            &reaction_.reactStr[1], eLab, &reaction_.header[1]);
                // '0' = blank line
                std::printf("0  %-4.4s ", labCm);
                for (int i = 1; i <= columnCount; i++)
                    std::printf("      %cT%d%d", colIFlag[i], colK[i], colQ[i]);
                std::printf("\n");
                std::printf("  ANGLE\n\n");
                lineCount = 8;
            }

            std::printf(" %6.2f  ", anglePrint);
            for (int ii = colStart; ii <= colEnd; ii += angleCount)
                std::printf("%10.5f", apow1Pointer[ii - colStart + angIndex - 1]);
            std::printf("\n");
            lineCount++;

            // Blank line after angleBlock-th angle
            double angleModulo = std::fmod(angle + angleStep + 1.e-7, angleBlock);
            if (angleModulo <= 1.e-5) {
                std::printf("\n");  // blank line
                lineCount++;
                if (lineCount + numBlocks > lineMax) lineCount = 1000;
            }
        } // angle loop

        colStart = colEnd + angleCount;
        kStart  = kEnd + 1;
    } // page group loop
}
