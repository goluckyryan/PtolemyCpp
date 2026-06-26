// interpolation.cpp — LINTRP: interpolation/extrapolation in L of reaction S-matrices.

#include "ptolemy_types.h"
#include "MathTables.h"
#include "Timing.h"
#include "l_extrapolation.h"
#include "math/angular_momentum_coeff.h"
#include "CrossSectionCalc.h"
#include "math/continued_fraction.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Reaction.h"
#include "Constants.h"

// interpCfrac — continued-fraction S-matrix interpolation; defined at end of
// file, used only here.
static void interpCfrac(int jProj, int jT, int lx, int lDelta,
            int nLValues, int* lIs, int lStep, int lxMax,
            double* sR, double* sI, int sDimension,
            int lStart, int lMax, float* sOut4,
            int lBase2,
            double& bigBoy, int& lPeak, double* complexL, double* complexS);

// Read the (jT, jProj, lx, lDelta) channel quad from tocsPointer at kOffset.
static void readChannelQuad(const int* tocsPointer, int kOffset,
                            int& jT, int& jProj, int& lx, int& lDelta) {
    jT     = tocsPointer[4 * kOffset];
    jProj  = tocsPointer[4 * kOffset - 1];
    lx     = tocsPointer[4 * kOffset - 2];
    lDelta = tocsPointer[4 * kOffset - 3];
}

void CrossSectionCalc::linterp()
{
    const auto PI     = Constants::PI;
    const auto BIGLOG = Constants::BIGLOG;
    auto& lMin   = reaction_.angMom.lMin;
    auto& lMax   = reaction_.angMom.lMax;
    auto& maxLExtrap = reaction_.integrationGrid.maxLExtrap;
    auto& printLevel = reaction_.flags.printLevel;
    auto& excitationType = reaction_.flags.excitationType;
    auto& lInMax = reaction_.internalState.lInMax;
    auto& undefValue  = reaction_.internalState.undefValue;
    auto& akIn    = reaction_.kin.akIn;
    auto& akOut    = reaction_.kin.akOut;
    auto& lOutMax = reaction_.kin.lOutMax;
    auto& lxMax  = reaction_.inelastic.lxMax;
    auto& nSpl  = reaction_.inelastic.nSpl;
    auto& nLValues = reaction_.inelastic.nLValues;
    auto& liFit  = reaction_.inelastic.liFit;
    auto& liFitIndex = reaction_.inelastic.liFitIndex;
    auto& nolFit = reaction_.inelastic.nolFit;

    // LINTRP's GRIDCM overlay (IBINDX/MBINDX/NBINDX/NMFFAC) was the only consumer
    // of the int* reinterpret-cast view above; all four are now class-owned vectors.

    // Local variables
    int mWork = 12;
    double tStart, aP, bP, cP;
    bool isTransferReaction, isExtrapolating, debugSwitch;
    int verbosity, liPrevious;
    std::vector<double> cfLValuesVector;    // CMPLXLS: continued-fraction complex L values (2*nLValues)
    std::vector<double> cfSValuesVector; // CMPLXIS: continued-fraction complex S(L) (2*nLValues)
    std::vector<double> workVector;  // LWORK: extrapolation fit workspace (mWork*nAspli)
    std::vector<double> fitValuesVector;  // FITLS: L values for extrapolation fit (2*nolFit)
    std::vector<double> fitSValuesVector; // FITIS: |S(L)| values for extrapolation fit (2*nolFit)
    double* workPointer  = nullptr;   // raw pointer into workVector (0-based)
    double* valsPointer  = nullptr;   // raw pointer into fitValuesVector (0-based)
    double* xintsPointer = nullptr;   // raw pointer into fitSValuesVector (0-based)
    int jMost, lMaxL;  // KSMAT/NCHNDF dropped: KSMAT was 0 always.
    int lnsMatrix, liloSize;
    int* lisPointer = nullptr;  // ptr to LIS array (set after LLIS assignment)
    int liStart;
    int lx, lDelta, jProj, jT, kOffset, lPeak;
    // factor / elPhaseIn init to 0 silences -Wmaybe-uninitialized for the
    // isTransferReaction / !isExtrapolating branches — factor's writer (~line 447) lives inside
    // `if (!isTransferReaction)`, every reader (~line 578) inside the same guard.
    // elPhaseIn's writer (~line 458) is in `if (isExtrapolating)`, reader (~line 526)
    // in the matching `if (isExtrapolating)` block. GCC can't see the guards line up.
    double bigBoy, lMaxDouble;
    double factor = 0, elPhaseIn = 0;
    bool isExtendedLRange;

    // (FOURSW=true at the sole caller made the SOUT branch dead).

    // Extrapolation type names
    const char* xWords[6][2] = {
        {"WOODS-SA", "XON"},
        {"POWER-LA", "W"},
        {"POWER-LA", "W 2"},
        {"WKB POWE", "R-LAW"},
        {"LEXTRAP4", " "},
        {"LEXTRAP5", " "}
    };

    tStart = second();

    // Define things that might never get defined
    aP = 0; bP = 0; cP = 0;

    isTransferReaction = (reaction_.internalState.stripPickup != 0);
    verbosity = printLevel % 10;
    debugSwitch = ((printLevel / 100) % 10) >= 4;

    isExtrapolating = (maxLExtrap > 0);

    // Print header
    std::printf("1%58sP T O L E M Y\n"
                "%47sINTERPOLATION AND EXTRAPOLATION IN L\n"
                "0%.45sELAB =%7.2f MEV     %.65s\n\n",
                "", "",
                &reaction_.reactStr[1], reaction_.energies.eLab, &reaction_.header[1]);

    int excitationIndex = excitationType - 1;
    if (excitationIndex < 0 || excitationIndex >= 6) excitationIndex = 0;
    std::printf("0USING %-8s%-8sEXTRAPOLATION FUNCTIONS.\n",
                xWords[excitationIndex][0], xWords[excitationIndex][1]);

    lMaxDouble = lMax;
    liPrevious = lInMax;
    if (isTransferReaction) lInMax = lMax;

    // These 2 work areas are used for continued-fraction interpolation.
    cfLValuesVector.assign(2 * nLValues, 0.0);
    cfSValuesVector.assign(2 * nLValues, 0.0);

    int* tocsPointer = nullptr;  // set after each LTOCS assignment
    if (isExtrapolating) {
    // Extrapolation work arrays
    workVector.assign(mWork * reaction_.inelastic.nAspli, 0.0);   // INIT8 equivalent: all zeros
    fitValuesVector.assign(2 * nolFit, 0.0);
    fitSValuesVector.assign(2 * nolFit, 0.0);
    workPointer  = workVector.data();
    valsPointer  = fitValuesVector.data();
    xintsPointer = fitSValuesVector.data();

    // Allocation is finished for now
    tocsPointer = reaction_.inelastic.tocsPointer;  // 1-based class-owned
    lisPointer = reaction_.inelastic.lisArr.data();  // 1-based (lisPointer[liIndex] valid for 1..nLValues)
    // valsPointer/xintsPointer indexing is inlined as `i - 1` (0-based) at the
    // direct-fit sites and `nolFit - 1 + i` for the +nolFit shift.

    std::printf("\nEXTRAPOLATION AND CONTINUED FRACTION OVERLAP"
                " REGION:%4d < LI <%4d\n\n", liFit, lMax);

    // ========================================================================
    // Possibly setup for multiple channels
    // ========================================================================
    }  // end if (isExtrapolating) — extrapolation setup

    jMost = 0;

    // ========================================================================
    // PASS 1: Determine lInMax and size requirements
    // ========================================================================
    {
        lMaxL = lMax;

        // skip extrapolation fitting if not enabled
        if (isExtrapolating) {

        // Pass 1 extrapolation fitting loop (isExtrapolating=true path)
        tocsPointer = reaction_.inelastic.tocsPointer;  // 1-based class-owned
        lisPointer = reaction_.inelastic.lisArr.data();  // refresh
        for (kOffset = 1; kOffset <= nSpl; kOffset++) {
            readChannelQuad(tocsPointer, kOffset, jT, jProj, lx, lDelta);
            if (jT < 0) continue;
            double* delsRPointer = reaction_.inelastic.smatrArr.data() + (kOffset - 1);
            double* delsIPointer = reaction_.inelastic.smatiArr.data() + (kOffset - 1);

            // Locate L-peak and check for decaying form
            int lPeakL = lMin;
            bigBoy = 0;
            bool cannotExtrapolate = false;
            double sizeLast = 0, signRealL = 0, signImagL = 0;

            for (int liIndex = 1; liIndex <= nLValues; liIndex++) {
                int li = lisPointer[liIndex];
                int lo = li + lDelta;
                if (lo + li < lx) continue;
                int lOffset = (liIndex - 1) * nSpl;
                double re = delsRPointer[lOffset];
                double him = delsIPointer[lOffset];
                double size = re * re + him * him;
                if (li >= liFit) {
                    double signReal = copysign(1.0, re);
                    double signImag = copysign(1.0, him);
                    if (li > liFit) {
                        if (sizeLast > size) {
                            // OK - decaying
                        } else {
                            cannotExtrapolate = true;
                            double sMag = sqrt(size);
                            double sMagLast = sqrt(sizeLast);
                            std::printf("\n**** ERROR:  CANNOT EXTRAPOLATE FOR "
                                        "CHANNEL%3d  "
                                        "(JP,JT,LX,LO-LI) = (%3d/2%3d/2%2d%3d) ****\n"
                                        " ****%9sSINCE FOR LI =%4d"
                                        " |S| EXCEEDS PREVIOUS:%11.4G%11.4G ****\n",
                                        2, jProj, jT, lx, lDelta, "", li, sMag, sMagLast);
                        }
                        if (signReal != signRealL || signImag != signImagL) {
                            std::printf("0**** WARNING:  PHASE OF S(L) FOR "
                                        "CHANNEL%3d  "
                                        "(JP,JT,LX,LO-LI,LI) = (%3d/2%3d/2%3d %3d %3d"
                                        ") HAS FLUCTUATING SIGN\n",
                                        2, jProj, jT, lx, lDelta, li);
                        }
                    }
                    sizeLast = size;
                    signRealL = signReal;
                    signImagL = signImag;
                }
                size = (2 * li + 1) * size;
                if (size >= bigBoy) {
                    lPeakL = li;
                    bigBoy = size;
                }
            }

            if (cannotExtrapolate) continue;

            // Find weeBoy
            double weeBoy = 1.0e+06;
            for (int liIndex = 1; liIndex <= nLValues; liIndex++) {
                int li = lisPointer[liIndex];
                int lo = li + lDelta;
                if (lo + li < lx) continue;
                if (li > lPeakL) continue;
                int lOffset = (liIndex - 1) * nSpl;
                double size = delsRPointer[lOffset] * delsRPointer[lOffset]
                            + delsIPointer[lOffset] * delsIPointer[lOffset];
                size = (2 * li + 1) * size;
                if (size < weeBoy) weeBoy = size;
            }

            weeBoy = 0.2 * sqrt(weeBoy);
            bigBoy = reaction_.integrationGrid.dwCutoff * sqrt(bigBoy);
            weeBoy = std::min(weeBoy, bigBoy);

            // Fit extrapolation function
            for (int liIndex = liFitIndex; liIndex <= nLValues; liIndex++) {
                int i = liIndex - liFitIndex + 1;
                valsPointer[i - 1] = lisPointer[liIndex];  // R6              // LLVALS+i = i-1
                valsPointer[nolFit - 1 + i] = lisPointer[liIndex];  // R6     // LLVAL2+i = nolFit-1+i
                int lOffset = (liIndex - 1) * nSpl;
                double F = sqrt(delsRPointer[lOffset] * delsRPointer[lOffset]
                              + delsIPointer[lOffset] * delsIPointer[lOffset]);
                xintsPointer[i - 1] = F;                              // LXINTS+i = i-1
                F = atan2(delsIPointer[lOffset], delsRPointer[lOffset]);
                xintsPointer[nolFit - 1 + i] = F;                    // LXINT2+i = nolFit-1+i
            }

            double aa, width, flCrit, barL, barA, b, barC, chiSq;
            int convergenceCode;
            lxTrp1(excitationType, nolFit, convergenceCode, verbosity, valsPointer - 1,     // 1-based: ptr[1]=fitValuesVector[0]
                   xintsPointer - 1, flCrit, aa, width,                // 1-based: ptr[1]=fitSValuesVector[0]
                   barL, barA, b, barC, lMaxDouble, chiSq, lx, lDelta, jProj, jT);

            if (convergenceCode < 0) {
                std::printf("\n  **** ERROR RETURN FROM L-EXTRAP SUBROUTINE\n"
                            "   CHANNEL%3d  "
                            "       EXTRAPOLATION SUPPRESSED FOR (JP,JT,LX,LO-LI)=("
                            "%3d/2%3d/2%4d%4d) \n", 2, jProj, jT, lx, lDelta);
                continue;
            }
            if (convergenceCode != 0) {
                std::printf("\n  SIMPLE EXP OR POWER  USED FOR |S(L)| "
                            "   CHANNEL%3d  "
                            "(JP,JT,LX,LO-LI)=(%3d/2%3d/2%4d%4d)\n",
                            2, jProj, jT, lx, lDelta);
            }

            // Phase extrapolation for power-law
            if (excitationType >= 2) {
                // LLVAL2+1 = nolFit, pointing to second half of fitValuesVector
                // 1-based: ptr[1]=fitValuesVector[nolFit]
                lxTrp1(3, nolFit, convergenceCode, verbosity, valsPointer + nolFit - 1,
                       xintsPointer + nolFit - 1, flCrit, aa, width,
                       barL, cP, bP, aP, lMaxDouble, chiSq, lx, lDelta, jProj, jT);
                if (convergenceCode < 0) {
                    std::printf("\n**** FOLLOWING REFERS TO PHASE EXTRAP.\n");
                    continue;
                }
            }

            // Find lInMax for this (lx, lDelta)
            int li = (int)lxTrpM(excitationType, barA, b, barL, lMaxDouble, weeBoy);
            if (debugSwitch) {
                std::printf(" FOR CHANNEL%3d    FOR 2*JP, 2*JT, LX, LO-LI =%3d%3d%3d%3d   "
                            "ESTIMATED LARGEST SIGNIFICANT LI IS%4d\n",
                            2, jProj, jT, lx, lDelta, li);
            }
            lMaxL = std::max(lMaxL, li);

            // Save extrapolation parameters into workPointer (0-based)
            int lArea = mWork * (kOffset - 1);  // workPointer base for this kOffset
            workPointer[lArea]      = aa;
            workPointer[lArea + 1]  = b;
            workPointer[lArea + 2]  = barL;
            workPointer[lArea + 3]  = flCrit;
            workPointer[lArea + 4]  = aP;
            workPointer[lArea + 5]  = cP;
            workPointer[lArea + 6]  = barA;
            workPointer[lArea + 7]  = barC;
            workPointer[lArea + 8]  = li;
            workPointer[lArea + 11] = bP;
        } // end kOffset loop

        } // end if (isExtrapolating) — Pass 1 extrapolation

        // Compute size requirements
        if (lInMax <= lMax) lInMax = std::min(lMaxL, lMax + maxLExtrap);
        jMost = std::max(jMost, lInMax);
        lnsMatrix = nSpl * (lInMax - lMin + 1);

    } // end pass-1 channel block

    liloSize = lnsMatrix;
    reaction_.inelastic.liloSize = liloSize;

    // Allocate space for reaction S-matrices to XSECTN
    // sMatrixArr: 0-based float vector of size 2*liloSize (interleaved real/imag
    // pairs). Reader convention: rdIntPointer = sMatrixArr.data() + offset.
    reaction_.inelastic.sMatrixArr.assign(2 * liloSize, 0.0f);
    // unitrArr: indexed by [li - lMin] for li = lMin..jMost.
    reaction_.inelastic.unitrArr.assign(jMost - lMin + 2, 0.0f);

    isExtendedLRange = isTransferReaction || (jMost > liPrevious);

    // Extend elastic S-matrix arrays if needed
    lOutMax = jMost + lxMax;
    if (jMost > lMax) {
        // Elastic S-matrix extension (Born approximation)

        for (int i = 1; i <= 2; i++) {
            if (isExtendedLRange) reaction_.distortedWave.channel[i].sigmaArr.resize(lOutMax + 1, 0.0);
            int zeroStart = 2 * reaction_.distortedWave.channel[i].nJStates * (lMax + 1) + 1;
            int smatSize = 2 * reaction_.distortedWave.channel[i].nJStates * (lOutMax + 1);
            reaction_.distortedWave.channel[i].smatArr.resize(smatSize, 0.0);
            double* smatPointer = reaction_.distortedWave.channel[i].smatArr.data();  // 0-based: smatPointer[K-1] = slot[K]
            double* sigmaPointer  = reaction_.distortedWave.channel[i].sigmaArr.data();  // 0-based: [lMax] = element index lMax
            for (int ii = zeroStart; ii <= smatSize; ii++) {
                smatPointer[ii - 1] = 0;  // smatPointer[ii - 1] = smatArr[ii - 1] = slot[ii]
            }

            double coulombPhase = sigmaPointer[lMax];
            double e2 = sqrt(reaction_.kin.etaCh[i] * reaction_.kin.etaCh[i] + lMaxDouble * (lMaxDouble + 1));
            double akChannel = (i == 1) ? akIn : akOut;
            double exrFactor = 1.0e10;
            if (reaction_.kin.aScts[i] != undefValue) exrFactor = 1.0 / (akChannel * reaction_.kin.aScts[i]);
            double exiFactor = 1.0e10;
            if (reaction_.distortedWave.channel[i].aI  != undefValue) exiFactor = 1.0 / (akChannel * reaction_.distortedWave.channel[i].aI);
            if (reaction_.distortedWave.channel[i].aSi != undefValue) exiFactor = 1.0 / (akChannel * reaction_.distortedWave.channel[i].aSi);
            bool cannotExtrapolate = (exrFactor == 1.0e10 || exiFactor == 1.0e10);
            if (cannotExtrapolate) {
                std::printf("0**** WARNING:  ELASTIC S MATRIX FOR CHANNEL"
                            "%2d NOT EXTRAPOLATED BECAUSE A OR AI IS"
                            " NOT DEFINED.\n", i);
            }
            double deltaReal = 0, deltaImag = 0, smatMagnitude;
            if (!cannotExtrapolate) {
                int lOff = 2 * lMax * reaction_.distortedWave.channel[i].nJStates;
                smatMagnitude = sqrt(smatPointer[lOff] * smatPointer[lOff] + smatPointer[lOff+1] * smatPointer[lOff+1]);
                deltaImag = -0.5 * log(smatMagnitude);
                deltaReal = 0.5 * atan2(smatPointer[lOff+1], smatPointer[lOff]);
            }

            double smatReal = 1, smatImag = 0;
            for (int L = liFit - lxMax; L <= lOutMax; L++) {
                int lOff = 2 * reaction_.distortedWave.channel[i].nJStates * L;  // smatPointer[lOff] = pool[LSMAT + lOff]
                bool writeSmat = cannotExtrapolate;  // if bad, always write
                if (!cannotExtrapolate) {
                    double e1 = sqrt(reaction_.kin.etaCh[i] * reaction_.kin.etaCh[i] + (double)L * (L + 1));
                    double deltaRho = e1 - e2;
                    double realDamping = 0;
                    double x = deltaRho * exrFactor;
                    if (x < BIGLOG) realDamping = exp(-x);
                    double imagDamping = 0;
                    x = deltaRho * exiFactor;
                    if (x < BIGLOG) imagDamping = exp(-x);
                    double dampedDeltaReal = realDamping * deltaReal;
                    double dampedDeltaImag = imagDamping * deltaImag;
                    smatMagnitude = exp(-2.0 * dampedDeltaImag);
                    smatReal = smatMagnitude * cos(2.0 * dampedDeltaReal);
                    smatImag = smatMagnitude * sin(2.0 * dampedDeltaReal);
                    if (L > lMax) writeSmat = true;  // extrapolated L: write
                }
                if (writeSmat) {
                    smatPointer[lOff] = smatReal;
                    smatPointer[lOff + 1] = smatImag;
                    if (L > lMax && isExtendedLRange) {
                        coulombPhase = coulombPhase + atan(reaction_.kin.etaCh[i] / (double)L);
                        sigmaPointer[L] = coulombPhase;
                    }
                }
            }
        } // end i=1,2 loop
    } // end if jMost > lMax

    // ========================================================================
    // PASS 2: Interpolation and extrapolation
    // ========================================================================
    setLog(2 * (lOutMax + lxMax));
    lisPointer = reaction_.inelastic.lisArr.data();  // refresh for PASS 2
    tocsPointer = reaction_.inelastic.tocsPointer;  // 1-based class-owned
    // smatArr class-owned vectors; smat1Pointer/smat2Pointer are 0-based slot ptrs
    double* smat1Pointer = reaction_.distortedWave.channel[1].smatArr.data();
    double* smat2Pointer = reaction_.distortedWave.channel[2].smatArr.data();
    float* unitrPointer = reaction_.inelastic.unitrArr.data() - lMin;  // unitrPointer[li] valid for li=lMin..jMost
    double* bratPointer = nullptr;  // set if !isTransferReaction

    if (!isTransferReaction) {
        factor = sqrt(akIn * akOut / (reaction_.distortedWave.channel[1].Ecm * reaction_.distortedWave.channel[2].Ecm * PI));
        // bratPointer base = betas vector start, offset by -lxMin/2 so that
        // bratPointer[lx/2] == betas[(lx-lxMin)/2].
        bratPointer = reaction_.inelastic.poolBetas[2]->data() - reaction_.inelastic.lxMin / 2;
    }

    if (isExtrapolating) {

    // Phase extrapolation setup
    {
        int lMaxOffset = 2 * reaction_.distortedWave.channel[1].nJStates * lMax;
        elPhaseIn = 0.5 * atan2(smat1Pointer[lMaxOffset + 1], smat1Pointer[lMaxOffset]);
    }
    } // end if (isExtrapolating) — L700 block


    // (Pass-2 KSMAT2=0 reset also dropped — was the only writer of the
    // rdIntPointer offset, now the offset is just kOffset - 1.)

    {
        // loop over (jT, jProj, lx, lo-li)
        for (kOffset = 1; kOffset <= nSpl; kOffset++) {
            readChannelQuad(tocsPointer, kOffset, jT, jProj, lx, lDelta);
            if (jT < 0) continue;
            double* delsRPointer = reaction_.inelastic.smatrArr.data() + (kOffset - 1);
            double* delsIPointer = reaction_.inelastic.smatiArr.data() + (kOffset - 1);
            float* rdIntPointer = reaction_.inelastic.sMatrixArr.data() + 2 * (kOffset - 1);
            liStart = std::max(lMin, (lx - lDelta + 1) / 2);

            // Continued fraction interpolation
            interpCfrac(jProj, jT, lx, lDelta,
                   nLValues, lisPointer, reaction_.integrationGrid.lStep, lxMax,  // was = lisPointer
                   delsRPointer, delsIPointer, nSpl,
                   liStart, lMax, rdIntPointer,
                   lMin - 1,                          // ISDIM2 dropped — == ISDIM1 == nSpl
                   bigBoy, lPeak, cfLValuesVector.data(), cfSValuesVector.data());

            if (isTransferReaction && ((lPeak == lMin && lMin > 0) || lPeak >= lMax)) {
                std::printf("0**** WARNING:  FOR (2*JP,2*JT,LX,LO-LI) = (%3d%3d%3d%3d"
                            ")    lPeak = %3d IS AT END OF L-RANGE (FISHY) ***\n",
                            jProj, jT, lx, lDelta, lPeak);
            }

            if (isExtrapolating) {

            // Extrapolation (isExtrapolating=true path)
            {
                int lArea = mWork * (kOffset - 1);  // workPointer base for this kOffset
                double b    = workPointer[lArea + 1];
                double barL = workPointer[lArea + 2];
                double barA = workPointer[lArea + 6];
                double barC = workPointer[lArea + 7];
                workPointer[lArea + 9]  = lPeak;
                workPointer[lArea + 10] = sqrt(bigBoy);

                if (barA != 0) {

                aP = workPointer[lArea + 4];
                cP = workPointer[lArea + 5];
                bP = workPointer[lArea + 11];

                if (excitationType < 2) {
                    // Phase extrapolation using elastic phases
                    int lo = lMax + lDelta;
                    int lOffset = (nLValues - 1) * nSpl;
                    double phase = atan2(delsIPointer[lOffset], delsRPointer[lOffset]);
                    cP = 0;
                    if (fabs(phase) > 0.25 * PI) cP = copysign(0.5 * PI, phase);
                    if (fabs(phase) > 0.75 * PI) cP = copysign(PI, phase);
                    phase = phase - cP;
                    int smat2Offset = 2 * reaction_.distortedWave.channel[2].nJStates * lo;
                    double elPhaseOut = 0.5 * atan2(smat2Pointer[smat2Offset + 1], smat2Pointer[smat2Offset]);
                    double elPhase = elPhaseOut + elPhaseIn;
                    aP = phase / elPhase;
                    workPointer[lArea + 4] = aP;
                    workPointer[lArea + 5] = cP;
                }

                if (debugSwitch) {
                    std::printf("\n  * EXTRAPOLATED S-MATRICES WITH  CHANNEL =%3d   "
                                "(JP,JT,LX,LO-LI) = (%3d/2%3d/2%3d %3d) *\n",
                                2, jProj, jT, lx, lDelta);
                }

                for (int li = liFit; li <= lInMax; li++) {
                    int lo = li + lDelta;
                    double size;
                    lxTrp2(excitationType, barA, b, barC, barL, lMaxDouble, li, size);
                    if (size < Constants::smlNum) break;  // exit li extrapolation loop

                    double phase, re, him;
                    int smat2Offset = 2 * reaction_.distortedWave.channel[2].nJStates * lo;   // SMAT2 offset for lo
                    int smat1Offset = 2 * reaction_.distortedWave.channel[1].nJStates * li;  // SMAT1 offset for li
                    if (excitationType >= 2) {
                        double tempPhase;
                        lxTrp2(3, aP, bP, cP, lMaxDouble, lMaxDouble, li, tempPhase);
                        phase = tempPhase;
                    } else {
                        double elPhase = 0.5 * (atan2(smat2Pointer[smat2Offset+1], smat2Pointer[smat2Offset])
                                           + atan2(smat1Pointer[smat1Offset+1], smat1Pointer[smat1Offset]));
                        phase = cP + aP * elPhase;
                    }
                    re = size * cos(phase);
                    him = size * sin(phase);
                    if (fabs(re) < 1.0e-30) re = 0;
                    if (fabs(him) < 1.0e-30) him = 0;

                    // In overlap region, don't overwrite
                    if (li <= lMax) continue;

                    int lOffset = 2 * (li - lMin) * nSpl;
                    rdIntPointer[lOffset] = (float)re;
                    rdIntPointer[lOffset + 1] = (float)him;
                }
            } // end if (barA != 0)
            } // end extrapolation block
            } // end if (isExtrapolating)

            // Add Coulomb excitation for inelastic scattering
            if (!isTransferReaction) {

            {
                double scaledFactor = factor;
                if (reaction_.inelastic.densitySwitch) scaledFactor = sqrt(2.0) * scaledFactor;
                if ((lx % 4) < 2) scaledFactor = -scaledFactor;
                bool isImaginary = (lx % 2) == 0;
                // cl2ffPointer is 0-based ptr into cl2ffArr at offset (kOffset-1).
                // New: cl2ffArr[kOffset - 1 + lOffset] (0-based) = same slot[kOffset + lOffset] = old result.
                double* cl2ffPointer = reaction_.inelastic.cl2ffArr.data() + (kOffset - 1);

                // ICL2FF was allocated by COULST with n_spl*(lInMax+1) elements,
                // where lInMax is the COULST asymptotic limit. Loop to lInMax is safe.
                for (int li = liStart; li <= lInMax; li++) {
                    int lo = li + lDelta;
                    int lOffset = (li - lMin) * nSpl;
                    double coefficient = scaledFactor * bratPointer[lx / 2]
                                * fabs(clebschGordan(2 * li, 2 * lx, 0, 0, 2 * lo, 0));
                    double ffI = coefficient * cl2ffPointer[lOffset];
                    if (!isImaginary)
                        rdIntPointer[2 * lOffset] += (float)ffI;
                    if (isImaginary)
                        rdIntPointer[2 * lOffset + 1] += (float)ffI;
                }
            }
            } // end if (!isTransferReaction) — Coulomb excitation


            // Unitarity contributions
            for (int li = liStart; li <= lInMax; li++) {
                int rOffset = 2 * nSpl * (li - lMin);  // offset into rdIntPointer for this li
                float v0 = rdIntPointer[rOffset], v1 = rdIntPointer[rOffset + 1];
                unitrPointer[li] += v0 * v0 + v1 * v1;
            }
        } // end kOffset loop (989)

    } // end pass-2 channel block

    // Pass 3 would go here for coupled channels (CC only)
    // ...

    // Extrapolation summary
    if (isExtrapolating) {

    // Print extrapolation summary
    if (!isTransferReaction) {  // IBRNSB==0 always (permanently 0)
        std::printf("\n\n%32sFOLLOWING REFER TO THE \"NUCLEAR\" PART OF THE INELASTIC "
                    "AMPLITUDES.\n", "");
    }
    std::printf("\n\n%48sSUMMARY OF EXTRAPOLATION PARAMETERS\n"
                "\n%19s|S| = A/( C + EXP( (LI-LCRIT)/LWIDTH ) )%23s"
                "PHASE = CP + AP*(SUM OF ELASTIC PHASE-SHIFTS)\n", "", "", "");

    tocsPointer = reaction_.inelastic.tocsPointer;  // 1-based class-owned

    {

        for (kOffset = 1; kOffset <= nSpl; kOffset++) {
            readChannelQuad(tocsPointer, kOffset, jT, jProj, lx, lDelta);
            if (jT < 0) continue;
            int lArea = mWork * (kOffset - 1);  // workPointer base for this kOffset at start
            double width = workPointer[lArea + 1];
            if (excitationType == 1 && width != 0) width = 1.0 / width;
            int I_v      = (int)workPointer[lArea + 7];
            int lMaxLV  = (int)workPointer[lArea + 8];
            int lPeakV  = (int)workPointer[lArea + 9];
            std::printf("%5d%6d/2%3d/2%3d%5d%6d%13.3G%6d"
                        "%17.5G%4d%10.3f%12.3G%6d%13.3f%10.3f%10.3f\n",
                        2, jProj, jT, lx, lDelta, lPeakV,
                        workPointer[lArea + 10], lInMax,
                        workPointer[lArea],      I_v, width, workPointer[lArea + 3],
                        lMaxLV, workPointer[lArea + 5], workPointer[lArea + 4],
                        workPointer[lArea + 11]);
        }

    } // end pass-3 (extrapolation summary) channel block

    workVector.clear(); fitValuesVector.clear(); fitSValuesVector.clear();

    } // end if (isExtrapolating) — extrapolation summary

    // cleanup
    cfLValuesVector.clear(); cfSValuesVector.clear();

    {
        double elapsedTime = second() - tStart;
        std::printf("0 TIME FOR L-INTERPOLATIONS : %7.3f SECONDS\n\n", elapsedTime);
    }
    return;
}

// ---------------------------------------------------------------------------
// Backward-compat free function wrapper — calls the class method
// ---------------------------------------------------------------------------

// ============================================================================
// folded in from source_misc.cpp: recordPeak (static) + interpCfrac
// ============================================================================
// Track the channel L with the largest |S|^2 (size) seen so far.
static void recordPeak(double size, double& bigBoy, int& lPeak, int li) {
    if (size >= bigBoy) { bigBoy = size; lPeak = li; }
}

void interpCfrac(int jProj, int jT, int lx, int lDelta,
            int nLValues, int* lIs, int lStep, int lxMax,
            double* sR, double* sI, int sDimension,
            int lStart, int lMax, float* sOut4,
            int lBase2,
            double& bigBoy, int& lPeak, double* complexL, double* complexS)  // definition
{
    // Uses continued fraction to interpolate S-matrices

    double re, him, f, g, size;
    complex16 cfLi, cTemp;

    bigBoy = 0;
    lPeak = lStart;

    if (lStep != 1) {
        // Fit continued fraction to S-matrices
        int cfPointCount = 0;
        for (int liIndex = 1; liIndex <= nLValues; liIndex++) {
            int li = lIs[liIndex];  // 1-based
            if (li < lxMax) continue;
            cfPointCount++;
            // complexL and complexS are REAL*8(2,n_L_values) treated as complex
            complexL[2*(cfPointCount-1)]   = (double)li;
            complexL[2*(cfPointCount-1)+1] = 0.0;
            int ii = liIndex;
            complexS[2*(cfPointCount-1)]   = sR[(ii-1)*sDimension];
            complexS[2*(cfPointCount-1)+1] = sI[(ii-1)*sDimension];
        }

        // cfracInit: complex continued fraction setup
        cfracInit(cfPointCount, reinterpret_cast<complex16*>(complexL),
                        reinterpret_cast<complex16*>(complexS));

        // Compute S for lStart to lMax using continued fraction
        int liIndex = 1;
        for (int li = lStart; li <= lMax; li++) {
            cfLi = complex16((double)li, 0.0);
            cfracEval(reinterpret_cast<complex16*>(complexL),
                            reinterpret_cast<complex16*>(complexS),
                            cfLi, cTemp);
            re = std::real(cTemp);
            him = std::imag(cTemp);
            size = re*re + him*him;

            // Check against exact values when available
            // Advance liIndex past any L-values below current li (capped at nLValues)
            while (lIs[liIndex] - li < 0 && liIndex != nLValues) liIndex++;
            // After loop:
            //   lIs[liIndex] - li == 0  → true match
            //   lIs[liIndex] - li <  0  → liIndex hit nLValues, force match (original L510→L550)
            //   lIs[liIndex] - li >  0  → no match, skip update (original goto L780)
            if (lIs[liIndex] - li <= 0) {
                int ii = liIndex;
                f = sR[(ii-1)*sDimension];
                g = sI[(ii-1)*sDimension];
                size = f*f + g*g;
                liIndex++;
                if (li >= lxMax &&
                    (f-re)*(f-re) + (g-him)*(g-him) > 1.0e-6 * size)
                    std::printf(" **** WARNING: FOR CHAN, 2*JP, 2*JT, LX, LO-LI, LI =%5d%5d%5d%5d%5d%5d\n"
                                "     CONT. FRAC. =%13.5G%13.5G     EXACT =%13.5G%13.5G\n",
                                2, jProj, jT, lx, lDelta, li, re, him, f, g);
                re = f;
                him = g;
            }
            sOut4[0 + 2*sDimension*(li-lBase2-1)] = (float)re;
            sOut4[1 + 2*sDimension*(li-lBase2-1)] = (float)him;
            recordPeak(size, bigBoy, lPeak, li);
        }
        return;
    }

    // No interpolation (lStep = 1), just transfer and scan
    {
        int lBase1 = lIs[1] - 1;
        for (int li = lStart; li <= lMax; li++) {
            f = sR[(li-lBase1-1)*sDimension];
            g = sI[(li-lBase1-1)*sDimension];
            sOut4[0 + 2*sDimension*(li-lBase2-1)] = (float)f;
            sOut4[1 + 2*sDimension*(li-lBase2-1)] = (float)g;
            size = f*f + g*g;
            recordPeak(size, bigBoy, lPeak, li);
        }
    }

    return;
}

