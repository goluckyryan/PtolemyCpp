// angular_coupling.cpp — A12: angular-momentum transformation function.

#include "Reaction.h"
#include "angular_coupling.h"
#include "math/angular_momentum_coeff.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

// The single internal use `jA12S[jA12M + 2] = LXTMP[...] + LHSM1` is now
// just `jA12S[jA12M + 2] = LXTMP[...]`.
void angularCoupling12(int li, int lxMin, int lxMax, int lMin, int lMax, double* xLamArg,
         int* nLamArg, double* a12VlArg, double* msvalArg,
         int* jA12SArg, int& jA12M, int& jA12N, int& jA12An, int* iIndexFlat,
         int& hCount, int& loMinMin, int& loMaxMax, double* dIntsArg, double* outTempArg,
         double* xLoTempArg, int* lxTempArg, int la12Vl, int printLevel,
         int* indxDwFlat, int& dwCount, int& numIi, int* iDwfiFlat, int* iDwfoFlat,
         Reaction& reaction)
{
    // All arrays are passed with Fortran 1-based convention.
    // Adjust base pointers so that arr[N] = Fortran arr(N).
    double* xLam   = xLamArg - 1;
    int*    nLam   = nLamArg - 1;
    double* a12Vl  = a12VlArg - 1;
    double* msval = msvalArg - 1;
    int*    jA12S  = jA12SArg - 1;
    double* dInts  = dIntsArg - 1;
    double* outTemp = outTempArg - 1;
    double* xLoTemp = xLoTempArg - 1;
    int*    lxTemp  = lxTempArg - 1;

    // Access iIndex, indxDw, iDwfi, iDwfo as 2D arrays (column-major, 1-based)
    #define iIndex(i,j) iIndexFlat[((j)-1)*4 + ((i)-1)]
    #define indxDw(i,j) indxDwFlat[((j)-1)*4 + ((i)-1)]
    #define iDwfi(i,j)  iDwfiFlat[((j)-1)*3 + ((i)-1)]
    #define iDwfo(i,j)  iDwfoFlat[((j)-1)*4 + ((i)-1)]

    auto& lBoundProj    = reaction.boundState.vertex[1].lBound;
    auto& lBoundTarg    = reaction.boundState.vertex[2].lBound;
    auto& jBt    = reaction.boundState.vertex[2].jB;
    auto& mStop  = reaction.inelastic.mStop;

    // SAVE variables (static)
    static bool isHalfInterval;
    static int mtMin, mpMin, lamIndex, loMost;
    static int llBt, llBp;
    static double xN;
    static bool lDebugSwitch, mDebugSwitch;

    // Local
    // lo init to 0 silences -Wmaybe-uninitialized for the mDebugSwitch "START"
    // debug printf (~line 101) — it dumps lo alongside lBoundTarg/lBoundProj etc.
    // *before* the per-li loop ever assigns lo. With mDebugSwitch disabled
    // (the default in production) the read is unreachable.
    int lo = 0;
    int lx, mu, mX, muPos, loStart, lxStart;
    int hIndex, ix, ll, mm, numM, llEven, mLow;
    int liPhase, loPhase, lTargPhase, lProjPhase, phaseSum;
    int jA12On, muLast, muOne;
    int zeroCount;
    double outer, temp, threeJValue;
    double tempMin, tempMax, tempAverage;
    bool isFirstValid;

    // ========================================================================
    // Lambda-table setup. INIT_arg parameter dropped — sole caller passed
    // 5 always, so this block always ran; gate removed.
    // ========================================================================
    {
    {
        int i = (printLevel / 10000) % 10;
        lDebugSwitch = (i >= 4) || ((printLevel % 10) >= 4);
        mDebugSwitch = (i >= 5);
    }

    // If lBoundTarg or lBoundProj is odd, halve that interval
    isHalfInterval = false;
    mtMin = -lBoundTarg;
    mpMin = -lBoundProj;
    if (lBoundTarg % 2 != 0) {
        isHalfInterval = true;
        mtMin = 1;
    } else if (lBoundProj % 2 != 0) {
        isHalfInterval = true;
        mpMin = 1;
    }
    lamIndex = 1;
    loMost = lMax + std::max({lxMax, lBoundProj, lBoundTarg});

    nLam[1] = 1;  // 1-based
    xLam[1] = 1.0;
    {
        int M = 0;
        outer = 1.0;
        if (loMost != 0) {

        if (mDebugSwitch) {
            std::printf("\nSTART%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%5s\n",
                        lBoundTarg, lBoundProj, lo, lxMin, lxMax, lMin, lMax,
                        loMost, mtMin, mpMin, mStop, reaction.inelastic.a12mSize, reaction.inelastic.nMloLx,
                        isHalfInterval ? " T" : " F");
        }

        if (!isHalfInterval) {

        // Generate LAMBDA(L, M)
        for (ll = 1; ll <= loMost; ll++) {
            M = 1 - M;
            outer = outer * sqrt(dInts[ll + M - 1] / dInts[ll + M]);
            if (M == 1) outer = -outer;
            lamIndex++;
            nLam[ll + 1] = lamIndex;
            xLam[lamIndex] = outer;
            mu = std::min(mStop, ll);
            int mmtMin = M + 2;
            if (mmtMin > mu) continue;
            for (mm = mmtMin; mm <= mu; mm += 2) {
                lamIndex++;
                xLam[lamIndex] = -xLam[lamIndex - 1] *
                    sqrt(dInts[ll - mm + 2] * dInts[ll + mm - 1] /
                        (dInts[ll + mm] * dInts[ll - mm + 1]));
            }
        }
        } else {
        // Generate LAMBDA(L, L-M) for halved interval
        for (ll = 1; ll <= loMost; ll++) {
            outer = -sqrt(dInts[2 * ll - 1] / dInts[2 * ll]) * outer;
            llEven = ll - (ll % 2);
            mLow = std::min(mStop, llEven);
            numM = mLow / 2 + 1;
            lamIndex = lamIndex + numM;
            xLam[lamIndex] = outer;
            nLam[ll + 1] = lamIndex - ll / 2;
            if (numM == 1) continue;
            mm = ll;
            int lamBackIndex = lamIndex;
            for (int i = 2; i <= numM; i++) {
                mm = mm - 2;
                lamBackIndex--;
                xLam[lamBackIndex] = -xLam[lamBackIndex + 1] * sqrt(dInts[ll - mm - 1]
                    * dInts[ll + mm + 2] / (dInts[ll - mm] * dInts[ll + mm + 1]));
            }
        }
        } // end isHalfInterval if/else
        } // end loMost != 0
    } // end M scope

    // INIT is modified in Fortran via SAVE — here we use static tracking
    llBt = 2 * lBoundTarg;
    llBp = 2 * lBoundProj;

    // ========================================================================
    // Main computation
    // ========================================================================
    } // end setup block (always runs now)
    {
        xN = dInts[2 * li + 1] * dInts[llBt + 1] * dInts[llBp + 1];
        xN = 0.5 * sqrt(xN);

        // Set up lx table: for each lx store (IH_start, loMin, loMax)
        loMinMin = li + 500;
        hIndex = 0;
        ix = 0;
        for (lx = lxMin; lx <= lxMax; lx++) {
            int loMin = abs(li - lx);
            int loMax = li + lx;
            loMin = loMin + (lBoundTarg + lBoundProj + li + loMin) % 2;
            loMax = loMax - (lBoundTarg + lBoundProj + li + loMax) % 2;
            lxTemp[ix + 1] = hIndex + 1;   // 1-based
            lxTemp[ix + 2] = loMin;
            lxTemp[ix + 3] = loMax;
            ix += 3;
            loMinMin = std::min(loMin, loMinMin);
            if (loMin > loMax) continue;
            hIndex = hIndex + (loMax - loMin) / 2 + 1;
        }

        hCount = hIndex;
        loMaxMax = lx > lxMin ? lxTemp[ix] : 0;  // loMax of last lx
        if (lxMax >= lxMin) loMaxMax = lxTemp[3 * (lxMax - lxMin) + 3];

        // Build indxDw array: pairs of (incident, outgoing) waves
        dwCount = 0;
        int jpoMin = std::max(0, 2 * li - reaction.distortedWave.channel[1].twoSpin - jBt);  // JBT
        for (int kWI = 1; kWI <= reaction.gridData.nWfi; kWI++) {
            int lasI = iDwfi(1, kWI) + li;
            if (lasI < lMin || lasI > lMax) continue;
            int jPi = iDwfi(2, kWI) + 2 * li;
            if (jPi < abs(2 * li - reaction.distortedWave.channel[1].twoSpin)) continue;
            int jpoMax = jPi + jBt;  // JBT
            if (reaction.distortedWave.channel[1].hasSpinorbit) jpoMin = abs(jPi - jBt);
            for (int kWO = 1; kWO <= reaction.gridData.nWfo; kWO++) {
                lo = iDwfo(1, kWO) + li;
                if (lo < loMinMin) continue;
                if (!reaction.distortedWave.channel[2].hasSpinorbit) {
                    // No spin-orbit on channel 2
                    if (2 * lo + reaction.distortedWave.channel[2].twoSpin < jpoMin
                        || 2 * lo - reaction.distortedWave.channel[2].twoSpin > jpoMax) continue;
                } else {
                    if (iDwfo(2, kWO) + li < 0) continue;
                    int jPo = iDwfo(3, kWO) + 2 * li;
                    if (jPo < abs(2 * lo - reaction.distortedWave.channel[2].twoSpin)) continue;
                    if (jPo < jpoMin || jPo > jpoMax) continue;
                }
                // Store pointers
                dwCount++;
                indxDw(1, dwCount) = kWI;
                indxDw(2, dwCount) = kWO;
                indxDw(3, dwCount) = iDwfi(3, kWI);
                indxDw(4, dwCount) = iDwfo(4, kWO);
            }
        }

        // Build iIndex array
        numIi = 0;
        ix = 0;
        for (lx = lxMin; lx <= lxMax; lx++) {
            hIndex = lxTemp[ix + 1] - 1;
            int loMin = lxTemp[ix + 2];
            int loMax = lxTemp[ix + 3];
            ix += 3;
            for (int kDw = 1; kDw <= dwCount; kDw++) {
                int kWO = indxDw(2, kDw);
                lo = iDwfo(1, kWO) + li;
                if (lo < loMin || lo > loMax) continue;
                numIi++;
                // smivlArr is class-owned; store 0-based offset (reader adds +uIndex-1 via smivlArr.data()[ITEMP]).
                iIndex(1, numIi) = (hIndex + (lo - loMin) / 2) * reaction.gridData.nInterpPoints - 1;
                // 1-based via reaction.gridData.dwPointer.
                // Store the dw_ subscript for the imaginary half of pair kDw (1-based);
                // the reader in DWBAGrid::inelDc (inelastic_dwba.cpp:659, numIi-indexed
                // accumulation loop) does dwPointer[kDw-1] (real) / dwPointer[kDw] (imag),
                // where dwPointer = reaction.gridData.dwPointer.
                iIndex(2, numIi) = 2 * kDw;
                iIndex(3, numIi) = kDw;
                iIndex(4, numIi) = lx;
            }
        }

        zeroCount = 0;

        // Prepare for the triple sum on mp, mt, mu
        int muStart = isHalfInterval ? -li : li % 2;

        jA12M = -4;
        jA12N = 0;
        jA12An = 0;
        tempMin = 10.0;
        tempMax = 0;
        tempAverage = 0;

        if (loMinMin > loMaxMax) return;

        // Main triple loop: mp, mt, mu
        for (int mp = mpMin; mp <= lBoundProj; mp += 2) {
            for (int mt = mtMin; mt <= lBoundTarg; mt += 2) {
                mX = mt + mp;
                if (abs(mX) > lxMax) continue;
                int lamLtIndex = nLam[lBoundTarg + 1] + abs(mt) / 2;
                int lamLpIndex = nLam[lBoundProj + 1] + abs(mp) / 2;
                lTargPhase = 0;
                lProjPhase = 0;
                if (mt < 0) lTargPhase = lBoundTarg;
                if (mp < 0) lProjPhase = lBoundProj;

                // Setup A12 factors for each lx
                isFirstValid = false;
                for (lx = lxMin; lx <= lxMax; lx++) {
                    outer = xLam[lamLtIndex] * xLam[lamLpIndex] * xN *
                        threeJ(llBt, llBp, 2 * lx, 2 * mt, 2 * mp, -2 * mX);
                    outTemp[lx - lxMin + 1] = outer;
                    if (outer != 0) isFirstValid = true;
                }
                if (!isFirstValid) continue;  // skip this (mt, mp) pair

                isFirstValid = true;
                for (mu = muStart; mu <= li; mu += 2) {
                    muPos = abs(mX - mu);
                    if (muPos > loMaxMax) continue;

                    if (isFirstValid) {
                        // First valid mu: initialize xLoTemp
                        for (lo = loMinMin; lo <= loMaxMax; lo++) {
                            xLoTemp[lo - loMinMin + 1] = 0;
                        }
                    }

                    loStart = std::max(loMinMin, muPos);

                    // Build lo-dependent factors
                    for (lo = loMinMin; lo <= loMaxMax; lo += 2) {
                        outer = 0;
                        if (lo >= loStart) {
                            outer = xLoTemp[lo - loMinMin + 1];
                            if (outer != 0) {
                                // Iterate to new lambda's
                                outer = outer * sqrt(dInts[li - mu + 2]
                                    * dInts[li + mu - 1] * dInts[lo - mX + mu - 1]
                                    * dInts[lo + mX - mu + 2] / (dInts[li + mu]
                                    * dInts[li - mu + 1] * dInts[lo - mX + mu]
                                    * dInts[lo + mX - mu + 1]));
                            } else {
                                // First mu valid for this lo, li
                                liPhase = 0;
                                loPhase = 0;
                                if (mu < 0) liPhase = li;
                                if ((mX - mu) < 0) loPhase = loMinMin;
                                phaseSum = lProjPhase + lTargPhase + loPhase + liPhase;
                                outer = xLam[nLam[li + 1] + abs(mu) / 2] *
                                    xLam[nLam[lo + 1] + muPos / 2] *
                                    sqrt(dInts[2 * lo + 1]);
                                if (phaseSum % 2 != 0) outer = -outer;
                            }
                        }
                        xLoTemp[lo - loMinMin + 1] = outer;
                    }

                    // Compute A12 terms for all valid (lx, lo)
                    lxStart = std::max(abs(mX), lxMin);
                    jA12On = jA12N + 1;
                    for (lx = lxStart; lx <= lxMax; lx++) {
                        int loMin = lxTemp[3 * (lx - lxMin) + 2];
                        int loMax = lxTemp[3 * (lx - lxMin) + 3];
                        if (loMin > loMax) continue;
                        for (lo = loMin; lo <= loMax; lo += 2) {
                            threeJValue = threeJ(2 * li, 2 * lo, 2 * lx,
                                            2 * mu, 2 * (mX - mu), -2 * mX);
                            temp = xLoTemp[lo - loMinMin + 1] *
                                   outTemp[lx - lxMin + 1] * threeJValue;

                            double absTemp = fabs(temp);
                            if (temp == 0) zeroCount++;
                            if (absTemp > tempMax) tempMax = absTemp;
                            if (absTemp < tempMin && absTemp != 0)
                                tempMin = absTemp;
                            if (isHalfInterval || mu != 0) temp = temp + temp;
                            tempAverage += fabs(temp);
                            jA12N++;
                            a12Vl[jA12N] = temp;
                        }
                    }
                    muLast = mu;

                    if (!isFirstValid) continue;
                    isFirstValid = false;

                    // First mu of this (mp, mt) pair — store info
                    muOne = mu;
                    jA12M += 5;
                    msval[jA12M]     = dInts[mt];
                    msval[jA12M + 1] = dInts[mp];
                    msval[jA12M + 2] = dInts[mu] - 2;
                    jA12S[jA12M + 2]  = lxTemp[3 * (lxStart - lxMin) + 1];  // + LHSM1 (=0) dropped
                    jA12S[jA12M]      = la12Vl + jA12On - jA12S[jA12M + 2];
                    jA12S[jA12M + 3]  = hCount - lxTemp[3 * (lxStart - lxMin) + 1] + 1;
                } // end mu loop

                if (!isFirstValid) {
                    jA12S[jA12M + 1] = (muLast - muOne) / 2 + 1;
                    jA12An += jA12S[jA12M + 1];
                }
            } // end mt loop
        } // end mp loop

        // Statistics
        if (jA12N > 0) tempAverage = tempAverage / jA12N;
        if (lDebugSwitch) {
            std::printf("\nLI, LOMNMN, LOMXMX, hCount, jA12M, jA12An, jA12N:%4d%4d%4d%4d%4d%8d%8d\n"
                        " NUM OF ZEROS; SMALLEST, LARGEST AND AVERAGE A12VAL:%6d%15.5G%15.5G%15.5G\n",
                        li, loMinMin, loMaxMax, hCount, jA12M, jA12An, jA12N,
                        zeroCount, tempMin, tempMax, tempAverage);
        }
    }

    #undef iIndex
    #undef indxDw
    #undef iDwfi
    #undef iDwfo

    return;
}
