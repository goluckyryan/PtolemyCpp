// CoulombWaveFunction_coulombIntegral.cpp — COULIN (Coulomb integrals by
// recursion), COULNG (penetrability factor), GETSCT (scattering Coulomb setup).

#include "CoulombWaveFunction.h"
#include "ptolemy_types.h"
#include "Timing.h"
#include "math/gauss_quadrature.h"
#include "math/special.h"
#include "Reaction.h"
#include "Constants.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>

// ============================================================================
// SECTION 3: CoulombWaveFunction::coulombIntegral — COULIN
//            CoulombWaveFunction::penetrability    — COULNG
//            CoulombWaveFunction::computeScatteringWaves — GETSCT
// ============================================================================


// ============================================================================
//
// COULOMB INTEGRALS BY RECURSION
//
//   FF = INTEGRAL( R TO INF ) DR F(lOut,ETAOUT,KOUT*R)
//        * F(lIn,ETAIN,KIN*R) / R**N
// AND SIMILARLY FOR FG, GF, AND GG  FOR
//   lMin =< (lOut+lIn)/2 =< lMax
//   -maxDel =< lIn-lOut =< maxDel
//
// 12/10/76 - FIRST VERSION - S. PIEPER
// ============================================================================
void CoulombWaveFunction::coulombIntegral(int rPower, int maxDel, int lMin, int lMax, double etaOut,
            double akOut, double* sigOut,
            double etaIn, double akIn, double* sigIn, double R, int includeIrregularG,
            double* ff, double* fg, double* gf, double* gg, int lDlDimension,
            double accuracy, int nPts,
            double* work, double* fIn, double* fOut, double* gIn, double* gOut,
            double* starts,
            int printLevel, int& returnCode, double& clTime,
            Reaction& reaction)
{
    constexpr int nTerms = 80;

    // Local variables
    double rTurns[3]; // 1-based
    double rMins[3];  // 1-based
    double accuracies[3];   // 1-based
    double a, b, c, d, e, x;
    // rTop/rChunk/pointIndex init to 0 silences -Wmaybe-uninitialized for the
    // R!=0 branch's returnCode-error printf (~line 364) — those vars are only
    // set in the R==0 branch (rChunk at line 215, rTop/pointIndex in the for-ii
    // loop at line 221+). The values are diagnostic-only; in the rare
    // R!=0 + returnCode!=0 failure they'll now print 0 instead of stack garbage.
    double rPowerDouble, acc, rBottom;
    // lInDouble/lOutDouble/lInLast/lOutLast/id2 init to 0 silences -Wmaybe-uninitialized for the
    // R!=0 trailing block (lines 421-438) that reads them after the lDelta inner
    // loop. In practice, lDelta=0 satisfies lSum=lIn+lOut on every first lSum
    // iteration (with all IABS guards), so the inner loop body always sets these
    // before the trailing block executes; the init is a no-op on the live path.
    double lInDouble = 0, lOutDouble = 0;
    double rValue = 0;
    double rChunk = 0, rTop = 0;
    double tt;
    int mxDel, lMaxExtended, lMinExtended, mxDelNeg;
    int termCount, stopBase, stopAt, lSumMin, lSumMax;
    int lSumIndex, lIndex, dCount, startSetIndex, dIndex;
    int id2 = 0;
    int lIn, lOut, lOutP, lSum, lDelta;
    int ii, i;
    int lInLast = 0, lOutLast = 0;
    int pointIndex = 0;
    int lMn, lMx;

    // STARTS is dimensioned STARTS(lDlDimension, 2, 4, 2) in Fortran
    // Access macro: STARTS(id, is, k, ii) => STARTS[ (ii-1)*lDlDimension*2*4 + (k-1)*lDlDimension*2 + (is-1)*lDlDimension + (id-1) ]
    // But we use a flat pointer, so define a macro
    #define starts4D(id,is,k,ii) starts[((ii)-1)*lDlDimension*2*4 + ((k)-1)*lDlDimension*2 + ((is)-1)*lDlDimension + ((id)-1)]

    // FF, FG, GF, GG are dimensioned (lDlDimension, *) — column-major
    #define ff2D(id,il) ff[((il)-1)*lDlDimension + ((id)-1)]
    #define fg2D(id,il) fg[((il)-1)*lDlDimension + ((id)-1)]
    #define gf2D(id,il) gf[((il)-1)*lDlDimension + ((id)-1)]
    #define gg2D(id,il) gg[((il)-1)*lDlDimension + ((id)-1)]

    // RCWFN-from-COULIN failure diagnostic, printed byte-identically before
    // each of the three `return;` bailouts on a nonzero computeFG returnCode.
    auto printCoulinError = [&] {
        std::printf("\n***RCWFN FROM COULIN IRET = %6d%15.5G%5d%5d%15.5G%15.5G%15.5G%15.5G\n $$$%20.10G%20.10G%20.10G%5d%5d\n",
            returnCode, rValue, lMn, lMx, etaIn, etaOut, akIn, akOut,
            rTop, rChunk, b, ii, pointIndex);
    };
    auto bailIfCoulinError = [&]() -> bool {
        if (returnCode != 0 && returnCode != 2) {
            printCoulinError();
            return true;
        }
        return false;
    };

    mxDel = std::max(2, maxDel);
    lMaxExtended = lMax + mxDel;
    lMinExtended = std::max(0, lMin - mxDel);
    if (lDlDimension <= mxDel) {
        // lDlDimension TOO SMALL
        returnCode = -10;
        return;
    }
    mxDelNeg = -mxDel;
    rPowerDouble = rPower;
    termCount = std::max(nTerms, lMax / 4);

    // FOLLOWING IS AN ATTEMPT TO ALLOW FOR LOSS OF ACCURACY

    acc = accuracy * std::exp(10 - 5 * rPower - 2 * std::fabs(etaIn - etaOut)
          - (lMax - lMin) / 100.);
    acc = std::max({acc, 1.0e-12, accuracy * .01});

    // IF THE INTEGRAL IS TO BE FROM 0, FIND A SUITABLE STARTING POINT

    rMins[1] = R;
    rMins[2] = R;
    accuracies[1] = acc;
    accuracies[2] = .01 * accuracy;
    if (R <= 0) {
        accuracies[1] = .1 * accuracy;
        accuracies[2] = .1 * acc;
        if (includeIrregularG) {
            returnCode = -12;
            return;
        }
        rTurns[1] = std::max((etaIn + std::sqrt(etaIn * etaIn + lMinExtended * (lMinExtended + 1))) / akIn,
            (etaOut + std::sqrt(etaOut * etaOut + lMinExtended * (lMinExtended + 1))) / akOut);
        rTurns[2] = std::max((etaIn + std::sqrt(etaIn * etaIn + lMaxExtended * (lMaxExtended + 1))) / akIn,
            (etaOut + std::sqrt(etaOut * etaOut + lMaxExtended * (lMaxExtended + 1))) / akOut);
        if (printLevel >= 1) std::printf(" RTURN'S = %12.4G%12.4G\n", rTurns[1], rTurns[2]);

        // WE DO THE NUMERIC INTEGRAL MORE EFFICENTLY HERE SINCE WE DO
        // ALL OF THEM AT ONCE.  SO FORCE CLINTS RIGHT TO THE PLACE WHERE
        // IT STARTS TO USE RCASYM.

        rMins[1] = 1.4 * rTurns[1];
        rMins[2] = 1.4 * rTurns[2];
    }

    // GET STARTING VALUES FOR RECURSIONS AND FINAL VALUES
    // TO CHECK THE STABILITY OF THE RECURSIONS.

    lSumIndex = 2 * lMin;
    if (printLevel >= 1) std::printf(" REQUESTED CLINTS ACCURACIES ARE %13.3G%13.3G\n", accuracies[1], accuracies[2]);

    // FOR RECURSION UPWARDS WE WILL NEED TO START WITH A
    // FULL WIDTH SET.  ALSO THE DOWNWARDS RECURSION CAN NOT GO
    // BENEATH  lIn+lOut = rPower-2  EVEN THOUGH FOR R > 0 THE RESULTS
    // WILL STILL BE FINITE.  THEREFORE WE GENERATE ALL THESE SMALL
    // L VALUES NUMERICALLY.

    stopBase = std::max(rPower - 2 * lMin, 2);
    if (includeIrregularG) stopBase = std::max(mxDel + 1 - 2 * lMin, stopBase);
    stopAt = stopBase;
    lSumMin = 2 * lMin + 1 + stopBase - 2;
    lSumMax = 2 * lMax - 1;

    tt = (float)dtime_();
    for (ii = 1; ii <= 2; ii++) {

        // lIn+lOut = lSumIndex+i-1

        for (i = 1; i <= stopAt; i++) {
            lIndex = lSumIndex - 2 * lMin + i;
            dCount = mxDel + ((mxDel + lIndex) % (2));
            startSetIndex = (i == stopAt) ? 2 : 1;
            for (dIndex = 1; dIndex <= dCount; dIndex++) {
                lIn = (lSumIndex + i + mxDel) / 2 + dIndex - 1 - mxDel;
                lOut = lSumIndex + i - 1 - lIn;
                if (lIn >= 0 && lOut >= 0) {
                    reaction.boundState.coulombIntegrals(rMins[ii], etaIn, etaOut, akIn, akOut,
                           sigIn[lIn], sigOut[lOut],
                           accuracies[ii],
                           starts4D(dIndex,startSetIndex,1,ii), starts4D(dIndex,startSetIndex,2,ii),
                           starts4D(dIndex,startSetIndex,3,ii), starts4D(dIndex,startSetIndex,4,ii),
                           work, work+nPts, fIn, fOut, gIn, gOut, work+2*nPts,
                           rPower, lIn, lOut, termCount, nPts, returnCode, printLevel);
                    if (returnCode != 0) return;
                }
            }

            // SAVE END OF THE FF RECURSION AND
            // STORE START OF FG, GF, GG STUFF

            for (dIndex = 1; dIndex <= dCount; dIndex++) {
                ff2D(dIndex, lIndex) = starts4D(dIndex, startSetIndex, 1, ii);
                if (includeIrregularG) {
                    fg2D(dIndex, lIndex) = starts4D(dIndex, startSetIndex, 2, ii);
                    gf2D(dIndex, lIndex) = starts4D(dIndex, startSetIndex, 3, ii);
                    gg2D(dIndex, lIndex) = starts4D(dIndex, startSetIndex, 4, ii);
                }
            }

        } // 69
        lSumIndex = 2 * lMax - 1;
        stopAt = 2;
    } // 79

    clTime = (float)dtime_() - tt;

    // Dump FF starting values before recursion
    for (int il=1; il<=4; il++) {
        for (int id=1; id<=lDlDimension; id++)
            std::fprintf(stderr, " FF(%d,%d)=%.8e", id, il, ff2D(id,il));
        std::fprintf(stderr, "\n");
    }
    if (R == 0) {

    // HERE WE DO THE INTEGRAL  0 < R < 1.4*RTURN  BACKWARDS UNTIL
    // WE COME TO INSIGNIFICANT CONTRIBUTIONS.

    rChunk = 3.1 / akIn;
    stopAt = stopBase;
    lMn = lMinExtended;
    lMx = lMin + (stopAt + mxDel + 1) / 2;
    lSumIndex = 2 * lMin;
    for (ii = 1; ii <= 2; ii++) {
        rTop = rMins[ii];
        while (true) {
        rBottom = std::max(0.0e0, rTop - rChunk);
        if (!(rBottom < rTop)) break;
        bool accuracyReached = false;
        b = rTop - rBottom;
        for (pointIndex = 1; pointIndex <= nPts; pointIndex++) {
            rValue = rTop - .5 * b * (1 + work[pointIndex]);
            CoulombWaveFunction::computeFG(akIn * rValue, etaIn, lMn, lMx, fIn, fOut,
                gIn, gOut, 1.0e-8, returnCode);
            if (bailIfCoulinError()) return;
            if (std::fabs(fIn[lMn]) < accuracy &&
                rValue < rTurns[ii]) { accuracyReached = true; break; }
            CoulombWaveFunction::computeFG(akOut * rValue, etaOut, lMn, lMx, fOut, gIn,
                gOut, &work[2 * nPts + 1], 1.0e-8, returnCode);
            if (bailIfCoulinError()) return;
            if (std::fabs(fOut[lMn]) < accuracy &&
                rValue < rTurns[ii]) { accuracyReached = true; break; }
            x = .5 * b * work[nPts + pointIndex] / std::pow(rValue, rPower);
            for (i = 1; i <= stopAt; i++) {
                lIndex = lSumIndex - 2 * lMin + i;
                dCount = mxDel + ((mxDel + lIndex) % (2));
                startSetIndex = (i == stopAt) ? 2 : 1;
                for (dIndex = 1; dIndex <= dCount; dIndex++) {
                    lIn = (lSumIndex + i + mxDel) / 2 + dIndex - 1 - mxDel;
                    lOut = lSumIndex + i - 1 - lIn;
                    if (lIn >= 0 && lOut >= 0) {
                        a = x * fIn[lIn] * fOut[lOut];
                        ff2D(dIndex, lIndex) = ff2D(dIndex, lIndex) + a;
                        starts4D(dIndex, startSetIndex, 1, ii) = ff2D(dIndex, lIndex);
                    }
                }
            }
        } // 159
        if (accuracyReached) break;
        rTop = rBottom;
        rValue = rBottom;
        }  // end while (rTop > 0, accuracy not reached)
        rMins[ii] = rValue;
        lSumIndex = 2 * lMax - 1;
        stopAt = 2;
        lMn = (lSumIndex - mxDel - 1) / 2;
        lMx = lMaxExtended;
    } // 169

    clTime = (float)dtime_() - tt;

    if (printLevel >= 1) std::printf(" RMIN'S =  %12.4G%12.4G\n", rMins[1], rMins[2]);
    dCount = mxDel + 1;
    if (printLevel >= 2) {
        std::printf("0FF STARTING VALUES:\n");
        for (startSetIndex = 1; startSetIndex <= 2; startSetIndex++) {
            for (dIndex = 1; dIndex <= dCount; dIndex++) {
                std::printf(" %20.10G", starts4D(dIndex, startSetIndex, 1, 2));
                if (dIndex % 3 == 0 || dIndex == dCount) std::printf("\n");
            }
        }
    }

    // RECURSE DOWNWARD FOR FF
    // WE GO FROM  lSum  TO  lSum-1  AND GENERATE ALL THE
    // lIn-lOut FOR EACH lSum.

    for (lSumIndex = lSumMin; lSumIndex <= lSumMax; lSumIndex++) {
        lSum = lSumMax + lSumMin - lSumIndex;
        for (lDelta = mxDelNeg; lDelta <= mxDel; lDelta++) {

            // lIn, lOut POINT TO THE CORNER WE WILL EXTRAPOLATE FROM

            lIn = (lSum - lDelta);
            if (lIn < 0) continue;
            lIn = lIn / 2;
            lOut = (lSum + lDelta) / 2;
            if (lSum != lIn + lOut) continue;
            lOutP = lOut - 1;
            if (std::abs(lIn - lOutP) > mxDel) continue;
            if (lOutP < 0) continue;
            if (std::abs(lOut + 1 - lIn) > mxDel) continue;

            // lIn, lOut-1  IS NEEDED AND WE HAVE ALL THE STUFF NEEDED TO
            // FIND IT. BY RECURSION DOWNWARD ON lOut.

            lInDouble = lIn;
            lOutDouble = lOut;
            dIndex = (lIn - lOut + mxDel) + 2;
            lIndex = lIn + lOut - 2 * lMin + 1;
            id2 = (dIndex + 1) / 2;
            dIndex = dIndex / 2;

            e = (2 * lOutDouble + 1) / ((lInDouble + lOutDouble - rPowerDouble + 2) * std::sqrt(lOutDouble * lOutDouble + etaOut * etaOut));
            a = etaOut * (lInDouble - rPowerDouble + 2) / (lOutDouble + 1)
                - (etaIn * (akIn / akOut)) * lOutDouble / (lInDouble + 1);
            b = (akIn / akOut) * (lOutDouble / (lInDouble + 1)) * std::sqrt((lInDouble + 1) * (lInDouble + 1) + etaIn * etaIn);
            c = lOutDouble * (lOutDouble - lInDouble + rPowerDouble - 1) * std::sqrt((lOutDouble + 1) * (lOutDouble + 1) + etaOut * etaOut)
                / ((lOutDouble + 1) * (2 * lOutDouble + 1));
            ff2D(id2, lIndex - 1) = e * (a * ff2D(dIndex, lIndex) + c * ff2D(id2 - 1, lIndex + 1)
                + b * ff2D(id2, lIndex + 1));
        }

        // lInDouble, lOutDouble, ETC, REFER TO THE LAST CASE FOR WHICH WE REDUCED
        // lOut.  NOW REDUCE lIn FOR THAT SAME CASE TO COMPLETE THE
        // GENERATION OF lSum-1 FROM lSum.

        if (lInDouble > 0) {
            e = (2 * lInDouble + 1) / ((lInDouble + lOutDouble - rPowerDouble + 2) * std::sqrt(lInDouble * lInDouble + etaIn * etaIn));
            a = etaIn * (lOutDouble - rPowerDouble + 2) / (lInDouble + 1)
                - (etaOut * (akOut / akIn)) * lInDouble / (lOutDouble + 1);
            b = (akOut / akIn) * (lInDouble / (lOutDouble + 1)) * std::sqrt((lOutDouble + 1) * (lOutDouble + 1) + etaOut * etaOut);
            c = lInDouble * (lInDouble - lOutDouble + rPowerDouble - 1) * std::sqrt((lInDouble + 1) * (lInDouble + 1) + etaIn * etaIn)
                / ((lInDouble + 1) * (2 * lInDouble + 1));
            ff2D(id2 - 1, lIndex - 1) = e * (a * ff2D(dIndex, lIndex) + b * ff2D(id2 - 1, lIndex + 1)
                + c * ff2D(id2, lIndex + 1));
        }
    } // 399

    } else {

    // GET THE INHOMO TERMS

    tt = (float)dtime_();

    CoulombWaveFunction::computeFG(akIn * R, etaIn, lMinExtended, lMaxExtended, fIn, work,
        gIn, &work[lMaxExtended], 1.0e-14, returnCode);
    if (returnCode == 0) {
        CoulombWaveFunction::computeFG(akOut * R, etaOut, lMinExtended, lMaxExtended, fOut, work,
            gOut, &work[lMaxExtended], 1.0e-14, returnCode);
    }
    if (returnCode != 0) {
        rValue = R;
        lMn = lMinExtended;
        lMx = lMaxExtended;
        b = 123456789.;
        printCoulinError();
        return;
    }

    clTime = clTime + (float)dtime_() - tt;

    // RECURSE UPWARDS ON FG, GF, AND GG. AND FF WHEN R > 0.

    for (lSum = lSumMin; lSum <= lSumMax; lSum++) {
        for (lDelta = mxDelNeg; lDelta <= mxDel; lDelta++) {

            // lIn, lOut POINT TO THE CORNER WE WILL EXTRAPOLATE FROM

            lIn = (lSum + lDelta) / 2;
            lOut = (lSum - lDelta) / 2;
            if (lSum != lIn + lOut) continue;
            lOutP = lOut + 1;
            if (std::abs(lIn - lOutP) > mxDel) continue;
            if (std::abs(lOut - 1 - lIn) > mxDel) continue;

            // lIn, lOut+1  IS NEEDED AND WE HAVE ALL THE STUFF NEEDED TO
            // FIND IT. BY RECURSION UPWARD ON lOut.

            lInDouble = lIn;
            lOutDouble = lOut;
            lInLast = lIn;
            lOutLast = lOut;
            dIndex = (lIn - lOut + mxDel) + 2;
            lIndex = lIn + lOut - 2 * lMin + 1;
            id2 = (dIndex + 1) / 2;
            dIndex = dIndex / 2;
            e = (2 * lOutDouble + 1) / ((lInDouble + lOutDouble + rPowerDouble) * std::sqrt((lOutDouble + 1) * (lOutDouble + 1) + etaOut * etaOut));
            a = etaOut * (lInDouble + rPowerDouble - 1) / lOutDouble
                - (etaIn * (akIn / akOut)) * (lOutDouble + 1) / lInDouble;
            b = (akIn / akOut) * ((lOutDouble + 1) / lInDouble) * std::sqrt(lInDouble * lInDouble + etaIn * etaIn);
            c = (lOutDouble + 1) * (lOutDouble - lInDouble - rPowerDouble + 1) * std::sqrt(lOutDouble * lOutDouble + etaOut * etaOut)
                / (lOutDouble * (2 * lOutDouble + 1));
            d = (lOutDouble + 1) / (akOut * std::pow(R, rPower));

            ff2D(id2 - 1, lIndex + 1) = e * (a * ff2D(dIndex, lIndex) + c * ff2D(id2, lIndex - 1) +
                b * ff2D(id2 - 1, lIndex - 1) + d * fOut[lOut] * fIn[lIn]);
            if (includeIrregularG) {
                fg2D(id2 - 1, lIndex + 1) = e * (a * fg2D(dIndex, lIndex) + c * fg2D(id2, lIndex - 1) +
                    b * fg2D(id2 - 1, lIndex - 1) + d * fOut[lOut] * gIn[lIn]);
                gf2D(id2 - 1, lIndex + 1) = e * (a * gf2D(dIndex, lIndex) + c * gf2D(id2, lIndex - 1) +
                    b * gf2D(id2 - 1, lIndex - 1) + d * gOut[lOut] * fIn[lIn]);
                gg2D(id2 - 1, lIndex + 1) = e * (a * gg2D(dIndex, lIndex) + c * gg2D(id2, lIndex - 1) +
                    b * gg2D(id2 - 1, lIndex - 1) + d * gOut[lOut] * gIn[lIn]);
            }
        }

        // lInDouble, lOutDouble, ETC, REFER TO THE LAST CASE FOR WHICH WE INCREASED
        // lOut.  NOW INCREASE lIn FOR THAT SAME CASE TO COMPLETE THE
        // GENERATION OF lSum+1 FROM lSum.

        e = (2 * lInDouble + 1) / ((lInDouble + lOutDouble + rPowerDouble) * std::sqrt((lInDouble + 1) * (lInDouble + 1) + etaIn * etaIn));
        a = etaIn * (lOutDouble + rPowerDouble - 1) / lInDouble
            - (etaOut * (akOut / akIn)) * (lInDouble + 1) / lOutDouble;
        b = (akOut / akIn) * ((lInDouble + 1) / lOutDouble) * std::sqrt(lOutDouble * lOutDouble + etaOut * etaOut);
        c = (lInDouble + 1) * (lInDouble - lOutDouble - rPowerDouble + 1) * std::sqrt(lInDouble * lInDouble + etaIn * etaIn)
            / (lInDouble * (2 * lInDouble + 1));
        d = (lInDouble + 1) / (akIn * std::pow(R, rPower));

        ff2D(id2, lIndex + 1) = e * (a * ff2D(dIndex, lIndex) + b * ff2D(id2, lIndex - 1)
            + c * ff2D(id2 - 1, lIndex - 1) + d * fOut[lOutLast] * fIn[lInLast]);
        if (includeIrregularG) {
            fg2D(id2, lIndex + 1) = e * (a * fg2D(dIndex, lIndex) + b * fg2D(id2, lIndex - 1)
                + c * fg2D(id2 - 1, lIndex - 1) + d * fOut[lOutLast] * gIn[lInLast]);
            gf2D(id2, lIndex + 1) = e * (a * gf2D(dIndex, lIndex) + b * gf2D(id2, lIndex - 1)
                + c * gf2D(id2 - 1, lIndex - 1) + d * gOut[lOutLast] * fIn[lInLast]);
            gg2D(id2, lIndex + 1) = e * (a * gg2D(dIndex, lIndex) + b * gg2D(id2, lIndex - 1)
                + c * gg2D(id2 - 1, lIndex - 1) + d * gOut[lOutLast] * gIn[lInLast]);
        }
    } // 599

    }

    // CHECK THE ACCURACY OF THE RECURSIONS

    for (i = 1; i <= 2; i++) {
        lIndex = i + stopBase - 2;
        dCount = mxDel + ((mxDel + lIndex) % (2));
        if (R > 0) {
            dCount = mxDel + ((mxDel + i + 1) % (2));
            for (dIndex = 1; dIndex <= dCount; dIndex++) {
                lIndex = 2 * (lMax - lMin) - 1 + i;
                lIn = (lIndex + 2 * lMin + mxDel) / 2 + dIndex - 1 - mxDel;
                lOut = lIndex + 2 * lMin - 1 - lIn;
                if (printLevel >= 1)
                    std::printf(" CHECKING %7d%7d%7d%7d%7d%7d\n          %20.10G%20.10G%20.10G%20.10G\n          %20.10G%20.10G%20.10G%20.10G\n",
                        i, dCount, dIndex, lIndex, lIn, lOut,
                        ff2D(dIndex, lIndex), fg2D(dIndex, lIndex), gf2D(dIndex, lIndex), gg2D(dIndex, lIndex),
                        starts4D(dIndex, i, 1, 2), starts4D(dIndex, i, 2, 2),
                        starts4D(dIndex, i, 3, 2), starts4D(dIndex, i, 4, 2));
                if (std::fabs(ff2D(dIndex, lIndex) / starts4D(dIndex, i, 1, 2) - 1) > 10 * accuracy)
                    std::printf("\n**** FOR LIN, LOUT =%4d%4d FF RECURSION IS POOR: %20.10G%20.10G\n",
                        lIn, lOut, ff2D(dIndex, lIndex), starts4D(dIndex, i, 1, 2));
                if (!includeIrregularG) continue;
                if (std::fabs(fg2D(dIndex, lIndex) / starts4D(dIndex, i, 2, 2) - 1) > 10 * accuracy)
                    std::printf("\n**** FOR LIN, LOUT =%4d%4d FG RECURSION IS POOR: %20.10G%20.10G\n",
                        lIn, lOut, fg2D(dIndex, lIndex), starts4D(dIndex, i, 2, 2));
                if (std::fabs(gf2D(dIndex, lIndex) / starts4D(dIndex, i, 3, 2) - 1) > 10 * accuracy)
                    std::printf("\n**** FOR LIN, LOUT =%4d%4d GF RECURSION IS POOR: %20.10G%20.10G\n",
                        lIn, lOut, gf2D(dIndex, lIndex), starts4D(dIndex, i, 3, 2));
                if (std::fabs(gg2D(dIndex, lIndex) / starts4D(dIndex, i, 4, 2) - 1) > 10 * accuracy)
                    std::printf("\n**** FOR LIN, LOUT =%4d%4d GG RECURSION IS POOR: %20.10G%20.10G\n",
                        lIn, lOut, gg2D(dIndex, lIndex), starts4D(dIndex, i, 4, 2));
            }
        } else {
            for (dIndex = 1; dIndex <= dCount; dIndex++) {
                lIn = (2 * lMin + lIndex + mxDel) / 2 + dIndex - 1 - mxDel;
                lOut = 2 * lMin + lIndex - 1 - lIn;
                if (lIn < 0 || lOut < 0) continue;
                if (printLevel >= 1)
                    std::printf(" CHECKING %7d%7d%7d%7d%7d%7d\n          %20.10G%20.10G\n",
                        i, dCount, dIndex, lIndex, lIn, lOut,
                        ff2D(dIndex, lIndex), starts4D(dIndex, i, 1, 1));
                if (std::fabs(ff2D(dIndex, lIndex) / starts4D(dIndex, i, 1, 1) - 1) > 10 * accuracy)
                    std::printf("\n**** FOR LIN, LOUT =%4d%4d FF RECURSION IS POOR: %20.10G%20.10G\n",
                        lIn, lOut, ff2D(dIndex, lIndex), starts4D(dIndex, i, 1, 1));
            }
        }
    } // 869

    returnCode = 0;
    return;

    #undef starts4D
    #undef ff2D
    #undef fg2D
    #undef gf2D
    #undef gg2D
}


// ============================================================================
// FUNCTION COULNG — lines 11133-11208
//
// COMPUTE COULOMB WAVEFUNCTION FOR NEGATIVE ENERGY
//
// THE RETURNED FUNCTION IS
//           * EXP(RHO) * (2*RHO)**ETA * W( -ETA, L+1/2, 2*RHO )
// WHERE W IS THE WHITTAKER FUNCTION OF ABRAMOWITZ AND STEGUN
//
// 4/23/75 - FIRST VERSION - S. PIEPER
// ============================================================================
double CoulombWaveFunction::penetrability(int L, double eta, double rho, double aNorm)
{
    // Local variables
    int pointCount = 28;
    double laguerrePoints[32];  // 1-based, extra space for LAGBC/LAGUER
    double laguerreWeights[32]; // 1-based
    // own.
    double result;
    double eta1, l1Double, gamma, xCutoff, temp;
    // CSX/CSW/TSX/TSW dropped — were LAGUER out-params, never read here.
    int l1, l2, etaInteger, nodeIndex;

    eta1 = std::fmod(eta, 1.0e0);
    l1 = (eta < 0) ? 0 : ((L) % (10));
    l2 = L - l1;
    l1Double = l1;
    eta1 = eta1 + std::max(std::min(10 - l1Double, eta - eta1), 0.0e0);
    etaInteger = (int)(eta - eta1);

    gamma = 1;
    if (eta + L != 0) gamma = std::exp(-dLogGamma(aNorm) / (eta + L));
    xCutoff = 0;
    if (l2 + etaInteger > 1) {
        xCutoff = std::exp(-Constants::BIGLOG / (l2 + etaInteger)) / gamma;
    }

    laguerre(pointCount, &laguerrePoints[1], &laguerreWeights[1], eta1 + l1Double);

    result = 0;
    for (nodeIndex = 1; nodeIndex <= pointCount; nodeIndex++) {
        if (laguerrePoints[nodeIndex] < xCutoff) continue;
        if (rho < 2) {
            temp = 1 + laguerrePoints[nodeIndex] / (2 * rho);
            result = result + laguerreWeights[nodeIndex] * std::pow(gamma * laguerrePoints[nodeIndex] / temp, l2 + etaInteger)
                * std::pow(temp, L + l2 - eta1);
        } else {
            result = result + laguerreWeights[nodeIndex] * std::pow(gamma * laguerrePoints[nodeIndex], l2 + etaInteger) *
                std::pow(1 + laguerrePoints[nodeIndex] / (2 * rho), L - eta);
        }
    }

    result = result * std::pow(gamma, eta1 + l1Double);
    if (eta + L == 0) result = result / dGamma(aNorm);
    return result;
}


// ============================================================================
//
// COMPUTES SCATTERING WAVES AT lMin AND lCrit FOR GRDSET
//
// 7/8/75 - FIRST VERSION - S. P.
// ============================================================================
void CoulombWaveFunction::computeScatteringWaves(int& returnCode, Reaction& reaction)
{

    auto& L       = reaction.angMom.L;
    auto& lMax    = reaction.angMom.lMax;
    auto& lMin    = reaction.angMom.lMin;



    auto& waveChannel  = reaction.internalState.waveChannel;

    auto& isStandalone  = reaction.distortedWave.scatteringSolver.isStandalone;

    auto& lCrit   = reaction.kin.lCrit;

    auto& rOfMax  = reaction.boundState.data.rOfMax;


    // Local scratch: WAVELJ may rewrite linkuleAddr while solving a single L for the
    // form factor — save it here and restore after the solve. Was scratchWork.LNKAD3.
    int savedLinkuleAddr[numLinkules + 1][7] = {};

    std::vector<double>* scrtsVectorPointers[2] = {
        &reaction.boundState.data.sctmnArr,
        &reaction.boundState.data.sctcrArr,
    };

    // Local variables
    float fgSeed[5]; // 1-based, REAL*4
    int savedSpinorbit;
    float dummy4[2]; // REAL*4 DUMMY4(1)
    int lCritUsed;
    int gridPointCount;
    int ii, i;
    // psi init to 0 silences -Wmaybe-uninitialized for the post-loop
    // verbosity>=4 "FINAL PSI" printf (~line 718) — set in the i=1..gridPointCount
    // loop body, read after the loop. Unreachable in practice (gridPointCount>0
    // always) but GCC can't see gridPointCount's lower bound.
    double rValue, psiMax, rPsiMax, rStep;
    double psi = 0;

    int verbosity = ((reaction.flags.printLevel) % (10));

    // INITIALIZE THE COULOMB FUNCTIONS, ETC.

    for (waveChannel = 1; waveChannel <= 2; waveChannel++) {
        bool ok = reaction.setupWavefunctionPotential();
        if (verbosity >= 5) std::printf(" ++++ setupWavefunctionPotential CALL, waveChannel, IRET = %10d%10d\n",
            waveChannel, ok ? 1 : 0);
        if (!ok) { returnCode = 0; return; }
    }

    lCritUsed = lCrit;
    if (lCritUsed < lMin || lCritUsed > lMax) lCritUsed = (lMin + lMax) / 2;

    if (verbosity >= 3) std::printf("\n L CRITICAL:  AVERAGE = %4d     USED FOR GRID SETUP = %4d\n",
        lCrit, lCritUsed);

    // COMPUTE THE WAVEFUNCTIONS, USING CENTRAL POTENTIALS ONLY.

    isStandalone = TRUE_F;
    savedSpinorbit = reaction.distortedWave.channel[1].hasSpinorbit;
    reaction.distortedWave.channel[1].hasSpinorbit = FALSE_F;

    L = std::max(0, lMin - reaction.inelastic.lxMax);
    gridPointCount = reaction.distortedWave.channel[1].nGridSteps + 1;
    rStep = reaction.distortedWave.channel[1].rStart;

    // WAVELJ WILL SET LNKBLK BACK TO THE STATUS FOR THE INCOMING
    // CHANNEL.  IT MAY HAVE BEEN CHANGED JUST FOR THE FORM factor.
    // THEREFORE WE SAVE ITS STATUS.

    for (ii = 1; ii <= numLinkules; ii++) {
        for (i = 1; i <= 6; i++) {
            savedLinkuleAddr[ii][i] = reaction.linkuleData.linkuleAddr[ii][i];
        }
    }

    for (ii = 1; ii <= 2; ii++) {
        // Dead Fortran-era locals (LWAVR/LWAVI/LVREAL/LVIMAG/LCENTR) dropped — never read.

        // ALWAYS FORCE WAVELJ TO START OVER AND STEP OUT SLOWLY LOOKING FOR
        // OVERFLOWS.  THIS IS NECESSARY SINCE WE MAKE A LARGE JUMP IN L THAT
        // WAVELJ IS NOT EXPECTING.

        reaction.distortedWave.channel[1].lastL = 999999;
        fgSeed[1] = (float)(reaction.distortedWave.channel[1].nF1sArr[L]);
        fgSeed[2] = (float)(reaction.distortedWave.channel[1].nG1sArr[L]);
        fgSeed[3] = (float)(reaction.distortedWave.channel[1].F_arr[L]);
        fgSeed[4] = (float)(reaction.distortedWave.channel[1].G_arr[L]);
        reaction.distortedWave.scatteringSolver.solvePartialWave(L, 2 * L + reaction.distortedWave.channel[1].twoSpin, 1, 0, dummy4, fgSeed,
            dummy4,
            reaction.distortedWave.scatteringSolver.wavRPointer, reaction.distortedWave.scatteringSolver.wavIPointer,   // use WAVCOM pointer fields
            reaction.distortedWave.channel[1].rlvPointer,
            reaction.distortedWave.channel[1].imvPointer,
            reaction.distortedWave.channel[1].centPointer,
            reaction);

        // JPTOLX HAS STORED AN S-MATRIX.  LATER ON WE WILL ACCUMULATE A
        // sum OF CONTRIBUTIONS IN THIS SPOT SO WE MUST ZERO IT AGAIN.

        reaction.distortedWave.channel[1].smatArr[2 * L]     = 0;
        reaction.distortedWave.channel[1].smatArr[2 * L + 1] = 0;

        // COMPUTE AND STORE   RPSI

        // iterations of the ii=1..2 loop now run unconditionally.
        psiMax = 0;
        rPsiMax = 0;
        // Size = gridPointCount (= nGridSteps + 1). 0-based access: scrts[i-1].
        std::vector<double>& scrtsVector = *scrtsVectorPointers[ii - 1];
        scrtsVector.assign(gridPointCount, 0.0);
        double* scrts = scrtsVector.data();  // 0-based: scrts[i-1]
        rValue = -rStep;

        for (i = 1; i <= gridPointCount; i++) {
            rValue = rValue + rStep;
            psi = std::fabs(reaction.distortedWave.scatteringSolver.wavRPointer[i]) +
                std::fabs(reaction.distortedWave.scatteringSolver.wavIPointer[i]);   // use WAVCOM pointer fields
            if (psi > psiMax) rOfMax = rValue;
            psiMax = std::max(psiMax, psi);
            rPsiMax = std::max(rPsiMax, rValue * psi);
            scrts[i-1] = rValue * psi;
        } // 389
        if (verbosity >= 4) std::printf(" FOR L = %4d, MAX(PSI) = %11.3G   IS AT R = %5.1f     MAX(R*PSI) = %11.3G     FINAL PSI = %11.3G\n",
            L, psiMax, rOfMax, rPsiMax, psi);

        L = lCritUsed;
    } // 399

    isStandalone = FALSE_F;
    reaction.distortedWave.channel[1].hasSpinorbit = savedSpinorbit;

    for (ii = 1; ii <= numLinkules; ii++) {
        for (i = 1; i <= 6; i++) {
            reaction.linkuleData.linkuleAddr[ii][i] = savedLinkuleAddr[ii][i];
        }
    }


    // SAVE STUFF FOR BSPROD

    reaction.boundState.data.scatPointCount = gridPointCount;
    reaction.boundState.data.scatRMax = reaction.distortedWave.channel[1].asymptopia;
    reaction.boundState.data.scatInvStep = 1 / rStep;
    return;
}
