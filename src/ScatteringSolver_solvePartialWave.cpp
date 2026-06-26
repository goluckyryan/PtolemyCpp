// ScatteringSolver_solvePartialWave.cpp — WAVELJ: Numerov integration of the radial
// Schrodinger equation for one partial wave (L, jProj) with optional spin-orbit;
// computes the outgoing wavefunction and S-matrix element.

#include "ptolemy_types.h"
#include "Timing.h"
#include "linkule.h"
#include "Reaction.h"
#include "Constants.h"
#include <cstdio>
#include <cmath>


// ============================================================================
// SECTION 4: wavefunction_solver.cpp — ScatteringSolver::solvePartialWave, WAVELJ wrapper
// ============================================================================

void ScatteringSolver::solvePartialWave(int L, int jProj, int channelIndex, int nPts, float* rGrid, float* waveR,
            float* waveI, double* waveReal, double* waveImag, double* vReal,
            double* vImag, double* vCent,
            Reaction& reaction)
{
//
//     COMPUTES SCATTERING WAVE FOR A GIVEN CHANNEL, L, jProj.
//
//     THE INPUT IS
//
//     L - L VALUE.
//     jProj - J OF PROJECTILE - USED ONLY FOR SPIN ORBIT FORCE.
//     channelIndex - CHANNEL NUMBER (1 OR 2).
//     nPts - NUMBER OF DESIRED OUTPUT (BY INTERPOLATION) POINTS.
//     rGrid - VALUES OF R AT WHICH OUTPUT IS TO BE FOUND.
//     waveR, waveI - THE COMPUTED WAVE FUNCTION AT rGrid IS STORED HERE.
//     vReal, vImag - THE REAL AND IMAGINARY PARTS OF THE POTENTIAL IN NUMEROV FORM.
//     vCent - 1/RHO**2 IN THE NUMEROV FORM.
//

    // Local variables
    double lDouble, lLp1, h;
    double thisR, thisI, thisRSaved, aA1, waveThreshold;
    double d, a1, a2, cReal, cImag, sJr, sJi;
    double F, G, f1, g1;
    double sDotL;
    double rT, rInv, p, pSquared, x1, x2, x3, x4, x5;
    double c1, c2, c3, c4, c5;
    double rByH, stepInverse;
    double rStart;
    double tt1, t1, t2;
    int i, ii, rGridIndex, n, n1, n2, stepCount, step2Count, startIndex, returnCode;
    int channelNameIndex, srcMaxIndex;
    // extraSteps init to 0 silences -Wmaybe-uninitialized for the useFullRange
    // branch — when L >= lastL the assignment block at line ~254 is skipped
    // AND every later extraSteps read is gated by the same (!useFullRange)
    // condition (the loop body always hits `if (n2 == stepCount) break;` on the
    // first iteration). GCC can't see the guards are loop-invariant.
    int extraSteps = 0;
    logical isBadRenorm;

    const char* channelName[4]; // 1-based
    channelName[1] = "INCOMING";
    channelName[2] = "OUTGOING";
    channelName[3] = "        ";

    double twelth = .0833333333333333333333e0;
    double over24 = .041666666666666666667e0;

//
//     SETUP THE PARAMETERS FOR THIS CHANNEL
//
    if (L < 0) return;
    i = std::abs(2 * L - reaction.distortedWave.channel[channelIndex].twoSpin);
    if (reaction.distortedWave.channel[channelIndex].hasSpinorbit && jProj < i) return;
//
    tt1 = second();
    lDouble = L;
    lLp1 = lDouble * (lDouble + 1);
    h = reaction.distortedWave.channel[channelIndex].stepSize;
    reaction.integrationGrid.stepSize = reaction.distortedWave.channel[channelIndex].rStart;
    stepCount = reaction.distortedWave.channel[channelIndex].nGridSteps;
//     STORE J-TOTAL FOR L-DEPENDENT LINKULES
    reaction.angMom.J = 2 * L;
//
//     IF LINKULE POTENTIALS OR WAVEFUNCTIONS ARE IN USE,
//     AND IF THIS IS A DWBA, THEN WE MUST SETUP THE PROPER /LNKBLK/
//
    if (reaction.flags.problemType >= 20) {
        for (i = 1; i <= numLinkules; i++) {
            for (ii = 1; ii <= 6; ii++) {
                reaction.linkuleData.linkuleAddr[i][ii] = reaction.linkuleData.lnkAd2[i][channelIndex][ii];
            }
        }
    }
//
//     START THE INTEGRATION AT THE POINT DETERMINED FOR THE PREVIOUS
//     SOLUTION IF WE HAVE NOT RESET TO A NEW, LOWER, L VALUE.
//
    reaction.distortedWave.scatteringSolver.nFirst = reaction.distortedWave.channel[channelIndex].lastNf;
    startIndex = reaction.distortedWave.scatteringSolver.nFirst - 2;
    if (L < reaction.distortedWave.channel[channelIndex].lastL) {
        reaction.distortedWave.scatteringSolver.nFirst = 2;
        startIndex = 0;
        reaction.distortedWave.scatteringSolver.lastZero = 1;
    }
//
//     STEP1I (=stepI) IS THE SMALLEST VALUE OF THE WAVEFUNCTION TO
//     ALLOW AND THUS CONTROLS THE CHOICE OF nFirst.
//
    if (reaction.distortedWave.scatteringSolver.stepI == 1) reaction.distortedWave.scatteringSolver.stepI = 1.e-6;
//
//     OUR INITIAL GUESS IS ALWAYS stepI (STEPR was permanently 1.0).
//
    thisR = reaction.distortedWave.scatteringSolver.stepI;
    thisI = thisR;
//
//     ASSUME   U(R) = R**(L+1)  TO GET RATIO OF FIRST TWO POINTS
//
    aA1 = startIndex;
//     ZERO THE FIRST FEW POINTS IF NECESSARY.
    if (startIndex >= reaction.distortedWave.scatteringSolver.lastZero) {
        for (i = reaction.distortedWave.scatteringSolver.lastZero; i <= startIndex; i++) {
            waveReal[i] = 0.;
            waveImag[i] = 0.;
        }
    }
    d = aA1 / (aA1 + 1);
    waveReal[startIndex + 1] = thisR * std::pow(d, (L + 1));
    waveImag[startIndex + 1] = thisI * std::pow(d, (L + 1));
    waveReal[startIndex + 2] = thisR;
    waveImag[startIndex + 2] = thisI;
    thisR = thisR + thisI;
    thisRSaved = thisR;
//
//     COMPUTE THE TOTAL EFFECTIVE POTENTIAL FOR THIS L AND TEMPORARILY
//     STORE IN THE WAVEFUNCTION ARRAYS.
//     MODIFY POTENTIAL FOR MULTI-CHANNELS
//     problemType == 24 CCH-offset else branch + the `if (!force_cch_loop)`
    for (i = startIndex; i <= stepCount; i++) {
        waveReal[i + 4] = vReal[i + 1] + lLp1 * vCent[i + 1];
        waveImag[i + 4] = vImag[i + 1];
    }
//
//        MAKE ONE EXTRA SAFE
//
    waveReal[stepCount + 5] = 1;
    waveImag[stepCount + 5] = 0;
//
//     ADD IN THE SPIN ORBIT FORCES IF THEY EXIST
//
    // use soRPointer/soIPointer != nullptr as guard (works for both pool and class-owned)
    if (reaction.distortedWave.channel[channelIndex].soRPointer != nullptr || reaction.distortedWave.channel[channelIndex].soIPointer != nullptr) {
//
//     sDotL IS REALLY  sigma.L
//
    sDotL = .25 * (jProj * (jProj + 2) - reaction.distortedWave.channel[channelIndex].twoSpin * (reaction.distortedWave.channel[channelIndex].twoSpin + 2)) - lLp1;
    sDotL = sDotL / reaction.distortedWave.channel[channelIndex].twoSpin;
    if (reaction.distortedWave.channel[channelIndex].soRPointer != nullptr) {
      const double* p = reaction.distortedWave.channel[channelIndex].soRPointer;
      for (i = startIndex; i <= stepCount; i++) waveReal[i + 4] += sDotL * p[i];
    }
    if (reaction.distortedWave.channel[channelIndex].soIPointer != nullptr) {
      const double* p = reaction.distortedWave.channel[channelIndex].soIPointer;
      for (i = startIndex; i <= stepCount; i++) waveImag[i + 4] += sDotL * p[i];
    }
    }
//
//
//
//     NOW PROCESS L-DEPENDENT POTENTIALS WHICH COME ONLY FROM LINKULES.
//
    a1 = -(reaction.integrationGrid.stepSize / Constants::hbar_c) * (reaction.integrationGrid.stepSize / Constants::hbar_c)
       * ((channelIndex == 1) ? reaction.kin.redMi : reaction.kin.redMo) / 6.0;
    n = stepCount - startIndex + 1;
    rStart = startIndex * reaction.integrationGrid.stepSize;
//
    // The two channel linkule() calls (slots [1]/[2]) are byte-identical except the
    // address-slot index and the work array (waveReal/waveImag); dedup via a lambda.
    auto callPartialWaveLinkule = [&](int idx, double* waveArray) {
        if (reaction.linkuleData.linkuleAddr[idx][3] != 0 && reaction.linkuleData.linkuleAddr[idx][4] == 1) {
            linkule(reaction.linkuleData.linkuleAddr[idx][3], *reinterpret_cast<char8*>(&reaction.linkuleData.linkuleAddr[idx][1]),
                    &reaction.linkuleData.linkuleAddr[idx][5],
                    idx, 4, returnCode,
                    L, (double)jProj, rStart, reaction.integrationGrid.stepSize, n, &waveArray[startIndex + 4], &a1, (char*)&channelIndex, reaction);
            if (returnCode < 0) std::exit(7777);
        }
    };
    callPartialWaveLinkule(1, waveReal);
    callPartialWaveLinkule(2, waveImag);
//
//     COMPUTE FIRST TWO XSI'S
//
    ii = startIndex + 1;
    for (i = startIndex; i <= ii; i++) {
        waveReal[i + 3] =
            twelth * (waveReal[i + 1] * waveReal[i + 4] - waveImag[i + 1] * waveImag[i + 4]);
        waveImag[i + 3] =
            twelth * (waveReal[i + 1] * waveImag[i + 4] + waveImag[i + 1] * waveReal[i + 4]);
    }
//
    if (nPts < 0) return;
//
//     READY TO DO THE INTEGRATION.
//     IS IT DONE HERE OR IN A LINKULE?
//
    t1 = second();
    reaction.timing.times[1] = reaction.timing.times[1] + (float)(t1 - tt1);
//
    if (reaction.linkuleData.linkuleAddr[6][3] != 0) {
//
    linkule(reaction.linkuleData.linkuleAddr[6][3], *reinterpret_cast<char8*>(&reaction.linkuleData.linkuleAddr[6][1]),
            &reaction.linkuleData.linkuleAddr[6][5],
            reaction.distortedWave.scatteringSolver.nFirst, 3, returnCode, L, (double)jProj, h * startIndex, h, stepCount + 1,
            waveReal, waveImag, (char*)&channelIndex, reaction);
    if (returnCode < 0) std::exit(7777);
    } else {
//
//     WILL DO INTEGRATION HERE
//
//
//      NUMEROV METHOD (A LA RAYNAL) IS USED IN THE FOLLOWING LOOP
//
    n1 = reaction.distortedWave.scatteringSolver.nFirst;
    n2 = stepCount;
    {
    bool useFullRange = (L >= reaction.distortedWave.channel[channelIndex].lastL);
    if (!useFullRange) {
        extraSteps = 20;
        if (L > 16) extraSteps = std::exp(50.0 / lDouble);
    }
//
    while (true) {
    if (!useFullRange) {
        n2 = n1 + extraSteps;
        n2 = std::min(n2, stepCount);
    }
//
    d = (12. / (waveReal[n1 + 4] * waveReal[n1 + 4] + waveImag[n1 + 4] * waveImag[n1 + 4]));
//
    for (i = n1; i <= n2; i++) {
//
        waveReal[i + 3] = ((-10.) * waveReal[i + 2]) + (waveReal[i] - waveReal[i + 1]);
        waveImag[i + 3] = ((-10.) * waveImag[i + 2]) + (waveImag[i] - waveImag[i + 1]);
//
        waveReal[i + 1] = (waveReal[i + 4] * waveReal[i + 3] + waveImag[i + 4] * waveImag[i + 3]) * d;
        waveImag[i + 1] = (waveReal[i + 4] * waveImag[i + 3] - waveImag[i + 4] * waveReal[i + 3]) * d;
//
        d = (12. / (waveReal[i + 5] * waveReal[i + 5] + waveImag[i + 5] * waveImag[i + 5]));
    }
    if (n2 == stepCount) break;
//
//     CHECK FOR NEARNESS OF OVERFLOW
//
    thisR = std::fabs(waveReal[n2 + 1]) + std::fabs(waveImag[n2 + 1]);
    if (reaction.distortedWave.scatteringSolver.pwBgSwitch) std::printf(" CHECKING SIZE: NWP, L, N1, N2, MAG =%5d%5d%5d%5d%13.3G\n",
        channelIndex, L, n1, n2, thisR);
//
    if (thisR > Constants::bigNum) {
//
//     MUST SCALE DOWN.  WE WILL MULTIPLY EVERYTHING BY 1/thisR.
//     FIRST WE FIND WHERE THIS WILL RESULT IN |U| = stepI.
//
    thisI = reaction.distortedWave.scatteringSolver.stepI * thisR;
    n1 = startIndex;
    for (startIndex = n1; startIndex <= n2; startIndex++) {
        if (std::fabs(waveReal[startIndex + 1]) >= thisI) break;
        waveReal[startIndex + 1] = 0;
        waveImag[startIndex + 1] = 0.;
    }
//
//
//     SHOULD WE REDUCE extraSteps
//
    if (thisR * 1.e-10 > Constants::bigNum) extraSteps = std::max(extraSteps / 2, 1);
    if (reaction.distortedWave.scatteringSolver.pwBgSwitch) std::printf("     RESCALING: NEW ISTRT =%6d%6d\n", startIndex, extraSteps);
//
    thisI = 1 / thisR;
    n1 = n2 + 2;
    for (i = startIndex; i <= n1; i++) {
        waveReal[i + 1] = thisI * waveReal[i + 1];
        waveImag[i + 1] = thisI * waveImag[i + 1];
    }
    thisR = 0;
    } else if (thisR < 10000 * thisRSaved) {
        extraSteps = extraSteps + 5;
    }
//
    n1 = n2 + 1;
    thisRSaved = thisR;
    }  // end while (Numerov integration)
    }  // end useFullRange scope
    }  // end else: Numerov integration branch
//
//     END OF THE INTEGRATION
//
    t2 = second();
    reaction.timing.times[2] = reaction.timing.times[2] + (float)(t2 - t1);
//
//     EXTRACT S-MATRIX AND NORMALIZATION FROM U AND U AT NBAKCM
//     STEPS BACK FROM END.
//
// F IS THE REGULAR COULOMB WAVE
//  G IS THE IRREGULAR COULOMB WAVE
//  NOTE THAT f1,g1 ARE REGULAR AND IRREGULAR COULOMB FUNTIONS AT NBACK
//
    i = L;
    F  = reaction.distortedWave.channel[channelIndex].F_arr[i];
    G  = reaction.distortedWave.channel[channelIndex].G_arr[i];
    f1 = reaction.distortedWave.channel[channelIndex].nF1sArr[i];
    g1 = reaction.distortedWave.channel[channelIndex].nG1sArr[i];
//
//  NOW FINDING REAL S AND IMAG S
//
    a1 = waveReal[stepCount + 1] * f1 + waveImag[stepCount + 1] * g1
         - waveReal[stepCount + 1 - 4] * F
         - waveImag[stepCount + 1 - 4] * G;
    a2 = -waveReal[stepCount + 1] * g1 + waveImag[stepCount + 1] * f1
         + waveReal[stepCount + 1 - 4] * G
         - waveImag[stepCount + 1 - 4] * F;
    cReal = -waveReal[stepCount + 1] * f1 + waveImag[stepCount + 1] * g1
         + waveReal[stepCount + 1 - 4] * F
         - waveImag[stepCount + 1 - 4] * G;
    cImag = -waveReal[stepCount + 1] * g1 - waveImag[stepCount + 1] * f1
         + waveReal[stepCount + 1 - 4] * G
         + waveImag[stepCount + 1 - 4] * F;
//
//     TRANSLATE THE ELASTIC S-MATRICES TO THE lx SYSTEM AND SAVE THEM.
//
    sJr = (cReal * a1 + cImag * a2) / (a1 * a1 + a2 * a2);
    sJi = (cImag * a1 - cReal * a2) / (a1 * a1 + a2 * a2);
    reaction.boundState.convertJtoL(L, jProj, channelIndex, sJr, sJi,
           reaction.distortedWave.channel[channelIndex].indxeArr.data(),  // 1-based
           reaction.distortedWave.channel[channelIndex].smatArr.data(),  // 0-based smatArr
           reaction);
//
//  NORMALIZATION HERE
//
    a1 = .5 * (F * (1. + sJr) + sJi * G);
    a2 = .5 * (G * (1. - sJr) + sJi * F);
    cReal = (waveReal[stepCount + 1] * a1 + waveImag[stepCount + 1] * a2) /
         (waveReal[stepCount + 1] * waveReal[stepCount + 1] + waveImag[stepCount + 1] * waveImag[stepCount + 1]);
    cImag = (waveReal[stepCount + 1] * a2 - waveImag[stepCount + 1] * a1) /
         (waveReal[stepCount + 1] * waveReal[stepCount + 1] + waveImag[stepCount + 1] * waveImag[stepCount + 1]);
//
//     FIND hasNextBlock STARTING POINT
//     FIND POINT WHERE ALL RENORMALIZED WAVEFUNCTION WILL BE < stepI
//
    aA1 = std::fabs(cReal) + std::fabs(cImag);
    waveThreshold = reaction.distortedWave.scatteringSolver.stepI / aA1;
    for (i = startIndex; i <= stepCount; i++) {
        if (std::fabs(waveReal[i + 1]) + std::fabs(waveImag[i + 1]) >= waveThreshold) break;
        waveReal[i + 1] = 0.;
        waveImag[i + 1] = 0.;
    }
//
    startIndex = i;
    reaction.distortedWave.scatteringSolver.lastZero = i + 1;
//
//     SET IT BACK SOME FOR LINKULES
//
    i = std::max(i + 1, 2);
    reaction.distortedWave.channel[channelIndex].lastNf = i;
    reaction.distortedWave.channel[channelIndex].lastL = L;
    isBadRenorm = ((aA1 > Constants::bigNum * 1.e-10 && reaction.distortedWave.scatteringSolver.nFirst > 2) ||
             aA1 < Constants::smlNum * 1.e-5);
    if (isBadRenorm) std::printf("0**** WARNING, LARGE WAVEFUNCTION RENORMALIZATION:\n");
    if (reaction.distortedWave.scatteringSolver.pwBgSwitch || isBadRenorm)
        std::printf(" NWP, L =%2d%4d  NFIRST =%4d  NEW NFIRST =%4d  RENORM =%11.3E   NUM PTS:%6d%6d%6d\n",
                    channelIndex, L, reaction.distortedWave.scatteringSolver.nFirst, i, aA1, stepCount, reaction.distortedWave.channel[channelIndex].nStp2s, nPts);
    reaction.distortedWave.scatteringSolver.nFirst = startIndex;
//
//
    if (reaction.distortedWave.scatteringSolver.pwAvSwitch) {
        channelNameIndex = reaction.distortedWave.scatteringSolver.isStandalone ? 3 : channelIndex;
        if (!reaction.distortedWave.channel[channelIndex].hasSpinorbit) {
            std::printf(" %8s ELASTIC S-MATRIX FOR L =%4d:%15.5G +%13.5G I\n",
                        channelName[channelNameIndex], L, sJr, sJi);
        } else {
            std::printf(" %8s ELASTIC S-MATRIX FOR L =%4d,   JP =%4d/2:%15.5G +%13.5G I\n",
                        channelName[channelNameIndex], L, jProj, sJr, sJi);
        }
    }
//
//  STORE NORMALIZED WAVE FUNCTIONS
//
    n2 = stepCount + 2;
    for (i = startIndex; i <= n2; i++) {
        a1 = cReal * waveReal[i + 1] - cImag * waveImag[i + 1];
        waveImag[i + 1] = cReal * waveImag[i + 1] + cImag * waveReal[i + 1];
        waveReal[i + 1] = a1;
    }
    ii = startIndex + 2;
    if (reaction.distortedWave.scatteringSolver.pwBgSwitch)
        std::printf(" STARTS:%5d%14.5G%14.5G%5d%14.5G%14.5G%5d%14.5G%14.5G\n",
                    startIndex, waveReal[startIndex], waveImag[startIndex],
                    startIndex + 1, waveReal[startIndex + 1], waveImag[startIndex + 1],
                    startIndex + 2, waveReal[startIndex + 2], waveImag[startIndex + 2]);
//
//
    step2Count = reaction.distortedWave.channel[channelIndex].nStp2s;
    if (step2Count > stepCount) {
        stepCount = stepCount + 1;
        rT = stepCount * reaction.integrationGrid.stepSize;
//
        for (i = stepCount; i <= step2Count; i++) {
            rInv = 1. / rT;
            waveReal[i + 3] = (-10.) * waveReal[i + 2] + (waveReal[i] - waveReal[i + 1]);
            waveImag[i + 3] = (-10.) * waveImag[i + 2] + (waveImag[i] - waveImag[i + 1]);
            d = 12. / (reaction.distortedWave.channel[channelIndex].xFacs[1] + rInv * (reaction.distortedWave.channel[channelIndex].xFacs[2] +
                (lLp1 * reaction.distortedWave.channel[channelIndex].xFacs[3]) * rInv));
            waveReal[i + 1] = d * waveReal[i + 3];
            waveImag[i + 1] = d * waveImag[i + 3];
            rT = rT + reaction.integrationGrid.stepSize;
        }
    }
//
    if (reaction.distortedWave.scatteringSolver.isStandalone || nPts == 0) {
        reaction.timing.times[3] = reaction.timing.times[3] + (float)(second() - t2);
        return;
    }
//
//     NOW INTERPOLATE TO THE DESIRED GRID
//
    stepInverse = 1 / reaction.distortedWave.channel[channelIndex].rStart;
    srcMaxIndex = step2Count - 2;
    if (reaction.distortedWave.scatteringSolver.pwBgSwitch) std::printf(" INTERPOLATING:%6d%6d%15.5G\n", nPts, srcMaxIndex, stepInverse);
    for (rGridIndex = 1; rGridIndex <= nPts; rGridIndex++) {
        rT = rGrid[rGridIndex];
//
//     FIND THE LOCATION OF R IN THE TABLE
//
        rByH = rT * stepInverse;
        i = (int)(rByH + 0.5);
        i = std::max(2, std::min(i, srcMaxIndex));
        p = rByH - i;
//     USE 5-POINT LAGRANGIAN INTERPOLATION (ABRAMOWITZ & STEGUN 25.2.15)
        pSquared = p * p;
        x1 = p * (pSquared - 1.) * over24;
        x2 = x1 + x1;
        x3 = x1 * p;
        x4 = x2 + x2 - 0.5 * p;
        x5 = x4 * p;
        c1 = x3 - x2;
        c5 = x3 + x2;
        c3 = x5 - x3;
        c2 = x5 - x4;
        c4 = x5 + x4;
        c3 = c3 + c3 + 1.;
        waveR[rGridIndex] = (float)(c1 * waveReal[i - 1] - c2 * waveReal[i] + c3 * waveReal[i + 1]
                            - c4 * waveReal[i + 2] + c5 * waveReal[i + 3]);
        waveI[rGridIndex] = (float)(c1 * waveImag[i - 1] - c2 * waveImag[i] + c3 * waveImag[i + 1]
                            - c4 * waveImag[i + 2] + c5 * waveImag[i + 3]);
//
    }
//
    reaction.timing.times[3] = reaction.timing.times[3] + (float)(second() - t2);
    return;
}
