// l_extrapolation.cpp — L-extrapolation helper functions for LINTRP.

#include "CoulombWaveFunction.h"
#include "l_extrapolation.h"
#include "math/numeric_utils.h"
#include <cstdio>
#include <cmath>

double lxTrpM(int excitationType, double barA, double b, double barL,
              double lMaxDouble, double weeBoy)
{
    // excitationType ∈ {1, 2} since defaults.cpp seeds 1 and parameters.cpp inelastic
    // path sets 2; no other writers. Cases 3 (phase-power) and 4 (WKB form)
    double lExtrap = lMaxDouble;
    switch (excitationType) {
    case 1: // Woods-Saxon
        lExtrap = CoulombWaveFunction::solveRTXLNX(b, -0.50, std::log(weeBoy/barA) - 0.5*std::log(2.0) - b*barL, 0.010);
        break;
    case 2: // Power law
        if (b > 0.5) lExtrap = std::exp((std::log(1.41*std::fabs(barA)/weeBoy) + b*std::log(lMaxDouble)) / (b - 0.5));
        break;
    }
    return lExtrap + 0.5;
}

// callers inlined), and eta_ch fed only the case-4 body that was dropped
// the same day along with the ETA local.
void lxTrp1(int excitationType, int N, int& convergenceCode, int printLevel, double* xVals, double* sVals,
            double& flCrit, double& aVal, double& width,
            double& barL, double& barA, double& b, double& barC, double lMaxDouble,
            double& chiSq, int lx, int lDelta, int jProj, int jT)
{
    // excitationType-as-field is bounded to {1, 2}; interpolation.cpp:300 also calls
    // with literal `3` (phase-power form). Case 4 (WKB) had no field or
    int debugSwitch = printLevel >= 4;
    int linPrintSwitch = printLevel >= 6;
    double del = 1.0e-5;
    double step = 0.30;
    int ixType2 = 2*excitationType;
    convergenceCode = 0;
    // multiplier" knob was never wired by any caller).
    double bMin = 0.001;

    // The "MINIMIZATION W.R.T. WIDTH HAS FAILED / EXTRAPOLATION MUST BE
    // SUPPRESSED" warning + convergenceCode=-5 marks two bailout points (the
    // initial-B out-of-bounds check and the L800 final-B check); share it.
    // Each caller prints its own preceding diagnostic and issues the return.
    auto suppressExtrapolation = [&] {
        std::printf("\n*** MINIMIZATION W.R.T. WIDTH HAS FAILED \n       EXTRAPOLATION MUST BE SUPPRESSED --\n");
        convergenceCode = -5;
    };

    double xBar = (xVals[1]+xVals[2]) / 2.0;
    if (excitationType >= 2) xBar = lMaxDouble;

    for (int i = 1; i <= N; i++) {
        switch (excitationType) {
        case 1: xVals[i] = xVals[i] - xBar; break;
        case 2: case 3: xVals[i] = xBar / xVals[i]; break;
        }
    }

    if (debugSwitch) std::printf("\nEXTRAPOLATION FOR MCHN =%5d%5d%5d%5d%5d   IEXTYP =%5d%5d\n",
        2, jProj, jT, lx, lDelta, excitationType, N);
    if (linPrintSwitch) { std::printf("\nX AND S:\n"); for(int i=1;i<=N;i++) std::printf("%4d%17.5G%17.5G\n",i,xVals[i],sVals[i]); }

    double cVal = (excitationType == 3) ? 3 : 1;
    double bExp=0, aExp=0, cExp=0, f=0, chiExp=0;
    linLsq(ixType2-1, N, xVals, sVals, bExp, aExp, cVal, f, linPrintSwitch);
    if (excitationType == 3) {
        cExp = bExp; ixType2 = 5; bExp = 3;
    } else if (bExp <= bMin) {
        std::printf("\n***  INITIAL B = %12.6G (OUT OF BOUNDS)--\n", bExp);
        suppressExtrapolation();
        return;
    }

    bool switchSucceeded = false;

    // The "FOR CHANNEL ... JP, JT, LX, LO-LI" provenance line follows several
    // warning prints byte-identically; share it via one helper.
    auto printChannelLine = [&] {
        std::printf("      FOR CHANNEL%3d, JP, JT, LX, LO-LI =%3d/2%3d/2%4d%4d\n", 2, jProj, jT, lx, lDelta);
    };

    // Shared epilogue. Each of the
    // 8 in-body `goto L700;` short-cuts becomes
    // `{ finishWithFallback(true); return; }`; the natural-success path
    // calls with `!switchSucceeded`.
    auto finishWithFallback = [&](bool useFallback) {
        if (useFallback) {
            // L700 body
            if (excitationType == 3) {
                aVal = aExp; cVal = cExp; b = bExp;
                barA = aVal; barC = cVal; flCrit = cVal; width = b;
            } else {
                barA = std::exp(aExp); b = bExp; chiSq = chiExp; barC = 0; convergenceCode = 5;
                if (excitationType != 1) {
                    width = bExp; flCrit = 0; barL = 0; aVal = barA;
                } else {
                    barL = xBar; aVal = barA; width = 1.0/bExp; flCrit = barL;
                }
            }
        }

        // L800 body
        if (b <= bMin) {
            std::printf("\n*** ERROR IN EXTRAP: B OUT OF BOUNDS:%15.5G%15.5G%15.5G%15.5G\n", b, barA, barC, barL);
            suppressExtrapolation();
            return;
        }
        if (debugSwitch) std::printf("\nLXTRP1 END%3d%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G\n\n",
            excitationType, barA, b, barC, barL, aVal, width, cVal, flCrit, chiSq);
        switch (excitationType) {
        case 1: if (b > 0.01) return; break;
        case 2: if (b > bMin+0.6) return; break;
        case 3: return;
        }
        std::printf("\n*** WARNING:  EXTRAPOLATION DECAYS VERY SLOWLY;  A, B, C, LCRIT = %15.5G%15.5G%15.5G%15.5G\n",
            aVal, b, cVal, flCrit);
        printChannelLine();
    };

    {
    if (debugSwitch) std::printf("\nB FOR EXP OR POWER FORM = %12.6G\n", bExp);
    step = step * std::fabs(bExp/10.0);
    if (excitationType == 3) step = 0.1;
    double bMax = 10*std::fabs(bExp);
    if (excitationType == 1) bMax = 1.2;

    double x = bExp;
    double c1=0, a1=0;
    linLsq(ixType2, N, xVals, sVals, c1, a1, x, f, linPrintSwitch);
    chiExp = f;

    { // Block for del-doubling loop
    double xPlus, xMinus, fPlus, fMinus, df, derivTest;
    while (true) {
        xPlus = x*(1+del);
        xMinus = x*(1-del);
        linLsq(ixType2, N, xVals, sVals, cVal, aVal, xPlus, fPlus, linPrintSwitch);
        linLsq(ixType2, N, xVals, sVals, cVal, aVal, xMinus, fMinus, linPrintSwitch);
        df = (fPlus-fMinus) / (2*x*del);
        derivTest = std::fabs(df/f);
        if (derivTest > 1.0e-10) break;
        if (del >= 0.5) { finishWithFallback(true); return; }
        del = 2*del;
    }

    { double sV = std::fabs(step);
    if (df > 0) sV = -sV;

    int iter = 1;
    double f1 = f, x1 = x;
    double x2, f2, x3, f3, c, a, c2=0, a2=0;

    while (true) {
        x2 = x1 + sV;
        if (excitationType == 1 && x2 > 1.2) {
            std::printf("\n**** ATTEMPTING TO FIT VERY SHARP WOODS-SAXON - ATTEMPT ABANDONED\n");
            printChannelLine();
            { finishWithFallback(true); return; }
        }
        if (excitationType >= 2 && x2 < bMin) {
            std::printf("\n*** ATTEMPTING TO FIT WITH POWER OF 0 -- DEFAULT STARTING VALUE USED.\n");
            printChannelLine();
            { finishWithFallback(true); return; }
        }
        if (excitationType == 1 && x2 < bMin) { finishWithFallback(true); return; }
        if (excitationType >= 2 && x2 > bMax) { finishWithFallback(true); return; }

        linLsq(ixType2, N, xVals, sVals, c, a, x2, f2, linPrintSwitch);
        if (f2 > f1) {
            // Passed minimum — bracket it
            sV = 0.5*sV;
            x3 = x2; f3 = f2;
            x2 = x3 - sV;
            linLsq(ixType2, N, xVals, sVals, c, a, x2, f, linPrintSwitch);
            if (f < f1) { f2=f; c2=c; a2=a; break; }
            x3 = x2; f3 = f;
            x2 = x1; f2 = f1; a2 = a1; c2 = c1;
            x1 = x2 - sV;
            linLsq(ixType2, N, xVals, sVals, c, a, x1, f1, linPrintSwitch);
            break;
        }

        x1=x2; f1=f2; a1=a; c1=c;
        iter = iter+1;
        if (iter > 1000) {
            std::printf("\n*** MORE THAN 1000 STEPS TO BOX MINIMUM IN B-- \n");
            { finishWithFallback(true); return; }
        }
    }

    { double f2Pivot = f1 - 2*f2 + f3;
    double sMin, xMin;
    if (f2Pivot != 0) {
        sMin = (0.5*sV*(f1-f3)) / f2Pivot;
        xMin = x2 + sMin;
        if (xMin <= bMin || xMin > bMax) { finishWithFallback(true); return; }
        linLsq(ixType2, N, xVals, sVals, cVal, aVal, xMin, f, linPrintSwitch);
        if (f > f2) { xMin = x2; aVal = a2; cVal = c2; f = f2; }
    } else {
        xMin = x2; aVal = a2; cVal = c2; f = f2;
    }
    b = xMin; chiSq = f;

    // Taylor refinement
    { double dx = del*xMin;
    double fMin = f;
    xPlus = xMin + dx;
    linLsq(ixType2, N, xVals, sVals, c, a, xPlus, fPlus, linPrintSwitch);
    xMinus = xMin - dx;
    linLsq(ixType2, N, xVals, sVals, c, a, xMinus, fMinus, linPrintSwitch);
    f2Pivot = (fPlus - 2*fMin + fMinus) / (dx*dx);
    if (f2Pivot != 0) {
        double xNew = xMin - ((fPlus-fMinus) / (2*dx)) / f2Pivot;
        linLsq(ixType2, N, xVals, sVals, c, a, xNew, fMin, linPrintSwitch);
        if (fMin < f) { aVal=a; cVal=c; b=xNew; chiSq=fMin; }
    }
    } // end Taylor block

    { double test = cVal + aVal;
    if (debugSwitch) std::printf("\n-*  (A,C) = (%12.6G %12.6G)\n     (B,B2,CH)=(%12.6G %12.6G %12.6G)\n",
        aVal, cVal, b, x2, chiSq);
    if (test <= 0 && excitationType != 3) {
        if (debugSwitch) std::printf("\n--- BEST-FIT IS AN EXP OR SIMPLE POWER --\n");
        { finishWithFallback(true); return; }
    }
    if (debugSwitch) std::printf("\n-COMPLETE IN %6d ITERATIONS\n", iter);
    }

    // Return values per excitationType
    switch (excitationType) {
    case 1: // Woods-Saxon
        barA = 1.0/std::fabs(aVal); width = 1.0/b;
        barC = std::copysign(1.0, aVal);
        flCrit = xBar - width*std::log(barA*cVal);
        barL = flCrit; aVal = barA;
        switchSucceeded = true;
        break;
    case 2: // Power fit
        width = b;
        if (aVal >= 0) {
            barA = aVal; barC = 1;
            barL = cVal*xBar/aVal; flCrit = barL;
            switchSucceeded = true;
        }
        break;
    case 3: // Phase power
        barA = aVal; barC = cVal; flCrit = cVal; width = b;
        switchSucceeded = true;
        break;
    }
    } // end L300 block
    } // end L200 block
    } // end L130 block
    } // end L110 scope

    finishWithFallback(!switchSucceeded);
}

//
// excitationType-as-field is bounded to {1, 2}; cases 1 and 2 fire from
// interpolation.cpp:542 with the field value. Case 3 stays live because
// interpolation.cpp:550 calls lxTrp2 with a hard-coded literal `3`
// (phase-power form used inside the excitationType>=2 branch). Case 4 (WKB form)
void lxTrp2(int excitationType, double barA, double b, double barC, double barL,
            double lMaxDouble, int li, double& size)
{
    double liDouble = li;
    double x;
    switch (excitationType) {
    case 1:
        x = b*(liDouble-barL);
        if (x < -46) x = -46;
        size = barA / (barC + std::exp(x));
        return;
    case 2:
        size = barA * std::pow(lMaxDouble/liDouble, b) * (1 + barL/liDouble);
        return;
    case 3:
        size = barC + barA * std::pow(lMaxDouble/liDouble, b);
        return;
    }
}

// MAKDER/SDERIV/SECDER/DERCHK/GENBNX/GETBFC/SETFIT/SETBFC/SETBRN —

