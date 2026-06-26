// BoundState_coulombIntegrals.cpp — CLINTS: evaluates INTEGRAL(rLower,INF) DR
// fCoul fCoul / R**N (ffIntegral, fgIntegral, gfIntegral, ggIntegral) for DWBA.

#include "ptolemy_types.h"
#include "math/gauss_quadrature.h"
#include "BoundState.h"
#include "CoulombWaveFunction.h"
#include <cstdio>
#include <cmath>

// ============================================================================
//
// Evaluates INTEGRAL(rLower, INF) DR fCoul fCoul / R**N
// Four integrals: ffIntegral, fgIntegral, gfIntegral, ggIntegral
// ============================================================================

void BoundState::coulombIntegrals(double rLower, double etaIn, double etaFinal, double fkIn, double fkFinal,
            double sigIn, double sigFinal, double accuracy,
            double& ffIntegral, double& fgIntegral, double& gfIntegral, double& ggIntegral,
            double* points, double* weights, double* fA, double* fpA,
            double* gA, double* gpA, double* work,
            int rPower, int li, int lf, int termCount, int nPts,
            int& returnCode, int printLevel)
{
    double asMult = 1.0;
    double PI = 3.14159265358979300;
    double accur = 1.0e-14;
    int isOddTerm;
    int isAsymptotic;

    double rhoIn, rhoFinal, zIn, zFinal, phaseIn, phaseFinal, dzIn, dzFinal;
    double fInit, fpIn, gInit, gpIn, fFinal, fpFinal, gFinal, gpFinal;
    double rho0In, rho0Final;
    double tIn, tFinal, chi, b, f, a0, r, eps, dB;
    double delta, rValue, rTurn;
    double ci, si, zFactor, rIn, rFinal, wIn, wFinal, nTerm, nIndex, term, term2;
    double cumulativeCos = 0, cumulativeSin = 0, ciPrevious = 0, siPrevious = 0;
    double dChi, d1, c1, dF, dtIn, dtFinal, d2zIn, d2zFinal, d2Chi, d12, d2, d2F;
    double invRChi, phaseAngle, sinPhase, cosPhase;
    double factorCoulomb, even, odd, sign, aPrevious, aPrevious2, a, bPrevious;
    double c2;
    double aRatio, ciLast, siLast;
    double epsilonInitial, epsilon;
    int phaseSign, pieceCount, retryCount;
    int asymMax, asymCountIn, asymCountFinal;
    int i, j, jj, iterCount;

    if ( printLevel >= 3 )  std::printf("\n CLINTS:%3d%4d%4d%3d%3d%15.6g%15.6g%15.6g%15.6g%15.6g%15.6g%15.6g\n accuracy=%15.5g\n",
      rPower, li, lf, termCount, nPts, rLower, etaIn, etaFinal, fkIn, fkFinal, sigIn, sigFinal, accuracy);

    // RCWFN-failure diagnostic, emitted byte-identically before each of the
    // four `return;` bailouts when computeFG sets a nonzero returnCode.
    auto printRcwfnError = [&] {
        std::printf("\n *** ERROR IRET = %5d IN RCWFN ***\n"
          " L's, ETA's, RHO's =%5d%5d%20.10g%20.10g%20.10g%20.10g\n",
          returnCode, li, lf, etaIn, etaFinal, rhoIn, rhoFinal);
    };

    // The computeFG calls below all bail out the same way on error: print the
    // RCWFN diagnostic and signal the caller to return. (return stays at the
    // call site — a return inside the lambda would only exit the lambda.)
    auto bailIfRcwfnError = [&]() -> bool {
        if (returnCode != 0) {
            printRcwfnError();
            return true;
        }
        return false;
    };

    // printLevel>=3 running-total summary, emitted byte-identically before
    // both `return;` exits of the integral-accumulation epilogue.
    auto printIntegralSummary = [&] {
        if ( printLevel >= 3 )  std::printf("\n C+ S+ C- S-:%15.8g%15.8g%15.8g%15.8g\n FF FG GF GG:%15.8g%15.8g%15.8g%15.8g\n",
          ciPrevious, siPrevious, cumulativeCos, cumulativeSin, ffIntegral, fgIntegral, gfIntegral, ggIntegral);
    };

    //
    // SETUP GAUSS POINTS
    //
    gaussL( nPts, points, weights );

    retryCount = 0;
    epsilonInitial = asMult * accuracy;
    epsilon = epsilonInitial;

    //
    // InItIALIZE BELLING METHOD FOR sum
    //
    phaseSign = 1;
    ffIntegral = 0;
    fgIntegral = 0;
    gfIntegral = 0;
    ggIntegral = 0;
    delta = PI / std::max(fkIn, fkFinal);
    pieceCount = 0;
    rValue = rLower;

    //
    // GET TURNING POINT
    //
    rTurn = std::max( (etaIn + std::sqrt(etaIn*etaIn + li*(li+1))) / fkIn,
      (etaFinal + std::sqrt(etaFinal*etaFinal + lf*(lf+1))) / fkFinal );
    isAsymptotic = FALSE_F;
    if ( printLevel >= 3 )  std::printf(" rTurn, epsilonInitial, DELTA =%15.6g%15.6g%15.6g\n",
      rTurn, epsilonInitial, delta);

    //
    // Outer iteration: each pass runs phaseSign=+1 then phaseSign=-1; L700-style
    // asMult-reduction retries also re-enter from here.
    //
    bool skipInnerSetup = false;
    while (true) {

    //
    // GET COULOMB WAVE FUNCTION AND BELLING PARAMETERS AT R = rValue
    //
    while (true) {
    if (!skipInnerSetup) {
    rhoIn = fkIn * rValue;
    rhoFinal = fkFinal * rValue;

    if ( isAsymptotic ) {
        //
        // FIND STUFF FROM THE ASYMPTOTIC EXPANSION
        //
        rIn = rho0In / rhoIn;
        rFinal = rho0Final / rhoFinal;
        zIn = 1;
        zFinal = 1;
        dzIn = 0;
        dzFinal = 0;
        wIn = rIn;
        wFinal = rFinal;
        phaseIn = 0;
        phaseFinal = 0;

        nTerm = 1;
        term2 = 1;
        for (i = 2; i <= asymCountIn; i++) {
            term = wIn * fA[i];
            wIn = wIn * rIn;
            zIn = zIn + term;
            dzIn = dzIn - nTerm * term;
            if ( i > 2 )  phaseIn = phaseIn + term / (nTerm - 1);
            nTerm = nTerm + 1;
            if ( std::fabs(term) + term2 < accuracy * zIn )  break;
            term2 = std::fabs(term);
        }
        if (i > asymCountIn) i = asymCountIn;
        asymCountIn = i;
        phaseIn = rhoIn - etaIn*std::log(2*rhoIn) + sigIn - .5*PI*li - rhoIn*phaseIn;
        dzIn = dzIn / rhoIn;

        nTerm = 1;
        term2 = 1;
        for (i = 2; i <= asymCountFinal; i++) {
            term = wFinal * fpA[i];
            wFinal = wFinal * rFinal;
            zFinal = zFinal + term;
            dzFinal = dzFinal - nTerm * term;
            if ( i > 2 )  phaseFinal = phaseFinal + term / (nTerm - 1);
            nTerm = nTerm + 1;
            if ( std::fabs(term) + term2 < accuracy * zFinal )  break;
            term2 = std::fabs(term);
        }
        if (i > asymCountFinal) i = asymCountFinal;
        asymCountFinal = i;
        phaseFinal = rhoFinal - etaFinal*std::log(2*rhoFinal) + sigFinal - .5*PI*lf - rhoFinal*phaseFinal;
        dzFinal = dzFinal / rhoFinal;
        fInit = std::sin(phaseIn) / std::sqrt(zIn);
        gInit = std::cos(phaseIn) / std::sqrt(zIn);
        fpIn = zIn*gInit - (.5*dzIn/zIn)*fInit;
        gpIn = -zIn*fInit - (.5*dzIn/zIn)*gInit;
        fFinal = std::sin(phaseFinal) / std::sqrt(zFinal);
        gFinal = std::cos(phaseFinal) / std::sqrt(zFinal);
        fpFinal = zFinal*gFinal - (.5*dzFinal/zFinal)*fFinal;
        gpFinal = -zFinal*fFinal - (.5*dzFinal/zFinal)*gFinal;
    } else {
        bool useRecursion = ( rValue < 1.2*rTurn );

        if (!useRecursion) {
            //
            // GET ASYMPTOTIC EXPANSION PARAMETERS
            //
            asymMax = std::max( li, 4*termCount );
            CoulombWaveFunction::asymptoticPhase( li, etaIn, rhoIn, printLevel-4, sigIn, &zIn, &phaseIn, &dzIn,
              &fInit, &fpIn, &gInit, &gpIn, fA, gA, gpA, work,
              accuracy, asymMax, asymCountIn, returnCode );
            if ( returnCode != 0 ) useRecursion = true;
        }
        if (!useRecursion) {
            CoulombWaveFunction::asymptoticPhase( lf, etaFinal, rhoFinal, printLevel-4, sigFinal, &zFinal, &phaseFinal, &dzFinal,
              &fFinal, &fpFinal, &gFinal, &gpFinal, fpA, gA, gpA, work,
              accuracy, asymMax, asymCountFinal, returnCode );
            if ( returnCode != 0 ) useRecursion = true;
        }

        if (!useRecursion) {
            //
            // fA NOW HAS EXPANSION PARAMETERS FOR INCOMING
            // fpA HAS PARAMETERS FOR FINAL.
            //
            rho0In = rhoIn;
            rho0Final = rhoFinal;
            isAsymptotic = TRUE_F;
        } else {
            CoulombWaveFunction::computeFG(rhoIn, etaIn, li, li, fA, fpA, gA, gpA, accur, returnCode);
            if (bailIfRcwfnError()) return;
            fInit = fA[li+1];
            gInit = gA[li+1];
            fpIn = fpA[li+1];
            gpIn = gpA[li+1];
            CoulombWaveFunction::computeFG(rhoFinal, etaFinal, lf, lf, fA, fpA, gA, gpA, accur, returnCode);
            if (bailIfRcwfnError()) return;
            fFinal = fA[lf+1];
            gFinal = gA[lf+1];
            fpFinal = fpA[lf+1];
            gpFinal = gpA[lf+1];
            zIn = 1 / (fInit*fInit + gInit*gInit);
            zFinal = 1 / (fFinal*fFinal + gFinal*gFinal);
        }
    }
    tIn = fInit*fpIn + gInit*gpIn;
    tFinal = fFinal*fpFinal + gFinal*gpFinal;
    }  // end if (!skipInnerSetup)
    skipInnerSetup = false;
    chi = fkIn*zIn + phaseSign*fkFinal*zFinal;

    //
    // WE DO THE WHOLE THING NUMERICALLY
    //
    b = 1.0e-8 * fkIn * zIn;
    if ( printLevel >= 4 )  std::printf(" R, Z's, CHI, TINY, T's:%12.4g%20.10g%20.10g%14.4g%12.4g%14.4g%12.4g\n",
      rValue, zIn, zFinal, chi, b, tIn, tFinal);

    if ( std::fabs(chi) >= b && std::fabs(gInit) <= 1000 ) {
        f = 1 / ( std::sqrt(zIn*zFinal) * std::pow(rValue, rPower) );
        a0 = f / chi;

        //
        // Precheck: can we use Belling at this rValue?
        //
        r = ( rPower/rValue + 2*(fkIn*zIn*std::fabs(tIn) + fkFinal*zFinal*std::fabs(tFinal)) ) / std::fabs(chi);
        eps = epsilon;
        if ( phaseSign == -1 )  eps = eps * ( 1 + std::fabs(cumulativeCos/a0) );
        dB = r*r*r;
        if ( dB*dB <= eps )  break;   // exit rValue-stepping integration

        //
        // Not yet — integrate one cycle numerically, then retry at larger rValue.
        //
        if ( phaseSign == -1 )  delta = PI / (std::fabs(chi) + rPower/rValue);
    } else if ( std::fabs(chi) < b ) {
        if ( phaseSign == -1 )  delta = PI / (std::fabs(chi) + rPower/rValue);
    }

    pieceCount = pieceCount + 1;
    ci = 0;
    si = 0;

    for (i = 1; i <= nPts; i++) {
        r = rValue + .5*delta + .5*delta*points[i];
        rhoIn = fkIn * r;
        rhoFinal = fkFinal * r;
        zFactor = .5*delta*weights[i] / std::pow(r, rPower);
        if ( !isAsymptotic ) {
            CoulombWaveFunction::computeFG(rhoIn, etaIn, li, li, fA, fpA, gA, gpA, accur, returnCode);
            if (bailIfRcwfnError()) return;
            fInit = fA[li+1];
            gInit = gA[li+1];
            CoulombWaveFunction::computeFG(rhoFinal, etaFinal, lf, lf, fA, fpA, gA, gpA, accur, returnCode);
            if (bailIfRcwfnError()) return;
            fFinal = fA[lf+1];
            gFinal = gA[lf+1];
            if ( phaseSign == -1 ) {
                ci = ci + zFactor*(gInit*gFinal + fInit*fFinal);
                si = si + zFactor*(fInit*gFinal - gInit*fFinal);
                continue;
            }
            ffIntegral = ffIntegral + zFactor*fFinal*fInit;
            fgIntegral = fgIntegral + zFactor*fFinal*gInit;
            gfIntegral = gfIntegral + zFactor*gFinal*fInit;
            ggIntegral = ggIntegral + zFactor*gFinal*gInit;
            continue;
        }

        //
        // HERE USE ASYMPTOTIC EXPANSION
        //
        rIn = rho0In / rhoIn;
        zIn = fA[asymCountIn] * rIn;
        phaseIn = 0;
        nIndex = asymCountIn - 2;
        for (jj = 3; jj <= asymCountIn; jj++) {
            j = asymCountIn + 2 - jj;
            zIn = rIn * (fA[j] + zIn);
            phaseIn = rIn * (fA[j+1]/nIndex + phaseIn);
            nIndex = nIndex - 1;
        }
        zIn = 1 + zIn;
        phaseIn = rhoIn - etaIn*std::log(2*rhoIn) + sigIn - .5*PI*li - rho0In*phaseIn;

        rFinal = rho0Final / rhoFinal;
        zFinal = fpA[asymCountFinal] * rFinal;
        phaseFinal = 0;
        nIndex = asymCountFinal - 2;
        for (jj = 3; jj <= asymCountFinal; jj++) {
            j = asymCountFinal + 2 - jj;
            zFinal = rFinal * (fpA[j] + zFinal);
            phaseFinal = rFinal * (fpA[j+1]/nIndex + phaseFinal);
            nIndex = nIndex - 1;
        }
        zFinal = 1 + zFinal;
        phaseFinal = rhoFinal - etaFinal*std::log(2*rhoFinal) + sigFinal - .5*PI*lf - rho0Final*phaseFinal;

        if ( phaseSign == -1 ) {
            ci = ci + (zFactor/sqrt(zIn*zFinal)) * cos(phaseIn - phaseFinal);
            si = si + (zFactor/sqrt(zIn*zFinal)) * sin(phaseIn - phaseFinal);
        } else {
            fInit = std::sin(phaseIn) / std::sqrt(zIn);
            gInit = std::cos(phaseIn) / std::sqrt(zIn);
            fFinal = std::sin(phaseFinal) / std::sqrt(zFinal);
            gFinal = std::cos(phaseFinal) / std::sqrt(zFinal);
            ffIntegral = ffIntegral + zFactor*fFinal*fInit;
            fgIntegral = fgIntegral + zFactor*fFinal*gInit;
            gfIntegral = gfIntegral + zFactor*gFinal*fInit;
            ggIntegral = ggIntegral + zFactor*gFinal*gInit;
        }
    }

    if ( printLevel >= 7 ) {
        if ( phaseSign == +1 )  std::printf(" INTS:%4d%6.1f%14.6g%14.6g%14.6g%14.6g\n",
          pieceCount, rValue, ffIntegral, fgIntegral, gfIntegral, ggIntegral);
        if ( phaseSign == -1 )  std::printf(" NUMERIC C-, S-:%4d%12.4g%14.6g%14.6g\n",
          pieceCount, rValue, ci, si);
    }

    rValue = rValue + delta;
    if ( phaseSign == +1 )  continue;   // was goto L100
    cumulativeCos = cumulativeCos + ci;
    cumulativeSin = cumulativeSin + si;

    //
    // HAVE WE GONE SO FAR IN DIFFERENCE THAT THE BELLING'S
    // PIECE WON'T MATTER.
    //
    if ( std::fabs(ci) < accuracy*std::fabs(ffIntegral + .5*(cumulativeCos - ciPrevious))  &&
         std::fabs(si) < accuracy*std::fabs(gfIntegral + .5*(cumulativeSin + siPrevious)) ) {
        ffIntegral = ffIntegral + .5*( cumulativeCos - ciPrevious );
        ggIntegral = ggIntegral + .5*( cumulativeCos + ciPrevious );
        fgIntegral = fgIntegral + .5*( siPrevious - cumulativeSin );
        gfIntegral = gfIntegral + .5*( siPrevious + cumulativeSin );
        printIntegralSummary();
        return;
    }
    }  // end while (rValue-stepping integration loop)

    //
    // ALL DONE WITH NUMERIC PART OF TOTALS
    //
    if ( printLevel >= 3  &&  pieceCount != 0 ) {
        if ( phaseSign == +1 )  std::printf("\n NUMERIC VALUES FOR%4d CYCLES:\n%15.8g%15.8g%15.8g%15.8g\n",
          pieceCount, ffIntegral, fgIntegral, gfIntegral, ggIntegral);
        if ( phaseSign == -1 )  std::printf(" NUMERIC C-, S- FROM%4d CYCLES:%15.8g%15.8g\n",
          pieceCount, cumulativeCos, cumulativeSin);
    }

    //
    // FOR sum WE USE THE FF INTEGRAL AS A GUIDE TO ERRORS
    //
    if ( phaseSign != -1 ) {
        cumulativeCos = ffIntegral;
        cumulativeSin = ffIntegral;
    }

    //
    // BELLINGS METHOD FOR sum OR DIF
    //
    phaseIn = std::atan2(fInit, gInit);
    phaseFinal = std::atan2(fFinal, gFinal);
    dzIn = -2*fkIn*zIn*zIn*tIn;
    dzFinal = -2*fkFinal*zFinal*zFinal*tFinal;
    dChi = fkIn*dzIn + phaseSign*fkFinal*dzFinal;
    d1 = dzIn/zIn + dzFinal/zFinal;
    c1 = static_cast<double>(rPower)/rValue + .5*d1;
    dF = -f*c1;
    wIn = 1 - (2*etaIn)/rhoIn - (li*(li+1))/(rhoIn*rhoIn);
    wFinal = 1 - (2*etaFinal)/rhoFinal - (lf*(lf+1))/(rhoFinal*rhoFinal);
    dtIn = fkIn*(fpIn*fpIn + gpIn*gpIn - wIn/zIn);
    dtFinal = fkFinal*(fpFinal*fpFinal + gpFinal*gpFinal - wFinal/zFinal);
    d2zIn = -2*fkIn*zIn*(2*dzIn*tIn + zIn*dtIn);
    d2zFinal = -2*fkFinal*zFinal*(2*dzFinal*tFinal + zFinal*dtFinal);
    d2Chi = fkIn*d2zIn + phaseSign*fkFinal*d2zFinal;
    d12 = (dzIn/zIn)*(dzIn/zIn) + (dzFinal/zFinal)*(dzFinal/zFinal);
    d2 = d2zIn/zIn + d2zFinal/zFinal;
    d2F = (f/rValue - dF)*c1 - (f*d1)/(2*rValue) - .5*f*(d2 - d12);
    if ( printLevel >= 4 )  std::printf("\n R =%12.4g   RHO's =%12.4g%12.4g   EPSILON =%10.2g\n",
      rValue, rhoIn, rhoFinal, r);
    dB = 0;
    invRChi = 1 / (rValue * chi);
    phaseAngle = phaseIn + phaseSign*phaseFinal;
    sinPhase = std::sin(phaseAngle);
    cosPhase = std::cos(phaseAngle);

    //
    // EACH ITERATION ON nItS BRINGS IN A HIGHER DERIVATIVE OF BETA.
    //
    if ( printLevel >= 6 )  std::printf("\n iterCount  K%6s B%9s DB%10s A%13s even%13s odd\n", "", "", "", "", "");

    bool needRetry = false;
    for (iterCount = 1; iterCount <= 2; iterCount++) {

    factorCoulomb = f / chi;

    //
    // BELLING'S EXPANSION
    //
        even = a0;
        aPrevious = a0;
        b = 1;
        odd = 0;
        sign = +1.0;
        isOddTerm = TRUE_F;
        bool converged = false;
        for (i = 1; i <= termCount; i++) {
            factorCoulomb = factorCoulomb * invRChi;
            c2 = dF/f - (i*dChi)/chi;
            bPrevious = b;
            b = rValue*dB + (c2*rValue - i + 1)*b;
            if ( iterCount > 1 )
              dB = (c2*rValue - i + 2)*dB + (c2 - rValue*((dF/f)*(dF/f)
                - i*(dChi/chi)*(dChi/chi)) + rValue*(d2F/f - i*d2Chi/chi))*bPrevious;
            a = b * factorCoulomb;

            if ( isOddTerm ) {
                odd = odd + sign*a;
                sign = -sign;
            } else {
                even = even + sign*a;
            }
            isOddTerm = !isOddTerm;

            if ( printLevel >= 6 )  std::printf("%3d%4d%11.4g%11.4g%11.4g%17.9g%17.9g\n",
              iterCount, i, b, dB, a, even, odd);
            if ( i != 1 ) {
              if ( std::fabs(a) > std::fabs(aPrevious2) )  { needRetry = true; break; }
              if ( std::fabs(a) < accuracy*(std::fabs(a0) + std::min(std::fabs(cumulativeCos),
                      std::fabs(cumulativeSin))) ) { converged = true; break; }
            }
            aPrevious2 = aPrevious;
            aPrevious = a;
        }
        if (needRetry) break;
        if (!converged) { needRetry = true; break; }

        ci = -even*sinPhase - odd*cosPhase;
        si = -odd*sinPhase + even*cosPhase;
        aRatio = a / a0;
        if ( printLevel >= 4 )  std::printf(" iterCount =%2d   A(%2d)/A(0) =%9.2g   COS, SIN =%15.8g%15.8g\n",
          iterCount, i, aRatio, ci, si);
        if ( iterCount != 1 ) {
            even = std::fabs( (ciLast - ci) / (cumulativeCos + ci) );
            odd  = std::fabs( (siLast - si) / (cumulativeSin + si) );
            if ( printLevel >= 4 )  std::printf("           REL. ERR. =%10.2g%10.2g\n\n", even, odd);
        }
        ciLast = ci;
        siLast = si;
    }

    //
    // WAS THE DERIVATIVE A BIG CORRECTION?
    //
    if ( !needRetry && std::max(even, odd)*std::max(even, odd) > accuracy )  needRetry = true;

    if (needRetry) {
        //
        // WE NEED TO MAKE asMult SMALLER
        //
        r = std::min( .5*epsilon, .75*std::pow(r, 6)*(epsilon/eps) );
        asMult = asMult * r / epsilon;
        epsilon = r;
        retryCount = retryCount + 1;
        if ( retryCount > 20 ) {
            std::printf("\n *** COULD NOT GET CONVERGENCE IN CLINTS *** \n");
            returnCode = -1;
            return;
        }
        if ( printLevel >= 4 )  std::printf("\n +++ REDOING CALCULATION WITH NEW asMult:%12.4g%12.4g\n",
          asMult, epsilon);
        skipInnerSetup = true;
        continue;   // — retry with new epsilon/asMult
    }

    retryCount = 0;
    if (phaseSign == -1) {
        cumulativeCos = cumulativeCos + ci; cumulativeSin = cumulativeSin + si;
        ffIntegral = ffIntegral + .5*( cumulativeCos - ciPrevious );
        ggIntegral = ggIntegral + .5*( cumulativeCos + ciPrevious );
        fgIntegral = fgIntegral + .5*( siPrevious - cumulativeSin );
        gfIntegral = gfIntegral + .5*( siPrevious + cumulativeSin );
        printIntegralSummary();
        return;
    }
    ciPrevious = ci;
    siPrevious = si;

    //
    // NOW DO  PH(IN)-PH(OUT)
    //
    pieceCount = 0;
    phaseSign = -1;
    cumulativeSin = 0;
    cumulativeCos = 0;
    epsilon = epsilonInitial;
    skipInnerSetup = true;
    continue;   // — second pass with phaseSign=-1

    }  // end outer while

}
