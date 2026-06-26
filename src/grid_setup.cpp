// grid_setup.cpp — GRDSET: sets up integration grids for DWBA transfer integrals.

#include "Timing.h"
#include "math/linear_algebra.h"
#include "math/numeric_utils.h"
#include "print_utils.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Reaction.h"
#include "Constants.h"

// The "INCREASE ASYMPTOPIA TO MORE THAN" advisory follows three separate
// "wavefunctions needed beyond asymptopia" warnings byte-identically (the
// source-level string was split at slightly different spots, but the
// concatenated format and args are the same); share it via one helper.
static void printAsymptopiaAdvice(double rPMax, double rTMax) {
    std::printf(" **** IN FUTURE RUNS INCREASE ASYMPTOPIA TO MORE THAN:"
                "%15.2f (PROJECTILE)%15.2f (TARGET)\n", rPMax, rTMax);
}

// The "VALUES SET TO MAXIMUM" diagnostic echoes the offending U/Vmin/Vmax
// then the clamp notice; it follows both the direct and interpolated V-range
// illegal-bounds checks below (only the leading message differs). Share it.
static void printValuesSetToMax(double U, double vMin, double vMax) {
    print_G(15, 5, U); print_G(15, 5, vMin); print_G(15, 5, vMax);
    std::printf("\n *** VALUES SET TO MAXIMUM ALLOWED.\n");
}

void DWBAGrid::gridSet(int& returnCode, Reaction& reaction) {
    // =========================================================================
    // Local variable declarations
    // =========================================================================

    int isFirstPass;
    int isAllZero;
    int printSwitch;
    int gridPrintSwitch;
    int skipPolynomials;   // true ⇒ V-polynomials are NOT used

    // Local scalars
    double tStart, tIntrp;
    int    interpCount;
    int    verbosity;
    double temp, time1;
    double fifo, rP, rT;
    double U, uMax = 0, vVal, vMax = 0, vMin, vMid;
    double d, syne;
    double wvwMax, rvrLimit;
    double uStep, uLimit;
    double sum0, sum1, sum2;
    double pvpMax;
    double vAtPeak = 0;
    double wow;
    int    inAdditional;
    int    n, i, ii, vIndex, uIndex;
    int    riRoIndex, vPtCount, xEnd;
    int    phiCount;
    int    phiBlockCount;
    int    formType, aitkenCount, iiMax;
    int    bumpCount;
    int    riRoBlockCount, riRoHCount, nRiRoInterp, polysCount, vTermCount;
    double vStepSize, xQuadStep;
    double x, x0, x0Min, x0Av, phi0 = 0, phi, dPhi;
    double rI, rO, area;
    double rTMax = 0, rPMax = 0;
    // uMax/VMAX/rPMax/rTMax/vAtPeak/phi0 inits silence pre-existing
    // -Wmaybe-uninitialized: each is set only inside a "first time the
    // observed max grew" conditional and read by a later debug printf;
    // gcc can't prove the conditional ever fires on the first pass.
    double vRangeMin, vRangeMax, vAbsMax;
    double vRange[4];   // 1-based, [1..3]

    std::vector<double> vPhiPointsVector, vPhiWeightsVector;

    std::vector<double> vDifPointsVector, vDifWeightsVector;   // DIFPTS/DIFWTS: size nPhiDifference
    std::vector<double> vPolysVector;            // VPOLYS: size polysCount = 3*vPolyDegree+3
    std::vector<double> vWorkVector;            // VPOLYWRK: size 4*nPhiSum+3+(vPolyDegree+1)^2
    std::vector<double> rIHVector, rOHVector;       // RIH/ROH: size nRiRoH_h+1 (1-based)
    // LOGIC replaced by local vector (int4 packed array, 1-based [1..riRoHCount])
    std::vector<int> logicVector;  // LOGIC: size riRoHCount+1 (index [0] unused)
    // Sub-array pointers into vWorkVector (set after resize)
    double *vWtsPointer = nullptr, *vResdPointer = nullptr, *vAmatPointer = nullptr;
    // Shifted pointers for 1-based element access
    double *vPolysPointer = nullptr;  // vPolysPointer[k] = vPolysVector[k-1], k>=1
    double *vDifptPointer = nullptr;  // vDifptPointer[k] = vDifPointsVector[k-1], k>=1
    double *vDifWtPointer = nullptr;  // vDifWtPointer[k] = vDifWeightsVector[k-1], k>=1

    double* rihPointer = nullptr;   // 0-based: rihPointer[k-1] = rIHVector[k-1]
    double* rohPointer = nullptr;   // 0-based: rohPointer[k-1] = rOHVector[k-1]
    // phi angle arrays + Stage5 arrays declared at function scope to avoid jump-over-init
    float*  phiTPointer  = nullptr;
    float*  phiPPointer  = nullptr;
    float*  phiPointer   = nullptr;
    float*  trapWeightPointer   = nullptr;
    int*    logicPointer  = nullptr;
    double* vPhiPtPointer  = nullptr;
    double* vPhiWtPointer  = nullptr;
    // Stage 5 arrays
    double* smivl5Pointer = nullptr;
    double* sumpt5Pointer = nullptr;
    float*  riPointer    = nullptr;
    float*  roPointer    = nullptr;
    float*  wioPointer   = nullptr;

    // Polynomial/interpolation counts
    int    nRiRoH_h;   // for H-search (= nPhiDifference*nPhiSum)

    // 1-based: xs[1]=1.0, xs[2]=0.0
    double xs[3] = { 0.0, 1.0, 0.0 };   // index 0 unused; xs[1]=1, xs[2]=0

    // VS(5) — work array for 5 test V values
    double vs[6] = {};   // 1-based

    double& undefValue      = reaction.internalState.undefValue;
    int&    stripPickup     = reaction.internalState.stripPickup;

    double& sumMax     = reaction.integrationGrid.sumMax;
    double& sumMid     = reaction.integrationGrid.sumMid;
    double& sumMin     = reaction.integrationGrid.sumMin;
    double& gammaSum     = reaction.rxn.gammaSum;
    double& gammaDif     = reaction.rxn.gammaDif;
    double& dwCutoff      = reaction.integrationGrid.dwCutoff;
    double& phiMid     = reaction.integrationGrid.phiMid;
    constexpr double gammaPhi = 1.0e-6;
    double& alphaP     = reaction.boundState.vertex[1].alpha;  // asymptotic decay constant for projectile
    double& alphaT     = reaction.boundState.vertex[2].alpha;  // asymptotic decay constant for target

    int&    printLevel     = reaction.flags.printLevel;
    int&    nPhiSum      = reaction.gridData.nPhiSum;
    int&    nPhiDifference      = reaction.gridData.nPhiDifference;
    int&    nPhiPoints      = reaction.gridData.nPhiPoints;
    int&    nPhiAdditional     = reaction.gridData.nPhiAdditional;
    constexpr int vPolyDegree = 3;
    constexpr int mapSum = 2;
    constexpr int mapDif = 1;
    constexpr int mapPhi = 2;


    // 1-based pointers into class-owned sctmnArr/sctcrArr.
    // Recompute these *just before each evaluateFormFactor call* — vector reallocations between
    // calls (boundState.setupFormFactors, allocateFormFactor, etc.) can invalidate cached ptrs.
    auto sctmn1b = [&reaction]() {
        auto& v = reaction.boundState.data.sctmnArr;
        return v.empty() ? nullptr : v.data() - 1;
    };
    auto sctcr1b = [&reaction]() {
        auto& v = reaction.boundState.data.sctcrArr;
        return v.empty() ? nullptr : v.data() - 1;
    };
    double riMax = 0;
    double roMax = 0;

    // Per-vertex: vertex[1].rLMax (projectile) / vertex[2].rLMax (target)
    double& rLMaxProj  = reaction.boundState.vertex[1].rLMax;
    double& rLMaxTarget   = reaction.boundState.vertex[2].rLMax;
    double& boundMxProj     = reaction.boundState.vertex[1].boundMx;
    double& boundMxTarget     = reaction.boundState.vertex[2].boundMx;
    double& rOfMax     = reaction.boundState.data.rOfMax;
    int     isInfoPrint;
    float&  phiSign     = reaction.boundState.data.phiSign;
    double& s1         = reaction.boundState.data.s1;
    double& s2         = reaction.boundState.data.s2;
    double& t1         = reaction.boundState.data.t1;
    double& t2         = reaction.boundState.data.t2;

    double& jacobian      = reaction.gridData.jacobian;
    int&    nInterpPoints     = reaction.gridData.nInterpPoints;
    // write-only self-assignments to gridData.iRi/IRO/IWIO (already 0
    // sentinels set by allocateRiRoWio); no reader anywhere.
    int&    nRiRoInterp_g   = reaction.gridData.nRiRoInterp;
    int&    maxCount     = reaction.gridData.maxCount;

    int&    nWaveF     = reaction.distortedWave.scatteringSolver.nWaveF;

    const double massRatios1 = reaction.boundState.vertex[1].massRatio;
    const double massRatios2 = reaction.boundState.vertex[2].massRatio;



    constexpr int facFr4 = 2;

    // st1, st2 — local temporaries for 2*s1*t1, 2*s2*t2
    double st1, st2;

    // =========================================================================
    //
    //     tIntrp = 0
    //     interpCount = 0
    //     returnCode = 0
    //
    // =========================================================================
    tStart  = second();
    tIntrp  = 0.0;
    interpCount  = 0;
    // explicitly sets returnCode = 0 before its early return (single case at
    // line 429: sumMin==undefValue) or reaches the end and sets returnCode = 1.

    verbosity  = printLevel % 10;
    gridPrintSwitch = (verbosity >= 5);
    isInfoPrint = (verbosity >= 4);
    printSwitch = (verbosity >= 1);

    // leading '0' is carriage control: blank line before content
    if (isInfoPrint) {
        std::printf("\n\n0SETUP OF THE INTEGRATION GRIDS:\n");
    }


    inAdditional = (reaction.flags.nuConL >= 2) ? 3 : 2;


    // =========================================================================
    //  THE CODE IS BASED ON THE ASSUMPTION THAT THE INCOMING PROJECTILE
    //  IS COMPOSITE ( A = B+X ).  ...
    //
    //      temp = 1 / ( mass_ratios(1) + mass_ratios(2) * (1+mass_ratios(1)) )
    //      S1 = ( 1 + mass_ratios(1) ) * ( (1 + mass_ratios(2) ) * temp )
    //      t1 = - ( 1 + mass_ratios(2) ) * temp
    //      S2 = ( 1 + mass_ratios(1) ) * temp
    //      T2 = -S1
    //      jacobian = S1**3
    //      phiSign = 1
    // =========================================================================
    temp = 1.0 / (massRatios1 + massRatios2 * (1.0 + massRatios1));
    s1   = (1.0 + massRatios1) * ((1.0 + massRatios2) * temp);
    t1   = -(1.0 + massRatios2) * temp;
    s2   = (1.0 + massRatios1) * temp;
    t2   = -s1;
    jacobian = s1 * s1 * s1;
    phiSign = 1.0f;
    if (stripPickup != 1) {
        // Pickup: B = A+X, flip coordinate mapping
        t2 =  s2;
        s2 = -s1;
        s1 =  t1;
        t1 = -s2;
        jacobian  = t1 * t1 * t1;
        phiSign = -1.0f;
    }
    if (isInfoPrint) {
        std::printf("\n0R TO RX MAPPING PARAMETERS (S1, T1, S2, T2) =%15.5G%15.5G%15.5G%15.5G     jacobian =%15.5G\n",
                    s1, t1, s2, t2, jacobian);
    }

    st1 = 2.0 * s1 * t1;
    st2 = 2.0 * s2 * t2;

    // =========================================================================
    //  INITIALIZE THE COMPUTATION OF  (PHI V PHI)
    //
    // =========================================================================
    reaction.boundState.setupFormFactors(reaction);


    constexpr int lookst = 250;
    xQuadStep   = 2.0 / ((double)lookst * (double)lookst);
    xs[2] = 1.0 - xQuadStep;

    // =========================================================================
    //  FIND THE MAXIMUM  ( R PSI PHI V PHI PSI R )   THAT WE WILL
    //  ENCOUNTER AND USE THIS TO SCALE dwCutoff
    // =========================================================================

    wvwMax = 0.0;
    U      = 0.5 * rOfMax;
    n      = (int)(1.5 * rOfMax / 0.20 + 1.0);
    vs[3]  = 0.0;
    time1  = second();

    for (uIndex = 1; uIndex <= n; uIndex++) {
        //  WE LOOK AT 5 DIFFERENT V'S FOR EACH U
        //         D = 2*(RLMAXS(2)-(S1+t1)*U)/(S1-t1)
        //         D = 2*(RLMAXS(1)-(S2+T2)*U)/(S2-T2)
        d = 2.0 * (rLMaxTarget - (s1 + t1) * U) / (s1 - t1);
        if (std::fabs(d) > 2.0 * U) d = std::copysign(2.0 * U, d);
        vs[1] = d;
        vs[2] = 0.5 * d;
        d = 2.0 * (rLMaxProj - (s2 + t2) * U) / (s2 - t2);
        if (std::fabs(d) > 2.0 * U) d = std::copysign(2.0 * U, d);
        vs[4] = 0.5 * d;
        vs[5] = d;

        for (vIndex = 1; vIndex <= 5; vIndex++) {
            vVal = vs[vIndex];
            for (i = 1; i <= 2; i++) {
                if (reaction.boundState.evaluateFormFactor(fifo, 3, U + 0.5 * vVal, U - 0.5 * vVal, xs[i],
                       sctmn1b(), 1, rP, rT, reaction) != 0) continue;
                if (std::fabs(fifo) <= wvwMax) continue;
                wvwMax = std::fabs(fifo);
                uMax   = U;
                vMax   = vVal;
                rPMax  = rP;
                rTMax  = rT;
            }
        }
        //         U = U + .20D0
        U = U + 0.20;
    }   // end DO 269

    tIntrp = tIntrp + second() - time1;
    interpCount = interpCount + 10 * n * (inAdditional + 2);

    //     1    ' sum, DIF, RP, RT =', 4F8.2, 5X, 'VALUE =', G13.3 )
    if (isInfoPrint) {
        std::printf("\n0MAX(RI PSI PHI V PHI PSI RO) OCCURS AT SUM, DIF, RP, RT =%8.2f%8.2f%8.2f%8.2f     VALUE =%13.3G\n",
                    uMax, vMax, rPMax, rTMax, wvwMax);
    }

    //      rvrLimit = dwCutoff*wvwMax
    rvrLimit = dwCutoff * wvwMax;

    // =========================================================================
    //  STEP OUT FROM  U = 0  LOOKING FOR sumMin.
    //  WE LOOK FOR U SUCH THAT
    //    PSI(lMin) PHI' V PHI' PSI(lMin)
    //  EXCEEDS VZERO.
    // =========================================================================

    //      uStep = .20D0
    //      U = 0
    //      time1 = second()
    uStep  = 0.20;
    U      = 0.0;
    time1  = second();

    while (true) {
        bool foundSumMin = false;
        for (i = 1; i <= 2; i++) {
            if (reaction.boundState.evaluateFormFactor(fifo, 4, U, U, xs[i], sctmn1b(), 1, rP, rT, reaction) != 0) {
                std::printf("\n0**** ASYMPTOPIA IS FAR TOO SMALL - CANNOT FIND SUMMIN.\n");
                U = U - uStep;
                if (sumMin == undefValue) { returnCode = 0; return; }
                foundSumMin = true; break;
            }
            if (std::fabs(fifo) >= rvrLimit) { foundSumMin = true; break; }
        }
        if (foundSumMin) break;
        interpCount = interpCount + 2 * (inAdditional + 2);
        U      = U + uStep;
    }

    tIntrp = tIntrp + second() - time1;
    U = std::fmax(0.0, U - uStep);
    if (sumMin == undefValue) sumMin = U;
    if (isInfoPrint) {
        std::printf("\n0SUMMIN:  PROPOSED =%14.5G     USED =%14.5G\n", U, sumMin);
    }

    // =========================================================================
    //  NOW FIND sumMax USING
    //    PSI(L CRIT) PHI V PHI PSI(L CRIT)
    //  AND COMPUTE MOMENTS OF THE WEIGHT FUNCTION.
    // =========================================================================

    //      sum0 = 0
    //      sum1 = 0
    //      sum2 = 0
    //      pvpMax = 0
    sum0   = 0.0;
    sum1   = 0.0;
    sum2   = 0.0;
    pvpMax = 0.0;

    //     1  5X, 'V OF MAX', 4X, 'RP OF MAX', 4X, 'rT OF MAX' / )
    if (verbosity >= 7) {
        std::printf("\n0      U         SUM(F'S)      MAX(F)     V OF MAX    RP OF MAX    RT OF MAX\n\n");
    }

    //      U = sumMin
    //      uLimit = std::max( rScts(1)+10*aScts(1), rOfMax+5 )
    U    = sumMin;
    uLimit = std::fmax(reaction.kin.rScts[1] + 10.0 * reaction.kin.aScts[1], rOfMax + 5.0);
    time1 = second();

    while (true) {
    temp   = 0.0;
    wvwMax = 0.0;
    isAllZero = true;

    //  WE LOOK AT 5 DIFFERENT V'S FOR EACH U
    //  AT EACH U, V PAIR WE LOOK AT 2 X'S
    //  FOR THIS SEARCH, VIOLATIONS OF BOUND-STATE ASYMPTOPIA ARE
    //  IGNORED (TREATED AS ZERO) UNLESS ALL TEST POINTS ARE BAD.
    {
        d = 2.0 * (rLMaxTarget - (s1 + t1) * U) / (s1 - t1);
        if (std::fabs(d) > 2.0 * U) d = std::copysign(2.0 * U, d);
        vs[1] = d;
        vs[2] = 0.5 * d;
        d = 2.0 * (rLMaxProj - (s2 + t2) * U) / (s2 - t2);
        if (std::fabs(d) > 2.0 * U) d = std::copysign(2.0 * U, d);
        vs[4] = 0.5 * d;
        vs[5] = d;

        for (vIndex = 1; vIndex <= 5; vIndex++) {
            vVal = vs[vIndex];
            for (i = 1; i <= 2; i++) {
                if (reaction.boundState.evaluateFormFactor(fifo, 3, U + 0.5 * vVal, U - 0.5 * vVal, xs[i],
                       sctcr1b(), 1, rP, rT, reaction) != 0) continue;
                temp = temp + std::fabs(fifo);
                if (wvwMax > std::fabs(fifo)) continue;
                wvwMax = std::fabs(fifo);
                vMax   = vVal;
                rPMax  = rP;
                rTMax  = rT;
                isAllZero = false;
            }
        }
    }

    if (!isAllZero) {
        //         interpCount = interpCount + 10*(inAdditional+2)
        interpCount = interpCount + 10 * (inAdditional + 2);
        sum0   = sum0 + temp;
        sum1   = sum1 + temp * U;
        sum2   = sum2 + temp * U * U;

        if (pvpMax <= temp) {
            pvpMax = temp;
            uMax   = U;
        }

        if (verbosity >= 7) {
            std::printf("%13.3G%13.3G%13.3G%13.3G%13.3G%13.3G\n",
                        U, temp, wvwMax, vMax, rPMax, rTMax);
        }

        if (!(wvwMax >= rvrLimit || U < uLimit) && rP > rLMaxProj && rT > rLMaxTarget) break;  // was goto L380
        U = U + uStep;
        continue;  // was goto L350
    }


    //     1   'BY THE BOUND STATE ASYMPTOPIA.' /
    //     2   '0**** IN FUTURE RUNS INCREASE ONE OR BOTH ASYMPTOPIA.' )
    std::printf("\n0**** WARNING:  SEARCH FOR SUMMAX WAS STOPPED"
                "BY THE BOUND STATE ASYMPTOPIA.\n"
                "0**** IN FUTURE RUNS INCREASE ONE OR BOTH ASYMPTOPIA.\n");
    U = U - uStep;
    break;
    }  // end while (sumMax search)

    tIntrp = tIntrp + second() - time1;
    if (sumMax == undefValue) sumMax = U;
    if (isInfoPrint) {
        std::printf(" SUMMAX:  PROPOSED =%14.5G     USED =%14.5G\n", U, sumMax);
    }

    // =========================================================================
    //  COMPUTE sumMid AS THE FIRST MOMENT OF U.
    // =========================================================================
    //      U = sum1/sum0
    //      sum2 = sum2/sum0
    //      sum2 = std::sqrt(sum2)
    U    = sum1 / sum0;
    sum2 = sum2 / sum0;
    sum2 = std::sqrt(sum2);
    if (sumMid == undefValue) sumMid = U * reaction.integrationGrid.midpointFactor;

    //     2  ' SUMMID:  PROPOSED =', G14.5, 5X, 'USED =', G14.5 )
    if (isInfoPrint) {
        std::printf("\n0<SUM> =%14.5G     SQRT<SUM**2> =%14.5G     SUM(MAX term) =%14.5G\n"
                    " SUMMID:  PROPOSED =%14.5G     USED =%14.5G\n",
                    U, sum2, uMax, U, sumMid);
    }

    // =========================================================================
    //  ADJUST sumMin, sumMid TO REQUIREMENTS OF CUBIC MAP.
    //
    //      sumMid = std::min( sumMid, .5*(sumMin+sumMax) )
    // =========================================================================
    sumMin = std::fmin(sumMin, 7.0 * (sumMid - sumMax / 7.0) / 6.0);
    sumMid = std::fmin(sumMid, 0.5 * (sumMin + sumMax));
    if (isInfoPrint) {
        std::printf(" FINAL CHOICE OF SUMMIN, SUMMID, SUMMAX =%15.5G%15.5G%15.5G\n",
                    sumMin, sumMid, sumMax);
    }

    // =========================================================================
    //  sumMin AND sumMax NOW KNOWN; PERFORM SINH MAP OF @sumMin,sumMax@
    //  TO @-1,1@ ; COMPUTE AND STORE GAUSS-LEGENDRE PTS
    //     FOR BOTH THE H'S AND THE INTERPOLATED INTEGRAL GRID
    // =========================================================================

    reaction.dwbaGrid.allocateSmhpts(nPhiSum, reaction);
    reaction.dwbaGrid.allocateSmhwk(nPhiSum, reaction);

    //  THE NUMBER OF POINTS FOR THE INTERPOLATED INTEGRAL GRID
    //  IS BASED ON THE WAVE-LENGTH
    //
    //      n_interp_pts = (sumMax-sumMin)*sumDensity*(akIn+akOut)/(4*PI)
    //      n_interp_pts = std::max( n_interp_pts, nPhiSum )
    {
        double& sumDensity = reaction.integrationGrid.sumDensity;
        nInterpPoints = (int)((sumMax - sumMin) * sumDensity * (reaction.kin.akIn + reaction.kin.akOut) / (4.0 * Constants::PI));
        nInterpPoints = std::max(nInterpPoints, nPhiSum);
    }

    skipPolynomials = (nInterpPoints == nPhiSum);

    //     1  6X, 'BUT DUE TO LIMITATIONS OF THE SAVEHS RUN',
    //     2    ' IT WILL BE LIMITED TO', I5 )
    //      n_interp_pts = nPhiSum
    if (skipPolynomials && nPhiSum != nInterpPoints) {
        std::printf("\n0**** WARNING:  n_interp_pts WAS COMPUTED TO BE%5d\n"
                    "      BUT DUE TO LIMITATIONS OF THE SAVEHS RUN IT WILL BE LIMITED TO%5d\n",
                    nInterpPoints, nPhiSum);
        nInterpPoints = nPhiSum;
    }
    reaction.dwbaGrid.allocateSmipts(nInterpPoints, reaction);
    reaction.dwbaGrid.allocateSmivl(nInterpPoints, reaction);

    //  GENERATE GAUSS POINTS FOR U
    //      CALL cubMap ( MAPSUM,  sumMin, sumMid, sumMax, gammaSum,
    cubMap(mapSum, sumMin, sumMid, sumMax, gammaSum,
           reaction.gridData.smiptsPointer, reaction.gridData.smivlPointer, nInterpPoints);

    nRiRoInterp  = nInterpPoints * nPhiDifference;
    nRiRoInterp_g = nRiRoInterp;
    riRoBlockCount   = (nRiRoInterp + facFr4 - 1) / facFr4;
    // iRi/IRO/IWIO → DWBAGrid class member float vectors (RIPTS/ROPTS/RIROWTS)
    reaction.dwbaGrid.allocateRiRoWio(riRoBlockCount, reaction);
    // gridData.IWIO/iRi/IRO are 0 sentinels; class-owned RIPTS/ROPTS/RIROWTS hold the storage.
    vDifPointsVector.assign(nPhiDifference, 0.0);
    vDifWeightsVector.assign(nPhiDifference, 0.0);
    polysCount = 3 * vPolyDegree + 3;
    vPolysVector.assign(polysCount, 0.0);
    vWorkVector.assign(4 * nPhiSum + 3 + (vPolyDegree + 1) * (vPolyDegree + 1), 0.0);

    // =========================================================================
    //  IF USING H'S, RESTORE GENERATED sum POINTS AND POLYS
    //
    //      READ ( 1, END=9950, ERR=9950 )
    //      ii = 3*nPhiSum
    // =========================================================================
    {

    cubMap(mapSum, sumMin, sumMid, sumMax, gammaSum,
           reaction.gridData.smhptsPointer, reaction.gridData.smhwkPointer, nPhiSum);

    // =========================================================================
    //  SETUP THE DIFFERENCE GRID.
    //
    //  THE STAGES ARE
    //  1)  DETERMINE THE RANGE OF V FOR EACH VALUE OF U.
    //  2)  INSIDE EACH V-RANGE, FIND THE MOMENTS OF V AND THUS
    //      THE "MID-POINT".
    //  3)  MAKE A LOW-DEGREE POLYNOMIAL FIT TO THE MID- AND END-
    //      POINTS.
    //  4)  GENERATE ALL THE rI, rO PAIRS FOR FINDING H'S
    //  5)  GENERATE ALL rI, rO PAIRS FOR THE INTERPOLATED INTEGRATION
    //      GRID.
    //
    //  IN ALL THESE OPERATIONS WE DO NOT INCLUDE THE SCATTERING
    //  WAVEFUNCTIONS SINCE  V << U  SO THAT rI, rO WILL VARY ONLY
    //  A LITTLE FOR FIXED U.
    //
    //  IN FINDING THE EXTREMA OF THE V'S, WE USE
    //    PHI' V PHI'
    // =========================================================================

    nRiRoH_h = nPhiDifference * nPhiSum;
    reaction.dwbaGrid.allocateRioEx(nRiRoH_h, reaction);
    rIHVector.assign(nRiRoH_h + 1, 0.0);
    rOHVector.assign(nRiRoH_h + 1, 0.0);
    // Store riRoHCount into reaction.gridData.nRiRoInterp (reused for H-grid count separately)
    // We track nRiRoH_h as the H-search grid size; nRiRoInterp is assigned after
    // the 50x H-grid is found.

    //  WORK AREAS FOR STAGE 3

    // set vector sub-array pointers (declared nullptr at top, sized above)
    // shifting a base pointer by -1 gives a 1-based view: ptr[1] = base[0]
    // Match: P_x[1] = V_x[0] ⇒ P_x = V_x.data() - 1
    vPolysPointer = vPolysVector.data();  // 0-based: vPolysPointer[k-1]=vPolysVector[k-1]
    vWtsPointer  = vWorkVector.data();                       // 0-based: vWtsPointer[k-1]=vWorkVector[k-1]
    vResdPointer = vWorkVector.data() + nPhiSum;               // 0-based: vResdPointer[k-1]=vWorkVector[nPhiSum+k-1]
    vAmatPointer = vWorkVector.data() + 4*nPhiSum + 3;         // 0-based: vAmatPointer[k-1]=vWorkVector[4*nPhiSum+3+k-1]
    vDifptPointer = vDifPointsVector.data() - 1;  // vDifptPointer[vIndex]=vDifPointsVector[vIndex-1] for vIndex>=1
    vDifWtPointer = vDifWeightsVector.data() - 1;  // vDifWtPointer[vIndex]=vDifWeightsVector[vIndex-1] for vIndex>=1
    // vMin/vMid/vMax live at reaction.gridData.smhwkPointer offsets (1-based);
    // the sum-grid points at reaction.gridData.smhptsPointer (1-based)
    // set 0-based pointers from vectors (accessed [riRoIndex-1])
    rihPointer   = rIHVector.data();
    rohPointer   = rOHVector.data();

    // =========================================================================
    //  STAGE 1)  END POINTS OF EACH V-RANGE.
    // =========================================================================
    //      vMin = 1
    //      bumpCount = 0
    //      rPMax = 0
    //      rTMax = 0
    vMax   = 1.0;
    vMin   = 1.0;
    bumpCount = 0;
    rPMax  = 0.0;
    rTMax  = 0.0;
    vStepSize     = 1.0 / 250.0;   // LOOKST fixed 250


    for (uIndex = 1; uIndex <= nPhiSum; uIndex++) {
        double vLength;
        U    = reaction.gridData.smhptsPointer[uIndex];
        vLength = 2.0 * U;
        if (U >= 1.0) {

        //         vVal=VMAX
        //         syne= +.5*vLength
        vVal = vMax;
        syne = +0.5 * vLength;

        time1 = second();

        while (true) {
        vVal = std::fmin(1.0, vVal + 3.0 * vStepSize);

        while (true) {
        if (vVal <= 0.5 * vStepSize) break;   // was goto L450
        //            rI=U+vVal*syne
        //            rO=U-vVal*syne
        //            uLimit = rvrLimit/std::max( 1.D-2, rI*rO )
        rI   = U + vVal * syne;
        rO   = U - vVal * syne;
        uLimit = rvrLimit / std::fmax(1.0e-2, rI * rO);

        bool fifoExceedsULimit = false;
        for (i = 1; i <= 2; i++) {
            if (reaction.boundState.evaluateFormFactor(fifo, 2, rI, rO, xs[i], sctmn1b(), 1, rP, rT, reaction) != 0) break;
            interpCount = interpCount + inAdditional;
            if (std::fabs(fifo) > uLimit) { fifoExceedsULimit = true; break; }
        }
        if (fifoExceedsULimit) break;   // was goto L450

        vVal = vVal - vStepSize;
        }  // end while (inner search; was L422 loop)

        vVal = std::fmin(1.0, vVal + vStepSize);

        //  MAKE SURE THIS DESIRED VALUE OF V DOES NOT VIOLATE ASYMPTOPIA
        //  BECAUSE  S1*t1 < 0  AND  S2*T2 < 0,  RP2 > RP1 AND
        //  RT2 > RT1 ALWAYS.
        //
        //         rI=U+vVal*syne
        //         rO=U-vVal*syne
        //         rTMax = std::max( rTMax, rT )
        //         rPMax = std::max( rPMax, rP )
        rI    = U + vVal * syne;
        rO    = U - vVal * syne;
        rT    = std::sqrt((s1 * rI) * (s1 * rI) + (t1 * rO) * (t1 * rO) + st1 * rI * rO * (1.0 - xQuadStep));
        rP    = std::sqrt((s2 * rI) * (s2 * rI) + (t2 * rO) * (t2 * rO) + st2 * rI * rO * (1.0 - xQuadStep));
        rTMax = std::fmax(rTMax, rT);
        rPMax = std::fmax(rPMax, rP);
        if (!(rT <= boundMxTarget && rP <= boundMxProj)) {
            //  CANNOT USE DESIRED VALUE
            //
            //         bumpCount = bumpCount + 1
            //         vVal = vVal - DV
            bumpCount = bumpCount + 1;
            vVal   = vVal - vStepSize;
        }


        if (syne > 0.0) vMax = vVal;
        if (syne < 0.0) {
            vMin = vVal;
            break;   // 2nd pass done
        }
        vVal = vMin;
        syne = -syne;
        }  // end while (outer search; was L420 loop)


        tIntrp = tIntrp + second() - time1;
        }

        if (gridPrintSwitch) {
            std::printf(" U VMIN VMAX =%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G\n",
                        U, vMin, vMax, rP, rT, fifo);
        }

        //  THESE ARE EXTREMA OF V=DIF FOR EACH POINT IN THE U=sum GRID.
        reaction.gridData.smhwkPointer[uIndex] = -vMin * vLength;
        reaction.gridData.smhwkPointer[2*nPhiSum + uIndex] =  vMax * vLength;
        //  CHECK FOR RANGES CONTRACTED TO ZERO BY ASYMPTOPIA LIMITATIONS
        //
        //         vVal = 2*U*(vMin+VMAX)
        {
            double vValCheck = 2.0 * U * (vMin + vMax);
            if (vValCheck < 0.1 * vStepSize * U) {
                std::printf("\n0**** ERROR:  THE DIF GRID RANGE FOR SUM =%8.2f"
                            " IS ESSENTIALLY ZERO DUE TO BOUND STATE ASYMPTOPIA "
                            "LIMITATIONS; RANGE =%13.3G\n"
                            " **** IN FUTURE RUNS ASYMPTOPIA MUST BE INCREASED.\n",
                            U, vValCheck);
            }
        }

    }   // end DO 489  uIndex=1,nPhiSum

    if (bumpCount != 0) {
        std::printf("\n0**** WARNING:  SEARCH FOR THE DIF GRIDS WAS STOPPED%4d"
                    " TIMES BY BOUND STATE ASYMPTOPIA.\n",
                    bumpCount);
        printAsymptopiaAdvice(rPMax, rTMax);
    }

    // =========================================================================
    // Stage 2) Find V-moments of the weight function for each U
    // =========================================================================
    if (gridPrintSwitch) std::printf("\n IU        U           VMIN        VMID        VMAX"
        "      V(MAX)      pvpMax         <V>      <V**2>     <V VAR>\n");

    for (uIndex = 1; uIndex <= nPhiSum; uIndex++) {
        U = reaction.gridData.smhptsPointer[uIndex];
        vMin = reaction.gridData.smhwkPointer[uIndex];
        vMax = reaction.gridData.smhwkPointer[2*nPhiSum + uIndex] ;
        vWtsPointer[uIndex - 1] = 1.0 / ((vMax - vMin) * (vMax - vMin)); // vWtsPointer 0-based, [uIndex-1]=vWorkVector[uIndex-1]
        vPtCount = 2 * nPhiDifference;
        vStepSize = (vMax - vMin) / (vPtCount + 1);
        vVal = vMin;
        pvpMax = 0; sum0 = 0; sum1 = 0; sum2 = 0;
        for (vIndex = 1; vIndex <= vPtCount; vIndex++) {
            vVal = vVal + vStepSize;
            rI = U + 0.5 * vVal;
            rO = U - 0.5 * vVal;
            temp = 0;
            bool skip539 = false;
            for (i = 1; i <= 2; i++) {
                if (reaction.boundState.evaluateFormFactor(fifo, 1, rI, rO, xs[i], sctcr1b(), 1, rP, rT, reaction) != 0)
                    { skip539 = true; break; }
                interpCount = interpCount + inAdditional;
                temp = std::fabs(fifo) + temp;
            }
            if (skip539) continue;
            sum0 = sum0 + temp;
            sum1 = sum1 + temp * vVal;
            sum2 = sum2 + temp * vVal * vVal;
            if (temp <= pvpMax) continue;
            pvpMax = temp;
            vAtPeak = vVal;
        }
        sum1 = sum1 / sum0;
        sum2 = sum2 / sum0;
        sum0 = std::sqrt(sum2 - sum1 * sum1);
        sum2 = std::sqrt(sum2);
        vMid = vAtPeak;
        temp = 0.3 * (vMax - vMin);
        vMid = std::min(std::max(vMid, vMin + temp), vMax - temp);
        if (gridPrintSwitch)
            std::printf("%4d%14.4G%14.4G%14.4G%14.4G%14.4G%14.4G%14.4G%14.4G%14.4G\n",
                        uIndex, U, vMin, vMid, vMax, vAtPeak, pvpMax, sum1, sum2, sum0);
        vMid = (vMid - vMin) / (vMax - vMin);
        reaction.gridData.smhwkPointer[nPhiSum + uIndex] = vMid;
    }

    // =========================================================================
    // Stage 3) Polynomial fits to V-ranges
    // =========================================================================
    vTermCount = vPolyDegree + 1;

    // lsqPol uses 1-based vector pointers (P_x[1]=V_x[0])
    // use class pointers directly (smhptsPointer[1]=smhpts[0], smhwkPointer[1]=smhwk[0])
    // P_VSMSQ dropped — sum out-param (per-poly sum of W*resid^2) was
    // written by lsqPol and never read; lsqPol no longer takes it.
    lsqPol(reaction.gridData.smhptsPointer, reaction.gridData.smhwkPointer, vWtsPointer,
           vResdPointer, nPhiSum, 3, vAmatPointer,
           vPolysPointer, vTermCount);

    if (isInfoPrint) {
        std::printf("\n\n%2d-DEGREE POLYNOMIAL FITS TO DIF RANGES:\n"
                    "\n      V LO           V MID          V HI\n", vPolyDegree);
        for (n = 1; n <= vTermCount; n++)
            std::printf("%15.5G%15.5G%15.5G\n",
                        vPolysPointer[n - 1], vPolysPointer[n + vTermCount - 1],
                        vPolysPointer[n + 2 * vTermCount - 1]);  // vPolysPointer 0-based, [k-1]=vPolysVector[k-1]
    }
    if (gridPrintSwitch) std::printf("\n       U             V LO     RESID"
        "               V MID     RESID               V HI     RESID\n\n");

    vRangeMax = 0; vRangeMin = Constants::bigNum; vAbsMax = 0;
    for (uIndex = 1; uIndex <= nPhiSum; uIndex++) {
        if (!skipPolynomials) {
            for (i = 1; i <= 3; i++)
                // vResdPointer is 0-indexed
                reaction.gridData.smhwkPointer[(i-1)*nPhiSum + uIndex] += vResdPointer[(i-1)*nPhiSum + uIndex - 1];
        }
        vMin = reaction.gridData.smhwkPointer[uIndex];
        vMax = reaction.gridData.smhwkPointer[2*nPhiSum + uIndex];
        vMid = (vMax - vMin) * reaction.gridData.smhwkPointer[nPhiSum + uIndex] + vMin;
        reaction.gridData.smhwkPointer[nPhiSum + uIndex] = vMid;
        if (gridPrintSwitch) {
            std::printf("%14.4G", reaction.gridData.smhptsPointer[uIndex]);
            for (i = 1; i <= 3; i++)
                std::printf("   %14.4G%12.2G",
                    reaction.gridData.smhwkPointer[(i-1)*nPhiSum + uIndex], vResdPointer[(i-1)*nPhiSum + uIndex - 1]);  // R4+R6
            std::printf("\n");
        }
        vRangeMax = std::max(vRangeMax, vMax - vMin);
        vRangeMin = std::min(vRangeMin, vMax - vMin);
        vAbsMax = std::max({vAbsMax, vMax, -vMin});
        temp = 0.2 * (vMax - vMin);
        if (vMin > vMid - temp || vMax < vMid + temp)
            std::printf(" **** INVALID VMIN, VMID, VMAX:%15.5G%15.5G%15.5G\n", vMin, vMid, vMax);
    }

    // =========================================================================
    // Stage 4) Generate the rI, rO points
    // =========================================================================

    for (uIndex = 1; uIndex <= nPhiSum; uIndex++) {
        U = reaction.gridData.smhptsPointer[uIndex];
        vMin = reaction.gridData.smhwkPointer[uIndex];
        vMax = reaction.gridData.smhwkPointer[2*nPhiSum + uIndex] ;
        vMid = reaction.gridData.smhwkPointer[nPhiSum + uIndex]   ;

        if (vMin < -2*U || vMax > 2*U) {
            std::printf("0*** ILLEGAL VMIN OR VMAX:%3d", uIndex);
            printValuesSetToMax(U, vMin, vMax);
        }
        vMin = std::max(vMin, -2.0*U);
        vMax = std::min(vMax, 2.0*U);
        temp = 0.3*(vMax - vMin);
        vMid = std::min(std::max(vMid, vMin + temp), vMax - temp);

        syne = 1;
        if (vMid > 0.5*(vMax + vMin)) {
            vMid = -vMid; temp = vMax; vMax = -vMin; vMin = -temp; syne = -1;
        }

        // use 1-based vDifptPointer/vDifWtPointer for cubMap (cubMap uses ptr[1..nPhiDifference])
        cubMap(mapDif, vMin, vMid, vMax, gammaDif, vDifptPointer, vDifWtPointer, nPhiDifference);

        for (vIndex = 1; vIndex <= nPhiDifference; vIndex++) {
            riRoIndex = (syne < 0) ? (nPhiDifference - vIndex) : (vIndex - 1);
            riRoIndex = nPhiSum * riRoIndex + uIndex;
            vVal = vDifptPointer[vIndex];  // vDifptPointer[vIndex]=vDifPointsVector[vIndex-1]
            rI = U + 0.5*vVal*syne;
            rO = U - 0.5*vVal*syne;
            rihPointer[riRoIndex - 1] = rI;
            rohPointer[riRoIndex - 1] = rO;
            rP = std::sqrt(1 + (s1*rI + t1*rO)*(s1*rI + t1*rO));
            rT = std::sqrt(1 + (s2*rI + t2*rO)*(s2*rI + t2*rO));
            reaction.gridData.rioExPointer[riRoIndex - 1] = std::exp(alphaP*rP + alphaT*rT);  // (0-based, matches rih/rohPointer[riRoIndex-1])
        }
    }


    } // end cubMap block

    // =========================================================================
    // 5) Generate points for the interpolated integral grid
    // =========================================================================
    riRoHCount = nRiRoH_h;  // H-grid size used for phi-grid passes
    area = 0;
    // Stage 5 pointer aliases: 1-based (uIndex=1 → first element)
    smivl5Pointer = reaction.dwbaGrid.smivl.data();  // 0-based: smivl5Pointer[uIndex-1]=smivl[uIndex-1]
    sumpt5Pointer = reaction.dwbaGrid.smipts.data(); // 0-based: sumpt5Pointer[uIndex-1]=smipts[uIndex-1]
    // iRi/IRO/IWIO class-owned; use reaction.gridData.riPointer/roPointer/wioPointer (0-based float*)
    riPointer  = reaction.gridData.riPointer;
    roPointer  = reaction.gridData.roPointer;
    wioPointer = reaction.gridData.wioPointer;

    for (uIndex = 1; uIndex <= nInterpPoints; uIndex++) {
        wow = smivl5Pointer[uIndex - 1];
        U = sumpt5Pointer[uIndex - 1];

        // Evaluate polynomial-determined ranges
        // ITEMP is now index into vPolysVector (0-based); LVPOLS=-1 → start at -1
        // Bug fix: only read from smhwkPointer when skipPolynomials=true (uIndex<=nPhiSum guaranteed);
        // when skipPolynomials=false, polynomial evaluation overwrites vRange anyway.
        int vPolyOffset = -1;
        for (i = 1; i <= 3; i++) {
            if (skipPolynomials) {
                vRange[i] = reaction.gridData.smhwkPointer[(i-1)*nPhiSum + uIndex];  // uIndex<=nPhiSum here
                continue;
            }
            vPolyOffset = vPolyOffset + vPolyDegree + 1;
            vRange[i] = vPolysVector[vPolyOffset];
            for (n = 1; n <= vPolyDegree; n++)
                vRange[i] = vPolysVector[vPolyOffset - n] + U * vRange[i];
        }
        vMin = vRange[1]; vMid = vRange[2]; vMax = vRange[3];
        if (!skipPolynomials) vMid = (vMax - vMin) * vMid + vMin;

        if (vMin < -2*U || vMax > 2*U) {
            std::printf("0*** ILLEGAL INTERPOLATED VMIN OR VMAX:%4d", uIndex);
            printValuesSetToMax(U, vMin, vMax);
        }

        area = area + wow * (vMax - vMin);
        vMin = std::max(vMin, -2.0*U);
        vMax = std::min(vMax, 2.0*U);
        temp = 0.3*(vMax - vMin);
        vMid = std::min(std::max(vMid, vMin + temp), vMax - temp);
        syne = 1;
        if (vMid > 0.5*(vMax + vMin)) {
            vMid = -vMid; temp = vMax; vMax = -vMin; vMin = -temp; syne = -1;
        }

        // use 1-based vDifptPointer/vDifWtPointer for cubMap (1-based, ptr[1..nPhiDifference])
        cubMap(mapDif, vMin, vMid, vMax, gammaDif, vDifptPointer, vDifWtPointer, nPhiDifference);

        for (vIndex = 1; vIndex <= nPhiDifference; vIndex++) {
            riRoIndex = (syne < 0) ? (nPhiDifference - vIndex) : (vIndex - 1);
            riRoIndex = nInterpPoints * riRoIndex + uIndex;
            vVal = vDifptPointer[vIndex];  // vDifptPointer[vIndex]=vDifPointsVector[vIndex-1]
            rI = U + 0.5*vVal*syne;
            rO = U - 0.5*vVal*syne;
            riMax = std::max(riMax, rI);
            roMax = std::max(roMax, rO);
            if (riRoIndex > nRiRoInterp_g) {
                std::printf(" *** IPLUNK OVERFLOW: %d > NRIROI=%d IU=%d IV=%d n_interp_pts=%d NPDIF=%d\n",
                    riRoIndex, nRiRoInterp_g, uIndex, vIndex, nInterpPoints, nPhiDifference);
                break;
            }
            riPointer[riRoIndex] = (float)rI;
            roPointer[riRoIndex] = (float)rO;

            rP = std::sqrt(1 + (s1*rI + t1*rO)*(s1*rI + t1*rO));
            rT = std::sqrt(1 + (s2*rI + t2*rO)*(s2*rI + t2*rO));
            temp = std::exp(-(alphaP*rP + alphaT*rT));
            wioPointer[riRoIndex] = (float)(jacobian * rI * rO * wow * vDifWtPointer[vIndex] * temp);  // R4+R6
        }
    }


    if (isInfoPrint)
        std::printf(" MAXIMUM RI AND RO:%15.5G%15.5G\n", riMax, roMax);

    // Expand wavefunction work area if needed
    reaction.distortedWave.channel[1].nStp2s = (int)(riMax / reaction.distortedWave.channel[1].rStart + 3.5);
    reaction.distortedWave.channel[2].nStp2s = (int)(roMax / reaction.distortedWave.channel[2].rStart + 3.5);
    i = std::max(reaction.distortedWave.channel[1].nStp2s, reaction.distortedWave.channel[2].nStp2s) + 6;
    if (i > nWaveF) {
        nWaveF = i;
        reaction.distortedWave.scatteringSolver.reallocate(i);
    }

    // =========================================================================
    // Phi grid setup (two passes)
    // =========================================================================

    logicVector.assign(riRoHCount + 2, 0);  // 1-based [1..riRoHCount], [0] unused

    iiMax = 250 + 1;   // LOOKST fixed 250
    xQuadStep = 2.0 / (250.0 * 250.0);
    aitkenCount = 2;
    formType = 2;
    isFirstPass = true;
    x0Min = 1; x0Av = 0;

    vPhiPointsVector.assign(nPhiPoints + 1, 0.0);
    vPhiWeightsVector.assign(nPhiPoints + 1, 0.0);
    vPhiPtPointer = vPhiPointsVector.data() - 1;  // 1-based
    vPhiWtPointer = vPhiWeightsVector.data() - 1;  // 1-based
    phiMid = std::max(phiMid, 0.10);
    phiMid = std::min(phiMid, 0.90);
    cubMap(mapPhi, 0.0, phiMid, 1.0, gammaPhi,
           vPhiPtPointer, vPhiWtPointer, nPhiPoints);

    // Two-pass loop: pass 1 finds x0, pass 2 stores phi grid
    while (true) {
    phiCount = 0; rTMax = 0; rPMax = 0; bumpCount = 0;
    // vPhiPtPointer/vPhiWtPointer already set from vPhiPointsVector/vPhiWeightsVector (no refresh needed)
    // LPHIPT/LPHIWT removed (use vPhiPtPointer/vPhiWtPointer directly)
    // rihPointer/rohPointer already point to rIHVector/rOHVector (no refresh needed)
    rihPointer = rIHVector.data();
    rohPointer = rOHVector.data();
    // logicPointer points to local logicVector (0-based: logicPointer[riRoIndex]=logicVector[riRoIndex])
    logicPointer = logicVector.data();

    for (riRoIndex = 1; riRoIndex <= riRoHCount; riRoIndex++) {
        bool earlyX0 = false;  // pass-1 inner-loop early exit
        rI = rihPointer[riRoIndex - 1];
        rO = rohPointer[riRoIndex - 1];
        uLimit = rvrLimit / (std::max(1.0, rI) * std::max(1.0, rO));

        if (!isFirstPass) {
            // Pass 2: set up end point
            xEnd = logicPointer[riRoIndex];
            xEnd = std::max(xEnd, 2);
            x0 = 1.0 - xQuadStep * (xEnd-1) * (xEnd-1);
            phi0 = std::acos(x0);
        }

        if (isFirstPass) time1 = second();

        wow = 0;
        for (ii = 1; ii <= iiMax; ii++) {
            x = 1.0 - xQuadStep * (ii-1) * (ii-1);
            if (!isFirstPass) {
                phi = phi0 * vPhiPtPointer[ii];
                dPhi = phi0 * vPhiWtPointer[ii] * std::sin(phi);
                x = std::cos(phi);
            }

            // Calculate bound-state product
            bumpCount++;
            {
                int bsReturn = reaction.boundState.evaluateFormFactor(fifo, formType, rI, rO, x, sctmn1b(), aitkenCount, rP, rT, reaction);
                bumpCount--;
                if (bsReturn == 0) interpCount += inAdditional;
            }
            rTMax = std::max(rTMax, rT);
            rPMax = std::max(rPMax, rP);
            if (!isFirstPass) {
                // Pass 2: store phi angles and weighted BS product
                phiCount++;
                temp = (t1*rO + s1*rI*x) / rT;
                phiTPointer[phiCount] = (float)(phiSign * std::acos(temp));
                temp = (t2*rO + s2*rI*x) / rP;
                phiPPointer[phiCount] = (float)(phiSign * std::acos(temp));
                phiPointer[phiCount]  = (float)phi;
                trapWeightPointer[phiCount]  = (float)(dPhi * fifo);
                continue;
            }

            // Pass 1: finding x0
            fifo = std::fabs(fifo);
            if (ii != 1 && fifo < uLimit && wow < uLimit) { earlyX0 = true; break; }
            wow = fifo;
        } // end ii loop

        if (!isFirstPass) {
            // Pass 2: dump if needed
            if (verbosity >= 9) {
                n = phiCount - nPhiPoints + 1;
                std::printf(" %7.2f%7.2f%7.4f", rI, rO, phi0);
                for (i = n; i <= phiCount; i++) std::printf("%9.1E", (double)trapWeightPointer[i]);
                std::printf("\n");
            }
            continue;
        }

        if (!earlyX0) ii = 250 + 2;   // LOOKST fixed 250
        tIntrp += second() - time1;
        xEnd = std::min(ii - 1 + nPhiAdditional, iiMax);
        logicPointer[riRoIndex] = xEnd;
        if (xEnd != 1) phiCount += nPhiPoints;
        x0 = 1.0 - xQuadStep * (xEnd-1) * (xEnd-1);
        x0Min = std::min(x0Min, x0);
        x0Av = x0Av + x0;

    } // end riRoIndex loop

    if (isFirstPass) {
        // End of pass 1 — allocate and setup for pass 2
        if (bumpCount != 0) {
            std::printf(" **** WARNING:  IN THE PHI-GRID SEARCHES%8d"
                        " BOUND STATE WAVEFUNCTIONS WERE NEEDED BEYOND THE"
                        " BOUND STATE ASYMPTOPIA.\n", bumpCount);
            printAsymptopiaAdvice(rPMax, rTMax);
        }

        // debug removed
        phiBlockCount = (phiCount + facFr4 - 1) / facFr4;
        reaction.dwbaGrid.allocatePhiArrays(phiBlockCount, reaction);
        maxCount = phiCount;
        x0Av = x0Av / riRoHCount;

        if (verbosity >= 9) std::printf("     RI     RO               PHI V PHI  AT EACH X\n\n");

        // function-scope float* aliases (use reaction.gridData.P_ set by allocatePhiArrays)
        phiTPointer = reaction.gridData.phiTPointer;
        phiPPointer = reaction.gridData.phiPPointer;
        phiPointer  = reaction.gridData.phiPointer;
        trapWeightPointer  = reaction.gridData.trapWeightPointer;
        iiMax = nPhiPoints;
        aitkenCount = 4;
        formType = 1;
        inAdditional = 0;
        isFirstPass = false;
        continue;  // start pass 2
    }
    break;  // pass 2 done — exit while(true)
    }  // end while (two-pass loop)

    // =========================================================================
    // End of pass 2 — check asymptopia overruns
    // =========================================================================
    // debug removed
    if (bumpCount != 0) {
        std::printf(" **** ERROR:%8d BOUND STATE WAVEFUNCTIONS WERE NEEDED"
                    " BEYOND THE BOUND STATE ASYMPTOPIA.\n", bumpCount);
        printAsymptopiaAdvice(rPMax, rTMax);
    }


    // =========================================================================
    // Print summary
    // =========================================================================
    if (printSwitch) {

    std::printf("\n\n0             SUMMARY OF THE INTEGRATION GRIDS\n"
        "0         NUM. PTS.   MAP TYPE   GAMMA    MINIMUM   \"MID. PT.\"   MAXIMUM\n"
        " (RI+RO)/2:%5d%11d%11.2f%10.2f%11.2f%11.2f"
        "     (H'S FOUND AT%4d (RI+RO)/2 VALUES)\n"
        "0                                        MIN. WIDTH   MAX. WIDTH\n"
        " RI-RO:%9d%11d%11.2f%11.3f%13.3f\n"
        "0                                          NPHIADD    MIN. START   AV. START     phiMid\n"
        " COS(PHI):%6d%11d%11.2f%9d%16.5f%12.5f%12.4f\n",
        nInterpPoints, mapSum, gammaSum, sumMin, sumMid, sumMax, nPhiSum,
        nPhiDifference, mapDif, gammaDif, vRangeMin, vRangeMax,
        nPhiPoints, mapPhi, gammaPhi, nPhiAdditional, x0Min, x0Av, phiMid);


    std::printf("0DWCUTOFF =%12.2E\n"
        " MAXIMUM R'S USED FOR BOUND STATES:  PROJECTILE =%5.1f     TARGET =%5.1f\n"
        " MAXIMUM R'S USED FOR SCATTERING STATES:  INCOMING =%5.1f     OUTGOING =%5.1f\n",
        dwCutoff, rPMax, rTMax, riMax, roMax);

    std::printf("0THERE ARE%6d POINTS IN THE (RI, RO) GRID ON WHICH THE H'S ARE FOUND,"
        "  AND%8d POINTS IN THE COMPLETE 3-D GRID.\n"
        " THERE ARE%6d POINTS IN THE (RI, RO) INTERPOLATED INTEGRATION GRID.\n"
        " AREA OF  SUM X DIF  INTEGRATION REGION =%9.3f FM**2.\n\n",
        riRoHCount, maxCount, nRiRoInterp_g, area);

    }  // end if (printSwitch)


    returnCode = 1;


    // Timing
    if (!isInfoPrint) return;
    wow = second() - tStart;
    temp = wow - tIntrp;
    U = 1000.0 * tIntrp / interpCount;
    std::printf("\n        END OF GRID INITIALIZATION\n"
        "\nROUGH INTERPOLATION TIME FOR GRID SEARCHING:%54s%7.3f SECONDS\n"
        " OTHER TIME (INCLUDES FINAL INTERPOLATIONS):%54s%7.3f SECONDS\n"
        " TOTAL TIME:%54s%7.3f SECONDS\n"
        "\n%8d ROUGH INTERPOLATIONS MADE; AVERAGE TIME =%54s%7.3f MILLISECONDS\n\n",
        "", tIntrp, "", temp, "", wow, interpCount, "", U);
    return;


} // end DWBAGrid::gridSet

// ---------------------------------------------------------------------------
// Backward-compat free function wrapper — calls the class method
// ---------------------------------------------------------------------------
