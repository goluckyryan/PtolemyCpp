// integration_grid.cpp — INGRST: initializes the grid for the inelastic-scattering
// radial integral (Gauss points/weights + nuclear/Coulomb potential form factors).

#include "ptolemy_types.h"
#include "linkule.h"
#include "math/numeric_utils.h"
#include "math/spline.h"
#include "print_utils.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include "Reaction.h"
#include "Constants.h"

void DWBAGrid::inelasticGridSet(int& returnCode, Reaction& reaction) {
    //
    // 1/77 - FIRST VERSION - S.P.
    // 11/22/77 - NEW SCTASYOPIA - S.P.
    // 12/4/77 - LINKULES - S.P.
    // 5/7/78 - STUFF TO ALLOW 2 COUPLED STATES - S.P.
    // 5/23/78 - FIX BUG IN 11/22/77 - ALWAYS DEFINE VWORK - S.P.
    //           BIGGER waveReal, waveImag ARRAYS FOR NEW WAVELJ - S.P.
    // 3/12/79 - CLEAN UP DEF OF sumMax, nSteps; DEFINE MORE FOR LINKULES
    // 12/20/79 - INELASTIC NAMES INPLACE OF TRANSFER NAMES - S.P.
    //            START OF FULL COUPLED CHANNELS CHANGES - S.P.
    //  1/2/80 PUT IN FULL COUPLED CHANNELS FOR ROTATIONAL MODEL
    //         UP TO 2 ORDER M.R.B.
    // 3/18/80 - FIX BUG IN FREEING ITEMP - S.P.
    // 4/27/80 - GET NPTMIN, FREE UNNEEDED ARRAYS FOR C.C. - S.P.
    // 4/28/80 - REASSIGN LRPTS AFTER GIVALL - S.P.
    // 6/24/80 - SAVE rC IN rcuEff - S.P.
    // 7/9/80 - FIX BUG IN FREEING COULH, DON'T OVERFLOW AREA
    // 8/25/80 - FIX TYPING ERROR AT STMNT 820 - S.P.
    // 12/9/80 - MAKPOT NOLONGER COMPUTES COULOMB; DON'T
    //           DEFINE COULOMB COUPLING FOR BETAC=0 - S.P.
    // 12/30/80 - FIX ERRORS IN 2ND ORDER NUCLEAR - S.P.
    // 1/11/81 - ALWAYS 1/R**(lx+1); NO MORE COUPL2 - S.P.
    // 7/13/81 - FIX BUG FROM 1/11/81 (JUMPSZ) - S.P.
    // 7/16/81 - SAVE LNKBLK TO ALLOW FORM-factor CHANGES - S.P.
    // 7/29/81 - FIRST CRUDE STUFF FOR DEFORMED POTENTIALS - S.P.
    // 8/28/81 - MORE FOR DEFORMED POTENTIALS - S.P.
    // 9/21/81 - DON'T LOOK WHEN THERE ARE NO POTENTIALS - S.P.
    // 10/16/81 - DONT' CONDISDER LXTOT FOR MUTUAL - S.P.
    // 12/14/18 - LET GETSCT GET WAVEFUNC ALSO FOR C.C. - S.P.
    // 4/14/82 - FORMASIS ALSO WORKS FOR DWBA - S.P.
    // 6/10/83 - FIX BUGS AT STMNT 235, 408 - S.P.
    // 7/24/84 - HOTFUDGE - S.P.
    // 10/16/84 - USE ABS(scatAsy) - S.P.
    // 12/19/84 - USE parameterPrint; DEFORMED SI - S.P.
    // 4/29/86 - FIX BUG FOR IMAG+DEFORMED+HOTFUDGE - S.P.
    // 11/16/01 - initialize IFIRST, ISECND at very start
    //

    double& undefValue   = reaction.internalState.undefValue;

    constexpr int mapSum = 2;
    // notDefSentinel/waveChannel/printLevel/nPhiSum/zProj/zTarget/L_int single-use aliases

    double& J_dbl   = reaction.angMom.J;       // reaction.angMom.J


    double& stepSize  = reaction.integrationGrid.stepSize;
    double& sumMax  = reaction.integrationGrid.sumMax;
    double& sumMid  = reaction.integrationGrid.sumMid;
    double& sumMin  = reaction.integrationGrid.sumMin;
    double& gammaSum  = reaction.rxn.gammaSum;
    double& scatAsy  = reaction.integrationGrid.scatAsy;
    double& r0Mass  = reaction.internalState.r0Mass;
    double& R       = reaction.integrationGrid.R;
    double& rC      = reaction.opticalPotentialParams.rC;
    double& rI      = reaction.opticalPotentialParams.rI;
    double& rSo     = reaction.opticalPotentialParams.rSo;
    double& rSoi    = reaction.opticalPotentialParams.rSoi;
    double& rSi     = reaction.opticalPotentialParams.rSi;
    double& tvReal     = reaction.distortedWave.scatteringSolver.potentialWork.tvReal;
    double& tvImag     = reaction.distortedWave.scatteringSolver.potentialWork.tvImag;
    double& taReal     = reaction.distortedWave.scatteringSolver.potentialWork.taReal;
    double& taImag     = reaction.distortedWave.scatteringSolver.potentialWork.taImag;
    double& E       = reaction.energies.E;
    double& vSi     = reaction.opticalPotentialParams.vSi;
    double& vSo     = reaction.opticalPotentialParams.vSo;
    double& vSoi    = reaction.opticalPotentialParams.vSoi;
    double& rcProj     = reaction.masses.rcProj;
    double& rcTarget     = reaction.masses.rcTarget;
    double& rc0Proj    = reaction.masses.rc0Proj;
    double& rc0Target    = reaction.masses.rc0Target;
    double& aM      = reaction.masses.aM;
    double& massProj     = reaction.masses.massProj;
    double& massTgt     = reaction.masses.massTgt;

    double& redMi   = reaction.kin.redMi;
    double& rcuEff  = reaction.gridData.jacobian;

    int&    lxMin   = reaction.inelastic.lxMin;
    // gridData.IWIO field stays for now (DWBAGrid sentinel).
    int&    numPoint   = reaction.gridData.nRiRoInterp;

    int*    zArray     = reaction.charges.zArray;      // 1-based

    // are already 0 sentinels set by allocateRiRoWio. CC repurposing of the
    // handles is gone (CC unreachable).


    int&    uniqueLinkuleId  = reaction.linkuleData.uniqueLinkuleId;

    // --- Local variables ---
    int    verbosity;

    int    i, ii;
    int    stepCount, gridStepsSave, nSteps;
    int    savedNumPoint;

    // Pointers into ALLOC
    double* rptsPointer = nullptr;  // direct pointer to RPTS pool
    double* rwtsPointer = nullptr;  // direct pointer to RWTS pool
    // LCOULH retired (LCOULH_loc local at the one use site is the only live one).
    std::vector<double> cubicSplineVector;  // allocated at L600 only if needed
    double *splineRPointer = nullptr, *splineBPointer = nullptr, *splineCPointer = nullptr, *splineDPointer = nullptr;


    double temp1, temp2, temp5, temp6;
    double eta;

    // Loop variables
    int    stepIndex;
    double rValue, wt, x, xx;
    double term;
    double coulombCoupling;

    // CC Coulomb-loop locals (LABEL/IORDER/IMULT/IPOINT/IMT/IPOW/TMPFAC/I1)

    // The Fortran ASSIGN-goto IGOTO_real/IGOTO_imag flags were collapsed

    // =========================================================================
    // Initialise at start (11/16/01 bug fix)
    // =========================================================================
    returnCode   = 0;

    verbosity = ((reaction.flags.printLevel) % (10));

    //
    // FIND WHERE lMin SCATTERING WAVE EXCEEDS 1E-15
    //
    nSteps   = reaction.distortedWave.channel[1].nGridSteps;
    {
    // read from boundState.data.sctmnArr directly. Preserves the
    // legacy off-by-one (writer fills slot[1..nSteps+1] via 1-based access;
    // reader here uses 0-based sctmnPointer so sctmnPointer[i] = sctmnArr[i] = slot[i+1]).
    const double* sctmnPointer = reaction.boundState.data.sctmnArr.data();
    for (i = 1; i <= nSteps; i++) {
        if (std::fabs(sctmnPointer[i]) > 1.0e-15) break;
    }
    } // sctmnPointer scope
    if (i > nSteps) {
        std::printf("0**** ERROR:  U(LMIN) < 1E-15 EVERYWHERE.\n");
        return;
    }

    if (sumMin == undefValue) sumMin = reaction.distortedWave.channel[1].rStart * (i - 1);
    if (sumMax == undefValue) sumMax = std::fabs(scatAsy);
    if (sumMid == undefValue) sumMid = 0.5 * (sumMin + sumMax);

    //
    // NUMBER OF GAUSS POINTS IS DETERMINED FROM SUMPTSPER.
    numPoint  = (int)((sumMax - sumMin) * (reaction.integrationGrid.sumDensity * (reaction.kin.akIn + reaction.kin.akOut) / (4.0 * Constants::PI)));
    numPoint  = std::max(numPoint, reaction.gridData.nPhiSum);
    savedNumPoint = numPoint;


    //
    // MAY HAVE TO REDEFINE WAVELJ WORK AREAS
    //
    i = (int)(sumMax / reaction.distortedWave.channel[1].rStart + 0.5);
    reaction.distortedWave.channel[1].nStp2s = std::max(i, reaction.distortedWave.channel[1].nStp2s);
    i = (int)(sumMax / reaction.distortedWave.channel[2].rStart + 0.5);
    reaction.distortedWave.channel[2].nStp2s = std::max(i, reaction.distortedWave.channel[2].nStp2s);
    i = std::max(reaction.distortedWave.channel[1].nStp2s, reaction.distortedWave.channel[2].nStp2s);
    if (i > std::max(reaction.distortedWave.channel[1].nGridSteps, reaction.distortedWave.channel[2].nGridSteps)) {
        reaction.distortedWave.scatteringSolver.reallocate(i + 6);
    }

    //
    // THE ARRAYS WILL BE:
    //   IRPTS, IRWTS - GAUSS POINTS AND WEIGHTS
    //   IRPTS4 - REAL*4 GAUSS POINTS FOR WAVELJ INTERPOLATOR
    //   INUCH - COMPLEX NUCLEAR POTENTIAL PART OF H WITH R'S AND GAUSS
    //           WEIGHT BUT NO BETA.
    //   ICOULH - COULOMB PART OF H WITH rC AND GAUSS WEIGHT BUT NO betas
    //
    // CCSW (CC) NBASCP form-factor enumeration + NNUCFF/NCOUFF nuch/coulh
    //
    // Reuse smiptsPointer/smivlPointer 1-based caches that downstream readers already use.
    reaction.gridData.rptsArr.assign(numPoint, 0.0);
    reaction.gridData.rwtsArr.assign(numPoint, 0.0);
    reaction.gridData.smiptsPointer = reaction.gridData.rptsArr.data() - 1;
    reaction.gridData.smivlPointer  = reaction.gridData.rwtsArr.data() - 1;
    reaction.gridData.rpts4Arr.assign(numPoint, 0.0f);
    reaction.gridData.rpts4Pointer = reaction.gridData.rpts4Arr.data() - 1;  // 1-based
    reaction.gridData.nuclearHArr.assign(2 * numPoint, 0.0);
    reaction.gridData.coulombHArr.assign(numPoint, 0.0);

    //
    // NOW SETUP THE POTENTIALS * GAUSS WEIGHTS * R'S
    //
    // FIRST WE GET THE POTENTIALS ON THE ARRAY ( 0, sumMax )
    // THEN WE DEFINE CUBICS AND USE THEM TO GET THE DERIVATIVES.
    //
    // TO GET THE POTENTIALS WE SETUP SOME WORK AREAS INPLACE OF
    // THE CHANNEL 1 ARRAYS.
    //

    reaction.internalState.waveChannel = 1;
    gridStepsSave  = reaction.distortedWave.channel[1].nGridSteps;
    {
        double* rlvSave  = reaction.distortedWave.channel[1].rlvPointer;
        double* imvSave  = reaction.distortedWave.channel[1].imvPointer;
        double* centSave = reaction.distortedWave.channel[1].centPointer;
        auto restoreChannelPotentials = [&]() {
            reaction.distortedWave.channel[1].rlvPointer  = rlvSave;
            reaction.distortedWave.channel[1].imvPointer  = imvSave;
            reaction.distortedWave.channel[1].centPointer = centSave;
            reaction.distortedWave.channel[1].nGridSteps = gridStepsSave;
        };
        // Error bailout shared by the two returnCode<0 checks below: clear the
        // error, restore the cached channel pointers, signal the caller to return.
        auto bailIfError = [&]() -> bool {
            if (returnCode < 0) {
                returnCode = 0;
                restoreChannelPotentials();
                return true;
            }
            return false;
        };

        //
        // WE USE THE STEPSIZE OF THE FIRST CHANNEL
        // WE GO SLIGHTLY FURTHER TO GIVE THE SPLINES TIME TO SETTLE DOWN
        //
        stepSize   = reaction.distortedWave.channel[1].rStart;
        E        = reaction.distortedWave.channel[1].Ecm;
        stepCount = (int)(sumMax / stepSize + 20.5);
        nSteps   = stepCount + 1;
        reaction.distortedWave.channel[1].nGridSteps = stepCount;
        std::vector<double> deformedRlv(nSteps, 0.0);
        std::vector<double> deformedImv(nSteps, 0.0);
        std::vector<double> deformedCentr(nSteps, 0.0);
        // allocateVWork sets vWork.size()=nSteps+1 and vWorkPointer = vWork.data()-1.
        reaction.distortedWave.scatteringSolver.allocateVWork(nSteps);
        reaction.distortedWave.channel[1].rlvPointer  = deformedRlv.data() - 1;
        reaction.distortedWave.channel[1].imvPointer  = deformedImv.data() - 1;
        reaction.distortedWave.channel[1].centPointer = deformedCentr.data() - 1;

        //
        // SETUP STUFF AS IF IN FIRST CHANNEL
        //
        massProj = reaction.masses.massesArr[0];
        massTgt = reaction.masses.massesArr[2];
        aM  = redMi;
        reaction.charges.zProj = zArray[1];
        reaction.charges.zTarget = zArray[3];

        //
        // INITIALIZE LINKULES IF NEEDED
        //
        // Linkule key slots scanned at grid setup: keys 1-5 plus key 13 (SIPOTENT).
        static const int linkuleSlots[7] = {0, 1, 2, 3, 4, 5, 13};  // [0] unused, 1-based
        for (ii = 1; ii <= 6; ii++) {
            i = linkuleSlots[ii];
            if (reaction.linkuleData.linkuleAddr[i][3] == 0) continue;

            //
            // GENERATE THE SPECIAL UNIQUE NAME
            //
            uniqueLinkuleId = uniqueLinkuleId + 1;
            // LNKNAM is a char8* overlay on linkuleAddr
            {
                // Build the 4-character LINKID as in Fortran:
                char linkId[5];
                linkId[0] = '*';
                linkId[1] = '3';
                linkId[2] = (char)('0' + ((uniqueLinkuleId / 10) % (10)));
                linkId[3] = (char)('0' + ((uniqueLinkuleId) % (10)));
                linkId[4] = '\0';

                //
                // MAKE INITIALIZING CALL
                //
                int linkulRet = 0;
                linkule(reaction.linkuleData.linkuleAddr[i][3],
                       reinterpret_cast<char8*>(&reaction.linkuleData.linkuleAddr[1][1])[i - 1],
                       &reaction.linkuleData.linkuleAddr[i][5], i, 1,
                       linkulRet, reaction.angMom.L, J_dbl, 0.0, stepSize, stepCount + 1,
                       nullptr, nullptr, linkId, reaction);
                returnCode = linkulRet;
                if (bailIfError()) return;
                reaction.linkuleData.linkuleAddr[i][4] = returnCode;
            }
        }

        //
        //
        // NOW COMPUTE THE POTENTIAL, channelIndex=3 IS CONVERT TO channelIndex=1
        // BY MAKPOT EXCEPT THE COULOMB IS NOT COMPUTED.
        //
        // line is `returnCode = 0;` so the 1-write is dead.
        J_dbl = reaction.internalState.notDefSentinel;
        reaction.makePotential(3, returnCode);
        if (bailIfError()) return;

        //
        // SUMMARIZE THE POTENTIALS TO BE USED
        //
        // FOLLOWING LIFTED FROM WAVSET
        //
        std::printf("-POTENTIALS WHOSE DERIVATIVES ARE THE EFFECTIVE INTERACTION OPERATOR:\n");

        eta   = reaction.kin.etaCh[1];

        //
        // IF POSSIBLE, RECOMPUTE R0'S FOR PRINTING
        //
        temp1 = 0.0;
        temp2 = 0.0;
        temp5 = 0.0;
        temp6 = 0.0;
        if (r0Mass != undefValue) {
            temp1 = R  / r0Mass;
            if (tvImag != 0.0)  temp2 = rI   / r0Mass;
            if (vSi != 0.0)  temp5 = rSi  / r0Mass;
            temp6 = rC / r0Mass;
        }

        std::printf("0POTENTIAL           COUPLING CONS.        RADIUS     DIFFUSENESS     RADIUS PARAMETER\n\n");
        parameterPrint(1,  "REAL CENTRAL       ", tvReal,  R,   taReal,  temp1, reaction);
        parameterPrint(2,  "VOLUME ABSORPTION   ", tvImag,  rI,  taImag,  temp2, reaction);
        parameterPrint(13, "SURFACE ABSORPTION   ", vSi,  rSi, 0.0,  temp5, reaction);
        parameterPrint(3,  "REAL SPIN-ORBIT   ", vSo,  rSo, 0.0,  0.0, reaction);
        parameterPrint(4,  "IMAG. SPIN-ORBIT   ", vSoi, rSoi, 0.0, 0.0, reaction);
        parameterPrint(5,  "COULOMB           ", eta,  rC,  0.0,  temp6, reaction);
        if ((rcProj != 0.0) && (rcTarget != 0.0)) {
            std::printf(" FOLDED COULOMB POTENTIALS - RCP = %7.4f  RCT = %7.4f  RC0P = %7.4f  RC0T = %7.4f\n",
                        rcProj, rcTarget, rc0Proj, rc0Target);
        }

        std::printf("0\n");

        //
        // PREPARE TO COMPUTE DERIVATIVE OF POTENTIAL
        // IF THE DERIVATIVE WAS DIRECTLY SUPPLIED, THEN WE
        // DON'T NEED TO COMPUTE IT UNLESS IF THERE IS 2ND ORDER.
        //

        // +4 extra for 1-based Fortran indexing in naturalCubicSpline (writes to array[nSteps])
        cubicSplineVector.assign(4 * nSteps + 4, 0.0);
        splineRPointer = cubicSplineVector.data();          // cubicSplineVector[0..nSteps-1]
        splineBPointer = splineRPointer + nSteps;           // cubicSplineVector[nSteps..2*nSteps-1]
        splineCPointer = splineBPointer + nSteps;           // cubicSplineVector[2*nSteps..3*nSteps-1]
        splineDPointer = splineCPointer + nSteps;           // cubicSplineVector[3*nSteps..4*nSteps-1]
        double* rlv1Pointer = reaction.distortedWave.channel[1].rlvPointer;
        double* imv1Pointer = reaction.distortedWave.channel[1].imvPointer;
        // that the inner loops expect (nuclearHPointer[2*i]/[2*i+1] for i=1..NUMPT; coulombHPointer[i]).
        double* nuclearHPointer  = reaction.gridData.nuclearHArr.data();  // 0-based: nuclearHPointer[2*i-2] = slot[2*i-2]
        double* coulombHPointer = reaction.gridData.coulombHArr.data();  // 0-based: coulombHPointer[i-1] = slot[i-1]
        //
        // Gauss-points DWBA branch — was wrapped in IFIRST+ISECND!=0 +
        // (rptsArr, rwtsArr) via the 1-based smiptsPointer/smivlPointer caches.
        rptsPointer = reaction.gridData.smiptsPointer;
        rwtsPointer = reaction.gridData.smivlPointer;
        cubMap(mapSum, sumMin, sumMid, sumMax, gammaSum,
               rptsPointer, rwtsPointer, numPoint);

        //
        // FIRST DO THE REAL AND COULOMB PARTS — GET THE SPLINES,
        // GENERATE INPUT GRID FOR THE SPLINES.
        //
        rValue = 0.0;
        for (i = 0; i < nSteps; i++) {
            splineRPointer[i] = rValue;
            rValue = rValue + stepSize;
        }

        naturalCubicSpline(nSteps, splineRPointer, rlv1Pointer + 1,
               splineBPointer, splineCPointer, splineDPointer);  // cubicSplineVector sub-arrays


        coulombCoupling = reaction.inelastic.r2s[4];  // r2s(4) — Coulomb VC value
        rcuEff = rC;

        for (i = 1; i <= numPoint; i++) {
            rValue = rptsPointer[i];
            wt   = rwtsPointer[i];
            reaction.gridData.rpts4Pointer[i] = (float)rValue;
            stepIndex = (int)(rValue / stepSize);
            x     = rValue - stepIndex * stepSize;

            // 652: derivative of W.S.
            xx = (splineBPointer[stepIndex] + x * (2.0*splineCPointer[stepIndex] + x*3.0*splineDPointer[stepIndex]));
            // term = + (12*E/H**2) * r2s(1) * wt * xx  — r2s(1) = real radius (inelastic)
            term = (12.0 * E / (reaction.distortedWave.channel[1].stepSize * reaction.distortedWave.channel[1].stepSize)) * reaction.inelastic.r2s[1] * wt * xx;
            nuclearHPointer[2 * i - 2] = term;

            if (rValue < rC) xx = std::pow(rValue / (rC * rC), (double)lxMin) / rC;
            if (rValue >= rC) xx = 1.0 / std::pow(rValue, (double)(lxMin + 1));
            coulombHPointer[i - 1] = -wt * coulombCoupling * xx;
        }

        //
        // NOW GO ON TO THE IMAGINARY PART
        //
        naturalCubicSpline(nSteps, splineRPointer,
               imv1Pointer + 1, splineBPointer, splineCPointer, splineDPointer);

        for (i = 1; i <= numPoint; i++) {
            rValue = rptsPointer[i];
            wt   = rwtsPointer[i];
            stepIndex = (int)(rValue / stepSize);
            x     = rValue - stepIndex * stepSize;

            // 732: derivative
            xx = (splineBPointer[stepIndex] + x*(2.0*splineCPointer[stepIndex] + x*3.0*splineDPointer[stepIndex]));
            // r2s(2) = imag radius (inelastic)
            term = (12.0 * E / (reaction.distortedWave.channel[1].stepSize * reaction.distortedWave.channel[1].stepSize)) * reaction.inelastic.r2s[2] * wt * xx;
            nuclearHPointer[2 * i - 1] = term;
        }

        // CCSW (problemType == 24) C.C. deformed-NUCH + COULOMB-form-factor block
        // Was ~170 lines of NNUCFF/NCOUFF radial-form generation that fed
        // into the CC NUCH/COULH pool slots.

        //
        //
        // FREE WORK AREAS
        //

        //
        // RESTORE FIRST CHANNEL
        //
        restoreChannelPotentials();

        if (verbosity != 0) {
            //
            // SUMMARIZE THINGS NICELY
            //
            std::printf("-ONE-DIMENSIONAL INTEGRATION GRID:\n"
                        "0NUM. PTS.   MAP TYPE   GAMMA    MINIMUM   \"MID. PT.\"   MAXIMUM\n");
            std::printf("%6d%11d%11.2f%10.2f%12.2f%14.2f\n",
                        savedNumPoint, mapSum, gammaSum, sumMin, sumMid, sumMax);
        }

        returnCode = 1;
        return;
    }  // end of deformed-potential scratch scope (auto-frees local vectors)
}

// ---------------------------------------------------------------------------
// Backward-compat free function wrapper — calls the class method
// ---------------------------------------------------------------------------
