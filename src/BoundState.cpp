// BoundState.cpp — BoundState class implementation: solve() (BOUND, Numerov
// eigenvalue search) and allocateFormFactor(). The form-factor / Coulomb-integral
// / JtoL methods live in the BoundState_*.cpp companion files.

#include "ptolemy_types.h"
#include "CoulombWaveFunction.h"
#include "Timing.h"
#include "linkule.h"
#include "print_utils.h"
#include "BoundState.h"
#include "Reaction.h"
#include "Constants.h"
#include "OpticalPotential.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>


// ============================================================================
// Part 1: Form-factor work array management
// ============================================================================
void BoundState::allocateFormFactor(int nSteps)
{
    if (nSteps <= 0) return;

    // Always reallocate fresh (VPHI is re-created each time BSSET is called)
    formFactor.assign(nSteps + 1, 0.0);

    // Note: vertex[1].jbdPointer / vertex[2].jbdPointer are set by the BSSET caller after this function
}

// zero callers remained.

// ============================================================================
// Part 2: Numerov bound-state solver (BOUND)
// ============================================================================
void BoundState::solve(int& returnCode, Reaction& reaction)
{
//
//     BOUND STATE PROGRAM
//
//  CODED 2/22/73  D.GLOECKNER
//     1/20/75 - FIX nStep CALCULATION - S. PIEPER
//
//     NOTE...  THE CHANNEL NUMBER (1 OR 2) IS INPUT VIA  boundChannel BELOW
//
//
//     ON OUTPUT:
//     returnCode IS 1 IF BOUND WAS SUCCESSFUL  OR
//          IS 0 IF BOUND WAS NOT SUCCESSFUL
//     IVOUT POINTS TO THE COMPLETE POTENTIAL
//     IWF2 POINTS TO THE WAVEFUNCTION
//

    // implicit real*8 (a-h, o-z)

    // --- Part 2 COMMON aliases (moved here to avoid goto-crosses-init) ---
    auto& accuracy = reaction.integrationGrid.accuracy;
    auto& Q      = reaction.energies.Q;
    constexpr int maxIteration = 10;
    auto& r0Mass = reaction.internalState.r0Mass;
    auto* eBnds  = reaction.internalState.eBnds;
    auto& J  = reaction.angMom.J;
    auto& spinTarget    = reaction.angMom.spinTarget;
    double rAsymp = 0, temp = 0, x = 0, xx = 0;
    double absPhi2 = 0, dPhi = 0, d2Phi = 0, xL = 0, xQ = 0, xP = 0, phiPrevious = 0;
    double centrifugal = 0;
    int loopBound = 0, trialIndex = 0, matchIndex = 0;

    const auto hbar_c   = Constants::hbar_c;
    const auto fine_structure_inv   = Constants::fine_structure_inv;

    auto& A       = reaction.opticalPotentialParams.A;
    auto& aSo     = reaction.opticalPotentialParams.aSo;
    auto& asymptopia  = reaction.integrationGrid.asymptopia;
    constexpr double delta_Vk = 0.05;
    auto& E       = reaction.energies.E;
    auto& aM      = reaction.masses.aM;
    auto& massProj     = reaction.masses.massProj;
    auto& massTgt     = reaction.masses.massTgt;
    auto& R       = reaction.integrationGrid.R;
    auto& rC      = reaction.opticalPotentialParams.rC;
    auto& rSo     = reaction.opticalPotentialParams.rSo;
    auto& stepSize  = reaction.integrationGrid.stepSize;
    auto& V       = reaction.opticalPotentialParams.V;
    auto& vSo     = reaction.opticalPotentialParams.vSo;
    auto& stepsPerUnit  = reaction.integrationGrid.stepsPerUnit;

    auto& L       = reaction.angMom.L;
    auto& nNodes   = reaction.angMom.nNodes;
    auto& printLevel  = reaction.flags.printLevel;
    auto& zProj     = reaction.charges.zProj;
    auto& zTarget     = reaction.charges.zTarget;
    auto& parity  = reaction.angMom.parity;
    auto* parityPt  = reaction.angMom.parityPt;  // 1-based

    auto& undefValue   = reaction.internalState.undefValue;
    int   notDefSentinel  = NOTDEF_INT;
    auto& boundChannel  = reaction.internalState.boundChannel;


    auto& jProj      = reaction.angMom.jProj;
    auto& spinProj     = reaction.angMom.spinProj;

    auto& uniqueLinkuleId  = reaction.linkuleData.uniqueLinkuleId;
    // linkuleAddr accessed as reaction.linkuleData.linkuleAddr[col][row] (1-based)

    // --- Local variables ---
    double tStart, uK, stepSize2, kinPrefactor, lLp1;
    double coulombInConst, coulombInR2Coef, coulombConst, vCoul, rValue;
    double lsCoupling, vSoSave;
    double rcSave;
    double delta, deltaV;
    double vStart;
    double eta;
    int nSteps, matchTopIndex, matchBotIndex;
    double* vOutPointer = nullptr;  // 0-based pointer into vertex[boundChannel].bsPotential
    int vIterCount, nodeIterCount;
    int actualNodeCount;
    int ii, i;
    // N1/N2 dropped — were woodsX out-params, never read here.
    int callStatus;


    // wavefunction & bsPotential now live on vertex[v] vectors, no pool slot to name.

    char channelWords[3][2][9];
    std::strncpy(channelWords[0][0], "PROJECTI", 9);
    std::strncpy(channelWords[0][1], "LE      ", 9);
    std::strncpy(channelWords[1][0], "    TARG", 9);
    std::strncpy(channelWords[1][1], "ET      ", 9);
    std::strncpy(channelWords[2][0], "        ", 9);
    std::strncpy(channelWords[2][1], "        ", 9);

    double uKs[4], etas[4], uK2s[4], vTrials[4], dRin[4];
    double dRout[4], phis[4], valueOut[4], sumOut[4], sumIn[4];
    int nodeCounts[4];
    // wfPointers: 1-based pointer array for WF1, WF2, WF3
    double* wfPointers[4] = {nullptr, nullptr, nullptr, nullptr};
    // Named vectors for the 4 internal arrays (allocated after nSteps known)
    std::vector<double> boundV1Vector, boundV2Vector, boundPhi1Vector, boundPhi3Vector;

    int ipTyps[5]; // 1-based
    ipTyps[1] = 1; ipTyps[2] = 3; ipTyps[3] = 5; ipTyps[4] = 6;

    char parityWord[3][9];
    std::strncpy(parityWord[0], "-1      ", 9);
    std::strncpy(parityWord[1], "UNKNOWN ", 9);
    std::strncpy(parityWord[2], "+1      ", 9);

    char linkId[5];
    linkId[0] = '*'; linkId[1] = '1'; linkId[2] = '0'; linkId[3] = '0';
    linkId[4] = '\0';
    char linkFmtOut[9];

    // LOGICAL  printSwitch, convergencePrintSwitch, hasVso
    bool printSwitch, convergencePrintSwitch, hasVso;

    tStart = (float)dtime_();
    convergencePrintSwitch = (printLevel / 100) % 10 >= 2;
    printSwitch = printLevel % 10 >= 1;
//
//   INITIALIZE PARAMETERS
//
    if (boundChannel > 9) boundChannel = 1;
    // NAMES_arr (pool-slot label suffix) deleted with the slots themselves.
//
    hasVso = (vSo != 0.0);
//
//
    if (V == 0.0) V = 60.0;
//
//     TEST THE INPUT AND RETURN  0  COMPLETION CODE IF BAD
//
    returnCode = 1;
    if (!(R > 0.0 && A > 0.0)
        && !(reaction.linkuleData.linkuleAddr[1][3] > 0 || reaction.linkuleData.linkuleAddr[6][3] > 0)) {
        std::printf("\n **** R OR A HAS INVALID VALUE: %15.5G%15.5G\n", R, A);
        returnCode = 0;
    }
    if (hasVso) {
    if (!(rSo > 0.0 && aSo > 0.0)
        && !(reaction.linkuleData.linkuleAddr[3][3] > 0)) {
        std::printf("\n **** RSO OR ASO HAS INVALID VALUE: %15.5G%15.5G\n", rSo, aSo);
        returnCode = 0;
    }
    if (jProj == reaction.internalState.notDefSentinel) {
        std::printf("\n **** JP MUST BE DEFINED FOR SPIN ORBIT FORCE.\n");
        returnCode = 0;
    }
    if (spinProj == reaction.internalState.notDefSentinel) {
        std::printf("\n **** WARNING:  SP WAS NOT DEFINED; "
                    "IT IS ASSUMED TO BE 1/2 FOR THE SPIN-ORBIT FORCE.\n");
        spinProj = 1;
    }
    if (!((int)jProj <= (int)spinProj + 2*L && (int)jProj >= std::abs((int)spinProj - 2*L)
          && ((int)jProj + (int)spinProj) % 2 == 0)) {
        std::printf("\n **** ERROR: L, SP AND JP DO NOT FORM A TRIANGLE OR "
                    "THEY ARE MIXED INTEGER AND HALF-INTEGER: %5d%5d/2%5d/2\n",
                    L, (int)spinProj, (int)jProj);
        returnCode = 0;
    }
    }
    if (V == 0.0) V = 50.0;
    if (!(E < 0.0 && V > 0.0 && aM > 0.0 && E != undefValue)) {
        std::printf("\n **** E OR V OR M (REDUCED MASS) HAS INVALID"
                    " VALUE: %15.5G%15.5G%15.5G\n", E, V, aM);
        returnCode = 0;
    }
    if (!(R < asymptopia || reaction.linkuleData.linkuleAddr[1][3] > 0)) {
        std::printf("\n **** R MUST BE LESS THAN ASYMTOPIA: %15.5G%15.5G\n", R, asymptopia);
        returnCode = 0;
    }
    rcSave = rC;
    if (rC == undefValue && zProj * zTarget == 0) rC = 1.0;
    if (!(rC > 0.0 || reaction.linkuleData.linkuleAddr[5][3] > 0)) {
        std::printf("\n **** RC IS INVALID: %15.5G\n", rC);
        returnCode = 0;
    }
    if (!(reaction.linkuleData.linkuleAddr[6][3] > 0)
        && !(L != notDefSentinel && nNodes != notDefSentinel)) {
        std::printf("\n **** BOTH L AND NODES MUST BE DEFINED.\n");
        returnCode = 0;
    }
    i = parity * parityPt[1] * parityPt[2];
    if (i != 0) {
        i = (i + 3) / 2 + L;
        if (i % 2 != 0) {
            std::printf("\n **** L AND PARITIES ARE INCOMPATABLE: %5d%5d%5d%5d\n",
                        L, parity, parityPt[1], parityPt[2]);
            returnCode = 0;
        }
    }
//
    vIterCount = 0;
    nodeIterCount = 0;
//
    // Shared solve()-exit epilogue: restore rC and report CPU time (printed only,
    // not curve data). Used by the returnCode==0 early return and the final return.
    auto finishSolve = [&]() {
        rC = rcSave;
        tStart = (float)dtime_() - tStart;
        if (convergencePrintSwitch) std::printf(" BOUND STATE CPU TIME =%7.3f SECONDS.\n", tStart);
    };
    // Identical post-linkule failure guard repeated after every linkule() call;
    // returns true (caller then returns) when the call reported callStatus < 0.
    auto bailIfCallFailed = [&]() {
        if (callStatus < 0) { returnCode = 0; return true; }
        return false;
    };
//
    if (returnCode == 0) {
        finishSolve();
        return;
    }
//
//  K IS
    uK = std::sqrt(-2.0 * aM * E) / hbar_c;
    if (stepsPerUnit != undefValue)
        stepSize = std::min(1.0 / uK, A) / stepsPerUnit;
    stepSize2 = stepSize * stepSize;
    kinPrefactor = 2.0 * aM / (hbar_c * hbar_c);
    lLp1 = (double)L * (L + 1);
//
    nSteps = (int)(asymptopia / stepSize + 1.5);
//
//     WE DO THE MATCHING BEYOND BOTH THE REAL RADIUS
//     AND THE LAST NODE.
//
    matchTopIndex = (int)(std::max(R, A) / stepSize + 0.5);
    matchTopIndex = std::max(matchTopIndex, 20);
    matchBotIndex = std::max(matchTopIndex - 20, 3);
//
    actualNodeCount = nNodes;
//
    boundV1Vector.assign(nSteps + 1, 0.0);   // BOUNDV1: indices 0..nSteps
    boundV2Vector.assign(nSteps + 1, 0.0);   // BOUNDV2: indices 0..nSteps
    boundPhi1Vector.assign(nSteps + 2, 0.0);  // BNDPHI1: indices 0..nSteps+1
    boundPhi3Vector.assign(nSteps + 2, 0.0);  // BNDPHI3: indices 0..nSteps+1
    vertex[boundChannel].bsPotential.assign(nSteps, 0.0);
    vertex[boundChannel].wavefunction.assign(nSteps, 0.0);
//
//     INITIALIZE THE LINKULES FOR CENTRAL, COULOMB, SPIN-ORBIT,
//     AND WAVE FUNCTION.
//
    for (ii = 1; ii <= 4; ii++) {
        i = ipTyps[ii];
        if (reaction.linkuleData.linkuleAddr[i][3] != 0) {
            uniqueLinkuleId = uniqueLinkuleId + 1;
            std::snprintf(linkFmtOut, sizeof(linkFmtOut), "*%03d", uniqueLinkuleId);
            std::memcpy(linkId, linkFmtOut, 4);
            double dummy;
            linkule(reaction.linkuleData.linkuleAddr[i][3],
                   *(char8*)&reaction.linkuleData.linkuleAddr[i][1],
                   &reaction.linkuleData.linkuleAddr[i][5],
                   i, 1,
                   callStatus, L, jProj, 0.0, stepSize, nSteps, &dummy, &dummy, linkId, reaction);
            if (bailIfCallFailed()) return;
            reaction.linkuleData.linkuleAddr[i][4] = returnCode;
        }
    }
//
//
    // Setup wfPointers: pointer array for wave functions (replaces LWFS int array)
    wfPointers[1] = boundPhi1Vector.data();   // BNDPHI1
    wfPointers[2] = vertex[boundChannel].wavefunction.data();   // PHIX (0-based view)
    wfPointers[3] = boundPhi3Vector.data();   // BNDPHI3
    // 0-based view into bsPotential (vOutPointer[i-1] = bsPotential[i-1])
    vOutPointer = vertex[boundChannel].bsPotential.data();
//
//
//
//  fine_structure_inv IS FINE STRUCTURE CONSTANT
// NEED coulombInConst,coulombInR2Coef FOR COULOMB POTENTIAL INSIDE NUCLEAR RADIUS
//
    coulombInConst = zProj * zTarget * hbar_c * 1.5 / (fine_structure_inv * rC);
    coulombInR2Coef = zProj * zTarget * hbar_c * 0.5 / (fine_structure_inv * rC * rC * rC);
//
//  C0 IS COULOMB OUTSIDE OF POTENTIAL
//
    coulombConst = zProj * zTarget * hbar_c / fine_structure_inv;
//
//
//  CALCULATE POTENTIAL AT ALL STEPS
//
//     THE ARRAYS ARE
//
//          D      AND vSo PART OF SPIN ORBIT.
//                  WELL AND TAU PART OF SPIN ORBIT.
//
//     TOTAL POTENTIAL IS
//
//     IN EACH CASE
//
//
//     GET THE COULOMB POTENTIAL.
//
    // Site 1: Coulomb potential. LINKUL or analytic coulombInConst/coulombInR2Coef/coulombConst.
    //
    // SCF potential sites 1/2/3 below share an identical LINKUL call (request
    // code 3), differing only in the linkuleAddr slot index and the destination
    // vector; deduped into one [&] lambda (the caller still issues the return).
    auto callBoundScfLinkule = [&](int idx, double* dataArr) -> bool {
        double neg1 = -1.0;
        linkule(reaction.linkuleData.linkuleAddr[idx][3],
               *(char8*)&reaction.linkuleData.linkuleAddr[idx][1],
               &reaction.linkuleData.linkuleAddr[idx][5],
               idx, 3, callStatus, L, jProj,
               0.0, stepSize, nSteps, dataArr, &neg1, (char*)nullptr, reaction);
        return bailIfCallFailed();
    };
    if (reaction.linkuleData.linkuleAddr[5][3] > 0) {
        if (callBoundScfLinkule(5, boundV1Vector.data() - 1)) return;
    } else {
        rValue = -stepSize;
        for (int stepIndex = 1; stepIndex <= nSteps; stepIndex++) {
            rValue = rValue + stepSize;
            if (rValue < rC) {
                vCoul = coulombInConst - coulombInR2Coef * (rValue * rValue);
            } else {
                vCoul = coulombConst / rValue;
            }
            boundV1Vector[stepIndex - 1] = -vCoul;
        }
    }
//
//     INITIALIZE THE FIXED PART OF THE EFFECTIVE TRANSITION OPERATOR.
//
    for (int stepIndex = 1; stepIndex <= nSteps; stepIndex++) {
        vOutPointer[stepIndex - 1] = boundV1Vector[stepIndex - 1];
    }

    // Site 2: Central potential. LINKUL or woodsX.
    if (reaction.linkuleData.linkuleAddr[1][3] > 0) {
        if (callBoundScfLinkule(1, boundV2Vector.data() - 1)) return;
    } else {
        // resize() zeroes the buffer so ADD-semantics fillWoodsSaxon == old SET.
        OpticalPotential pot;
        pot.resize(nSteps, 0.0, stepSize);
        pot.fillWoodsSaxon(-1.0, R, A);
        std::copy(pot.values.begin(), pot.values.end(), boundV2Vector.begin());
    }

    // Site 3: Spin-orbit potential (if any). LINKUL or woodsX, then split.
    if (hasVso) {
    // sigma.L.  Note that jProj, two_spin_proj are doubled.
    lsCoupling = jProj * (jProj + 2) - spinProj * (spinProj + 2);
    lsCoupling = (0.25 * lsCoupling - lLp1) / spinProj;
    lsCoupling = vSo * lsCoupling;
    vSoSave = vSo;
    vSo = 1.0;
    if (reaction.linkuleData.linkuleAddr[3][3] > 0) {
        if (callBoundScfLinkule(3, boundPhi1Vector.data() - 1)) return;
    } else {
        // resize() zeroes the buffer so ADD-semantics fillSpinOrbit == old SET;
        // the lsCoupling multiply below stays AFTER the fill, as before.
        OpticalPotential pot;
        pot.resize(nSteps, 0.0, stepSize);
        pot.fillSpinOrbit(-1.0, rSo, aSo);
        std::copy(pot.values.begin(), pot.values.end(), boundPhi1Vector.begin());
    }
    for (int stepIndex = 1; stepIndex <= nSteps; stepIndex++) {
        boundV1Vector[stepIndex - 1] += lsCoupling * boundPhi1Vector[stepIndex - 1];
        vOutPointer[stepIndex - 1] += lsCoupling * boundPhi1Vector[stepIndex - 1];
    }
    vSo = vSoSave;
    }  // end if (hasVso) — spin-orbit setup
//
//
//     NOW CALCULATE THE WAVE FUNCTION.
//
//     IF A LINKULE DOES IT, SKIP MOST OF THE FOLLOWING.
//
    // The two [6]-LINKULE invocations (wavefunction-compute below at requestCode
    // 3, print-self in the finalization at requestCode 2) are byte-identical
    // except the request code — share via a lambda taking it as a param.
    auto callBoundLinkule = [&](int requestCode) {
        linkule(reaction.linkuleData.linkuleAddr[6][3],
               *(char8*)&reaction.linkuleData.linkuleAddr[6][1],
               &reaction.linkuleData.linkuleAddr[6][5],
               6, requestCode,
               callStatus, L, jProj, 0.0, stepSize, nSteps,
               wfPointers[2] - 1, boundV2Vector.data() - 1, (char*)nullptr, reaction);
    };

    // do-while-false wraps the LINKUL/SCF section so the 4 former
    // `goto L850;` early-exits become `break;` to the finalization below.
    do {
    if (reaction.linkuleData.linkuleAddr[6][3] > 0) {
        callBoundLinkule(3);
        if (bailIfCallFailed()) return;
        sumIn[2] = 1.0;
        sumOut[2] = 1.0;
        matchIndex = nSteps / 2;
        actualNodeCount = nNodes;
        break;
    }
//
//     WE FIND SOLUTIONS FOR THREE VALUES OF V OR K.  DELTA IS THE
//     RELATIVE SPACING TO USE.  AS APPROPRIATE, DELTAV OR DELTAK ARE
//     THE ABSOLUTE SPACINGS TO USE.  SOLUTIONS ARE FOUND FOR
//     V-DELTAV, V, V+DELTAV   OR   K-DELTAK, K, K+DELTAK
//
    // 3-level SCF iteration as nested while-true:
    //   outer = was L3020 (each starting V/K trial; L6666 retry feeds here)
    //   middle = was L3030 (was: each K within current V; IFIT==0 retry path
    //                       — now permanently 1-trip since IFIT is always 1)
    //   inner = was L500   (each Numerov pass with new V)
    bool exitDoWhileScf = false;
    while (true) {
    delta = delta_Vk;
    if (convergencePrintSwitch) std::printf(" %5d%16.8G%16.8G%16.8G\n", vIterCount, V, uK, delta);
    vStart = V;
//
//     SET K+, K, K- AND CORRESPONDING eta_ch EACH TIME E CHANGES
//     (THIS IS FOR WHEN E AND NOT V IS BEING VARIED)
//
    while (true) {
    deltaV = V * delta;
    uKs[1] = uK;
    uKs[2] = uK;
    uKs[3] = uK;
    for (int trialIndex = 1; trialIndex <= 3; trialIndex++) {
        uK2s[trialIndex] = uKs[trialIndex] * uKs[trialIndex];
        etas[trialIndex] = zProj * zTarget * aM / (hbar_c * fine_structure_inv * uKs[trialIndex]);
    }
    eta = etas[2];
//
//     FIND ASYMPTOTIC FORM AS A WHITTAKER FUNCTION.
//


    // --- BOUND PART 2 ---

    // (Part 2 aliases and locals moved to top of function)

    // wfPointers already set up above (WF1/WF3 as vectors, WF2 as ALLOC)
    // (wfPointers replaces old LWFS int array)

//
//     FIND ASYMPTOTIC FORM AS A WHITTAKER FUNCTION.
//     The IASYMP=1 (SKIPASYMP) shortcut — use 2 and 1 as the last two
//     function values and extend asymptopia by ~5 fm — was dropped
//
    rAsymp = asymptopia - stepSize;
//
//     WE NEED A NORM term THAT IS SMALL FOR SMALL RHO AND
//     ETA+L+1 FOR LARGE RHO.  WHEN THE WHITTAKER IS DIVIDED BY
//     GAMMA(THIS NORM) IT WILL GENERALLY NOT BE TOO LARGE OR SMALL
//
    temp = rAsymp * uK;
    temp = (L + eta + 1) * std::pow(temp / (2.0 + temp), 2);
//
    // IFIT permanently 1: only i==1 evaluates penetrability; i==2,3 fall
    // through to the wfPointers write with the x from i==1.
    for (ii = 1; ii <= 2; ii++) {
        for (i = 1; i <= 3; i++) {
            if (i == 1) {
                x = CoulombWaveFunction::penetrability(L, etas[i], uKs[i] * rAsymp, temp);
//
//     COULNG RETURNS THE WHITTAKER times A factor.  HERE WE
//     PATCH UP THE RATIO BETWEEN W(asymptopia-stepSize) AND W(asymptopia)
//
                if (ii == 1) x = x *
                    std::exp(uKs[i] * stepSize + etas[i] * std::log(asymptopia / rAsymp));
            }
            wfPointers[i][nSteps - 3 + ii] = x;
        }
        rAsymp = rAsymp + stepSize;
    }
    if (convergencePrintSwitch) std::printf(" ASYMPTOTIC FORM AT LAST TWO STEPS: %18.8G%18.8G\n",
        wfPointers[2][nSteps - 2], wfPointers[2][nSteps - 1]);
//
//  NOW HAVE ASYMPTOTIC WAVEFUNCTIONS
//
//
//
//
//
//  STEP INTO MATCHING RADIUS
//
//     FIND THE MAXIMUM VALUE OF THE MINIMUM OVER V-DELTA, V, V+DELTA,
//     OF THE WAVEFUNCTION IN THE INNER 20 STEPS AND USE AS THE MATCHING
//     RADIUS.
//
//
//     WE COME HERE FOR EACH NEW VALUE OF V
//
    while (true) {
    vTrials[1] = V - deltaV;
    vTrials[2] = V;
    vTrials[3] = V + deltaV;
//
    xx = 0.0;
    loopBound = nSteps - matchBotIndex + 1;
    for (i = 3; i <= loopBound; i++) {
        ii = nSteps - i + 1;
        centrifugal = lLp1 / std::pow(stepSize * ii, 2);
        temp = Constants::bigNum;
        for (trialIndex = 1; trialIndex <= 3; trialIndex++) {
            x = vTrials[trialIndex] * boundV2Vector[ii] + boundV1Vector[ii];
            x = 2.0 + stepSize2 * (uK2s[trialIndex] - kinPrefactor * x + centrifugal);
            wfPointers[trialIndex][ii - 1] = x * wfPointers[trialIndex][ii] - wfPointers[trialIndex][ii + 1]; // K-1,K,K+1
            temp = std::min(temp, std::fabs(wfPointers[trialIndex][ii]));
        }
        if (ii <= matchTopIndex && temp > xx) {
            matchIndex = ii;
            xx = temp;
        }
    }
//
//     COMPUTE OUTER MATCHING PARAMETERS -- DERIVATIVES, NORM
//     OF OUTER PART, VALUE AT MATCHING POINT AND NUMBER OF nNodes.
//
    for (trialIndex = 1; trialIndex <= 3; trialIndex++) {
        // was: K = LWFS[trialIndex]+matchIndex
        dRout[trialIndex] = (0.5 / stepSize) * (wfPointers[trialIndex][matchIndex + 1] - wfPointers[trialIndex][matchIndex - 1]);
        valueOut[trialIndex] = wfPointers[trialIndex][matchIndex];
        nodeCounts[trialIndex] = 0;
        sumIn[trialIndex] = 0.0;
        sumOut[trialIndex] = 0.5 * (std::pow(wfPointers[trialIndex][matchIndex], 2) + std::pow(wfPointers[trialIndex][nSteps - 1], 2));
        loopBound = matchIndex + 3;
        for (ii = loopBound; ii <= nSteps; ii++) {
            // was: K = LWFS[trialIndex]-2+ii → index ii-2
            sumOut[trialIndex] = sumOut[trialIndex] + std::pow(wfPointers[trialIndex][ii - 2], 2);
            if (wfPointers[trialIndex][ii - 2] * wfPointers[trialIndex][ii - 3] <= 0.0) nodeCounts[trialIndex] = nodeCounts[trialIndex] + 1;
        }
    }
//
//  FIND POWER SERIES SOLUTION (KEEP FIRST TWO TERMS) IN ORDER
//  TO FIX FUNCTION AT FIRST STEP.  THEN STEP TO MATCHING RADIUS
//
//     INCLUDE IT.
//
    for (trialIndex = 1; trialIndex <= 3; trialIndex++) {
        wfPointers[trialIndex][0] = 0.0;
//
//     POTENTIAL + K**2  AT R = stepSize.
//
        temp = -kinPrefactor * (boundV1Vector[1] + vTrials[trialIndex] * boundV2Vector[1])
            + uK2s[trialIndex];
        wfPointers[trialIndex][1] = stepSize * (1.0 + stepSize2 * temp / 6.0);
    }
//
//  STEP OUT TO MATCHING RADIUS
//
    for (ii = 1; ii <= matchIndex; ii++) {
        centrifugal = lLp1 / std::pow(stepSize * ii, 2);
        for (trialIndex = 1; trialIndex <= 3; trialIndex++) {
            // was: K = LWFS[trialIndex]+ii
            x = vTrials[trialIndex] * boundV2Vector[ii] + boundV1Vector[ii];
            x = 2.0 + stepSize2 * (uK2s[trialIndex] - kinPrefactor * x + centrifugal);
            wfPointers[trialIndex][ii + 1] = x * wfPointers[trialIndex][ii] - wfPointers[trialIndex][ii - 1];
            sumIn[trialIndex] = sumIn[trialIndex] + std::pow(wfPointers[trialIndex][ii], 2);
            if (wfPointers[trialIndex][ii] * wfPointers[trialIndex][ii - 1] < 0.0) nodeCounts[trialIndex] = nodeCounts[trialIndex] + 1;
        }
    }
//
//     MATCH SOLUTIONS;  DERIVATIVES FOR INNER SIDE.
//
    for (trialIndex = 1; trialIndex <= 3; trialIndex++) {
        // was: K = LWFS[trialIndex]+matchIndex
        sumIn[trialIndex] = sumIn[trialIndex] - 0.5 * std::pow(wfPointers[trialIndex][matchIndex], 2);
//
//     sumOut MULTIPLIES THE OUTER SOLUTION TO NORMALIZE IT
//     sumIn  MULTIPLIES THE INNER SOLUTION TO NORMALIZE IT
//
        xx = wfPointers[trialIndex][matchIndex] / valueOut[trialIndex];
        x = 1.0 / std::sqrt(stepSize * (sumIn[trialIndex] + xx * xx * sumOut[trialIndex]));
        sumIn[trialIndex] = x;
        sumOut[trialIndex] = xx * x;
        valueOut[trialIndex] = (xx * x) * valueOut[trialIndex];
        dRout[trialIndex] = (xx * x) * dRout[trialIndex];
        dRin[trialIndex] = x * (0.5 / stepSize) * (wfPointers[trialIndex][matchIndex + 1] - wfPointers[trialIndex][matchIndex - 1]);
        phis[trialIndex] = dRout[trialIndex] - dRin[trialIndex];
    }
//
//     WE LOOK FOR A ZERO OF THE DIFFERENCE OF THE INNER AND OUTER
//     DERIVATIVES.
//
    actualNodeCount = nodeCounts[2];
    x = matchIndex * stepSize;
    if (convergencePrintSwitch) std::printf(" DERIVATIVES AT R =%7.2f FM"
        "%57sV =%15.6G%15.6G%15.6G\n"
        "     INNER =%16.6G%16.6G%16.6G"
        "%35sK =%15.6G%15.6G%15.6G\n"
        "    OUTER =%16.6G%16.6G%16.6G"
        "%32sU(R) =%15.6G%15.6G%15.6G\n"
        "       DIF =%16.6G%16.6G%16.6G"
        "%31sNODES =%9d%6s%9d%6s%9d%6s\n",
        x, "", vTrials[1], vTrials[2], vTrials[3],
        dRin[1], dRin[2], dRin[3],
        "", uKs[1], uKs[2], uKs[3],
        dRout[1], dRout[2], dRout[3],
        "", valueOut[1], valueOut[2], valueOut[3],
        phis[1], phis[2], phis[3],
        "", nodeCounts[1], "", nodeCounts[2], "", nodeCounts[3], "");
    if (!(std::fabs(phis[2]) < accuracy)) {
        absPhi2 = std::fabs(phis[2]);
        if (convergencePrintSwitch) std::printf(" FOR vIterCount=%3d THE  FUNCTION IS%11.3G >%11.3G\n",
            vIterCount, absPhi2, accuracy);
        if (vIterCount > maxIteration) {
            if (nNodes != actualNodeCount) break;  // — fall through to L6666 block below
            std::printf(" ******** NO SOLUTION FOUND *********\n");
            returnCode = 0;
            exitDoWhileScf = true;
            break;
        }
        dPhi = (phis[3] - phis[1]) / (2.0 * delta);
        d2Phi = (phis[3] - 2.0 * phis[2] + phis[1]) / (delta * delta);
        xL = -phis[2] / dPhi;
        xQ = 0.0;
        x = xL;
    //
    //  IF d2Phi IS ZERO DO LINEAR APPROXIMATION
    //
        if (!(std::fabs(d2Phi) < 1.0e-5 * std::fabs(dPhi))) {
            double discriminant = (dPhi * dPhi) - 2.0 * phis[2] * d2Phi;
            if (!(discriminant < 0.0)) {
    //
                xQ = (-dPhi + (dPhi / std::fabs(dPhi)) * std::sqrt(discriminant)) / d2Phi;
                x = xQ;  // (LINEAR != 1) gate dropped — permanently 0.
            }
        }
    //
        if (std::fabs(x) > 0.5) x = std::copysign(0.50, x);
        if (convergencePrintSwitch) std::printf(" dPhi,d2Phi =%14.5G%14.5G\n"
            " LINEAR, QUADRATIC METHODS =%12.3G%12.3G"
            "     REL. STEP OF V OR K USED =%14.5G\n",
            dPhi, d2Phi, xL, xQ, x);
    //
    //     FIND A NEW DELTAVK.  WE USE THE STEP LENGTH USED IN THIS
    //     ITERATION UNLESS IT STAYED WITHIN THE THE PREVIOUS DELTA,
    //     INWHICH CASE WE USE .5 OF IT.
    //
        xx = delta;
        delta = std::fabs(x);
        if (delta < xx) delta = 0.5 * delta;
        if (delta < 1.0e-7) delta = 1.0e-7;
        delta = std::min(delta, 0.20);
        // IFIT permanently 1: take only the FIT-V0 branch.
        xx = V * (1.0 + x);
        if (xx > 0.0) V = xx;
        if (xx <= 0.0) V = 0.5 * V;
        vIterCount = vIterCount + 1;
        if (convergencePrintSwitch) std::printf("\nITERATION%3d     V =%16.8G  K =%16.8G"
            "     DELTAVK =%12.4G\n", vIterCount, V, uK, delta);
        deltaV = delta * V;
        continue;   // — retry Numerov with new V
    }
//
//     HAVE CONVERGED ON A V OR UK - CHECK FOR PROPER NUMBER OF nNodes
//
    absPhi2 = std::fabs(phis[2]);
    if (convergencePrintSwitch) std::printf(" ON ITERATION %3d THE FUNCTION IS%11.3G <%11.3G\n",
        vIterCount, absPhi2, accuracy);
    E = -uK2s[2] * hbar_c * hbar_c / (2.0 * aM);
    if (convergencePrintSwitch) std::printf("\n   FIT LOCATED WITH V=%14.8f AND E=%14.8f\n\n",
        V, E);
//
//
//  RENORMALIZE WAVEFUNCTION AND CHECK nNodes
//
    if (actualNodeCount == nNodes) {
        if (convergencePrintSwitch) std::printf("  SOLUTION HAS %1d NODE(S). L=%2d JP=%2d/2\n",
            nNodes, L, (int)jProj);
        exitDoWhileScf = true;
        break;   // exit L500
    }
//
//     SOLUTION HAS WRONG NUMBER OF nNodes - TRY ANOTHER V OR UK
//
    break;   // exit L500 — fall through to L6666 block below (wrong-nodes branch)

    }  // end while

    // L3030 retry path was IFIT==0 — gone since IFIT permanently 1.
    // The middle while-true now executes exactly one trip.
    break;

    }  // end while

    if (exitDoWhileScf) break;                    // exit L3020 → exit do-while-false

    // — wrong-nodes branch; reached when L500 breaks without setting exitDoWhileScf
    std::printf(" FOR V =%8.2f,  E =%7.2f, WAVEFUNCTION HAS"
        "%3d NODES BUT%3d ARE DESIRED.\n", V, E, actualNodeCount, nNodes);
    if (nodeIterCount > maxIteration) {
        std::printf("\n**** COULD NOT FIND A SOLUTION WITH%3d NODES.\n", nNodes);
        returnCode = 0;
        break;   // exit L3020 → exit do-while-false
    }
    vIterCount = 0;
    nodeIterCount = nodeIterCount + 1;
//
//     WE USE
//
//       X = DESIRED NUMBER OF nNodes
//       XQ = NUMBER OF nNodes WE JUST GOT
//       XP = PREVIOUS NUMBER OF nNodes
//
    x = nNodes;
    xQ = actualNodeCount;
//
//     WE USE THE STARTING OR FINAL V OR K, WHICH EVER IS CLOSER
//     TO WHERE WE WANT TO GO.
//
    if (nNodes < actualNodeCount) {
        V = std::min(V, vStart);
    } else {
        V = std::max(V, vStart);
    }
    phis[1] = V;
    if (nodeIterCount > 1 && xP != xQ) {
//
//     HAVE TWO POINTS - USE STRAIGHT LINE
//
        phis[2] = phis[1] + (x - xQ) * (phis[1] - phiPrevious) / (xQ - xP);
    } else {
//
//     WE HAVE ONLY ONE POINT TO WORK WITH
//
        // IFIT permanently 1: take only the FIT-V0 branch.
        phis[2] = (x + 1.0) * phis[1] / (xQ + 1.0);
    }
//
//     DO NOT ALLOW TO GO NEGATIVE
//
    if (phis[2] <= 0.0) phis[2] = 0.5 * phis[1];
//
    phiPrevious = phis[1];
    xP = xQ;
    // IFIT permanently 1: take only the FIT-V0 branch.
    V = phis[2];
    std::printf(" WILL TRY A NEW V:  %14.5G\n", V);
    continue;   // — try a new starting V/UK
    }  // end while
    } while (false);  // end do-while wrap
//
//
    rValue = stepSize;
    loopBound = matchIndex + 1;
    for (ii = 1; ii <= loopBound; ii++) {
        wfPointers[2][ii] = sumIn[2] * wfPointers[2][ii] / rValue;
        rValue = rValue + stepSize;
    }
    loopBound = loopBound + 2;
    for (ii = loopBound; ii <= nSteps; ii++) {
        wfPointers[2][ii - 1] = sumOut[2] * wfPointers[2][ii - 1] / rValue;
        rValue = rValue + stepSize;
    }
//
//     FOR L = 0 GET WAVEFUNCTION AT ORIGIN BASED ON STARTING
//     VALUE DEFINED ABOVE
//
    if (L == 0) wfPointers[2][0] = sumIn[2];
//
//     FOR L=0 WAVEFUNCTIONS FROM LINKULES, WE EXTRAPOLATE, USING THE
//
    if (L == 0 && reaction.linkuleData.linkuleAddr[6][3] > 0)
        wfPointers[2][0] = (4.0 * wfPointers[2][1] - wfPointers[2][2]) / 3.0;
//
//
//     PRINT RESULTS NICELY
//

    if (printSwitch) {
    i = (boundChannel > 2) ? 3 : boundChannel;
    uK = std::sqrt(-2.0 * aM * E) / hbar_c;
    eta = zProj * zTarget * aM / (hbar_c * fine_structure_inv * uK);
    std::printf("\n0%8s%.8s%.2s BOUND STATE PARAMETERS\n",
        "", channelWords[i-1][0], channelWords[i-1][1]);
    std::printf("0E =%10.4f MEV     KAPPA =%8.5f\n", E, uK);
    if (massProj != undefValue && massTgt != undefValue) {
        std::printf(" PROJECTILE MASS =%7.2f AMU     "
            "TARGET MASS =%7.2f AMU     "
            "REDUCED MASS =%10.2f MEV/C**2\n", massProj, massTgt, aM);
    } else {
        std::printf(" REDUCED MASS =%10.2f MEV/C**2\n", aM);
    }
    std::printf(" L =%3d%9d NODES\n", L, actualNodeCount);
    if ((int)spinProj != NOTDEF_INT && (int)spinTarget != NOTDEF_INT)
        std::printf(" PROJECTILE SPIN =%3d/2     "
            "TARGET SPIN =%3d/2\n", (int)spinProj, (int)spinTarget);
    if ((int)jProj != NOTDEF_INT)
        std::printf(" J PROJECTILE =%3d/2\n", (int)jProj);
    if ((int)J != NOTDEF_INT)
        std::printf(" TOTAL J =%3d/2\n", (int)J);
    std::printf(" PROJECTILE PARITY IS %.7s"
        "   TARGET PARITY IS %.7s"
        "   TOTAL PARITY IS %.7s\n",
        parityWord[parityPt[1] + 1], parityWord[parityPt[2] + 1], parityWord[parity + 1]);
    std::printf(" Z PROJECTILE =%4d     Z TARGET =%4d\n"
        "0POTENTIAL         COUPLING CONS.    RADIUS    DIFFUSENESS    RADIUS PARAM.\n",
        zProj, zTarget);
    parameterPrint(1, "REAL CENTRAL      ", V, R, A, R/r0Mass, reaction);
    parameterPrint(3, "REAL SPIN-ORBIT   ", vSo, rSo, aSo, rSo/r0Mass, reaction);
    parameterPrint(5, "PT.&SPHERE COULOMB", eta, rC, 0.0, rC/r0Mass, reaction);
//
    temp = std::min(1.0 / uK, A) / stepSize;
    std::printf("0ASYMPTOPIA = %8.3f FM\n"
        " STEP SIZE =  %8.3f FM     %8.3f STEPS PER \"ASYMPTOTIC RANGE\"\n\n",
        asymptopia, stepSize, temp);
//
//     PRINT STUFF ABOUT ITSELF.
//
    if (reaction.linkuleData.linkuleAddr[6][3] != 0) {
        std::printf("0 THE BOUND WAVEFUNCTION WAS COMPUTED BY THE %.8s LINKULE:\n",
            (char*)&reaction.linkuleData.linkuleAddr[6][1]);
        callBoundLinkule(2);
        if (bailIfCallFailed()) return;
        std::printf("\n");
    }
    }
//
//
//     SAVE THE POTENTIAL AND OTHER STUFF FOR GRDSET
//
    if (boundChannel >= 1 && boundChannel <= 2) {
    vertex[boundChannel].lBound = L;
    vertex[boundChannel].nodeCount = nNodes;
    vertex[boundChannel].jB = (int)jProj;
    vertex[boundChannel].boundMx = asymptopia;
    vertex[boundChannel].bsVstep = stepSize;
    vertex[boundChannel].nSpBd = nSteps;
    vertex[boundChannel].bsMass = aM;
    eBnds[boundChannel] = E;
    vertex[boundChannel].massRatio = reaction.internalState.ratMass;
    // ff.LNKADB[i][boundChannel][ii] = reaction.linkuleData.linkuleAddr[i][ii] copy loop
//
//     CAN WE COMPUTE THE Q OF THE REACTION
//
    if (eBnds[3 - boundChannel] != undefValue) {
        Q = eBnds[1] - eBnds[2];
        if (reaction.internalState.stripPickup == -1) Q = -Q;
    }
    }
//
//
//     CONSTRUCT THE COMPLETE POTENTIAL FOR V PHI
//
    for (i = 1; i <= nSteps; i++) {
        vOutPointer[i - 1] = -(vOutPointer[i - 1]
            + V * boundV2Vector[i - 1]);
    }
//
//     IF REQUESTED, PRINT THE WAVEFUNCTION
//
    rValue = 0.0;
//
    finishSolve();
    return;
//
//
//
}

// ---------------------------------------------------------------------------
// Backward-compat free function wrapper — calls the class method
// ---------------------------------------------------------------------------

// ============================================================================
// Part 3: Form-factor eval, Coulomb integrals, J->L coupling
// ============================================================================







