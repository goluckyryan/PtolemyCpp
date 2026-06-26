// ScatteringSolver_setupScatteringWaves.cpp — WAVSET: sets up scattering-wave
// parameters and dispatches the Numerov solve over all partial waves (L, jProj).

#include "ptolemy_types.h"
#include "ScatteringSolver.h"
#include "coulomb_utils.h"
#include "linkule.h"
#include "print_utils.h"
#include "Reaction.h"
#include "Constants.h"
#include <cstdio>
#include <cstring>
#include <cmath>



// ============================================================================
// SECTION 3: wavefunction_setup.cpp — ScatteringSolver::setupScatteringWaves, WAVSET wrapper
// ============================================================================

// setupScatteringWaves: setup scattering wave parameters.

void ScatteringSolver::setupScatteringWaves(int& returnCode, int isStandalone, Reaction& reaction)
{
    static const char channelName[2][9] = { "INCOMING", "OUTGOING" };
    static const char statNames[2][9] = { "BOSONS  ", "FERMIONS" };
    static const char soTNames[2][9] = { "SPIN-ORB", "TENSOR  " };
    char linkId[4];
    char linkFmtOut[9];

    linkId[0] = '*';
    linkId[1] = '1';
    linkId[2] = '0';
    linkId[3] = '0';

    int hasVso, hasVsoi;
    double eLab, tvReal, tvImag, taReal, taImag;
    double kWave, waveLength, eta, H, rMax, temp;
    int verbosity, stepCount, lSkip, statsCode, gridPointCount;
    int localRc, lC1, lC2;
    double dummy;
//
    // caller passed 0 / FALSE_F). The `if (!FITSW)` guard always fired.
    verbosity = reaction.flags.printLevel % 10;
    reaction.distortedWave.scatteringSolver.pwAvSwitch = ((reaction.flags.printLevel / 1000) % 10) >= 1;
    reaction.distortedWave.scatteringSolver.pwBgSwitch = ((reaction.flags.printLevel / 1000) % 10) >= 4;
//
    reaction.distortedWave.scatteringSolver.isStandalone = isStandalone;
//
//     SOME OF THE PARAMETERS CAN HAVE AN ENERGY DEPENDENCE WHICH
//     IS EXPRESSED IN TERMS OF THE LABORATORY ENERGY
//
    eLab = (1 + reaction.internalState.ratMass) * reaction.energies.E;
    tvReal = reaction.opticalPotentialParams.V;
    tvImag = reaction.opticalPotentialParams.vI;
    taReal = reaction.opticalPotentialParams.A;
    taImag = reaction.opticalPotentialParams.aI;
    // Store in potentialWork member for use by MAKPOT and downstream readers.
    reaction.distortedWave.scatteringSolver.potentialWork.tvReal = tvReal;
    reaction.distortedWave.scatteringSolver.potentialWork.tvImag = tvImag;
    reaction.distortedWave.scatteringSolver.potentialWork.taReal = taReal;
    reaction.distortedWave.scatteringSolver.potentialWork.taImag = taImag;
//
//     IS THERE A SPIN-ORBIT POTENTIAL
//
    hasVso = reaction.opticalPotentialParams.vSo != 0;
    hasVsoi = reaction.opticalPotentialParams.vSoi != 0;
    reaction.distortedWave.scatteringSolver.hasAnySpinorbit = hasVso || hasVsoi;
//
//     IS THERE A TENSOR POTENTIAL?
//
    reaction.distortedWave.channel[reaction.internalState.waveChannel].hasSpinorbit = reaction.distortedWave.scatteringSolver.hasAnySpinorbit;
//
//     TEST THE INPUT AND RETURN  0  COMPLETION CODE IF BAD
//
    localRc = 1;
    // Identical post-linkule failure guard repeated after three linkule() calls;
    // returns true (caller then returns) when the call reported localRc < 0.
    auto bailIfCallFailed = [&]() {
        if (localRc < 0) { localRc = 0; return true; }
        return false;
    };
    if (reaction.linkuleData.linkuleAddr[1][3] == 0 && !(reaction.integrationGrid.R > 0 && taReal > 0)) {
        std::printf("0**** R OR A HAS INVALID VALUE:%15.5G%15.5G\n", reaction.integrationGrid.R, taReal);
        localRc = 0;
    }
    if (tvImag != 0 && reaction.linkuleData.linkuleAddr[2][3] == 0 && !(reaction.opticalPotentialParams.rI > 0 && taImag > 0)) {
        std::printf("0**** RI OR AI HAS INVALID VALUE:%15.5G%15.5G\n", reaction.opticalPotentialParams.rI, taImag);
        localRc = 0;
    }
    if (hasVso && reaction.linkuleData.linkuleAddr[3][3] == 0 && !(reaction.opticalPotentialParams.rSo > 0 && reaction.opticalPotentialParams.aSo > 0)) {
        std::printf("0**** RSO OR ASO HAS INVALID VALUE:%15.5G%15.5G\n", reaction.opticalPotentialParams.rSo, reaction.opticalPotentialParams.aSo);
        localRc = 0;
    }
//
//     CHECK  SP  FOR EITHER REAL OR IMAG SPIN ORBIT.
//
    if (reaction.angMom.spinProj == reaction.internalState.notDefSentinel && reaction.distortedWave.channel[reaction.internalState.waveChannel].hasSpinorbit) {
        reaction.angMom.spinProj = 1;
        std::printf("0**** WARNING:  SP WAS NOT DEFINED; "
            "IT IS ASSUMED TO BE %1d/2 FOR THE %.8s FORCE.\n",
            (int)reaction.angMom.spinProj, soTNames[(int)reaction.angMom.spinProj - 1]);
    }
    if (!(reaction.energies.E > 0 && reaction.masses.aM > 0)) {
        std::printf("0**** ENERGY OR M (REDUCED MASS) HAS INVALID"
            " VALUE:%15.5G%15.5G%15.5G\n", reaction.energies.E, reaction.masses.aM, 0.0);
        localRc = 0;
    }
    if (!(reaction.masses.massProj != reaction.internalState.undefValue && reaction.masses.massTgt != reaction.internalState.undefValue)) {
        std::printf("0**** BOTH MP AND MT MUST BE DEFINED.\n");
        localRc = 0;
    }
    if (!(reaction.opticalPotentialParams.rC < reaction.integrationGrid.asymptopia)) {
        std::printf("0**** RC MUST BE LESS THAN ASYMPTOPIA:%15.5G%15.5G\n",
            reaction.opticalPotentialParams.rC, reaction.integrationGrid.asymptopia);
        localRc = 0;
    }
    {
        double rcc = reaction.opticalPotentialParams.rC;
        if (reaction.opticalPotentialParams.rC == reaction.internalState.undefValue) {
            reaction.internalState.wasSet[9] = 1;
            rcc = 1;
        }
        if (reaction.charges.zProj * reaction.charges.zTarget != 0 && !(reaction.opticalPotentialParams.rC > 0)) {
            std::printf("0**** RC IS INVALID:%15.5G\n", reaction.opticalPotentialParams.rC);
            localRc = 0;
        }
        reaction.opticalPotentialParams.rC = rcc;
    }
    if (reaction.opticalPotentialParams.vSi != 0 && reaction.linkuleData.linkuleAddr[2][3] == 0 && !(reaction.opticalPotentialParams.rSi > 0 && reaction.opticalPotentialParams.aSi > 0)) {
        std::printf("0**** RSI OR ASI IS INVALID:%15.5G%15.5G\n", reaction.opticalPotentialParams.rSi, reaction.opticalPotentialParams.aSi);
        localRc = 0;
    }
    if (hasVsoi && reaction.linkuleData.linkuleAddr[4][3] == 0 && !(reaction.opticalPotentialParams.rSoi > 0 && reaction.opticalPotentialParams.aSoi > 0)) {
        std::printf("0**** RSOI OR ASOI HAS INVALID VALUE:%15.5G%15.5G\n",
            reaction.opticalPotentialParams.rSoi, reaction.opticalPotentialParams.aSoi);
    }
//
    if (localRc == 0) return;
//
//
//
//  SOME DERIVED DATA
//
//  K VALUE
//
    kWave = std::sqrt(2 * reaction.masses.aM * reaction.energies.E) / Constants::hbar_c;
//
    waveLength = 2 * Constants::PI / kWave;
//
//  COULOMB PARAMETER AND SIMPLE TERMS USING IT
//     INITIALIZE THE COULOMB POTENTIAL
//
    eta = (reaction.charges.zProj * reaction.charges.zTarget) * std::sqrt(reaction.masses.aM / (2.0 * reaction.energies.E)) / Constants::fine_structure_inv;
    // setVsq is now called AFTER MAKPOT (below) where rcTarget/rcProj have valid values.
    // Calling it here with potentially undefValue values (reset by CLRCHN) would corrupt the arrays.
//
//
    if (reaction.integrationGrid.stepsPerUnit != reaction.internalState.undefValue)
        reaction.integrationGrid.stepSize = std::min(waveLength, 1.0) / reaction.integrationGrid.stepsPerUnit;
//
    // Ensure undefValue has the correct value (may have been corrupted by raw array writes)
    {
        uint64_t undefPattern = 0xF0F0F0F0F0F0F0F0ull;
        double undefValue;
        std::memcpy(&undefValue, &undefPattern, sizeof(double));
        if (reaction.internalState.undefValue != undefValue) reaction.internalState.undefValue = undefValue;
    }
    rMax = reaction.integrationGrid.asymptopia;
//
//     ADJUST ASYMPTOPIA IF THERE WILL BE VERY LARGE L'S.
//     MATCHING MUST BE AT OR OUTSIDE THE CLASSICAL TURNING POINT.
//     HOWEVER, FITTER DOES NOT ALLOW VARIABLE ASYMPTOPIA RIGHT NOW
//     ALSO IF SCATASYM WAS READ IN, WE DO NOT OVERRIDE IT
//
    if (reaction.angMom.lMax != NOTDEF_INT && reaction.integrationGrid.scatAsy < 0)
        rMax = std::max(rMax, (eta + std::sqrt(eta * eta + (double)reaction.angMom.lMax * (reaction.angMom.lMax + 1))) / kWave);
//
//     NOTE - stepCount IS COMPUTED HERE SEPARATLY FOR EACH CHANNEL
//     AND IS STORED IN /WAVCOM/.  IT NEED NOT BE RELATED TO THE
//     "nSteps" USED IN THE BOUND STATE CALCULATION.  HOWEVER, IF
//     asymptopia AND stepSize  ARE THE SAME FOR THE ENTIRE JOB, THEN
//        stepCount = nSteps-1.
//
    stepCount = (int)(rMax / reaction.integrationGrid.stepSize + 0.5);
    rMax = stepCount * reaction.integrationGrid.stepSize;
    H = kWave * reaction.integrationGrid.stepSize;
//
//     TEST FOR IDENTICAL PARTICLES
//
    lSkip = 1;
    statsCode = 3;
    if (reaction.charges.zProj == reaction.charges.zTarget && reaction.masses.massProj == reaction.masses.massTgt
        && reaction.angMom.spinProj == reaction.angMom.spinTarget && reaction.energies.exsPt[1] == reaction.energies.exsPt[2]) {
        // Identical particles
        if (reaction.angMom.spinProj == reaction.internalState.notDefSentinel) {
            std::printf("0*** WARNING:  THE PARTICLES APPEAR TO BE IDENTICAL"
                " BUT THEIR SPINS HAVE NOT BEEN DEFINED;\n"
                "%15s%s\n", "", "NON-IDENTICAL PARTICLE SCATTERING WILL BE ASSUMED.");
        } else {
            if ((int)reaction.angMom.spinProj == 0) lSkip = 2;
            statsCode = ((int)reaction.angMom.spinProj & 1) + 1;
        }
    }
//
//     SAVE THINGS FOR WAVPOT AND WAVELJ
//
    reaction.distortedWave.channel[reaction.internalState.waveChannel].stepSize = H;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].v0R   = tvReal;
    reaction.kin.rScts[reaction.internalState.waveChannel] = reaction.integrationGrid.R;
    reaction.kin.aScts[reaction.internalState.waveChannel] = taReal;
    // rcScts reader is channel_clear's inelastic L800 block — channel[1] only.
    if (reaction.internalState.waveChannel == 1)
        reaction.kin.rcScts[1] = reaction.opticalPotentialParams.rC;
    reaction.kin.rcSctP[reaction.internalState.waveChannel] = reaction.masses.rcProj;
    reaction.kin.rcSctT[reaction.internalState.waveChannel] = reaction.masses.rcTarget;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].aI = taImag;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].aSi = reaction.opticalPotentialParams.aSi;
    // v0I/v0Si/rI/RSI restored to channel[1] only — sole reader is the
    // inelastic L800 block in channel_clear which restores from [1].
    if (reaction.internalState.waveChannel == 1) {
        reaction.distortedWave.channel[1].v0I  = tvImag;
        reaction.distortedWave.channel[1].v0Si = reaction.opticalPotentialParams.vSi;
        reaction.distortedWave.channel[1].rI   = reaction.opticalPotentialParams.rI;
        reaction.distortedWave.channel[1].rSi  = reaction.opticalPotentialParams.rSi;
    }
    reaction.kin.etaCh[reaction.internalState.waveChannel] = eta;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].asymptopia   = rMax;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].Ecm          = reaction.energies.E;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].rStart      = reaction.integrationGrid.stepSize;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].nGridSteps = stepCount;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].nStp2s       = stepCount;
    if (reaction.internalState.waveChannel == 1) {
        reaction.kin.akIn   = kWave;
        reaction.kin.redMi = reaction.masses.aM;
    } else {
        reaction.kin.akOut   = kWave;
        reaction.kin.redMo = reaction.masses.aM;
    }
    reaction.kin.tauRatio[reaction.internalState.waveChannel] = reaction.internalState.ratMass;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].twoSpin = (int)reaction.angMom.spinProj;
    // STEP1R / STEP1I were ReactionParams fields bulk-defaulted to 1.0 by
    reaction.distortedWave.scatteringSolver.stepI = 1.0;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].lSkips = lSkip;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].statsCode = statsCode;
    // wavR_arr/wavI_arr are reset by .assign() in angular_setup when TCSWS is set.
//
//     QUANTITIES FOR INTEGRATION PAST ASYMPTOPIA
//
    reaction.distortedWave.channel[reaction.internalState.waveChannel].xFacs[1] = 1 + H * H / 12;
    reaction.distortedWave.channel[reaction.internalState.waveChannel].xFacs[2] = -(H * H / 12) * reaction.charges.zProj * reaction.charges.zTarget * Constants::hbar_c / (reaction.energies.E * Constants::fine_structure_inv);
    reaction.distortedWave.channel[reaction.internalState.waveChannel].xFacs[3] = -(H * H / 12) / (kWave * kWave);
    reaction.timing.times[1] = 0;
    reaction.timing.times[2] = 0;
    reaction.timing.times[3] = 0;
//
//     INITIALIZE LINKULES IF NEEDED
//
    for (int i = 1; i <= numLinkules; i++) {
        if (reaction.linkuleData.linkuleAddr[i][3] != 0) {
//
//     GENERATE THE SPECIAL UNIQUE NAME
//
            reaction.linkuleData.uniqueLinkuleId = reaction.linkuleData.uniqueLinkuleId + 1;
            std::snprintf(linkFmtOut, sizeof(linkFmtOut), "*%03d", reaction.linkuleData.uniqueLinkuleId);
            std::memcpy(linkId, linkFmtOut, 4);
//
//     MAKE INITIALIZING CALL
//
            linkule(reaction.linkuleData.linkuleAddr[i][3], *(char8*)&reaction.linkuleData.linkuleAddr[i][1], &reaction.linkuleData.linkuleAddr[i][5], i, 1,
                localRc, reaction.angMom.L, reaction.angMom.J, 0.0, reaction.integrationGrid.stepSize, stepCount + 1, &dummy, &dummy, linkId, reaction);
            if (bailIfCallFailed()) return;
            reaction.linkuleData.linkuleAddr[i][4] = localRc;
        }
    } // 259
//
    if (verbosity >= 1) {
//
//     PRINT OUT INPUT NICELY
//
    std::printf("\n0        OPTICAL MODEL SCATTERING FOR THE %.8s CHANNEL\n",
        channelName[reaction.internalState.waveChannel - 1]);
    std::printf("0E LAB =%9.3f MEV,    E CM =%9.3f MEV,     K =", eLab, reaction.energies.E);
    print_G(12, 5, kWave);
    std::printf("     WAVELENGTH =%8.4f FM\n"
        " PROJECTILE MASS =%7.2f AMU,%5sTARGET MASS =%7.2f AMU,%5sREDUCED MASS =%10.2f MEV/C**2\n",
        waveLength, reaction.masses.massProj, "", reaction.masses.massTgt, "", reaction.masses.aM);
    if (statsCode == 3) std::printf(" THIS IS NON-IDENTICAL PARTICLE SCATTERING\n");
    if (statsCode != 3) std::printf(" THIS IS SCATTERING OF IDENTICAL %.8s\n", statNames[statsCode - 1]);
    if (reaction.angMom.spinProj != reaction.internalState.notDefSentinel && reaction.angMom.spinTarget != reaction.internalState.notDefSentinel)
        std::printf(" PROJECTILE SPIN =%3d/2%5sTARGET SPIN =%3d/2\n", (int)reaction.angMom.spinProj, "", (int)reaction.angMom.spinTarget);
    std::printf(" Z PROJECTILE =%4d%5sZ TARGET =%4d\n"
        "0POTENTIAL         COUPLING CONS.    RADIUS    DIFFUSENESS    RADIUS PARAMETER\n",
        reaction.charges.zProj, "", reaction.charges.zTarget);
//
    localRc = 1;
    {
        double r0M = reaction.internalState.r0Mass;
        if (r0M == reaction.internalState.undefValue || r0M == 0.0) r0M = 1.0;
        parameterPrint(1,  "REAL CENTRAL       ", tvReal,              reaction.integrationGrid.R,    taReal,             reaction.integrationGrid.R    / r0M, reaction);
        parameterPrint(2,  "VOLUME ABSORPTION  ", tvImag,              reaction.opticalPotentialParams.rI,   taImag,             reaction.opticalPotentialParams.rI   / r0M, reaction);
        parameterPrint(13, "SURFACE ABSORPTION ", reaction.opticalPotentialParams.vSi,   reaction.opticalPotentialParams.rSi,  reaction.opticalPotentialParams.aSi,  reaction.opticalPotentialParams.rSi  / r0M, reaction);
        parameterPrint(3,  "REAL SPIN-ORBIT    ", reaction.opticalPotentialParams.vSo,   reaction.opticalPotentialParams.rSo,  reaction.opticalPotentialParams.aSo,  reaction.opticalPotentialParams.rSo  / r0M, reaction);
        parameterPrint(4,  "IMAG. SPIN-ORBIT   ", reaction.opticalPotentialParams.vSoi,  reaction.opticalPotentialParams.rSoi, reaction.opticalPotentialParams.aSoi, reaction.opticalPotentialParams.rSoi / r0M, reaction);
        parameterPrint(5,  "PT&SPHERE  COULOMB ", eta,                reaction.opticalPotentialParams.rC,   0.0,               reaction.opticalPotentialParams.rC   / r0M, reaction);
    }

    if (!((reaction.masses.rcProj == 0) || (reaction.masses.rcTarget == 0) ||
        (reaction.masses.rcProj == reaction.internalState.undefValue) || (reaction.masses.rcTarget == reaction.internalState.undefValue))) {
        std::printf(" FOLDED COULOMB POTENTIALS - RCP = %7.4f  RCT = %7.4f  RC0P = %7.4f  RC0T = %7.4f\n",
            reaction.masses.rcProj, reaction.masses.rcTarget,
            reaction.masses.rcProj / std::pow(reaction.masses.massProj, 0.333333333333333333333330),
            reaction.masses.rcTarget / std::pow(reaction.masses.massTgt, 0.333333333333333333333330));
    }
//
    temp = std::min(waveLength, 1.0) / reaction.integrationGrid.stepSize;
    std::printf("0ASYMPTOPIA = %8.3f FM\n"
        " STEP SIZE =  %8.3f FM     %8.3f STEPS PER \"WAVELENGTH\"\n\n",
        rMax, reaction.integrationGrid.stepSize, temp);
//
//     WILL THE WAVE FUNCTIONS BE FOUND BY A LINKULE
//
    if (reaction.linkuleData.linkuleAddr[6][3] != 0) {
        std::printf("0 THE SCATTERING WAVEFUNCTION IS BEING"
            " COMPUTED BY THE %.8s LINKULE:\n", (char*)&reaction.linkuleData.linkuleAddr[6][1]);
        {
            int i = 6;
//
//
            linkule(reaction.linkuleData.linkuleAddr[i][3], *(char8*)&reaction.linkuleData.linkuleAddr[i][1], &reaction.linkuleData.linkuleAddr[i][5],
                i, 2, localRc,
                reaction.angMom.L, reaction.angMom.J, 0.0, reaction.integrationGrid.stepSize, stepCount + 1, &dummy, &dummy, (char*)&reaction.internalState.waveChannel, reaction);
            std::printf(" \n");
            if (bailIfCallFailed()) return;
        }
    }
    } // end if (verbosity >= 1)
//
//
//
//     SET ASIDE AREAS FOR THE POTENTIAL AND FOR V' USED BY WAVELJ
//     THESE ARE STORED ONLY TO ASYMPTOPIA
//
//     THIS IS DONE AS FOLLOWS:
//
//     2) FOR ELASTIC FITS: MAXIMUM NUMBER OF STEPS IS ACCUMULATED IN
//        nGridSteps(2); AS THE SIZE GROWS, THEY ARE ALLOCATED.
//     3) FOR DWBA: POTENTIAL AREAS ARE ALLOCATED AT BOTH CALLS;
//        WAVEFUNCTION AREAS ARE ALLOCATED ONLY ON 2ND CALL.
//
    gridPointCount = stepCount + 1;
//
    {
//
    reaction.distortedWave.scatteringSolver.allocateVWork(gridPointCount);
    {
        auto& dwc = reaction.distortedWave.channel[reaction.internalState.waveChannel];
        dwc.rlvsArr.assign(gridPointCount, 0.0);
        dwc.imvsArr.assign(gridPointCount + 1, 0.0);  // +1 guard cell (legacy)
        dwc.centrArr.assign(gridPointCount, 0.0);
        dwc.rlvPointer  = dwc.rlvsArr.data() - 1;   // 1-based ptr
        dwc.imvPointer  = dwc.imvsArr.data() - 1;
        dwc.centPointer = dwc.centrArr.data() - 1;
    }
    // solvePartialWave guards on soRPointer / soIPointer != nullptr.
    {
        int channelIndex = reaction.internalState.waveChannel;
        reaction.distortedWave.channel[channelIndex].soRPointer = nullptr;
        reaction.distortedWave.channel[channelIndex].soIPointer = nullptr;
        if (hasVso)
            reaction.distortedWave.scatteringSolver.allocateSors(channelIndex, gridPointCount, reaction);
        if (hasVsoi)
            reaction.distortedWave.scatteringSolver.allocateSois(channelIndex, gridPointCount, reaction);
    }
    } // end pool alloc else
//
//     AND COMPUTE THE POTENTIAL
//
    // Initialize vcsq12 Coulomb coefficients for this channel BEFORE MAKPOT.
    // reaction.opticalPotentialParams.rC was computed by SETPOT. rcTarget/rcProj may be undefValue if SETFIT is not
    // fully translated, so fall back to rcTarget=rC (full Coulomb radius), rcProj=0 (point projectile).
    {
        double rcTarget = reaction.masses.rcTarget;
        double rcProj = reaction.masses.rcProj;
        if (!std::isfinite(rcTarget) || rcTarget <= 0.0 ||
            !std::isfinite(rcProj) || rcProj < 0.0) {
            rcTarget = reaction.opticalPotentialParams.rC;
            rcProj = 0.0;
        }
        setVsq(rcTarget, rcProj, reaction.charges.zTarget, reaction.charges.zProj, reaction.internalState.waveChannel);
        // Also fix up rcSctP/rcSctT which were saved above with undefValue values
        if (!std::isfinite(reaction.kin.rcSctP[reaction.internalState.waveChannel]) || reaction.kin.rcSctP[reaction.internalState.waveChannel] < 0.0)
            reaction.kin.rcSctP[reaction.internalState.waveChannel] = rcProj;
        if (!std::isfinite(reaction.kin.rcSctT[reaction.internalState.waveChannel]) || reaction.kin.rcSctT[reaction.internalState.waveChannel] <= 0.0)
            reaction.kin.rcSctT[reaction.internalState.waveChannel] = rcTarget;
    }
    reaction.makePotential(reaction.internalState.waveChannel, localRc);
//
//     FOR DWBA MUST SAVE TWO SETS OF LINKULE INDICATORS
//
    for (int i = 1; i <= numLinkules; i++) {
        for (int ii = 1; ii <= 6; ii++) {
            reaction.linkuleData.lnkAd2[i][reaction.internalState.waveChannel][ii] = reaction.linkuleData.linkuleAddr[i][ii];
        }
    }
//
//     WORK AREA NOT NEEDED AFTER MAKPOT
//
    // vWork vector lifetime owned by ScatteringSolver — pool free no longer needed
//
    if (bailIfCallFailed()) return;
//
//     DETERMINE LCRITICAL.
//     IF BOTH lMin AND lMax ARE DEFINED WE BYPASS THIS AND JUST USE
//     THE AVERAGE.
//
    if (reaction.angMom.lMin == NOTDEF_INT || reaction.angMom.lMax == NOTDEF_INT) {
        lCritL(kWave, eta, reaction.opticalPotentialParams.rC, stepCount + 1,
            reaction.integrationGrid.stepSize,
            reaction.distortedWave.channel[reaction.internalState.waveChannel].rlvsArr.data(),
            reaction.distortedWave.channel[reaction.internalState.waveChannel].imvsArr.data(),
            (reaction.flags.printLevel / 10) % 10,
            lC1, lC2, reaction.kin.lCrits[reaction.internalState.waveChannel]);
        if (verbosity != 0)
            std::printf(" ESTIMATED CRITICAL L:%4d%5sBY |S|:%4d%5sBY DEFLECTION FUNC.:%4d\n",
                reaction.kin.lCrits[reaction.internalState.waveChannel], "", lC1, "", lC2);
        if (lC1 == 0 && reaction.angMom.lMin == NOTDEF_INT) reaction.angMom.lMin = 0;
    } else {
        reaction.kin.lCrits[reaction.internalState.waveChannel] = (reaction.angMom.lMin + reaction.angMom.lMax) / 2;
    }
//
//
    if (reaction.distortedWave.scatteringSolver.isStandalone || reaction.internalState.waveChannel != 1) {
        if (!reaction.distortedWave.scatteringSolver.isStandalone) gridPointCount = std::max(reaction.distortedWave.channel[1].nGridSteps + 1, gridPointCount);
        {
            int i = gridPointCount + 6;
            reaction.distortedWave.scatteringSolver.nWaveF = i;
//
            reaction.distortedWave.scatteringSolver.reallocate(i);
        }
    }
    localRc = 1;
    return;
//
//
//
}
