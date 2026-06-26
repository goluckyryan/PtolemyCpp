// Elastic.cpp — Elastic scattering: Elastic::calculate (L-loop + S-matrix + DCS)
// and dispatchPartialWave (forwards to solvePartialWave).

#include "Elastic.h"

// ============================================================================
// SHARED INCLUDES (deduplicated)
// ============================================================================

#include "ptolemy_types.h"
#include "CrossSectionCalc.h"
#include "ScatteringSolver.h"
#include <cstdio>
#include <vector>
#include "Reaction.h"
#include "Constants.h"

Elastic::Elastic(Reaction& reaction) : reaction_(reaction) {}


// ============================================================================
// SECTION 1: wavefunction_io.cpp — calculate, dispatchPartialWave
// ============================================================================

// ============================================================================
//
//  calculate CALCULATES RADIAL WAVEFUNCTION FOR ALL PARTIAL WAVES
//  USING WOOD-SAXON POTENTIAL
//
//     THIS ROUTINE IS USED ONLY FOR STAND-ALONE CALCULATIONS OF
//     ELASTIC SCATTERING.  IT IS AN INTERFACE TO WAVSET, setupWavefunctionPotential
//     AND WAVELJ WITH A COUPLE OF BELLS AND WHISTLES ADDED.
//
//     IF L IS DEFINED THEN THE CALC IS DONE FOR ONLY ONE VALUE OF
//     L; ELSE THE RANGE (lMin, lMax) IS USED.
//     IF THERE IS A SPIN-ORBIT FORCE THEN jProj IS CHECKED.  IF IT IS
//     DEFINED THEN ONLY THAT jProj (AND THE ASSOCIATED L'S) ARE DONE.
//     IF IT IS NOT DEFINED THEN ALL jProj'S FOR EACH L ARE DONE.
//
// ============================================================================
bool Elastic::calculate()
{
    int returnCode;

    // implicit real*8 ( a-h, o-z )
    // Local variables (i-N are integer, A-H,O-Z are double)

        // All from globals (shared mutable state)
    auto& L      = reaction_.angMom.L;
    auto& lMax   = reaction_.angMom.lMax;
    auto& lMin   = reaction_.angMom.lMin;
    int lx;
    auto& printLevel = reaction_.flags.printLevel;
    auto& jProj     = reaction_.angMom.jProj;
    auto& spinProj    = reaction_.angMom.spinProj;


    auto& waveChannel = reaction_.internalState.waveChannel;

    auto& angleMin = reaction_.rxn.angleMin;
    auto& angleMax = reaction_.rxn.angleMax;
    auto& angleStep = reaction_.rxn.angleStep;
    auto& eLab   = reaction_.energies.eLab;

    // DEGREE was only used by the deleted IRLWAV rotation block.

    int notDefSentinel = NOTDEF_INT;

    // LOGICAL locals
    int hasJProj, hasSpinorbit, keepFAmplitude;

    // Other locals
    int verbosity, llMax, llMin, lSave;
    double jpSave;  // must be double: preserves notDefSentinel bit pattern (int cast loses it)
    int lSkip, statsCode, nSpline;
    // WRITES / ICHECK blocks.
    // asymptotic-check block (ICHECK permanently 0).
    // pointer aliases (declared at function scope to avoid jump-over-init)
    double* smatsPointer = nullptr;
    int*    tocePointer   = nullptr;
    int lines, twoMsStart, twoMsEnd, twoMs;
    int i, kOffset, lL;
    // TORUT/ruth/CROSS/F_amp from elasticDcs are caller-owned std::vector.
    // Other vectors are write-only by this caller.
    std::vector<double> torutVector, ruthVector, crossSectionVector, fAmpVector;
    double tL, sReal, sImag, temp, phase;
    double sigReaction;

    //
    verbosity = ((printLevel) % (10));
    if (verbosity > 0) {
        std::printf("1%47s%s\n", "", "P T O L E M Y");
        std::printf("0");
        for (int k = 1; k <= 45; k++) std::printf("%c", reaction_.reactStr[k]);
        std::printf("%8.2f MEV     ", eLab);
        for (int k = 1; k <= 65; k++) std::printf("%c", reaction_.header[k]);
        std::printf("\n");
    }

    //
    //     TEST THE INPUT AND RETURN  0  COMPLETION CODE IF BAD
    //
    // sets returnCode = 1 at its first executable line and writes it on every
    // path; no reader observes the init.
    llMax = (int)lMax;
    llMin = (int)lMin;
    lSave = L;
    jpSave = jProj;  // save as double to preserve notDefSentinel bit pattern

    //
    if (L != notDefSentinel) {
        llMax = L;
        llMin = L;
    }
    if (waveChannel > 2) waveChannel = 1;

    //
    //     SET UP ARRAYS AND PRINT PARAMETERS.
    //
    reaction_.distortedWave.scatteringSolver.setupScatteringWaves(returnCode, TRUE_F, reaction_);
    if (returnCode == 0) { L = lSave; jProj = jpSave; return false; }
    // overwrites returnCode unconditionally and no reader between here and there.
    lSkip = reaction_.distortedWave.channel[waveChannel].lSkips;
    statsCode = reaction_.distortedWave.channel[waveChannel].statsCode;

    //
    if (llMin != notDefSentinel && llMax != notDefSentinel) {
        if (lSkip == 2) llMin = std::abs(llMin - ((llMin) % (2)));
        if (lSkip == 2) llMax = llMax + ((llMax) % (2));
    } else {
        reaction_.kin.lCrit = reaction_.kin.lCrits[waveChannel];
        if (llMin == notDefSentinel) {
            llMin = (int)(reaction_.kin.lCrit * reaction_.opticalPotentialParams.alMnMt);
            llMin = std::min(llMin, reaction_.kin.lCrit - reaction_.integrationGrid.lMinSub);
            llMin = std::max(llMin, 0);
            if (lSkip == 2) llMin = std::abs(llMin - ((llMin) % (2)));
        }
        if (llMax == notDefSentinel)
            llMax = std::max(reaction_.kin.lCrit + reaction_.angMom.lMaxAdditional, (int)(reaction_.opticalPotentialParams.alMxMt * reaction_.kin.lCrit));
        if (lSkip == 2) llMax = llMax + ((llMax) % (2));
    }
    std::printf("0%10d =< L =<%4d\n", llMin, llMax);

    //     ARE WE DOING ONLY ONE jProj
    //
    hasSpinorbit = reaction_.distortedWave.channel[waveChannel].hasSpinorbit;
    hasJProj = (jProj != reaction_.internalState.notDefSentinel);
    twoMsStart = 0;
    twoMsEnd = 0;
    if (ftobool(hasSpinorbit)) {
        if (!ftobool(hasJProj)) {
            twoMsStart = -(int)spinProj;
            twoMsEnd = (int)spinProj;
        } else {
            if ((((int)jProj + (int)spinProj) % (2)) != 0) {
                std::printf("0**** JP AND SP ARE AN INVALID HALF-INTEGER AND "
                            "INTEGER COMBINATION:%6d/2%6d/2\n", (int)jProj, (int)spinProj);
                L = lSave; jProj = jpSave; return false;
            }
            llMin = std::max(llMin, ((int)jProj - (int)spinProj) / 2);
            llMax = std::min(llMax, ((int)jProj + (int)spinProj) / 2);
            if (llMin > llMax) {
                std::printf("0**** THERE ARE NO VALID COMBINATIONS OF JP, SP AND L:"
                            "   JP, SP =%5d/2%5d/2\n LMIN, LMAX =%5d%5d  OR  L =%5d\n",
                            (int)jProj, (int)spinProj, lMin, lMax, L);
                L = lSave; jProj = jpSave; return false;
            }
        }
    }
    reaction_.kin.lOutMax = llMax;

    //
    //     SET UP S-MATRIX AND COULOMB WF ARRAYS.
    //
    if (!reaction_.setupWavefunctionPotential()) { L = lSave; jProj = jpSave; return false; }
    nSpline = reaction_.distortedWave.channel[waveChannel].nJStates;



    //
    //     ALL ALLOCATIONS DONE, GET LOCATIONS
    //
    // smatArr class-owned vector. smatsPointer is a 0-based ptr.
    // Accessed as [i-1]/[i] (shifted from the former 1-based [i]/[i+1]).
    smatsPointer = reaction_.distortedWave.channel[waveChannel].smatArr.data();
    // toceArr replaces pool slot; tocePointer is 0-based (accessed [4*kOffset-3]/[4*kOffset-4],
    // matching the writer's [4*ii-4 .. 4*ii-1] layout in source_potentials).
    tocePointer   = reaction_.distortedWave.channel[waveChannel].toceArr.data();
    std::fill(reaction_.distortedWave.channel[waveChannel].smatArr.begin(),
              reaction_.distortedWave.channel[waveChannel].smatArr.end(), 0.0);

    //
    if (llMin == llMax) verbosity = 1;

    //
    //     PRINT HEADINGS IF NECESSARY.
    //
    lines = 35;
    if (verbosity > 0) {
        std::printf("\n0 L    L'   LX%19sS%24s|S|%10sPHASE SHIFT%11sFRACTION\n",
                    "", "", "", "");
        std::printf("%70sDEGREES%11sABSORBED\n", "", "");
        std::printf(" \n");
    }

    //
    ;

    // DUMMY1() (Fortran overlay-loader hint for RACAH tables) deleted

    //
    //     LOOP THROUGH INCIDENT L.
    //
    for (L = llMin; L <= llMax; L += lSkip) {

        //
        //        LOOP THROUGH jProj = 2*J(PROJECTILE).  SPIN OF TARGET IS IGNORED.
        //        S-MATRIX ELEMENTS FOR ALL lx ARE ACCUMULATED AS SUMS OVER jProj.
        //
        for (twoMs = twoMsStart; twoMs <= twoMsEnd; twoMs += 2) {
            if (!ftobool(hasJProj)) jProj = 2 * L + twoMs;

            //
            //           (TENSOR COUPLING) TO DO THE ACTUAL CALCULATION.
            //
            // nPts=0 → solvePartialWave/dispatchPartialWave never touch
            // rGrid/WAVER/WAVEI; pass nullptr instead of pool sentinels.
            dispatchPartialWave(reaction_, L, (int)jProj, waveChannel, 0, nullptr, nullptr, nullptr,
                  reaction_.distortedWave.scatteringSolver.wavRPointer, reaction_.distortedWave.scatteringSolver.wavIPointer,   // WAVCOM pointer fields
                  reaction_.distortedWave.channel[waveChannel].rlvsArr.data() - 1,
                  reaction_.distortedWave.channel[waveChannel].imvsArr.data() - 1,
                  reaction_.distortedWave.channel[waveChannel].centrArr.data() - 1);

            //
        } // end twoMs loop (259)

        //
        //        END OF CALCULATION LOOP THROUGH jProj.
        //        IF NECESSARY, RETRIEVE AND PRINT THE S-MATRIX ELEMENTS.
        //
        if (verbosity != 0) {

        //
        //        PRINT HEADINGS IF NECESSARY.
        //
        if (lines + nSpline >= 59) {
        std::printf("1%47s%s\n", "", "P T O L E M Y");
        std::printf("0");
        for (int k = 1; k <= 45; k++) std::printf("%c", reaction_.reactStr[k]);
        std::printf("%8.2f MEV     ", eLab);
        for (int k = 1; k <= 65; k++) std::printf("%c", reaction_.header[k]);
        std::printf("\n");
        std::printf("\n0 L    L'   LX%19sS%24s|S|%10sPHASE SHIFT%11sFRACTION\n",
                    "", "", "", "");
        std::printf("%70sDEGREES%11sABSORBED\n", "", "");
        std::printf(" \n");
        lines = 7;

        //
        //        CALCULATE THE FRACTION ABSORBED FROM THIS L.
        } // end page header
        tL = 1.0;
        i = 2 * nSpline * L - 1;
        for (kOffset = 1; kOffset <= nSpline; kOffset++) {
            i = i + 2;
            tL = tL - smatsPointer[i-1] * smatsPointer[i-1] - smatsPointer[i] * smatsPointer[i];
        } // end kOffset loop (279)

        //
        //        LOOP THROUGH kOffset = lx, LAS.  QUANTUM NUMBERS ARE
        //        IN THE TABLE-OF CONTENTS ARRAY.
        //
        for (kOffset = 1; kOffset <= nSpline; kOffset++) {
            lx = tocePointer[4 * kOffset - 3];
            lL = L + tocePointer[4 * kOffset - 4];
            if (lx > L + lL) continue;

            //
            //           RETRIEVE THE S-MATRIX ELEMENT, AND CALCULATE THE
            //           MAGNITUDE, PHASE, AND FRACTION ABSORBED.
            //
            i = 2 * (nSpline * L + kOffset) - 1;
            sReal = smatsPointer[i-1];
            sImag = smatsPointer[i];
            temp = std::sqrt(sReal * sReal + sImag * sImag);
            if (temp == 0) continue;
            phase = 0.5 * Constants::RADIAN * std::atan2(sImag, sReal);

            //
            //           PRINT IT.
            //
            if (kOffset == 1)
                std::printf(" %3d%5d%5d%17.5G +%12.5G I%17.5G%12.2f%16.5G\n",
                            L, lL, lx, sReal, sImag, temp, phase, tL);
            else
                std::printf(" %3d%5d%5d%17.5G +%12.5G I%17.5G%12.2f\n",
                            L, lL, lx, sReal, sImag, temp, phase);
            lines = lines + 1;
        } // end kOffset loop (299)

        //
        //        END OF kOffset LOOP FOR PRINTING.  INSERT A BLANK LINE IF NEEDED.
        //
        if (nSpline != 1) {
            std::printf(" \n");
            lines = lines + 1;
        }
        } // end if (verbosity != 0)
    } // end L loop (399)

    //
    //     END OF LOOP OVER ALL L'S AND J'S
    //
    L = llMax;

    // FOR SINGLE WAVEFUNCTION (LLMAX==LLMIN, non-tensor-coupled): the
    // Fortran path here called IREDEF on the pool slots to shrink the
    // wavefn buffers to nStep+1. With wavR/wavI as std::vectors owned
    // by ScatteringSolver, pool_wavr/pool_wavi are permanent-0 sentinels
    // and the IREDEF would have been a no-op (or UB) on a 0 handle.

    //     IRLWAV "make scattering wave real" rotation block deleted
    //     was only ever consumed by the dead ICHECK asymptotic-check
    //     block; both locals dropped with that block.
    //
    // ICHECK asymptotic-check + WRITES wavefunction-print blocks dropped

    //
    //     IF DESIRED, COMPUTE THE ELASTIC DIFFERENTIAL CROSS SECTIONS
    //
    if (reaction_.flags.isElastic != 0) {
    keepFAmplitude = ftobool(hasSpinorbit);
    CrossSectionCalc(reaction_).elasticDcs(reaction_.distortedWave.channel[waveChannel].stepSize/reaction_.integrationGrid.stepSize,
          reaction_.kin.etaCh[waveChannel], angleMin, angleMax, angleStep,
          llMin, llMax, lSkip, statsCode, (int)spinProj,
          reaction_.distortedWave.channel[waveChannel].smatArr.data(),
          reaction_.distortedWave.channel[waveChannel].toceArr.data(),
          nSpline, reaction_.distortedWave.channel[waveChannel].sigmaArr.data(), eLab, eLab/reaction_.energies.E - 1.0, TRUE_F,
          torutVector, keepFAmplitude,
          fAmpVector, sigReaction,
          ruthVector, crossSectionVector);

    //
    //     NOW COMPUTE THE ANALYZING POWERS (IF NONZERO).
    //
    if (ftobool(hasSpinorbit)) {
    reaction_.distortedWave.scatteringSolver.pwBgSwitch = (((printLevel) % (10)) >= 3);
    CrossSectionCalc(reaction_).analyzingPower(angleMin, angleMax, angleStep, (int)spinProj, (int)spinProj, (int)reaction_.angMom.spinTarget, (int)reaction_.angMom.spinTarget,
           nSpline, 1, reaction_.distortedWave.scatteringSolver.pwBgSwitch,
           eLab, "ELASTIC ",
           fAmpVector, reaction_.distortedWave.channel[waveChannel].toceArr.data() - 1);

    //
    //     A L L   D O N E
    //
    } // end if (hasSpinorbit) analyzing powers
    } // end if (reaction_.isElastic != 0)

    L = lSave;
    jProj = jpSave;
    return true;
}


// ============================================================================
//     LAS parameter dropped — sole consumer was the tensor-coupled gate
// ============================================================================
void dispatchPartialWave(Reaction& reaction, int L, int jProj, int channelIndex, int nPts,
           float* rGrid, float* waveR, float* waveI,
           double* waveReal, double* waveImag, double* vReal, double* vImag,
           double* vCent)
{
    // implicit real*8 ( a-h, o-z )

    //
    if (reaction.distortedWave.scatteringSolver.pwBgSwitch)
        std::printf(" dispatchPartialWave:  L,JP,NWP,NUMPTS=%8d%8d%8d%8d\n",
                    L, jProj, channelIndex, nPts);

    //
    //     CHECK FOR LEGALITY AND COUPLING.
    //
    if (L < 0) {
        // Illegal L
        if (nPts != 0) for (int i = 1; i <= nPts; i++) { waveR[i] = 0; waveI[i] = 0; }
        return;
    }

    if (ftobool(reaction.distortedWave.channel[channelIndex].hasSpinorbit)) {
        if (jProj < std::abs(2 * L - reaction.distortedWave.channel[channelIndex].twoSpin)) {
            // Invalid jProj — set to zero
            if (nPts != 0) for (int i = 1; i <= nPts; i++) { waveR[i] = 0; waveI[i] = 0; }
            return;
        }
    }

    reaction.distortedWave.scatteringSolver.solvePartialWave(L, jProj, channelIndex, nPts, rGrid, waveR, waveI,
            waveReal, waveImag, vReal, vImag, vCent, reaction);
    return;
}

