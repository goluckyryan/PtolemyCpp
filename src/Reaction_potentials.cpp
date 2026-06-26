// Reaction_potentials.cpp — Reaction optical/bound-state potential assembly:
// makePotential, setupOpticalPotential, setupWavefunctionPotential,
// setupInelasticAngMomTable, and the sFromI S-matrix renormalization helper.

#include "ptolemy_types.h"
#include "CoulombWaveFunction.h"
#include "Reaction_potentials.h"
#include "linkule.h"
#include "math/angular_momentum_coeff.h"
#include "math/coulomb_sigma.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include "Reaction.h"
#include "Constants.h"
#include "OpticalPotential.h"

// The legacy Potential hierarchy (WoodsSaxonPotential/SpinOrbitPotential::fill


// ============================================================================
// MAKPOT — Computes potentials for WAVELJ
// ============================================================================
// makePotential — assemble the optical potential for one channel directly from
// flat OpticalPotential grids (Phase E: the virtual Potential hierarchy and the
// CompositeOpticalPotential walk/assemble factory are gone). Each contribution is
// a plain OpticalPotential whose fill* owns the bit-identical woodsX/Coulomb math;
// the 5-target fan-out into the channel buffers is explicit below. This is the
// only reader of reaction.opticalPotentialParams.* for the assembly.
void Reaction::makePotential(int channelIndexIn, int& returnCode)
{
    Reaction& reaction = *this;
    int channelIndex = (channelIndexIn == 3) ? 1 : channelIndexIn;
    returnCode = 0;
    double rStart = reaction.distortedWave.channel[channelIndex].rStart;
    double eInv = 1.0 / reaction.energies.E;
    double h = reaction.distortedWave.channel[channelIndex].stepSize;
    double h2 = h * h;
    int gridPointCount = reaction.distortedWave.channel[channelIndex].nGridSteps + 1;
    // cached 1-based pointers set by setupScatteringWaves.
    double* realVPointer = reaction.distortedWave.channel[channelIndex].rlvPointer;
    double* imagVPointer = reaction.distortedWave.channel[channelIndex].imvPointer;
    double* centrPointer = reaction.distortedWave.channel[channelIndex].centPointer;
    double* soRPointer = reaction.distortedWave.channel[channelIndex].soRPointer;
    double* soIPointer = reaction.distortedWave.channel[channelIndex].soIPointer;
    auto& pw = reaction.distortedWave.scatteringSolver.potentialWork;
    auto& params = reaction.opticalPotentialParams;  // user-input optical params

    // Indicate WAVELJ must reestablish starting point
    reaction.distortedWave.channel[channelIndex].lastL = 999999;
    reaction.distortedWave.scatteringSolver.lastZero = 1;

    // Flat-class assembly (Phase D). The virtual Potential hierarchy + composite
    // walk is gone: each contribution is a plain OpticalPotential whose fill*
    // method owns the (bit-identical) woodsX/Coulomb math. The 5-target fan-out
    // is now explicit per channel buffer below. The five fill grids share the
    // grid spec (origin 0.0, step = the channel rStart) — same args the old
    // fill(out, nPts, 0.0, rStart) calls used.
    //
    // OPTICAL-FITTER LINKULE PATH (out of scope — stays LinkulePlugin): when a
    // fitter is registered for a slot, fill that grid via linkule() run() (rc 3)
    // instead of the built-in fill*. fillLinkule marshals the slot's name/ints
    // into the run() call writing the grid's 1-based buffer; returns false on a
    // linkule error (caller returns early). Coulomb (slot 5) additionally needs
    // channelIndexIn != 3. The fitter writes a slot's own grid, so Coulomb and
    // the real central stay in SEPARATE grids (added below) exactly as before —
    // a fitter must not clobber the Coulomb grid.
    auto fillLinkule = [&](int slot, double* buf) -> bool {
        char8 lname;
        std::memcpy(lname.data, reinterpret_cast<char*>(&reaction.linkuleData.linkuleAddr[slot][1]), 8);
        int myInts[3];
        myInts[1] = reaction.linkuleData.linkuleAddr[slot][5];
        myInts[2] = reaction.linkuleData.linkuleAddr[slot][6];
        double dummy = 0.0;
        linkule(reaction.linkuleData.linkuleAddr[slot][3], lname, myInts, slot, 3, returnCode,
                0, dummy, 0.0, rStart, gridPointCount,
                buf, &dummy, (char*)"", reaction);
        return returnCode >= 0;
    };
    auto hasFitter = [&](int slot) { return reaction.linkuleData.linkuleAddr[slot][3] != 0; };

    int i;

    // --- Coulomb (slot 5) + real central (slot 1) -> realV, centr ---
    OpticalPotential coul;   coul.resize(gridPointCount, 0.0, rStart);
    if (reaction.linkuleData.linkuleAddr[5][3] > 0 && channelIndexIn != 3) {
        if (!fillLinkule(5, coul.data1Based())) return;
    } else {
        coul.fillCoulomb(channelIndex, channelIndexIn);
    }
    OpticalPotential realCentral; realCentral.resize(gridPointCount, 0.0, rStart);
    if (hasFitter(1)) { if (!fillLinkule(1, realCentral.data1Based())) return; }
    else realCentral.fillWoodsSaxon(pw.tvReal, reaction.integrationGrid.R, pw.taReal);
    coul.add(realCentral);   // coul[i] = Coulomb + real
    {
        double rho = 1.0e-20;
        // problemType == 24 dead (CC mode unreachable; InputParser sets 6/20/22).
        double aj = (channelIndexIn == 3) ? 0.0 : 1.0;
        double temp = (channelIndexIn == 3) ? 0.0 : 1.0;
        double* combined = coul.data1Based();
        for (i = 1; i <= gridPointCount; i++) {
            double ak1 = eInv * combined[i] - aj;
            realVPointer[i] = temp - ak1 * (h2 / 12.0);
            centrPointer[i] = -(h2 / 12.0) / (rho * rho);
            rho = rho + h;
        }
    }

    // --- Imaginary central (slot 2) -> imagV ---
    OpticalPotential imagCentral; imagCentral.resize(gridPointCount, 0.0, rStart);
    if (hasFitter(2)) { if (!fillLinkule(2, imagCentral.data1Based())) return; }
    else imagCentral.fillWoodsSaxon(pw.tvImag, params.rI, pw.taImag);
    {
        double* v = imagCentral.data1Based();
        for (i = 1; i <= gridPointCount; i++)
            imagVPointer[i] = -(h2 / (12.0 * reaction.energies.E)) * v[i];
    }

    // --- Surface imaginary (slot 13) -> imagV, only if active ---
    if (params.vSi != 0.0) {
        OpticalPotential surface; surface.resize(gridPointCount, 0.0, rStart);
        if (hasFitter(13)) { if (!fillLinkule(13, surface.data1Based())) return; }
        else surface.fillSurface(params.vSi, params.rSi, params.aSi);
        double* v = surface.data1Based();
        for (i = 1; i <= gridPointCount; i++)
            imagVPointer[i] = imagVPointer[i] - (h2 / (12.0 * reaction.energies.E)) * v[i];
    }

    // --- Spin-orbit real (slot 3) -> soR[i-1], 0-based, if active ---
    if (soRPointer != nullptr) {
        OpticalPotential soReal; soReal.resize(gridPointCount, 0.0, rStart);
        if (hasFitter(3)) { if (!fillLinkule(3, soReal.data1Based())) return; }
        else soReal.fillSpinOrbit(params.vSo, params.rSo, params.aSo);
        double* v = soReal.data1Based();
        for (i = 1; i <= gridPointCount; i++)
            soRPointer[i-1] = (-h2 / (12.0 * reaction.energies.E)) * v[i];
    }

    // --- Spin-orbit imaginary (slot 4) -> soI[i-1], 0-based, if active ---
    if (soIPointer != nullptr) {
        OpticalPotential soImag; soImag.resize(gridPointCount, 0.0, rStart);
        if (hasFitter(4)) { if (!fillLinkule(4, soImag.data1Based())) return; }
        else soImag.fillSpinOrbit(params.vSoi, params.rSoi, params.aSoi);
        double* v = soImag.data1Based();
        for (i = 1; i <= gridPointCount; i++)
            soIPointer[i-1] = (-h2 / (12.0 * reaction.energies.E)) * v[i];
    }

}


// ============================================================================
// SETPOT — Sets up potential parameters for a single channel
// (large routine: 447 lines of Fortran)
// ============================================================================
bool Reaction::setupOpticalPotential()
{
    Reaction& reaction = *this;
    bool ok = true;

    // Setup R0 mass factor
    double aLighter = std::min(reaction.masses.massProj, reaction.masses.massTgt);
    double aHeavier = reaction.masses.massProj + reaction.masses.massTgt - aLighter;
    double oneThird = 1.0 / 3.0;

    // r0Type ∈ {0, 1} (parser sets 1 only on the r0target keyword; defaults
    // leaves 0). Cases 3 (R0SUM) and 4 (R0MTOT) had no parser plumbing —
    if (reaction.flags.r0Type == 0) {
        reaction.internalState.r0Mass = std::pow(aHeavier, oneThird);
        if (aLighter > 2.5) reaction.internalState.r0Mass += std::pow(aLighter, oneThird);
    } else {
        // R0TARGET (r0Type == 1)
        reaction.internalState.r0Mass = std::pow(aHeavier, oneThird);
    }


    if (reaction.opticalPotentialParams.V == reaction.internalState.undefValue) reaction.opticalPotentialParams.V = 0.0;
    if (reaction.opticalPotentialParams.vSo == reaction.internalState.undefValue) reaction.opticalPotentialParams.vSo = 0.0;

    // V real params must always be defined
    if (reaction.opticalPotentialParams.R0 != reaction.internalState.undefValue) {
        reaction.integrationGrid.R = reaction.opticalPotentialParams.R0 * reaction.internalState.r0Mass;
    } else if (reaction.integrationGrid.R == reaction.internalState.undefValue) {
        if (reaction.linkuleData.linkuleAddr[1][3] == 0) {
            std::printf("\n**** R OR R0 MUST BE DEFINED.\n");
            ok = false;
        }
    }
    if (reaction.opticalPotentialParams.A == reaction.internalState.undefValue) {
        std::printf("\n**** A MUST BE DEFINED.\n");
        ok = false;
    }

    if (reaction.opticalPotentialParams.vI != 0.0) {
        if (reaction.opticalPotentialParams.rI0 != reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.rI = reaction.opticalPotentialParams.rI0 * reaction.internalState.r0Mass;
        } else if (reaction.opticalPotentialParams.rI == reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.rI = reaction.integrationGrid.R;
            reaction.opticalPotentialParams.rI0 = reaction.opticalPotentialParams.R0;
            reaction.internalState.wasSet[1] = TRUE_F;
        }
        if (reaction.opticalPotentialParams.aI == reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.aI = reaction.opticalPotentialParams.A;
            reaction.internalState.wasSet[2] = TRUE_F;
        }
    }

    if (reaction.opticalPotentialParams.vSo != 0.0) {
        if (reaction.opticalPotentialParams.rSo0 != reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.rSo = reaction.opticalPotentialParams.rSo0 * reaction.internalState.r0Mass;
        } else if (reaction.opticalPotentialParams.rSo == reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.rSo = reaction.integrationGrid.R;
            reaction.internalState.wasSet[3] = TRUE_F;
        }
        if (reaction.opticalPotentialParams.aSo == reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.aSo = reaction.opticalPotentialParams.A;
            reaction.internalState.wasSet[4] = TRUE_F;
        }
    }

    if (reaction.opticalPotentialParams.vSoi != 0.0) {
        if (reaction.opticalPotentialParams.rSoi0 != reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.rSoi = reaction.opticalPotentialParams.rSoi0 * reaction.internalState.r0Mass;
        } else if (reaction.opticalPotentialParams.rSoi == reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.rSoi = reaction.opticalPotentialParams.rSo;
            reaction.internalState.wasSet[5] = TRUE_F;
        }
        reaction.internalState.wasSet[6] = (reaction.opticalPotentialParams.aSoi == reaction.internalState.undefValue);
        if (reaction.opticalPotentialParams.aSoi == reaction.internalState.undefValue) reaction.opticalPotentialParams.aSoi = reaction.opticalPotentialParams.aSo;
    }

    // Surface absorption
    if (reaction.opticalPotentialParams.vSi != 0.0) {
        if (reaction.opticalPotentialParams.rSi0 != reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.rSi = reaction.opticalPotentialParams.rSi0 * reaction.internalState.r0Mass;
        } else if (reaction.opticalPotentialParams.rSi == reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.rSi = reaction.opticalPotentialParams.rI;
            reaction.internalState.wasSet[7] = TRUE_F;
        }
        reaction.internalState.wasSet[8] = (reaction.opticalPotentialParams.aSi == reaction.internalState.undefValue);
        if (reaction.opticalPotentialParams.aSi == reaction.internalState.undefValue) reaction.opticalPotentialParams.aSi = reaction.opticalPotentialParams.aI;
    }

    // Coulomb
    if (reaction.charges.zProj != NOTDEF_INT && reaction.charges.zTarget != NOTDEF_INT) {
        if (reaction.opticalPotentialParams.rC0 != reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.rC = reaction.opticalPotentialParams.rC0 * reaction.internalState.r0Mass;
        } else if (reaction.opticalPotentialParams.rC == reaction.internalState.undefValue) {
            reaction.opticalPotentialParams.rC = 1.0;
            reaction.internalState.wasSet[9] = TRUE_F;
            // Check if folded Coulomb needed
            // hasNextBlock∈{0,1,2,3,4,6} per InputParser — the `!= 21` guard from
            // the Fortran days is dead. Folded-Coulomb requires hasNextBlock to
            // name a bound-state phase (3/4/6), not a projectile/target
            // optical setup (0/1/2).
            if (reaction.charges.zProj * reaction.charges.zTarget != 0 &&
                reaction.flags.hasNextBlock != 3 && reaction.flags.hasNextBlock != 4 && reaction.flags.hasNextBlock != 6) {
                std::printf("\n**** RC OR RC0 MUST BE DEFINED FOR BS.\n");
                ok = false;
            }
        }
    } else {
        std::printf("\n**** ZP AND ZT MUST ALWAYS BE DEFINED.\n");
        ok = false;
    }

    return ok;
}





// ============================================================================
// setupWavefunctionPotential — Gets L-dependent stuff for WAVELJ
// ============================================================================
bool Reaction::setupWavefunctionPotential()
{
    Reaction& reaction = *this;
    // NAMES[2][8] pool-slot label table (SIN/SOUT/COULFA/.../TOCEB) dropped

    int waveChannel = reaction.internalState.waveChannel;
    double eta = reaction.kin.etaCh[waveChannel];
    reaction.energies.E = reaction.distortedWave.channel[waveChannel].Ecm;
    int stepCount = reaction.distortedWave.channel[waveChannel].nGridSteps;
    int lOutMax = reaction.kin.lOutMax;
    int i = (lOutMax + 2);  // +2: RCWFN writes FC[lMax+1]
    int l1 = 0, l2 = lOutMax;

    // F_arr/G_arr/nF1s/nG1s: caller-owned per-channel vectors.
    reaction.distortedWave.channel[waveChannel].F_arr.assign(i, 0.0);
    reaction.distortedWave.channel[waveChannel].G_arr.assign(i, 0.0);
    reaction.distortedWave.channel[waveChannel].nF1sArr.assign(i, 0.0);
    reaction.distortedWave.channel[waveChannel].nG1sArr.assign(i, 0.0);

    // Compute Coulomb wave functions
    double rho = stepCount * reaction.distortedWave.channel[waveChannel].stepSize;
    std::vector<double> workVector(2 * (lOutMax + 2) + 1, 0.0);
    int returnCode;

    // Identical RCWFN error-return guard at both computeFG call sites (forward
    // rho and NBACK-steps-back rho1); returns true (caller then returns false)
    // on a nonzero return code. Only the printed rho value varies.
    auto bailIfRcwfnReturnError = [&](double rhoVal) {
        if (returnCode != 0) {
            std::printf("-***** ERROR RETURN FROM RCWFN:%8d%18.8g%18.8g%8d%8d\n",
                returnCode, eta, rhoVal, l1, l2);
            return true;
        }
        return false;
    };

    CoulombWaveFunction::computeFG(rho, eta, l1, l2,
        reaction.distortedWave.channel[waveChannel].F_arr.data(),
        workVector.data(),
        reaction.distortedWave.channel[waveChannel].G_arr.data(),
        workVector.data() + (lOutMax + 1),
        1.0e-14, returnCode);
    if (bailIfRcwfnReturnError(rho)) return false;

    // Calculate Coulomb function at NBACK steps back
    double rho1 = (stepCount - 4) * reaction.distortedWave.channel[waveChannel].stepSize;
    CoulombWaveFunction::computeFG(rho1, eta, l1, l2,
        reaction.distortedWave.channel[waveChannel].nF1sArr.data(),
        workVector.data(),
        reaction.distortedWave.channel[waveChannel].nG1sArr.data(),
        workVector.data() + (lOutMax + 1),
        1.0e-14, returnCode);
    if (bailIfRcwfnReturnError(rho1)) return false;


    reaction.distortedWave.channel[waveChannel].sigmaArr.assign(i, 0.0);
    coulombSigmaL(eta, lOutMax, reaction.distortedWave.channel[waveChannel].sigmaArr.data());

    // Allocate S-matrix space
    {
        int lxECount = reaction.distortedWave.channel[waveChannel].hasSpinorbit
                           ? (int)reaction.distortedWave.channel[waveChannel].twoSpin + 1
                           : 1;
        // indxeArr: 0-based (matches convertJtoL's indxePointer convention).
        reaction.distortedWave.channel[waveChannel].indxeArr.assign(3 * lxECount * lxECount + 2, 0);
        int* indxePointer = reaction.distortedWave.channel[waveChannel].indxeArr.data();  // 0-based
        i = (lxECount * lxECount + 1) / 2;
        // toceArr: 4-int records, 1-based via .data()-1.
        reaction.distortedWave.channel[waveChannel].toceArr.assign(4 * i, 0);
        int* tocePointer = reaction.distortedWave.channel[waveChannel].toceArr.data();  // 0-based (accessed [4*ii-4 .. 4*ii-1])
        int ii = 1;
        for (int lxP = 1; lxP <= lxECount; lxP++) {
            int lxTemp = lxP - 1;
            int lDeltaCount = (lxTemp + lxP + (lxP % 2)) / 2;
            int lDelta = 1 - lDeltaCount;
            int k = 3 * lxTemp * ((int)reaction.distortedWave.channel[waveChannel].twoSpin + 2);  // 0-based offset
            indxePointer[k]     = ii;
            indxePointer[k + 1] = lDelta;
            indxePointer[k + 2] = lDeltaCount;
            for (int deltaIndex = 1; deltaIndex <= lDeltaCount; deltaIndex++) {
                tocePointer[4 * ii - 4] = lDelta;
                tocePointer[4 * ii - 3] = lxTemp;
                tocePointer[4 * ii - 2] = 2 * lxTemp;
                tocePointer[4 * ii - 1] = 0;
                ii++;
                lDelta += 2;
            }
        }
        reaction.distortedWave.channel[waveChannel].nJStates = ii - 1;
        ii = 2 * reaction.distortedWave.channel[waveChannel].nJStates * (lOutMax + 1);
        reaction.distortedWave.channel[waveChannel].smatArr.assign(ii, 0.0);
    }
    return true;
}


// ============================================================================
// SETSPT — Sets up S-matrix pointers
// indxsBase and tocsBase are 1-based int* (reaction.inelastic.indxsPointer /
// reaction.inelastic.tocsPointer, or +per_ch_off for CC). counter is the size accumulator
// for passes 1/2 and unused in pass 3.
// ============================================================================
void Reaction::setupInelasticAngMomTable(int& counter, int* indxsBase, int* tocsBase, int pass)
{
    Reaction& reaction = *this;
    // problemType == 24 (CC) jpMin/jpMax/jtMin/jtMax/vertex.lBound seeding dropped

    reaction.inelastic.jtBase = std::abs((int)reaction.angMom.js[3] - (int)reaction.angMom.js[4]);
    reaction.inelastic.nJt = ((int)reaction.angMom.js[3] + (int)reaction.angMom.js[4] - reaction.inelastic.jtBase) / 2 + 1;
    reaction.inelastic.jpBase = std::abs((int)reaction.angMom.js[1] - (int)reaction.angMom.js[2]);
    reaction.inelastic.nJp = ((int)reaction.angMom.js[1] + (int)reaction.angMom.js[2] - reaction.inelastic.jpBase) / 2 + 1;
    reaction.inelastic.nLx = ((int)reaction.angMom.js[1] + (int)reaction.angMom.js[2] + (int)reaction.angMom.js[3] + (int)reaction.angMom.js[4]) / 2 + 1;

    switch (pass) {
        case 1: // Return required INDXS size
            counter += 3 * reaction.inelastic.nLx * reaction.inelastic.nJp * reaction.inelastic.nJt;
            return;

        case 2: // Fill in INDXS and compute nSpl
        {
            reaction.inelastic.nSpl = 0;
            int* indxsPointer = indxsBase;  // 1-based; indxsPointer[3k-2..3k]
            bool useExtendedJpRange = reaction.distortedWave.channel[1].hasSpinorbit || reaction.distortedWave.channel[2].hasSpinorbit;
            int jpMinExtended = reaction.inelastic.jpMin, jpMaxExtended = reaction.inelastic.jpMax;
            if (useExtendedJpRange) { jpMinExtended = reaction.inelastic.jpBase; jpMaxExtended = (int)reaction.angMom.js[1] + (int)reaction.angMom.js[2]; }

            for (int jTt = reaction.inelastic.jtMin; jTt <= reaction.inelastic.jtMax; jTt += 2) {
                for (int jTp = jpMinExtended; jTp <= jpMaxExtended; jTp += 2) {
                    int lxMin = std::abs(jTt - jTp) / 2;
                    int lxMax = (jTt + jTp) / 2;
                    if (!useExtendedJpRange) { lxMin = std::max(lxMin, reaction.inelastic.lxMin); lxMax = std::min(lxMax, reaction.inelastic.lxMax); }
                    if (lxMax < lxMin) continue;
                    for (int lx = lxMin; lx <= lxMax; lx++) {
                        int deltaCount = lx + 1 - ((lx + reaction.boundState.vertex[1].lBound + reaction.boundState.vertex[2].lBound) % 2);
                        if (deltaCount == 0) continue;
                        int k = (jTp - reaction.inelastic.jpBase + reaction.inelastic.nJp * (jTt - reaction.inelastic.jtBase)) / 2;
                        k = 1 + lx + reaction.inelastic.nLx * k;
                        int lDeltaMin = 1 - deltaCount;
                        indxsPointer[3 * k - 3] = reaction.inelastic.nSpl + 1;
                        indxsPointer[3 * k - 2] = lDeltaMin;
                        indxsPointer[3 * k - 1] = deltaCount;
                        reaction.inelastic.nSpl += deltaCount;
                    }
                }
            }
            counter += reaction.inelastic.nSpl;
            return;
        }

        case 3: // Fill in TOCS
        {
            int* tocsPointer  = tocsBase;   // 1-based write of tocsArr
            int* indxsPointer = indxsBase;  // 0-based read of indxsArr
            reaction.distortedWave.scatteringSolver.hasAnySpinorbit = reaction.distortedWave.channel[1].hasSpinorbit || reaction.distortedWave.channel[2].hasSpinorbit;
            int jpMinExtended = reaction.inelastic.jpMin, jpMaxExtended = reaction.inelastic.jpMax;
            if (reaction.distortedWave.scatteringSolver.hasAnySpinorbit) { jpMinExtended = reaction.inelastic.jpBase; jpMaxExtended = (int)reaction.angMom.js[1] + (int)reaction.angMom.js[2]; }
            int lxJpCount = reaction.inelastic.nLx * reaction.inelastic.nJp;
            int slotCount = lxJpCount * reaction.inelastic.nJt;

            for (int k = 1; k <= slotCount; k++) {
                int deltaCount = indxsPointer[3 * k - 1];
                if (deltaCount == 0) continue;
                int jTt = ((k - 1) / lxJpCount) * 2 + reaction.inelastic.jtBase;
                if (jTt < reaction.inelastic.jtMin || jTt > reaction.inelastic.jtMax) continue;
                int jTp = ((k - 1) / reaction.inelastic.nLx % reaction.inelastic.nJp) * 2 + reaction.inelastic.jpBase;
                if (jTp < jpMinExtended || jTp > jpMaxExtended) continue;
                int lx = (k - 1) % reaction.inelastic.nLx;
                if (!reaction.distortedWave.scatteringSolver.hasAnySpinorbit && (lx < reaction.inelastic.lxMin || lx > reaction.inelastic.lxMax)) continue;
                int lDelta = indxsPointer[3 * k - 2];
                int kOffset = indxsPointer[3 * k - 3];
                for (int deltaIndex = 1; deltaIndex <= deltaCount; deltaIndex++) {
                    tocsPointer[4 * kOffset - 3] = lDelta;
                    tocsPointer[4 * kOffset - 2] = lx;
                    tocsPointer[4 * kOffset - 1] = jTp;
                    tocsPointer[4 * kOffset]     = jTt;
                    lDelta += 2;
                    kOffset++;
                }
            }
            return;
        }
    }
}


// ============================================================================
// sFromI — Computes S matrix elements from transfer radial integrals
// ============================================================================
// accumulator promoted to a function-local static; the 1000 sentinel still
// suppresses the first "PTOLEMY" banner per program run (the rest of the
// pipeline is single-shot post-Phase-8, so persistence is not an issue).
void sFromI(int li, int liIndex, double* sMatR, double* sMatI, int* indxsFlat,
            double* xiReal, double* xiImag, int* iIndexFlat, int numIi,
            int* indxDwFlat, int* iDwfiFlat, int* iDwfoFlat,
            double* abs1, double* aTerm, double factor,
            int isInfoPrint, Reaction& reaction)
{
    static int lineCount = 1000;

    // Multi-dimensional index macros
    #define indxs(i,k) indxsFlat[((k)-1)*3 + (i) - 1]
    #define iIndex(i,k) iIndexFlat[((k)-1)*4 + (i) - 1]
    #define indxDw(i,k) indxDwFlat[((k)-1)*4 + (i) - 1]
    #define iDwfi(i,k) iDwfiFlat[((k)-1)*3 + (i) - 1]
    #define iDwfo(i,k) iDwfoFlat[((k)-1)*4 + (i) - 1]

    int jA = (int)reaction.distortedWave.channel[1].twoSpin;
    int jB = (int)reaction.distortedWave.channel[2].twoSpin;
    double amplitude = 0.0;

    // Print headings if necessary
    if (isInfoPrint) {
        lineCount += numIi + 1;
        if (lineCount > 58 && lineCount < 1000) {
            std::printf("1%58sP T O L E M Y\n", "");
            std::printf("%10sCOMPUTATION OF TRANSFER S-MATRIX ELEMENTS\n", "");
        }
        if (!reaction.distortedWave.scatteringSolver.hasAnySpinorbit) {
            std::printf("   L   L  LX%20c S-MATRIX%17c CANCELLATION\n", ' ', ' ');
            std::printf("   IN OUT%11c REAL%8c IMAG%7c MAGNITUDE%9c RI, RO\n", ' ', ' ', ' ', ' ');
        }
        lineCount = 7 + numIi;
        std::printf("\n");
    }

    // Loop through radial integrals
    int jpMax = jA + jB;
    for (int ii = 1; ii <= numIi; ii++) {
        int lxP = iIndex(4, ii);
        int kDw = iIndex(3, ii);
        int kWI = indxDw(1, kDw);
        int kWO = indxDw(2, kDw);
        int lasI = iDwfi(1, kWI) + li;
        int jPi = iDwfi(2, kWI) + 2 * li;
        int lo = iDwfo(1, kWO) + li;
        int lasO = iDwfo(2, kWO) + li;
        int jPo = iDwfo(3, kWO) + 2 * li;

        int ksBase = liIndex - 1;
        ksBase = reaction.inelastic.nSpl * ksBase;

        double temp = factor * aTerm[lxP] / std::sqrt(2.0 * lasI + 1.0);  // aTerm = atermArr (0-based)
        int phaseExponent = lasI + lasO + 2 * lxP + 1;
        if ((phaseExponent % 4) >= 2) temp = -temp;
        // Guard against NaN radial integrals (forbidden transitions)
        double xr = xiReal[ii - 1], xi = xiImag[ii - 1];  // xiReal/xiImag = liloR/liloI (0-based)
        if (std::isnan(xr)) xr = 0.0;
        if (std::isnan(xi)) xi = 0.0;
        double tempReal = temp * xr;
        double tempImag = temp * xi;
        double sR = tempReal;
        double sI = tempImag;
        if ((phaseExponent % 2) == 0) {
            // already correct
        } else {
            sR = -tempImag;
            sI = tempReal;
        }
        if (isInfoPrint) {
            amplitude = std::sqrt(sR * sR + sI * sI);
            if (amplitude == 0.0) amplitude = 1.0e-30;
            temp = std::fabs(temp) * abs1[ii - 1] / amplitude;  // abs1 = abs1Pointer (0-based)
        }

        if (!reaction.distortedWave.scatteringSolver.hasAnySpinorbit) {
            int k = ((int)reaction.boundState.vertex[1].jB - reaction.inelastic.jpBase + reaction.inelastic.nJp * ((int)reaction.boundState.vertex[2].jB - reaction.inelastic.jtBase)) / 2;
            k = 1 + lxP + reaction.inelastic.nLx * k;
            int kOffset = indxs(1, k);
            int lDeltaMin = indxs(2, k);
            int ks = ksBase + kOffset + (lasO - lasI - lDeltaMin) / 2;
            sMatR[ks - 1] += sR;  // sMatR/sMatI = smatRPointer/smatIPointer (0-based)
            sMatI[ks - 1] += sI;
            if (isInfoPrint)
                std::printf(" %3d%4d%3d%16.4g%12.4g%14.4g%12.2f\n", li, lo, lxP, sR, sI, amplitude, temp);
        } else {
            // S.O. / tensor case - 9-J symbols
            if (amplitude == 1.0e-30) continue;
            int jpiMin = jPi, jpiMax = jPi;
            if (!reaction.distortedWave.channel[1].hasSpinorbit) jpiMin = std::abs(2 * li - (int)reaction.distortedWave.channel[1].twoSpin);
            int jpoMin = jPo, jpoMax = jPo;
            if (!reaction.distortedWave.channel[2].hasSpinorbit) jpoMin = std::abs(2 * lo - (int)reaction.distortedWave.channel[2].twoSpin);
            for (int jpi = jpiMin; jpi <= jpiMax; jpi += 2) {
                for (int jpo = jpoMin; jpo <= jpoMax; jpo += 2) {
                    double sav9J = (jpi + 1) * (jpo + 1) * (2 * lxP + 1) * ((int)reaction.boundState.vertex[1].jB + 1);
                    sav9J = std::sqrt(sav9J) * wig9J((int)reaction.boundState.vertex[2].jB, 2 * lxP, (int)reaction.boundState.vertex[1].jB,
                        jpi, 2 * li, jA, jpo, 2 * lo, jB);
                    if (sav9J == 0.0) continue;
                    tempReal = sav9J * sR;
                    tempImag = sav9J * sI;
                    for (int jp = reaction.inelastic.jpBase; jp <= jpMax; jp += 2) {
                        int lxMin = std::max(std::abs((int)reaction.boundState.vertex[2].jB - jp) / 2, std::abs(lasO - lasI));
                        int lxMax = std::min(((int)reaction.boundState.vertex[2].jB + jp) / 2, lasO + lasI);
                        if (lxMin > lxMax) continue;
                        for (int lx = lxMin; lx <= lxMax; lx++) {
                            int k = (jp - reaction.inelastic.jpBase + reaction.inelastic.nJp * ((int)reaction.boundState.vertex[2].jB - reaction.inelastic.jtBase)) / 2;
                            k = 1 + lx + reaction.inelastic.nLx * k;
                            int kOffset = indxs(1, k);
                            int lDeltaMin = indxs(2, k);
                            double temp2 = sav9J;
                            if (!(lx == lxP && jp == (int)reaction.boundState.vertex[1].jB && lasI == li && lasO == lo)) {
                                temp2 = (jpi + 1) * (jpo + 1) * (2 * lx + 1) * (jp + 1);
                                temp2 = std::sqrt(temp2) * wig9J((int)reaction.boundState.vertex[2].jB, 2 * lx, jp,
                                    jpi, 2 * lasI, jA, jpo, 2 * lasO, jB);
                                if (temp2 == 0.0) continue;
                            }
                            if ((lx + lxP) % 2 != 0) temp2 = -temp2;
                            int ks = ksBase + kOffset + (lasO - lasI - lDeltaMin) / 2;
                            sMatR[ks - 1] += temp2 * tempReal;
                            sMatI[ks - 1] += temp2 * tempImag;
                        }
                    }
                }
            }
        }
    }

    #undef indxs
    #undef iIndex
    #undef indxDw
    #undef iDwfi
    #undef iDwfo
}
