// BoundState_formFactor.cpp — BSPROD/BSSET. BSPROD computes the bound-state
// product PHI(rT)*VEFF*PHI(rP); BSSET initializes the form-factor calculation.
// Uses enum dispatch + _impl pattern for the Fortran ENTRY.

#include "ptolemy_types.h"
#include "coulomb_utils.h"
#include "BoundState.h"
#include "Reaction.h"
#include "Constants.h"
#include <cstdio>
#include <cmath>

// Aitken-Lagrange interpolation for DWBA form-factor lookups. Definition
//   table[i+1] = F(i/stepInverse), table[1] = F(0).
//   aitkenOrder      = order (aitkenOrder+1 points used).
static double aitkenLagrange(double x, double stepInverse, const double* table, int tableSize, int aitkenOrder)
{
    static const double reciprocals[17] = { 0,
        1.0, 0.5, 1.0/3.0, 0.25, 0.2,
        1.0/6.0, 1.0/7.0, 0.125, 1.0/9.0, 0.1,
        1.0/11.0, 1.0/12.0, 1.0/13.0, 1.0/14.0, 1.0/15.0, 0.0625
    };

    double fs[17], deltas[17];
    // f init to 0 silences -Wmaybe-uninitialized for the (unreachable in
    // practice) aitkenOrder==0 empty-loop edge — every caller passes aitkenOrder ∈ {1, 2, 4}
    // so the for-i loop always runs at least once and assigns f before
    // return.
    double f = 0;
    double del1, fracIndex;
    int tableIndex, aitkenOrder2, startIndex, endIndex;

    if (x < 0) {
        std::printf("\n****** AITLAG - X < 0:%15.5G\n", x);
        return table[1];
    }

    fracIndex = x * stepInverse;
    tableIndex = (int)fracIndex;

    if (tableIndex >= tableSize) {
        std::printf("\n***** AITLAG - X TOO LARGE:%15.5G%10d%10d\n", x, tableIndex, tableSize);
        return table[tableSize];
    }

    aitkenOrder2 = aitkenOrder / 2;
    startIndex = tableIndex - aitkenOrder2;
    if (startIndex < 0) startIndex = 0;
    endIndex = startIndex + aitkenOrder;
    if (endIndex >= tableSize) endIndex = tableSize - 1;
    startIndex = endIndex - aitkenOrder;

    fs[1] = table[startIndex + 1];
    del1 = fracIndex - startIndex;
    deltas[1] = del1;

    for (int i = 1; i <= aitkenOrder; i++) {
        f = table[startIndex + 1 + i];
        del1 = del1 - 1;
        for (int j = 1; j <= i; j++) {
            f = (f * deltas[j] - fs[j] * del1) * reciprocals[i + 1 - j];
        }
        deltas[i + 1] = del1;
        fs[i + 1] = f;
    }

    return f;
}

// ============================================================================
//
// BSPROD: Computes bound-state product PHI(rT)*VEFF*PHI(rP)
// BSSET: entry point — initializes form-factor calculation
//
// Uses enum dispatch + _impl pattern for the BSSET entry point
// Returns: 0 = normal return, 1 = alternate return (rP or rT exceeds max)
// ============================================================================

// SAVE variables (persist across calls, shared between evaluateFormFactorImpl
namespace {
int scChannel;
int verbosity;
int isOuterRadius;   // file-scope persistent flag (single TU).
double rnScattering, rnCore, vOpt, aOpt;
}

static int evaluateFormFactorImpl(
    double& fpFt, int ffKind, double rA, double rB, double cosTheta,
    const double* scatPointer, int aitkenOrder, double& rP, double& rT,
    Reaction& reaction)
{
    // Locals
    double fP, fT, factor, delVNu, rScattering, rCore;
    double pot, delVCore, delVScattering, delV, xCore, fCore, xScattering, fScattering;
    double rBound, vEffective;
    logical forceInterpolation;

    // ---------------------------------------------------------------
    // BSPROD entry point
    // ---------------------------------------------------------------
    //
    // STEP 1 - COMPUTE  PHI V PHI   OR   PHI' V PHI'
    //
    rP = (reaction.boundState.data.s2*rA)*(reaction.boundState.data.s2*rA) + (reaction.boundState.data.t2*rB)*(reaction.boundState.data.t2*rB) +
       2*(reaction.boundState.data.s2*rA)*(reaction.boundState.data.t2*rB)*cosTheta;
    rT = (reaction.boundState.data.s1*rA)*(reaction.boundState.data.s1*rA) + (reaction.boundState.data.t1*rB)*(reaction.boundState.data.t1*rB) +
       2*(reaction.boundState.data.s1*rA)*(reaction.boundState.data.t1*rB)*cosTheta;
    if ( std::fabs(rP) < 1.0e-8 )  rP = 0;
    if ( std::fabs(rT) < 1.0e-8 )  rT = 0;
    rP = std::sqrt(rP);
    rT = std::sqrt(rT);

    fpFt = 0;
    if ( rP > reaction.boundState.vertex[1].boundMx  ||  rT > reaction.boundState.vertex[2].boundMx )  return 1;

    fP = reaction.boundState.vertex[1].vMax;
    fT = reaction.boundState.vertex[2].vMax;
    forceInterpolation = (( ffKind) % (2 )) != 0;
    delVNu = 0;

    if ( forceInterpolation || rT > reaction.boundState.vertex[2].rLMax )
      fT = aitkenLagrange(rT, reaction.boundState.vertex[2].boundSp, reaction.boundState.vertex[2].jbdPointer, reaction.boundState.vertex[2].nSpBd, aitkenOrder);  // per-vertex
    if ( forceInterpolation || rP > reaction.boundState.vertex[1].rLMax )
      fP = aitkenLagrange(rP, reaction.boundState.vertex[1].boundSp, reaction.boundState.vertex[1].jbdPointer, reaction.boundState.vertex[1].nSpBd, aitkenOrder);  // per-vertex
    factor = 1;
    fpFt = fP*fT;
    // nuConL∈{2,3} always (defaults seeds 3, linkule clash demotes to 2);
    // the `nuConL == 2 || nuConL == 3` disjunct is permanently true.
    if ( std::fabs(fpFt) >= Constants::smlNum ) {

    //
    // MULTIPLY V*PHI BY (V+DV)/V ; THIS GIVES (V+DV)*PHI
    // DV IS CORE CORRECTION TO VEFF
    //
    rScattering = isOuterRadius ? rB : rA;
    rCore = std::sqrt( ((reaction.boundState.data.s1-reaction.boundState.data.s2)*rA)*((reaction.boundState.data.s1-reaction.boundState.data.s2)*rA)
      + ((reaction.boundState.data.t1-reaction.boundState.data.t2)*rB)*((reaction.boundState.data.t1-reaction.boundState.data.t2)*rB)
      + 2*((reaction.boundState.data.s1-reaction.boundState.data.s2)*rA)*((reaction.boundState.data.t1-reaction.boundState.data.t2)*rB)*cosTheta );
    vcsq12(rCore, pot, 3);
    delVCore = pot;
    vcsq12(rScattering, pot, scChannel);
    delVScattering = pot;
    delV = delVCore - delVScattering;
    if ((reaction.flags.nuConL == 3) && (aitkenOrder != 2)) {
        xCore = (rCore - rnCore) / aOpt;
        if (xCore >= Constants::BIGLOG) {
            fCore = 0;
        } else if (xCore <= -Constants::BIGLOG) {
            fCore = 1;
        } else {
            fCore = 1 / (1 + std::exp(xCore));
        }
        xScattering = (rScattering - rnScattering) / aOpt;
        if (xScattering >= 174) {
            fScattering = 0;
        } else if (xScattering <= -175) {
            fScattering = 1;
        } else {
            fScattering = 1 / (1 + std::exp(xScattering));
        }
        delVNu = vOpt * (fCore - fScattering);
        delV = delV + delVNu;
    }
    // IVRTEX permanently 1: rBound = rT*0 + rP*1 = rP.
    rBound = rP;
    vEffective = aitkenLagrange(rBound, reaction.boundState.vertex[1].boundSp, reaction.boundState.vertex[1].getPotential(), reaction.boundState.vertex[1].nSpBd, aitkenOrder);
    factor = 1 + delV / vEffective;
    fpFt = factor * fpFt;
    if ( verbosity >= 9 )
      std::printf(" BSP%8.4f%8.4f%8.4f%8.4f%8.4f%12.4g%12.4g%12.4g%12.4g%12.4g%12.4g%12.4g\n",
        rA, rB, rP, rT, rCore, vEffective, delVScattering, delVCore, delVNu, delV, factor, fpFt);
    }

    //
    // NOW WE MULTIPLY IN THE PSI'S IF DESIRED.
    // IN FACT WE MULTIPLY IN R*PSI
    //
    if ( ffKind <= 2 )  return 0;

    if ( rA > reaction.boundState.data.scatRMax ) {
        fpFt = fpFt * rA;
    } else {
        fpFt = fpFt *
          aitkenLagrange( rA, reaction.boundState.data.scatInvStep, scatPointer, reaction.boundState.data.scatPointCount, aitkenOrder );
        // AVOID UNDERFLOWS
        if ( std::fabs(fpFt) <= 1.0e-34 ) {
            fpFt = 0;
            return 0;
        }
    }

    if ( rB > reaction.boundState.data.scatRMax ) {
        fpFt = fpFt * rB;
    } else {
        fpFt = fpFt *
          aitkenLagrange( rB, reaction.boundState.data.scatInvStep, scatPointer, reaction.boundState.data.scatPointCount, aitkenOrder );
    }
    return 0;
}

// ---------------------------------------------------------------
// BSSET entry point
// INITIALIZE THE FORM-factor CALCULATION
// ---------------------------------------------------------------

static void setupFormFactorsImpl(Reaction& reaction)
{
    int i;
    int zS1, zS2, zC1, zC2;
    int needsCoreCorrection;
    double temp, rRp, rRt;
    double factor;
    double amA3, amB3, amBgA3, amBgB3, power;
    double rCcoP, rCcoT, rcpSct, rctSct;

    //
    // COMPUTE APPROPRIATE  V*PHI  AND STORE IT IN "VPHI".
    //
    verbosity = ((reaction.flags.printLevel / 10000) % (10));
    // IVRTEX permanently 1: "active vertex" = vertex[1], "other vertex" = vertex[2].
    reaction.boundState.allocateFormFactor(reaction.boundState.vertex[1].nSpBd);
    double* vphiPointer = reaction.boundState.getFormFactor();  // 1-based
    // 1-based views into per-vertex bsPotential / wavefunction.
    double* bd1Pointer  = reaction.boundState.vertex[1].getWavefunction();           // 1-based
    double* bd2Pointer  = reaction.boundState.vertex[2].getWavefunction();           // 1-based
    double* potPointer  = reaction.boundState.vertex[1].getPotential();              // 1-based
    reaction.boundState.vertex[1].vMax = 0;
    reaction.boundState.vertex[2].vMax = 0;

    rRp = -reaction.boundState.vertex[1].bsVstep;
    for (i = 1; i <= reaction.boundState.vertex[1].nSpBd; i++) {
        rRp = rRp + reaction.boundState.vertex[1].bsVstep;
        temp = bd1Pointer[i] * potPointer[i];
        double absTemp = std::fabs(temp);
        if ( absTemp > reaction.boundState.vertex[1].vMax ) {
            reaction.boundState.vertex[1].rLMax = rRp;
            reaction.boundState.vertex[1].vMax = absTemp;
        }
        vphiPointer[i] = temp;
    }
    rRt = -reaction.boundState.vertex[2].bsVstep;
    for (i = 1; i <= reaction.boundState.vertex[2].nSpBd; i++) {
        temp = bd2Pointer[i];
        rRt = rRt + reaction.boundState.vertex[2].bsVstep;
        double absTemp = std::fabs(temp);
        if ( absTemp > reaction.boundState.vertex[2].vMax ) {
            reaction.boundState.vertex[2].rLMax = rRt;
            reaction.boundState.vertex[2].vMax = absTemp;
        }
    }

    if (reaction.flags.printLevel % 10 >= 4)  std::printf(" PROJECTILE AND TARGET FORMFACTOR MAXIMA:%15.5g%15.5g\n"
      " LOCATION OF PROJECTILE AND TARGET FORM-factor MAXIMA:%15.5g%15.5g\n",
      reaction.boundState.vertex[1].vMax, reaction.boundState.vertex[2].vMax,
      reaction.boundState.vertex[1].rLMax, reaction.boundState.vertex[2].rLMax);


    // wavefunction already populated by solve() directly into vertex[vv].wavefunction
    // (no pool→vector copy needed; pool slot no longer exists).
    // Set jbdPointer pointers: active vertex → formFactor, other → vertex wavefunction
    // IVRTEX permanently 1: active vertex (gets formFactor) is vertex[1].
    reaction.boundState.vertex[1].jbdPointer = reaction.boundState.getFormFactor();          // class-owned VPHI
    reaction.boundState.vertex[2].jbdPointer = reaction.boundState.vertex[2].getWavefunction(); // vertex-owned wavefunction
    reaction.boundState.vertex[1].boundSp = 1 / reaction.boundState.vertex[1].bsVstep;
    reaction.boundState.vertex[2].boundSp = 1 / reaction.boundState.vertex[2].bsVstep;
    // (vertex[1].nSpBd / vertex[2].nSpBd already written at solve loop L990;

    //
    // EXPONENTIAL TAIL FACTORS FOR H-INTERPOLATION
    //
    reaction.boundState.vertex[1].alpha = sqrt( -2*reaction.boundState.vertex[1].bsMass*reaction.internalState.eBnds[1] ) / Constants::hbar_c;
    reaction.boundState.vertex[2].alpha = sqrt( -2*reaction.boundState.vertex[2].bsMass*reaction.internalState.eBnds[2] ) / Constants::hbar_c;
    if ( verbosity >= 4 )  std::printf("\n EXPONENTIAL DECAYS OF FORM factor:  RP =%13.5g     RT =%13.5g\n",
      reaction.boundState.vertex[1].alpha, reaction.boundState.vertex[2].alpha);

    // Sync data paired fields into vertex[1] (projectile) and vertex[2] (target)
    // (vertex[1].rLMax/VMAX already set above; old data.RLPMAX/VPMAX sync deleted)
    // (vertex[1].boundSp already set above; old data.BNDSPP sync deleted)
    // (vertex[1].nSpBd already set above; old data.NSPBDP sync deleted)


    //
    // COMPUTE POTENTIAL PARAMETERS FOR CORE TERMS IN VEFF
    //
    power = 1.0 / 3.0;
    amA3 = std::pow(reaction.masses.massesArr[0], power);
    amB3 = std::pow(reaction.masses.massesArr[1], power);
    amBgA3 = std::pow(reaction.masses.massesArr[2], power);
    amBgB3 = std::pow(reaction.masses.massesArr[3], power);
    factor = 1;

    //
    // FOLDED COULOMB POTENTIALS IN BOTH CHANNELS
    //
    needsCoreCorrection = 0;
    for (i = 1; i <= 2; i++) {
        if ((reaction.kin.rcSctP[i] == 0) || (reaction.kin.rcSctT[i] == 0)) needsCoreCorrection = 1;
    }

    //
    // phiSign = +1 ; STRIPPING
    // phiSign = -1 ; PICK-UP
    //
    if (reaction.boundState.data.phiSign != -1) {

        // STRIPPING ; PROJECTILE VERTEX
        scChannel = 2;
        zS1 = reaction.charges.zArray[2];
        zS2 = reaction.charges.zArray[4];
        zC1 = reaction.charges.zArray[2];
        zC2 = reaction.charges.zArray[3];
        factor = (needsCoreCorrection == 0) ? (amBgA3 / amBgB3)
                                            : (amBgA3 + amB3) / (amBgB3 + amB3);
        rcpSct = reaction.kin.rcSctP[2];
        rctSct = reaction.kin.rcSctT[2];
        rCcoP = rcpSct;
        rCcoT = rctSct * factor;
        isOuterRadius = TRUE_F;
        if (reaction.flags.nuConL == 3) {
            rnScattering = reaction.kin.rScts[2];
            rnCore = rnScattering * ((amBgA3 + amB3) / (amBgB3 + amB3));
            vOpt = -reaction.distortedWave.channel[2].v0R;
            aOpt = reaction.kin.aScts[2];
        }
    } else {

    //
    // PICK-UP ; PROJECTILE VERTEX
    //
        scChannel = 1;
        zS1 = reaction.charges.zArray[1];
        zS2 = reaction.charges.zArray[3];
        zC1 = reaction.charges.zArray[1];
        zC2 = reaction.charges.zArray[4];
        factor = (needsCoreCorrection == 0) ? (amBgB3 / amBgA3)
                                            : (amBgB3 + amA3) / (amBgA3 + amA3);
        rcpSct = reaction.kin.rcSctP[1];
        rctSct = reaction.kin.rcSctT[1];
        rCcoP = rcpSct;
        rCcoT = rctSct * factor;
        isOuterRadius = FALSE_F;
        if (reaction.flags.nuConL == 3) {
            rnScattering = reaction.kin.rScts[1];
            rnCore = rnScattering * ((amBgB3 + amA3) / (amBgA3 + amA3));
            vOpt = -reaction.distortedWave.channel[1].v0R;
            aOpt = reaction.kin.aScts[1];
        }
    }

    //
    //
    setVsq(rCcoP, rCcoT, zC1, zC2, 3);

    if (verbosity <= 4) return;
    std::printf("\n  ++++ SETRS ++++\n   CHANNEL %2d(Z1,Z2) = ( %2d %2d )\n"
      "     CORE  (Z1,Z2) = (%2d %2d) NUCONL = %2d\n",
      scChannel, zS1, zS2, zC1, zC2, reaction.flags.nuConL);
    std::printf("   SCATTERING :(RT,RP) = %12.6g %12.6g(COULOMB)\n"
      "   CORE :(RT,RP) = %12.6g %12.6g(COULOMB)\n",
      rcpSct, rctSct, rCcoT, rCcoP);
    if (reaction.flags.nuConL != 3) return;
    std::printf("\n   (rnScattering,rnCore)=%12.6g %12.6g\n"
      "   (-V,A)=(%12.6g %12.6g %12.6g)\n",
      rnScattering, rnCore, vOpt, aOpt, aOpt);
}


// BoundState::evaluateFormFactor method
// Returns 0 = normal, 1 = alternate return (rP or rT out of range)
int BoundState::evaluateFormFactor(double& fpFt, int ffKind, double rA, double rB, double cosTheta,
           const double* scatPointer, int aitkenOrder, double& rP, double& rT,
           Reaction& reaction)
{
    return evaluateFormFactorImpl(fpFt, ffKind, rA, rB, cosTheta,
                          scatPointer, aitkenOrder, rP, rT, reaction);
}

// BoundState::setupFormFactors method
void BoundState::setupFormFactors(Reaction& reaction)
{
    setupFormFactorsImpl(reaction);
}
