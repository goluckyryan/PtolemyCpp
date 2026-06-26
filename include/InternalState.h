#pragma once

#include "ptolemy_types.h"   // logical1

// ============================================================================
// ============================================================================
struct InternalState {
    double   undefValue;          // sentinel: undefined floating-point value (set to a recognizable pattern)
    double   notDefSentinel;         // integer sentinel stored as double bits: 0xF0F0F0F0 = NOTDEF_INT (-252645136)
    int      boundChannel;         // index of current bound-state channel (0-based)
    int      waveChannel;         // index of current wave-function channel (0-based)
    double   eBnds[3];       // 1-based: energy bounds for search [1]=lower, [2]=upper (keV)
    int      iDone;          // 1=current sub-calculation is complete
    int      stripPickup;         // number of nucleons being transferred (stripping sign: +1=stripping, -1=pickup)
    int      lInMax;         // L_max for incoming channel (most L-values)
    int      iExcit;         // 1=excitation reaction (inelastic), 0=transfer
    double   r0Mass;         // R0 * A^(1/3) — radius parameter times mass^1/3 (fm)
    logical1 wasSet[30];     // 1-based [1..29]: flag array — wasSet[k]=true if parameter k was explicitly set.
                             //   Live indices: 1..9 (radii/diffuseness from source_potentials), 12..16, 29
    double   ratMass;         // mass ratio for reduced-mass calculation: ma*mA/(ma+mA)
    int      lSpecs[5];      // 1-based: L values of spectroscopic factors (projectile bound states)
    int      nodesP[5];      // 1-based: node counts for projectile bound-state solutions
    int      lSpcPt2;        // L value for target bound state
    int      nodePt2;        // node count for target bound-state solution
};

// Integer notDefSentinel value: the Fortran stores 0xF0F0F0F0 as an int bit pattern in notDefSentinel
// (int)internalState.notDefSentinel gives WRONG result (casts tiny double to 0)
// Use this constant for integer comparisons instead
constexpr int NOTDEF_INT = static_cast<int>(0xF0F0F0F0u); // = -252645136
