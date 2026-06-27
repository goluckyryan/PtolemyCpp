#pragma once
// OpticalPotentialLibrary.h — named optical-model potential parameterizations.
//
// Ported verbatim from Cleopatra's potentials.h (Ryan Tang, 2018) — the same
// physics Ben Kay's globals_beta_v5 uses. The only structural change is that
// each parameterization returns its 16 numbers by value in an OMPset struct
// instead of writing global variables, so the library has no global state.
//
// Each parameterization gives a Woods-Saxon optical potential for a light
// projectile (n, p, d, t, 3He, alpha) on a target nucleus A(Z) at lab energy
// E (MeV). The single-letter codes match the DWBA-input potential codes:
//
//   deuteron : A H B D C L Q Z
//   proton   : K V M G P
//   A=3      : x l p c t h b
//   alpha    : s a f
//   custom   : X Y
//
// See potentialRef() for the citation of each code.

#include <string>

// 16 Woods-Saxon optical-model parameters (depths in MeV, radii/diffuseness in fm).
//   v,r0,a        : real volume
//   vi,ri0,ai     : imaginary volume
//   vsi,rsi0,asi  : imaginary surface
//   vso,rso0,aso  : real spin-orbit
//   vsoi,rsoi0,asoi: imaginary spin-orbit
//   rc0           : Coulomb radius
struct OMPset {
    double v   = 0, r0    = 0, a    = 0;
    double vi  = 0, ri0   = 0, ai   = 0;
    double vsi = 0, rsi0  = 0, asi  = 0;
    double vso = 0, rso0  = 0, aso  = 0;
    double vsoi= 0, rsoi0 = 0, asoi = 0;
    double rc0 = 0;
    bool   ok  = false;   // false if the potential code was not recognized
};

// Top-level dispatcher. name is a one-character potential code (see above).
// A,Z = target mass/charge, E = lab energy (MeV), Zproj = projectile charge.
// Returns an OMPset with ok=false if name is not a known code.
OMPset callPotential(const std::string& name, int A, int Z, double E, int Zproj);

// Human-readable citation / validity-range string for a potential code.
// Returns "" if the code is unknown.
std::string potentialRef(const std::string& name);
