#pragma once
// distorted_wave.h — DistortedWave class
// Everything about solving the scattering Schrödinger equation:
// wave numbers, Coulomb parameters, S-matrix, phase shifts, L-space,
// matching radii, wavefunctions, integration grid.

#include "ScatteringSolver.h"
#include <vector>

class DistortedWave {
public:
    // -----------------------------------------------------------------------
    // Channel — per-channel scattering data
    // Channel 1 = incoming (a+A), Channel 2 = outgoing (b+B).
    // -----------------------------------------------------------------------
    struct Channel {
        // --- Grid ---
        double stepSize     = 0.0;   // Numerov step size h (fm)
        double rStart      = 0.0;   // radial start for integration
        int    nGridSteps  = 0;     // number of Numerov steps
        double asymptopia   = 0.0;   // max integration radius (fm)

        // --- Optical potential ---
        double v0R          = 0.0;   // real central depth (MeV)
        double v0I          = 0.0;   // imaginary central depth
        double v0Si         = 0.0;   // surface imaginary depth
        double rI           = 0.0;   // imaginary radius
        double aI           = 0.0;   // imaginary diffuseness
        // Readers use reaction.rxn.{rSo,ASO,RSOI,ASOI} directly.
        double rSi          = 0.0;   // surface imaginary radius
        double aSi          = 0.0;   // surface imaginary diffuseness

        // --- Channel state ---
        double Ecm          = 0.0;   // channel cm energy (MeV)
        int    twoSpin     = 0;     // 2*spin+1 of projectile
        int    hasSpinorbit = 0;    // has SO potential
        int    nJStates        = 0;     // number of J values

        // CC writer in CWF_scattering FANDG path: 2*lMaxMax+2 doubles (dead in tests).
        // IREDEF→resize at L-extrapolation (interpolation.cpp).
        std::vector<double> smatArr;
        // Sized by source_potentials (non-CC: (lOutMax+2)*MAX(NUMFIT,1)) or
        // CWF_scattering CC path (lMaxMax+1). Resized down by IREDEF→resize at
        // L-extrapolation time (Born approximation bigRSw).
        std::vector<double> sigmaArr;
        // Non-CC writer in source_potentials.cpp: size = (lOutMax+2)*std::max(NUMFIT,1) doubles.
        // (CC writer at CWF_scattering FANDG path also assigns F_arr with packed
        //  size 4*NBASDF*n_L_values/FACFR4; CC mode is currently dead in tests.)
        std::vector<double> F_arr;
        std::vector<double> G_arr;
        // No CC alternate; size = (lOutMax+2)*std::max(NUMFIT,1) doubles.
        std::vector<double> nF1sArr;
        std::vector<double> nG1sArr;

        // --- Solver state ---
        int    lastL        = 0;     // last L computed
        int    lastNf       = 0;     // last NF
        int    lSkips       = 0;     // L values skipped
        int    statsCode    = 0;     // particle statistics: 1/2 = identical (even/odd 2*spinProj), 3 = non-identical
        int    nStp2s       = 0;     // 2nd-order steps

        // Sized gridPointCount in setupScatteringWaves (imvs gets gridPointCount+1 — extra guard cell).
        // rlvPointer/imvPointer/centPointer cache 1-based pointers (.data()-1) into these vectors.
        // integration_grid INGRST swaps in local scratch via the P_* pointers and
        // restores them after the deformed-potential pass.
        std::vector<double> rlvsArr;
        std::vector<double> imvsArr;
        std::vector<double> centrArr;
        std::vector<int> indxeArr; // L-index array, 0-based
        std::vector<int> toceArr;  // TOC array (4-int records, 1-based access via .data()-1)

        // --- Tensor potential ---
        double xFacs[4] = {};       // cross-section normalization (1-based [1..3])

        // --- Pointer caches into pool / class-owned vectors (set by allocators) ---
        double* soRPointer = nullptr;       // 0-based spin-orbit real ptr
        double* soIPointer = nullptr;       // 0-based spin-orbit imag ptr
        double* rlvPointer  = nullptr;       // 1-based real local potential ptr
        double* imvPointer  = nullptr;       // 1-based imag local potential ptr
        double* centPointer = nullptr;       // 1-based centrifugal term ptr (ICENTR cache)
    };

    // -----------------------------------------------------------------------
    // Per-channel data: channel[1]=incoming, channel[2]=outgoing
    // -----------------------------------------------------------------------
    Channel channel[3] = {};  // 1-based: [1]=incoming, [2]=outgoing

    // -----------------------------------------------------------------------
    // Scattering solver (HAS-A): owns shared work arrays, solver scalars,
    // ScatteringSolver lives inside DistortedWave instead of being a sibling
    // global. Future: flatten members into DistortedWave directly.
    // -----------------------------------------------------------------------
    ScatteringSolver scatteringSolver = {};

};

