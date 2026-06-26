#pragma once
// kinematics.h — Pure particle kinematics
//
// NOTE: The original KANDM struct layout is preserved here for binary
// compatibility with linkule (which casts the whole struct to double*).
// Fields that belong to DistortedWave are kept here but accessed through
// reaction.distortedWave backward-compat references. When linkule is modernized,
// these fields will physically move.

struct Kinematics {
    // === Distorted-wave fields (accessed via reaction.distortedWave refs) ===
    double akIn;          // → reaction.distortedWave: incoming wave number k_i (fm^-1)
    double akOut;          // → reaction.distortedWave: outgoing wave number k_o (fm^-1)
    // === Pure kinematics ===
    double redMi;        // incoming reduced mass mu_i (MeV)
    double redMo;        // outgoing reduced mass mu_o (MeV)
    // === Distorted-wave fields ===
    double rScts[3];     // → reaction.distortedWave: scattering radii (fm)
    double aScts[3];     // → reaction.distortedWave: asymptotic matching radii (fm)
    double rcScts[3];    // → reaction.distortedWave: Coulomb scattering radii (fm)
    int    lOutMax;       // → reaction.distortedWave: L_max for outgoing channel (= lInMost + lxMax; see CWF_scattering.cpp:275, probe_print.cpp:1019)
    int    lCrit;        // → reaction.distortedWave: critical L
    double etaCh[3];    // → reaction.distortedWave: Sommerfeld parameter eta
    double rcSctP[3];    // → reaction.distortedWave: projectile Coulomb radii (fm)
    double rcSctT[3];    // → reaction.distortedWave: target Coulomb radii (fm)
    int    lCrits[3];    // → reaction.distortedWave: critical L per channel
    double tauRatio[3]; // → reaction.distortedWave: mass ratios for Coulomb
    // === Pure kinematics ===
    double aBar;         // average A for grid setup
};
// Access via reaction.kin (declared in reaction.h)
// Inline helpers (k_ch, REDMS, RMAXS) historical note — were in reaction.h
