#pragma once
// Energies.h — energy / Q-value / excitation block
// Extracted from Reaction.h. Member: reaction.energies.

struct Energies {
    double E;        // kinetic energy (MeV) for Schrödinger equation (cm energy or binding energy)
    double eCm;      // center-of-mass energy (MeV)
    double eLab;     // lab-frame energy (MeV)
    double Q;        // Q-value of reaction (MeV): Q = E_out - E_in = E_b + E_B - E_a - E_A
    double exs[6];   // 1-based [1..5]: excitation energies for each state (MeV); exs[0] unused (padding)
    double exsPt[3]; // 1-based [1..2]: excitation energies for projectile/target vertices (MeV)
};
