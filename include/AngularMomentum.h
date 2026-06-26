#pragma once
// AngularMomentum.h — angular-momentum / parity quantum numbers for a channel.
// Holds the full Reaction group-5 block: parities, L-range, J-values, spins, and
// the bound-state radial node count (the n companion to L/jProj).

struct AngularMomentum {
    int nNodes = 0;        // bound-state radial node count (n companion to L/jProj)
    int parities[6] = {};  // 1-based per-particle parities [1..5]
    int parityPt[3] = {};  // 1-based projectile/target parities [1..2]
    int parity = 0;        // total channel parity (+1 / -1 / 0 = undefined)
    int L = 0;             // orbital angular momentum of current channel
    int lMax = 0;          // maximum partial-wave L for this channel
    int lMin = 0;          // minimum partial-wave L for this channel
    int lMaxAdditional = 0; // extra L added past lCrit when sizing lMax
    double J = 0;          // total angular momentum of current channel
    double js[6] = {};     // 1-based [1..5]: J values for each channel
    double jProj  = 0;     // J of projectile bound state (transfer)
    // spinProj/spinTarget MUST stay adjacent (channel_setup jSpts pointer alias).
    double spinProj = 0;   // spin of projectile particle a
    double spinTarget = 0; // spin of target particle A
};
