#pragma once
// Flags.h — orthogonal control flags for a reaction calculation.
// These select output framing and dispatch; none participate in the physics math.

struct Flags {
    int outputInLab = 0;     // 0 = c.m. output angles, !=0 = lab frame
    int problemType = 0;     // calc dispatch code (InputParser sets 6/20/22)
    int r0Type = 0;          // 0 = default radius convention, 1 = R0TARGET
    int excitationType = 0;  // extrapolation form (defaults.cpp seeds 1)
    int printLevel = 0;      // packed verbosity digits (defaults.cpp seeds 10001)
    int isElastic = 0;       // !=0 = elastic-scattering dispatch
    int nuConL = 0;          // continuum-L level selector (defaults.cpp seeds 3; ∈{2,3})
    int hasNextBlock = 0;    // channel-block dispatch code (InputParser sets 0..6)
};
