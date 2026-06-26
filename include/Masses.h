#pragma once
// Masses.h — reaction masses / Coulomb radii block
// Extracted from Reaction.h. Member: reaction.masses.
// aMs(reaction) returns masses.massesArr.data() - 1 (1-based);
// AMPTS(reaction) = &masses.massProj (2 contiguous doubles).

#include <array>

struct Masses {
    std::array<double, 5> massesArr;  // [0]=mass_a, [1]=mass_b, [2]=mass_A, [3]=mass_B, [4]=mass_x; aMs()[i] = massesArr[i-1]
    double massProj;     // mass of projectile (AMPTS[0])
    double massTgt;      // mass of target (AMPTS[1])
    double aM;            // reduced mass of scattering pair (AMU)
    double amxcs[6];      // 1-based [1..5]: cross section mass numbers, channels 1..5
    double amxcgs[6];     // 1-based [1..5]: center-of-mass corrected mass numbers, channels 1..5
    double amxgPt[3];     // 1-based [1..2]: mass numbers for gamma-phi transfer vertices
    double rcProj;           // Coulomb radius for projectile (fm; rC * massProj^(1/3))
    double rcTarget;           // Coulomb radius for target (fm; rC * massTgt^(1/3))
    double rc0Proj;          // Coulomb radius parameter for projectile (fm; bare)
    double rc0Target;          // Coulomb radius parameter for target (fm; bare)
};
