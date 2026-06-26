#pragma once
// Spectroscopy.h — spectroscopic amplitude / factor block
// Extracted from Reaction.h. Member: reaction.spec.

struct Spectroscopy {
    double spam;             // spectroscopic amplitude product specAmpProj*specAmpTgt (combined)
    double specAmpProj;    // projectile spectroscopic amplitude (dimensionless; default 1)
    double specAmpTgt;     // target spectroscopic amplitude (dimensionless; default 1)
    double specFactorProj; // projectile spectroscopic factor (specAmpProj^2; input keyword SPFACP)
    double specFactorTgt;  // target spectroscopic factor (specAmpTgt^2; input keyword SPFACT)
};
