#pragma once
// OpticalPotentialParams.h — Woods-Saxon optical-potential INPUT parameters
// (V, A, R0, vSo, ...). The flat user-input struct; the computed grid V(r) is the
// separate OpticalPotential class (OpticalPotential.h). Member: reaction.opticalPotentialParams.
// (Renamed from OpticalPotential / Reaction::opticalPotential, ponytail-flat Phase A.)

struct OpticalPotentialParams {
    // Diffusenesses (fm)
    double A;       // real central
    double aI;      // imaginary central
    double aSo;     // real spin-orbit
    double aSoi;    // imaginary spin-orbit
    double aSi;     // surface imaginary
    // Depths (MeV; negative = attractive for real)
    double V;       // real central
    double vI;      // imaginary central
    double vSo;     // real spin-orbit
    double vSoi;    // imaginary spin-orbit
    double vSi;     // surface imaginary
    // Radii (fm)
    double R0;      // real central radius parameter (bare, before A^1/3)
    double rI;      // imaginary central radius (scaled)
    double rI0;     // imaginary central radius parameter (bare)
    double rSo;     // real spin-orbit radius
    double rSo0;    // real spin-orbit radius parameter
    double rSoi;    // imaginary spin-orbit radius
    double rSoi0;   // imaginary spin-orbit radius parameter
    double rSi;     // surface imaginary radius
    double rSi0;    // surface imaginary radius parameter
    double rC;      // Coulomb radius (scaled = rC * A^1/3)
    double rC0;     // Coulomb radius parameter (bare)
    // L-value limits for amplitude sum
    double alMnMt;  // minimum L (as float)
    double alMxMt;  // maximum L (as float)
};
