#pragma once
// inelastic.h — Inelastic scattering class (DWBA)
//
// Handles: (p,p'), (α,α'), (d,d'), etc.
//
// Takes: Reaction + DistortedWave(×2) + β deformation parameters
// Does:  Integration grid, inelastic radial integrals, angular coupling
// Output: dσ/dΩ(θ), analyzing powers, spin observables
//
// Absorbs: INGRST, INRDIN, LINTRP, XSECTN (for inelastic)

class Reaction;

class InelasticReaction {
public:
    explicit InelasticReaction(Reaction& reaction);
    // Full inelastic calculation
    bool calculate();

private:
    Reaction& reaction_;
};
