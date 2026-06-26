#pragma once
// IntegrationGrid.h — Numerov / DWBA integration-grid parameters
// Extracted from Reaction.h. Member: reaction.integrationGrid.

struct IntegrationGrid {
    double stepSize;       // Numerov step size (fm)
    double accuracy;        // accuracy parameter for Numerov integration (convergence tolerance)
    double accuracyInel;    // accuracy for inelastic convergence
    double asymptopia;      // asymptotic matching radius (fm) — wavefunction integrated out to here
    double dwCutoff;       // DW cutoff parameter (absorption cutoff for deep imaginary wells)
    double sumMax;          // maximum sum radius for DWBA integration (fm)
    double sumMid;          // middle radius of DWBA sum grid (fm)
    double sumMin;          // minimum radius for DWBA sum (fm)
    double sumDensity;     // number of sum points per unit length (for DWBA grid density)
    double midpointFactor; // multipole amplitude multiplier
    double stepsPerUnit;  // step-size spread (ratio for adaptive step)
    double R;               // current integration radius (fm) — running variable in Numerov loop
    double boundAsy;          // bound-state asymptotic normalization (Whittaker W coefficient)
    double scatAsy;          // scattering asymptotic amplitude
    double phiMid;          // phi_mid: midpoint angle (radians) for phi integration grid

    // Partial-wave (L) grid parameters (integers; zero-init via `= {}`).
    int lStep;             // partial-wave step between computed L values
    int maxLExtrap;        // extra L's to extrapolate beyond lMax
    int lMinSub;           // lMin subtraction below lCrit
    int nCoulombPoints;    // number of Coulomb-integral quadrature points
};
