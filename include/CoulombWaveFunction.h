#pragma once
// coulomb.h — Coulomb wavefunction and scattering class
//
// Provides a unified interface for all Coulomb-related computations:
// - Regular F_L(η,ρ) and irregular G_L(η,ρ) Coulomb wavefunctions
// - Asymptotic Coulomb phase shifts σ_L
// - Coulomb scattering amplitudes
// - Coulomb integrals for DWBA transfer matrix elements
//
// All methods are static — use CoulombWaveFunction::methodName(...) directly.
// Thin backward-compat free-function wrappers (RCWFN, RCASYM, etc.) are
// provided in coulomb.cpp and delegate to these static methods.

class Reaction;

class CoulombWaveFunction {
public:
    // -----------------------------------------------------------------------
    // computeFG — RCWFN
    // Computes regular (F) and irregular (G) Coulomb wavefunctions for
    // L = minL..maxL using continued fractions, Maclaurin series, and
    // Taylor expansion.  Original: Barnett, Feng, Steed, Goldfarb (1974).
    // -----------------------------------------------------------------------
    static void computeFG(double rho, double eta, int minL, int maxL,
                          double* fCArg, double* fCpArg,
                          double* gCArg, double* gCpArg,
                          double accuracy, int& returnCode);

    // -----------------------------------------------------------------------
    // asymptoticPhase — RCASYM
    // Asymptotic series for Coulomb functions F, G and phase φ.
    // -----------------------------------------------------------------------
    static void asymptoticPhase(int L, double eta, double rho, int printLevel, double sigL,
                                double* zetaPointer, double* phiPointer, double* dZetaPointer,
                                double* fPointer, double* fpPointer, double* gPointer, double* gpPointer,
                                double* z, double* dzSquared, double* s, double* zInv,
                                double eps, int nMax, int& nTz, int& convergenceCode);

    // -----------------------------------------------------------------------
    // coulombIntegral — COULIN
    // Coulomb integrals by recursion: ∫(R→∞) F·F/rⁿ, F·G/rⁿ, etc.
    // -----------------------------------------------------------------------
    static void coulombIntegral(int rPower, int maxDel, int lMin, int lMax,
                                double etaOut, double akOut, double* sigOut,
                                double etaIn,  double akIn,  double* sigIn,
                                double R, int includeIrregularG,
                                double* ff, double* fg, double* gf, double* gg,
                                int lDlDimension,
                                double accuracy, int nPts,
                                double* work, double* fIn, double* fOut,
                                double* gIn, double* gOut, double* starts,
                                int printLevel, int& returnCode, double& clTime,
                                Reaction& reaction);

    // -----------------------------------------------------------------------
    // penetrability — COULNG
    // Coulomb wavefunction for negative energy (Whittaker function).
    // -----------------------------------------------------------------------
    static double penetrability(int L, double eta, double rho, double aNorm);

    // -----------------------------------------------------------------------
    // computeScatteringWaves — GETSCT
    // Computes scattering waves at lMin and lCrit for GRDSET.
    // -----------------------------------------------------------------------
    static void computeScatteringWaves(int& returnCode, Reaction& reaction);

    // -----------------------------------------------------------------------
    // scattering — COULST
    // Coulomb integrals for inelastic DWBA and coupled-channels calculations.
    // -----------------------------------------------------------------------
    static void scattering(int& returnCode, Reaction& reaction);

    // -----------------------------------------------------------------------
    // solveRTXLNX — RTXLNX (4-arg)
    // Finds real solution of a·x + b·ln(x) + c = 0 by Newton iteration.
    // -----------------------------------------------------------------------
    static double solveRTXLNX(double a, double b, double c, double acc);

    // setupFG / generateBasisIndex / setBasisFactors (SETFG/GENBNX/SETBFC) —

};
