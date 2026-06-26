#pragma once
// NumerovSolver.h — Numerov method radial Schrödinger equation solver
//
// Solves: u''(r) + [k² - V_eff(r)] u(r) = 0
// where u(r) = r * R(r) is the reduced radial wavefunction.
//
// Shared by BoundState and DistortedWave (composition).
// Single copy of the Numerov integration code.
//
// Aligned with Raphael: ~/PtolemyGUI/Raphael/solveSE.py
// Note: Raphael uses RK4, Ptolemy uses Numerov (3-point recursion).
//
// Currently delegates to existing solve()/wavelj() implementations.
// The methods below provide a clean interface for future direct implementation.

#include <vector>

class NumerovSolver {
public:
    // -----------------------------------------------------------------------
    // Grid configuration
    // -----------------------------------------------------------------------
    double stepSize  = 0.0;     // integration step h (fm)

    // -----------------------------------------------------------------------
    // Results
    // -----------------------------------------------------------------------
    std::vector<double> solution;   // u(r) on grid [0..N-1]
    int    nodesFound = 0;          // number of nodes in solution

    // -----------------------------------------------------------------------
    // Core Numerov integration: solve u''(r) = f(r) * u(r)
    //
    // Uses the three-point recursion:
    //   u_{n+1} = (2 - 5h²f_n/6) * u_n - (1 + h²f_{n-1}/12) * u_{n-1}
    //             / (1 + h²f_{n+1}/12)
    //
    // f(r) = V_eff(r) - k² (the effective potential minus the energy)
    // -----------------------------------------------------------------------

    // Integrate outward from rStart
    // fOnGrid: precomputed f(r) values [0..nSteps-1]
    // u0, u1: starting values for u[0] and u[1]
    void integrateOutward(const std::vector<double>& fOnGrid,
                          double u0, double u1);

    // Count nodes in solution
    int countNodes() const;
};
