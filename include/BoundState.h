#pragma once
// BoundState.h — BoundState class.
//
// Owns the form-factor work array (φ(r)·V(r) product) and the per-vertex
// bound-state data. After allocation, sets vertex[v].jbdPointer to point at the
// 1-based form-factor storage (single source of truth, no sync).
//
// Thread safety: not thread-safe (ptolemy is single-threaded).

#include <vector>

// Forward decl so BoundState method decls can name `Reaction&` as a
// parameter without including Reaction.h (which itself includes this header).
// A `Reaction&` *member* on BoundState would shift boundState's interior layout
// in Reaction and break the FLOAT_arr / aMs-style pointer arithmetic, so the
// reference is threaded per method instead.
class Reaction;

class BoundState {
public:
    // -----------------------------------------------------------------------
    // Vertex — per-vertex bound-state data.
    // 1-based: vertex[1]=projectile, vertex[2]=target; vertex[0] unused.
    // -----------------------------------------------------------------------
    struct Vertex {
        double bsVstep = 0.0;     // step size for BS Numerov integration, fm
        double bsMass    = 0.0;    // reduced mass, AMU
        double massRatio = 0.0;    // mass ratio m_x / m_composite

        // --- Quantum numbers ---
        int    lBound = 0;             // orbital angular momentum L
        int    nodeCount = 0;          // number of radial nodes
        int    jB = 0;                 // 2*J+1

        double rLMax  = 0.0;       // radius of wavefunction maximum, fm
        double vMax   = 0.0;       // potential at wavefunction maximum, MeV
        double boundSp  = 0.0;     // sampling radius, fm
        int    nSpBd  = 0;         // number of sample points
        double alpha  = 0.0;       // decay constant, fm^-1
        double boundMx  = 0.0;     // max radius for BS integration, fm
        double* jbdPointer = nullptr;   // direct pointer for wavefunction

        std::vector<double> wavefunction;   // φ(r) on radial grid, 0-based: [0..nSteps-1]
        std::vector<double> bsPotential;    // V(r) on radial grid, 0-based: [0..nSteps-1]

        // 1-based pointer to data (ptr[1] = first element); matches pool convention.
        double* getWavefunction() {
            return wavefunction.empty() ? nullptr : wavefunction.data() - 1;
        }
        double* getPotential() {
            return bsPotential.empty() ? nullptr : bsPotential.data() - 1;
        }
    };

    // -----------------------------------------------------------------------
    // Data — bound-state search and shape parameters (reaction-wide).
    // -----------------------------------------------------------------------
    struct Data {
        // --- R → RX Jacobi-transformation coefficients ---
        // Set once in DWBAGrid::gridSet from mass ratios; read by
        // BoundState::evaluateFormFactor to compute rP / rT / rCore for the
        // form-factor coordinate change. Stripping vs pickup flips the mapping.
        double s1;
        double s2;
        double t1;
        double t2;
        double scatRMax;       // scattering wavefunction maximum
        double scatInvStep;    // scattering grid: 1/rStep for aitkenLagrange fractional-index lookup
        int    scatPointCount; // number of scattering sample points
        // --- Scattering wavefunction sample tables ---
        // Written by CWF::coulombIntegral (rValue * |PSI|), read by
        // BoundState::evaluateFormFactor (aitkenLagrange interpolation) and by
        // integration_grid (lMin-wave search).
        // Size = nGridSteps + 1 elements; element [0] = legacy slot[1], etc.
        std::vector<double> sctmnArr;   // minimum-L scattering integral
        std::vector<double> sctcrArr;   // critical-L scattering integral
        float  phiSign;   // sign of wavefunction at phi=0 (REAL*4; +1 or -1)
        double rOfMax;    // radius of maximum |wavefunction| found during search (fm)
    };

    // 1-based: vertex[1]=projectile, vertex[2]=target; vertex[0] unused.
    Vertex vertex[3] = {};

    Data data = {};

    std::vector<double> formFactor;

    // Resize formFactor to nSteps+1 doubles (1-based: indices 1..nSteps valid).
    void allocateFormFactor(int nSteps);

    // Returns 1-based pointer (same convention as pool pointers).
    double* getFormFactor() { return formFactor.empty() ? nullptr : formFactor.data() - 1; }

    void solve(int& returnCode, Reaction& reaction);
    void setupFormFactors(Reaction& reaction);
    int  evaluateFormFactor(double& fpFt, int ffKind, double rA, double rB, double cosTheta,
                const double* scatPointer, int aitkenOrder, double& rP, double& rT,
                Reaction& reaction);
    void coulombIntegrals(double rLower, double etaIn, double etaFinal, double fkIn, double fkFinal,
                double sigIn, double sigFinal, double accuracy,
                double& ffIntegral, double& fgIntegral, double& gfIntegral, double& ggIntegral,
                double* points, double* weights, double* fA, double* fpA,
                double* gA, double* gpA, double* work,
                int rPower, int li, int lf, int termCount, int nPts,
                int& returnCode, int printLevel);
    void convertJtoL(int L, int jProj, int channelIndex, double sJr, double sJi,
                int* indxePointer, double* sLx, Reaction& reaction);
};
