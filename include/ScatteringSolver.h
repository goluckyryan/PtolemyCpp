#pragma once
// scattering_solver.h — ScatteringSolver class
// Encapsulates wavefunction work arrays + local potential / spin-orbit /
// tensor vectors that used to live in the Fortran pool.
//
// Design:
//   - Owns wavefunction work arrays (wavR, wavI) as std::vector<double>
//   - reallocate(size) sets wavRPointer/wavIPointer to point at the vector data
//   - Callers go through wavRPointer / wavIPointer (1-based pointers)
//
// Usage:
//   // Callers use wavRPointer / wavIPointer (1-based pointers) directly.
//
// Thread safety: not thread-safe (ptolemy is single-threaded).

#include <vector>

class Reaction;

class ScatteringSolver {
public:
    // Per-radius scratch values written by MAKPOT and read by Numerov inner
    // loop. Lives on ScatteringSolver because both writer (setupScatteringWaves)
    // and readers (Numerov, integration_grid, source_potentials, probe_print,
    // linkulesfitters_potentials) operate on a single shared instance per call
    struct PotentialWork {
        double tvReal;    // temp real central potential V_R at current r (MeV)
        double tvImag;    // temp imag central potential V_I at current r (MeV)
        double taReal;    // temp real potential asymptotic amplitude (MeV)
        double taImag;    // temp imag potential asymptotic amplitude (MeV)
    };

    // -----------------------------------------------------------------------
    // Wavefunction work arrays (1-based: data()[0] unused, data()[1..N] valid)
    // These replace the pool-allocated WAVER/WAVEI arrays.
    // For tensor-coupled channels, size = 4*(nStep+6); otherwise (nStep+6).
    // -----------------------------------------------------------------------
    std::vector<double> wavR;   // Re(chi) work array
    std::vector<double> wavI;   // Im(chi) work array

    // -----------------------------------------------------------------------
    // 1-based: data()[0] unused, data()[1..N] valid.
    // -----------------------------------------------------------------------
    std::vector<double> vWork;     // work array for MAKPOT

    // class-owned spin-orbit and tensor potential vectors.
    // setupScatteringWaves sets P_ISORS[k]=sors[k].data()-1 (1-based double*);
    // solvePartialWave uses P_ISORS[channelIndex] != nullptr as the allocation guard.
    std::vector<double> sors[3];   // spin-orbit real [per channel]
    std::vector<double> sois[3];   // spin-orbit imag [per channel]
    // callers after the tensor block (VTR..VTPI etc.) was removed.


    // -----------------------------------------------------------------------
    // Solver-wide step-size state (not per-channel — both channels share
    // these during the integration of a single partial wave).
    // -----------------------------------------------------------------------
    double stepI = 0.0;   // imag wavefunction step size (fm)
    int    isStandalone = 0;    // LOGICAL: 1=use standard (non-adaptive) step
    int    nFirst = 0;    // L index of first partial wave computed
    int    pwBgSwitch = 0;    // LOGICAL: 1=plane-wave background subtracted
    int    pwAvSwitch = 0;    // LOGICAL: 1=partial-wave (plane wave basis)
    int    lastZero = 0;    // L index of last zero-crossing
    int    hasAnySpinorbit  = 0;    // LOGICAL: 1=spin-orbit weight correction
    int    nWaveF = 0;    // number of wavefunction files written so far

    // Per-radius MAKPOT scratch.
    PotentialWork potentialWork = {};

    // -----------------------------------------------------------------------
    // Pointer caches for wavefunction work arrays. 1-based: data()-1 so
    // wavRPointer[1] == wavR[1]. Set by reallocate() from wavR.data() / wavI.data().
    // -----------------------------------------------------------------------
    double* wavRPointer = nullptr;   // 1-based wavefn-real ptr
    double* wavIPointer = nullptr;   // 1-based wavefn-imag ptr
    double* vWorkPointer = nullptr;  // 1-based vWork ptr

    // -----------------------------------------------------------------------
    // reallocate(size)
    //
    // Resizes wavR/wavI to `size` doubles (1-based: indices 1..size-1 valid)
    // and updates wavRPointer / wavIPointer to point at the vector data.
    // -----------------------------------------------------------------------
    void reallocate(int size);

    // -----------------------------------------------------------------------
    // allocators for class-owned scratch vectors.
    // k = channel index (1 or 2). n = number of elements (1-based: [1..n]).
    // -----------------------------------------------------------------------
    void allocateVWork(int n);         // vWork
    void allocateSors(int k, int n, Reaction& reaction);    // sors[k] (when V0SORS!=0)
    void allocateSois(int k, int n, Reaction& reaction);    // sois[k] (when V0SOIS!=0)


    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------
    void setupScatteringWaves(int& returnCode, int isStandalone, Reaction& reaction);
    void solvePartialWave(int L, int jProj, int channelIndex, int nPts, float* rGrid, float* waveR,
                float* waveI, double* waveReal, double* waveImag, double* vReal,
                double* vImag, double* vCent,
                Reaction& reaction);
};

// ScatteringSolver is owned by reaction.distortedWave (composition).
// Access via reaction.distortedWave.scatteringSolver.

// Scattering free-function helpers (Reaction fwd-declared above).
void dispatchPartialWave(Reaction& reaction, int L, int jProj, int channelIndex, int nPts, float* rGrid,
           float* waveR, float* waveI,
           double* waveReal, double* waveImag, double* vReal, double* vImag,
           double* vCent);
void lCritL(double fK, double eta, double rC, int nSteps,
            double stepSize, const double* vInPointer, const double* wInPointer,
            int printLevel, int& lC1, int& lC2, int& lC);
