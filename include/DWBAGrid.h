#pragma once
// DWBAGrid.h — DWBAGrid class
//
// Centralizes DWBA integration grid array storage for
// gridSet/inelDc/inelasticRadialIntegrals/inelasticGridSet. After allocation,
// updates gridData.*Pointer pointers to point here.
//
// Thread safety: not thread-safe (ptolemy is single-threaded).

class Reaction;

#include <vector>

class DWBAGrid {
public:
    // RIOEX: exponential decay factors exp(alpha_p*r_i + alpha_t*r_o) for each
    // (r_i, r_o) grid point pair. 1-based: valid range [1..nRiRoH].
    std::vector<double> rioEx;

    // SMHPTS: H-integration sum grid points (nPhiSum elements). 1-based.
    std::vector<double> smhpts;

    // SMHWK: H-integration work area (3*nPhiSum: min/mid/max per point). 1-based.
    std::vector<double> smhwk;

    // SMIPTS: Interpolation grid points (NPSUMI elements). 1-based.
    std::vector<double> smipts;

    // SMIVL: Interpolation grid values (NPSUMI elements). 1-based.
    std::vector<double> smivl;

    // SMHVL: H-integral values (nPhiSum * NUMHS elements).
    // 0-based: ptr[0..nPhiSum*NUMHS-1] valid (gridData.smhvlPointer = data()-1).
    std::vector<double> smhvl;

    // IHINT: H integral array (NMLOLX elements).
    // IHABS: H absolute-value integrals (NMLOLX elements). 1-based.
    std::vector<double> hint;
    std::vector<double> habs;

    // LILOR/LILOI: real/imag 3D integrals I1(lx,li). 1-based [1..liloSize].
    std::vector<double> liloR;
    std::vector<double> liloI;

    // IABS1: absorptive term integral (liloSize elements). 1-based [1..liloSize].
    std::vector<double> abs1;

    // iIndex: H/DW index table (int4 packed, 0-based). iiindxPointer=iiindx.data().
    std::vector<int> iiindx;

    // DW: distorted wave products (nWfi*NWFO*2 doubles). dwPointer=dw_.data() (0-based).
    std::vector<double> dw_;

    // DW index table (int4 packed, 0-based).
    std::vector<int> dwi;

    // IWFII: incoming WF index table (int4 packed, 0-based).
    std::vector<int> iwfii;

    // IWFIO: outgoing WF index table (int4 packed, 0-based).
    std::vector<int> iwfio;

    // NLAMBDA: per-L count of lambda-multipole entries (lOutMax+1 ints, 0-based).
    // Filled by A12. Caller passes nlam.data(); A12 does nLam=arg-1, so A12's
    // nLam[1..lOutMax+1] = nlam[0..lOutMax].
    std::vector<int> nlam;

    // XLAMBDA: λ(L,M) angular coefficients (INEED doubles, 0-based, A12 sees
    // 1-based). Caller passes xlam.data(); A12 does xLam=arg-1.
    std::vector<double> xlam;

    // A12TEMPS: A12 scratch buffer. Three logical sub-arrays packed into one
    // std::vector<double>:
    //   [0 .. numLx-1]                   : OUTTMP (REAL*8) — A12 sees 1-based
    //   [numLx .. numLx+2*lxMax+1]       : xLoTemp (REAL*8) — A12 sees 1-based
    //   ints at &a12tm_[numLx+2*lxMax+2] : LXTMP  (INTEGER, 3*numLx values)
    // Sized 4*numLx+2*lxMax+2.
    std::vector<double> a12tm_;

    // rI/rO/WIO: r_i, r_o grid points and integration weights (float, N4RIO
    // each). 0-based; riPointer = ri.data() etc.
    std::vector<float> ri_;   // RIPTS   — incoming r grid points
    std::vector<float> ro_;   // ROPTS   — outgoing r grid points
    std::vector<float> wio_;  // RIROWTS — integration weights

    // LIR/LII: incoming wavefunction real/imag on (r_i,r_o) grid.
    // LOR/LOI: outgoing wavefunction real/imag. 0-based floats.
    std::vector<float> lir;  // ILIR
    std::vector<float> lii;  // ILII
    std::vector<float> lor_; // ILOR (trailing _ to avoid log conflict)
    std::vector<float> loi;  // ILOI

    // Phi-angle arrays (RIRO integration grid).
    // Accessed 1-based as float: phiT[0] unused, phiT[1..phiCount] valid.
    std::vector<float> phiT;   // PHIANG1 — theta angles
    std::vector<float> phiP;   // PHIANG2 — phi angles
    std::vector<float> phi_;   // PHIANG  — phi values
    std::vector<float> trapWeight_;   // FIFODX — trapezoidal rule weights

    // cos table (ICOSSP=NCOSIN+1 elements, 0-based).
    std::vector<double> cosin;

    // Allocators (set gridData.*Pointer caches).
    void allocateRioEx(int size, Reaction& reaction);
    void allocateSmhpts(int nPhiSum, Reaction& reaction);
    void allocateSmhwk(int nPhiSum, Reaction& reaction);    // size = 3*nPhiSum
    void allocateSmipts(int nPhiSumI, Reaction& reaction);
    void allocateSmivl(int nPhiSumI, Reaction& reaction);   // gridSet version
    void allocateSmhvl(int nPhiSumHCount, Reaction& reaction); // angleSet version: nPhiSum * hCount
    void allocateHint(int nMloLx, Reaction& reaction);    // IHINT — H integral values
    void allocateHabs(int nMloLx, Reaction& reaction);    // IHABS — H absolute-value integrals
    void allocateLiloR(int size, Reaction& reaction);     // LILOR — real 3D integrals
    void allocateLiloI(int size, Reaction& reaction);     // LILOI — imag 3D integrals
    void allocateAbs1(int size, Reaction& reaction);      // IABS1 — absorptive term integral
    void allocateIiindx(int n, Reaction& reaction);        // H/DW index (int4, n doubles → FACFR4*n ints)
    void allocateDw(int n, Reaction& reaction);             // distorted wave products
    void allocateDwi(int n, Reaction& reaction);            // DW index (int4, n doubles → FACFR4*n ints)
    void allocateIwfii(int n, Reaction& reaction);          // IWFII  — incoming WF index (int4, n doubles)
    void allocateIwfio(int n, Reaction& reaction);          // IWFIO  — outgoing WF index (int4, n doubles)
    void allocateNlam(int size);        // NLAMBDA — per-L lambda counts (size ints)
    void allocateXlam(int size);        // XLAMBDA — λ(L,M) angular coefficients (size doubles)
    void resizeIwfio(int n, Reaction& reaction);            // IWFIO  — shrink
    void allocateLir(int nRiRoH, Reaction& reaction);     // ILIR — incoming wavefunction real (2*nRiRoH floats)
    void allocateLii(int nRiRoH, Reaction& reaction);     // ILII — incoming wavefunction imag
    void allocateLor(int nRiRoH, Reaction& reaction);     // ILOR — outgoing wavefunction real
    void allocateLoi(int nRiRoH, Reaction& reaction);     // ILOI — outgoing wavefunction imag
    void allocateCosin(int cosinSize, Reaction& reaction);   // cos table (0-based)
    void allocatePhiArrays(int jSize, Reaction& reaction);  // PHIANG1/2/PHIANG/FIFODX (FACFR4=2 stride baked in)
    void allocateRiRoWio(int n4rio, Reaction& reaction);    // r_i/r_o/weight float arrays (FACFR4=2 stride baked in)

    void gridSet(int& returnCode, Reaction& reaction);
    void inelDc(Reaction& reaction);
    void inelasticRadialIntegrals(Reaction& reaction);
    void inelasticGridSet(int& returnCode, Reaction& reaction);
};
