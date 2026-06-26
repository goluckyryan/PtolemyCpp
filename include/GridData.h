#pragma once
// GridData.h — split out from Reaction so Reaction can hold a GridData member.

#include <vector>

// ============================================================================
// DWBA integration grid and phi-summation control
// ============================================================================
struct GridData {
    double jacobian;    // Jacobian of the (r_i, r_o, phi) → (r, s, phi) coordinate transformation
    double stepInverse;   // 1/stepSize for the sum grid (inverse of sum step in fm^-1)
    double cosStep;   // cos(theta) step: cos-spacing for the phi-integration Gaussian quadrature
    int    intsOffset;   // offset into intsArr for current channel
    int    noFlo;    // number of floating-point operations counted (profiling)
    int    nRiRoInterp;   // total number of (r_i, r_o) grid pairs
    int    nCrit;    // critical grid index (boundary between nuclear and Coulomb)
    int    a12nSize;    // size of the A12 * N product array storage
    int    maxCount;   // maximum integration count (grid points per r_i shell)
    int    cosinQuarter;   // nCosin / 4 — quarter-period index offset for sin-from-cos table lookup (cos(θ - π/2) = sin θ)
    int    nInterpPoints;   // number of phi-sum integration points (actual count used)
    int    nWfi;   // number of incoming wavefunction grid points stored
    int    nWfo;   // number of outgoing wavefunction grid points stored

    // Phi-integration grid sizes (input/default counts).
    int    nPhiSum;          // sum-direction phi points
    int    nPhiDifference;   // difference-direction phi points
    int    nPhiPoints;       // phi-quadrature points per shell
    int    nPhiAdditional;   // extra phi points appended per shell

    // Class-owned pointer fields: set inside DWBAGrid::allocate*; used
    // throughout grid/angular code as the canonical access pointer.
    double* smhptsPointer = nullptr;   // SMHPTS — sum grid points (1-based)
    double* smhwkPointer  = nullptr;   // SMHWORK — sum grid work area (1-based)
    double* smiptsPointer = nullptr;   // SMIPTS — interp grid points (1-based)
    double* smivlPointer  = nullptr;   // SMIVL — interp grid values (1-based)
    double* dwPointer    = nullptr;   // DW — distorted wave products (0-based double*)
    double* smhvlPointer  = nullptr;   // SMHVL — H-integral values (1-based)
    double* rioExPointer  = nullptr;   // RIROEXPS — exponential factors (0-based)
    double* hintPointer  = nullptr;   // HINTEGRL — H integral values (0-based, accessed [hIndex-1])
    double* habsPointer  = nullptr;   // HABSINT — H absolute-value integrals (0-based, accessed [hIndex-1])
    double* abs1Pointer  = nullptr;   // RIROABS — absorptive term integral (0-based)
    float*  riPointer    = nullptr;   // RIPTS — r_i grid (float, 0-based)
    float*  roPointer    = nullptr;   // ROPTS — r_o grid (float, 0-based)
    float*  wioPointer   = nullptr;   // RIROWTS — r_i*r_o weights (float, 0-based)
    float*  lirPointer   = nullptr;   // WAVEAR — incoming wavefunction real (float, 0-based)
    float*  liiPointer   = nullptr;   // WAVEAI — incoming wavefunction imag (float, 0-based)
    float*  lorPointer   = nullptr;   // WAVEBR — outgoing wavefunction real (float, 0-based)
    float*  loiPointer   = nullptr;   // WAVEBI — outgoing wavefunction imag (float, 0-based)
    double* cosinPointer = nullptr;   // COSINES — cos table (0-based)
    float*  phiTPointer  = nullptr;   // PHIANG1 — phi-angle theta array (0-based)
    float*  phiPPointer  = nullptr;   // PHIANG2 — phi-angle phi array (0-based)
    float*  phiPointer   = nullptr;   // PHIANG — phi-angle array (0-based)
    float*  trapWeightPointer = nullptr;   // FIFODX — trapezoidal integration points (0-based)
    int*    iDwfiPointer = nullptr;   // incoming WF index table (0-based)
    int*    iDwfoPointer = nullptr;   // outgoing WF index table (0-based)
    int*    dwiPointer   = nullptr;   // indxDw — DW index table (0-based)
    int*    iiindxPointer = nullptr;   // iIndex — H/DW index table (0-based)

    // INTEGERS table: doubles holding successive integer values for fast
    // A12-side index→double conversion.
    std::vector<double> intsArr;

    // RPTS4 — REAL*4 Gauss points for WAVELJ interpolator. rpts4Pointer is the
    // 1-based access pointer (rpts4Pointer[1] = first elt).
    std::vector<float> rpts4Arr;
    float*  rpts4Pointer = nullptr;

    // RPTS / RWTS — Gauss-quadrature points/weights for the inelastic radial
    // integral. Reuse smiptsPointer/smivlPointer 1-based caches.
    std::vector<double> rptsArr;
    std::vector<double> rwtsArr;
    // NUCH/COULH — H-array (V * R * Gauss weights) for the inelastic radial
    // integral. Readers use .data() - 2 / .data() - 1 (1-based double view).
    // Size: 2*NUMPT doubles (NUCH: NUMPT complex pairs), NUMPT doubles (COULH).
    std::vector<double> nuclearHArr;
    std::vector<double> coulombHArr;
    // SUMIVALS. Sized nInterpPoints * NUMHS doubles. Writer stores a 0-based
    // offset into iIndex(1, numIi); reader reads
    // smivlArr[iIndex(1, numIi) + uIndex - 1] (the legacy "+iU" is 1-based).
    std::vector<double> smivlArr;

    // A12MSVAL / JA12S — msval doubles and JA12S ints walked by the kA12m loop
    // in inelastic_dwba's H sum.
    std::vector<double> msvalArr;
    std::vector<int>    ja12sArr;

    // H + A12VL unified into one vector. hsA12Arr holds H at indices
    // [0..NMLOLX-1] and A12 at [NMLOLX..NMLOLX+a12nSize-1], so the absolute
    // pool offsets stored in JA12S address one contiguous range:
    // [lA12Of+LH]=NMLOLX+jA12On reproduces the legacy pool arithmetic.
    std::vector<double> hsA12Arr;
};
