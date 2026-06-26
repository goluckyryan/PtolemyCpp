#pragma once
// InelasticData.h — split out from Reaction so Reaction can hold an
// InelasticData member.

#include <vector>

// ============================================================================
// Inelastic and collective DWBA scattering control
// ============================================================================
struct InelasticData {
    int    lxMin;   // minimum transferred L_x for inelastic sum
    int    lxMax;   // maximum transferred L_x for inelastic sum
    int    mStop;   // 0=continue; nonzero=stop after current L (convergence reached)
    int    a12mSize;   // size of the A12 * M (multipole) product array storage
    int    nMloLx;  // number of (L_i, L_o, L_x) triplets = numLx * lOutMax^2
    int    numLx;   // number of L_x values in sum
    // S reaction S-matrix scratch. 0-based float vector of size 2*liloSize
    // (interleaved real/imag pairs). Written by interpolation, read by
    // phase_shift_print / xSection:non-CC.
    std::vector<float> sMatrixArr;
    int    liloSize;  // size of L_i * L_o block (for index arithmetic)
    int    nLValues;  // number of L_i values in incoming-channel sum
    // 1-based: indices [1..nLValues] valid; [0] unused.
    std::vector<int> lisArr;
    // 0-based [0..2]: pointers to the beta deformation-parameter arrays
    // (BETA / BETACOUL / BETARATS) owned by reaction.named. Set by probe_print;
    // read by input_reader, interpolation, CoulombWaveFunction_scattering.
    // nullptr = "no data".
    std::vector<double>* poolBetas[3] = {nullptr, nullptr, nullptr};
    // ATERM (transfer mode). Size lxMax+1, 0-based; atermArr[LXv] for LXv in
    // [lxMin..lxMax]; passed to sFromI as atermArr.data() and read 0-based as aTerm[LXP].
    std::vector<double> atermArr;
    // BETANRAT. Size numLxI = (lxMax - lxMin)/2 + 1, 0-based; element k =
    // (LXv - lxMin)/2. Writer (probe_print) and reader (input_reader /
    // DWBAGrid::inelasticRadialIntegrals) both use (betnrArr.data() - lxMin/2)[LXv/2].
    std::vector<double> betnrArr;
    // FF2INTS — Coulomb-to-form-factor coupling array. Sized
    // nSpl*(lInMax-lMin+1) doubles in CWF_scattering coulombInel.
    std::vector<double> cl2ffArr;
    int    lxStep;  // step in L_x (usually 1; 2 for even parity)
    int    densitySwitch;  // LOGICAL: 1=use density-dependent (Bohr-Mottelson) form factor
    int    lSkip;   // number of L values to skip at low end of sum
    // --- J projection loop bounds for coupled channels ---
    int    jpMin;   // minimum m_p (projectile J projection) in CC sum
    int    jpMax;   // maximum m_p
    int    jpBase;  // base offset for m_p index into CC matrix
    int    nJp;     // number of m_p values = jpMax - jpMin + 1
    int    jtMin;   // minimum m_t (target J projection)
    int    jtMax;   // maximum m_t
    int    jtBase;  // base offset for m_t index
    int    nJt;     // number of m_t values
    // --- S-matrix and TOC ---
    // INDXS — L-index array for S-matrix storage, produced by SETSPT pass 2.
    // 0-based via indxsPointer = indxsArr.data(). Sized ISIZE entries
    // (3*nLx*nJp*nJt summed over CC channels).
    std::vector<int> indxsArr;
    int*   indxsPointer = nullptr;  // 0-based: indxsPointer[K] for K in [0..ISIZE-1]
    // SMATR/SMATI — real/imag parts of the inelastic S-matrix. Sized nSmatPerL
    // (nAspli*nLValues doubles) + 1-pad. 0-based via .data(); callers subscript
    // [i-1]; interpolation uses .data()+(kOffset-1) directly.
    std::vector<double> smatrArr;
    std::vector<double> smatiArr;
    double* smatRPointer = nullptr;  // 0-based: smatRPointer[i-1] for i in [1..nSmatPerL]
    double* smatIPointer = nullptr;
    int    nSmatPerL;   // total size of S-matrix storage
    int    nLx;     // number of L_x values with nonzero contribution
    int    nSpl;   // number of amplitude splittings (spin projections)
    // TOCS — 4 ints per amplitude slot (lDelta, lx, jProj, jT) produced by
    // SETSPT pass 3. Sized 4*nAspli + pad. 1-based via tocsPointer =
    // tocsArr.data() - 1; CC adds per-channel LTOCOF offset. -1 sentinel marks
    // skipped (kOffset, ITOC[4*kOffset] < 0).
    std::vector<int> tocsArr;
    int*   tocsPointer = nullptr;  // 1-based: tocsPointer[4*kOffset + j] for j in [-3..0]
    int    nAspli;  // number of (m_p, m_t) angular-momentum amplitude pairs
    // 0-based, indexed by [li - lMin] for li=lMin..jMost.
    std::vector<float> unitrArr;
    // sMag/sPhase reaction-channel S-matrix magnitudes/phases. Written by
    // CrossSectionCalc::phasePrint for chanNumber=1 (reaction), consumed by
    // xSection. 0-based float vector of size liloSize. Index layout:
    // (li-lBase)*nSpl + kOffset, accessed via smagArr.data() + LSOFF
    // (LSOFF=0 for non-CC, nonzero for CC).
    std::vector<float> smagArr;
    std::vector<float> sphaseArr;
    // Non-CC Coulomb integrals (R→∞) for the FF/FG/GF/GG products. Written by
    // CoulombWaveFunction::scattering, consumed by
    // DWBAGrid::inelasticRadialIntegrals (input_reader). Sized nSmatPerL;
    // 0-based; index = nSpl*(liIndex-1) + kOffset - 1.
    std::vector<double> cl1ffArr;
    std::vector<double> cl1fgArr;
    std::vector<double> cl1gfArr;
    std::vector<double> cl1ggArr;
    int    liFit;   // L_i_max for fit (last L used in L-interpolation)
    int    liFitIndex;  // index of liFit in the L_i list
    int    nolFit;  // number of L_i values above liFit (not fitted)
    // Direct pointer fields for class-owned arrays.
    double* liloRPointer = nullptr;  // real part of 3D integral (0-based, class-owned)
    double* liloIPointer = nullptr;  // imag part of 3D integral (0-based, class-owned)
    // r2s[1..4] — inelastic form-factor radii + Coulomb VC.
    //   r2s[1] = real radius    (= R0  * r2Mass)
    //   r2s[2] = imag radius    (= rI0 * r2Mass)
    //   r2s[3] = Coulomb radius (= rC0 * r2Mass)
    //   r2s[4] = Coulomb VC value (= -3 Z1 Z2 hbar_c / fine_structure_inv)
    double r2s[5] = {0,0,0,0,0};  // 1-based, [0] unused
};
