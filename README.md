# PtolemyCpp — Modern C++ Ptolemy

PtolemyCpp is a modern C++ port of the Ptolemy DWBA nuclear-reaction code. It
produces bit-identical angular differential cross sections to the original
Fortran Ptolemy while running ~4.8× faster, with a clean object-oriented
architecture.

## History — Standing on shoulders

- **Original Fortran Ptolemy** by **Stephen C. Pieper** and collaborators at
  Argonne National Laboratory (development beginning ~1977, ~40k LOC accumulated
  over decades).
- **Faithful 1:1 Fortran-to-C++ translation: Ptolemy-f2c** by Ryan Tang — a
  1:1 transliteration preserving every COMMON block and GOTO.
  Repo: https://github.com/goluckyryan/Ptolemy-f2c
- **This project (PtolemyCpp)** — a deep refactor of Ptolemy-f2c into modern C++17,
  eliminating all globals, all GOTOs, all COMMON blocks, and the NALLOC memory
  pool. Bit-identical to Cleopatra-verified Maple (where Maple is the
  bit-identical-to-32-bit-Cleopatra reference build of the original Fortran).

## Quick start

```bash
make                    # builds ./ptolemy
./ptolemy < test_inputs/transfer_12C_dp.in    # run one reaction
./test_ptolemy.sh       # full 35-test bit-identical regression vs Maple
make test               # builds + runs the 58 unit tests
make distclean          # clean build + binaries
```

Requires `g++` with C++17 support. No external dependencies.

The regression suite compares PtolemyCpp output line-by-line to a reference
Maple build (the 1:1 Fortran-faithful predecessor). To override the Maple
path or skip the diff when unavailable:

```bash
MAPLE=/path/to/your/ptolemy ./test_ptolemy.sh    # custom reference
MAPLE=missing ./test_ptolemy.sh                  # PtolemyCpp-only, skip diff
```

## Stats

- **28,668 LOC** across 54 `.cpp` / 49 `.h` (down from 38k Fortran baseline, −25%)
- **35/35** reaction tests bit-identical to Cleopatra (via Maple oracle)
- **58/58** unit tests pass
- **~0.89s** wall vs Maple ~4.2s on the full 35-test suite at `-O2` (~4.7× faster)

## Architecture

One `Reaction` instance owned by `main()` holds all physics state, composed
of purpose-named sub-structs:

```
class Reaction {
    Kinematics       kin;             // momenta, eta, LCRIT
    Energies         energies;        // E, Q, excitation
    Masses           masses;          // projectile/target/residual
    Charges          charges;         // Z values
    AngularMomentum  angMom;          // J, L, parities
    OpticalPotential opticalPotential; // V(r) on the radial grid
    BoundState       boundState;      // vertex wavefunctions
    DistortedWave    distortedWave;   // per-channel scattering wavefunctions
    DWBAGrid         dwbaGrid;        // radial-integral buffers
    InelasticData    inelastic;       // S-matrix scratch, λ-coefficients
    // ... 6 more focused sub-structs ...
};
```

Three orchestrators drive the calculation, each taking `Reaction&`:

- `Elastic`     — elastic scattering, prints DCS / analyzing powers
- `Transfer`    — DWBA transfer, computes radial integrals + cross sections
- `Inelastic`   — inelastic / collective DWBA

The **`OpticalPotential`** class is the extension point for new V(r) forms:

```cpp
class OpticalPotential {
    std::vector<double> values;       // grid values, 0-based
    void resize(int nPts, double rStart, double step);
    void fillWoodsSaxon(double V, double R, double a);
    void fillSpinOrbit (double V, double R, double a);
    void fillSurface   (double V, double R, double a);
    void fillCoulomb   (int channelIdx, int channelIdxIn);
    void add(const OpticalPotential&);   // composition
    void scale(double k);                  // e.g. Numerov −h²/12E
};
```

To add a new potential form (CD-Bonn, Koning-Delaroche, Becchetti-Greenlees,
custom): add one method to `OpticalPotential`. No factory, no inheritance,
no dispatch.

Pure math (Clebsch-Gordan, 3j/6j/9j, Legendre, Gauss-Legendre quadrature,
Numerov integration) lives in `MathAngular` / `MathFunctions` / `NumerovSolver`
with no state.

Plugin-style bound-state solvers (AV18 deuteron interior, Phiffer Weinberg
expansion) live behind the `LinkulePlugin` interface.

## Status

Stable. Bit-identical regression-tested. Open for contribution.

### Known limitations — collective inelastic model has a narrow domain

PtolemyCpp's inelastic scattering uses the **vibrational / collective
excitation model** inherited from the original Fortran Ptolemy: the
excitation is treated as a one-step surface deformation

  V(r, θ) = V₀(r) - β_L R₀ ∂V/∂R · Y_LM(θ)

parametrized by a single deformation amplitude `BELX` = β_L (input from
`PARAMETERSET INELOCA*`). This is the right physics for low-lying
collective 2⁺ / 3⁻ states of **even-even spherical heavy nuclei**
(e.g. 90Zr, 208Pb, Sn isotopes, 40Ca) at moderate excitation. It is
**physically inappropriate** for:

- **Light targets** (12C, 16O, 20Ne, ...) where the excited states are
  single-particle or cluster in character, not collective surface modes.
- **Excitation energies above the particle-emission threshold** (~7 MeV
  for 12C, ~12 MeV for 16O) where the daughter nucleus is unbound and
  the DWBA assumption of a stable bound final state is violated.
- **Large deformation amplitudes** outside the small-β linearization the
  one-step model assumes.

When the model is applied outside its domain — most visibly in
`(d,d')` on light targets with E_x ≳ 4 MeV at 20–60 MeV beam energy —
the code fails numerically downstream of the physics mismatch. The
L-extrapolation breaks (`**** ERROR: CANNOT EXTRAPOLATE FOR CHANNEL`),
|S(L)| diverges at high L, and the reported DCS becomes meaningless.
All four codes in the lineage exhibit this in different ways:

| Code | DCS @ 0° for 12C(d,d') 2⁺ 5MeV, Ed=30 | Behavior past failure |
|---|---|---|
| Cleopatra (Fortran ifort 32-bit) | 9 × 10¹⁶ mb/sr | extrapolation kept → garbage |
| Maple (Cleopatra-faithful C++) | 9 × 10¹⁶ mb/sr | matches Cleopatra bit-identically |
| Ptolemy-f2c (1:1 translation) | 9 × 10¹⁶ mb/sr | matches Cleopatra bit-identically |
| **PtolemyCpp** (this code) | 35 mb/sr | partial sum, also wrong |

For a random sweep of 800 reactions, the affected band was ~60 cases,
all in the high-E_x light-target `(d,d')` corner described above. The
DCS column flag `% FROM L>LMAX` is the canonical diagnostic: any value
far from 0.0 means the answer is dominated by extrapolated
(non-converged) partial waves and should not be trusted. Increasing
`LMAX` / `LMAXADD` / `ASYMPTOPIA` typically makes things worse, not
better, because the broken extrapolation simply runs longer.

**For collective inelastic scattering outside the model's domain** —
use FRESCO or ECIS with an explicit microscopic form factor, or
compare against experimental data.

**Reactions inside the model's domain are unaffected.** The 35-test
regression suite (PtolemyCpp vs Cleopatra) passes 35/35 bit-identical;
the 800-case random sweep shows 294/300 transfer + 199/200 elastic +
238/300 inelastic PASS, with all 62 inelastic FAILs concentrated in the
problematic regime.

## Contributing

Run `./test_ptolemy.sh` to verify any change is bit-identical against the
Maple oracle. Run `make test` for the unit tests. PRs welcome.

The single hard rule: **commits must preserve bit-identical 35/35 + 58/58.**
The exp recurrences in `OpticalPotential::fillWoodsSaxon` etc. are not
equivalent to per-point `std::exp()`, so refactors there require care.

## License

Derived from the original Argonne Ptolemy code. No formal license file at
this time — please contact the maintainer for use beyond personal/research
purposes.
