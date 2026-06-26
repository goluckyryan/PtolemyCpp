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
