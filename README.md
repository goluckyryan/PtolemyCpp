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
./test_ptolemy.sh       # full 36-test bit-identical regression vs Maple
make test               # builds + runs the 117 unit tests
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

## DWBA input format

Besides a full Ptolemy input deck, `ptolemy` accepts a compact, human-readable
**DWBA reaction description** — the workflow formerly handled by Cleopatra's
standalone `InFileCreator`, now built in. One reaction per line:

```
target(in,out)residual   gs-Jpi   orbital   Jpi(Ex)   Ex   ELab   Potentials [beta]
```

```bash
echo '208Pb(d,p)209Pb  0+  0g9/2  9/2+  0.000  7MeV/u  AK' | ./ptolemy
./ptolemy test_inputs/dwba/transfer_208Pb_dp.dwba
./ptolemy --dwba my_reactions.txt        # force DWBA mode
```

On a DWBA input, `ptolemy` prints `DWBA input detected, expanding...` to stderr,
expands the description into a full Ptolemy deck (optical-model potentials,
bound-state and projectile vertices, angular grid), and runs it.

| Field        | Meaning                                                                 |
|--------------|------------------------------------------------------------------------|
| `target(in,out)residual` | e.g. `208Pb(d,p)209Pb`; light particles: `n p d t 3He a` |
| `gs-Jpi`     | ground-state spin-parity of the target, e.g. `0+`, `3/2+`              |
| `orbital`    | `none` (elastic/inelastic); `1g9/2` = node-l-j for single-nucleon transfer; `0L=2` = node + transferred L for two-nucleon transfer |
| `Jpi(Ex)`    | spin-parity of the populated state, e.g. `9/2+`                        |
| `Ex`         | excitation energy [MeV]                                                |
| `ELab`       | beam energy: `30MeV` (total) or `7.39MeV/u` (per nucleon)             |
| `Potentials` | 1 code for elastic/inelastic, 2 codes (incoming,outgoing) for transfer |
| `beta`       | optional deformation length for inelastic scattering                   |

Lines beginning with `#` and lines shorter than 5 characters are ignored.
Parity and angular-momentum consistency (`σ_gs·σ_state = (−1)^l`, triangle rule)
are checked; inconsistent lines are skipped with a diagnostic.

**Detection.** A file argument is auto-detected by content (handles leading
comments and an optional leading `DWBA` keyword line). On stdin, a line starting
with a mass-number digit (`206Hg(...`) is treated as DWBA and any other start as
a native deck; for a piped DWBA stream with leading comments, use `--dwba`.

**Optical-potential codes** (`OpticalPotentialLibrary`, ported 1:1 from
Cleopatra/Kay's `globals_beta_v5`):

| Code | Reference | | Code | Reference |
|------|-----------|-|------|-----------|
| `A` | An & Cai (2006), d        | | `K` | Koning & Delaroche (2009), p |
| `H` | Han, Shi, Shen (2006), d  | | `V` | Varner CH89 (1991), p |
| `B` | Bojowald (1988), d        | | `M` | Menet (1971), p |
| `D` | Daehnick rel. (1980), d   | | `G` | Becchetti & Greenlees (1969), p |
| `C` | Daehnick non-rel. (1980), d (n/i) | | `P` | Perey (1963), p |
| `L` | Lohr & Haeberli (1974), d | | `x` | Xu et al. (2011), 3He |
| `Q` | Perey & Perey (1963), d   | | `l` | Liang, Li, Cai (2009), 3He |
| `Z` | Zhang, Pang, Lou (2016), d (6,7Li) | | `p` | Pang et al. (2009), 3He/t |
| `s` | Su & Han (2015), α        | | `c` | Li, Liang, Cai (2007), t |
| `a` | Avrigeanu et al. (2009), α | | `t` | Trost et al. (1987), t |
| `f` | Bassani & Picard (1969), α (fixed) | | `h` | Hyakutake et al. (1980), t |
| `X` `Y` | Bardayan (2008), custom | | `b` | Becchetti & Greenlees (1971), 3He/t |

(Code `C` is not implemented in the library yet.) `potentialRef(code)` returns the
full citation; `callPotential(code, A, Z, E, Zproj)` returns the 16 Woods-Saxon
parameters as an `OMPset`. The library is verified bit-identical to the original
(see the `OpticalPotentialLibrary` unit tests).

Examples live in `test_inputs/dwba/`. Single-nucleon `(d,p)` transfer, elastic,
and two-nucleon `(t,p)` transfer run end-to-end. Two cases expand faithfully but
do not yet run in PtolemyCpp's engine: collective `(p,p')`/`(d,d')` **inelastic**
(the engine keys off `BELX` = B(EL), not the Cleopatra `BETA` deformation length
— see [Known limitations](#known-limitations--collective-inelastic-model-has-a-narrow-domain)),
and A=3 single-nucleon transfer using the `phiffer` projectile fit (e.g.
`(d,3He)`, `(d,t)`); `(d,p)`-type `av18` transfers are unaffected.

## Stats

- **28,668 LOC** across 54 `.cpp` / 49 `.h` (down from 38k Fortran baseline, −25%)
- **36/36** reaction tests bit-identical to Cleopatra (via Maple oracle)
- **117/117** unit tests pass
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

### Spin-orbit convention — `--fixedLS` flag

Ptolemy uses a non-standard spin-orbit coupling factor inherited from
Cleopatra:

```
sDotL = [J(J+1) - L(L+1) - S(S+1)] / (2S)     ← default
```

For **spin-1/2 projectiles (p, n, t)** this equals `2*<L*S>`, which when
combined with `woodsX`'s built-in factor of 2 in the radial form gives
the standard `Vso * <L*S> * (1/r) dV/dr` physics. So proton optical-model
tables (Koning-Delaroche, Becchetti-Greenlees, etc.) work without
adjustment — they were calibrated against this convention.

For **spin > 1/2 projectiles (deuteron, 6Li, 7Li, alpha bound states)**
the `2S` divisor gives the wrong scaling. The effective spin-orbit
strength is `Vso / S` instead of `Vso`. To use spin > 1/2 optical-model
parameters as-published, the simple workaround is to multiply input `Vso`
by `S` (e.g. ×2 for deuteron, ×3/2 for 7Li).

Either the `--fixedLS` command-line flag OR a bare `FIXEDLS` keyword in
the input deck selects the physics-standard convention instead:

```
sDotL = (1/2) * [J(J+1) - L(L+1) - S(S+1)] = <L*S>   ← with --fixedLS / FIXEDLS
```

CLI flag:

```bash
./ptolemy --fixedLS < my_input.in     # standard physics
./ptolemy           < my_input.in     # Cleopatra-compatible (default)
```

Or in the input deck (any top-level position):

```
HEADER: ...
REACTION: ...
FIXEDLS                                 # ← keyword equivalent of --fixedLS
PARAMETERSET dpsb r0target
...
```

**Default mode is bit-identical to Cleopatra/Maple/Ptolemy-f2c** — all
36 regression tests pass bit-identically.

**Equivalence between modes (verified):**

The two changes (kernel divisor `2S → 2`, and adjusting input `Vso`)
compensate exactly. For a given physical input:

| projectile | default mode `Vso` to use | `--fixedLS` mode `Vso` to use |
|---|---|---|
| p, n, t (spin-1/2) | `Vso_published` (no change) | `Vso_published / 2` |
| d (spin-1) | `Vso_published` (no change) | `Vso_published` (no change — modes coincide) |
| ³He (spin-1/2) | `Vso_published` (no change) | `Vso_published / 2` |
| ⁶Li (spin-1) | `Vso_published` (no change) | `Vso_published` (no change — modes coincide) |
| ⁷Li (spin-3/2) | `Vso_published × (3/2)` | `Vso_published × 3` |
| α (spin-0) | irrelevant (no spin-orbit kernel) | irrelevant |

Verified bit-identical: `./ptolemy --fixedLS` with `Vso=6.0` on
208Pb(p,p) Ep=30 produces the same DCS as `./ptolemy` (default) with
`Vso=3.0` on the same input — 0 curve-diff lines.

For spin-1 deuteron, the divisor in default mode is already `2S = 2`,
so `--fixedLS` is a no-op (also verified bit-identical on the
`elastic_40Ca_d.in` test).

**Use `--fixedLS` if your input `Vso` was tabulated for a code using
standard `<L*S>` coupling** (FRESCO, ECIS, textbook DWBA), and you want
to avoid manually adjusting `Vso` by the projectile spin. Default mode
is preferable when reproducing legacy Ptolemy / Cleopatra results
verbatim.

## License

Derived from the original Argonne Ptolemy code. No formal license file at
this time — please contact the maintainer for use beyond personal/research
purposes.
