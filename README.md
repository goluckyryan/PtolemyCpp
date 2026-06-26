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

```
make
./test_ptolemy.sh
```

## Stats (current)

- 28,668 LOC across 54 `.cpp` / 49 `.h`
- 35/35 reaction tests bit-identical to Cleopatra (via Maple)
- 58/58 unit tests pass
- PtolemyCpp ~0.89s wall vs Maple ~4.2s (~4.7× faster) on the full 35-test suite at `-O2`

## Architecture

See `flow.md` for the class diagram. (A developer who clones with full history
can browse it; it documents how `Reaction`, `OpticalPotential`, `BoundState`,
the scattering solver, and the cross-section pipeline fit together.)

## Status

Stable. Bit-identical regression-tested. Open for contribution.

## Contributing

Run `./test_ptolemy.sh` to verify any change is bit-identical against the Maple
oracle. PRs welcome.

## License

Derived from the original Argonne Ptolemy code. No formal license file at this
time — derived from internally-distributed source.
