#pragma once
// OpticalPotential.h — the COMPUTED distorted-wave optical potential V(r) on a
// channel grid. The flat USER-INPUT parameters (V, A, R0, vSo, ...) are the
// separate OpticalPotentialParams struct (OpticalPotentialParams.h).
//
// This is the flat-class redesign that replaces the virtual `Potential`
// hierarchy (include/Potential.h): instead of a base + WoodsSaxon/SpinOrbit/
// Surface/Coulomb subclasses dispatched through a composite, one fat class owns
// one `fill*` method per V(r) form. Adding a form = adding a method, no factory,
// no virtual dispatch.
//
// GRID / INDEXING. `values` is 0-based with size nPts; legacy Numerov code wants
// a 1-based array, so data1Based() hands back values.data()-1 (valid subscripts
// 1..nPts). The fill* methods write the 1..nPts range via that 1-based view, so
// their bodies stay byte-for-byte identical to the former woodsX()/Potential::fill
// loops (region-bracket + multiplicative exp recurrence — see Potential.h on why
// a per-point std::exp() is NOT bit-identical).
//
// NAMING of the grid spec: to keep the fill* bodies verbatim from the old fill()
// (whose signature was fill(vRay, nPts, rStart, stepSize) and which the composite
// called as fill(out, nPts, 0.0, channelRStart)), this class stores the same two
// quantities under the same names: rStart == the fill() `rStart` arg (0.0 at the
// grid origin) and stepSize == the fill() `stepSize` arg (the channel rStart).

#include <vector>

class OpticalPotential {
public:
    std::vector<double> values;   // grid values, 0-based, size = nPts
    int    nPts     = 0;
    double rStart   = 0.0;        // fill()'s rStart arg (grid origin, 0.0)
    double stepSize = 0.0;        // fill()'s stepSize arg (channel rStart)

    // Size the grid + record its spec; zeroes values.
    void resize(int nPts, double rStart, double stepSize);
    // Zero values (keeps size + grid spec).
    void clear();

    // COMPOSITION. add() accumulates another grid of identical size; scale()
    // multiplies every value by k (e.g. the Numerov -h2/12E factor).
    void add(const OpticalPotential& other);
    void scale(double k);

    // 1-based view for legacy Numerov code and the fill* loop bodies.
    double* data1Based() { return values.data() - 1; }

    // BUILDERS (one per V(r) form). Each ADDS its contribution to values[1..nPts]
    // (via data1Based()), so they compose onto a pre-filled buffer (e.g. Coulomb
    // then Woods-Saxon for the real-central channel). Bodies are the former
    // Potential subclass fill() bodies, with region-1 bounded at n1-1 instead of
    // n1 to drop the (formerly harmless overwrite) overlap at n1 that ADD would
    // otherwise double-count — region-2 always owns n1 (n2 = max(n2,n1)), so the
    // final values are byte-for-byte identical.
    void fillWoodsSaxon(double V, double R, double a);   // central WS (potForm 1)
    void fillSpinOrbit (double V, double R, double a);   // derivative-WS/r (potForm 2)
    void fillSurface   (double V, double R, double a);   // derivative-WS (potForm 3)
    void fillCoulomb   (int channelIndex, int channelIndexIn);  // vcsq12 sphere/point
};
