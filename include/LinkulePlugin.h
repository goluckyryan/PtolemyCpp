#pragma once
// LinkulePlugin.h — OOP hierarchy for the "linkule" (linked-potential) plugins,
//
// WHY A SEPARATE BASE FROM Potential (deviation from the Phase-4 spec wording):
// the spec asked for each linkule to become a `Potential` subclass overriding
// fill(vRay,nPts,rStart,stepSize). That fits the array-fill optical fitters, but
// NOT the two linkules any test actually exercises — av18 / phiffer are bound-
// state wavefunction SOLVERS dispatched from BoundState.cpp across the full
// linkule lifecycle: requestCode 1 (setup: writes jp + reaction fields),
// 2 (print), and 3 (compute), and they need the `jp` out-param and the second
// output buffer (array2). A single fill() cannot carry that. So the linkule
// plugins get their own base whose virtual run() matches the real linkule()
// calling convention; each subclass owns one plugin's body verbatim, and a
// factory replaces the switch. The optical fitters share the same base (they
// are dispatched identically); Reaction::makePotential adapts them to a flat
// OpticalPotential grid by marshalling the fill args into run() (requestCode 3)
// inline in its fillLinkule helper.
//
// PERSISTENT STATE: phiffer keeps state across calls (rc 1 sets it, rc 3 reads
// it). Plugin objects are constructed per linkule() call by the factory, so that
// state stays in `static` locals inside the run() body (shared, exactly as the
// former file-statics were) — never in plugin members.

#include "ptolemy_types.h"   // char8
#include <memory>

class Reaction;

// One linkule plugin. run() reproduces the former linkule() dispatch for one
// slot; the arg list is the union the former trimmed free-fn signatures needed
// (each body ignores the params it does not use, as in the original Fortran).
struct LinkulePlugin {
    virtual ~LinkulePlugin() = default;
    virtual void run(char8 alias, int* linkuleInts, int potType, int requestCode,
                     int& callStatus, int L, double& J, double rStart,
                     double stepSize, int nPts, double* array1, double* array2,
                     Reaction& reaction) = 0;
};

// Factory: build the plugin for a 1-based linkule index (the value stored in
// linkuleAddr[k][3] by loadLinkule). Returns nullptr for an unavailable slot
// (caller emits "NOT AVAILABLE" + FSTOP, matching the old default branch).
std::unique_ptr<LinkulePlugin> makeLinkulePlugin(int linkuleIndex);

// Per-plugin factories, each defined in its own TU where the subclass is
// complete; makeLinkulePlugin() dispatches to these by index. (Added one per
std::unique_ptr<LinkulePlugin> makeFixedWoodsSaxonPlugin();          // linkule 2
std::unique_ptr<LinkulePlugin> makeGaussianPlugin();                 // linkule 3
std::unique_ptr<LinkulePlugin> makeLagrangePlugin();                 // linkule 4
std::unique_ptr<LinkulePlugin> makeShapePlugin();                    // linkule 8
std::unique_ptr<LinkulePlugin> makeSplinePlugin();                   // linkule 9
std::unique_ptr<LinkulePlugin> makeJDependentWoodsSaxonPlugin();     // linkule 12
std::unique_ptr<LinkulePlugin> makeJDependentWoodsSaxonFermiPlugin();// linkule 13
std::unique_ptr<LinkulePlugin> makeParityWoodsSaxonPlugin();         // linkule 15
std::unique_ptr<LinkulePlugin> makeAV18Plugin();                     // linkule 16
std::unique_ptr<LinkulePlugin> makePhifferPlugin();                  // linkule 17
