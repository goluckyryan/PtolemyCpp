#pragma once
// DwbaInputExpander.h — expand a human-readable DWBA reaction description into
// a full Ptolemy input deck.
//
// This is the InFileCreator workflow from digios_master/analysis/Cleopatra,
// absorbed into PtolemyCpp so users can feed a one-line-per-reaction file like
//
//   206Hg(d,p)207Hg   0+   1g9/2   9/2+   0.000   7.39MeV/u   AK
//
// directly to ptolemy. Each line has the form
//
//   target(in,out)residual  gsJpi  orbital  Jpi(Ex)  Ex  ELab  Potentials [beta]
//
//   gsJpi      ground-state spin-parity of the target, e.g. 0+
//   orbital    transferred orbital: "none" (elastic/inelastic),
//              "1g9/2" (single-nucleon, node-l-j), or "0L=2" (two-nucleon, L=)
//   Jpi(Ex)    spin-parity of the populated state
//   Ex         excitation energy [MeV]
//   ELab       beam energy: "30MeV" (total) or "7.39MeV/u" (per nucleon)
//   Potentials 1 char for elastic/inelastic, 2 chars (in,out) for transfer.
//              See OpticalPotentialLibrary potentialRef() for the codes.
//   beta       optional deformation length for inelastic scattering.
//
// Lines beginning with '#' and lines shorter than 5 characters are ignored.

#include <string>

struct DwbaExpander {
    // Expand a multi-line DWBA description into a Ptolemy input deck.
    // dwbaInput is the full file content; the returned string is the deck.
    // Unrecognized / inconsistent reaction lines are skipped with a diagnostic
    // on stderr (mirroring the original InFileCreator behavior).
    static std::string expand(const std::string& dwbaInput,
                              double angMin = 0.0, double angMax = 180.0,
                              double angStep = 1.0);

    // Heuristic: does the first non-blank, non-comment line look like a DWBA
    // reaction description (target(in,out)residual followed by >=6 tokens)?
    // Used by InputParser to auto-detect DWBA vs. native Ptolemy decks.
    static bool looksLikeDwba(const std::string& input);
};
