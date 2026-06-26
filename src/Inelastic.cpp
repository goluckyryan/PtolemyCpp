// inelastic_reaction.cpp — Inelastic scattering class implementation.
// Composes the per-phase class methods (DWBAGrid::inelasticGridSet/inelasticRadialIntegrals,
// CrossSectionCalc::linterp/xSection, CoulombWaveFunction::scattering).

#include "Inelastic.h"
#include "Elastic.h"
#include "Reaction.h"
#include "CoulombWaveFunction.h"
#include "CrossSectionCalc.h"

InelasticReaction::InelasticReaction(Reaction& reaction) : reaction_(reaction) {}

bool InelasticReaction::calculate()
{
    int returnCode = 1;

    // Elastic scattering must run first — it produces the S-matrix
    // and distorted waves used by the inelastic calculation.
    // Save/restore rxnParams because computeElasticScattering modifies them.
    auto savedRxn = reaction_.rxn;
    Elastic(reaction_).calculate();
    reaction_.rxn = savedRxn;

    // Setup grid: PRBPRT + COULST + INGRST + ANGSET + scattering
    reaction_.probePrint(returnCode);
    if (returnCode == 0) return false;
    CoulombWaveFunction::computeScatteringWaves(returnCode, reaction_);
    if (returnCode == 0) return false;
    reaction_.dwbaGrid.inelasticGridSet(returnCode, reaction_);
    if (returnCode == 0) return false;
    CrossSectionCalc(reaction_).angleSet();
    CoulombWaveFunction::scattering(returnCode, reaction_);
    if (returnCode == 0) return false;

    // Inelastic radial integrals + L-interpolation + cross section.
    reaction_.dwbaGrid.inelasticRadialIntegrals(reaction_);
    CrossSectionCalc(reaction_).linterp();
    CrossSectionCalc(reaction_).xSection();
    return true;
}
