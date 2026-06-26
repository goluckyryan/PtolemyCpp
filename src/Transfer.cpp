// transfer.cpp — Transfer reaction class implementation.
// Composes the per-phase class methods (DWBAGrid::gridSet/inelDc,
// CrossSectionCalc::linterp/xSection).

#include "Transfer.h"
#include "Reaction.h"
#include "CoulombWaveFunction.h"
#include "CrossSectionCalc.h"

Transfer::Transfer(Reaction& reaction) : reaction_(reaction) {}

bool Transfer::calculate()
{
    int returnCode = 1;

    // Setup grid: PRBPRT + COULST + GRDSET + ANGSET.
    reaction_.probePrint(returnCode);
    if (returnCode == 0) return false;
    CoulombWaveFunction::computeScatteringWaves(returnCode, reaction_);
    if (returnCode == 0) return false;
    reaction_.dwbaGrid.gridSet(returnCode, reaction_);
    if (returnCode == 0) return false;
    CrossSectionCalc(reaction_).angleSet();

    // Transfer DWBA radial integrals + L-interpolation + cross section.
    reaction_.dwbaGrid.inelDc(reaction_);
    CrossSectionCalc(reaction_).linterp();
    CrossSectionCalc(reaction_).xSection();
    return true;
}
