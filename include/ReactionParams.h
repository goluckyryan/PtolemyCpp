#pragma once
// ReactionParams.h — reaction output/Gamow scalar parameters
// Extracted from Reaction.h. Member: reaction.rxn.

struct ReactionParams {
    double angleStep;   // angular step size for cross section output (degrees)
    double angleMin;   // minimum output angle (degrees; cm or lab per outputInLab)
    double angleMax;   // maximum output angle (degrees)
    double gammaDif;    // gamma_diff: difference of absorption widths (in Gamow factor calculation)
    double gammaSum;   // gammaSum: sum of absorption widths
};
