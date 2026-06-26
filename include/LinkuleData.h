#pragma once
// LinkuleData.h — linkule (linked-potential) address block
// Extracted to its own header so Reaction.h can own a LinkuleData member
//
// A "linkule" is a linked group of potentials that share a common
// scaling parameter (e.g. all components of one optical potential).
// The 1-based [0..6] inner dimension overlays char8 names on the first
// two columns (linkuleAddr[k][1..2] = 8 chars).

constexpr int numLinkules = 13;

#include <vector>
#include <array>

// Linkule user-parameter block. The flat
// std::array<double, 21> PARAM_arr is accessed by the params(reaction)
// positional accessor and the 2x10-pair SPLINE_/LAGRANGE_linkule indexing:
//   PARAM_arr[0..4]  = PARAM1..PARAM5 (named input keywords)
//   PARAM_arr[5]     = dead padding (value-initialized 0.0, never written)
//   PARAM_arr[6..20] = PARAM6..PARAM20
struct LinkuleParams {
    std::array<double, 21> PARAM_arr;
};

struct LinkuleData {
    int uniqueLinkuleId;                         // serial counter for unique linkule ID tags ("*101", "*102", ...)
    int linkuleAddr[numLinkules + 1][7];              // 1-based: linkuleAddr(6,numLinkules)
    int lnkAd2[numLinkules + 1][3][7];           // per-channel linkule addresses

    // in JDEPENWS / JDEPEN / PARITWOO / SPLINE_linkule / LAGRANGE. Each fitter
    // stores a 1-based index into this vector in linkuleInts[0]; the requestCode==3,4
    // readers fetch the array via workarrs[linkuleInts[0] - 1].data() - 1.
    //
    // Append-only across the calculation. Memory bounded by the number of
    // fitter linkule activations (≤ numLinkules × max_components per problem).
    // The regression suite doesn't exercise these fitters, but the vector replaces the
    // NAMLOC-tagged pool slot cleanly: no test uses REDEF to pre-define a
    // worknm name (the names are derived from each fitter's iD), so the
    // NAMLOC probe was dead in this slot.
    std::vector<std::vector<double>> workarrs;
};

// Access via reaction.linkuleData (Reaction.h owns the storage).
