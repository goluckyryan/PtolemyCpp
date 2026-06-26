//
// GETLNK - processes LINKULE keywords and loads linkule
// LINKUL - linkule caller for Ptolemy (dispatches to specific linkule routines)
//
// 11/24/77 - first version
// 3/15/03 - add AV18
// 2/15/07 - add PHIFFER
//

#include "ptolemy_types.h"
#include "InputParser.h"
#include "linkule.h"
#include "Reaction.h"
#include "LinkulePlugin.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

// ============================================================================
//
// Processes linkule keywords and loads linkule.
// Determines what (real potential, imag potential, wave functions) is to be
// set by the linkule and also the linkule name.
// Then loads the linkule and saves its address.
//
// read back. The error/success paths now just return early; the printf
// already named the failure mode at the user.
// ============================================================================

void Reaction::loadLinkule(char8 linkuleKey)
{
    // implicit real*8 (a-h, o-z)

    auto* linkuleAddr = linkuleData.linkuleAddr;    // 1-based: linkuleAddr[i][j]

    // LNKNAM(1,i) <=> linkuleAddr(1,i) as character*8

    char8 word;

    // Key list for linkule usage types
    static const char8 keyNames[14] = {
        char8(),              // [0] unused
        char8("REALPOTE"),    // [1]
        char8("IMAGPOTE"),    // [2]
        char8("REALSOPO"),    // [3]
        char8("IMAGSOPO"),    // [4]
        char8("COULOMBP"),    // [5]
        char8("WAVEFUNC"),    // [6]
        char8("REALTRPO"),    // [7]
        char8("IMAGTRPO"),    // [8]
        char8("REALTLPO"),    // [9]
        char8("IMAGTLPO"),    // [10]
        char8("REALTPRP"),    // [11]
        char8("IMAGTPRP"),    // [12]
        char8("SIPOTENT")     // [13]
    };

    // List of built-in linkule names
    constexpr int numNames = 17;
    static const char8 names[18] = {
        char8(),              // [0] unused
        char8("BKGPTELP"),    // [1]
        char8("FIXEDWOO"),    // [2]
        char8("GAUSSIAN"),    // [3]
        char8("LAGRANGE"),    // [4]
        char8("LTSTELP "),    // [5]
        char8("RAWITSCH"),    // [6]
        char8("REID"),        // [7]
        char8("SHAPE"),       // [8]
        char8("SPLINE"),      // [9]
        char8("TWOSHAPE"),    // [10]
        char8("DEFORMED"),    // [11]
        char8("JDEPEN"),      // [12]
        char8("JDEPENWS"),    // [13]
        char8("OHTA"),        // [14]
        char8("PARITWOO"),    // [15]
        char8("AV18"),        // [16]
        char8("PHIFFER")      // [17]
    };

    static const char8 internal("INTERNAL");

    int keyIndex, linkuleIndex;

    // Find the key
    for (keyIndex = 1; keyIndex <= numLinkules; keyIndex++) {
        if (linkuleKey == keyNames[keyIndex]) break;
    }

    // Get the linkule name
    if (nxWord(word.data) != 0) {
        printf("\n**** A LINKULE NAME MUST FOLLOW THE %.8s KEYWORD.\n", linkuleKey.data);
        return;
    }

    // Store linkule name: first 8 bytes of linkuleAddr[keyIndex]
    std::memcpy(&linkuleAddr[keyIndex][1], &word, 8);

    linkuleIndex = 0;
    if (word != internal) {
        // In the non-IBM version, linkules are built-in
        for (int i = 1; i <= numNames; i++) {
            if (word == names[i]) { linkuleIndex = i; break; }
        }
        if (linkuleIndex == 0) {
            printf("\n*** %.8s IS NOT IN THIS VERSION OF PTOLEMY.\n", word.data);
            return;
        }
    }

    linkuleAddr[keyIndex][3] = linkuleIndex;
    printf(" LINKULE %.8s IS LINKULE NUMBER/location%8d\n", word.data, linkuleIndex);
}


// ============================================================================
// Maps a 1-based linkule index (the value loadLinkule stored in
// linkuleAddr[k][3]) to its LinkulePlugin subclass. Returns nullptr for the
// deleted/never-implemented slots (1/5/6/7/10/11/14) so the caller can emit
// "NOT AVAILABLE" exactly as the old default branch did.
// ============================================================================

std::unique_ptr<LinkulePlugin> makeLinkulePlugin(int linkuleIndex)
{
    switch (linkuleIndex) {
        case 2:  return makeFixedWoodsSaxonPlugin();
        case 3:  return makeGaussianPlugin();
        case 4:  return makeLagrangePlugin();
        case 8:  return makeShapePlugin();
        case 9:  return makeSplinePlugin();
        case 12: return makeJDependentWoodsSaxonPlugin();
        case 13: return makeJDependentWoodsSaxonFermiPlugin();
        case 15: return makeParityWoodsSaxonPlugin();
        case 16: return makeAV18Plugin();
        case 17: return makePhifferPlugin();
        default: return nullptr;
    }
}


// ============================================================================
//
// Linkule caller for Ptolemy.
// Adds the fixed part of a linkule argument list and calls the linkule.
// Non-IBM version: linkules are part of the processor and we just call them.
// ============================================================================

void linkule(int linkuleIndex, char8 alias, int* linkuleInts, int potType, int requestCode,
             int& callStatus, int L, double J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, char* iD,
             Reaction& reaction)
{
    // implicit real*8 (a-h, o-z)


    bool debugSwitch = (reaction.flags.printLevel % 10) == 9;

    if (debugSwitch) {
        printf(" CALLING LINKULE AT%8d %.8s%8d%8d%4d%4d %.4s\n",
               linkuleIndex, alias.data, linkuleInts[0], linkuleInts[1], potType, requestCode, iD);
    }

    // NOTE: F(112) is beginning of PARAM's
    //       i(46) is beginning of IPARAM's

    // All linkule routines take the same long argument list.
    // We dispatch based on linkuleIndex (1-17) to the appropriate built-in linkule.

    // Arguments common to all fitter-linkule calls (ALIAS + J dropped
    // out-param and phiffer reads jp, so the J parameter stays on those
    // two cleaner signatures. ALIAS kept on av18/phiffer for error printf.
    // (linkuleInts, potType, requestCode, callStatus,
    //  L, rStart, stepSize, nPts, array1, array2,
    //  F, T, F(112), i, JB, IS, IN, IN, C, WAV, KM,

    // For the dispatch we pass the raw int* pointer.

    // makeLinkulePlugin() (below) maps the index to the right LinkulePlugin
    // subclass; we call its uniform run(). A null result means an
    // unavailable/deleted slot. Each plugin body uses only the args it needs.
    std::unique_ptr<LinkulePlugin> plugin = makeLinkulePlugin(linkuleIndex);
    if (!plugin) {
        printf("\n**** LINKULE%10d NOT AVAILABLE.\n", linkuleIndex);
        std::exit(1122);  // was FSTOP(1122) from ptolemy_io.h
    }
    plugin->run(alias, linkuleInts, potType, requestCode, callStatus,
                L, J, rStart, stepSize, nPts, array1, array2, reaction);

    if (debugSwitch) {
        printf(" LINKULE RETURNS%10d%10d%10d\n",
               callStatus, linkuleInts[0], linkuleInts[1]);
    }
    return;

}
