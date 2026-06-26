// defaults.cpp — DEFALT: seeds default control / parameter values.

#include "Timing.h"
#include <cstdio>
#include <cstring>
#include <ctime>
#include <unistd.h>
#include "Reaction.h"

// Machine/host banner info. Sole caller is the prologue below, so it lives here
// as a file-static.
static void getMachineInfo(char* mach, char* hostname) {
    std::strcpy(mach, "Intel Pentium, Redhat Linux");
    // Fortran hostnm intrinsic — use POSIX gethostname.
    char hBuf[256];
    if (gethostname(hBuf, sizeof(hBuf)) == 0) {
        std::strcpy(hostname, hBuf);
    } else {
        hostname[0] = '\0';
    }
}

// Positional Fortran-COMMON array accessors SWITCH_arr / FLOAT_arr /
// DEFALT bulk-init loops they served. Sentinel-init is now per-field
// (see undefValue / NOTDEF_INT explicit-write blocks in DEFALT below).

// ============================================================================
//     INITIALIZATION OF THE DWBA CODE
//
//     1/18/75 - SEPARATEDED FROM CONTRL - S.P.
//     1/18/75 - REARRANGED TO BE OVERLAYED - S. P.
//     2/10/75 - MOVE  NUCONLY AND JX
//     2/16/75 - LOAD IEWSZOVR, PRINT DATE
//     2/20/75 - MOST DEFAULTS SET TO 0 - S. P.
//     3/5/75 - NEW GRID DEFAULTS
//     4/15/75 - MAJOR REVISION FOR NEW COMMON BLOCKS
//     5/5/75 - STUFF FOR SETCHN
//     5/18/75 - RESTSW ADDED.
//     5/22/75 - vSi ADDED
//     5/30/75 - FIX FOR RESET
//     6/11/75 - SAVEHS/USEHS
//     8/7/75  - ADD R0E,VE,ETC TO COMMON
//     8/28/75 - LMINMULT, LMAXMULT, specFactorProj specFactorTgt
//     9/3/75 - MAXLEXTRAP ADDED; USECOULOMB IS DEFAULT
//     9/25/75 - CHANGE dwCutoff; REMOVE VCUT
//     9/30/75 - LMAXMULT BECOMES LMAXADD
//     10/1/75 - FITTING KEYWORDS
//     10/19/75 - ESQ PARAMETERS, TBEGIN
//     11/6/75 - DERIVSTEP, MAXFUNCTIONS;  LMCHOL IS DEFAULT FITTER
//     1/29/76 - FITMULTI, FITRATIO, FITMODE, NUMRANDO
//     2/26/76 - LMAXADD INCREASED FOR ELASTIC STUFF
//     6/1/76 - NVPOLY, sumDensity ADDED
//     6/17/76 - LARGE FAKE midpointFactor
//     9/16/76 - SEPT. 76 VERSION; FITMUL FIX
//     12/26/76 - INELASTIC SCATTERING
//     4/15/77 - lMinSub, ALMXAD ADDED
//     5/8/77 - PRINT=10001; MAY 77 VERSION
//     9/24/77 - USECORE DEFAULT; SEPT 77 VERSION
//     10/16/77 - POWERS OF WOODS-SAXONS
//     11/20/77 - SEPARATE BOUND & SCAT ASYMP; DEC 77 VERSION
//     11/25/77 - LINKULES
//     3/19/78 - MARCH, 1978 VERSION
//     4/24/78 - APRIL VERSION FOR li, lo REVERSALS - S.P.
//     4/28/78 - BIG AND SMALL STUFF IN /CNSTNT/ - S.P.
//     5/23/78 - JUNE 1978 VERSION - S.P.
//     12/7/78 - ADD TENSOR POTENTIALS - RPG
//     3/21/79 - APRIL 1979 VERSION - S.P.
//     5/6/79 - MAY 1979 VERSION; PARITIES AND LEVEL SPEC ARRAYS - S.P.
//     6/10/79 - CMS VERSION - S.P.
//     6/8/79 - JUNE 1979 VERSION; STANDARD PREFIXES, CHARACTER TYPE,
//        UNIVAC STUFF - RPG
//     6/18/79 - PROPERLY RESTORE IECHO FOR RESET - S.P.
//     9/10/79 - SEPTEMBER 1979 VERSION - RPG
//     11/3/79 - OCTOBER 1979 VERSION; TIME AND CPU VERSION STAMP - S.P.
//     11/15/79 - NOVEMBER 1979 VERSION; CMPUTR SUBROUTINE - S.P.
//     12/31/79 - DECEMBER 1979 VERSION; VAX STUFF - RPG
//     1/14/80 - FIX CVA IUNDEF/AUNDEF DATA; CHARACTER*1 BLANK1 - RPG
//     1/25/80 - CND CHARACTER*1 THEDAT, WRITE(6,-) TO PRINT - RPG
//     2/19/80 - INITIALIZE MAPPHI, phiMid, GAMPHI; FEB. 1980 VERSION - R
//     3/15/80 - MARCH 1980 VERSION FOR RPG MERGED BACK TO ANL - S.P.
//     5/2/80 - APRIL 1980 VERSION - C.C. WORKING FAIRLY WELL! S.P.
//     6/13/80 - CHAR*8 FOR CDC THEDAT, SHIFT FOR NOS - S.P.
//     7/10/80 - CWD, CHR PREFIXES, JULY 1980 VERSION - S.P.
//     7/17/80 - FIX CRAY notDefSentinel, undefValue VALUES - S.P.
//     8/25/80 - JULY 1980 VERSION, SOMEHOW LOST ON CDC - S.P.
//     9/16/80 - SEPT. 1980 VERSION; ZERO CCBLK; L-EXTRAP DEFAULT - S.P.
//     10/8/80 - OCTOBER 1980 VERSION FOR NEW BASCPL, ETC. - S.P.
//     12/17/80 - DECEMBER 1980 VERSION FOR COMPLETION OF ABOVE - S.P.
//     1/8/81 - JANUARY 1981 VERSION FOR NEW COUPLING SCHEME  - S.P.
//     6/24/81 - JUNE 1981 VERSION FOR REWORKS OF COUPLN - S.P.
//     8/28/81 - AUGUST 1981 VERSION FOR DEFORMED POTENS - S.P.
//     9/17/81 - SEPT 1981 VERSION; LEBACK DEFAULT - S.P.
//     10/12/81 - OCT 1981 VERSION FOR GIESSEN CHANGES - S.P.
//     2/11/82 - JANUARY 1982 VERSION
//     3/12/82 - MARCH 1982 VERSION
//     4/10/82 - APRIL 1982 VERSION
//     5/3/83 - MAY 1983 VERSION - VARIOUS SMALL FIXES - S.P.
//     7/20/83 - USE ISIZE=-500 FOR NOS - S.P
//     8/17/83 - AUGUST 1983 VERSION = FINAL AT MUNICH - S.P.
//     10/16/84 - OCT 1984 VERSION; scatAsy < 0 - S.P.
//     12/19/84 - DEC 1984 VERSION; ROOM FOR SIPOTENT - S.P.
//     5/10/85 - INITIALIZE PLM WORK ARRAY FOR CRA - S.P.
//     6/10/85 - JUNE, 1985 VERSION - S.P.
//     11/16/01 - Nov, 2001 version; new date routine
//
// ============================================================================

void Reaction::applyDefaults()
{
    //
    // replaced with 9 explicit NOTDEF_INT writes (L/lMax/lMin/nNodes/zArray[1..5]/
    // zProj/zTarget). The remaining 20 INTGER fields are either explicitly set to a
    // non-notDefSentinel value or zeroed by DEFALT below.
    // replaced with explicit notDefSentinel writes for all 9 J-block doubles
    // (J/js[1..5]/jProj/spinProj/spinTarget). Every J-block field has a notDefSentinel-guard reader,
    // so dropping any was not safe.

    int undefInteger = NOTDEF_INT;  // 0xF0F0F0F0 as signed int
    double undefDouble;
    {
        uint64_t tmp = 0xF0F0F0F0F0F0F0F0ull;
        std::memcpy(&undefDouble, &tmp, sizeof(double));
    }

    char dateString[10];    // CHARACTER*9 + NUL — fully written by getDate()
    char timeString[11];    // CHARACTER*10 + NUL — fully written by strftime() below

    char cpuId[41];
    std::memset(cpuId, ' ', sizeof(cpuId));
    cpuId[40] = '\0';

    char hostname[41];
    std::memset(hostname, ' ', sizeof(hostname));
    hostname[40] = '\0';

    //
    //
    //  INITIALIZE INPUT
    //
    //     FIRST SET EVERYTHING TO UNDEFINED STATUS
    {
        // notDefSentinel = undefInteger
        // In Fortran, notDefSentinel is an integer-sized bit pattern stored in a double-word.
        // We replicate: write the bit pattern of undefInteger into the first 4 bytes of notDefSentinel.
        double notDefScratch = 0.0;
        std::memcpy(&notDefScratch, &undefInteger, sizeof(int));
        internalState.notDefSentinel = notDefScratch;
    }
    internalState.undefValue = undefDouble;
    //
    // ReactionParams (FLOAT-block) fields that act as undefValue sentinels — read
    // somewhere via `if (X == undefValue)` / `if (X != undefValue)` guards. Fields not
    // listed here are either explicitly written below (V/vI/vSo/VSOI/vSi/exs/
    // accuracy/stepSize/angles/spec_*/gamma_*/midpointFactor/phiMid/
    // sumDensity/accuracyInel/alMnMt/alMxMt/boundAsy/scatAsy/dwCutoff) or have
    // no undefValue-guard reader at all. Replaces the former
    // `for i=1..NUMFLT: FLOAT(i) = undefValue` bulk loop + FLOAT_arr() positional
    {
        const double undefValue = internalState.undefValue;
        // Potential geometry — radii / diffuseness (real, imag, SO, surf, Coul)
        opticalPotentialParams.A     = undefValue;
        opticalPotentialParams.aI    = undefValue;
        opticalPotentialParams.aSo   = undefValue;
        opticalPotentialParams.aSoi  = undefValue;
        opticalPotentialParams.aSi   = undefValue;
        integrationGrid.R = undefValue;
        opticalPotentialParams.R0    = undefValue;
        opticalPotentialParams.rI    = undefValue;
        opticalPotentialParams.rI0   = undefValue;
        opticalPotentialParams.rSo   = undefValue;
        opticalPotentialParams.rSo0  = undefValue;
        opticalPotentialParams.rSoi  = undefValue;
        opticalPotentialParams.rSoi0 = undefValue;
        opticalPotentialParams.rSi   = undefValue;
        opticalPotentialParams.rSi0  = undefValue;
        opticalPotentialParams.rC    = undefValue;
        opticalPotentialParams.rC0   = undefValue;
        masses.rcProj   = undefValue;
        masses.rcTarget   = undefValue;
        masses.rc0Proj  = undefValue;
        masses.rc0Target  = undefValue;
        // Masses — massesArr[0..4] = mass_a, mass_b, mass_A, mass_B, mass_x
        // (aMs()[1..5] positional accessor over Reaction::Masses::massesArr).
        masses.aM        = undefValue;
        for (int i = 0; i < 5; i++) masses.massesArr[i] = undefValue;
        masses.massProj = undefValue;
        masses.massTgt  = undefValue;
        // Energies / Q-value
        energies.E    = undefValue;
        energies.eCm  = undefValue;
        energies.eLab = undefValue;
        energies.Q    = undefValue;
        // Spectroscopic / asymptote / step-spread
        spec.spam          = undefValue;
        integrationGrid.asymptopia     = undefValue;
        integrationGrid.stepsPerUnit = undefValue;
        // Sum-grid radii (computed in integration_grid / grid_setup if undefValue)
        integrationGrid.sumMax = undefValue;
        integrationGrid.sumMid = undefValue;
        integrationGrid.sumMin = undefValue;
        // Per-channel mass-excess arrays (1-based [1..5])
        for (int i = 1; i <= 5; i++) {
            masses.amxcs[i]  = undefValue;
            masses.amxcgs[i] = undefValue;
        }
        // Projectile/target vertex excitation + mass-number (1-based [1..2])
        for (int i = 1; i <= 2; i++) {
            energies.exsPt[i] = undefValue;
            masses.amxgPt[i]     = undefValue;
        }
        // User PARAM keywords — PARAM_arr[0..4] = PARAM1..PARAM5 (named input
        // keywords); PARAM_arr[6..20] = PARAM6..PARAM20 (former PAR620[1..15],
        // accessed positionally via params(reaction)[i]). PARAM_arr[5] is the
        // dead padding slot (former PAR620[0]) — left at 0.0 (value-init) to
        // preserve the bit-identical pre-extraction memory state.
        linkuleParams.PARAM_arr[0] = undefValue;
        linkuleParams.PARAM_arr[1] = undefValue;
        linkuleParams.PARAM_arr[2] = undefValue;
        linkuleParams.PARAM_arr[3] = undefValue;
        linkuleParams.PARAM_arr[4] = undefValue;
        for (int i = 6; i <= 20; i++) linkuleParams.PARAM_arr[i] = undefValue;
    }
    // INTGER-block fields used as NOTDEF_INT sentinels — explicit init replaces
    // the former `for i=1..NUMINT: INTGER(i) = notDefSentinel` bulk loop. Fields not
    // listed here are either explicitly written below to a non-notDefSentinel value
    // (lStep=1, nPhiPoints=10, nPhiSum=15, nPhiDifference=10, nPhiAdditional=4, printLevel=10001, nCoulombPoints=8,
    // lMaxAdditional=30, lMinSub=20, maxLExtrap=100) or are explicitly written to 0 below
    // (parities[1..5], parityPt[1..2], parity).
    angMom.L     = NOTDEF_INT;
    angMom.lMax  = NOTDEF_INT;
    angMom.lMin  = NOTDEF_INT;
    angMom.nNodes = NOTDEF_INT;
    for (int i = 1; i <= 5; i++) charges.zArray[i] = NOTDEF_INT;
    charges.zProj   = NOTDEF_INT;
    charges.zTarget   = NOTDEF_INT;
    // JBLOCK is a double-word array; the Fortran assignment stores the 4-byte
    // notDefSentinel bit pattern into its first 4 bytes (rest zero) — preserved here
    // via internalState.notDefSentinel (the partial-bit-pattern double).
    {
        const double notDefSentinel = internalState.notDefSentinel;
        angMom.J   = notDefSentinel;
        for (int i = 1; i <= 5; i++) angMom.js[i] = notDefSentinel;
        angMom.jProj  = notDefSentinel;
        angMom.spinProj = notDefSentinel;
        angMom.spinTarget = notDefSentinel;
    }
    //
    //     SET ALL OPTIONS OFF
    //
    // SWITCH-block fields are value-initialized to 0; the bulk-zero loop is redundant.
    // RESTSW always 0 via the single InputParser caller). Echo guard
    // in nxintReadCard was always-true; both went the same day.
    //
    for (int i = 1; i <= 45; i++) {
        reactStr[i] = ' ';
    }
    for (int i = 1; i <= 65; i++) {
        header[i] = ' ';
    }
    //
    //
    //     FIND OUT ABOUT THE MACHINE AND OPERATING SYSTEM
    //
    getMachineInfo(cpuId, hostname);
    //
    std::printf("%49s%s\n", "", "P T O L E M Y");
    //
    //     GET CURRENT DATE
    //
    getDate(dateString);
    //
    {
        time_t now = time(nullptr);
        struct tm* tmInfo = localtime(&now);
        strftime(timeString, sizeof(timeString), "%H%M%S.000", tmInfo);
    }
    //
    //     **********************************************************
    //     *                                                        *
    //     *                STOP     THINK                          *
    //     *                                                        *
    //     *    IS IT TIME TO CHANGE THE DATE OF THE VERSION ??     *
    //     *                                                        *
    //     **********************************************************
    //
    // Trim cpuId and hostname for printing (remove trailing spaces)
    {
        // trim cpuId
        int len = 40;
        while (len > 0 && cpuId[len - 1] == ' ') len--;
        cpuId[len] = '\0';
        // trim hostname
        len = 40;
        while (len > 0 && hostname[len - 1] == ' ') len--;
        hostname[len] = '\0';
    }

    std::printf(" April 2007 version     Computer: %s  %s          %s  %s\n\n",
                cpuId, hostname, dateString, timeString);
    //
    //
    //
    //     SETUP DEFAULT VALUES
    //
    integrationGrid.dwCutoff = 1.0E-3;
    integrationGrid.accuracy = 1.0E-12;
    integrationGrid.boundAsy = 20.0;
    integrationGrid.scatAsy = -20.0;
    // numbered-store pool size). Inlined at the three callers
    // (ISTART + DEFINE × 2). angular_setup was using the same field
    // as a local-in-disguise counter; it has its own `int ISIZE` now.
    integrationGrid.stepSize = 0.10;
    //  DEFAULT VALUES FOR ANGLES
    rxn.angleMin = 0;
    rxn.angleMax = 90;
    rxn.angleStep = 1;
    integrationGrid.lStep = 1;
    gridData.nPhiPoints = 10;
    gridData.nPhiSum = 15;
    gridData.nPhiDifference = 10;
    integrationGrid.sumDensity = 6;
    spec.specAmpProj = 1.0;
    spec.specAmpTgt = 1.0;
    spec.specFactorProj = 1;
    spec.specFactorTgt = 1;
    flags.printLevel = 10001;
    // internalState boundChannel/waveChannel/iDone/stripPickup zero-inits
    // that bulk-write notDefSentinel/undefValue garbage above, internalState is its own
    // struct outside those alias-array regions; the struct value-init covers
    // these single-shot defaults.
    rxn.gammaSum = 1;
    rxn.gammaDif = 5;
    integrationGrid.midpointFactor = 2;
    integrationGrid.phiMid = 0.5;
    gridData.nPhiAdditional = 4;
    flags.nuConL = 3;
    opticalPotentialParams.alMxMt = 1.6;
    angMom.lMaxAdditional = 30;
    opticalPotentialParams.alMnMt = 0.6;
    integrationGrid.lMinSub = 20;
    integrationGrid.maxLExtrap = 100;
    //
    //     EXCITATION ENERGIES ARE 0 BY DEFAULT
    //     parity = 0 MEANS UNDEFINED
    //
    for (int i = 1; i <= 5; i++) {
        angMom.parities[i] = 0;
        energies.exs[i] = 0;
    }
    angMom.parity = 0;
    //
    //
    // Use NOTDEF_INT constexpr instead of type-punning notDefSentinel double —
    // same 0xF0F0F0F0 bit pattern, no strict-aliasing UB.
    for (int i = 1; i <= 4; i++) {
        internalState.lSpecs[i] = NOTDEF_INT;
        internalState.nodesP[i] = NOTDEF_INT;
    }
    for (int i = 1; i <= 2; i++) {
        angMom.parityPt[i] = 0;
        internalState.eBnds[i] = internalState.undefValue;
    }
    internalState.lSpcPt2 = NOTDEF_INT;
    internalState.nodePt2 = NOTDEF_INT;
    //
    //     POTENTIALS ARE 0 BY DEFAULT
    //
    opticalPotentialParams.V = 0;
    opticalPotentialParams.vI = 0;
    opticalPotentialParams.vSo = 0;
    opticalPotentialParams.vSoi = 0;
    opticalPotentialParams.vSi = 0;
    //
    //  ENERGY DEPENDENCE ZERO BY DEFAULT
    //
    //   R0E / RI0E / R0ESQ / RI0ESQ — radius coeffs
    //   AE / AIE / AESQ / AIESQ — diffuseness coeffs
    //   VE / VIE / VESQ / VIESQ — depth coeffs
    // All 12 fields were DEFALT-only 0 with no input plumbing.
    //
    //     POTENTIALS ARE LINEAR WOODS SAXONS
    //
    //
    //     INDICATE NO LINKULES IN USE.
    //
    linkuleData.uniqueLinkuleId = 100;
    for (int i = 1; i <= numLinkules; i++) {
        linkuleData.linkuleAddr[i][3] = 0;
    }
    //
    //     FITTING KEYWORDS — only MAXFUN survives; the rest of the legacy
    //     fitter scalars (FITACC, FITMUL, FITRAT, MODEFT, NUMRAN, IREINI,
    //
    //
    //     INELASTIC SCATTERING
    //
    integrationGrid.accuracyInel = 1.0E-5;
    integrationGrid.nCoulombPoints = 8;
    //
    //     COUPLED CHANNELS
    //
    flags.excitationType = 1;
    //
    return;
}
