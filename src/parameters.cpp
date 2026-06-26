// parameters.cpp — PARAM, REACTN, GSINFO: parameter / reaction / ground-state setup.

#include "ptolemy_types.h"
#include "masstable.h"
#include <cstdio>
#include <cstring>
#include "Reaction.h"
#include "InputParser.h"

// File-scope alias to the InputParser-owned scanner buffer. Was an
namespace { auto& inputBuffer = InputParser::buffer(); }

// aMs(r) — 0-based pointer into the massesArr std::array<double, 5>
// Last of three TUs with positional aMs(reaction)[i] callers (others:
// set_channels.cpp, probe_print.cpp) — all three anon-namespace copies are now
// 0-based; callers shift their 1-based index by -1.
namespace {
double* aMs(Reaction& r) { return r.masses.massesArr.data(); }  // 0-based

// Extract particle i's integer mass number (A), charge (Z), and neutron number
// (N = A - Z) from the reaction arrays. The same A/Z/N read precedes both the
// mass-excess scan and the ground-state-spin scan in parseReactionString().
void readParticleAZN(Reaction& reaction, int i, int& ia, int& iz, int& n) {
    ia = static_cast<int>(aMs(reaction)[i - 1]);
    iz = reaction.charges.zArray[i];
    n  = ia - iz;
}

// Normalize the "explicitly-zero" excitation sentinel (1e-30) back to 0 for
// particle i. The 1e-30 marker means E* was entered as exactly zero; this is
// repeated verbatim across several exit branches of the per-particle J loop in
// parseReactionString() (and assignGroundStateSpin), so it lives here once.
void clearExcitationSentinel(Reaction& reaction, int i) {
    if (reaction.energies.exs[i] == 1.e-30) reaction.energies.exs[i] = 0;
}

// Finalize the ground-state spin/parity for particle i and normalize the
// "explicitly-zero" excitation sentinel (1e-30) back to 0. Shared verbatim by
// the dineutron/diproton (J=0) and known-G.S.-spin exit branches of the
// per-particle J loop in parseReactionString().
void assignGroundStateSpin(Reaction& reaction, int i, int jVal, int gsParity) {
    reaction.angMom.js[i] = jVal;
    if (gsParity != 0) reaction.angMom.parities[i] = gsParity;
    clearExcitationSentinel(reaction, i);
}
}

// ============================================================================
//
// Returns ground state properties of nuclei.
// Data from the mass-card by Ajenberg-Selove.
//
// kZ      - charge Z of the nucleus
// kA      - atomic number A of the nucleus
// jVal    - set to 2*J of ground state (output)
// parity  - set to +1, -1, or 0 if not known (output)
// notTabulated  - 0 = found, +1 = not tabulated (output)
//
// Fortran's CHARACTER*4 JPI out-param dropped — both callers passed an
// ALPHA[5] scratch buffer they never read. The 4-byte string is now an
// internal local since the J / parity decoder still needs it.
// ============================================================================

void groundStateInfo(int kZ, int kA, int& jVal, int& parity, int& notTabulated)
{
    // INTEGER*2 IE, IA arrays (270 elements each)
    static const integer2 kZLoArr[271] = { 0, // index 0 unused (1-based)
    //  IE1 (1..90)
           0,   0,   1,   1,   1,   1,   2,   3,   2,   3,
           4,   4,   5,   5,   6,   6,   6,   7,   8,   8,
           8,   9,   9,  10,  10,  11,  11,  12,  12,  13,
          13,  14,  14,  15,  15,  16,  16,  16,  16,  17,
          17,  18,  18,  19,  19,  19,  19,  19,  20,  20,
          20,  21,  22,  23,  23,  24,  24,  25,  25,  26,
          26,  26,  27,  27,  28,  28,  28,  28,  29,  30,
          30,  30,  30,  31,  31,  32,  32,  32,  32,  33,
          33,  33,  33,  34,  34,  35,  35,  35,  35,  35,
    //  IE2 (91..180)
          36,  36,  36,  36,  36,  37,  38,  39,  40,  40,
          41,  41,  42,  42,  42,  42,  43,  44,  44,  45,
          45,  46,  46,  46,  46,  47,  47,  48,  48,  48,
          49,  49,  49,  49,  49,  50,  50,  50,  50,  50,
          51,  51,  51,  51,  52,  53,  53,  53,  54,  54,
          54,  54,  54,  54,  57,  58,  58,  59,  59,  59,
          60,  60,  61,  61,  61,  62,  62,  62,  63,  63,
          63,  64,  64,  65,  65,  66,  66,  67,  67,  67,
          67,  68,  68,  68,  69,  69,  69,  70,  70,  71,
    //  IE3 (181..270)
          71,  72,  72,  72,  73,  73,  73,  74,  74,  74,
          75,  75,  75,  76,  76,  76,  77,  77,  77,  78,
          78,  78,  79,  79,  80,  80,  80,  81,  81,  81,
          81,  82,  82,  82,  82,  83,  83,  84,  84,  85,
          85,  86,  86,  86,  87,  88,  88,  88,  88,  89,
          89,  89,  90,  90,  90,  91,  91,  91,  92,  92,
          92,  93,  93,  94,  94,  94,  94,  95,  96,  96,
          96,  97,  98,  98,  98,  99,  99, 100, 101, 104,
         105, 105, 106, 106, 107, 107, 108, 108, 109, 200
    };

    static const integer2 ajBaseArr[271] = { 0, // index 0 unused (1-based)
    //  IA1 (1..90)
           0,   0,   2,   3,   5,   7,  10,  13,  15,  19,
          23,  26,  29,  32,  35,  38,  41,  45,  49,  52,
          55,  59,  63,  67,  71,  75,  79,  83,  86,  90,
          93,  97, 100, 104, 107, 111, 114, 118, 123, 128,
         132, 137, 141, 146, 150, 154, 158, 163, 168, 173,
         178, 184, 189, 194, 198, 203, 207, 212, 216, 221,
         225, 229, 234, 238, 243, 247, 252, 257, 262, 266,
         270, 275, 280, 285, 290, 296, 301, 306, 311, 317,
         322, 327, 332, 338, 343, 349, 354, 360, 366, 373,
    //  IA2 (91..180)
         380, 387, 394, 402, 410, 419, 427, 434, 441, 447,
         455, 462, 469, 476, 483, 490, 497, 504, 510, 517,
         523, 529, 534, 540, 546, 553, 559, 565, 571, 578,
         585, 591, 597, 603, 609, 615, 622, 629, 636, 644,
         652, 660, 668, 676, 684, 692, 700, 708, 716, 722,
         729, 736, 744, 754, 764, 771, 778, 785, 792, 799,
         806, 813, 821, 828, 836, 844, 851, 857, 864, 870,
         876, 882, 888, 894, 899, 905, 910, 915, 920, 926,
         932, 938, 943, 949, 955, 961, 967, 977, 983, 989,
    //  IA3 (181..270)
         994,1000,1005,1011,1017,1025,1033,1041,1048,1055,
        1062,1068,1074,1081,1087,1093,1100,1106,1112,1119,
        1125,1132,1139,1146,1153,1160,1167,1175,1182,1189,
        1196,1203,1210,1217,1225,1233,1239,1245,1250,1255,
        1259,1263,1267,1271,1276,1280,1284,1288,1292,1297,
        1301,1306,1311,1316,1321,1326,1330,1335,1340,1345,
        1350,1355,1360,1365,1370,1375,1381,1387,1392,1397,
        1402,1407,1412,1416,1421,1426,1430,1435,1440,1444,
        1445,1446,1447,1448,1449,1450,1451,1452,1453,1454
    };

    // AJ — 1485 entries of CHARACTER*4 ground state J-pi strings
    // AJ1..AJ33, 45 entries each = 1485 total
    static const char ajData[1486][5] = { "", // index 0 unused (1-based)
    // AJ1 (1..45)
    /*   1 */ "1/2+", "1/2+", " 1+ ", "1/2+", "1/2+",
    /*   6 */ "X   ", " 0+ ", "X   ", "3/2-", "3/2-",
    /*  11 */ " 0+ ", " 1+ ", " 0+ ", "3/2-", "3/2-",
    /*  16 */ " 0+ ", " 2+ ", " 0+ ", " 2+ ", "3/2-",
    /*  21 */ "3/2-", "3/2-", "3/2-", " 0+ ", " 3+ ",
    /*  26 */ " 0+ ", "1/2+", "3/2-", "3/2-", " 1+ ",
    /*  31 */ " 0+ ", " 1+ ", "3/2-", "1/2-", "1/2-",
    /*  36 */ " 0+ ", " 1+ ", " 0+ ", "1/2+", "1/2-",
    /*  41 */ "1/2-", " 0+ ", " 2- ", " 0+ ", " 0- ",
    // AJ2 (46..90)
    /*  46 */ "1/2-", "5/2+", "5/2+", "1/2-", " 0+ ",
    /*  51 */ " 1+ ", " 0+ ", "5/2+", "1/2+", "1/2+",
    /*  56 */ " 0+ ", " 2+ ", " 0+ ", " 2+ ", "5/2+",
    /*  61 */ "3/2+", "3/2+", "5/2+", "    ", " 0+ ",
    /*  66 */ " 3+ ", " 0+ ", "5/2+", "3/2+", "3/2+",
    /*  71 */ "5/2+", " 0+ ", " 4+ ", " 0+ ", " 4+ ",
    /*  76 */ "5/2+", "5/2+", "5/2+", "5/2+", ",3+ ",
    /*  81 */ " 0+ ", " 5+ ", " 0+ ", "1/2+", "5/2+",
    /*  86 */ "5/2+", " 0+ ", " 3+ ", " 0+ ", " 3+ ",
    // AJ3 (91..135)
    /*  91 */ "5/2+", "1/2+", "1/2+", "2,3+", " 0+ ",
    /*  96 */ " 1+ ", " 0+ ", "3/2+", "1/2+", "1/2+",
    /* 101 */ " 0+ ", " 1+ ", " 0+ ", " 2+ ", "1/2+",
    /* 106 */ "3/2+", "3/2+", " 1+ ", " 0+ ", " 0+ ",
    /* 111 */ " 0+ ", "3/2+", "3/2+", "3/2+", " 0+ ",
    /* 116 */ " 2+ ", " 0+ ", "X   ", "7/2-", "3/2+",
    /* 121 */ "3/2+", "3/2+", "3/2+", " 0+ ", " 2- ",
    /* 126 */ " 0+ ", " 3+ ", " 0+ ", "3/2+", "7/2-",
    /* 131 */ "3/2+", "3/2+", " 2- ", " 0+ ", " 4- ",
    // AJ4 (136..180)
    /* 136 */ " 0+ ", " 4- ", "7/2-", "3/2+", "7/2-",
    /* 141 */ "7/2-", " 0+ ", " 2- ", " 0+ ", " 0+ ",
    /* 146 */ " 0+ ", "3/2+", "7/2-", "7/2-", "7/2-",
    /* 151 */ " 2- ", " 0+ ", " 2+ ", " 0+ ", "3/2+",
    /* 156 */ "7/2-", "7/2-", "7/2-", " 2- ", " 0+ ",
    /* 161 */ " 4+ ", " 0+ ", " 0+ ", "1/2+", "7/2-",
    /* 166 */ "7/2-", "5/2-", "3/2-", " 0+ ", " 6+ ",
    /* 171 */ " 0+ ", " 4+ ", " 0+ ", "3/2-", "7/2-",
    /* 176 */ "7/2-", "7/2-", "5/2-", " 0+ ", " 5+ ",
    // AJ5 (181..225)
    /* 181 */ " 0+ ", " 6+ ", " 0+ ", " 0+ ", "7/2-",
    /* 186 */ "3/2-", "7/2-", "7/2-", "5/2-", " 0+ ",
    /* 191 */ " 3+ ", " 0+ ", " 6+ ", " 0+ ", "7/2-",
    /* 196 */ "3/2-", "7/2-", "7/2-", " 5+ ", " 0+ ",
    /* 201 */ " 3+ ", " 0+ ", " 0+ ", "3/2-", "5/2-",
    /* 206 */ "3/2-", "7/2-", " 0+ ", " 3+ ", " 0+ ",
    /* 211 */ " 4+ ", " 0+ ", "5/2-", "1/2-", "7/2-",
    /* 216 */ "3/2-", "    ", " 0+ ", " 2+ ", " 0+ ",
    /* 221 */ " 1+ ", "3/2-", "7/2-", "3/2-", "3/2-",
    // AJ6 (226..270)
    /* 226 */ " 0+ ", " 5+ ", " 0+ ", " 2+ ", "3/2-",
    /* 231 */ "7/2-", "3/2-", "3/2-", "3/2-", " 5+ ",
    /* 236 */ " 0+ ", " 1+ ", " 0+ ", "    ", "1/2-",
    /* 241 */ "3/2-", "3/2-", "3/2-", " 0+ ", " 1+ ",
    /* 246 */ " 0+ ", " 0+ ", "5/2-", "3/2-", "5/2-",
    /* 251 */ "3/2-", "    ", " 0+ ", " 1+ ", " 0+ ",
    /* 256 */ " 0+ ", " 0+ ", "    ", "3/2-", "5/2-",
    /* 261 */ "3/2-", "    ", " 1+ ", " 0+ ", " 1+ ",
    /* 266 */ " 0+ ", "1/2-", "3/2-", "5/2-", "    ",
    // AJ7 (271..315)
    /* 271 */ " 0+ ", " 1+ ", " 0+ ", " 4  ", " 0+ ",
    /* 276 */ "1/2-", "3/2-", "1/2-", "5/2-", "5/2-",
    /* 281 */ " 0+ ", " 3- ", " 0+ ", " 2- ", " 0+ ",
    /* 286 */ "3/2-", "9/2+", "3/2-", "9/2+", "X   ",
    /* 291 */ " 3- ", " 0+ ", " 2- ", " 0+ ", " 1+ ",
    /* 296 */ " 0+ ", "1/2-", "3/2-", "5/2+", "3/2-",
    /* 301 */ "    ", " 0+ ", " 2- ", " 0+ ", " 1  ",
    /* 306 */ " 0+ ", "7/2+", "3/2-", "1/2-", "3/2-",
    /* 311 */ "7/2+", " 0+ ", " 2- ", " 0+ ", " 1+ ",
    // AJ8 (316..360)
    /* 316 */ " 0+ ", "    ", "3/2-", "7/2+", "3/2-",
    /* 321 */ "1/2-", "3/2-", " 1+ ", " 0+ ", " 1+ ",
    /* 326 */ " 0+ ", " 1+ ", "3/2-", "1/2-", "3/2-",
    /* 331 */ "7/2+", "3/2-", " 5- ", " 0+ ", " 5- ",
    /* 336 */ " 0+ ", " 1+ ", " 0+ ", "9/2+", "3/2-",
    /* 341 */ "9/2+", "5/2-", "7/2+", " 0+ ", " 2- ",
    /* 346 */ " 0+ ", " 2- ", " 0+ ", " 4- ", "3/2-",
    /* 351 */ "9/2+", "5/2-", "9/2+", "9/2+", "1,2 ",
    /* 356 */ " 0+ ", " 2- ", " 0+ ", " 4- ", " 0+ ",
    // AJ9 (361..405)
    /* 361 */ "    ", "5/2+", "3/2-", "9/2+", "1/2-",
    /* 366 */ "    ", "    ", " 0+ ", " 2- ", " 0+ ",
    /* 371 */ " 4- ", " 0+ ", " 8+ ", "    ", "5/2+",
    /* 376 */ "3/2-", "5/2+", "1/2-", "9/2+", "9/2+",
    /* 381 */ " 0+ ", " 1- ", " 0+ ", " 2- ", " 0+ ",
    /* 386 */ " 8+ ", " 0+ ", "5/2+", " -  ", "5/2+",
    /* 391 */ "1/2-", "5/2+", "9/2+", "9/2+", " 0+ ",
    /* 396 */ "    ", " 0+ ", " 2- ", " 0+ ", " 7+ ",
    /* 401 */ " 0+ ", " 8+ ", "    ", "    ", "    ",
    // AJ10 (406..450)
    /* 406 */ "1/2-", "5/2+", "9/2+", "5/2+", "9/2+",
    /* 411 */ " 0+ ", "    ", " 0+ ", " 2- ", " 0+ ",
    /* 416 */ " 6+ ", " 0+ ", "6,7+", " 0+ ", "    ",
    /* 421 */ "    ", " -  ", "5/2+", "9/2+", "5/2+",
    /* 426 */ "9/2+", "5/2+", " 0+ ", "    ", " 0+ ",
    /* 431 */ " 5+ ", " 0+ ", " 7+ ", " 0+ ", "1/2-",
    /* 436 */ "1/2+", "9/2+", "5/2+", "9/2+", "5/2+",
    /* 441 */ "9/2+", " 0+ ", " 1+ ", " 0+ ", "7,6+",
    /* 446 */ " 0+ ", " 2+ ", "    ", "9/2+", "1/2+",
    // AJ11 (451..495)
    /* 451 */ "9/2+", "5/2+", "1/2-", "5/2+", "    ",
    /* 456 */ "    ", " 0+ ", " 1+ ", " 0+ ", " 1- ",
    /* 461 */ " 0+ ", " 5+ ", "    ", "1/2+", "9/2+",
    /* 466 */ "5/2+", "1/2-", "5/2+", "9/2+", " 0+ ",
    /* 471 */ " 1+ ", " 0+ ", " 1- ", " 0+ ", " 5+ ",
    /* 476 */ " 0+ ", "    ", "    ", "5/2+", "1/2-",
    /* 481 */ "5/2+", "7/2+", "    ", " 0+ ", "    ",
    /* 486 */ " 0+ ", " 1+ ", " 0+ ", " 5+ ", " 0+ ",
    /* 491 */ "    ", "    ", "5/2+", "7/2+", "5/2+",
    // AJ12 (496..540)
    /* 496 */ "1/2-", "5/2+", "    ", " 0+ ", " 1+ ",
    /* 501 */ " 0+ ", " 1+ ", " 0+ ", "    ", "    ",
    /* 506 */ "5/2+", "5/2+", "1/2-", "5/2+", "9/2+",
    /* 511 */ " 0+ ", " 1+ ", " 0+ ", " 1+ ", " 0+ ",
    /* 516 */ " 2+ ", " 0+ ", "    ", "5/2+", "1/2-",
    /* 521 */ "5/2+", "9/2+", "    ", " 1+ ", " 0+ ",
    /* 526 */ " 1+ ", " 0+ ", " 2+ ", " 0+ ", "5/2+",
    /* 531 */ "1/2-", "1/2+", "9/2+", "7/2+", " 0+ ",
    /* 536 */ " 2- ", " 0+ ", " 1+ ", " 0+ ", " 3+ ",
    // AJ13 (541..585)
    /* 541 */ "    ", "1/2-", "1/2+", "9/2+", "1/2+",
    /* 546 */ "5/2+", " 0+ ", "    ", " 0+ ", " 1+ ",
    /* 551 */ " 0+ ", " 3+ ", " 0+ ", "1/2-", "1/2+",
    /* 556 */ "9/2+", "1/2+", "5/2+", "1/2+", "    ",
    /* 561 */ " 0+ ", " 1+ ", " 0+ ", " 3+ ", " 0+ ",
    /* 566 */ "1/2+", "9/2+", "1/2+", "5/2+", "1/2+",
    /* 571 */ "5/2+", " 0+ ", " 1+ ", " 0+ ", " 1+ ",
    /* 576 */ " 0+ ", "    ", " 0+ ", "    ", "9/2+",
    /* 581 */ "1/2+", "5/2+", "1/2+", "5/2+", "    ",
    // AJ14 (586..630)
    /* 586 */ " 1+ ", " 0+ ", " 1+ ", " 0+ ", "    ",
    /* 591 */ " 0+ ", "9/2+", "3/2+", "5/2+", "1/2+",
    /* 596 */ "5/2+", "    ", "    ", " 0+ ", " 2- ",
    /* 601 */ " 0+ ", " 1+ ", " 0+ ", "9/2+", "11/-",
    /* 606 */ "7/2+", "1/2+", "5/2+", "1/2+", "    ",
    /* 611 */ " 0+ ", " 3- ", " 0+ ", " 2- ", " 0+ ",
    /* 616 */ "11/-", "7/2+", "1/2+", "5/2+", "1/2+",
    /* 621 */ "1/2+", "    ", " 0+ ", " 8- ", " 0+ ",
    /* 626 */ " 2- ", " 0+ ", " 1+ ", " 0+ ", "11/-",
    // AJ15 (631..675)
    /* 631 */ "7/2+", "3/2+", "5/2+", "1/2+", "1/2+",
    /* 636 */ "    ", " 0+ ", " 8- ", " 0+ ", " 1+ ",
    /* 641 */ " 0+ ", " 1+ ", " 0+ ", "    ", "    ",
    /* 646 */ "7/2+", "3/2+", "7/2+", "1/2+", "1/2+",
    /* 651 */ "    ", "3/2+", " 5+ ", " 0+ ", " 5- ",
    /* 656 */ " 0+ ", " 1+ ", " 0+ ", " 3+ ", " 0+ ",
    /* 661 */ "7/2+", "3/2+", "7/2+", "3/2+", "5/2+",
    /* 666 */ "1/2+", "1/2+", "    ", "7,8-", " 0+ ",
    /* 671 */ " 4+ ", " 0+ ", " 2- ", " 0+ ", " 2- ",
    // AJ16 (676..720)
    /* 676 */ " 0+ ", "7/2+", "3/2+", "7/2+", "3/2+",
    /* 681 */ "7/2+", "1/2+", "5/2+", "5/2+", " 0+ ",
    /* 686 */ " 4+ ", " 0+ ", " 4+ ", " 0+ ", " 1+ ",
    /* 691 */ " 0+ ", "    ", "7/2+", "3/2+", "7/2+",
    /* 696 */ "3/2+", "5/2+", "1/2+", "5/2+", "    ",
    /* 701 */ " 2- ", " 0+ ", " 5+ ", " 0+ ", " 1+ ",
    /* 706 */ " 0+ ", " 2+ ", " 0+ ", "    ", "7/2-",
    /* 711 */ "7/2+", "3/2+", "7/2+", "3/2+", "5/2+",
    /* 716 */ "1/2+", " 0+ ", " 3- ", " 0+ ", " 5- ",
    // AJ17 (721..765)
    /* 721 */ " 0+ ", " 1+ ", "    ", "    ", "7/2-",
    /* 726 */ "7/2+", "3/2+", "5/2+", "3/2+", " 0+ ",
    /* 731 */ "    ", " 0+ ", " 3- ", " 0+ ", " 1+ ",
    /* 736 */ " 0+ ", "    ", "    ", "    ", "    ",
    /* 741 */ "7/2-", "5/2+", "3/2+", "5/2+", " 0+ ",
    /* 746 */ "    ", " 0+ ", " 2- ", " 0+ ", " 2- ",
    /* 751 */ " 0+ ", " 1+ ", " 0+ ", " 1+ ", "    ",
    /* 756 */ "    ", "    ", "    ", "3/2-", "7/2+",
    /* 761 */ "7/2-", "5/2+", "3/2+", "5/2+", "X   ",
    // AJ18 (766..810)
    /* 766 */ " 0+ ", " 0- ", " 0+ ", " 6- ", " 0+ ",
    /* 771 */ " 1+ ", "    ", "7/2+", "7/2-", "5/2+",
    /* 776 */ "7/2-", "5/2+", "1/2+", " 0+ ", " 3- ",
    /* 781 */ " 0+ ", "3,4-", " 0+ ", " 4- ", " 0+ ",
    /* 786 */ "    ", "5/2-", "7/2+", "7/2-", "5/2+",
    /* 791 */ "7/2-", "5/2-", "    ", " 0+ ", " 1- ",
    /* 796 */ " 0+ ", " 5- ", " 0+ ", "    ", "    ",
    /* 801 */ "5/2-", "7/2+", "7/2-", "5/2+", "7/2-",
    /* 806 */ "    ", " 0+ ", " 1  ", " 0+ ", "0,1-",
    // AJ19 (811..855)
    /* 811 */ " 0+ ", "    ", " 0+ ", "    ", "5/2+",
    /* 816 */ "3/2-", "5/2+", "7/2-", "1/2 ", "7/2 ",
    /* 821 */ "    ", "    ", " 0+ ", " 3- ", " 0+ ",
    /* 826 */ " 2- ", " 0+ ", "    ", "    ", "3/2+",
    /* 831 */ "5/2+", "3/2+", "5/2-", "7/2 ", "    ",
    /* 836 */ "    ", "    ", " 0+ ", " 3- ", " 0+ ",
    /* 841 */ "    ", " 0+ ", " 1  ", " 0+ ", "3/2-",
    /* 846 */ "5/2+", "3/2-", "3/2+", "3/2-", "5/2 ",
    /* 851 */ "    ", " 0+ ", " 0+ ", " 0+ ", " 3- ",
    // AJ20 (856..900)
    /* 856 */ " 0+ ", " 1  ", "    ", "    ", "3/2-",
    /* 861 */ "3/2+", "3/2-", "7/2 ", "3/2 ", "    ",
    /* 866 */ " 0+ ", " 3- ", " 0+ ", " 5+ ", " 0+ ",
    /* 871 */ "5/2+", "3/2-", "3/2+", "3/2-", "7/2-",
    /* 876 */ "3/2-", "X   ", " 0+ ", " 3- ", " 0+ ",
    /* 881 */ " 5+ ", " 0+ ", "5/2-", "3/2+", "5/2+",
    /* 886 */ "7/2-", "3/2-", "7/2 ", "X   ", " 1- ",
    /* 891 */ " 0+ ", " 1+ ", " 0+ ", " 1- ", "    ",
    /* 896 */ "5/2-", "7/2-", "5/2-", "1/2+", "    ",
    // AJ21 (901..945)
    /* 901 */ " 0+ ", " 1+ ", " 0+ ", " 1+ ", " 0+ ",
    /* 906 */ "7/2+", "7/2-", "5/2-", "1/2+", "5/2-",
    /* 911 */ " 0+ ", " 0- ", " 0+ ", " 2+ ", " 0+ ",
    /* 916 */ "7/2-", "7/2+", "1/2+", "5/2-", "    ",
    /* 921 */ "    ", " 0+ ", " 3+ ", " 0+ ", " 1- ",
    /* 926 */ " 0+ ", "    ", "1/2-", "1/2+", "7/2+",
    /* 931 */ "7/2+", "5/2-", "X   ", " 0+ ", " 1- ",
    /* 936 */ " 0+ ", " 0+ ", " 0+ ", "5/2-", "1/2+",
    /* 941 */ "1/2-", "7/2+", "    ", " 0+ ", " 2- ",
    // AJ22 (946..990)
    /* 946 */ " 0+ ", " 4- ", " 0+ ", "    ", "    ",
    /* 951 */ "1/2+", "5/2-", "7/2+", "1/2-", "7/2+",
    /* 956 */ " 4- ", " 0+ ", " 1- ", " 0+ ", " 4- ",
    /* 961 */ " 0+ ", "1/2+", "7/2-", "7/2+", "5/2-",
    /* 966 */ "7/2+", "1/2-", "    ", " 0+ ", " 7- ",
    /* 971 */ " 0+ ", " 1- ", " 0+ ", "X   ", "X   ",
    /* 976 */ "X   ", " 0+ ", "9/2+", "7/2+", "7/2-",
    /* 981 */ "7/2+", "    ", "    ", "X   ", " 1+ ",
    /* 986 */ " 0+ ", " 1+ ", " 0+ ", " 5- ", "7/2+",
    // AJ23 (991..1035)
    /* 991 */ "9/2+", "7/2+", "7/2-", "5/2+", " 3- ",
    /* 996 */ " 0+ ", " 8+ ", " 0+ ", " 1- ", " 0+ ",
    /*1001 */ "1/2-", "7/2+", "9/2+", "5/2+", "    ",
    /*1006 */ " 0+ ", " 3- ", " 0+ ", " 7+ ", " 0+ ",
    /*1011 */ " 5- ", "3/2-", "7/2+", "1/2-", "5/2+",
    /*1016 */ "9/2+", "    ", " 5- ", " 0+ ", " 3- ",
    /*1021 */ " 0+ ", "    ", " 0+ ", "    ", " 0+ ",
    /*1026 */ "7/2+", "3/2-", "5/2+", "1/2-", "3/2+",
    /*1031 */ "    ", "    ", "    ", "    ", " 0+ ",
    // AJ24 (1036..1080)
    /*1036 */ " 1- ", " 0+ ", " 6- ", " 0+ ", "    ",
    /*1041 */ " 0+ ", "3/2-", "5/2+", "1/2-", "3/2+",
    /*1046 */ "    ", "    ", "    ", " 0+ ", " 1- ",
    /*1051 */ " 0+ ", " 2- ", " 0+ ", "    ", " 0+ ",
    /*1056 */ "    ", "5/2+", "3/2-", "3/2+", "    ",
    /*1061 */ "    ", "    ", "    ", " 0+ ", " 4+ ",
    /*1066 */ " 0+ ", " 1- ", " 0+ ", "    ", "9/2-",
    /*1071 */ "3/2+", "3/2-", "3/2+", "    ", "    ",
    /*1076 */ " 0+ ", " 4- ", " 0+ ", " 1- ", " 0+ ",
    // AJ25 (1081..1125)
    /*1081 */ " 2- ", "3/2-", "3/2+", "1/2-", "3/2+",
    /*1086 */ "3/2-", "1/2+", " 0+ ", " 1- ", " 0+ ",
    /*1091 */ " 1- ", " 0+ ", "(2-)", "X   ", "11/-",
    /*1096 */ "1/2-", "3/2+", "1/2-", "1/2+", "    ",
    /*1101 */ "    ", " 0+ ", " 2- ", " 0+ ", " 2- ",
    /*1106 */ " 0+ ", "    ", "1/2-", "3/2+", "1/2-",
    /*1111 */ "1/2+", "3/2-", "X   ", " 0+ ", " 2- ",
    /*1116 */ " 0+ ", " 2- ", " 0+ ", "    ", "5/2-",
    /*1121 */ "3/2+", "1/2-", "1/2+", "5/2-", "9/2 ",
    // AJ26 (1126..1170)
    /*1126 */ " 0+ ", " 1- ", " 0+ ", " 2- ", " 0+ ",
    /*1131 */ " 7+ ", " 0+ ", "X   ", "    ", "3/2-",
    /*1136 */ "1/2+", "5/2-", "9/2-", "3/2-", " 1- ",
    /*1141 */ " 0+ ", " 2- ", " 0+ ", " 5+ ", " 0+ ",
    /*1146 */ "X   ", "    ", "5/2-", "1/2+", "5/2-",
    /*1151 */ "9/2-", "5/2-", "X   ", " 0+ ", " 2- ",
    /*1156 */ " 0+ ", " 6+ ", " 0+ ", "X   ", "X   ",
    /*1161 */ "1/2-", "1/2+", "5/2-", "9/2-", "5/2-",
    /*1166 */ "X   ", "X   ", "X   ", " 0- ", " 0+ ",
    // AJ27 (1171..1215)
    /*1171 */ " 6+ ", " 0+ ", "    ", "X   ", "X   ",
    /*1176 */ "1/2+", "1/2-", "9/2-", "5/2-", "    ",
    /*1181 */ "X   ", "X   ", " 5+ ", " 0+ ", " 5+ ",
    /*1186 */ " 0+ ", " 7+ ", "X   ", "X   ", "1/2+",
    /*1191 */ "9/2+", "9/2-", "1/2-", "9/2-", "X   ",
    /*1196 */ "X   ", "4,5 ", " 0+ ", " 1- ", " 0+ ",
    /*1201 */ " 5+ ", " 0+ ", "X   ", "9/2+", "9/2-",
    /*1206 */ "9/2+", "9/2-", "1/2-", "X   ", "X   ",
    /*1211 */ " 0+ ", " 1- ", " 0+ ", " 1- ", " 0+ ",
    // AJ28 (1216..1260)
    /*1216 */ "X   ", "X   ", "    ", "9/2-", "9/2+",
    /*1221 */ "    ", "    ", "X   ", "X   ", "X   ",
    /*1226 */ " 0+ ", " 1- ", " 0+ ", "    ", " 0+ ",
    /*1231 */ " 1- ", "X   ", "X   ", "9/2-", "9/2+",
    /*1236 */ "    ", "    ", "X   ", "X   ", "X   ",
    /*1241 */ " 0+ ", " 1- ", " 0+ ", "    ", "X   ",
    /*1246 */ "X   ", "    ", "9/2+", "    ", "X   ",
    /*1251 */ " 0+ ", "    ", " 0+ ", "    ", "X   ",
    /*1256 */ "    ", "3/2+", "    ", "    ", "X   ",
    // AJ29 (1261..1305)
    /*1261 */ " 0+ ", "    ", " 0+ ", "    ", "    ",
    /*1266 */ "    ", "    ", " 0+ ", "    ", " 0+ ",
    /*1271 */ "    ", "    ", "3/2+", "1/2+", "5/2-",
    /*1276 */ "X   ", "    ", " 0+ ", "    ", " 0+ ",
    /*1281 */ "5/2-", "3/2-", "3/2+", "X   ", " 0+ ",
    /*1286 */ "    ", " 0+ ", "    ", "    ", "3/2-",
    /*1291 */ "3/2+", "5/2-", " 0+ ", " 3+ ", " 0+ ",
    /*1296 */ " 3+ ", " 0+ ", "X   ", "5/2+", "5/2-",
    /*1301 */ "3/2+", "X   ", " 0+ ", " 2- ", " 0+ ",
    // AJ30 (1306..1350)
    /*1306 */ "X   ", "X   ", "5/2+", "3/2-", "5/2-",
    /*1311 */ "X   ", " 0+ ", "2,3 ", " 0+ ", "X   ",
    /*1316 */ "X   ", "1/2+", "3/2-", "5/2+", "    ",
    /*1321 */ "X   ", " 0+ ", " 4+ ", " 0+ ", " 0+ ",
    /*1326 */ "X   ", "3/2-", "7/2-", "5/2+", "5/2+",
    /*1331 */ " 1- ", " 0+ ", " 6- ", " 0+ ", "X   ",
    /*1336 */ "3/2-", "1/2+", "5/2+", "7/2-", "X   ",
    /*1341 */ " 0+ ", " 2+ ", " 0+ ", "    ", "X   ",
    /*1346 */ "5/2+", "5/2+", "1/2+", "5/2-", "X   ",
    // AJ31 (1351..1395)
    /*1351 */ " 0+ ", " 5+ ", " 0+ ", " 3- ", " 0+ ",
    /*1356 */ "5/2+", "5/2+", "5/2-", "1/2+", "X   ",
    /*1361 */ "X   ", " 0+ ", " 1- ", " 0+ ", "X   ",
    /*1366 */ "7/2+", "5/2-", "5/2+", "3/2-", "X   ",
    /*1371 */ " 0+ ", " 6- ", " 0+ ", " 4- ", " 0+ ",
    /*1376 */ "9/2-", "5/2+", "7/2+", "3/2-", "1/2+",
    /*1381 */ "X   ", " 0+ ", " 2+ ", " 0+ ", " 2- ",
    /*1386 */ " 0+ ", "X   ", "5/2 ", "9/2-", "3/2-",
    /*1391 */ "7/2+", "X   ", " 0+ ", " 8- ", " 0+ ",
    // AJ32 (1396..1440)
    /*1396 */ "    ", " 0+ ", "1/2+", "7/2+", "9/2-",
    /*1401 */ "7/2+", "7/2+", " 0+ ", " 2- ", " 0+ ",
    /*1406 */ "    ", " 0+ ", "7/2+", "1/2+", "3/2-",
    /*1411 */ "9/2-", "X   ", " 0+ ", " 7+ ", " 0+ ",
    /*1416 */ "X   ", "7/2+", "7/2+", "1/2+", "X   ",
    /*1421 */ "9/2-", " 0+ ", " 7+ ", " 0+ ", "X   ",
    /*1426 */ "    ", "    ", "7/2+", "7/2-", "1/2+",
    /*1431 */ "    ", " 0+ ", " 0- ", " 0+ ", "    ",
    /*1436 */ "9/2+", "    ", "    ", "    ", "    ",
    // AJ33 (1441..1485)
    /*1441 */ "    ", "X   ", "X   ", "X   ", "X   ",
    /*1446 */ "X   ", "X   ", "X   ", "X   ", "X   ",
    /*1451 */ "X   ", "X   ", "X   ", "X   ", "X   ",
    /*1456 */ "X   ", "X   ", "X   ", "X   ", "X   ",
    /*1461 */ "X   ", "X   ", "X   ", "X   ", "X   ",
    /*1466 */ "X   ", "X   ", "X   ", "X   ", "X   ",
    /*1471 */ "X   ", "X   ", "X   ", "X   ", "X   ",
    /*1476 */ "X   ", "X   ", "X   ", "X   ", "X   ",
    /*1481 */ "X   ", "X   ", "X   ", "X   ", "X   "
    };

    char qMark[5] = "X   ";

    // Local variables
    int kZLo, kZHi, ajIndex;
    int jVl, jFound, kChar;
    char jPi[5];        // was an out-param; J/parity decoder still needs it locally
    char jPiChars[5]; // scratch buffer for the J^pi label

    notTabulated = 0;
    jVal = NOTDEF_INT;
    parity = 0;

    // Greatest A is 268, least is 1
    if (kA <= 0 || kA > 268) {
        notTabulated = 1;
        return;
    }
    kZLo = kZLoArr[kA + 1];
    kZHi = ajBaseArr[kA + 2] - ajBaseArr[kA + 1] + kZLo - 1;

    // Check for valid Z
    if (kZ > kZHi || kZ < kZLo) {
        notTabulated = 1;
        return;
    }

    // Calculate position in AJ array
    ajIndex = kZ - kZLo + ajBaseArr[kA + 1];
    std::memcpy(jPi, ajData[ajIndex + 1], 4);

    // Check for missing entry
    if (std::memcmp(jPi, qMark, 4) == 0) {
        notTabulated = 1;
        return;
    }

    // Decode the J value
    jVl = 0;
    std::memcpy(jPiChars, jPi, 4); jPiChars[4] = '\0';
    jFound = 0;

    bool halfIntFound = false;
    for (int i = 0; i < 4; i++) {
        if (jPiChars[i] == ' ') continue;
        if (jPiChars[i] == 'X') return;
        // Slash found so it is half-integer and is defined
        if (jPiChars[i] == '/') { jVal = jVl; halfIntFound = true; break; }
        kChar = static_cast<int>(jPiChars[i]) - static_cast<int>('0');
        if (kChar < 0 || kChar > 9) break;
        jFound = 1;
        jVl = 10 * jVl + kChar;
    }

    if (!halfIntFound) {
        // Done with field and no slash found - either nothing or integer
        if (jFound != 1) return;
        jVal = 2 * jVl;
    }

    // Now scan for the parity
    for (int i = 0; i < 4; i++) {
        if (jPiChars[i] == '+') parity = +1;
        if (jPiChars[i] == '-') parity = -1;
    }
}


// ============================================================================
//
// Sets up canned packages of parameters.
// ============================================================================

bool Reaction::applyParameterSet()
{
    Reaction& reaction = *this;
    int dataCount = 16;

    static const char8 dataSetNames[17] = { char8(""),  // index 0 unused (1-based)
        char8("CA60A"),    char8("CA60B"),    char8("PB100A"),   char8("PB100B"),   char8("PB100C"),
        char8("PBO1A"),    char8("PBO1B"),    char8("PBO1C"),    char8("PBO1D"),
        char8("ALPHA1"),   char8("ALPHA2"),   char8("ALPHA3"),
        char8("DPSA"),     char8("DPDA"),     char8("DPSB"),     char8("DPDB")
    };

    static const int gridIndices[17] = { 0, // index 0 unused
        1, 2, 3, 4, 5,
        3, 4, 5, 6,
        7, 8, 9,
        10, 11, 12, 13
    };

    // tranRGrids: column-major, row i=1..11, col J=1..13
    // Access: tranRGrids[J][i], 1-based
    static const double tranRGrids[14][12] = {
        { 0 }, // row 0 unused
        //       1        2         3        4       5        6        7        8        9       10       11
        { 0, 1.e-2,  6.e0,   12.e0,   5.e0,   5.e0,  20.e0,   .60e0,  1.50e0,  15.e0,   .50e0,  2.e0 },  // J=1
        { 0, 1.e-3,  6.e0,   20.e0,   5.e0,   5.e0,  25.e0,   .50e0,  1.70e0,  20.e0,   .50e0,  2.e0 },  // J=2
        { 0, 1.e-3,  2.50e0,  8.e0,  1.e-2,   3.e0,  45.e0,   .86e0,  1.60e0,  20.e0,   .50e0,  2.e0 },  // J=3
        { 0, 3.e-4,  3.00e0, 12.e0,  1.e-2,   3.e0,  50.e0,   .78e0,  1.80e0,  23.e0,   .50e0,  2.e0 },  // J=4
        { 0, 1.e-4,  3.50e0, 16.e0,  1.e-2,   3.e0,  55.e0,   .70e0,  2.00e0,  26.e0,   .50e0,  2.e0 },  // J=5
        { 0, 1.e-5,  4.00e0, 20.e0,  1.e-2,   3.e0,  60.e0,   .60e0,  2.20e0,  30.e0,   .50e0,  2.e0 },  // J=6
        { 0, 1.e-3,  6.00e0, 12.e0,   3.e0,   5.e0,  20.e0,   .60e0,  1.50e0,  15.e0,   .50e0,  2.e0 },  // J=7
        { 0, 1.e-4,  7.00e0, 20.e0,   3.e0,   5.e0,  25.e0,   .50e0,  1.70e0,  20.e0,   .50e0,  2.e0 },  // J=8
        { 0, 1.e-5,  8.00e0, 25.e0,   3.e0,   5.e0,  30.e0,   .40e0,  2.00e0,  24.e0,   .50e0,  2.e0 },  // J=9
        { 0, 2.e-5,  6.00e0,  5.e0,   2.e0,  12.e0,  20.e0,   0.0e0,  1.60e0,  20.e0,   .20e0,   .90e0 }, // J=10
        { 0, 2.e-6,  6.00e0,  5.e0,   2.e0,  20.e0,  20.e0,   0.0e0,  1.60e0,  20.e0,   .15e0,   .90e0 }, // J=11
        { 0, 2.e-6,  8.00e0,  8.e0,   2.e0,  12.e0,  20.e0,   0.0e0,  2.00e0,  20.e0,   .20e0,   .90e0 }, // J=12
        { 0, 5.e-7,  8.00e0,  8.e0,   2.e0,  20.e0,  20.e0,   0.0e0,  2.00e0,  20.e0,   .15e0,   .90e0 }  // J=13
    };

    // tranIGrids: column-major, [J][i], 1-based
    static const int tranIGrids[14][8] = {
        { 0 }, // row 0 unused
        //   1   2   3   4   5   6   7
        { 0, 10, 10, 10,  4,  3, 15, 10 },  // J=1
        { 0, 15, 13, 12,  4,  2, 20, 15 },  // J=2
        { 0, 10,  8, 10,  4,  5, 30, 12 },  // J=3
        { 0, 15, 10, 12,  4,  4, 35, 16 },  // J=4
        { 0, 20, 14, 16,  4,  3, 40, 20 },  // J=5
        { 0, 30, 18, 20,  4,  2, 45, 25 },  // J=6
        { 0, 40, 25, 12,  4,  2, 15, 10 },  // J=7
        { 0, 50, 30, 14,  4,  1, 20, 15 },  // J=8
        { 0, 70, 35, 16,  4,  1, 25, 20 },  // J=9
        { 0, 20, 20, 10,  4,  1, 10,  8 },  // J=10
        { 0, 20, 30, 15,  4,  1, 10,  8 },  // J=11
        { 0, 40, 40, 20,  4,  1, 12, 10 },  // J=12
        { 0, 40, 60, 30,  4,  1, 12, 10 }   // J=13
    };

    // Inelastic parameter sets
    int inelCount = 3;
    static const char8 inelSetNames[4] = { char8(""),  // index 0 unused
        char8("INELOCA1"), char8("INELOCA2"), char8("INELOCA3")
    };

    static const double inelRGrids[4][8] = {
        { 0 },
        { 0, 1.e-3,  1.e-3,  15.e0,   6.e0,  20.e0,  1.60e0,  1.e-2 },  // N=1
        { 0, 1.e-4,  1.e-4,  20.e0,   8.e0,  25.e0,  2.00e0,  5.e-3 },  // N=2
        { 0, 1.e-6,  1.e-5,  25.e0,  12.e0,  30.e0,  2.60e0,  1.e-3 }   // N=3
    };

    static const int inelIGrids[4][3] = {
        { 0 },
        { 0, 20,  6 },  // N=1
        { 0, 30,  8 },  // N=2
        { 0, 50, 10 }   // N=3
    };

    // Elastic parameter sets
    int elasCount = 3;
    static const char8 elasSetNames[4] = { char8(""),  // index 0 unused
        char8("EL1"), char8("EL2"), char8("EL3")
    };

    static const double elasRGrids[4][5] = {
        { 0 },
        { 0, 12.e0,  15.e0,  .70e0,  1.60e0 },  // N=1
        { 0, 15.e0,  20.e0,  .60e0,  1.80e0 },  // N=2
        { 0, 20.e0,  25.e0,  .50e0,  2.00e0 }   // N=3
    };

    static const int elasIGrids[4][3] = {
        { 0 },
        { 0, 15, 15 },  // N=1
        { 0, 20, 20 },  // N=2
        { 0, 25, 25 }   // N=3
    };

    // Local
    char8 dataName;
    int setIndex, i;

    // CALL nxWord ( dataName, *80, *80, *80 )
    if (nxWord(dataName.data) != 0) {
        std::printf("\n0**** DATA MUST BE FOLLOWED BY A DATA NAME\n");
        return false;
    }

    for (setIndex = 1; setIndex <= dataCount; setIndex++) {
        if (dataName == dataSetNames[setIndex]) {
            // Transfer grid words
            i = gridIndices[setIndex];

            reaction.integrationGrid.dwCutoff  = tranRGrids[i][1];
            reaction.integrationGrid.sumDensity = tranRGrids[i][2];
            reaction.integrationGrid.stepsPerUnit = tranRGrids[i][3];
            reaction.rxn.gammaSum = tranRGrids[i][4];
            reaction.rxn.gammaDif = tranRGrids[i][5];
            reaction.integrationGrid.boundAsy = tranRGrids[i][6];
            reaction.opticalPotentialParams.alMnMt = tranRGrids[i][7];
            reaction.opticalPotentialParams.alMxMt = tranRGrids[i][8];
            reaction.integrationGrid.scatAsy = -tranRGrids[i][9];
            reaction.integrationGrid.phiMid = tranRGrids[i][10];
            reaction.integrationGrid.midpointFactor = tranRGrids[i][11];

            reaction.gridData.nPhiSum  = tranIGrids[i][1];
            reaction.gridData.nPhiDifference  = tranIGrids[i][2];
            reaction.gridData.nPhiPoints  = tranIGrids[i][3];
            reaction.gridData.nPhiAdditional = tranIGrids[i][4];
            reaction.integrationGrid.lStep  = tranIGrids[i][5];
            reaction.angMom.lMaxAdditional = tranIGrids[i][6];
            reaction.integrationGrid.lMinSub = tranIGrids[i][7];

            return true;
        }
    }

    for (setIndex = 1; setIndex <= inelCount; setIndex++) {
        if (dataName == inelSetNames[setIndex]) {
            // Inelastic parameters
            reaction.opticalPotentialParams.alMnMt = 0;
            reaction.rxn.gammaSum = 5;
            reaction.integrationGrid.lStep = 1;
            reaction.integrationGrid.maxLExtrap = 5000;
            reaction.flags.excitationType = 2;
            reaction.integrationGrid.accuracyInel = inelRGrids[setIndex][1];
            reaction.integrationGrid.dwCutoff  = inelRGrids[setIndex][2];
            reaction.integrationGrid.stepsPerUnit = inelRGrids[setIndex][3];
            reaction.integrationGrid.sumDensity = inelRGrids[setIndex][4];
            reaction.integrationGrid.scatAsy = -inelRGrids[setIndex][5];
            reaction.opticalPotentialParams.alMxMt = inelRGrids[setIndex][6];
            reaction.integrationGrid.accuracy = inelRGrids[setIndex][7];
            reaction.angMom.lMaxAdditional = inelIGrids[setIndex][1];
            reaction.integrationGrid.nCoulombPoints = inelIGrids[setIndex][2];

            return true;
        }
    }

    for (setIndex = 1; setIndex <= elasCount; setIndex++) {
        if (dataName == elasSetNames[setIndex]) {
            // Elastic parameters
            reaction.integrationGrid.stepsPerUnit = elasRGrids[setIndex][1];
            reaction.integrationGrid.scatAsy = -elasRGrids[setIndex][2];
            reaction.opticalPotentialParams.alMnMt = elasRGrids[setIndex][3];
            reaction.opticalPotentialParams.alMxMt = elasRGrids[setIndex][4];
            reaction.integrationGrid.lMinSub = elasIGrids[setIndex][1];
            reaction.angMom.lMaxAdditional = elasIGrids[setIndex][2];
            return true;
        }
    }

    // Invalid data name
    std::printf("\n0**** INVALID DATA NAME: %.8s.  ALLOWABLE DATA NAMES ARE:\n", dataName.data);
    std::printf("  ");
    for (int k = 1; k <= dataCount; k++) std::printf("%.8s  ", dataSetNames[k].data);
    for (int k = 1; k <= inelCount; k++) std::printf("%.8s  ", inelSetNames[k].data);
    for (int k = 1; k <= elasCount; k++) std::printf("%.8s  ", elasSetNames[k].data);
    std::printf("\n");
    return false;
}


// ============================================================================
//
// Scans a reaction of the general form
//    40CA(O16, 12C(2+ 3.26))44TI(7.23 3/2)
// and stores Z's, M's, J's, E*'s, and mass excesses of the 5 particles.
// ============================================================================

bool Reaction::parseReactionString()
{
    Reaction& reaction = *this;
    // Local arrays
    char8 guy[5];        // guy(4), 1-based
    double eStars[5];    // eStars(4), 1-based
    int nodeVals[5];       // nodeVals(4), 1-based
    int lVals[5];        // lVals(4), 1-based
    int jVals[5];        // jVals(4), 1-based
    int iParities[5];       // iParities(4), 1-based


    static const int particleIndices[5] = { 0, 3, 1, 2, 4 };  // 1-based

    static const char names[6][5] = { "", "A   ", "B   ", "BIGA", "BIGB", "X   " };

    int inChStart, returnCode, partIndex, ia, iz, n, numParticles;
    int jVal, gsParity;
    double atomicMassExcess;

    // Get past any initial delimiters (assume the reaction is in upper case)
    while (true) {
        char myChar = inputBuffer.iBuf[inputBuffer.inCh];
        if ((myChar >= 'A' && myChar <= 'Z') ||
            (myChar >= 'a' && myChar <= 'z') ||
            (myChar >= '0' && myChar <= '9')) break;
        inputBuffer.inCh = inputBuffer.inCh + 1;
        if (inputBuffer.inCh > inputBuffer.nOch) {
            std::printf("\n0**** REACTION KEYWORD AND COMPLETE REACTION MUST BE ON ONE LINE.\n");
            return false;
        }
    }

    // Start of reaction found. Break it into 4 pieces.
    inChStart = inputBuffer.inCh;
    qvScan(guy, eStars, nodeVals, lVals, jVals, iParities, returnCode, reaction);
    if (returnCode != 0) return false;

    // Now get the A and Z values
    for (int i = 1; i <= 4; i++) {
        azCode(guy[i].data, iz, ia, returnCode);
        if (returnCode != 0) {
            if (returnCode != -2) std::printf("\n0**** A SYMBOL HAS INCORRECT SYNTAX:  SYMBOL = %.8s\n", guy[i].data);
            if (returnCode == -2) std::printf("\n0**** THE SYMBOL %.6s IS NOT A KNOWN ELEMENT.\n", guy[i].data);
            return false;
        }
        partIndex = particleIndices[i];
        reaction.charges.zArray[partIndex] = iz;
        aMs(reaction)[partIndex - 1] = ia;
        if (eStars[i] != reaction.internalState.undefValue) reaction.energies.exs[partIndex] = eStars[i] + 1.e-30;
        if (nodeVals[i] != NOTDEF_INT) reaction.internalState.nodesP[partIndex] = nodeVals[i];
        if (lVals[i] != NOTDEF_INT) reaction.internalState.lSpecs[partIndex] = lVals[i];
        if (jVals[i] != NOTDEF_INT) reaction.angMom.js[partIndex] = jVals[i];
        if (iParities[i] != 0) reaction.angMom.parities[partIndex] = iParities[i];
    }

    // Find X
    reaction.charges.zArray[5] = std::abs(reaction.charges.zArray[1] - reaction.charges.zArray[2]);
    aMs(reaction)[4] = std::fabs(aMs(reaction)[0] - aMs(reaction)[1]);
    numParticles = (aMs(reaction)[4] == 0) ? 4 : 5;

    // Now get ground state mass excess and determine if the nuclei exist
    for (int i = 1; i <= numParticles; i++) {
        readParticleAZN(reaction, i, ia, iz, n);
        atomicMassExcess = excess(iz, ia, returnCode);
        if (returnCode != 0) {
            if (i == 5) {
                std::printf("\n0**** THE EXCHANGED CLUSTER HAS A =%4d, Z =%4d,  N =%4d  AND IS NOT BOUND.\n"
                            "      MERELY TO CONCEIVE OF SUCH THINGS MAKES THEM APPEAR RIDICULOUS.\n",
                            ia, iz, n);
                return false;
            }
            std::printf("\n0**** THE NUCLEUS WITH A =%4d,  Z =%4d,  N =%4d  IS NOT BOUND "
                        "ACCORDING TO PTOLEMY\"S MASS TABLE.\n", ia, iz, n);
            continue;
        }

        // Print warning for absurd exchanges (only when i==5)
        if (i == 5 && !(ia <= 2 || iz != 0 || iz != ia)) {
            std::printf("\n0**** WARNING:  THE EXCHANGED CLUSTER IS ABSURB,"
                        " BUT IF THAT'S WHAT YOU WANT, SO BE IT.\n");
        }
        reaction.masses.amxcgs[i] = atomicMassExcess;
    }

    // Now get J for each of the 5
    for (int i = 1; i <= numParticles; i++) {
        // If E* has been entered then G.S. J is meaningless
        if (reaction.energies.exs[i] > 1.e-5) continue;

        readParticleAZN(reaction, i, ia, iz, n);
        groundStateInfo(iz, ia, jVal, gsParity, returnCode);
        if (!(returnCode == 0 && jVal != NOTDEF_INT)) {
            // No error if J is already defined or is explicitly given in the reaction
            if (reaction.angMom.js[i] != reaction.internalState.notDefSentinel) {
                clearExcitationSentinel(reaction, i);
                continue;
            }

            // Dineutron or diproton is J = 0
            if (ia == 2 && iz != 1) {
                jVal = 0;
                assignGroundStateSpin(reaction, i, jVal, gsParity);
                continue;
            }
            std::printf("\n0**** WARNING:  THE GROUND STATE SPIN OF THE NUCLEUS WITH A =%4d,"
                        "  Z =%4d,  N =%4d  IS NOT KNOWN TO PTOLEMY.\n", ia, iz, n);
            clearExcitationSentinel(reaction, i);
            continue;
        }

        // If J has been defined to be different from G.S. then need E*
        if (reaction.angMom.js[i] == reaction.internalState.notDefSentinel) {
            assignGroundStateSpin(reaction, i, jVal, gsParity);
            continue;
        }
        if (static_cast<int>(reaction.angMom.js[i]) == jVal || reaction.energies.exs[i] != 0) {
            clearExcitationSentinel(reaction, i);
            continue;
        }
        std::printf("\n0**** WARNING:  GROUND STATE SPIN FOR PARTICLE %.4s IS%3d/2."
                    "  YOU HAVE SPECIFIED J%.4s =%3d/2, BUT HAVE NOT YET SPECIFIED E*%.4s.\n",
                    names[i], jVal, names[i], static_cast<int>(reaction.angMom.js[i]), names[i]);
    }

    // It was successful, copy reaction for header
    for (int i = 1; i <= 45; i++) {
        reaction.reactStr[i] = ' ';
    }
    n = std::min(45, inputBuffer.inCh - inChStart);
    for (int i = 1; i <= n; i++) {
        reaction.reactStr[i] = inputBuffer.iBuf[inChStart - 1 + i];
    }
    return true;
}
