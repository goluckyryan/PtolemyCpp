#pragma once
// constants.h — Physical and numerical constants for Ptolemy
// All constexpr; no runtime initialization needed.

namespace Constants {
    constexpr double PI       = 3.14159265358979320;
    constexpr double RADIAN   = 180.0 / PI;          // degrees per radian
    constexpr double DEGREE   = PI / 180.0;          // radians per degree
    constexpr double hbar_c   = 197.328580;           // MeV·fm
    constexpr double amu_MeV  = 931.50160;            // atomic mass unit in MeV/c²
    constexpr double fine_structure_inv = 137.036040;  // 1/α (Ptolemy convention: AFINE=137.036)
    constexpr double bigNum   = 1.0E+30;              // sentinel infinity
    constexpr double smlNum   = 1.0 / bigNum;         // near-zero sentinel
    constexpr double BIGLOG   = 69.07755278982137;    // log(bigNum)
}
