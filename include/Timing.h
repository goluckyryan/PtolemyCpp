#pragma once
// timing.h — CPU timing and performance counters
//
// TODO: Replace float times[] with std::chrono. Currently kept because
//       Ptolemy prints timing breakdowns as part of standard output
//       (input_reader.cpp L861-871). Changing format would break
//       bit-identical output tests. Modernize when output format is decoupled.

struct Timing {
    float  times[9];   // 1-based [1..8]: elapsed CPU times for sub-phases (seconds)
};

// Timing/date free functions.
double dtime_();
double second();        // timing stub returning 0.0 — Fortran-era CPU-time placeholder
void   getDate(char* date);
