#pragma once
// MathTables.h — sqrt-factorial and log-factorial tables (accessor functions)
// Used by CG, 6J, 9J, Racah coefficient routines.
//
// Storage is hidden inside MathTables.cpp as static-local-init objects to
// avoid header-extern global mutable state. Callers reach the tables via
// the factorialTable() / logFactorialTable() accessor functions below
// (same names as the former extern variables, just with () appended).

// sqrt(n!) table for Clebsch-Gordan coefficients.
// maxFactorial = 96; factTable[1..97] = sqrt((n-1)!) for n=1..97 (1-based);
// factTable[0] is the unused 0th slot.
struct FactorialTable {
    int    maxFactorial;
    double factTable[98];
};

// log(n!) table for Racah / CG coefficients.
// maxLf tracks how far lf[] has been populated by setLog();
// lf[1] = log(0!) = 0 is the seed.
struct LogFactorialTable {
    int    maxLf;
    double lf[2100];
};

// Accessor functions — call once per use. Static-local-init guarantees a
// single instance per program, initialized on first call. No header-extern.
const FactorialTable& factorialTable();
LogFactorialTable&    logFactorialTable();

// GaussLegendreRoots / gaussRoots / ROOTI sqrt-integer table dropped

// Extend the log-factorial table up to index maxIndex (fills logFactorialTable()).
void setLog(int maxIndex);
