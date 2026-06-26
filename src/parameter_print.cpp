// parameter_print.cpp — parameterPrint: prints one potential-term row
// (suppressed when vPar=0, except i=3,4 spin-orbit).

#include "Reaction.h"
#include "print_utils.h"
#include <cstdio>
#include <cmath>

// Fortran Gw.d format: prints val in a field of width w with d significant digits.
// When 0.1 <= |val| < 10^d: uses F(w-4).(d-e) followed by 4 trailing spaces,
// where e = number of digits before the decimal point (0 if |val| < 1).
// Otherwise: uses E(w).(d-1) format.
void print_G(int w, int d, double val)
{
    double absVal = std::fabs(val);
    double limit  = std::pow(10.0, (double)d);
    if (absVal == 0.0 || (absVal >= 0.1 && absVal < limit)) {
        int e = 0;
        if (absVal >= 1.0)
            e = (int)std::floor(std::log10(absVal)) + 1;
        int decimals = d - e;
        if (decimals < 0) decimals = 0;
        // F(w-4).decimals + 4 trailing spaces = total w chars
        std::printf("%*.*f    ", w - 4, decimals, val);
    } else {
        // Fortran E notation: 0.dddddE±ee (exponent is 1 higher than C's)
        int exponent = (int)std::floor(std::log10(absVal)) + 1;
        double mantissa = val / std::pow(10.0, (double)exponent);
        std::printf("%*.*fE%+03d", w - 4, d, mantissa, exponent);
    }
}

// parameterPrint: print parameters for one potential term.
//
// i        term index in linkuleAddr (0 = no linkule possible).
// name     18-char label for the term (right-padded with spaces).
// vPar     coupling constant.
// rPar     radius.
// aPar     diffuseness.
// rPar0    radius parameter (R/r0Mass).
//  P5=1.0,"POWER" (i=1,2: REAL CENTRAL / VOLUME ABSORPTION) or
//  P5=0.0,"" (i=3,4,5,13: SO/Coulomb/surface). With POWRL/POWIM gone
//  DEFALT-only-1.0 and TAU/TAUI gone DEFALT-only-0, the trailing
//
// Suppression rule (no LINKULE case): if vPar == 0, return.

void parameterPrint(int i, const char* name, double vPar, double rPar, double aPar,
            double rPar0, Reaction& reaction)
{
    // If a linkule exists for this term, call it to print its own info.
    if (i > 0 && i <= numLinkules && reaction.linkuleData.linkuleAddr[i][3] != 0) {
        // Blank line before (FORMAT '0'), then header
        std::printf("\n%.18s POTENTIAL IS BEING COMPUTED BY THE %.8s LINKULE:\n",
                    name, (char*)&reaction.linkuleData.linkuleAddr[i][1]);
        // Full LINKUL call with mode=2 (print) would go here.
        // Linkule printing not yet fully implemented — fall through to standard print.
        std::printf(" \n");
        return;
    }

    // No linkule: print 4 parameters, unless vPar is zero.
    if (vPar == 0.0) return;

    // Leading space (Fortran carriage control ' ') + name (a18)
    std::printf(" %-18.18s", name);
    // vPar in G14.5
    print_G(14, 5, vPar);
    // 3 fill spaces to reach T37, then rPar in F7.4
    std::printf("   %7.4f", rPar);
    // 5 fill spaces to reach T49, then aPar in F7.4
    std::printf("     %7.4f", aPar);
    // 11 fill spaces to reach T67, then rPar0 in F7.4
    std::printf("           %7.4f\n", rPar0);
}
