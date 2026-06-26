// math_special.cpp — special mathematical functions.
//   dGamma     — complete Gamma function Γ(x).
//   dLogGamma  — natural log of the Gamma function, ln Γ(x).

#include <cstdio>
#include <cmath>
#include "MathTables.h"

// ============================================================================
// FUNCTION DGAMMA + DLGAMA — complete Gamma function and its natural log.
//
// Two entry points share the log-gamma front half; the selector below
// picks the raw vs log return.
// ============================================================================

static double dgammaImpl(double x, bool logOnly) {
    bool lNeg, lExp;
    double z, g;

    const char* name = logOnly ? "DLGAMA" : "DGAMMA";  // must be per-call: 'static' would freeze the name to the first caller's flag
    static const double PI    = 3.141592653589793240e0;
    static const double R2PILN = 0.918938533204672420e0;
    static const double expLimit = 695.0e0;  // biggest x s.t. exp(x) doesn't overflow

    // Log Gamma front half.
    z = x;
    lExp = false;
    lNeg = z < 0.50e0;
    if (lNeg) z = 1.0e0 - z;

    // If z < 12, transform to range 0.5 <= z < 1.5 (A&S 6.1.15)
    // and use 19-term series expansion (A&S 6.1.34).
    if (z < 12.0e0) {
        double h = 1.0e0;
        int n = (int)(z - 0.50e0);
        for (int i = 1; i <= n; i++) {
            z = z - 1.0e0;
            h = h * z;
        }
        double y = z - 1.0e0;
        g = -.000020134854780700e0;
        if (std::fabs(y) > 0.10e0)
            g = g
                + y * (-.00000125049348210e0
                + y * ( 0.00000113302723200e0
                + y * (-.00000020563384170e0
                + y * ( 0.00000000611609500e0
                + y * ( 0.00000000500200750e0
                + y * (-.00000000118127460e0
                + y * ( 0.00000000010434270e0
                )))))));
        g = 1.0e0
            + y * ( 0.577215664901532900e0
            + y * (-.655878071520253800e0
            + y * (-.042002635034095200e0
            + y * ( 0.166538611382291500e0
            + y * (-.042197734555544300e0
            + y * (-.009621971527877000e0
            + y * ( 0.007218943246663000e0
            + y * (-.001165167591859100e0
            + y * (-.000215241674114900e0
            + y * ( 0.000128050282388200e0
            + y * g ))))))))));
        g = h / g;
        lExp = true;
    } else {
        // label_200: Calculate log gamma using Bernoulli expansion.
        double zI = 1.0e0 / z;
        double zis = zI * zI;
        g = (-1.0e0 / 360.0e0);
        if (z < 164.0e0)
            g = g
                + zis * ((1.0e0 / 1260.0e0)
                + zis * ((-1.0e0 / 1680.0e0)
                + zis * ((1.0e0 / 1188.0e0)
                + zis * ((-691.0e0 / 312840.0e0)
                ))));
        g = R2PILN - z + (z - 0.50e0) * std::log(z)
            + zI * ((1.0e0 / 12.0e0)
            + zis * g);
    }

    // label_250: If x < 0.5, transform using A&S 6.1.17.
    if (lNeg) {
        double t = std::sin(PI * z);
        if (t == 0.0e0) {
            printf("\n **** ERROR IN %6s:  ARGUMENT = %25.16G ****\n"
                   " **** ARGUMENT IS A NEGATIVE INTEGER OR ZERO ****\n", name, x);
        }
        if (lExp) {
            g = PI / (g * t);
            lNeg = false;
        } else {
            g = std::log(PI / std::fabs(t)) - g;
            lNeg = t < 0;
        }
    }

    if (logOnly) {
        if (lExp) g = std::log(g);
        return g;
    }
    if (!lExp) {
        if (g > expLimit) {
            printf("\n **** ERROR IN %6s:  ARGUMENT = %25.16G ****\n"
                   " **** EXPONENT OVERFLOW ****\n", name, x);
        }
        g = std::exp(g);
        if (lNeg) g = -g;
    }
    return g;
}

double dGamma(double x) { return dgammaImpl(x, /*logOnly=*/false); }
double dLogGamma(double x) { return dgammaImpl(x, /*logOnly=*/true);  }  // was DLGAMA

// ============================================================================
// folded in from source_misc.cpp: setLog, second
// ============================================================================
void setLog(int maxIndex)
{
    // Set up log factorial table
    // lf[i+1] = lf[i] + log(i) for i=startIndex..maxIndex; max index = maxIndex+1.
    // lf array has 2100 slots (indices 0..2099); guard against overflow.
    static const int lfMax = 2098;  // max maxIndex: lf[2099] is last valid write
    if (maxIndex > lfMax) {
        std::printf("setLog: N=%d exceeds LF table limit %d; clamping.\n", maxIndex, lfMax);
        maxIndex = lfMax;
    }
    if (maxIndex <= logFactorialTable().maxLf) return;
    int startIndex = logFactorialTable().maxLf + 1;
    logFactorialTable().maxLf = maxIndex;
    for (int i = startIndex; i <= maxIndex; i++) {
        logFactorialTable().lf[i + 1] = logFactorialTable().lf[i] + std::log((double)i);
    }
}

double second() { return 0.0; } // timing stub

