// MathAngular.cpp — angular-momentum coefficients (CG, 3J, 6J, Racah, 9J),
// continued fractions, Coulomb sigma, Gauss-Legendre, getDate.

#include "Timing.h"
#include "math/special.h"
#include "math/angular_momentum_coeff.h"
#include "math/continued_fraction.h"
#include "math/coulomb_sigma.h"
#include "math/gauss_quadrature.h"
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <algorithm>
#include "MathTables.h"


// the table was populated at static init but no caller ever read it. The
// Legendre recursions that originally consumed it must have been refactored
// or deleted upstream of this cleanup.

// ============================================================================
// cfracInit / cfracEval — complex continued fraction interpolation
// Fortran lines 1-174
// ============================================================================

// File-static saved variable (Fortran SAVE nMax)
static int cfracMaxIndex = 0;

void cfracInit(int nPts, complex16* xs, complex16* ys) {
    // 0-based version (matches Maple). Callers pass plain .data() pointers.
    // WARNING: the earlier 1-based version read ys[nPts] one past the end of a
    // nPts-sized array (ASan OOB); keep this 0-based to match the callers.
    complex16 d, xJ, yJ;
    double comp, temp;
    int k = 0, i, j;

    cfracMaxIndex = nPts - 1;

    // Search for the largest Y
    comp = -1.0;
    for (i = 0; i < nPts; i++) {
        d = ys[i];
        temp = std::norm(d);
        if (temp <= comp) continue;
        comp = temp;
        k = i;
    }
    if (cfracMaxIndex <= 0) return;

    // The first coefficient is the biggest Y — swap it to position 0.
    yJ = ys[k];
    if (k != 0) {
        ys[k] = ys[0];
        ys[0] = yJ;
        xJ = xs[k];
        xs[k] = xs[0];
        xs[0] = xJ;
    }

    // Candidates for the second coefficient.
    comp = -1.0;
    for (i = 1; i < nPts; i++) {
        d = 1.0 - yJ / ys[i];
        ys[i] = d;
        temp = std::norm(d);
        if (temp <= comp) continue;
        comp = temp;
        k = i;
    }

    // Loop through the rest of the coefficients.
    for (j = 1; j < nPts; j++) {
        xJ = xs[k];
        yJ = ys[k];
        if (k != j) {
            ys[k] = ys[j];
            xs[k] = xs[j];
            xs[j] = xJ;
        }
        yJ = yJ / (xs[j - 1] - xJ);
        ys[j] = yJ;
        if (j == nPts - 1) break;

        // Quit if this coefficient is tiny.
        if (comp < 1.0e-14) {
            cfracMaxIndex = j - 1;
            std::printf("\n**** WARNING:  CONTINUED FRACTION USED ONLY%3d OUT OF%3d POINTS.\n",
                       cfracMaxIndex, nPts);
            break;
        }

        // Candidates for the next coefficient.
        {
            comp = -1.0;
            xJ = xs[j - 1];
            for (i = j + 1; i < nPts; i++) {
                d = 1.0 + yJ * (xs[i] - xJ) / ys[i];
                ys[i] = d;
                temp = std::norm(d);
                if (temp <= comp) continue;
                comp = temp;
                k = i;
            }
        }
    }
    return;
}

void cfracEval(complex16* xs, complex16* ys, complex16 x, complex16& y) {
    // 0-based version (matches Maple). y = ys[0] / (1 + y) at the end.
    y = complex16(0.0, 0.0);
    if (cfracMaxIndex >= 1) {
        for (int j = 1; j <= cfracMaxIndex; j++) {
            int k = cfracMaxIndex - j;
            y = ys[k + 1] * (x - xs[k]) / (1.0 + y);
        }
    }
    y = ys[0] / (1.0 + y);
    return;
}


// ============================================================================
// Begin clebschGordan/THRJ/THREEJ and SIXJ/RACAH translations
// ============================================================================


// ============================================================================
// clebschGordan / THRJ / THREEJ
// ============================================================================
//
// IMPLICIT INTEGER*4 (A-C, i-Z), REAL*8 (D-H)
//   Variables A-C, i-Z are int; D-H are double.
//   sum, answer are explicitly REAL*8.
//
//   clebschGordan(A,B,X,Y,cIn,zIn)  — Clebsch-Gordan coefficient
//   THRJ(A,B,cIn,X,Y,zIn)    — 3-J symbol (z = -zIn)
//   THREEJ same as THRJ
//
// Handled via static helper with isThreeJ flag.
// DBLSW is always true (double precision only).

static double clebschGordanImpl(int a, int b, int x, int y, int cIn, int zIn,
                           bool isThreeJ)
{
    // ---- All local variable declarations before any gotos ----

    // IMPLICIT INTEGER*4 (A-C, i-Z)
    int c, z;
    int x1, y1, z1, t30;
    int x2, y2, z2;
    int r;
    int a1, a2, a3;
    int b1, b2;
    int kMin, kMax, kMaxM2;
    int n;
    int jSum;
    int temp;
    int i;

    // IMPLICIT REAL*8 (D-H)
    double b1Double, b2Double;
    double kMaxDouble, kMinDouble, cDouble;
    double eLog, e1, e2, e3, e4, e5, e6;
    double f1, f2, f3, f4;
    double g1, g4;
    double h1, h4;
    double jSumDouble;

    // REAL*8 (explicit)
    double sum, answer;

    //   dRay(1)=a1Double, dRay(2)=a2Double, dRay(3)=a3Double,
    //   dRay(4)=y1Double, dRay(5)=x2Double, dRay(6)=y2Double,
    //   dRay(7)=x1Double, dRay(8)=z1Double, dRay(9)=z2Double
    double dRay[10]; // 1-based: dRay[1]..dRay[9]
    double& a1Double = dRay[1];
    double& a2Double = dRay[2];
    double& a3Double = dRay[3];
    double& y1Double = dRay[4];
    double& x2Double = dRay[5];
    double& y2Double = dRay[6];
    double& x1Double = dRay[7];
    double& z1Double = dRay[8];
    double& z2Double = dRay[9];

    // ---- Begin translated code ----


    if (isThreeJ) {
        z = -zIn;
    } else {
        z = zIn;
    }
    c = cIn;

    if (z != x+y) return 0.0;

    x1 = ((a+x) >> 1);
    y1 = ((b+y) >> 1);
    z1 = ((c+z) >> 1);

    t30 = 30;

    if ((((x1) >> (t30)) & 1) || (((y1) >> (t30)) & 1) || (((z1) >> (t30)) & 1)) return 0.0;

//      x2 = x1 - X
    x2 = x1 - x;
    y2 = y1 - y;
    z2 = z1 - z;

//      r = x1 + y1 + z2 + 1
    r = x1 + y1 + z2 + 1;

//      a1 = r-1 - A
    a1 = r-1 - a;
    a2 = r-1 - b;
    a3 = r-1 - c;

//      b1 = a1 - y1
    b1 = a1 - y1;
    b2 = a2 - x2;

//      a1Double = a1
    a1Double = a1;
    a2Double = a2;
    a3Double = a3;
    y1Double = y1;
    x2Double = x2;

//      b1Double = a1Double - y1Double
    b1Double = a1Double - y1Double;
    b2Double = a2Double - x2Double;

//      sum = 1
    sum = 1;

//      kMin = MAX0 (0, -b1, -b2)
    kMin = std::max({0, -b1, -b2});
//      kMax = MIN0 (a3, x2, y1)
    kMax = std::min({a3, x2, y1});

    if (kMax - kMin < 0) return 0.0;

    if (kMax - kMin > 0 && std::abs(x)+std::abs(y) == 0) {
        // Special case: M1 = M2 = M3 = 0
        // RESULT IS 0 IF J1+J2+J3 IS ODD
        if (!(((r) >> (0)) & 1)) return 0.0;
        jSum = r-1;
        if (r <= factorialTable().maxFactorial) {
            answer = ((factorialTable().factTable[1+a1]*factorialTable().factTable[1+a2])/factorialTable().factTable[1+r])
              * ((factorialTable().factTable[2+c]*factorialTable().factTable[1+a3])/factorialTable().factTable[1+c])
              * (factorialTable().factTable[1+jSum/2]
                 / (factorialTable().factTable[1+a1/2]*factorialTable().factTable[1+a2/2]*factorialTable().factTable[1+a3/2]))
              * (factorialTable().factTable[1+jSum/2]
                 / (factorialTable().factTable[1+a1/2]*factorialTable().factTable[1+a2/2]*factorialTable().factTable[1+a3/2]));
        } else if (r <= logFactorialTable().maxLf) {
            eLog = .5 * (logFactorialTable().lf[1+a1] + logFactorialTable().lf[1+a2]
                      + logFactorialTable().lf[1+a3]
                + logFactorialTable().lf[2+c] - logFactorialTable().lf[1+c]
                      - logFactorialTable().lf[1+r] )
              - logFactorialTable().lf[1+((a1) >> 1)]
                      - logFactorialTable().lf[1+((a2) >> 1)]
              - logFactorialTable().lf[1+((a3) >> 1)]
                      + logFactorialTable().lf[1+((jSum) >> 1)];
            answer = std::exp(eLog);
        } else {
            jSumDouble = a1Double + a2Double + a3Double;
            eLog = std::log(static_cast<double>(c+1)) - dLogGamma(jSumDouble+2);
            e2 = -dLogGamma(.5*jSumDouble+1);
            for (i = 1; i <= 3; i++) {
                eLog = eLog + dLogGamma(1+dRay[i]);
                e2 = e2 + dLogGamma(1 + .5*dRay[i]);
            }
            eLog = .5*eLog - e2;
            answer = std::exp(eLog);
        }
        if ((((a3) >> (1)) & 1)) answer = -answer;
    } else {
        // General case (M1,M2 != 0, or kMax == kMin)
        if (kMax - kMin > 0) {
            kMaxDouble = kMax;
            e1 = -kMaxDouble;
            e2 = -b1Double-kMaxDouble;
            e3 = -b2Double - kMaxDouble;
            e4 = a3Double - (kMaxDouble-1);
            e5 = x2Double - (kMaxDouble-1);
            e6 = y1Double - (kMaxDouble-1);

            f1 = e1*e2*e3;
            f4 = e4*e5*e6;

            kMaxM2 = kMax - 2;
            sum = 1 + f4/f1;

            if (kMaxM2 >= kMin) {
                g1 = (e1+1)*(e2+1)*(e3+1) - f1;
                g4 = (e4+1)*(e5+1)*(e6+1) - f4;
                h1 = (e1+2)*(e2+2)*(e3+2) - g1 - g1 - f1;
                h4 = (e4+2)*(e5+2)*(e6+2) - g4 - g4 - f4;

                for (n = kMin; n <= kMaxM2; n++) {
                    f1 = f1 + g1;
                    f4 = f4 + g4;
                    g1 = g1 + h1;
                    g4 = g4 + h4;
                    h1 = h1 + 6;
                    h4 = h4 + 6;
                    sum = 1 + sum*(f4/f1);
                }
            }
        }

        // Compute answer from factorial tables, log-factorial, or DLGAMA
        if (r <= factorialTable().maxFactorial) {
            f1 = factorialTable().factTable[1+a1];
            f2 = factorialTable().factTable[1+a2];
            f3 = factorialTable().factTable[1+a3];
            eLog = (f3/(factorialTable().factTable[1+kMin]*factorialTable().factTable[1+a3-kMin]))
              * (f2/(factorialTable().factTable[1+x2-kMin]*factorialTable().factTable[1+b2+kMin]))
              * (f1/(factorialTable().factTable[1+y1-kMin]*factorialTable().factTable[1+b1+kMin]));
            sum = std::round(eLog*sum*eLog);
            answer = sum * (factorialTable().factTable[1+x1]*factorialTable().factTable[1+x2]/(f3*f2))
              * (factorialTable().factTable[1+y1]*factorialTable().factTable[1+y2]/f1)
              * (factorialTable().factTable[2+c]/factorialTable().factTable[1+c])
              * (factorialTable().factTable[1+z1]*factorialTable().factTable[1+z2]/factorialTable().factTable[1+r]);
        } else if (r <= logFactorialTable().maxLf) {
            eLog = logFactorialTable().lf[1+a3] - logFactorialTable().lf[1+kMin]
                  - logFactorialTable().lf[1+a3-kMin]
              + logFactorialTable().lf[1+a1] - logFactorialTable().lf[1+b1+kMin]
                  - logFactorialTable().lf[1+y1-kMin]
              + logFactorialTable().lf[1+a2] - logFactorialTable().lf[1+b2+kMin]
                  - logFactorialTable().lf[1+x2-kMin];
            e2 = .5*(-logFactorialTable().lf[1+a1] - logFactorialTable().lf[1+a2]
                     - logFactorialTable().lf[1+a3]
               + logFactorialTable().lf[1+x1] + logFactorialTable().lf[1+x2]
               + logFactorialTable().lf[1+y1] + logFactorialTable().lf[1+y2]
                     + logFactorialTable().lf[1+z1]
               + logFactorialTable().lf[1+z2] - logFactorialTable().lf[1+r]
               + logFactorialTable().lf[2+c] - logFactorialTable().lf[1+c] );
            if (std::fabs(eLog) > 40) {
                eLog = eLog + e2;
                answer = std::exp(eLog)*sum;
            } else {
                sum = std::round(sum*std::exp(eLog));
                answer = std::exp(e2)*sum;
            }
        } else {
            // Large J: DLGAMA
            x1Double = a2Double + a3Double - x2Double;
            cDouble = c;
            kMinDouble = kMin;
            z1Double = x1Double + y1Double - a3Double;
            y2Double = a1Double + a3Double - y1Double;
            z2Double = cDouble - z1Double;
            eLog = std::log(cDouble+1) - dLogGamma(2+x1Double+y1Double+z2Double);
            for (i = 1; i <= 9; i++) {
                eLog = eLog + dLogGamma(1+dRay[i]);
            }
            eLog = .5*eLog;
            a1Double = b1Double + (kMinDouble+1);
            a2Double = b2Double + (kMinDouble+1);
            a3Double = a3Double - (kMinDouble-1);
            y1Double = y1Double - (kMinDouble-1);
            x2Double = x2Double - (kMinDouble-1);
            y2Double = (kMinDouble+1);
            for (i = 1; i <= 6; i++) {
                eLog = eLog - dLogGamma(dRay[i]);
            }
            answer = std::exp(eLog) * sum;
        }
        if ((((kMin) >> (0)) & 1)) answer = -answer;
    }

    // Apply 3-J factor if needed
    if (isThreeJ) {
        temp = a2-z2;
        if ((((temp) >> (0)) & 1)) answer = -answer;
        if (c > factorialTable().maxFactorial-1) {
            answer = answer / std::sqrt(c+1.0);
        } else {
            answer = answer * (factorialTable().factTable[1+c]/factorialTable().factTable[2+c]);
        }
    }

    return answer;
}

// ---- Public entry points ----

// clebschGordan(A,B,X,Y,cIn,zIn)
double clebschGordan(int a, int b, int x, int y, int cIn, int zIn)
{
    return clebschGordanImpl(a, b, x, y, cIn, zIn, false);
}

// Fortran ENTRY THRJ(A,B,cIn,X,Y,zIn) reordered args vs clebschGordan; the impl's
// isThreeJ path takes (A, B, X, Y, cIn, zIn, true). threeJ is the live caller.
double threeJ(int a, int b, int cIn, int x, int y, int zIn)
{
    return clebschGordanImpl(a, b, x, y, cIn, zIn, true);
}


// ============================================================================
// SIXJ / RACAH
// ============================================================================
//
// IMPLICIT INTEGER*4 (A-C, i-Z), REAL*8 (D-H)
//   Variables A-C, i-Z are int; D-H are double.
//   sum, answer, SIXJ, RACAH are explicitly REAL*8.
//
//   SIXJ(A,B,C,X,Y,Z)    — 6-J coefficient
//   RACAH(A,B,Y,X,C,Z)   — Racah W coefficient
//     Fortran: ENTRY RACAH(A,B,Y,X,C,Z)
//     The parameter names in the ENTRY match the FUNCTION's names,
//     so RACAH's 1st arg -> A, 2nd -> B, 3rd -> Y, 4th -> X, 5th -> C, 6th -> Z
//
// DBLSW is always true (double precision only).

static double sixjImpl(int a, int b, int c, int x, int y, int z,
                          bool isRacah)
{
    // ---- All local variable declarations before any gotos ----

    // IMPLICIT INTEGER*4 (A-C, i-Z)
    int a1, a2, a3, a4;
    int b1, b2;
    int r;
    int kMin, kMax, kMaxM2;
    int n;
    int temp;
    int i;

    // IMPLICIT REAL*8 (D-H)
    double b1Double, b2Double;
    double rDouble, kMinDouble, kMaxDouble, yDouble;
    double d1, d2, d3, d4, d5, d6, d7, d8;
    double e, e1, e2, e5;
    double f1, f5;
    double g1, g5;
    double h1, h5;
    double denom;

    // REAL*8 (explicit)
    double sum, answer;

    //   dRay(1)=a1Double, dRay(2)=a2Double, dRay(3)=a3Double, dRay(4)=a4Double,
    //   dRay(5)=cDouble,  dRay(7)=zDouble
    //   dRay(6) and dRay(8) are not named but used directly
    double dRay[9]; // 1-based: dRay[1]..dRay[8]
    double& a1Double = dRay[1];
    double& a2Double = dRay[2];
    double& a3Double = dRay[3];
    double& a4Double = dRay[4];
    double& cDouble  = dRay[5];
    double& zDouble  = dRay[7];
    // dRay[6] and dRay[8] are used directly (not named)

    // ---- Begin translated code ----


// 100  sum=1
    sum = 1;

//      a2 = ((X+Y-C) >> 1)
    a2 = ((x+y-c) >> 1);
    a3 = ((a+y-z) >> 1);
    a4 = ((b+x-z) >> 1);
    a1 = ((a+b-c) >> 1);

    if ((((a1) >> (30)) & 1) || (((a2) >> (30)) & 1) ||
        (((a3) >> (30)) & 1) || (((a4) >> (30)) & 1)) return 0.0;

//      b1 = Y - a2 - a3
    b1 = y - a2 - a3;
//      b2 = a3 - a1 + Z - Y
    b2 = a3 - a1 + z - y;

//      r = a1 + a2 + C + 1
    r = a1 + a2 + c + 1;

//      kMax = MIN0 (a1, a2, a3, a4)
    kMax = std::min({a1, a2, a3, a4});
//      kMin = MAX0 (0, -b1, -b2)
    kMin = std::max({0, -b1, -b2});

//      a1Double = a1
    a1Double = a1;
    a2Double = a2;
    a3Double = a3;
    a4Double = a4;
    yDouble = y;
    zDouble = z;
    cDouble = c;
    b1Double = yDouble - a2Double - a3Double;
    b2Double = a3Double - a1Double + zDouble - yDouble;
    rDouble = a1Double + a2Double + cDouble + 1;
    kMinDouble = kMin;
    kMaxDouble = kMax;

    if (kMax - kMin < 0) return 0.0;

    if (kMax - kMin > 0) {
        d1 = a1Double - (kMaxDouble-1);
        d2 = a2Double - (kMaxDouble-1);
        d3 = a3Double - (kMaxDouble-1);
        d4 = a4Double - (kMaxDouble-1);
        d5 = -kMaxDouble;
        d6 = -b1Double - kMaxDouble;
        d7 = -b2Double - kMaxDouble;
        d8 = rDouble - (kMaxDouble-1);

        e1 = d1*d2*d3*d4;
        e5 = d5*d6*d7*d8;

        kMaxM2 = kMax-2;
        sum = 1 + e1/e5;

        if (kMaxM2 >= kMin) {
            f1 = (d1+1)*(d2+1)*(d3+1)*(d4+1) - e1;
            f5 = (d5+1)*(d6+1)*(d7+1)*(d8+1) - e5;
            g1 = (d1+2)*(d2+2)*(d3+2)*(d4+2) - f1 - f1 - e1;
            g5 = (d5+2)*(d6+2)*(d7+2)*(d8+2) - f5 - f5 - e5;
            h1 = (d1+3)*(d2+3)*(d3+3)*(d4+3) - 3*(g1+f1) - e1;
            h5 = (d5+3)*(d6+3)*(d7+3)*(d8+3) - 3*(g5+f5) - e5;

            for (n = kMin; n <= kMaxM2; n++) {
                e1 = e1 + f1;
                e5 = e5 + f5;
                f1 = f1 + g1;
                f5 = f5 + g5;
                g1 = g1 + h1;
                g5 = g5 + h5;
                h1 = h1 + 24;
                h5 = h5 + 24;
                sum = 1 + sum*(e1/e5);
            }
        }
    }

    // Compute answer: factorial table / log-factorial / DLGAMA
    if (r <= factorialTable().maxFactorial) {

//      E = factTable(1+r-kMin) / (factTable(1+kMin) * ...)
    e = factorialTable().factTable[1+r-kMin] / (factorialTable().factTable[1+kMin]
      * factorialTable().factTable[1+b1+kMin] * factorialTable().factTable[1+b2+kMin]
      * factorialTable().factTable[1+a1-kMin]
      * factorialTable().factTable[1+a2-kMin] * factorialTable().factTable[1+a3-kMin]
      * factorialTable().factTable[1+a4-kMin]);

//      sum = std::round(E*sum*E)
    sum = std::round(e*sum*e);
        denom = (factorialTable().factTable[2+a1+c]/(factorialTable().factTable[1+a1]*factorialTable().factTable[1+a4+b1]
                                        *factorialTable().factTable[1+a3+b2]))
          * (factorialTable().factTable[2+a2+c]/(factorialTable().factTable[1+a2]*factorialTable().factTable[1+a3+b1]
                                        *factorialTable().factTable[1+a4+b2]))
          * (factorialTable().factTable[2+a3+z]/(factorialTable().factTable[1+a3]*factorialTable().factTable[1+a2+b1]
                                        *factorialTable().factTable[1+a1+b2]))
          * (factorialTable().factTable[2+a4+z]/(factorialTable().factTable[1+a4]*factorialTable().factTable[1+a1+b1]
                                        *factorialTable().factTable[1+a2+b2]));
        answer = sum/(denom);
    } else if (r <= logFactorialTable().maxLf) {
        // Medium J: log-factorial table
        e = logFactorialTable().lf[1+r-kMin] - logFactorialTable().lf[1+kMin]
          - logFactorialTable().lf[1+b1+kMin] - logFactorialTable().lf[1+b2+kMin]
          - logFactorialTable().lf[1+a1-kMin]
          - logFactorialTable().lf[1+a2-kMin] - logFactorialTable().lf[1+a3-kMin]
          - logFactorialTable().lf[1+a4-kMin];
        e2 = .5 * ( logFactorialTable().lf[1+a1] + logFactorialTable().lf[1+a4+b1]
          + logFactorialTable().lf[1+a3+b2]
            - logFactorialTable().lf[2+a1+c]
          + logFactorialTable().lf[1+a2] + logFactorialTable().lf[1+a3+b1]
          + logFactorialTable().lf[1+a4+b2]
            - logFactorialTable().lf[2+a2+c]
          + logFactorialTable().lf[1+a3] + logFactorialTable().lf[1+a2+b1]
          + logFactorialTable().lf[1+a1+b2]
            - logFactorialTable().lf[2+a3+z]
          + logFactorialTable().lf[1+a4] + logFactorialTable().lf[1+a1+b1]
          + logFactorialTable().lf[1+a2+b2]
            - logFactorialTable().lf[2+a4+z] );
        answer = std::exp(e+e2) * sum;
    } else {
        // Large J: dLogGamma
        dRay[6] = cDouble;
        dRay[8] = zDouble;
        e2 = 0;
        e = dLogGamma(1-kMinDouble+rDouble) - dLogGamma(1+kMinDouble)
          - dLogGamma(1+kMinDouble+b1Double) - dLogGamma(1+kMinDouble+b2Double);
        for (i = 1; i <= 4; i++) {
            e2 = e2 + dLogGamma(1+dRay[i]) + dLogGamma(1+b1Double+dRay[i])
              + dLogGamma(1+b2Double+dRay[i]) - dLogGamma(2+dRay[i]+dRay[4+i]);
            e = e - dLogGamma(1-kMinDouble+dRay[i]);
        }
        e = e + .5*e2;
        answer = std::exp(e) * sum;
    }

    temp = kMin + r;
    if (!(((temp) >> (0)) & 1)) answer = -answer;

    if (isRacah) {
        if (!(((r) >> (0)) & 1)) answer = -answer;
    }

    return answer;
}

// ---- Public entry points ----

// sixJ(A,B,C,X,Y,Z)
double sixJ(int a, int b, int c, int x, int y, int z)
{
    return sixjImpl(a, b, c, x, y, z, false);
}

// racah(A,B,Y,X,C,Z)
// Fortran: ENTRY RACAH(A,B,Y,X,C,Z)
// The ENTRY statement reuses the FUNCTION's parameter names.
// So RACAH's 1st arg maps to A, 2nd to B, 3rd to Y, 4th to X, 5th to C, 6th to Z.
// The caller calls RACAH(j1, j2, j5, j4, j3, j6) to compute W(j1,j2,j5,j4;j3,j6).
double racah(int a, int b, int y, int x, int c, int z)
{
    return sixjImpl(a, b, c, x, y, z, true);
}

// ============================================================================
// Begin WIG9J, coulombSigmaL, getDate, gaussL translations
// ============================================================================


// ============================================================================
// FUNCTION WIG9J — Wigner 9-J coefficient (lines 1070-1518)
//
// COMPUTES 9-J COEFFICIENTS.
// Input: 9 integers that are TWICE the corresponding J's.
// Returns double precision 9-J coefficient.
//
// S. PIEPER
// 9/6/74 - FIRST VERSION
// 4/11/88 - VAX/VMS VERSION - R.OSBORN AND G.L.GOODMAN
// 7/3/91 - RS6000 - s.p.
// ============================================================================
double wig9J(int j1, int j2, int j3, int j4, int j5, int j6,
             int j7, int j8, int j9)
{
    // IMPLICIT INTEGER*4 (A-C, i-R, T-Z), REAL*8 (D-H, S)
    // Vars A-C, i-R, T-Z are int; D-H, S are double.

    auto dlf = [](int n) -> double { return dLogGamma(static_cast<double>(n)); };

    // zeroOne[1]=0.0, zeroOne[2]=1.0 (1-based; index 0 unused)
    double zeroOne[3] = {0, 0.0, 1.0};

    // Local integer variables
    int j1S[6];   // DIMENSION j1S(5), 1-based
    int j2S[6];   // DIMENSION j2S(5)
    int j3S[6];   // DIMENSION j3S(5)
    int as[5][4]; // DIMENSION AS(4,3) -> AS[5][4], 1-based
    int bs[3][4]; // DIMENSION BS(2,3) -> BS[3][4], 1-based
    int mx2s[4];  // DIMENSION MX2S(3), 1-based
    int mn2s[4];  // DIMENSION MN2S(3), 1-based
    int mn1s[4];  // DIMENSION MN1S(3), 1-based
    int rs[4];    // DIMENSION RS(3), 1-based
    int ids[3][4]; // DIMENSION IDS(2,3) -> IDS[3][4], 1-based

    // Local double variables
    double asDouble[5][4]; // DIMENSION asDouble(4,3) -> asDouble[5][4], 1-based
    double bsDouble[3][4]; // DIMENSION bsDouble(2,3) -> bsDouble[3][4], 1-based
    double rsDouble[4];    // DIMENSION DRS(3), 1-based
    double mn1sDouble[4];  // DIMENSION DMN1S(3), 1-based
    double mn2sDouble[4];  // DIMENSION DMN2S(3), 1-based
    double mx2sDouble[4];  // DIMENSION DMX2S(3), 1-based
    double sums[4];   // DIMENSION SUMS(3), 1-based

    // Scalars — integer (A-C, i-R, T-Z)
    int i, j, n, x;
    int xStart, xEnd, xSum, xHalf, xPiece;
    int mx2, mn2, kMax, kMin, kMaxM2;
    int signExponent;

    // Scalars — double (D-H, S)
    double xDouble, xHalfDouble;
    double sumOut, sumLog;
    double sum, sumP;
    double d1, d2, d3, d4, d5, d6, d7, d8;
    double e1, e5;
    double f1, f5;
    double g1, g5;
    double h1, h5;
    double kMaxDouble, kMinDouble;

    double wig9jResult;

    //
    // PICK UP THE J'S IN THE APPROPRIATE ORDER
    // NOTE THAT THE JIS ARRAYS REALLY CONTAIN 2*J
    //
    j1S[1] = j1;
    j1S[2] = j6;
    j1S[3] = j8;
    j2S[1] = j2;
    j2S[2] = j4;
    j2S[3] = j9;
    j3S[1] = j3;
    j3S[2] = j5;
    j3S[3] = j7;
    //
    // EXPAND THE ARRAYS TO ALLOW CYCLIC INDICES
    //
    for (i = 1; i <= 2; i++) {
        j1S[i+3] = j1S[i];
        j2S[i+3] = j2S[i];
        j3S[i+3] = j3S[i];
    } // 59

    //
    // DETERMINE THE RANGE OF 2*X
    //
    xStart = std::abs(j2S[1] - j1S[2]);
    xEnd = j2S[1] + j1S[2];
    xSum = xEnd;
    for (i = 2; i <= 3; i++) {
        xStart = std::max(xStart, std::abs(j2S[i]-j1S[i+1]));
        xPiece = j2S[i] + j1S[i+1];
        xEnd = std::min(xEnd, xPiece);
        xSum = xSum + xPiece;
    } // 129

    if (xEnd < xStart) return 0.0;

    xSum = xSum - xEnd;

    //
    // EXTRACT 1/2 INTEGER PART OF X AND CONVERT 2*X TO INTPART(X)
    //
    xHalf = 0;
    if ((((xStart) >> (0)) & 1)) xHalf = 1;
    xStart = ((xStart) >> 1);
    xEnd = ((xEnd) >> 1);
    xSum = ((xSum) >> 1);
    xHalfDouble = zeroOne[xHalf+1];
    xDouble = xStart - 1;
    xSum = xSum + 1;

    //
    // SETUP THE J1+J2-J3, ETC
    //
    for (i = 1; i <= 3; i++) {
        //
        // THE AS, FOR EACH 6-J THEY ARE
        //   a1 = J1+J2-J3
        //   A2 = J4+J5-J3
        //   A3 = J1+J5-J6
        //   A4 = J4+J2-J6
        // IN A3 AND A4 THE J6 = X IS LEFT OUT.
        //
        as[1][i] = ((j1S[i]+j2S[i]-j3S[i]) >> 1);
        as[2][i] = ((j1S[i+1]+j2S[i+2]-j3S[i]) >> 1);
        as[3][i] = ((j1S[i]+j2S[i+2]) >> 1);
        as[4][i] = ((j2S[i]+j1S[i+1]) >> 1);
        //
        // NOW MAKE SURE NONE OF THEM ARE NEGATIVE WHICH WOULD INDICATE
        // A VIOLATION OF THE TRIANGLE RULES.
        //
        if ((((as[1][i]) >> (30)) & 1) || (((as[2][i]) >> (30)) & 1)) return 0.0;
        //
        // THE BS, FOR EACH 6-J THEY ARE
        //    B1 = J3 + J6 - J1 - J4
        //    B2 = J3 + J6 - J2 - J5
        // THE J6 = X IS LEFT OUT
        //
        bs[1][i] = j2S[i+2] - as[2][i] - as[3][i];
        bs[2][i] = j1S[i] - as[1][i] - as[3][i];
    } // 159

    //
    // SETUP THE QUANTITIES FOR EACH OF THE 3 6-J SUMS.
    //
    for (i = 1; i <= 3; i++) {
        mn1s[i] = std::min(as[1][i], as[2][i]);
        mn2s[i] = std::min(as[3][i], as[4][i]);
        mx2s[i] = -std::min(bs[1][i], bs[2][i]);
        rs[i] = 1 + xHalf + as[3][i] + as[4][i];
        for (j = 1; j <= 4; j++) {
            asDouble[j][i] = as[j][i];
        } // 239
        bsDouble[1][i] = bs[1][i];
        bsDouble[2][i] = bs[2][i];
        rsDouble[i] = 1 + xHalfDouble + asDouble[3][i] + asDouble[4][i];
        mn1sDouble[i] = mn1s[i];
        mn2sDouble[i] = mn2s[i];
        mx2sDouble[i] = mx2s[i];
        //
        // THESE ARE INDICES FOR THE 4 FACTORIALS THAT DEPEND ON X
        //
        ids[1][i] = as[2][i] + bs[2][i];
        ids[2][i] = as[1][i] + bs[1][i];
    } // 299

    //
    // FIND THE LOG OF THE OUTSIDE FACTORS NOW
    //
    sumOut = 0;

    //
    // xSum IS >= ANY OF THE REQUIRED FACTORIAL ARGUMENTS
    //
    if (xSum <= logFactorialTable().maxLf) {
        // WE CAN USE THE LOG(FACTORIAL) TABLE
        for (i = 1; i <= 3; i++) {
            for (j = 1; j <= 2; j++) {
                sumOut = sumOut +
                    logFactorialTable().lf[1+as[j][i]] - logFactorialTable().lf[2+j3S[i]+as[j][i]]
                    + logFactorialTable().lf[1-as[j][i]+j1S[i-1+j]]
                    + logFactorialTable().lf[1-as[j][i]+j2S[i-2+2*j]];
            }
        }
    } else {
        // WE MUST USE LOG GAMMA
        for (i = 1; i <= 3; i++) {
            for (j = 1; j <= 2; j++) {
                sumOut = sumOut +
                    dLogGamma(1+asDouble[j][i]) - dlf(2+j3S[i]+as[j][i])
                    + dlf(1-as[j][i]+j1S[i-1+j])
                    + dlf(1-as[j][i]+j2S[i-2+2*j]);
            }
        }
    }

    sumOut = .5*sumOut;

    //
    // WE ARE NOW READY TO DO THE 9-J LOOP
    //
    wig9jResult = 0;

    for (x = xStart; x <= xEnd; x++) {
        xDouble = xDouble + 1;
        //
        // DO THE PARTS OF EACH 6-J THAT DEPEND ON X
        //
        sumLog = sumOut;
        signExponent = 0;

        for (i = 1; i <= 3; i++) {
            //
            // GET THE COEFFICIENTS FOR THIS 6-J
            //
            mx2 = mx2s[i] - x;
            mn2 = mn2s[i] - x;
            kMax = mn1s[i];
            kMaxDouble = mn1sDouble[i];
            kMinDouble = 0;
            kMin = 0;
            //
            // GET THE RANGE OF THE 6-J sum
            //
            if (mn2 < kMax) {
                kMax = mn2;
                kMaxDouble = mn2sDouble[i] - xDouble;
            }

            if (mx2 > 0) {
                kMin = mx2;
                kMinDouble = mx2sDouble[i] - xDouble;
            }

            sum = 1;
            signExponent = signExponent + kMin;

            //
            // THIS ONE TEST CHECKS ALL 12 INEQUALITIES IMPLIED BY THE
            //
            if (kMax - kMin < 0) return 0.0;
            if (kMax - kMin == 0) { /* skip sum, go to label_700 */ }
            else {
            // else (kMax - kMin > 0) fall through to label_600

            // label_600:
            d1 = asDouble[1][i] - (kMaxDouble-1);
            d2 = asDouble[2][i] - (kMaxDouble-1);
            d8 = rsDouble[i] - (kMaxDouble-1);
            d3 = asDouble[3][i] - (kMaxDouble-1+xDouble);
            d4 = asDouble[4][i] - (kMaxDouble-1+xDouble);
            d5 = -kMaxDouble;
            d6 = -bsDouble[1][i] - (xDouble + kMaxDouble);
            d7 = -bsDouble[2][i] - (xDouble + kMaxDouble);

            e1 = d1*d2*d3*d4;
            e5 = d5*d6*d7*d8;

            kMaxM2 = kMax-2;
            sum = 1 + e1/e5;
            if (kMaxM2 >= kMin) {

            f1 = (d1+1)*(d2+1)*(d3+1)*(d4+1) - e1;
            f5 = (d5+1)*(d6+1)*(d7+1)*(d8+1) - e5;

            g1 = (d1+2)*(d2+2)*(d3+2)*(d4+2) - f1 - f1 - e1;
            g5 = (d5+2)*(d6+2)*(d7+2)*(d8+2) - f5 - f5 - e5;

            h1 = (d1+3)*(d2+3)*(d3+3)*(d4+3) - 3*(g1+f1) - e1;
            h5 = (d5+3)*(d6+3)*(d7+3)*(d8+3) - 3*(g5+f5) - e5;

            //
            // COMPUTE THE sum IN THE RACAH FORMULA
            //
            for (n = kMin; n <= kMaxM2; n++) {
                e1 = e1 + f1;
                e5 = e5 + f5;
                f1 = f1 + g1;
                f5 = f5 + g5;
                g1 = g1 + h1;
                g5 = g5 + h5;
                h1 = h1 + 24;
                h5 = h5 + 24;
                sum = 1 + sum*(e1/e5);
            } // 659
            } // kMaxM2 >= kMin
            } // kMax - kMin > 0

            sums[i] = sum;

            //
            // HERE ARE THE PARTS OF THE FACTORIALS FOR kMin AND ALSO THOSE
            // FACTORIALS THAT DEPEND ON X
            //
            if (rs[i] <= logFactorialTable().maxLf) {
                // WE CAN USE THE LOG FACTORIAL TABLE
                sumLog = sumLog +
                    logFactorialTable().lf[1-x+as[4][i]] + logFactorialTable().lf[1+x+ids[1][i]]
                    + logFactorialTable().lf[1+x+ids[2][i]] - logFactorialTable().lf[2+xHalf+x+as[4][i]]
                    + logFactorialTable().lf[1-kMin+rs[i]] - logFactorialTable().lf[1+kMin]
                    - logFactorialTable().lf[1+kMin+x+bs[1][i]]
                    - logFactorialTable().lf[1+kMin+x+bs[2][i]]
                    - logFactorialTable().lf[1-kMin+as[1][i]] - logFactorialTable().lf[1-kMin+as[2][i]]
                    - logFactorialTable().lf[1-kMin-x+as[3][i]]
                    - logFactorialTable().lf[1-kMin-x+as[4][i]];
            } else {
                // WE MUST USE LOG GAMMA
                sumLog = sumLog +
                    dLogGamma(1-xDouble+asDouble[4][i]) + dLogGamma(1+xDouble+asDouble[2][i]+bsDouble[2][i])
                    + dLogGamma(1+xDouble+asDouble[1][i]+bsDouble[1][i])
                    - dLogGamma(2+xHalfDouble+xDouble+asDouble[4][i])
                    + dLogGamma(1-kMinDouble+rsDouble[i]) - dLogGamma(1+kMinDouble)
                    - dLogGamma(1+kMinDouble+xDouble+bsDouble[1][i])
                    - dLogGamma(1+kMinDouble+xDouble+bsDouble[2][i])
                    - dLogGamma(1-kMinDouble+asDouble[1][i]) - dLogGamma(1-kMinDouble+asDouble[2][i])
                    - dLogGamma(1-kMinDouble-xDouble+asDouble[3][i])
                    - dLogGamma(1-kMinDouble-xDouble+asDouble[4][i]);
            }
        } // DO 759 i = 1, 3

        //
        // NOW MULTIPLY THE 3 6-J PARTS TOGETHER AND sum
        //
        sumP = sums[1] * sums[2] * std::exp(sumLog) * sums[3]
            * (1+xHalfDouble+(xDouble+xDouble));

        if ((((signExponent) >> (0)) & 1)) sumP = -sumP;
        wig9jResult = wig9jResult + sumP;
    } // 799 DO X = xStart, xEnd

    return wig9jResult;
}

// ============================================================================
//
// COULOMB PHASES
// ?/?/? - ORIGINAL VERSION
// 12/28/79 - CDC TO CNI, TITLE COMMENT - RPG
// 7/9/80 - USE CDM PREFIX - S.P.
// 7/3/91 - RS6000 version - s.p.
// ============================================================================
void coulombSigmaL(double eta, int lMax, double* dsg)
{
    // IMPLICIT REAL*8 ( A-H, O-Z )
    // i-N are int, rest are double

    // Local integer variables
    int i, j, k, l, m, n;

    // Local double variables
    double x, xSq, rSq, qSum, pSum, r;

    // DIMENSION P(26), Q(23) — 1-based
    static const double p[27] = {0,
        0.3101657810129948870e-08,
        0.3015947478999109100e-05,
        0.2056442871538789580e-03,
        0.4445588614720562960e-02,
        0.4241372519181223210e-01,
        0.2063788681972964780e+00,
        0.5468922695201207380e+00,
        0.7949480657792201980e+00,
        0.5931315608708119800e+00,
        0.1770600075369600650e+00,
        0.2410878569735943570e-08,
        0.9326569949959554900e-06,
        0.7585813793142706060e-04,
        0.2007067188606569890e-02,
        0.1889244677027970660e-01,
        0.4443159218008918380e-01,
       -0.1596016183850622740e+00,
       -0.6851417732759696180e+00,
       -0.5772019207036128280e+00,
        0.1395410517881128990e+03,
        0.1090522693580003650e+04,
        0.4853585241219512340e+03,
       -0.6587580539038850700e+03,
       -0.4130231383850326670e+02,
        0.9229335782342382280e+01,
       -0.9999999999999999440e+00
    };

    static const double q[24] = {0,
        0.1626320933946765800e-05,
        0.1808585994795038710e-03,
        0.5771439575139208560e-02,
        0.7876590720775496210e-01,
        0.5451579810646677990e+00,
        0.2085536101102752230e+01,
        0.4573132992231460970e+01,
        0.5697178591946184810e+01,
        0.3737311340573166380e+01,
        0.8418712796550826150e-09,
        0.5251990683341128460e-06,
        0.6387938953567890980e-04,
        0.2628228899692266520e-02,
        0.4482759373310784680e-01,
        0.3453892205334134070e+00,
        0.1223260295405861340e+01,
        0.1881032630842315980e+01,
       -0.1411156525262207530e+03,
       -0.9965362519625687700e+03,
       -0.5358167734466346130e+03,
        0.6544271625884707650e+03,
        0.4215892515370163450e+02,
       -0.9312669115675720870e+01
    };

    static const double ez  = 0.1805547071605106970e+01;
    static const double ez1 = 0.1805547714233398440e+01;
    static const double ez2 = 0.6426282915176236330e-06;

    x = std::fabs(eta);
    xSq = x*x;

    if (x <= 2.0) {
        // ABS(ETA) IN (0,2) RANGE
        i = 7; j = 1; k = 1; m = 1;
        rSq = xSq;
    } else if (x <= 4.0) {
        // ABS(ETA) IN (2,4) RANGE
        i = 6; j = 11; k = 10; m = 2;
        rSq = xSq;
    } else {
        // ABS(ETA) IN (4,inf) RANGE
        i = 4; j = 20; k = 18; m = 3;
        rSq = 1.0e0/xSq;
    }

    qSum = q[k];
    n = k+i;
    for (l = k; l <= n; l++) {
        qSum = qSum*rSq + q[l+1];
    }
    pSum = p[j];
    n = j+i+1;
    for (l = j; l <= n; l++) {
        pSum = pSum*rSq + p[l+1];
    }
    r = pSum/(qSum*rSq + 1.0e0);

    // dsg is 0-based: dsg[L] = sigma_L for L = 0..lMax.
    if (m == 1) {
        dsg[0] = x*r*(x+ez)*((x-ez1)+ez2);
    } else if (m == 2) {
        dsg[0] = x*r;
    } else {
        dsg[0] = std::atan(x)/2.0e0 + x*(std::log(1.0e0+xSq)/2.0e0 + r);
    }

    // SL(ETA)=SO(ETA)+sum(i=1TOL)std::atan(ETA/i)
    if (lMax != 0) {
        r = 1.0e0;
        for (i = 1; i <= lMax; i++) {
            rSq = x/r;
            dsg[i] = dsg[i-1] + std::atan(rSq);
            r = r + 1.0e0;
        }
        if (eta < 0.0e0) {
            for (i = 1; i <= lMax; i++) {
                dsg[i] = -dsg[i];
            }
        }
    }

    if (eta < 0.0e0) dsg[0] = -dsg[0];

    return;
}

// ============================================================================
//
// RETURNS THE DATE IN THE FORM  19 Jan 01
// date - a character*9
//
// 11/16/01 - all new
// ============================================================================
void getDate(char* date)
{
    static const char* months[13] = {"",
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };

    // Use C++ time functions in place of Fortran date_and_time
    std::time_t now = std::time(nullptr);
    struct std::tm* tmInfo = std::localtime(&now);

    int day   = tmInfo->tm_mday;        // 1-31
    int month = tmInfo->tm_mon + 1;     // 1-12
    int year  = tmInfo->tm_year % 100;  // 2-digit year

    // Format: "DD Mon YY" (9 characters)
    // date(1:2) = day, date(4:6) = month name, date(8:9) = year
    std::memset(date, ' ', 9);
    date[9] = '\0';

    // Write day (positions 0-1, i.e. Fortran 1:2)
    date[0] = (day / 10 == 0) ? ' ' : ('0' + day / 10);
    date[1] = '0' + day % 10;

    // date(3) = ' '  (space separator, already set)

    // Write month name (positions 3-5, i.e. Fortran 4:6)
    date[3] = months[month][0];
    date[4] = months[month][1];
    date[5] = months[month][2];

    // date(7) = ' '  (space separator, already set)

    // Write year (positions 7-8, i.e. Fortran 8:9)
    date[7] = '0' + year / 10;
    date[8] = '0' + year % 10;

    return;
}

// ============================================================================
//
// POINTS AND WEIGHTS FOR GAUSS-LEGENDRE INTEGRATION OVER (-1 1)
// COMPUTED BY METHOD GIVEN ON PP.88 AND 89 OF
// 'METHODS OF NUMERICAL INTEGRATION' DAVIS AND RABINOWITZ
// ACADEMIC PRESS (1975)
//
// P.369 (APP 2) OF DAVIS AND RABINOWITZ
// MHM  JUNE 29, 1975
// 7/5/91 - RS6000 - s.p.
// ============================================================================
void gaussL(int n, double* x, double* w)
{
    int i, k, m;
    double x0, e1, e2, e3, t, t1;
    double pk, pkM1, pkP1, fk, den;
    double d1, dPn, d2Pn, d3Pn, d4Pn;
    double u, v, h, h1, h2;
    double p, dp, fx;

    if (n <= 0) {
        std::printf("\n  0000 NUMBER OF GAUSS POINTS = %5d 0000\n"
                    "   0000             ABSURD             0000\n", n);
        return;
    }

    if (n == 1) {
        x[1] = 0;
        w[1] = 2;
        return;
    }

    // BEGIN COMPUTATION
    m = (n+1)/2;
    e1 = n*(n+1);
    e2 = 3.141592653589793200e0 / (4*n+2);
    e3 = 1.-(1.-1.0e0/n)/(8.0e0*n*n);
    for (i = 1; i <= m; i++) {
        t = (4*i-1)*e2;
        x0 = e3*std::cos(t);
        pkM1 = 1;
        pk = x0;
        fk = 1;
        for (k = 2; k <= n; k++) {
            fk = fk + 1;
            t1 = x0*pk;
            pkP1 = t1 - pkM1 - (t1-pkM1)/fk + t1;
            pkM1 = pk;
            pk = pkP1;
        }
        den = 1. - x0*x0;
        d1 = n*(pkM1-x0*pk);
        dPn = d1/den;
        d2Pn = (2.*x0*dPn - e1*pk)/den;
        d3Pn = (4.*x0*d2Pn + (2.-e1)*dPn)/den;
        d4Pn = (6.*x0*d3Pn + (6.-e1)*d2Pn)/den;
        u = pk/dPn;
        v = d2Pn/dPn;
        h1 = (1.+.5*u*(v+u*(v*v-u*d3Pn/(3.*dPn))));
        h2 = -u*h1;
        h = h2;
        p = pk + h*(dPn + .5*h*(d2Pn + h/3.*(d3Pn + .25*h*d4Pn)));
        dp = dPn + h*(d2Pn + .5*h*(d3Pn + h*d4Pn/3.));
        h = h - p/dp;
        x[n+1-i] = x0 + h;
        x[i] = -x[n+1-i];
        fx = d1 - h*e1*(pk + .5*h*(dPn + h/3.*(d2Pn + .25*h*(d3Pn + .2*h*d4Pn))));
        w[n+1-i] = 2.*(1-x[i]*x[i])/(fx*fx);
        w[i] = w[n+1-i];
    }
    if (m+m > n) x[m] = 0;

    return;
}
