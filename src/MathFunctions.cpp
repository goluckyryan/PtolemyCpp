// math_functions.cpp — second half of fortlib.f translation
// cubicSplineInterp, LAGBC, LGRECR, LGroot, LAGUER, lsqPol, MATINV,

#include "math/numeric_utils.h"
#include "math/special.h"
#include "math/spline.h"
#include "math/linear_algebra.h"
#include "math/legendre.h"
#include "math/gauss_quadrature.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include <vector>
#include "Constants.h"

// File-scope scratch for matInverse (pivot/index) and lsqPol (xPower reuses the
// bytes as a 1-based double array). Was global matrixWork in ptolemy_commons;
// no other translation unit ever read or wrote it.
namespace {
struct MatrixWork {
    double pivot[51];  // 1-based [1..50]: pivot elements
    int    index[51];  // 1-based [1..50]: pivot index array
};
MatrixWork matrixWork = {};
}

// matInverse is file-static (only lsqPol calls it) — forward decl for use-before-def.
// determinant out-param dropped — sole caller (lsqPol) declared it, passed by
// ref, never read it. MATINV still tracks determinant internally (singularity
// check + scaling).
static void matInverse(double* a, int n, double* b, int m);

// ==========================================================================
// cubicSplineInterp — cubic spline interpolation (lines 2312-2419)
// ==========================================================================
void cubicSplineInterp(int cubicCount, double* xCubes, double* aCoef, double* bCoef, double* cCoef,
            double* dCoef, int nPts, double* xPts, double* yOut)
{
    double inf = 1.0e300;

    double xSign = std::copysign(1.0, xCubes[2] - xCubes[1]);

    // Fake starting conditions that will cause reset on first point
    double xPrevious = inf;
    double xNext = inf;
    double xBase = 0.0;
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
    int n = 0;

    for (int i = 1; i <= nPts; i++) {
        double x = xPts[i];
        double xc = xSign * x;

        // Find the correct spline segment for this x
        for (;;) {
            if (xc < xNext) {
                // x is below the upper bound of current segment
                if (xc >= xPrevious) {
                    // x is in range [xPrevious, xNext) — evaluate
                    break;
                }
                // x is below xPrevious — reset to start
                n = 1;
                xPrevious = -inf;
                xNext = xCubes[2] * xSign;
                xBase = xCubes[1];
                // Load coefficients and re-check
                a = aCoef[n];
                b = bCoef[n];
                c = cCoef[n];
                d = dCoef[n];
                continue;
            }

            // x >= xNext — advance to next spline segment
            xPrevious = xNext;
            xBase = xNext * xSign;
            n = n + 1;
            xNext = xCubes[n + 1] * xSign;
            if (n == cubicCount - 1) xNext = inf;

            a = aCoef[n];
            b = bCoef[n];
            c = cCoef[n];
            d = dCoef[n];
            // Re-check if x is inside this new spline
        }

        // Evaluate the spline
        double delta = x - xBase;
        yOut[i] = a + delta * (b + delta * (c + delta * d));
    }

    return;
}

// ==========================================================================
// LAGBC — Laguerre recurrence coefficients (lines 2564-2580)
// ==========================================================================
static void laguerreRecurCoeffs(int nn, double alpha, double* a, double* b, double* c)
{
    for (int n = 1; n <= nn; n++) {
        double en = n;
        a[n] = 1.0 / en;
        b[n] = alpha + en + en - 1.0;
        c[n] = (alpha + en - 1.0) / en;
    }
    return;
}

// ==========================================================================
// LGRECR — Laguerre polynomial recurrence (lines 2539-2563)
// ==========================================================================
static void laguerreRecur(double& pN, double& dPn, double& pN1, double x,
            int nn, double alpha, double* a, double* b, double* c)
{
    double pPrevious = 1.0;
    double pCurr = alpha + 1.0 - x;
    double pPreviousPrime = 0.0;
    double pCurrPrime = -1.0;

    for (int j = 2; j <= nn; j++) {
        double pNext = a[j] * (b[j] - x) * pCurr - c[j] * pPrevious;
        double pNextPrime = a[j] * (b[j] - x) * pCurrPrime - c[j] * pPreviousPrime - a[j] * pCurr;
        pPrevious = pCurr;
        pCurr = pNext;
        pPreviousPrime = pCurrPrime;
        pCurrPrime = pNextPrime;
    }

    pN = pCurr;
    dPn = pCurrPrime;
    pN1 = pPrevious;
    return;
}

// ==========================================================================
// LGroot — Newton iteration for Laguerre root (lines 2515-2538)
// ==========================================================================
static void laguerreRoot(double& x, int nn, double alpha, double& dPn, double& pN1,
            double* a, double* b, double* c)
{
    constexpr double eps = 1.0e-20;
    double polyValue, polyDeriv, newtonStep;

    for (int iter = 1; iter <= 10; iter++) {
        laguerreRecur(polyValue, polyDeriv, pN1, x, nn, alpha, a, b, c);
        newtonStep = polyValue / polyDeriv;
        x = x - newtonStep;
        if (std::fabs(newtonStep / x) < eps) break;
    }

    dPn = polyDeriv;
    return;
}

// ==========================================================================
// LAGUER — Gauss-Laguerre quadrature (lines 2421-2514)
// ==========================================================================
// 32-element 1-based scratch arrays only to hand them through, and never
// read them back. Allocated as locals here now.
void laguerre(int nn, double* x, double* w, double alpha)
{
    // CSX/CSW/TSX/TSW out-params dropped — sole caller (CWF coulombIntegral)
    // declared them, passed them in, and never read them back. CSX/CSW
    // accumulated the sum of x and w; TSX/TSW computed FN*(FN+alpha) and CC
    // (look like vestigial sanity-check anchors). All 4 assignments removed.
    double fn = nn;
    // xTrial init to 0 silences -Wmaybe-uninitialized for the (unreachable in
    // practice) nn==0 early-loop-exit edge — first iteration (i=1, i-2<0)
    // unconditionally assigns xTrial before any read.
    double xTrial = 0;
    double r1, r2, ratio, fi, dPn, pN1;

    // 1-based scratch; nn comes in <= 28 from the sole caller's nPts=28.
    double a[32], b[32], c[32];
    laguerreRecurCoeffs(nn, alpha, a, b, c);

    double weightNorm = dGamma(alpha + 1.0);

    for (int j = 2; j <= nn; j++) {
        weightNorm = weightNorm * c[j];
    }

    for (int i = 1; i <= nn; i++) {
        // Arithmetic IF: IF (i - 2) 2, 4, 5
        if (i - 2 < 0) {
            // label_2: smallest zero
            xTrial = (1.0 + alpha) * (3.0 + 0.920 * alpha) /
                 (1.0 + 2.40 * fn + 1.80 * alpha);
        } else if (i - 2 == 0) {
            // label_4: second zero
            xTrial = xTrial + (15.0 + 6.250 * alpha) /
                 (1.0 + 0.90 * alpha + 2.50 * fn);
        } else {
            // label_5: all other zeros
            fi = i - 2;
            r1 = (1.0 + 2.550 * fi) / (1.90 * fi);
            r2 = 1.260 * fi * alpha / (1.0 + 3.50 * fi);
            ratio = (r1 + r2) / (1.0 + 0.30 * alpha);
            xTrial = xTrial + ratio * (xTrial - x[i - 2]);
        }

        // label_6:
        laguerreRoot(xTrial, nn, alpha, dPn, pN1, a, b, c);
        x[i] = xTrial;
        w[i] = -weightNorm / dPn / pN1;
    }

    return;
}

// ==========================================================================
// lsqPol — Least squares polynomial fit (lines 2712-2824)
// ==========================================================================
// sum out-param dropped — sole caller (grid_setup poly-fit stage) never
// read it. Was computing sum of W*residual^2 per polynomial as a vestigial
// goodness-of-fit anchor; the residual array carries everything
// the caller needs (it adds residuals back into the work block).
void lsqPol(double* X, double* Y, double* W, double* residual, int nSub,
            int lSub, double* A, double* B, int mSub)
{
    // xPower overlays matrixWork block (100 doubles of scratch space)
    // matrixWork is reused as raw double workspace here
    double* xPower = reinterpret_cast<double*>(&matrixWork) - 1; // 1-based access

    double term, poly;
    int k1, k2, i2;

    int n = nSub;
    int l = lSub;
    int m = mSub;
    const int nMax = n;
    const int mMax = m;
    int m1 = m + 1;
    int m3 = m + m + m;
    int m31 = m3 - 1;
    int m41 = m31 + m;

    // Scale X into (-1,1) to prevent overflow or underflow
    double xMax = 0.0;
    for (int k = 1; k <= n; k++) {
        xMax = std::max(xMax, std::fabs(X[k]));
    }
    double xScale = 1.0 / xMax;
    for (int k = 1; k <= n; k++) {
        X[k] = xScale * X[k];
    }

    // Formation and inversion of system of normal equations
    std::fill(xPower + m1, xPower + m41 + 1, 0.0);

    for (k1 = 1; k1 <= n; k1++) {
        term = W[k1 - 1];  // W (vWtsPointer) is 0-based; X/Y stay 1-based
        for (k2 = m1; k2 <= m31; k2++) {
            xPower[k2] = term + xPower[k2];
            term = X[k1] * term;
        }
    }

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= m; j++) {
            k2 = i + j + m - 1;
            // A(i,j) = xPower(k2)
            A[i + (j - 1) * mMax - 1] = xPower[k2];  // A (vAmatPointer) 0-based
        }
    }

    for (int j = 1; j <= l; j++) {
        for (int k = 1; k <= n; k++) {
            // term = W(k)*Y(k,j)
            term = W[k - 1] * Y[k + (j - 1) * nMax];  // W 0-based; Y 1-based
            for (k2 = m3; k2 <= m41; k2++) {
                xPower[k2] = term + xPower[k2];
                term = X[k] * term;
            }
        }
        for (int i = 1; i <= m; i++) {
            k2 = i + m31;
            // B(i,j) = xPower(k2)
            B[i + (j - 1) * mMax - 1] = xPower[k2];
            if (j != l) xPower[k2] = 0.0;
        }
    }

    // Call matInverse
    matInverse(A, m, B, l);

    // Evaluation of residuals
    for (int j = 1; j <= l; j++) {
        for (int k = 1; k <= n; k++) {
            poly = 0.0;
            for (i2 = 1; i2 <= m; i2++) {
                int i = m1 - i2;
                // poly = X(k)*poly + B(i,j)
                poly = X[k] * poly + B[i + (j - 1) * mMax - 1];
            }
            // RESID(k,j) = poly - Y(k,j)
            residual[k - 1 + (j - 1) * nMax] = poly - Y[k + (j - 1) * nMax];  // residual 0-based; Y 1-based
        }
    }

    // Un-scale the coefficients and X
    for (int j = 1; j <= l; j++) {
        for (i2 = 2; i2 <= m; i2++) {
            for (int i = i2; i <= m; i++) {
                // B(i,j) = xScale*B(i,j)
                B[i + (j - 1) * mMax - 1] = xScale * B[i + (j - 1) * mMax - 1];
            }
        }
    }
    for (int k = 1; k <= n; k++) {
        X[k] = xMax * X[k];
    }

    return;
}

// ==========================================================================
// MATINV — Matrix inversion with linear equation solution (lines 2826-2948)
// Gauss-Jordan elimination with full pivot searching
// ==========================================================================
static void matInverse(double* a, int n, double* b, int m)
{
    const int nMax = n;
    double* pivot = matrixWork.pivot; // 1-based, [1..50]
    int*    index = matrixWork.index; // 1-based, [1..50]

    double detMax = 1.0e300;
    double detMin = 1.0e-300;

    int pivotRow = 0, pivotCol = 0;
    double aMax, temp, swap, t;

    // Initialize determinant and pivot element array.
    // determinant is now a local; only used for singularity break + scaling.
    // iDet (sign-flip counter into pivot[1]) dropped — was never read after
    // matrix inversion completed.
    double determinant = 1.0;
    for (int i = 1; i <= n; i++) {
        pivot[i] = 0.0;
        index[i] = 0;
    }

    // Perform successive pivot operations (grand loop)
    for (int i = 1; i <= n; i++) {

        // Search for pivot element
        aMax = 0.0;
        for (int j = 1; j <= n; j++) {
            if (pivot[j] != 0.0) continue;
            for (int k = 1; k <= n; k++) {
                if (pivot[k] != 0.0) continue;
                temp = std::fabs(a[j + (k - 1) * nMax - 1]);
                if (temp < aMax) continue;
                pivotRow = j;
                pivotCol = k;
                aMax = temp;
            }
        }

        index[i] = 4096 * pivotRow + pivotCol;
        int j = pivotRow;
        aMax = a[j + (pivotCol - 1) * nMax - 1];
        determinant = aMax * determinant;

        // Determinant scaling to avoid overflow/underflow
        if (std::fabs(determinant) >= detMax) {
            determinant = determinant * detMin;
        } else if (std::fabs(determinant) <= detMin) {
            determinant = determinant * detMax;
        }

        // Return if matrix is singular (zero pivot)
        if (determinant == 0.0) break;

        pivot[pivotCol] = aMax;

        // Interchange rows to put pivot element on diagonal
        if (pivotRow != pivotCol) {
            determinant = -determinant;
            for (int k = 1; k <= n; k++) {
                swap = a[j + (k - 1) * nMax - 1];
                a[j + (k - 1) * nMax - 1] = a[pivotCol + (k - 1) * nMax - 1];
                a[pivotCol + (k - 1) * nMax - 1] = swap;
            }
            if (m > 0) {
                for (int k = 1; k <= m; k++) {
                    swap = b[j + (k - 1) * nMax - 1];
                    b[j + (k - 1) * nMax - 1] = b[pivotCol + (k - 1) * nMax - 1];
                    b[pivotCol + (k - 1) * nMax - 1] = swap;
                }
            }
        }

        // Divide pivot row by pivot element
        a[pivotCol + (pivotCol - 1) * nMax - 1] = 1.0;
        for (int k = 1; k <= n; k++) {
            a[pivotCol + (k - 1) * nMax - 1] = a[pivotCol + (k - 1) * nMax - 1] / aMax;
        }
        if (m > 0) {
            for (int k = 1; k <= m; k++) {
                b[pivotCol + (k - 1) * nMax - 1] = b[pivotCol + (k - 1) * nMax - 1] / aMax;
            }
        }

        // Reduce non-pivot rows
        for (int j2 = 1; j2 <= n; j2++) {
            if (j2 == pivotCol) continue;
            t = a[j2 + (pivotCol - 1) * nMax - 1];
            a[j2 + (pivotCol - 1) * nMax - 1] = 0.0;
            for (int k = 1; k <= n; k++) {
                a[j2 + (k - 1) * nMax - 1] = a[j2 + (k - 1) * nMax - 1]
                                        - a[pivotCol + (k - 1) * nMax - 1] * t;
            }
            if (m > 0) {
                for (int k = 1; k <= m; k++) {
                    b[j2 + (k - 1) * nMax - 1] = b[j2 + (k - 1) * nMax - 1]
                                            - b[pivotCol + (k - 1) * nMax - 1] * t;
                }
            }
        }
    } // end grand loop on i

    // Interchange columns after all pivot operations
    for (int i = 1; i <= n; i++) {
        int i1 = n + 1 - i;
        int k = index[i1] / 4096;
        pivotCol = index[i1] - 4096 * k;
        if (k != pivotCol) {
            for (int j = 1; j <= n; j++) {
                swap = a[j + (k - 1) * nMax - 1];
                a[j + (k - 1) * nMax - 1] = a[j + (pivotCol - 1) * nMax - 1];
                a[j + (pivotCol - 1) * nMax - 1] = swap;
            }
        }
    }

    // iDet sign-flip count stored to pivot[1] dropped — written to a
    // file-scope matrixWork slot that no caller ever reads back.
    return;
}

// ==========================================================================
// plmSub — Associated Legendre functions P(L,M)(X) (lines 2949-3072)
// ==========================================================================
void plmSub(int lMax, int mMax, double X, double* plmPointer)
{
    double root = 0.0;

    if (mMax != 0) {
        root = std::sqrt(std::fabs(1.0 - X * X));
        if (std::fabs(X) > 1.0) root = -root;
    }

    plmPointer[1] = 1.0;

    int i = 1;
    int j = 1;
    int recurrenceStart, recurrenceEnd;

    for (int m = 0; m <= mMax; m++) {
        if (m >= lMax) return;
        j = i;

        // Generate the P(L=M+1, M) for this M
        i = i + 1;
        plmPointer[i] = (2 * m + 1) * X * plmPointer[i - 1];
        int mP2 = m + 2;
        if (lMax >= mP2) {
            // Now do L = M+2, M+3, ..., lMax for this M
            double mP2Double = mP2;
            double temp1 = (2 * mP2Double - 1) * X;
            double temp2 = m - 1 + mP2Double;
            double temp3 = mP2Double - m;
            recurrenceStart = i + 1;
            recurrenceEnd = i + 1 + lMax - mP2;

            for (i = recurrenceStart; i <= recurrenceEnd; i++) {
                plmPointer[i] = (temp1 * plmPointer[i - 1] - temp2 * plmPointer[i - 2]) / temp3;
                temp1 = temp1 + (2.0 * X);
                temp2 = temp2 + 1.0;
                temp3 = temp3 + 1.0;
            }

            i = recurrenceEnd;
        }

        if (m == mMax) return;

        // Generate the next P(L=M+1, M+1)
        i = i + 1;
        plmPointer[i] = -(2 * m + 1) * root * plmPointer[j];
    }

    return;
}

// wrapper around plmSub(lMax, 0, X, plmPointer) with no callers.

// ==========================================================================
// naturalCubicSpline — Natural cubic spline (lines 4946-5133)
// ==========================================================================
void naturalCubicSpline(int pointCount, double* xPts, double* yIn, double* bCoef, double* cCoef, double* dCoef)
{
    if (pointCount <= 1) return;

    int lastIndex = pointCount - 1;
    int baseOffset;
    int i, i1, j;
    double segWidth, segSlope, priorSegWidth, priorSegSlope, priorSegWidthThird, segWidthSum, segSlopeDiff;

    // Long strings of constant values can cause underflow problems;
    // remove such strings from the start.
    i = lastIndex; // default if no variation found
    for (int ii = 1; ii <= lastIndex; ii++) {
        bCoef[ii] = 0.0;
        cCoef[ii] = 0.0;
        dCoef[ii] = 0.0;
        if (1.0e15 * std::fabs(yIn[ii + 1] - yIn[ii]) > std::fabs(yIn[ii])) {
            i = ii;
            break;
        }
    }

    baseOffset = i - 1;
    j = 0; // will be set in loop
    for (i1 = i; i1 <= lastIndex; i1++) {
        j = pointCount + i - i1;
        bCoef[j] = 0.0;
        cCoef[j] = 0.0;
        dCoef[j] = 0.0;
        if (1.0e15 * std::fabs(yIn[j - 1] - yIn[j]) > std::fabs(yIn[j])) break;
    }

    int useCount = j - baseOffset;
    if (useCount <= 1) return;
    lastIndex = useCount - 1;

    segWidth = xPts[baseOffset + 2] - xPts[baseOffset + 1];
    segSlope = (yIn[baseOffset + 2] - yIn[baseOffset + 1]) / segWidth;

    if (useCount == 2) {
        // Special case for pointCount = 2
        dCoef[baseOffset + 1] = 0.0;
        bCoef[baseOffset + 1] = segSlope;
        return;
    }

    for (i = 2; i <= lastIndex; i++) {
        priorSegWidth = segWidth;
        segWidth = xPts[baseOffset + i + 1] - xPts[baseOffset + i];
        priorSegSlope = segSlope;
        segSlope = (yIn[baseOffset + i + 1] - yIn[baseOffset + i]) / segWidth;
        priorSegWidthThird = priorSegWidth / 3.0;
        dCoef[baseOffset + i - 1] = priorSegWidthThird * bCoef[baseOffset + i - 1];
        segWidthSum = priorSegWidth + segWidth;
        segSlopeDiff = segSlope - priorSegSlope;

        bCoef[baseOffset + i] = 1.0 / ((2.0 / 3.0) * segWidthSum - priorSegWidthThird * dCoef[baseOffset + i - 1]);
        cCoef[baseOffset + i] = segSlopeDiff - dCoef[baseOffset + i - 1] * cCoef[baseOffset + i - 1];
    }
    dCoef[baseOffset + lastIndex] = 0.0;

    // Back substitution
    for (i1 = 2; i1 <= lastIndex; i1++) {
        i = lastIndex + 2 - i1;
        cCoef[baseOffset + i] = bCoef[baseOffset + i] * cCoef[baseOffset + i]
                     - dCoef[baseOffset + i] * cCoef[baseOffset + i + 1];
    }

    // Compute bCoef, cCoef, dCoef from the tridiagonal solution
    for (i = 1; i <= lastIndex; i++) {
        segWidth = xPts[baseOffset + i + 1] - xPts[baseOffset + i];
        dCoef[baseOffset + i] = (cCoef[baseOffset + i + 1] - cCoef[baseOffset + i]) / (3.0 * segWidth);
        bCoef[baseOffset + i] = (yIn[baseOffset + i + 1] - yIn[baseOffset + i]) / segWidth
                     - (segWidth * dCoef[baseOffset + i] + cCoef[baseOffset + i]) * segWidth;
    }

    // Define the pointCount'th cubic as the (pointCount-1)'th cubic
    dCoef[baseOffset + useCount] = dCoef[baseOffset + lastIndex];
    bCoef[baseOffset + useCount] = bCoef[baseOffset + lastIndex]
                    + (2 * cCoef[baseOffset + lastIndex] + 3 * dCoef[baseOffset + lastIndex] * segWidth) * segWidth;

    return;
}

// ==========================================================================
// SYSERR — Fatal error handler (lines 5135-5146)
// ==========================================================================
void fatalError()
{
    std::printf("\n\n\n$*$*$*$*$*$  SYSERR CALLED     "
                "$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*\n");
    std::exit(9876);
}


// ============================================================================
// folded in from source_misc.cpp: cubMap, epsLon, linLsq
// ============================================================================
void cubMap(int mapType, double xLow, double xMidIn, double xHigh, double gamma,
            double* args, double* weights, int nPts)
{
    // Cubic-sinh mapping of [-1,1] into [xLow,xHigh]
    // mapType ∈ {1=cubic-sinh, 2=rational-sinh} — cases 0 (linear) and 3
    // Arrays args, weights are 1-based.

    double tau, xLen, xSum, a, b, c, d, tu, xi, sh, denom;
    double xMid = xMidIn;

    // Compute arcsinh(gamma)
    if (gamma > 1.0e-6)
        tau = std::log(gamma + std::sqrt(gamma*gamma + 1.0));
    else
        tau = gamma * (1.0 - gamma*gamma / 6.0);

    gaussL(nPts, args, weights);
    xLen = xHigh - xLow;
    xSum = xLow + xHigh;

    switch (mapType) {
    case 1: // Cubic-sinh mapping
        xMid = std::max(xMid, xLow + xLen/7.0);
        xMid = std::min(xMid, 0.5*xSum);
        a = 0.5*xSum - xMid;
        b = 0.5*xLen;
        c = 0.5*xSum;
        for (int i = 1; i <= nPts; i++) {
            tu = tau * args[i];
            xi = std::sinh(tu) / gamma;
            args[i] = a*(xi*xi - 1.0)*(xi + 1.0) + b*xi + c;
            weights[i] = weights[i] * (tau/gamma) * std::cosh(tu)
                   * ((3.0*xi - 1.0)*(xi + 1.0)*a + b);
        }
        return;

    case 2: // Rational-sinh mapping
        a = -xMid * xLen;
        b = xLen;
        c = xMid*xSum - 2.0*xLow*xHigh;
        d = xSum - 2.0*xMid;
        for (int i = 1; i <= nPts; i++) {
            tu = tau * args[i];
            sh = std::sinh(tu);
            denom = b - (d/gamma)*sh;
            args[i] = (-a + (c/gamma)*sh) / denom;
            weights[i] = weights[i] * (tau/gamma) * std::cosh(tu) * (b*c - a*d) / (denom*denom);
        }
        return;
    }
}

// epsLon — Wynn epsilon algorithm for complex sequences.
// xIn(2,pointCount): input partial sums as real*8 pairs (real, imag).
// fRet(2): output accelerated estimate (real, imag).
void epsLon(double* xIn, int pointCount, double* fRet, double& relativeError)
{
    using c16 = std::complex<double>;
    const double big = 1.0e38;
    const double eps = 1.0e-5;

    if (pointCount > 2000) {
        std::fprintf(stderr, "epsLon: N > 2000 NOT ALLOWED %d\n", pointCount);
        std::exit(9876);
    }

    std::vector<c16> x(pointCount + 1);
    for (int i = 1; i <= pointCount; i++)
        x[i] = c16(xIn[2*(i-1)], xIn[2*(i-1)+1]);

    auto approxAbs = [](c16 z) { return std::abs(std::real(z)) + std::abs(std::imag(z)); };
    auto nonZero = [&](c16 z) { return approxAbs(z) != 0.0; };
    auto isZero = [&](c16 z) { return approxAbs(z) == 0.0; };

    if (pointCount < 5) {
        if (pointCount > 0) { fRet[0] = std::real(x[pointCount]); fRet[1] = std::imag(x[pointCount]); }
        return;
    }

    double acc = std::max(1.0e-8, eps * eps);
    int n = pointCount;
    c16 fIn, t, w1, w2, w3, w4, w5, w6, w7;
    int windowStart = 4, outIndex = 0;
    bool isAlternateMode = false, isAlternateModeNext = false;

    // Outer convergence-acceleration loop
    while (true) {
    if (n <= 0) return;
    fIn = x[n];
    fRet[0] = std::real(fIn);
    fRet[1] = std::imag(fIn);
    if (n == 1) return;
    relativeError = 0;
    for (int i = 1; i <= n; i++) {
        double d = approxAbs(fIn - x[i]);
        if (d > relativeError) relativeError = d;
    }
    relativeError = relativeError / (approxAbs(fIn) + Constants::smlNum);
    if (relativeError <= acc) return;
    acc = eps;
    if (n < 6) return;

    // L2: initial epsilon values
    isAlternateMode = false; isAlternateModeNext = false;
    w1 = c16(big, 0);
    w7 = x[4] - x[3];
    if (nonZero(w7)) w1 = c16(1.0, 0) / w7;
    w5 = c16(big, 0);
    w7 = x[2] - x[1];
    if (nonZero(w7)) w5 = c16(1.0, 0) / w7;
    w4 = x[3] - x[2];
    if (!nonZero(w4)) {
        w4 = c16(big, 0);
        t  = x[2];
        w2 = x[3];
        w3 = c16(big, 0); // L8
    } else {
        w4 = c16(1.0, 0) / w4;
        t  = c16(big, 0);
        w7 = w4 - w5;
        if (nonZero(w7)) t = x[2] + c16(1.0, 0) / w7;
        w2 = w1 - w4;
        if (!nonZero(w2)) {
            w2 = c16(big, 0);
            isAlternateModeNext = (std::real(t) != big);
            w3 = w4;
        } else {
            w2 = x[3] + c16(1.0, 0) / w2;
            w7 = w2 - t;
            if (isZero(w7))
                w3 = c16(big, 0);
            else
                w3 = w4 + c16(1.0, 0) / w7;
        }
    }

    isAlternateMode = isAlternateModeNext;
    isAlternateModeNext = false;
    windowStart = 4;

    bool restartOuter = false;

    // compute x[outIndex] via the w6/w2 update rule with two
    // sub-modes selected by isAlternateMode. All paths inside set x[outIndex] (and
    // possibly windowStart, n, restartOuter) and return to the caller, which
    // then breaks out of the do-while body.
    auto computeViaL28 = [&](int i, int outIndex) {
        if (isAlternateMode) {
            if (std::real(w2) == big) { x[outIndex] = w5 + t - x[i-2]; return; }
            w7 = w5/(w2-w5) + t/(w2-t) + x[i-2]/(x[i-2]-w2);
            // divergent denominator → mark x[outIndex] as BIG
            if (isZero(w7 + c16(1.0, 0))) { x[outIndex] = c16(big, 0); windowStart = i; return; }
            x[outIndex] = w7 * w2 / (c16(1.0, 0) + w7);
            return;
        }
        w7 = w6 - w3;
        // divergent denominator → mark x[outIndex] as BIG
        if (isZero(w7)) { x[outIndex] = c16(big, 0); windowStart = i; return; }
        x[outIndex] = w2 + c16(1.0, 0) / w7;
        if (approxAbs(x[outIndex]) < 1.0e-10 * approxAbs(w2) && nonZero(x[outIndex])) { n = outIndex - 1; restartOuter = true; return; }
        if (std::real(w2) == big) { x[outIndex] = c16(big, 0); windowStart = i; }
    };

    restartOuter = false;
    for (int i = 5; i <= n; i++) {
        // per-iteration body wrapped in do-while-false so the 5
        // former `goto L39` sites become `break` to fall through to the
        // shared epilogue below.
        do {
            outIndex = i - windowStart;
            w4 = c16(big, 0);
            w5 = x[i-1];
            w7 = x[i] - x[i-1];
            if (!isZero(w7)) {
                w4 = c16(1.0, 0) / w7;
                // w1 == BIG → can't combine slopes; mark x[outIndex]=w2
                if (std::real(w1) == big) { w6 = c16(big, 0); isAlternateModeNext = false; x[outIndex] = w2; break; }
                w6 = w4 - w1;
                // split into "large slope change" + "small change w/ nonZero w6";
                // both fall through naturally to the L22 body below. The else-of-else (w6==0
                // but isAlternateModeNext set) takes the special w5=BIG path and breaks before L22.
                if (!(approxAbs(w6) > 1.0e-12 * approxAbs(w4))) {
                    isAlternateModeNext = true;
                    if (!nonZero(w6)) {
                        w5 = c16(big, 0);
                        w6 = w1;
                        if (std::real(w2) != big) { computeViaL28(i, outIndex); break; }
                        isAlternateModeNext = false;
                        x[outIndex] = w2;
                        break;
                    }
                }
                // Compute new w5 via (x[i-1] + 1/w6) refinement.
                w5 = x[i-1] + c16(1.0, 0) / w6;
                if (approxAbs(w5) < 1.0e-10 * approxAbs(x[i-1]) && nonZero(w5)) { n = outIndex - 1; restartOuter = true; break; }
            }
            w7 = w5 - w2;
            if (nonZero(w7)) {
                w6 = w1 + c16(1.0, 0) / w7;
                computeViaL28(i, outIndex);
            } else {
                // w7 zero → can't update; mark x[outIndex]=w2
                w6 = c16(big, 0); isAlternateModeNext = false; x[outIndex] = w2;
            }
        } while (false);
        if (restartOuter) break;
        // iteration epilogue
        w1 = w4;
        t  = w2;
        w2 = w5;
        w3 = w6;
        isAlternateMode = isAlternateModeNext;
        isAlternateModeNext = false;
    } // end DO 40

    if (restartOuter) continue;  // was: 2 internal goto L_OUTER sites
    n -= windowStart;
    // was: goto L_OUTER — implicit continue at end of while-true
    }  // end while (outer convergence-acceleration loop)
}

void linLsq(int fitForm, int pointCount, double* xVals, double* sVals, double& cVal,
             double& aVal, double& b, double& chiSq, int debugSwitch)
{
    double g1=0, g2=0, g11=0, g12=0, g22=0;
    // f/f2 init to 0 silences -Wmaybe-uninitialized for the switch
    // default-fallthrough that no caller exercises (fitForm ∈ {1..5} at
    // every callsite). The G-accumulators below would produce 1/0 = inf
    // in the fallthrough path; init keeps the failure mode bounded.
    double scale=1, x, s, wt, delta, a, c;
    double f2 = 0, f = 0;

    if (fitForm==2 || fitForm==4) scale = 0.5*(sVals[1]+sVals[2]);
    if (fitForm==4) scale = std::pow(0.5*(xVals[1]+xVals[2]), b) / scale;

    for (int i=1; i<=pointCount; i++) {
        x = xVals[i];
        s = sVals[i];
        switch (fitForm) {
        case 1: f2 = -x;             f = std::log(std::fabs(s));       break;
        case 2: f2 = std::exp(b*x);      f = 1.0/s;              break;
        case 3: f2 = std::log(std::fabs(x));   f = std::log(std::fabs(s));      break;
        case 4: f2 = x;               f = s/std::pow(std::fabs(x), b); break;
        case 5: f2 = std::pow(std::fabs(x), b); f = s;             break;
        }
        f = scale * f;
        wt = 1.0 / (f*f);
        g1  = g1  + wt*f;
        g2  = g2  + wt*f*f2;
        g11 = g11 + wt;
        g12 = g12 + wt*f2;
        g22 = g22 + wt*f2*f2;
    }

    delta = g11*g22 - g12*g12;
    a = (g22*g1 - g12*g2) * (1.0/(scale*delta));
    c = (g11*g2 - g12*g1) * (1.0/(scale*delta));

    chiSq = pointCount - (a*g1 + c*g2) * scale;
    aVal = a;
    cVal = c;
    if (debugSwitch) std::printf(" linLsq:%16.9G%14.6G%14.6G%16.9G%11.2G\n", b, aVal, cVal, chiSq, scale);
}

