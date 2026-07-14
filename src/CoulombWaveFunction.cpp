// CoulombWaveFunction.cpp — Coulomb wavefunctions F_L/G_L, phases, and integrals,
// implemented as CoulombWaveFunction static methods.

#include "CoulombWaveFunction.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

// ============================================================================
// SECTION 1: CoulombWaveFunction::computeFG — RCWFN
// Computes regular (F) and irregular (G) Coulomb wavefunctions
// Original: Barnett, Feng, Steed, Goldfarb — cpC 8 (1974) 377-395
// ============================================================================
// Original: A.R. Barnett, D.H. Feng, J.W. Steed, L.J.B. Goldfarb
//           Computer Physics Communications 8 (1974) 377-395


// label_70: zero lMin/lMin1 and, for large positive eta, bump the region
// code to 3 (Maclaurin). Shared by the four region-routing branches below;
// pure integer/region-code logic — no wavefunction math.
static void applyLabel70(int& region, int& lMin, int& lMin1, double eta, double rho)
{
    lMin = 0;
    lMin1 = 1;
    if (!(eta < 10 || rho <= eta)) region = 3;
}

void CoulombWaveFunction::computeFG(double rho, double eta, int minL, int maxL,
           double* fCArg, double* fCpArg, double* gCArg, double* gCpArg,
           double accuracy, int& returnCode)
{
    // Fortran arrays are 1-based: FC(1)..FC(lMax+1). The caller passes
    // a 0-based pointer where [0]=L_min position. Adjust to 1-based
    // so that fC[l+1] in C++ = Fortran FC(L+1).
    double* fC  = fCArg  - 1;
    double* fCp = fCpArg - 1;
    double* gC  = gCArg  - 1;
    double* gCp = gCpArg - 1;
    // implicit real*8 (a-h, o-z)
    double acc, ar, ai, b1, b2, b3, br, bi, c, d, del, delInverse;
    double gpFactor, di, dp, dq, dr, eta2, f, fMag, g, gp, h;
    double p, pL, pLSave, q, r, rhoUse, s, stop, sum, t, turn;
    double w, twoEta, x, xLL1, accHalf;
    int i, region, l, lMax, lMin, lMin1, ll, lp;
    bool isFirstIteration;

    // Machine-dependent constants (IEEE / RS6000)
    double veryBig = 1.79769e+308;
    double big    = 1.0e+300;
    double small  = 1.0e-300;
    double smallN = -690.775527898214;
    double precisionRt3 = 2.7e-5;
    double precision = 2.0e-14;
    double precisionLn = 32.0e0;

    double PI = 3.14159265358979323800e+00;

    // Initialize variables that may be printed in error path before assignment
    rhoUse = rho;
    p = 0.0; q = 0.0; r = 0.0; t = 0.0; x = 0.0; sum = 0.0; i = 0;

    // Error reporting helper (replaces goto error_exit / goto range_error)
    auto reportError = [&](int errorCode) {
        returnCode = errorCode;
        printf("0***RCWFIN IRET = %5d      INPUT = %20.10G%20.10G%10d%10d%15.5G\n",
               returnCode, rho, eta, minL, maxL, accuracy);
        printf(" ***%18.8G%18.8G%18.8G%18.8G%18.8G%18.8G\n",
               rhoUse, p, q, r, t, x);
        printf(" %7d%18.8G\n", i, sum);
    };

    // Here we limit accuracy to reasonable values for the machine

    acc = accuracy;
    acc = std::max(acc, precision);
    acc = std::min(acc, precisionRt3);

    lMax = maxL;
    lMin = minL;
    lMin1 = lMin + 1;
    xLL1 = (double)lMin * (double)lMin1;
    eta2 = eta * eta;

    // Determine which region we are in
    // For RHO < .45, Q of P+IQ is poorly determined so we don't use it.
    // Except that for large negative eta the Maclauren series also has problems

    // Region routing: sets region and adjusts lMin/lMin1/xLL1 before label_100
    bool goTo100 = false;  // whether to skip directly to the CF section

    if (rho > 0.45e0) {
        // label_20 block
        turn = eta + std::sqrt(eta2 + xLL1);
        region = 1;
        if (rho >= turn - 1.0e-4) {
            goTo100 = true;
        }
        // We are inside the turning point for minL, can we get outside
        // of it by reducing minL. (This is always possible for eta < 0).
        else if (rho < eta + std::fabs(eta)) {
            // Must use a different method to supplement the bad IQ value.
            // Always start with lMin = 0 for simplicity.
            // Note only eta > 0 gets to here (except when RHO < .45)
            // label_60 → label_70 → label_80
            region = 2;
            applyLabel70(region, lMin, lMin1, eta, rho);
        } else {
            // Yes, reduce lMin to get outside turning point
            lMin = (int)(0.5 * (std::sqrt(1 + 4 * ((rho - eta) * (rho - eta) - eta2)) - 1));
            lMin1 = lMin + 1;
            // fall through to label_80
        }
    } else {
        // RHO <= 0.45
        bool goLabel60 = false;
        if (eta >= 0) {
            // label_10
            if (rho > 0.005e0) {
                goLabel60 = true;
            } else {
                region = 5;
                // label_70
                applyLabel70(region, lMin, lMin1, eta, rho);
            }
        }
        // eta < 0
        else if (-eta * rho > 7) {
            // label_20
            turn = eta + std::sqrt(eta2 + xLL1);
            region = 1;
            if (rho >= turn - 1.0e-4) {
                goTo100 = true;
            } else if (rho < eta + std::fabs(eta)) {
                goLabel60 = true;
            } else {
                lMin = (int)(0.5 * (std::sqrt(1 + 4 * ((rho - eta) * (rho - eta) - eta2)) - 1));
                lMin1 = lMin + 1;
            }
        } else if (rho > 0.005e0) {  // label_10
            goLabel60 = true;
        } else {
            region = 5;
            // label_70
            applyLabel70(region, lMin, lMin1, eta, rho);
        }
        if (goLabel60 && !goTo100) {
            // label_60: Must use a different method
            region = 2;
            applyLabel70(region, lMin, lMin1, eta, rho);
        }
    }

    // label_80: update xLL1
    if (!goTo100) {
        xLL1 = (double)lMin * (double)lMin1;
    }

    // Section routing: 0=continue to P+IQ, 1=maclaurin, 2=normalize
    int nextSection = 0;

    // Here we compute F'/F for L = maxL
    // We then recurse down to lMin to generate the unnormalized F's
    // This section is used for all RHO.
    // label_100 / label_105 entry point

    // The label_500 block needs to loop back to label_105 with different rhoUse.
    // We model this with a loop that runs at most twice.
    bool doLabel500 = false;  // will be set true if region==3 after label_200

    for (int outerPass = 0; outerPass < 2; outerPass++) {

        if (outerPass == 0) {
            pL = lMax + 1;
            rhoUse = rho;
        } else {
            // label_500: eta > 15 and eta < RHO < 2*eta
            // Find G and G' for lMin at RHO=2*eta using the CF method.
            rhoUse = eta + eta;
            pL = lMin + 1;
            region = 4;
        }

        // label_105
        pLSave = pL;

        // label_110: Continued fraction for R = FP(lMax)/F(lMax)
        // Outer retry loop: if d passes near zero, increment pL and retry
        bool cfOk = false;
        while (!cfOk) {
            isFirstIteration = true;
            r  = eta / pL + pL / rhoUse;
            dq = (eta * rhoUse) * 2.0 + 6 * pL * pL;
            dr = 12 * pL + 6;
            del = 0.0;
            d   = 0.0;
            f   = 1.0;
            x   = (pL * pL - pL + (eta * rhoUse)) * (2.0 * pL - 1.0);
            ai  = rhoUse * pL * pL;
            di  = (2 * pL + 1) * rhoUse;

            bool converged = false;
            bool nearZero = false;
            for (i = 1; i <= 100000; i++) {
                h = (ai + rhoUse * eta2) * (rhoUse - ai);
                x = x + dq;
                d = d * h + x;

                // If we pass near a zero of the divisor, start over at larger pL
                if (std::fabs(d) > precisionRt3 * std::fabs(dr)) {
                    // label_130
                    d = 1 / d;
                    dq = dq + dr;
                    dr = dr + 12;
                    ai = ai + di;
                    di = di + 2 * rhoUse;
                    del = del * (d * x - 1.0);
                    if (isFirstIteration) del = -rhoUse * (pL * pL + eta2) * (pL + 1.0) * d / pL;
                    isFirstIteration = false;
                    r = r + del;
                    if (d < 0.0) f = -f;
                    if (std::fabs(del) < std::fabs(r * acc)) {
                        converged = true;
                        break;  // goto label_140
                    }
                } else {
                    // Near zero: bump pL and restart CF
                    pL = pL + 1;
                    if (pL < pLSave + 10) {
                        nearZero = true;
                        break;  // restart while loop (label_110)
                    }
                    reportError(5); return;
                }
            }

            if (nearZero) continue;  // restart while loop

            if (!converged) {
                reportError(6); return;
            }

            cfOk = true;
        }
        // label_140: R has converged; did we increase lMax (pL > pLSave)?

        if (pL != pLSave) {
            // Recurse down on R to lMax
            // Here the only part of F that is of interest is the sign
            pL = pL - 1;
            // label_150 loop
            while (pL > pLSave) {
                d = eta / pL + pL / rhoUse;
                f = (r + d) * f;
                r = d - (1 + eta2 / (pL * pL)) / (r + d);
                pL = pL - 1;
            }
        }

        // label_160: Now have R(lMax, RHO) or if region=4, R(lMin, 2*eta)
        if (region == 4) {
            // Skip straight to label_210 (P+IQ computation for region=4)
            // (falls through to P+IQ section below, outside outerPass loop)
            break;
        }

        fC [lMax + 1] = f;
        fCp[lMax + 1] = f * r;

        if (lMax != lMin) {
            // Downward recursion to lMin for F and FP, arrays GC,gCp are storage
            l  = lMax;
            pL = lMax;
            ar = 1 / rho;
            for (lp = lMin1; lp <= lMax; lp++) {
                gC [l + 1] = eta / pL + pL * ar;
                gCp[l + 1] = std::sqrt((eta / pL) * (eta / pL) + 1);
                fC [l]     = (gC[l + 1] * fC[l + 1] + fCp[l + 1]) / gCp[l + 1];
                fCp[l]     =  gC[l + 1] * fC[l]     - gCp[l + 1] * fC[l + 1];
                pL = pL - 1;
                l  = l - 1;

                // If we are getting near an overflow, renormalize everything down
                if (std::fabs(fC[l + 1]) >= big) {
                    for (ll = l; ll <= lMax; ll++) {
                        fC [ll + 1] = small * fC [ll + 1];
                        fCp[ll + 1] = small * fCp[ll + 1];
                    }
                }
            }
            f = fC [lMin1];
            r = fCp[lMin1] / f;
        }

        // label_200: Here we find P + IQ = (G'+IF')/(G+IF)
        // This section is used in all cases except when 15 < eta < RHO < 2*eta

        if (region == 3) {
            // label_500 path: go to second outerPass
            doLabel500 = true;
            break;  // break outerPass loop, will re-enter as pass 1
        }

        if (region == 5) {
            nextSection = 1; // maclaurin
        }

        // Fall through to P+IQ computation (label_210)
        break;  // done with outerPass loop for region != 3,4,5
    }

    if (doLabel500) {
        // Second pass was requested (region=3 → label_500 → label_105 with rhoUse=2*eta, region=4)
        // Run the CF again for the second pass
        rhoUse = eta + eta;
        pL = lMin + 1;
        region = 4;
        pLSave = pL;

        bool cfOk2 = false;
        while (!cfOk2) {
            isFirstIteration = true;
            r  = eta / pL + pL / rhoUse;
            dq = (eta * rhoUse) * 2.0 + 6 * pL * pL;
            dr = 12 * pL + 6;
            del = 0.0;
            d   = 0.0;
            f   = 1.0;
            x   = (pL * pL - pL + (eta * rhoUse)) * (2.0 * pL - 1.0);
            ai  = rhoUse * pL * pL;
            di  = (2 * pL + 1) * rhoUse;

            bool converged2 = false;
            bool nearZero2 = false;
            for (i = 1; i <= 100000; i++) {
                h = (ai + rhoUse * eta2) * (rhoUse - ai);
                x = x + dq;
                d = d * h + x;

                if (std::fabs(d) > precisionRt3 * std::fabs(dr)) {
                    d = 1 / d;
                    dq = dq + dr;
                    dr = dr + 12;
                    ai = ai + di;
                    di = di + 2 * rhoUse;
                    del = del * (d * x - 1.0);
                    if (isFirstIteration) del = -rhoUse * (pL * pL + eta2) * (pL + 1.0) * d / pL;
                    isFirstIteration = false;
                    r = r + del;
                    if (d < 0.0) f = -f;
                    if (std::fabs(del) < std::fabs(r * acc)) {
                        converged2 = true;
                        break;
                    }
                } else {
                    pL = pL + 1;
                    if (pL < pLSave + 10) {
                        nearZero2 = true;
                        break;
                    }
                    reportError(5); return;
                }
            }
            if (nearZero2) continue;
            if (!converged2) { reportError(6); return; }
            cfOk2 = true;
        }

        // Recurse down R from pL to pLSave if needed
        if (pL != pLSave) {
            pL = pL - 1;
            while (pL > pLSave) {
                d = eta / pL + pL / rhoUse;
                f = (r + d) * f;
                r = d - (1 + eta2 / (pL * pL)) / (r + d);
                pL = pL - 1;
            }
        }
        // region == 4: fall through to label_210 (P+IQ)
    }

    // label_210: Now obtain P + i.Q for lMin from continued fraction (32)
    // Real arithmetic to facilitate conversion to IBM using REAL*8
    {
        p  = 0.0;
        q  = rhoUse - eta;
        pL = 0.0;
        ar = -(eta2 + xLL1);
        ai = eta;
        br = q + q;
        bi = 2.0;
        twoEta = eta + eta;
        dr =  br / (br * br + bi * bi);
        di = -bi / (br * br + bi * bi);
        dp = -(ar * di + ai * dr);
        dq =  (ar * dr - ai * di);

        // label_230: Loop and converge on P + IQ
        while (true) {
            p  = p + dp;
            q  = q + dq;
            pL = pL + 2.0;
            ar = ar + pL;
            ai = ai + twoEta;
            bi = bi + 2.0;
            d  = ar * dr - ai * di + br;
            di = ai * dr + ar * di + bi;
            t  = 1.0 / (d * d + di * di);
            dr =  t * d;
            di = -t * di;
            h  = br * dr - bi * di - 1.0;
            x  = bi * dr + br * di;
            t  = dp * h  - dq * x;
            dq = dp * x  + dq * h;
            dp = t;
            if (pL > 46000.0) { reportError(7); return; }
            if (std::fabs(dp) + std::fabs(dq) < (std::fabs(p) + std::fabs(q)) * acc) break;
        }
        p = p / rhoUse;
        q = q / rhoUse;

        // We now have R and P+IQ, is this enough
        if (region == 2) {
            nextSection = 1; // maclaurin
        }

      if (nextSection == 0) {
        // Solve for FP, G, gp and normalise F at L=lMin
        // Since this is for RHO > RHO(turn), F and G are reasonable numbers
        x = (r - p) / q;
        fMag = std::sqrt(1 / (q * (1 + x * x)));
        w = fMag / std::fabs(f);
        f = w * f;
        g = f * x;
        gp = r * g - 1 / f;

        if (region == 4) {
            // label_600: Taylor series from turning point to RHO
            del = rhoUse - rho;
            b1 = g;
            b2 = -del * gp;
            b3 = 0;
            g = b1 + b2;
            accHalf = acc / 2;
            delInverse = -1 / del;
            gpFactor = 3 * delInverse;
            x = del / rhoUse;
            ai = x + x;
            di = ai + ai;
            ar = 6;
            dr = 6;
            bool taylorOk = false;
            for (i = 1; i <= 10000; i++) {
                s = (ai * b3 + (x * del * del) * b1) / ar;
                ar = ar + dr;
                dr = dr + 2;
                ai = ai + di;
                di = di + 2 * x;
                g = g + s;
                gp = gp + gpFactor * s;
                if (g >= veryBig) { reportError(4); return; }
                gpFactor = gpFactor + delInverse;
                b1 = b2;
                b2 = b3;
                b3 = s;
                if (s < accHalf * g) { taylorOk = true; break; }
            }
            if (!taylorOk) { reportError(9); return; }

            // label_650: Here we have R = F'/F, G, G'
            // Use Wronskian as the 4th condition
            f = fC[lMin1];
            r = fCp[lMin1] / f;
            sum = 1 / (r * g - gp);
            w = sum / f;
            f = sum;
            nextSection = 2;
        }

        // region == 1: direct result
        if (nextSection == 0) nextSection = 2;
      } // end nextSection == 0 (P+IQ path)
    }

  if (nextSection == 1) {
    // maclaurin_series:
    // label_400: Here RHO < eta or RHO < 2*eta < 20 or RHO < .45
    // We use the Maclauren series to get F(L=0, eta, RHO)
    // First compute RHO*c(L=0, eta)
    {
        c = 2 * PI * eta;
        x = 0.0; t = 1.0;  // defaults for eta > 0 case
        if (std::fabs(c) > 0.5) {
            // label_410
            if (eta > 0) {
                // label_420
                x = -smallN - PI * eta;
                t = small;
            } else {
                c = -c;
                x = 0;
                t = 1;
            }
            // label_425
            if (c < precisionLn) c = c / (1 - std::exp(-c));
        } else {
            // Use Maclaurin expansion of X / (EXP(X)-1)
            x = 0;
            t = 1;
            ar = 1;
            br = c;
            ai = 1;
            c = 1;
            // label_405 loop
            while (true) {
                ai = ai + 1;
                ar = ar * br / ai;
                c = c + ar;
                if (std::fabs(ar) < acc * c) break;
            }
            c = 1 / c;
        }

        // label_430
        c = rho * std::sqrt(c);
        b1 = 1;
        b2 = eta * rho;
        sum = b1 + b2;
        ai = 6;
        di = 6;
        bool seriesOk = false;
        for (i = 1; i <= 10000; i++) {
            b3 = ((2 * eta * rho) * b2 - (rho * rho) * b1) / ai;
            ai = ai + di;
            di = di + 2;
            sum = sum + b3;
            stop = std::fabs(b1) + std::fabs(b2) + std::fabs(b3);
            b1 = b2;
            b2 = b3;
            if (std::fabs(sum) >= big) {
                x = x - smallN;
                sum = sum * small;
                b1 = b1 * small;
                b2 = b2 * small;
            }
            if (stop < acc * std::fabs(sum)) { seriesOk = true; break; }
        }
        if (!seriesOk) { reportError(8); return; }

        // label_450
        sum = (c * std::exp(x) * sum) * t;

        // Did it underflow?
        if (sum == 0) { reportError(4); return; }

        // We now have F (=sum), R, and P (P only if RHO > .005)
        // Use the Wronskian as the 4th condition
        w = sum / f;
        f = sum;

        if (region == 5) {
            // label_850: RHO < .005; cannot find P or Q; return only F, F'.
            fC [lMin1] = f;
            fCp[lMin1] = r * f;
            returnCode = 2;
            if (lMax == lMin) return;
            for (l = lMin1; l <= lMax; l++) {
                fC [l + 1] = w * fC [l + 1];
                fCp[l + 1] = w * fCp[l + 1];
            }
            // label_840
            if (std::fabs(fC[lMax + 1]) + std::fabs(fCp[lMax + 1]) == 0) returnCode = returnCode + 1;
            return;
        }

        x = (r - p) * f;
        if (std::fabs(x) > precisionRt3) {
            // label_480: Must include F**3, F**4 terms
            // Determine which sign of the root applies (G > F if eta >= 0)
            b1 = 0.5 / x;
            b2 = b1 * std::sqrt(1 - 4 * (x * f) * (x * f));
            g = b1 + b2;
            if (eta < 0) {
                sum = 1 / q - f * f;
                gp = b1 - b2;
                if (std::fabs(g * g - sum) > std::fabs(gp * gp - sum)) g = gp;
            }
            // label_490
            gp = p * g - x * f / g;
        } else {
            // label_340: F**3 and F**4 terms are less than machine precision
            g = 1 / x;
            gp = p * g;
        }
        nextSection = 2;
    }

  } // end maclaurin block

    // normalize_and_recur:
    // label_800: We now have F, R = F'/F, G, G' at lMin
    // Upward recursion from GC(lMin) and gCp(lMin), stored values are RHO
    // Renormalise FC,Fcp for each L-value
    gC [lMin1] = g;
    gCp[lMin1] = gp;
    fC [lMin1] = f;
    fCp[lMin1] = r * f;
    returnCode = 0;
    if (lMax == lMin) return;
    for (l = lMin1; l <= lMax; l++) {
        t         = gC[l + 1];
        gC [l + 1] = (gC[l] * gC [l + 1] - gCp[l]) / gCp[l + 1];
        gCp[l + 1] =  gC[l] * gCp[l + 1] - gC[l + 1] * t;
        fC [l + 1] = w * fC [l + 1];
        fCp[l + 1] = w * fCp[l + 1];
    }
    // label_840
    if (std::fabs(fC[lMax + 1]) + std::fabs(fCp[lMax + 1]) == 0) returnCode = returnCode + 1;
    return;

}

// ============================================================================
// SECTION 2: CoulombWaveFunction::asymptoticPhase — RCASYM
// Asymptotic series for Coulomb functions
// ============================================================================
// Asymptotic series for Coulomb functions


// caller (BoundState_coulombIntegrals) sets it once before each invocation
// and never reads the resulting value back.
void CoulombWaveFunction::asymptoticPhase(int L, double eta, double rho, int printLevel, double sigL,
            double* zetaPointer, double* phiPointer, double* dZetaPointer,
            double* fPointer, double* fpPointer, double* gPointer, double* gpPointer,
            double* z, double* dzSquared, double* s, double* zInv,
            double eps, int nMax, int& nTz, int& convergenceCode)
{
    static const double PI = 3.14159265358979300;

    // RHOT init to 0 silences -Wmaybe-uninitialized — the three if-chain
    // arms (eta>0, eta==0, eta<0) are exhaustive but GCC can't fold them.
    double rhoT = 0;
    double zeta, dZeta, phi, fLl, etaSquared, rhoSquared, bf1, bf2;
    double cf, bf, tf, d2F, factor, sqrtZeta;
    int i, j, jTop, jHigh, jjTop, zInvIndex, kSer;

    convergenceCode = 0;
    fLl = (double)L * (L + 1);
    etaSquared = eta * eta;
    rhoSquared = rho * rho;
    bf1 = 2.0 * (eta / rho);
    bf2 = fLl / rhoSquared;
    if (eta > 0) rhoT = eta * (1.0 + std::sqrt(1.0 + fLl / etaSquared));
    if (eta == 0) rhoT = std::sqrt(fLl);
    if (eta == 0 && L == 0) rhoT = 0.5;
    if (eta < 0) rhoT = std::fabs(eta) * (std::sqrt(1.0 + fLl / etaSquared) - 1.0);
    if (rho <= rhoT) {
        if (printLevel > -4) std::printf("\n -- RCASYM. ASYMPTOTIC EXPANSION INSIDE TP --\n");
        convergenceCode = -5;
        return;
    }

    z[1] = 1.0;
    s[1] = 1.0;
    dzSquared[1] = 0.0;
    zInv[1] = 1.0;
    zeta = z[1];
    dZeta = 0.0;
    phi = rho - eta * std::log(2.0 * rho) + sigL - 0.5 * L * PI;
    kSer = 2;
    nTz = 1;
    i = 2;

    // Main loop
    bool seriesDone = false;
    bool seriesError = false;
    bool seriesDiverge = false;
    while (true) {
        // L500 body
        z[i] = 0.0;
        jTop = (i + 1) / 2;
        s[i] = 0.0;
        dzSquared[i] = 0.0;
        cf = 0.0;
        // (i==2 skips the s[i] accumulation loop)
        if (i != 2) {
            for (j = 2; j <= jTop; j++) {
                factor = (j == (i - j + 1)) ? 1.0 : 2.0;
                s[i] = s[i] + factor * z[j] * z[i - j + 1];
            }
        }
        jHigh = i - 1;
        for (j = 1; j <= jHigh; j++)
            cf = cf + z[j] * s[i - j + 1];
        bf = -bf1 * z[i - 1];
        if (i > 2) bf = bf - bf2 * z[i - 2];
        // (i<=5 skips zInv computation)
        if (i > 5) {
            zInvIndex = i - 4;
            jjTop = zInvIndex - 1;
            zInv[zInvIndex] = 0.0;
            for (j = 1; j <= jjTop; j++)
                zInv[zInvIndex] = zInv[zInvIndex] - zInv[j] * z[zInvIndex - j + 1];
        }
        tf = 0.0;
        // (i<=4 skips dzSquared and tf computation)
        if (i > 4) {
            for (j = 3; j <= jTop; j++) {
                factor = 2.0 / rhoSquared;
                if (j == (i - j + 1)) factor = 1.0 / rhoSquared;
                dzSquared[i] = dzSquared[i] + factor * (j - 2) * (i - j - 1) * z[j - 1] * z[i - j];
            }
            jjTop = i - 4;
            for (j = 1; j <= jjTop; j++)
                tf = tf + zInv[j] * dzSquared[i - j + 1];
        }
        d2F = 0.0;
        if (i > 3) d2F = (double)(i - 2) * (i - 3) * (z[i - 2] / rhoSquared);
        z[i] = (-cf + bf + 0.75 * tf - 0.5 * d2F) * 0.5;
        if (printLevel >= 5) std::printf("%5d Z = %13.6G%13.6G%13.6G%13.6G%13.6G\n",
            i, z[i], cf, bf, tf, d2F);
        nTz = nTz + 1;
        s[i] = s[i] + 2.0 * z[i];
        zeta = zeta + z[i];
        if (std::fabs(z[i]) < eps) kSer = kSer - 1;
        dZeta = dZeta - (double)(i - 2) * z[i - 1] / rho;
        if (i != 2 && std::fabs((double)(i - 2) * z[i - 1] / rho) < eps) kSer = kSer - 1;
        if (i != 2) phi = phi - (z[i] * rho) / (double)(i - 2);
        if (kSer <= 0) {
            seriesDone = true;
            break;
        }
        i = i + 1;
        if (i > nMax) {
            seriesError = true;
            break;
        }
        if (i < 8) continue;
        if (std::fabs(z[i - 1]) > std::fabs(z[i - 7])) {
            seriesDiverge = true;
            break;
        }
        // else: continue loop
    }

    if (seriesDone) {
        if (printLevel > 0) std::printf("\n WITH %2d TERMS,  ZETA = %22.14G\n DZETA = %22.14G PHI = %22.14G\n",
            nTz, zeta, dZeta, phi);

        sqrtZeta = std::sqrt(zeta);
        *fPointer = std::sin(phi) / sqrtZeta;
        *gPointer = std::cos(phi) / sqrtZeta;
        factor = (0.5 * dZeta) / (zeta * zeta);
        *fpPointer = zeta * ((*gPointer) - factor * (*fPointer));
        *gpPointer = -zeta * ((*fPointer) + factor * (*gPointer));

        *zetaPointer = zeta;
        *phiPointer = phi;
        *dZetaPointer = dZeta;
        return;
    }

    // Error/divergence cases (L1000, L1020, L1050)
    if (seriesError && printLevel > -4)
        std::printf(" **** MORE THAN %4d TERMS NEEDED IN SERIES FOR ZETA ****\n", nMax);
    if (seriesDiverge && printLevel > -4)
        std::printf(" **** COULOMB SERIES turnS AT%6d TERMS ****\n", i);
    if (printLevel > -4) std::printf(" **** L, ETA, RHO, ZETA, LAST = %4d%15.5G%15.5G%15.5G%15.5G\n",
        L, eta, rho, zeta, z[i - 1]);
    convergenceCode = -5;
    return;
}





// ---------------------------------------------------------------------------
// Weak implementations for the extern routines that COULST calls but that
// have not yet been translated.  These forward to the stubs already in the
// code base where the name matches; for genuinely new signatures we provide
// minimal stubs here so that this translation unit compiles and links cleanly.
// The actual subroutines are expected to be translated in later phases.
// ---------------------------------------------------------------------------

// RTXLNX — finds real solution of a*X + b*LN(X) + c = 0
double CoulombWaveFunction::solveRTXLNX(double a, double b, double c, double acc) {
    if (a == 0.0) {
        if (b == 0.0) return 0.0;
        return std::exp(-c / b);
    }
    double bp = b / a;
    double cp = c / a;
    if (b == 0.0) return -cp;
    double x = -cp;
    if (x < 0.1) x = 0.1;
    for (int i = 1; i <= 100; i++) {
        double xOld  = x;
        double delta = -(x + bp * std::log(x) + cp) / (1.0 + bp / x);
        x = x + delta;
        if (x <= 0.0) x = 0.1 * xOld;
        if (std::fabs(delta) <= acc) return x;
    }
    std::printf("0**** COULD NOT CONVERGE IN RTXLNX%20.10G%20.10G%20.10G%20.10G\n",
                a, b, c, x);
    return x;
}


// CoulombWaveFunction::generateBasisIndex, ::setBasisFactors, ::setupFG —
// CC-only helpers (GENBNX/SETBFC/SETFG family). All three lost their last
// callers when the COUPSW CC dispatch block in CoulombWaveFunction_scattering

