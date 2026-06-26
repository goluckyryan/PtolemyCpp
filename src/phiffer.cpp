//
// PHIFFER - linkule for phiffer-computed bound states
// PHIFER  - generates single-particle radial function
// STOP_F  - print error message and stop (renamed to avoid clash with C stop)
//
// 2/13/07 - new linkule based on av18 somewhat
//

#include "ptolemy_types.h"
#include "Reaction.h"
#include "linkule.h"
#include "math/spline.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include "Constants.h"
#include "LinkulePlugin.h"

// phifer + stopF are file-static (called only from this TU) —
// forward decls below for use-before-def.
static void phifer(int itype, int ievOpt, int nodes, int l, int jTwo,
            double h2o2mu, double* params, bool prntSwitch, int nPts,
            double rMax, double* rGrid, double* phiByRl,
            int nVGrid, double* vOnGrid, double& vRGrid);
static void stopF(const char* msg, int code);

// ============================================================================
// ============================================================================

// PhifferLinkulePlugin — Phiffer single-particle bound-state linkule.
// linkuleInts is unused; persistent state lives in static locals.
struct PhifferLinkulePlugin : LinkulePlugin {
    void run(char8 alias, int* /*linkuleInts*/, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, Reaction& reaction) override;
};

void PhifferLinkulePlugin::run(char8 alias, int* /*linkuleInts*/, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart,
             double stepSize, int nPts, double* array1, double* array2,
             Reaction& reaction)
{
    // implicit real*8 (a-h, o-z)
    int l = L;          // run() param name -> former phiffer param name
    double& jp = J;     // bound-state j out-param

    // Local persistent variables (SAVE)
    static double params[15];  // 1-based: params(14) -> [15]
    static int nodes, zProj, zTarget;
    static double v, rv, rc, av, e, am;
    static double rhowb, alphawb;
    static double vso, rso, aso;

    // Local
    double undef, rt, rEnd, hbarc, hb2o2m, alphaCoulomb;
    int notDef;
    int itype, ievOpt;
    double x;

    callStatus = 0;

    undef = reaction.internalState.undefValue;
    notDef = NOTDEF_INT;

    // =========================================================================
    // Setup (requestCode = 1): check for errors
    // =========================================================================
    if (requestCode < 2) {
        if (potType != 6) {
            printf(" **** bad potType in %.8s %d\n", alias.data, potType);
            callStatus = -1;
            return;
        }

        // set jp, nodes, v, r, a, e, j, jsp, jst, spam
        if (l == notDef || (int)jp == notDef) {
            printf(" L or JP not defined: %d %d\n", l, (int)jp);
            callStatus = -1;
            return;
        }

        nodes = reaction.angMom.nNodes;
        v = reaction.opticalPotentialParams.V;
        if (v == undef) v = 60.0;
        rv = reaction.integrationGrid.R;
        rc = reaction.opticalPotentialParams.rC;
        av = reaction.opticalPotentialParams.A;
        e = reaction.energies.E;
        am = reaction.masses.aM;
        if (e == undef || am == undef) {
            printf(" E or red. mass undefined: %g %g\n", e, am);
            callStatus = -1;
            return;
        }
        if (rv == undef || av == undef) {
            printf(" R or A undefined: %g %g\n", rv, av);
            callStatus = -1;
            return;
        }
        if (rc == undef) rc = rv;
        rhowb = reaction.linkuleParams.PARAM_arr[1];
        alphawb = reaction.linkuleParams.PARAM_arr[2];
        if (alphawb == undef) alphawb = 0.0;
        if (alphawb == 0.0) rhowb = 1.0;
        if (rhowb == undef) {
            printf(" rho(wb) undefined, rho, alpha: %g %g\n", reaction.linkuleParams.PARAM_arr[1], reaction.linkuleParams.PARAM_arr[2]);
            callStatus = -1;
            return;
        }

        zProj = reaction.charges.zProj;
        zTarget = reaction.charges.zTarget;
        if (zProj == notDef || zTarget == notDef) {
            printf(" zp or zt not defined %d %d\n", zProj, zTarget);
            callStatus = -1;
            return;
        }

        vso = reaction.opticalPotentialParams.vSo;
        rso = reaction.opticalPotentialParams.rSo;
        aso = reaction.opticalPotentialParams.aSo;
        if (vso == undef) vso = 0.0;
        if (vso != 0.0) {
            if (rso == undef || aso == undef) {
                printf(" rso or vso undef\" vso, rso, aso = %g %g %g\n", vso, rso, aso);
                callStatus = -1;
                return;
            }
        } else {
            if (rso == undef) rso = 1.0;
            if (aso == undef) aso = 1.0;
        }

        return;
    }

    // =========================================================================
    // Printing (requestCode = 2)
    // =========================================================================
    if (requestCode == 2) {
        if (potType == 6) {
            // Printout for wave function calculation.
            printf(" Phiffer calculation of wave function\n");
            printf("L, nodes, jp =%4d%4d%4d/2\n", l, nodes, (int)jp);
            printf("E, mu =%10.3f%10.3f\n", e, am);
            printf("V(guess), V(convrg) =%10.3f%10.3f\n", v, -params[10]);
            printf("R, A =%10.3f%10.3f\n", rv, av);
            printf("w.b. rho, alpha =%10.3f%10.3f\n", rhowb, alphawb);
            printf("Vso, Rso, Aso =%10.3f%10.3f%10.3f\n", vso, rso, aso);
            printf("Zp, Zt, Rc =%4d%4d%10.3f\n", zProj, zTarget, rc);
        } else {
            printf(" ******* we should not be here****\n");
        }

        return;
    }

    // =========================================================================
    // Calculation (requestCode = 3)
    // =========================================================================
    am = reaction.masses.aM;
    e = reaction.energies.E;
    hbarc = Constants::hbar_c;

    hb2o2m = hbarc * hbarc / (2.0 * am);
    alphaCoulomb = zProj * zTarget * hbarc / Constants::fine_structure_inv;


    rt = rStart;
    rEnd = rStart + (nPts - 1) * stepSize;

    if (rt != 0.0) {
        printf(" rStart not 0: %g\n", rt);
        callStatus = -1;
        return;
    }

    for (int i = 1; i <= 14; i++) params[i] = 0.0;
    params[1] = rv;
    params[2] = av;
    params[3] = e;
    params[4] = rhowb;
    params[5] = alphawb;
    params[6] = rso;
    params[7] = aso;
    params[8] = vso;
    params[10] = v;
    params[11] = alphaCoulomb;
    params[12] = rc;

    // allocate rGrid
    double* rGrid = new double[nPts + 2];  // 1-based

    itype = 11;
    ievOpt = 1;

    phifer(itype, ievOpt, nodes, l, (int)jp, hb2o2m,
           params, true, nPts, rEnd, rGrid,
           array1, nPts, array2, x);

    // phifer returns phi/r^L. We want u = r^(L+1) * phi/r^L
    for (int i = 1; i <= nPts; i++) {
        array1[i] = std::pow(rGrid[i], l + 1) * array1[i];
    }

    reaction.opticalPotentialParams.V = -params[10];

    delete[] rGrid;
    return;
}

std::unique_ptr<LinkulePlugin> makePhifferPlugin() {
    return std::make_unique<PhifferLinkulePlugin>();
}


// ============================================================================
//
// Generates single-particle radial function.
// Returns phiByRl = phi / r**L
// Integral r^2 phi^2 deltaR = 1
// ============================================================================

static void phifer(int itype, int ievOpt, int nodes, int l, int jTwo,
            double h2o2mu, double* params, bool prntSwitch, int nPts,
            double rMax, double* rGrid, double* phiByRl,
            int nVGrid, double* vOnGrid, double& vRGrid)
{
    // implicit real*8 (a-h, o-z)

    bool infwellSwitch, varyE, scatSwitch, eminSwitch, emaxSwitch, gridSwitch;
    int nodeVal[3];    // 1-based
    double xVal[3];   // 1-based

    // Work arrays (1-based via vec.data() - 1). std::vector replaces the
    std::vector<double> workSpaceVector, psiRealVector, centrifugalVector, splineBVector, splineCVector, splineDVector;
    double *workSpace = nullptr, *psiReal = nullptr, *centrifugal = nullptr;
    double *splineB = nullptr, *splineC = nullptr, *splineD = nullptr;

    double wsRadius, wsAlpha, wsEnergy, rho, alpha, rSoL, aSoL, vSoL;
    double alphaFermi, wsVolume, alphaCoulomb, rc, rb, xa, xb;
    double small, deltaR, dx, dl, consSo, derivative;
    double eStart, vStart, vme;
    double ak2, aKappa, gg1, gg2, ff1, ff2, logDerivDiff;
    // logDerivDiff1/wsEnergy1/wsVolume1 set on loop==2, logDerivDiff2/wsEnergy2/wsVolume2 set on loop==1
    // (which then `continue`s before any read). The wsEnergy/wsVolume-update block
    // (~line 548/572) only runs on loop>=2 so both sets are populated by
    // then — but GCC can't follow the loop-counter discriminator. Init
    // to 0 to silence -Wmaybe-uninitialized.
    double logDerivDiff1 = 0, logDerivDiff2 = 0, wsEnergy1 = 0, wsEnergy2 = 0, wsVolume1 = 0, wsVolume2 = 0;
    double efac, efac0, emin, emax;
    double sum1, sum2, anorm, rms, fac;
    double r, vCoul, so, x, priorPsir, pi;
    int iu, it, ih, nPointsMinusOne, nPointsMinusTwo;
    int match, matchPoint, nodeLoop, nodeLast, nodeHave, nodeCount;
    int maxGrid, loop;

    if (prntSwitch) {
        printf("phifer input:%5d%3d%3d%3d%4d/2%6d%6.1f%10.4f\n",
               itype, ievOpt, nodes, l, jTwo, nPts, rMax, h2o2mu);
        printf(" %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
               params[1], params[2], params[3], params[4],
               params[5], params[6], params[7], params[8]);
        printf(" %10.4f%10.4f%10.4f%10.4f%10.4f\n",
               params[9], params[10], params[11], params[12], params[13]);
    }

    maxGrid = nPts;
    // gridSwitch not yet set, initialize
    iu = itype % 10;
    gridSwitch = (iu == 9);
    if (gridSwitch) maxGrid = std::max(maxGrid, nVGrid);

    // Resize work arrays and cache 1-based pointers (data() - 1 → ptr[1] is element 0).
    workSpaceVector.resize(maxGrid + 1);     workSpace     = workSpaceVector.data();
    psiRealVector.resize(maxGrid + 1);   psiReal   = psiRealVector.data();
    centrifugalVector.resize(maxGrid + 1); centrifugal = centrifugalVector.data();
    splineBVector.resize(maxGrid + 1);   splineB   = splineBVector.data();
    splineCVector.resize(maxGrid + 1);   splineC   = splineCVector.data();
    splineDVector.resize(maxGrid + 1);   splineD   = splineDVector.data();

    it = (itype / 10) % 10;
    ih = (itype / 100) % 10;
    infwellSwitch = (ih == 1);
    scatSwitch = (ih == 3);
    varyE = (ievOpt == 0);

    small = 1.e-10;
    deltaR = rMax / (nPts - 1);
    dx = deltaR * deltaR / (12.0 * h2o2mu);
    nPointsMinusOne = nPts - 1;
    nPointsMinusTwo = nPts - 2;
    for (int i = 1; i <= nPts; i++) {
        rGrid[i] = deltaR * (i - 1);
    }

    wsRadius = params[1];
    wsAlpha = params[2];
    wsEnergy = params[3];
    rho = params[4];
    alpha = params[5];
    rSoL = params[6];
    aSoL = params[7];
    vSoL = params[8];
    alphaFermi = params[9];
    wsVolume = params[10];
    alphaCoulomb = params[11];
    rc = params[12];
    rb = params[13];
    xa = rc / std::sqrt(12.0);
    xb = rb / std::sqrt(12.0);

    if (iu == 8) {
        // Return a h.o. solution
        if (nodes != 0) {
            std::exit(9876);
        }
        pi = std::acos(-1.0);
        anorm = std::pow(2.0 / wsAlpha, 2 * l + 3) * std::sqrt(2.0 / pi);
        for (int ll = 1; ll <= 2 * l + 1; ll += 2) {
            anorm = anorm / ll;
        }
        anorm = std::sqrt(anorm);
        printf(" H.O.; L= %d  A= %g  Norm= %g\n", l, wsAlpha, anorm);
        for (int i = 1; i <= nPts; i++) {
            x = rGrid[i] / wsAlpha;
            // We leave out the r**l because phifer returns phi/r**l
            phiByRl[i] = anorm * std::exp(-(x * x));
        }
        return;
    }

    if (scatSwitch) derivative = params[14];
    if (rho == 0.0) rho = 1.0;
    dl = l * (l + 1) * deltaR * deltaR / 12.0;
    consSo = 0.0;
    if (jTwo != -1) {
        consSo = -2.0 * vSoL * (jTwo * (jTwo + 2) / 4.0 - l * (l + 1) - 0.75) / aSoL;
    }

    if (gridSwitch) {
        naturalCubicSpline(nVGrid, &vRGrid, vOnGrid, splineB + 1, splineC + 1, splineD + 1);
        cubicSplineInterp(nVGrid, &vRGrid, vOnGrid, splineB + 1, splineC + 1, splineD + 1,
               nPts, rGrid + 1, workSpace + 1);
    }

    // Build potential and centrifugal arrays
    for (int i = 1; i <= nPts; i++) {
        r = deltaR * (i - 1);

        switch (iu) {
            case 1:
                workSpace[i] = 1.0 / (1.0 + std::exp((r - wsRadius) / wsAlpha));
                if (wsRadius > 0.0) workSpace[i] = workSpace[i] * (1.0 + (alphaFermi / (wsRadius * wsRadius)) * r * r);
                workSpace[i] = workSpace[i] - alpha * std::exp(-(r / rho) * (r / rho));
                break;
            case 2:
                workSpace[i] = std::exp(-((r - wsRadius) / wsAlpha) * ((r - wsRadius) / wsAlpha));
                break;
            case 9:
                // potential already in workSpace from spline interpolation
                break;
            default:
                stopF("PHIFER: invalid itype", itype);
                break;
        }

        workSpace[i] = dx * workSpace[i];
        centrifugal[i] = dl / ((r + 1.e-20) * (r + 1.e-20));
        if (consSo != 0.0) {
            x = std::exp((r - rSoL) / aSoL);
            so = consSo * x / ((r + 1.e-20) * (1.0 + x) * (1.0 + x));
            centrifugal[i] = centrifugal[i] + dx * so;
        }

        switch (it) {
            case 0:
                vCoul = 0.0;
                break;
            case 1:
                if (r > rc) {
                    vCoul = alphaCoulomb / r;
                } else {
                    vCoul = (alphaCoulomb / (2.0 * rc)) * (3.0 - (r / rc) * (r / rc));
                }
                break;
            case 2:
                if (r > 1.e-3) {
                    vCoul = alphaCoulomb * (1.0
                         - (0.5 / ((1.0 - (xa / xb) * (xa / xb)) * (1.0 - (xa / xb) * (xa / xb))))
                           * std::exp(-r / xb) * (2.0 + (r / xb) + 4.0 / (1.0 - (xb / xa) * (xb / xa)))
                         - (0.5 / ((1.0 - (xb / xa) * (xb / xa)) * (1.0 - (xb / xa) * (xb / xa))))
                           * std::exp(-r / xa) * (2.0 + (r / xa) + 4.0 / (1.0 - (xa / xb) * (xa / xb)))
                         ) / r;
                } else {
                    vCoul = (alphaCoulomb / (2.0 * std::pow(xb * xb - xa * xa, 3)))
                       * (xb * xb * xb * (xb * xb - 5.0 * xa * xa)
                        + xa * xa * xa * (5.0 * xb * xb - xa * xa));
                }
                break;
            default:
                stopF("PHIFER: invalid itype", itype);
                break;
        }
        centrifugal[i] = centrifugal[i] + dx * vCoul;
    }

    if (ievOpt == 9) {
        for (int i = 1; i <= nPts; i++) {
            phiByRl[i] = (wsVolume * workSpace[i] + centrifugal[i]) / dx;
        }
        return;
    }

    nodeLoop = 1;
    efac0 = 1.3;
    efac = efac0;
    nodeLast = 9999;
    nodeHave = 0;
    eminSwitch = false;
    emaxSwitch = false;
    emin = 0.0;
    emax = 0.0;
    if (!infwellSwitch) {
        if (!gridSwitch) wsVolume = -std::fabs(wsVolume);
        wsEnergy = std::fabs(wsEnergy);
    }

    // ===== Eigenvalue-search outer loop =====
    while (true) {
    eStart = wsEnergy;
    vStart = wsVolume;
    vme = (wsVolume * workSpace[1] / dx + wsEnergy) / h2o2mu;
    match = (int)(std::sqrt(-l * (l + 1) / vme) / deltaR + 0.5);
    match = std::max(match, nPts / 20);
    matchPoint = match + 1;

    // Converge on the eigenvalue
    for (loop = 1; loop <= 100; loop++) {

        // Integrate from outside in
        ak2 = -dx * wsEnergy;
        aKappa = std::sqrt(+wsEnergy / h2o2mu);
        if (infwellSwitch) {
            psiReal[nPts] = 0.0;
            psiReal[nPointsMinusOne] = 0.01;
        } else if (scatSwitch) {
            ak2 = -ak2;
            if (derivative == 0.0) {
                psiReal[nPts] = 0.0;
                psiReal[nPointsMinusOne] = 0.01;
            } else {
                psiReal[nPts] = 1.0;
                double d1 = derivative + 1.0 / rMax;
                double d2 = -aKappa * aKappa + l * (l + 1) / (rMax * rMax);
                psiReal[nPointsMinusOne] = 1.0 - d1 * deltaR + 0.5 * d2 * deltaR * deltaR;
            }
        } else {
            psiReal[nPts] = std::exp(-aKappa * rGrid[nPts]);
            psiReal[nPointsMinusOne] = std::exp(-aKappa * rGrid[nPointsMinusOne]);
        }

        double xaL = psiReal[nPts];
        double xbL = psiReal[nPointsMinusOne];
        double xd = wsVolume * workSpace[nPts] + centrifugal[nPts] - ak2;
        double xe = wsVolume * workSpace[nPointsMinusOne] + centrifugal[nPointsMinusOne] - ak2;
        for (int j = nPointsMinusTwo; j >= match; j--) {
            double xf = wsVolume * workSpace[j] + centrifugal[j] - ak2;
            psiReal[j] = ((2.0 + 10.0 * xe) * xbL - (1.0 - xd) * xaL) / (1.0 - xf);
            xaL = xbL;
            xbL = psiReal[j];
            xd = xe;
            xe = xf;
        }
        gg1 = xbL;
        gg2 = xaL;

        // Integrate from inside out
        psiReal[1] = 0.0;
        psiReal[2] = (1.0 + 12.0 * (wsVolume * workSpace[2] - ak2) / (4 * l + 6)) * std::pow(deltaR, l + 1);
        xbL = psiReal[2];
        xe = wsVolume * workSpace[2] + centrifugal[2] - ak2;
        double xf = wsVolume * workSpace[3] + centrifugal[3] - ak2;
        psiReal[3] = (2.0 + 10.0 * xe) * xbL / (1.0 - xf);
        if (l == 1) psiReal[3] = psiReal[3] + 2.0 * h2o2mu * dx / (1.0 - xf);
        xaL = xbL;
        xbL = psiReal[3];
        xd = xe;
        xe = xf;
        for (int j = 4; j <= matchPoint; j++) {
            xf = wsVolume * workSpace[j] + centrifugal[j] - ak2;
            psiReal[j] = ((2.0 + 10.0 * xe) * xbL - (1.0 - xd) * xaL) / (1.0 - xf);
            xaL = xbL;
            xbL = psiReal[j];
            xd = xe;
            xe = xf;
        }
        ff1 = xaL;
        ff2 = xbL;
        logDerivDiff = (ff2 * gg1 - ff1 * gg2) / (std::fabs(ff1 * gg1) + std::fabs(ff2 * gg2));

        if (loop == 1) {
            logDerivDiff2 = logDerivDiff;
            wsEnergy2 = wsEnergy;
            wsVolume2 = wsVolume;
            if (varyE) {
                wsEnergy = wsEnergy2 * 1.2;
            } else {
                wsVolume = wsVolume2 * 1.2;
            }
            continue;  // go to 80
        } else if (loop == 2) {
            logDerivDiff1 = logDerivDiff;
            wsEnergy1 = wsEnergy;
            wsVolume1 = wsVolume;
        } else {
            if (std::fabs(logDerivDiff) <= small) break;
            // Arithmetic IF replacements (three-way branch)
            bool toL65 = false;
            if (logDerivDiff * logDerivDiff1 <= 0.0) {
                // L54
                if (logDerivDiff * logDerivDiff2 <= 0.0) {
                    // L58
                    if (std::fabs(logDerivDiff1) - std::fabs(logDerivDiff2) > 0.0) {
                        // L65
                        toL65 = true;
                    }
                }
                // else: L60 — no skip
            }
            // L56
            else if (logDerivDiff * logDerivDiff2 <= 0.0) {
                // L65
                toL65 = true;
            }
            // L58
            else if (std::fabs(logDerivDiff1) - std::fabs(logDerivDiff2) > 0.0) {
                // L65
                toL65 = true;
            }
            // else: L60 — no skip

            if (!toL65) {
                logDerivDiff2 = logDerivDiff1;
                wsEnergy2 = wsEnergy1;
                wsVolume2 = wsVolume1;
            }
            logDerivDiff1 = logDerivDiff;
            wsEnergy1 = wsEnergy;
            wsVolume1 = wsVolume;
        }


        if (varyE) {
            wsEnergy = (logDerivDiff2 * wsEnergy1 - logDerivDiff1 * wsEnergy2) / (logDerivDiff2 - logDerivDiff1);
            if (wsEnergy < 0.0 && !infwellSwitch) {
                wsEnergy = 0.5 * std::min(wsEnergy1, wsEnergy2);
                if (wsEnergy < 0.001) {
                    if (prntSwitch)
                        printf(" wsEnergy tried for 0, start= %g\n", eStart);
                    wsEnergy = 2.0 * eStart;
                    continue;  // restart outer
                }
            }
            if (std::min(std::fabs(wsEnergy - wsEnergy1), std::fabs(wsEnergy - wsEnergy2))
                / std::fabs(wsEnergy1 - wsEnergy2) < 0.1)
                wsEnergy = 0.5 * (wsEnergy1 + wsEnergy2);
            if (eminSwitch && wsEnergy < emin) {
                wsEnergy = 0.5 * (emin + std::min(wsEnergy1, wsEnergy2));
            }
            if (emaxSwitch && wsEnergy < emax) {
                wsEnergy = 0.5 * (emax + std::max(wsEnergy1, wsEnergy2));
            }
        } else {
            wsVolume = (logDerivDiff2 * wsVolume1 - logDerivDiff1 * wsVolume2) / (logDerivDiff2 - logDerivDiff1);
            // Don't let wsVolume change sign, but it could be either sign
            if (wsVolume * vStart <= 0.0) {
                wsVolume = std::copysign(0.5 * std::min(std::fabs(wsVolume1), std::fabs(wsVolume2)), vStart);
            }
            if (std::min(std::fabs(wsVolume - wsVolume1), std::fabs(wsVolume - wsVolume2))
                / std::fabs(wsVolume1 - wsVolume2) < 0.1)
                wsVolume = 0.5 * (wsVolume1 + wsVolume2);
        }
    }  // end loop 80

    // Did not converge (only diagnose if the for-loop exited by completing all 100 iterations)
    if (loop > 100) {
        if (std::fabs(logDerivDiff) > 10000.0 * small) {
            printf("\n*********** did not converge on wavefunction ********************\n");
            printf("%25.15g%25.15g\n", wsEnergy, wsVolume);
            printf("%25.15g%25.15g%25.15g\n", wsEnergy1, wsVolume1, logDerivDiff1);
            printf("%25.15g%25.15g%25.15g\n", wsEnergy2, wsVolume2, logDerivDiff2);
            printf("%25.15g%25.15g%25.15g%25.15g\n", ff1, gg1, ff2, gg2);
            printf("******** phiffer could not converge ********\n");
        } else if (prntSwitch) {
            printf("phifer accepts poor converg%15.5g%15.5g\n", logDerivDiff, small);
        }
    }

    // Normalize inside r(matchPoint)
    fac = gg2 / psiReal[matchPoint];
    for (int i = 2; i <= matchPoint; i++) {
        psiReal[i] = fac * psiReal[i];
    }

    // Normalize and count nodes
    sum1 = 0.0;
    sum2 = 0.0;
    priorPsir = psiReal[4];
    nodeCount = 0;
    for (int i = 2; i <= nPts; i++) {
        sum1 = sum1 + psiReal[i] * psiReal[i];
        sum2 = sum2 + psiReal[i] * psiReal[i] * rGrid[i] * rGrid[i];
        phiByRl[i] = psiReal[i] / std::pow(rGrid[i], l + 1);
        if (psiReal[i] * priorPsir <= 0.0 && i > 4 && i < nPts) {
            nodeCount = nodeCount + 1;
            priorPsir = psiReal[i + 1];
        }
    }
    if (nodeCount == nodes) break;  // success — exit outer

    // Wrong number of nodes; try again
    if (prntSwitch) {
        printf(" got%3d nodes:  e =%15.5g%15.5g efac=%15.5g     v =%10.5f%10.5f%5d matching loops\n",
               nodeCount, eStart, wsEnergy, efac, vStart, wsVolume, loop);
    }
    nodeLoop = nodeLoop + 1;
    if (nodeLoop > 40) {
        printf(" *** node loop limit reached: %d %g %g %g %g %g\n",
               nodeCount, eStart, wsEnergy, efac, vStart, wsVolume);
        break;  // give up — exit outer
    }

    if (varyE) {
        // Varying e
        if (nodeCount == nodeLast) {
            efac = 1.1 * efac;
        } else {
            efac = efac0;
        }
        nodeLast = nodeCount;

        // Set barriers on e based on what we have found
        double xDelta = 0.1 * std::fabs(wsEnergy);
        if (nodeCount > nodes) {
            if (!eminSwitch || (eminSwitch && wsEnergy - xDelta > emin)) {
                eminSwitch = true;
                emin = wsEnergy - xDelta;
            }
        } else if (!emaxSwitch || (emaxSwitch && wsEnergy + xDelta < emax)) {
            emaxSwitch = true;
            emax = wsEnergy + xDelta;
        }

        // Accumulate eigenvalues and number of nodes.
        // doPredict toggles between "have enough data to extrapolate wsEnergy"
        // (run prediction + restart outer) and "still gathering"
        // (skip prediction, fall through to wsEnergy adjustment).
        bool doPredict = true;
        if (nodeHave == 0) {
            nodeHave = 1;
            nodeVal[1] = nodeCount;
            xVal[1] = wsEnergy;
            doPredict = false;
        } else if (nodeHave == 1) {
            if (nodeCount == nodeVal[1]) { doPredict = false; }
            else {
                nodeHave = 2;
                nodeVal[2] = nodeCount;
                xVal[2] = wsEnergy;
            }
        } else if (nodeCount == nodeVal[1] || nodeCount == nodeVal[2]) { doPredict = false; }
        else {
            // Replace the worse number of nodes
            int idx;
            if (std::abs(nodeVal[1] - nodes) > std::abs(nodeVal[2] - nodes)) {
                idx = 1;
            } else {
                idx = 2;
            }
            if (std::abs(nodeVal[idx] - nodes) > std::abs(nodeCount - nodes)) {
                nodeVal[idx] = nodeCount;
                xVal[idx] = wsEnergy;
            } else {
                doPredict = false;
            }
        }

        if (doPredict) {
            // Predict the new eigenval based on number of nodes found
            x = xVal[1] + ((xVal[2] - xVal[1]) / (double)(nodeVal[2] - nodeVal[1]))
                  * (nodes - nodeVal[1]);
            wsEnergy = x;
            continue;  // restart outer
        }

        // wsEnergy adjustment when not yet ready to predict
        if (nodeCount < nodes) {
            wsEnergy = std::min(wsEnergy, eStart);
            if (wsEnergy > 0.0) {
                wsEnergy = wsEnergy / efac;
            } else {
                wsEnergy = wsEnergy * efac;
            }
        } else {
            wsEnergy = std::max(wsEnergy, eStart);
            if (wsEnergy > 0.0) {
                wsEnergy = wsEnergy * efac;
            } else {
                wsEnergy = wsEnergy / efac;
            }
        }
    }
    // Varying v
    else if (nodeCount < nodes) {
        wsVolume = std::copysign(1.3 * std::max(std::fabs(wsVolume), std::fabs(vStart)), vStart);
    } else {
        wsVolume = std::copysign(std::min(std::fabs(wsVolume), std::fabs(vStart)) / 1.2, vStart);
    }
    continue;  // restart outer
    }  // end while (eigenvalue-search outer loop)

    // Have correct number of nodes
    phiByRl[1] = (rGrid[3] * rGrid[3] * phiByRl[2]
                   - rGrid[2] * rGrid[2] * phiByRl[3])
                 / (rGrid[3] * rGrid[3] - rGrid[2] * rGrid[2]);
    anorm = std::sqrt(sum1 * deltaR);
    for (int i = 1; i <= nPts; i++) {
        phiByRl[i] = phiByRl[i] / anorm;
    }
    rms = std::sqrt(sum2 / sum1);
    if (prntSwitch) {
        printf(" woods-saxon strength v =%15.10f     e(sep) =%15.10f;   rms radius =%8.5f%5d iters;\n",
               wsVolume, wsEnergy, rms, loop);
        printf(" convrg =%12.3g     phi/r**l:%18.9g%18.9g%18.9g%18.9g%18.9g\n",
               logDerivDiff, phiByRl[1], phiByRl[2], phiByRl[3], phiByRl[4], phiByRl[5]);
    }

    // If we completely failed to converge, stop
    if (std::fabs(logDerivDiff) > 1000.0 * small) {
        printf(" could not converge: %g %g\n", logDerivDiff, 1000.0 * small);
    }
    params[3] = wsEnergy;
    params[10] = wsVolume;

    // PTOLEMY: return potential shape
    for (int i = 1; i <= nVGrid; i++) {
        vOnGrid[i] = workSpace[i] / dx;
    }
    // std::vector destructors free workSpaceVector/psiRealVector/centrifugalVector/splineBVector/
    // splineCVector/splineDVector on return — no manual delete[] needed.
}


// ============================================================================
//
// Print an error message and stop.
// ============================================================================

static void stopF(const char* msg, int code)
{
    printf("\n %s\n", std::string(78, '*').c_str());
    printf(" *%20sStopping because of an error  !!%24s*\n", "", "");
    printf(" * %s : %10d%*s*\n", msg, code, (int)(74 - strlen(msg) - 13), "");
    printf(" %s\n\n", std::string(78, '*').c_str());

    fprintf(stderr, "\n");
    fprintf(stderr, " %s\n", std::string(78, '*').c_str());
    fprintf(stderr, " *%20sStopping because of an error  !!%24s*\n", "", "");
    fprintf(stderr, " * %s : %10d%*s*\n", msg, code, (int)(74 - strlen(msg) - 13), "");
    fprintf(stderr, " %s\n\n", std::string(78, '*').c_str());

    std::exit(9999);
}
