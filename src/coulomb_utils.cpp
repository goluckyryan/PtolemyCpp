// coulomb_utils.cpp — Coulomb potential between two uniform spheres
// (Poling et al. Phys.Rev.C13,648,1976): setVsq precomputes the per-sphere-pair
// coefficients, vcsq12 evaluates V(r) for a pair. DO NOT alter the math —
// the sphere-split kernel is bit-identical-critical.

#include "coulomb_utils.h"
#include "Constants.h"
#include <cstdio>
#include <cmath>

// File-scope arrays shared between setVsq and vcsq12
static double vcsq12R1[4] = {}, vcsq12R2[4] = {}, vcsq12Vc0[4] = {};
static double vcsq12A[4] = {}, vcsq12B[4] = {}, vcsq12C[4] = {};
static double vcsq12D[4] = {}, vcsq12E[4] = {}, vcsq12F[4] = {};
static double vcsq12X[4] = {}, vcsq12Y[4] = {};

void vcsq12(double rValue, double& x, int k)
{
    // Coulomb potential between two uniform spheres (Poling et al. Phys.Rev.C13,648,1976)
    // Uses coefficients precomputed by setVsq for sphere pair k (1-indexed).
    // Three regions: outside (point-charge), fully inside (quadratic), partial overlap (polynomial/r).
    double r1 = vcsq12R1[k];  // larger radius (enforced by setVsq)
    double r2 = vcsq12R2[k];  // smaller radius
    if (rValue >= r1 + r2) {
        // Outside both spheres: point-charge Coulomb
        x = vcsq12Vc0[k] / rValue;
    } else if (rValue <= r1 - r2) {
        // Fully inside larger sphere: quadratic formula
        x = vcsq12X[k] + vcsq12Y[k] * rValue * rValue;
    } else {
        // Partial overlap: polynomial / r
        double r = rValue;
        x = vcsq12A[k] + r*(vcsq12B[k] + r*(vcsq12C[k] + r*(vcsq12D[k]
              + r*(vcsq12E[k] + r*r*vcsq12F[k]))));
        x = x / r;
    }
    if (k <= 2 && !std::isfinite(x))
        std::fprintf(stderr, "VCSQ12_NAN K=%d rValue=%.5e R1=%.5e R2=%.5e VC0=%.5e X=%.5e\n",
            k, rValue, r1, r2, vcsq12Vc0[k], x);
}

void setVsq(double rr1, double rr2, int iz1, int iz2, int k)
{
    // Computes coefficients for Coulomb potential between two uniform spheres
    vcsq12R1[k] = rr1;
    vcsq12R2[k] = rr2;
    vcsq12Vc0[k] = iz1 * iz2 * (Constants::hbar_c / Constants::fine_structure_inv);
    if (vcsq12R1[k] < vcsq12R2[k]) {
        double swapTemp = vcsq12R2[k]; vcsq12R2[k] = vcsq12R1[k]; vcsq12R1[k] = swapTemp;
    }
    double r1Sq = vcsq12R1[k] * vcsq12R1[k];
    double r2Sq = vcsq12R2[k] * vcsq12R2[k];
    double r13 = r1Sq * vcsq12R1[k];
    if (vcsq12R2[k] != 0.0) {
        double r14 = r13 * vcsq12R1[k];
        double r15 = r14 * vcsq12R1[k];
        double r16 = r15 * vcsq12R1[k];
        double r23 = r2Sq * vcsq12R2[k];
        double r24 = r23 * vcsq12R2[k];
        double r25 = r24 * vcsq12R2[k];
        double r26 = r25 * vcsq12R2[k];
        double r12C = r13 * r23;
        double v = vcsq12Vc0[k] / r12C;
        vcsq12A[k] = (v * (r16 - 9*r14*r2Sq + 16*r12C - 9*r1Sq*r24 + r26)) / 32;
        vcsq12B[k] = -(3*v * (r15 - 5*r13*r2Sq - 5*r1Sq*r23 + r25)) / 20;
        vcsq12C[k] = (9*v * (r14 - 2*r1Sq*r2Sq + r24)) / 32;
        vcsq12D[k] = -(v * (r13 + r23)) / 4;
        vcsq12E[k] = (3*v * (r1Sq + r2Sq)) / 32;
        vcsq12F[k] = -v / 160;
    }
    if (vcsq12R1[k] == 0.0) return;
    vcsq12X[k] = 0.3 * (5*r1Sq - r2Sq) * (vcsq12Vc0[k] / r13);
    vcsq12Y[k] = -0.5 * (vcsq12Vc0[k] / r13);
}
