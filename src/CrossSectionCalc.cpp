// CrossSectionCalc.cpp — CrossSectionCalc class: shared cross-section machinery.
// Individual method groups live in the CrossSection_*.cpp companion files.

#include "CrossSectionCalc.h"
#include "math/angular_momentum_coeff.h"
#include "math/legendre.h"
#include "math/numeric_utils.h"
#include <cstdio>
#include <cmath>
#include "Constants.h"


// ============================================================================
// Merged from amplitude_calc.cpp
// ============================================================================

// Helper: Fortran 1-based access for raw pointer arrays.
// f1(ptr, i) accesses Fortran index i (1-based) from a 0-based C pointer.
#define f1(ptr, i) (ptr)[(i) - 1]

// Helper: Fortran column-major 2D access (1-based) for raw pointer.
// f2(ptr, dim1, r, c) accesses element (r,c) in array with first dim = dim1.
#define f2(ptr, dim1, r, c) (ptr)[((c) - 1) * (dim1) + (r) - 1]

// Helper: Fortran column-major 3D access (1-based) for raw pointer.
// f3(ptr, d1, d2, r, c, k) accesses element (r,c,k) with dims (d1,d2,*).
#define f3(ptr, d1, d2, r, c, k) (ptr)[((k) - 1) * (d1) * (d2) + ((c) - 1) * (d1) + (r) - 1]

// ============================================================================
//
// CALCULATES THE TRANSITION AMPLITUDES F AT ONE angle.
//
// INPUT SCALARS:
//   angle = C.M. angle (DEGREES).
//   nSpline = NUMBER OF AMPLITUDES F.
//   lMn, lMx, lSkp = MINIMUM, MAXIMUM, INCREMENT IN INCOMING L.
//   lxMax = MAXIMUM lx.
//   identicalParticles = 0 IF THERE ARE NO IDENTICAL PARTICLES
//          = 1 IF ONLY ONE PAIR (IN OR OUT) IS IDENTICAL
//          = 2 IF BOTH PAIRS ARE IDENTICAL
//   eta = COULOMB PARAMETER (ELASTIC ONLY).
//   kWave = INCIDENT WAVENUMBER (ELASTIC ONLY).
//   sigZero = L=0 COULOMB PHASE (ELASTIC ONLY).
//   lHigh = MAXIMUM L TO BE INCLUDED IN fHigh.
//
// OUTPUT SCALAR:
//   flopCount = APPROX NUMBER OF FLOATING OPS IN THIS CALL
//
// INPUT ARRAYS:
//   jtocs(4,nSpline)
//   betas(2,nSpline,lMx-lBase+1) — REAL*4
//   aLowFc(2,2)
//
// OUTPUT ARRAYS:
//   F(2,nSpline), fLow(2,nSpline), fHigh(2,nSpline), fError(nSpline), fCoul(2,3)
//   plm(lMx+1+((2*lMx+1-lxMax)*lxMax)/2)
//   contR(lMx+1), contI(lMx+1), fEpsLow(2,LEBACK)
// ============================================================================
// LEBACK > 0 is always true; the epsLon-acceleration branch always runs.
void CrossSectionCalc::ampCalc(double angle, int isElastic, int nSpline, int lMn, int lMx, int lSkp,
            int lxMax, int identicalParticles, double eta, double kWave, double sigZero,
            int returnFLowHigh, int lHigh,
            int* jtocs, float* betas, double* aLowFc,
            double* F, double* fLow, double* fHigh, double* fError,
            double* fCoul, double* plm, double* contR, double* contI,
            double* fEpsLow, int& flopCount)
{
    const double bigNum  = Constants::bigNum;

    // Local variables
    int densitySwitch;
    int densityElasticSwitch;
    double fT[3];    // 1-based: fT[1], fT[2]
    double fts[3];   // 1-based: fts[1], fts[2]
    double ftC[3];   // 1-based: ftC[1], ftC[2]
    double ftSc[3];  // 1-based: ftSc[1], ftSc[2]
    double fOut[3];  // 1-based: fOut[1], fOut[2]

    int loopFlop, parityShift, loMin, loMax, loHigh, lBase;
    int n, k, kOffset, lx, mX, loMinMax, plmBase, lTemp, lo, lParit;
    int riIndex, epsIndex, epsSize;
    double angleRadians, cosAngle, temp, phase, plmFactor, factorMbl;
    double conReal, conImag, dummy;

    //
    //
    flopCount = 0;
    loopFlop = 0;
    parityShift = (std::abs(f2(jtocs, 4, 1, 1)) % 2);
    loMin = (lSkp == 2) ? lMn + parityShift : lMn;
    loMax = lMx;
    if (lSkp == 2) loMax = loMax - ((loMax + parityShift) % (2));
    loHigh = lHigh;
    if (lSkp == 2) loHigh = loHigh - ((loHigh + parityShift) % (2));
    lBase = lMn;
    densitySwitch = (identicalParticles != 0);
    factorMbl = 1;
    if (densitySwitch && !isElastic) factorMbl = std::sqrt(2.0);
    if (identicalParticles == 2 && !isElastic) factorMbl = 2;
    densityElasticSwitch = (densitySwitch && isElastic);
    n = densityElasticSwitch ? 6 : 2 * nSpline;
    for (k = 1; k <= n; k++) {
        f1(F, k) = 0;
        if (returnFLowHigh) {
            f1(fLow, k) = 0;
            f1(fHigh, k) = 0;
        }
    }

    //
    // CALCULATE THE LEGENDRE FUNCTIONS.
    //
    angleRadians = Constants::DEGREE * angle;
    cosAngle = std::cos(angleRadians);
    plmSub(lMx, lxMax, cosAngle, plm-1);
    flopCount = 16 + 7 * (lMx + 1) * (lxMax + 1);

    //
    // ELASTIC: CALCULATE THE COULOMB AMPLITUDE.
    // IDENTICAL PARTICLES NEED IT FOR BOTH angle AND PI-angle.
    //
    if (isElastic) {
    loMin = 0;
    lBase = 0;
    temp = std::max(std::sin(0.5 * angleRadians), 1.0e-10);
    phase = 2 * (sigZero - eta * std::log(temp));
    temp = -0.5 * eta / (kWave * temp * temp);
    if (std::fabs(temp) > bigNum) temp = std::copysign(bigNum, temp);
    ftC[1] = temp * std::cos(phase);
    ftC[2] = temp * std::sin(phase);
    ftSc[1] = 0;
    ftSc[2] = 0;
    flopCount = flopCount + 72;
    if (densitySwitch) {
    angleRadians = Constants::PI - angleRadians;
    temp = std::max(std::sin(0.5 * angleRadians), 1.0e-10);
    phase = 2 * (sigZero - eta * std::log(temp));
    temp = -0.5 * eta / (kWave * temp * temp);
    if (std::fabs(temp) > bigNum) temp = std::copysign(bigNum, temp);
    ftSc[1] = temp * std::cos(phase);
    ftSc[2] = temp * std::sin(phase);
    flopCount = flopCount + 72;
    fts[1] = 0;
    fts[2] = 0;

    //
    // STILL ELASTIC: STORE SYMMETRIC, ANTISYMMETRIC, AND
    // UNSYMMETRIZED COULOMB AMPLITUDES.
    } // end if (densitySwitch)

    for (riIndex = 1; riIndex <= 2; riIndex++) {
        f2(fCoul, 2, riIndex, 1) = ftC[riIndex] + ftSc[riIndex];
        f2(fCoul, 2, riIndex, 2) = ftC[riIndex] - ftSc[riIndex];
        f2(fCoul, 2, riIndex, 3) = ftC[riIndex];
    }
    flopCount = flopCount + 4;
    } // end if (isElastic)

    //
    // LOOP THROUGH kOffset = jT, jProj, lx, mX.
    //
    for (kOffset = 1; kOffset <= nSpline; kOffset++) {
        if (f2(jtocs, 4, 4, kOffset) < 0) continue;
        lx = f2(jtocs, 4, 2, kOffset);
        mX = (f2(jtocs, 4, 1, kOffset) + lx + 1) / 2;
        loMinMax = std::max(loMin, mX);
        if (lSkp == 2) loMinMax = loMinMax + ((loMinMax + parityShift) % (2));
        fT[1] = 0;
        fT[2] = 0;
        plmBase = mX * (2 * lMx + 1 - mX) / 2 + 1;

        //
        // sum OVER lo.
        //
        // fLow = AMPLITUDE FROM lo = loMin, loMin+1.
        // FOR ELASTIC AND lx=0, IT DOES NOT INCLUDE THE 1 IN 1-S(L)
        //
        if (lx == 0 && isElastic && returnFLowHigh) {
        lTemp = loMin + 1;
        for (lo = loMin; lo <= lTemp; lo += lSkp) {
            plmFactor = factorMbl * f1(plm, lo + 1);
            lParit = ((lo) % (2)) + 1;
            f2(fLow, 2, 1, 1) = f2(fLow, 2, 1, 1) + plmFactor * f2(aLowFc, 2, 1, lParit);
            f2(fLow, 2, 2, 1) = f2(fLow, 2, 2, 1) + plmFactor * f2(aLowFc, 2, 2, lParit);
        }
        loopFlop = 5 * (2 + lSkp);
        } // end if (lx == 0 && isElastic && returnFLowHigh)

        //
        // COMPUTE ALL CONTRIBUTIONS FROM NON-ZERO S-MATRICES
        //
        for (lo = loMinMax; lo <= loMax; lo += lSkp) {
            plmFactor = factorMbl * f1(plm, lo + plmBase);
            f1(contR, lo + 1) = plmFactor * f3(betas, 2, nSpline, 1, kOffset, lo - lBase + 1);
            f1(contI, lo + 1) = plmFactor * f3(betas, 2, nSpline, 2, kOffset, lo - lBase + 1);
        }
        loopFlop = loopFlop + 5 * (loMax - loMinMax + lSkp);

        //
        // FIRST sum AMPLITUDE UP TO lMin+1 TO GET fLow
        //
        lTemp = std::min(loMin + 1, loMax);
        if (lTemp >= loMinMax) {
        for (lo = loMinMax; lo <= lTemp; lo += lSkp) {
            fT[1] = fT[1] + f1(contR, lo + 1);
            fT[2] = fT[2] + f1(contI, lo + 1);
        }

        //
        } // end if (lTemp >= loMinMax)

        if (returnFLowHigh && !(isElastic && lx == 0)) {
        f2(fLow, 2, 1, kOffset) = fT[1];
        f2(fLow, 2, 2, kOffset) = fT[2];

        //
        // sum CONTRIBUTION FOR lMin+2 TO lMax (lMax-2 FOR ELASTIC)
        //
        } // end if (returnFLowHigh && ...)

        lTemp = std::max(loMin + 2, loMinMax);
        if (lTemp <= loHigh) {
        for (lo = lTemp; lo <= loHigh; lo += lSkp) {
            fT[1] = fT[1] + f1(contR, lo + 1);
            fT[2] = fT[2] + f1(contI, lo + 1);
        }

        //
        // fHigh = AMPLITUDE FROM lo <= loHigh.
        //
        } // end if (lTemp <= loHigh)

        f2(fHigh, 2, 1, kOffset) = fT[1];
        f2(fHigh, 2, 2, kOffset) = fT[2];

        //
        // NOW THE LAST TERMS > lMax OR lMax-1
        //
        lTemp = loHigh + lSkp;
        if (lTemp <= loMax) {
        for (lo = lTemp; lo <= loMax; lo += lSkp) {
            fT[1] = fT[1] + f1(contR, lo + 1);
            fT[2] = fT[2] + f1(contI, lo + 1);
        }

        //
        // IF NECESSARY APPLY EPSiLoN ALGORITHM
        //
        } // end if (lTemp <= loMax)

        f1(fError, kOffset) = (std::fabs(f1(contR, loMax + 1)) + std::fabs(f1(contI, loMax + 1)))
            / (std::fabs(fT[1]) + std::fabs(fT[2]) + Constants::smlNum);
        constexpr int leBack = 15;
        n = std::min((loMax - loMinMax) / lSkp, leBack - 1);
        if (n > 5) {
        epsSize = n + 1;
        f2(fEpsLow, 2, 1, epsSize) = fT[1];
        f2(fEpsLow, 2, 2, epsSize) = fT[2];
        for (int i = 1; i <= n; i++) {
            epsIndex = epsSize - i;
            f2(fEpsLow, 2, 1, epsIndex) = f2(fEpsLow, 2, 1, epsIndex + 1) - f1(contR, loMax + 1 - lSkp * (i - 1));
            f2(fEpsLow, 2, 2, epsIndex) = f2(fEpsLow, 2, 2, epsIndex + 1) - f1(contI, loMax + 1 - lSkp * (i - 1));
        }
        epsLon(fEpsLow, epsSize, &fT[1], f1(fError, kOffset));

        n = std::min((loHigh - loMinMax) / lSkp, leBack - 1);
        if (n > 5) {
        epsSize = n + 1;
        f2(fEpsLow, 2, 1, epsSize) = f2(fHigh, 2, 1, kOffset);
        f2(fEpsLow, 2, 2, epsSize) = f2(fHigh, 2, 2, kOffset);
        for (int i = 1; i <= n; i++) {
            epsIndex = epsSize - i;
            f2(fEpsLow, 2, 1, epsIndex) = f2(fEpsLow, 2, 1, epsIndex + 1) - f1(contR, loHigh + 1 - lSkp * (i - 1));
            f2(fEpsLow, 2, 2, epsIndex) = f2(fEpsLow, 2, 2, epsIndex + 1) - f1(contI, loHigh + 1 - lSkp * (i - 1));
        }
        epsLon(fEpsLow, epsSize, &fOut[1], dummy);
        f2(fHigh, 2, 1, kOffset) = fOut[1];
        f2(fHigh, 2, 2, kOffset) = fOut[2];

        //
        // FOR ELASTIC SCATTERING OF IDENTICAL PARTICLES,
        // DO IT AGAIN FOR PI - angle.
        //
        } // end if n > 5 (second epsilon)
        } // end if n > 5 (first epsilon)

        if (densityElasticSwitch) {
        for (lo = loMinMax; lo <= loMax; lo += lSkp) {
            conReal = f1(contR, lo + 1);
            conImag = f1(contI, lo + 1);
            if (((lo) % (2)) != 0) {
                conReal = -conReal;
                conImag = -conImag;
            }
            fts[1] = fts[1] + conReal;
            fts[2] = fts[2] + conImag;
            if (returnFLowHigh && lo <= loMin + 1) {
            lParit = ((lo) % (2)) + 1;
            plmFactor = factorMbl * f1(plm, lo + plmBase);
            if (lParit == 2) plmFactor = -plmFactor;
            f2(fLow, 2, 1, 2) = f2(fLow, 2, 1, 2) + plmFactor * f2(aLowFc, 2, 1, lParit);
            f2(fLow, 2, 2, 2) = f2(fLow, 2, 2, 2) + plmFactor * f2(aLowFc, 2, 2, lParit);
            } // end if (returnFLowHigh && lo <= loMin + 1)
            if (returnFLowHigh && lo == loHigh) {
                f2(fHigh, 2, 1, 2) = fts[1];
                f2(fHigh, 2, 2, 2) = fts[2];
            }
        }
        loopFlop = loopFlop + 5 * (loMax - loMinMax + lSkp);

        //
        // END OF lo sum. TRANSFER TO OUTPUT ARRAYS.
        // fLow IS ALL EXCEPT THE FIRST 2 L'S.
        //
        } // end if (densityElasticSwitch)

        f2(F, 2, 1, kOffset) = fT[1];
        f2(F, 2, 2, kOffset) = fT[2];

        if (returnFLowHigh) {
        f2(fLow, 2, 1, kOffset) = fT[1] - f2(fLow, 2, 1, kOffset);
        f2(fLow, 2, 2, kOffset) = fT[2] - f2(fLow, 2, 2, kOffset);
        } // end if (returnFLowHigh)
    } // end kOffset loop

    //
    // END OF kOffset LOOP.
    //
    flopCount = flopCount + loopFlop / lSkp;

    //
    // FOR ELASTIC SCATTERING, ADD THE COULOMB AMPLITUDES.
    //
    if (!isElastic) return;
    for (riIndex = 1; riIndex <= 2; riIndex++) {
        f2(F, 2, riIndex, 1) = f2(F, 2, riIndex, 1) + ftC[riIndex];
        if (returnFLowHigh) {
        f2(fHigh, 2, riIndex, 1) = f2(fHigh, 2, riIndex, 1) + ftC[riIndex];
        f2(fLow, 2, riIndex, 1) = f2(fLow, 2, riIndex, 1) + ftC[riIndex];

        //
        // FOR ELASTIC SCATTERING OF IDENTICAL PARTICLES,
        // FORM SYMMETRIC, ANTISYMMETRIC, AND UNSYMMETRIZED COMBINATIONS.
        //
        } // end if (returnFLowHigh) — fHigh/fLow Coulomb

        if (densityElasticSwitch) {

        fts[riIndex] = fts[riIndex] + ftSc[riIndex];
        f2(F, 2, riIndex, 3) = f2(F, 2, riIndex, 1);
        f2(F, 2, riIndex, 1) = f2(F, 2, riIndex, 3) + fts[riIndex];
        f2(F, 2, riIndex, 2) = f2(F, 2, riIndex, 3) - fts[riIndex];
        if (returnFLowHigh) {
        f2(fLow, 2, riIndex, 2) = fts[riIndex] - f2(fLow, 2, riIndex, 2);
        f2(fLow, 2, riIndex, 3) = f2(fLow, 2, riIndex, 1);
        f2(fLow, 2, riIndex, 1) = f2(fLow, 2, riIndex, 3) + f2(fLow, 2, riIndex, 2);
        f2(fLow, 2, riIndex, 2) = f2(fLow, 2, riIndex, 3) - f2(fLow, 2, riIndex, 2);
        f2(fHigh, 2, riIndex, 2) = f2(fHigh, 2, riIndex, 2) + ftSc[riIndex];
        f2(fHigh, 2, riIndex, 3) = f2(fHigh, 2, riIndex, 1);
        f2(fHigh, 2, riIndex, 1) = f2(fHigh, 2, riIndex, 3) + f2(fHigh, 2, riIndex, 2);
        f2(fHigh, 2, riIndex, 2) = f2(fHigh, 2, riIndex, 3) - f2(fHigh, 2, riIndex, 2);
        } // end if (returnFLowHigh)
        } // end if (densityElasticSwitch)
    }
    flopCount = flopCount + 6;
    return;
}


// ============================================================================
// Merged from angular_distributions.cpp
// ============================================================================
// ============================================================================
//
// CALCULATES THE angle-INDEPENDENT PARTS (BETA) OF THE AMPLITUDES F.
// ============================================================================
void CrossSectionCalc::betCalc(int isElastic, double kWave, int spinProj, int nSpline, int lMn, int lMx, int lSkp,
            int lxMin, int lxMax, int lDeltaMax, int tempsCount, int statsCode, int printSwitch,
            double& sigTotal, double& sigReaction,
            int* jtocs, double* S, float* sMag, float* sPhase, double* sigIn, double* sigOut,
            float* betas, double* totLx, double* totMx, double* aLowFc, double* temps)
{
    // jtocs(4,nSpline), S(2,nSpline,lMx+1), sMag(nSpline,*), sPhase(nSpline,*)
    // sigIn(lMx+1), sigOut(lMx+1)
    // betas(2,nSpline,lMx-lBase+1) — REAL*4
    // totLx(lxMax+1), totMx(nSpline), aLowFc(2,2), temps(tempsCount)

    const double PI = Constants::PI;

    double fOpt[3];   // 1-based: fOpt[1], fOpt[2]
    double reaction[3]; // 1-based: reaction[1], reaction[2]

    //
    //
    int parityShift = ((lDeltaMax) % (2));
    double factor = 0.5 / kWave;
    if (parityShift != 0) factor = -factor;
    int mXCount = lxMax + 1;
    int lxCount = lxMax - lxMin + 1;
    int lBase, loMin;
    if (!isElastic) {
        // REACTION ONLY
        for (int kOffset = 1; kOffset <= nSpline; kOffset++) {
            f1(totMx, kOffset) = 0;
        }
        for (int lx = 0; lx <= lxMax; lx++) {
            f1(totLx, lx + 1) = 0;
        }
        sigTotal = 0;
        lBase = lMn;
        loMin = (lSkp == 2) ? lMn + parityShift : lMn;
    } else {
        // ELASTIC SCATTERING ONLY
        lBase = 0;
        loMin = 0;
        fOpt[1] = 0;
        fOpt[2] = 0;
        reaction[1] = 0;
        reaction[2] = 0;
        for (int i = 1; i <= 4; i++) {
            // aLowFc(i,1) = 0  — linear fill of (2,2) array, elements 1..4
            f1(aLowFc, i) = 0;
        }
    }
    // BOTH ELASTIC AND REACTION: (L40 merge)
    {
        int n = 2 * nSpline * (lMx - lBase + 1);
        for (int i = 1; i <= n; i++) {
            // betas(i,1,1) = 0 — linear fill
            f1(betas, i) = 0;
        }
    }

    // dummy1() Fortran overlay-hint (force clebschGordan COMMON load) deleted

    //
    // THE OUTER LOOP IS OVER lo.
    //
    for (int lo = loMin; lo <= lMx; lo += lSkp) {
        double loDouble = lo;
        double d2L1 = 2 * loDouble + 1.;
        int loIndex = lo - lBase + 1;
        int lParit = ((lo) % (2)) + 1;
        int liMin = lo - lDeltaMax;
        int liMax = lo + lDeltaMax;
        int mXMax = std::min(lxMax, lo);

        //
        // INSERT CLEBSCH-GORDANS ETC. INTO THE temp ARRAY.
        //
        for (int i = 1; i <= tempsCount; i++) {
            f1(temps, i) = 0;
        }
        for (int li = liMin; li <= liMax; li += 2) {
            if (li >= lBase) {
                int lxStart = std::max(std::abs(li - lo), lxMin);
                for (int lx = lxStart; lx <= lxMax; lx++) {
                    int mXz = ((lx + li - lo) % (2));
                    int i = 1 + mXCount * (lx - lxMin + lxCount * (li - liMin) / 2);
                    for (int mX = mXz; mX <= lx; mX++) {
                        f1(temps, i + mX) = factor * (2 * li + 1)
                            * clebschGordan(2 * li, 2 * lx, 0, 2 * mX, 2 * lo, 2 * mX);
                    }
                }
            }
        }

        //
        // START THE MAIN LOOP OVER jT, jProj, lx, li-lo.
        //
        int liPrevious = -100000;
        int mXz = 0, kOffZ = 0;
        for (int kOffset = 1; kOffset <= nSpline; kOffset++) {
            if (f2(jtocs, 4, 4, kOffset) < 0) continue;
            { int lx = f2(jtocs, 4, 2, kOffset);
            int li = lo - f2(jtocs, 4, 1, kOffset);

            if (li >= liPrevious) {
                mXz = lx + f2(jtocs, 4, 1, kOffset);
                kOffZ = kOffset - mXz;
            }
            liPrevious = li;
            if (li < lBase || li > lMx) continue;

            {
                int i = li - lBase + 1;
                double sMatR, sMatI, phase, aMag, cosPhase, sinPhase, sR, sI;

                if (!isElastic) {
                    // Reaction: S from sMag and sPhase
                    aMag = f2(sMag, nSpline, kOffset, i);
                    if (aMag == 0) continue;
                    phase = f2(sPhase, nSpline, kOffset, i) + f1(sigIn, li + 1) + f1(sigOut, lo + 1);
                    sMatR = aMag * std::sin(phase);
                    sMatI = -aMag * std::cos(phase);
                } else {
                    // Elastic: S from S-matrix
                    phase = f1(sigIn, li + 1) + f1(sigIn, lo + 1);
                    cosPhase = std::cos(phase);
                    sinPhase = std::sin(phase);
                    if (lo < lMn) {
                        // Below lMn: S=0, but may need 1-S subtraction for lx=0
                        sMatR = 0;
                        sMatI = 0;
                        if (lx == 0) {
                            // 1-S subtraction
                            sMatR = sMatR - sinPhase;
                            sMatI = sMatI + cosPhase;
                            reaction[lParit] = reaction[lParit] + d2L1;
                        }
                        // lx < 0 or lx > 0: skip to betas with sMatR=sMatI=0
                    } else {
                        // lo >= lMn, have S-matrix
                        sR = f3(S, 2, nSpline, 1, kOffset, i);
                        sI = f3(S, 2, nSpline, 2, kOffset, i);
                        sMatR = sR * sinPhase + sI * cosPhase;
                        sMatI = sI * sinPhase - sR * cosPhase;
                        reaction[lParit] = reaction[lParit] - d2L1 * (sR * sR + sI * sI);
                        if (lx == 0) {
                            if (lo <= lMn + 1) {
                                f2(aLowFc, 2, 1, lParit) = factor * d2L1 * sMatR;
                                f2(aLowFc, 2, 2, lParit) = factor * d2L1 * sMatI;
                            }
                            // 1-S subtraction
                            sMatR = sMatR - sinPhase;
                            sMatI = sMatI + cosPhase;
                            reaction[lParit] = reaction[lParit] + d2L1;
                        }
                        // lx != 0: skip L240, go directly to betas
                    }
                }

                // Multiply S by CG coefficients and sum into betas
                i = 1 + mXCount * (lx - lxMin + lxCount * (li - liMin) / 2);
                for (int mX = mXz; mX <= lx; mX++) {
                    // betas(1,kOffZ+mX,loIndex) and betas(2,kOffZ+mX,loIndex)
                    f3(betas, 2, nSpline, 1, kOffZ + mX, loIndex) =
                        f3(betas, 2, nSpline, 1, kOffZ + mX, loIndex)
                        + (float)(f1(temps, i + mX) * sMatR);
                    f3(betas, 2, nSpline, 2, kOffZ + mX, loIndex) =
                        f3(betas, 2, nSpline, 2, kOffZ + mX, loIndex)
                        + (float)(f1(temps, i + mX) * sMatI);
                }
            }
            } // end scope for lx, li
        } // end kOffset loop

        //
        // NOW CALCULATE THE SQRT(FACTORIALS) INTO temps.
        //
        f1(temps, 1) = 1;
        if (mXMax != 0) {
            for (int mX = 1; mX <= mXMax; mX++) {
                f1(temps, mX + 1) = f1(temps, mX) / std::sqrt((loDouble + mX) * (loDouble - mX + 1));
            }
        }

        //
        // CROSS SECTIONS MUST BE DOUBLED FOR IDENTICAL PARTICLES.
        //
        {
            double temp = 40 * PI / (2 * loDouble + 1);
            if (statsCode != 3) temp = 2 * temp;

            //
            // SECOND PASS THROUGH kOffset LOOP:
            // INCREMENT THE TOTAL AND PARTIAL CROSS SECTIONS,
            // AND PUT THE SQRT(FACTORIALS) INTO BETA.
            //
            for (int kOffset = 1; kOffset <= nSpline; kOffset++) {
                if (f2(jtocs, 4, 4, kOffset) >= 0) {
                int lx = f2(jtocs, 4, 2, kOffset);
                int mX = (f2(jtocs, 4, 1, kOffset) + lx + 1) / 2;
                float betaReal = f3(betas, 2, nSpline, 1, kOffset, loIndex);
                float betaImag = f3(betas, 2, nSpline, 2, kOffset, loIndex);
                if (!isElastic) {
                    // REACTION ONLY.
                    double term = temp * ((double)betaReal * betaReal + (double)betaImag * betaImag);
                    f1(totMx, kOffset) = f1(totMx, kOffset) + term;
                    if (mX > 0) term = 2 * term;
                    sigTotal = sigTotal + term;
                    f1(totLx, lx + 1) = f1(totLx, lx + 1) + term;
                } else if (mX == 0) {
                    // ELASTICS ONLY
                    fOpt[lParit] = fOpt[lParit] + (double)betaImag;
                }
                betaReal = (float)(f1(temps, mX + 1) * betaReal);
                betaImag = (float)(f1(temps, mX + 1) * betaImag);
                f3(betas, 2, nSpline, 1, kOffset, loIndex) = betaReal;
                f3(betas, 2, nSpline, 2, kOffset, loIndex) = betaImag;

                if (printSwitch) std::printf(" JP, JT, LX, MX, LO, B =%3d/2%3d/2%3d%3d%4d%17.5G +%13.5GI\n",
                    f2(jtocs, 4, 3, kOffset), f2(jtocs, 4, 4, kOffset),
                    lx, mX, lo, (double)betaReal, (double)betaImag);
                } // end if (jtocs >= 0)
            }
        }
    }

    //
    // FOR ELASTICS, FINISH CALCULATING THE REACTION AND
    // NUCLEAR CROSS SECTIONS.
    //
    if (!isElastic) return;
    sigReaction = reaction[1] + reaction[2];
    sigTotal = fOpt[1] + fOpt[2];
    if (statsCode != 3) {
        double temp = spinProj;
        double term = temp + 1;
        sigReaction = (temp * sigReaction + 2 * reaction[statsCode]) / term;
        sigTotal = (temp * sigTotal + 2 * fOpt[statsCode]) / term;
    }
    sigReaction = 10 * PI * sigReaction / (kWave * kWave);
    sigTotal = 40 * PI * sigTotal / kWave;
    return;
}















