// CrossSection_muelco.cpp — MUELCO: angle-independent parts of the Mueller matrix for fixed k.

#include "CrossSectionCalc.h"
#include "math/angular_momentum_coeff.h"
#include <cstdio>
#include <cmath>

// ============================================================================
// Merged from analyzing_powers.cpp
// ============================================================================
// ============================================================================
//
// COEF(qaCount,qbCount,qBigACount,qBigBCount,nSpline,nSpline) — 6D Fortran column-major array,
// passed as double* (0-based flattened). Use macro below.
// COEP(qaCount,qbCount), COET(qBigACount,qBigBCount) — 2D scratch arrays, column-major.
// coEx(kMax+1), termK(kMax+1) — 1D scratch, 0-based index k (Fortran K+1).
// ITOC(4,nSpline) — column-major integer array, passed as int* 0-based.
// ============================================================================

// 6-dimensional COEF access macro: Fortran COEF(qaIndex,qbIndex,qBigAIndex,qBigBIndex,kOff,kOffP)
// dims = (qaCount,qbCount,qBigACount,qBigBCount,nSpline,nSpline), 1-based all indices
#define coef6(qaIndex,qbIndex,qBigAIndex,qBigBIndex,kOff,kOffP) \
    coefficient[ (int)(kOffP-1)*(int)nSpline*(int)qBigBCount*(int)qBigACount*(int)qbCount*(int)qaCount \
        + (int)(kOff-1)*(int)qBigBCount*(int)qBigACount*(int)qbCount*(int)qaCount \
        + (int)(qBigBIndex-1)*(int)qBigACount*(int)qbCount*(int)qaCount \
        + (int)(qBigAIndex-1)*(int)qbCount*(int)qaCount \
        + (int)(qbIndex-1)*(int)qaCount \
        + (int)(qaIndex-1) ]

// ITOC(i,J) with Fortran dims (4,nSpline), passed as 0-based int*
#define tocsf(i,j) tocsPointer[(int)(j-1)*4 + (int)(i-1)]

void CrossSectionCalc::muElCoupling(int kA, int jA, int jB,
            int jResidual, int jBigB, int qaCount,
            int nSpline, int lxParity, int* tocsPointer, double* coefficient, double* coEp,
            double* coEt, double* coEx, double* termK, int debugSwitch)
{
    const int kB = 0, kBigA = 0, kBigB = 0;
    const int qbCount = 1, qBigACount = 1, qBigBCount = 1;

    // Fortran: COEP(qaCount,qbCount), COET(qBigACount,qBigBCount) column-major
    // COEP[col*qaCount + row], COET[col*qBigACount + row] (0-based)
    #define coEpf(qaIndex,qbIndex) coEp[(int)(qbIndex-1)*(int)qaCount + (int)(qaIndex-1)]
    #define coEtf(qBigAIndex,qBigBIndex) coEt[(int)(qBigBIndex-1)*(int)qBigACount + (int)(qBigAIndex-1)]

    int qAs = kA - qaCount + 1;
    int qBs = kB - qbCount + 1;
    int qBigAs = kBigA - qBigACount + 1;
    int qBigBs = kBigB - qBigBCount + 1;
    double parityDouble = (double)lxParity;

    // Clear coefficient (coefCount elements)
    int coefCount = nSpline * nSpline * qaCount * qbCount * qBigACount * qBigBCount;
    if (debugSwitch)
        std::printf("\nMUELCO... %12d%12d%12d%12d%12d%12d%12d%12d%12d%12d%12d%12d%12d%12d%12d\n",
                    kA, kB, kBigA, kBigB, jA, jB, jResidual, jBigB,
                    qaCount, qbCount, qBigACount, qBigBCount, nSpline, lxParity, coefCount);
    for (int i = 0; i < coefCount; i++) coefficient[i] = 0.0;

    // DUMMY2() Fortran overlay-loader hint (factorialTable load) deleted

    // Outer loops over jProj and jT groups (grouped by jProj,jT in ITOC)
    int jTPrevious = -1, jTpPrevious = -1, ktPrevious = -1, jProjPrevious = -1, jPpPrevious = -1, kpPrevious = -1;
    double c9jT = 0.0, c9jP = 0.0;

    int kOfe = 0;
    while (true) {
        // Find next valid kOfs
        int kOfs = kOfe + 1;
        while (true) {
            kOfe = kOfs;
            if (kOfs > nSpline) return;
            if (tocsf(4,kOfs) >= 0) break;
            kOfs++;
        }
        int jT = tocsf(4,kOfs);
        int jProj = tocsf(3,kOfs);
        int lxMin = tocsf(2,kOfs);

        // Find last element with this jProj,jT
        kOfe = kOfs;
        while (true) {
            kOfe++;
            if (kOfe > nSpline) break;
            if (tocsf(4,kOfe) != jT || tocsf(3,kOfe) != jProj) break;
        }
        kOfe--;
        int lxMax = tocsf(2,kOfe);

        // Inner: loop over jProj', jT' groups
        int kOfpe = 0;
        bool outerDone = false;
        while (!outerDone) {
            // Find next valid kOfps
            int kOfps = kOfpe + 1;
            while (true) {
                kOfpe = kOfps;
                if (kOfps > nSpline) { outerDone = true; break; }
                if (tocsf(4,kOfps) >= 0) break;
                kOfps++;
            }
            if (outerDone) break;
            int jTp = tocsf(4,kOfps);
            int jPp = tocsf(3,kOfps);
            int lxpMin = tocsf(2,kOfps);
            kOfpe = kOfps;
            while (true) {
                kOfpe++;
                if (kOfpe > nSpline) break;
                if (tocsf(4,kOfpe) != jTp || tocsf(3,kOfpe) != jPp) break;
            }
            kOfpe--;
            int lxpMax = tocsf(2,kOfpe);

            // kT loop
            int ktMin = std::max((double)std::abs(kBigA-kBigB), (double)std::abs(jT-jTp)/2);
            int ktMax = std::min((double)(kBigA+kBigB), (double)(jT+jTp)/2);
            if (ktMin > ktMax) continue;

            for (int kT = ktMin; kT <= ktMax; kT++) {
                int qtMin = std::max(-(double)kT, (double)(qBigBs-kBigA));
                int qtMax = std::min((double)kT, (double)(kBigB-qBigAs));
                if (qtMin > qtMax) continue;

                // Compute COET if kT/jT/jTp changed
                if (kT != ktPrevious || jTp != jTpPrevious || jT != jTPrevious) {
                    double fac = (double)(jResidual+1)*(jBigB+1)*(2*kBigA+1)*(jT+1)*(jTp+1)*(2*kT+1);
                    c9jT = std::sqrt(fac) * wig9J(jResidual, jT, jBigB,
                                               2*kBigA, 2*kT, 2*kBigB, jResidual, jTp, jBigB);
                    if (debugSwitch)
                        std::printf(" 9J AT 203:%8d%8d%8d%16.6g\n", jT, kT, jTp, c9jT);
                    ktPrevious = kT; jTpPrevious = jTp; jTPrevious = jT;
                }
                if (c9jT == 0.0) continue;

                // Fill COET
                for (int qBigAIndex = 1; qBigAIndex <= qBigACount; qBigAIndex++) {
                    int qBigA = qBigAIndex + qBigAs - 1;
                    for (int qBigBIndex = 1; qBigBIndex <= qBigBCount; qBigBIndex++) {
                        int qBigB = qBigBIndex + qBigBs - 1;
                        coEtf(qBigAIndex,qBigBIndex) = 0.0;
                        if (std::abs(qBigB - qBigA) <= kT)
                            coEtf(qBigAIndex,qBigBIndex) = c9jT * clebschGordan(2*kBigA, 2*kT,
                                2*qBigA, 2*(qBigB-qBigA), 2*kBigB, 2*qBigB);
                    }
                }

                // kP loop
                int kpMin = std::max((double)std::abs(kB-kA), (double)std::abs(jProj-jPp)/2);
                int kpMax = std::min((double)(kB+kA), (double)(jProj+jPp)/2);
                if (kpMin > kpMax) continue;

                for (int kP = kpMin; kP <= kpMax; kP++) {
                    int qpMin = std::max(-(double)kP, (double)(qAs-kB));
                    int qpMax = std::min((double)kP, (double)(kA-qBs));
                    if (qpMin > qpMax) continue;

                    // Compute COEP if kP/jProj/jPp changed
                    if (kP != kpPrevious || jPp != jPpPrevious || jProj != jProjPrevious) {
                        double fac2 = (double)(jA+1)*(jB+1)*(2*kB+1)*(jProj+1)*(jPp+1);
                        c9jP = (double)(2*kP+1) * std::sqrt(fac2)
                             * wig9J(jB, jProj, jA, 2*kB, 2*kP, 2*kA, jB, jPp, jA);
                        if (debugSwitch)
                            std::printf(" 9J AT 253:%8d%8d%8d%16.6g\n", jProj, kP, jPp, c9jP);
                        kpPrevious = kP; jPpPrevious = jPp; jProjPrevious = jProj;
                    }
                    if (c9jP == 0.0) continue;

                    // Fill COEP
                    for (int qaIndex = 1; qaIndex <= qaCount; qaIndex++) {
                        int qA = qaIndex + qAs - 1;
                        for (int qbIndex = 1; qbIndex <= qbCount; qbIndex++) {
                            int qB = qbIndex + qBs - 1;
                            coEpf(qaIndex,qbIndex) = 0.0;
                            if (std::abs(qA - qB) <= kP)
                                coEpf(qaIndex,qbIndex) = c9jP * clebschGordan(2*kB, 2*kP,
                                    2*qB, 2*(qA-qB), 2*kA, 2*qA);
                        }
                    }

                    // lx, lxP loops with kOffset/kOffsetP tracking
                    int kOffset = kOfs;
                    for (int lx = lxMin; lx <= lxMax; lx++) {
                        int mXMin = -lx;
                        int mXz = (tocsf(1,kOffset) + lx + 1) / 2;

                        int kOffsetP = kOfps;
                        for (int lxP = lxpMin; lxP <= lxpMax; lxP++) {
                            int mXpZ = (tocsf(1,kOffsetP) + lxP + 1) / 2;

                            // Compute 3rd 9-J for all k
                            int kMin = std::max(std::max((double)std::abs(lx-lxP), (double)std::abs(kP-kT)),
                                            (double)(qtMin-qpMax));
                            int kMax = std::min((double)(lx+lxP), (double)(kP+kT));
                            if (kMin > kMax) { kOffsetP = kOffsetP + lxP - mXpZ + 1; continue; }

                            {
                                int qMin = std::max(-(double)kMax, (double)(qtMin-qpMax));
                                int qMax = std::min((double)kMax, (double)(qtMax-qpMin));
                                if (qMin > qMax) { kOffsetP = kOffsetP + lxP - mXpZ + 1; continue; }

                                double tLx = (double)(2*lx + 1);
                                for (int k = kMin; k <= kMax; k++) {
                                    coEx[k] = (double)(2*k+1) * std::sqrt(tLx)
                                        * wig9J(2*lx, jProj, jT, 2*k, 2*kP, 2*kT, 2*lxP, jPp, jTp);
                                    if (debugSwitch)
                                        std::printf(" 9J:%5d%5d%5d%5d%5d%5d%5d%5d%5d%16.6g\n",
                                                    kT, kP, lx, lxP, k, jProj, jT, jPp, jTp, coEx[k]);
                                }

                                // mX, mXp loops
                                for (int mX = mXMin; mX <= lx; mX++) {
                                    int kOff = kOffset + std::abs(mX) - mXz;
                                    if (kOff < kOffset) continue;
                                    double fact;
                                    if (mX >= 0) {
                                        fact = 1.0;
                                    } else {
                                        fact = parityDouble;
                                        if ((lx + mX) % 2 != 0) fact = -fact;
                                    }

                                    int mXpMin = std::max(-(double)lxP, (double)(mX+qMin));
                                    int mXpMax = std::min((double)lxP, (double)(mX+qMax));
                                    if (mXpMin > mXpMax) continue;

                                    for (int mXp = mXpMin; mXp <= mXpMax; mXp++) {
                                        int kOffP = kOffsetP + std::abs(mXp) - mXpZ;
                                        if (kOffP < kOffsetP) continue;
                                        double fac = fact;
                                        if (mXp < 0) {
                                            fac *= parityDouble;
                                            if ((lxP + mXp) % 2 != 0) fac = -fac;
                                        }

                                        int q = mXp - mX;

                                        // Compute termK for all k
                                        int kStart = std::max((double)kMin, (double)std::abs(q));
                                        for (int k = kStart; k <= kMax; k++) {
                                            termK[k] = fac * coEx[k]
                                                * clebschGordan(2*lx, 2*k, 2*mX, 2*q, 2*lxP, 2*mXp);
                                        }

                                        // qP loop
                                        int qpMn = std::max((double)qpMin, (double)(qtMin-q));
                                        int qpMx = std::min((double)qpMax, (double)(qtMax-q));
                                        for (int qP = qpMn; qP <= qpMx; qP++) {
                                            int qT = q + qP;

                                            // Sum over k
                                            double sumK = 0.0;
                                            for (int k = kStart; k <= kMax; k++)
                                                sumK += termK[k] * clebschGordan(2*k, 2*kP, 2*q, 2*qP, 2*kT, 2*qT);

                                            // QBIGA, QB loops
                                            int qBigAMin = std::max((double)qBigAs, (double)(qBigBs-qT));
                                            int qBigAMax = std::min((double)kBigA, (double)(kBigB-qT));
                                            int qBMin = std::max((double)qBs, (double)(qAs-qP));
                                            int qBMax = std::min((double)kB, (double)(kA-qP));

                                            for (int qBigA = qBigAMin; qBigA <= qBigAMax; qBigA++) {
                                                int qBigAIndex = qBigA - qBigAs + 1;
                                                int qBigB = qBigA + qT;
                                                int qBigBIndex = qBigB - qBigBs + 1;
                                                double termT = sumK * coEtf(qBigAIndex,qBigBIndex);
                                                if (termT == 0.0) continue;

                                                for (int qB = qBMin; qB <= qBMax; qB++) {
                                                    int qbIndex = qB - qBs + 1;
                                                    int qA = qB + qP;
                                                    int qaIndex = qA - qAs + 1;
                                                    double term = termT * coEpf(qaIndex,qbIndex);
                                                    coef6(qaIndex,qbIndex,qBigAIndex,qBigBIndex,kOff,kOffP)
                                                        += term;
                                                    if (debugSwitch)
                                                        std::printf(
                                                            " COEF:%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%15.5g%15.5g%15.5g\n",
                                                            qB, qBigA, qP, mXp, mX,
                                                            qaIndex, qbIndex, kOffP, kOff, qBigAIndex, qBigBIndex,
                                                            sumK, termT, term);
                                                }
                                            }
                                        } // qP loop
                                    } // mXp loop
                                } // mX loop
                            } // block for qMin/qMax

                            kOffsetP = kOffsetP + lxP - mXpZ + 1;
                        } // lxP loop

                        kOffset = kOffset + lx - mXz + 1;
                    } // lx loop
                } // kP loop
            } // kT loop
        } // inner jProj'/jT' while
    } // outer jProj/jT while

    #undef coEpf
    #undef coEtf
}

#undef tocsf
#undef coef6
