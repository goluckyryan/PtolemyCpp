// CrossSection_eldcs.cpp — ELDCS: differential cross sections for elastic scattering.

#include "CrossSectionCalc.h"
#include "ptolemy_types.h"
#include "MathTables.h"
#include "Reaction.h"
#include "Constants.h"
#include <cstdio>
#include <cmath>
#include <vector>

// Fortran 1-based access helpers (duplicated from CrossSectionCalc.cpp).
#define f1(ptr, i) (ptr)[(i) - 1]
#define f2(ptr, dim1, r, c) (ptr)[((c) - 1) * (dim1) + (r) - 1]
#define f3(ptr, d1, d2, r, c, k) (ptr)[((k) - 1) * (d1) * (d2) + ((c) - 1) * (d1) + (r) - 1]

// Print the ELDCS page header + column titles. Identical block appeared twice
// (initial header at top, reprinted when a page fills); extracted verbatim.
static void printElasticDcsHeader(const Reaction& reaction_, double eLab)
{
    std::printf("1");
    for (int i = 0; i < 47; i++) std::printf(" ");
    std::printf("P T O L E M Y\n");
    std::printf(" ");
    for (int i = 1; i <= 45; i++) std::printf("%c", reaction_.reactStr[i]);
    std::printf("%8.2f MEV     ", eLab);
    for (int i = 1; i <= 65; i++) std::printf("%c", reaction_.header[i]);
    std::printf("\n");
    std::printf("0    ANGLE");
    std::printf("%*sSIGMA/", 11, "");
    std::printf("%*sSIGMA", 19, "");
    std::printf("%*sRUTHERFORD", 19, "");
    std::printf("%*s%% PER", 19, "");
    std::printf("%*s%% PER\n", 9, "");
    std::printf("  C.M.    LAB     RUTHERFORD");
    std::printf("%*sC.M.", 6, "");
    std::printf("%*sLAB", 13, "");
    std::printf("%*sC.M.", 15, "");
    std::printf("%*sLOW L", 16, "");
    std::printf("%*sHIGH L\n", 8, "");
    std::printf("\n");
}

// ============================================================================
//
// CALCULATES DIFFERENTIAL CROSS SECTIONS FOR ELASTIC SCATTERING
// ============================================================================
void CrossSectionCalc::elasticDcs(double kWave, double eta, double angMinArg, double angMaxArg, double angStepArg,
           int lMn, int lMx, int lSkp, int statsCode, int spinProj,
           const double* smatBase, const int* tocsBase, int nSpline, const double* sigBase, double eLab, double tau,
           int printSwitch, std::vector<double>& torutOut,
           int keepFAmplitude, std::vector<double>& F_out, double& sigReaction,
           std::vector<double>& ruthOut, std::vector<double>& crossOut)
{
    const auto RADIAN = Constants::RADIAN;
    const auto DEGREE = Constants::DEGREE;
    auto& outputInLab = reaction_.flags.outputInLab;
    constexpr int leBack = 15;
    // function now use the caller-supplied F_out / sigmaArr pointers.
    auto& aBar  = reaction_.kin.aBar;

    int betPrintSwitch;
    int isIdenticalParticles;

    double fCoul[7];     // flat 1-based: Fortran fCoul(2,3), access via f2(fCoul,2,r,c)
    double aLowFc[5]; // flat 1-based: Fortran aLowFc(2,2), 4 elements
    float dummy4[2];     // DUM4(1,1) scratch
    double dummy8[2];    // DUM8(1) scratch

    double angleStep = angStepArg;

    betPrintSwitch = (((reaction_.flags.printLevel) % (10)) >= 4);

    //
    // PRINT header
    //
    if (printSwitch) {
        printElasticDcsHeader(reaction_, eLab);
    }

    double angleMin = angMinArg;
    double angleMax = angMaxArg;
    if (printSwitch) aBar = 180;

    //
    // FOR LAB ANGLES WE HAVE FUNNY STUFF IF TAU > 1
    // FOR XSECTN CALLS, IT HAS ALL BEEN WORKED OUT ALREADY
    //
    if (outputInLab != 0 && printSwitch && tau > 1) {
        //
        // OUR CONVENTION IS THAT NEGATIVE VALUES FOR ANGLEMIN OR ANGLEMAX
        // CORRESPOND TO THE 2ND BRANCH.
        //
        double angBar = RADIAN * std::asin(1.0 / tau);
        angleStep = std::fabs(angleStep);
        aBar = angleStep * (std::trunc((angBar - angleMin) / angleStep + 1.e-10)
            + angleMin);
    }
    // L50:
    {
        double twoABar = 2 * aBar;
        angleMax = std::min(angleMax, twoABar);
        angleStep = std::copysign(angleStep, angleMax - angleMin);
        int angleCount = (int)((angleMax - angleMin) / angleStep + 1.5);

        //
        // MAKE A SPACE EVERY ONCE AND A WHILE.
        //
        double angleBlock = 52;
        if (angleCount > 50) {
            if (angleStep < .5)     angleBlock = 1;
            else if (angleStep < 1) angleBlock = 5;
            else                    angleBlock = 10;
        }
        int lineCount = 1;  // L80:
        int numBlocks = (int)(angleBlock / angleStep + .5);
        int lineMax = 53;

        if (!(nSpline == 1 || statsCode == 3)) {
            std::printf("0**** CANNOT EVALUATE ELASTIC SCATTERING FOR "
                        "IDENTICAL PARTICLES WITH S.O. FORCE\n");
            return;
        }
        //
        // ALLOCATE SPACE FOR THETA-INDEPENDENT ARRAYS
        //
        int lxMax = (nSpline > 1) ? spinProj : 0;
        int lDeltaMax = lxMax - ((lxMax) % (2));
        int tempsCount = (lxMax + 1) * (lxMax + 1) * (lDeltaMax + 1);
        // tempsVector: local BETCAL scratch (angle-independent temp storage)
        std::vector<double> tempsVector(tempsCount, 0.0);
        // betasVector: local BETCAL/AMPCAL float scratch — not exported
        std::vector<float> betasVector(2 * nSpline * (lMx + 1), 0.0f);

        setLog(2 * (lMx + lxMax));

        //
        // CALCULATE "BETAELAS", THE ANGLE-INDEPENDENT PARTS OF F.
        // ALSO CALCULATE TOTAL REACTION AND NUCLEAR CROSS SECTIONS.
        //
        double sigNuc;
        // tocsBase is 0-based (tocsBase[0] = first int); betCalc expects ITOC as int*
        // where ITOC[0] = first int of slot (0-based), so pass tocsBase directly.
        betCalc(TRUE_F, kWave, spinProj, nSpline, lMn, lMx, lSkp,
            0, lxMax, lDeltaMax, tempsCount, statsCode, betPrintSwitch,
            sigNuc, sigReaction,
            const_cast<int*>(tocsBase), const_cast<double*>(smatBase), dummy4, dummy4, const_cast<double*>(sigBase),
            dummy8, betasVector.data(), dummy8, dummy8, aLowFc, tempsVector.data());

        // tempsVector and betasVector auto-destruct (no manual NFREE needed)

        //
        // ALLOCATE SPACE FOR LEGENDRE FUNCTIONS, AMPLITUDES, CROSS SECTIONS.
        //
        int plmSize = lMx + 1 + ((2 * lMx + 1 - lxMax) * lxMax) / 2;
        std::vector<double> plmVector(plmSize, 0.0);
        std::vector<double> contRIVector(2 * (lMx + 1), 0.0);
        isIdenticalParticles = (statsCode != 3);
        int identicalParticles = isIdenticalParticles ? 2 : 0;
        int fSize = isIdenticalParticles ? 6 : 2 * nSpline;
        std::vector<double> fTempVector(fSize, 0.0);
        std::vector<double> fLowVector(fSize, 0.0);
        std::vector<double> fHighVector(fSize, 0.0);
        std::vector<double> fErrVector(fSize, 0.0);
        std::vector<double> fEpsLowVector(2 * leBack + 1, 0.0);
        // F_out: sized fSize*angleCount when keepFAmplitude, written via 1-based fOutBase[2*angIndex + 2*angleCount*K] alias below.
        //        F_out feeds analyzingPower() — caller passes the same vector through.
        // torutOut/ruthOut/crossOut: sized angleCount; written via torutBase/ruthBase/crossBase aliases below.
        if (keepFAmplitude) F_out.assign(fSize * angleCount, 0.0);
        torutOut.assign(angleCount, 0.0);
        ruthOut.assign(angleCount, 0.0);
        crossOut.assign(angleCount, 0.0);

        //
        // ALLOCATIONS ARE DONE, GET THE BASE ADDRESSES.
        //
        // tocsBase is 0-based (tocsBase[0] = first int of slot)
        const int* tocsPointer = tocsBase;  // tocsPointer[4*K+j-1] addresses slot K's j-th int
        // fOutBase into caller-owned F_out.
        // Original Fortran offset preserved: fOutBase[2*angIndex + 2*angleCount*K] addresses F_out[2*angIndex + 2*angleCount*(K-1) - 2].
        double* fOutBase = keepFAmplitude ? F_out.data() - 2 * angleCount - 2 : nullptr;
        // 0-based pointers for cross section arrays (angIndex runs 1..angleCount; accessed [angIndex-1])
        double* crossBase = crossOut.data();        // crossBase[angIndex-1] = crossOut[angIndex-1]
        double* torutBase = torutOut.data();        // torutBase[angIndex-1] = torutOut[angIndex-1]
        double* ruthBase  = ruthOut.data();         // ruthBase[angIndex-1]  = ruthOut[angIndex-1]

        double statFactor = 0;
        if (statsCode < 3) statFactor = (statsCode - 1.50e0) * spinProj / (spinProj + 1.e0);
        int lHigh = lMx - 2 - ((lMx) % (lSkp));
        double sigZero = sigBase[0];  // first element of sigma array

        //
        // LOOP OVER ALL ANGLES
        //
        for (int angIndex = 1; angIndex <= angleCount; angIndex++) {

            double angle = angleMin + (angIndex - 1) * angleStep;

            if (outputInLab != 0) {
                double an = DEGREE * angle;
                if (angle > aBar) an = DEGREE * (twoABar - angle);
                double term2 = 1 - (tau * std::sin(an)) * (tau * std::sin(an));
                if (term2 < 0) term2 = 0;
                term2 = std::cos(an) * std::sqrt(term2);
                if (angle > aBar) term2 = -term2;
                angle = RADIAN * std::acos(-tau * (std::sin(an)) * (std::sin(an)) + term2);
            }

            {
                double cosAngle = std::cos(angle * DEGREE);
                double angleLab = RADIAN * std::atan2(std::sin(angle * DEGREE), cosAngle + tau);

                double term2 = 1 + tau * (tau + 2 * cosAngle);
                double temp = 1 + tau * cosAngle;
                double aJacob = 0;
                if (term2 < 0) term2 = 0;
                if (temp != 0) aJacob = std::fabs(term2 * std::sqrt(term2) / temp);

                //
                // CALCULATE THE AMPLITUDES FOR THIS ANGLE.
                //
                int flopCount = 0;
                ampCalc(angle, TRUE_F, nSpline, lMn, lMx, lSkp,
                       lxMax, identicalParticles, eta, kWave, sigZero, printSwitch, lHigh,
                       const_cast<int*>(tocsBase), betasVector.data(), aLowFc,
                       fTempVector.data(), fLowVector.data(), fHighVector.data(),
                       fErrVector.data(), fCoul,
                       plmVector.data(), contRIVector.data(), contRIVector.data() + lMx + 1,
                       fEpsLowVector.data(), flopCount);

                //
                // CALCULATE THE CROSS SECTIONS.
                //
                double sigma = 0;
                double sigLow = 0;
                double sigHigh = 0;
                double ruth;
                if (!isIdenticalParticles) {
                    for (int kOffset = 1; kOffset <= nSpline; kOffset++) {
                        if (tocsPointer[4 * kOffset - 1] >= 0) {  // (0-based)
                        int mX = (tocsPointer[4 * kOffset - 4] + tocsPointer[4 * kOffset - 3] + 1) / 2;  // (0-based)
                        temp = (mX != 0) ? 20. : 10.;
                        sigma = sigma + temp *
                            (fTempVector[2*(kOffset-1)] * fTempVector[2*(kOffset-1)]
                           + fTempVector[2*(kOffset-1)+1] * fTempVector[2*(kOffset-1)+1]);
                        if (printSwitch) {
                            sigLow = sigLow + temp *
                                (fLowVector[2*(kOffset-1)] * fLowVector[2*(kOffset-1)]
                               + fLowVector[2*(kOffset-1)+1] * fLowVector[2*(kOffset-1)+1]);
                            sigHigh = sigHigh + temp *
                                (fHighVector[2*(kOffset-1)] * fHighVector[2*(kOffset-1)]
                               + fHighVector[2*(kOffset-1)+1] * fHighVector[2*(kOffset-1)+1]);
                        }
                        } // end scope
                    }
                    ruth = (f2(fCoul, 2, 1, 3) * f2(fCoul, 2, 1, 3)
                          + f2(fCoul, 2, 2, 3) * f2(fCoul, 2, 2, 3)) * 10.;
                } else {
                //
                // IDENTICAL PARTICLES (NO SPIN-DEPENDENT FORCES).
                //
                // IDENTICAL PARTICLES: fCoul(r,c) via f2(fCoul,2,r,c)
                ruth = (f2(fCoul, 2, 1, statsCode) * f2(fCoul, 2, 1, statsCode)
                    + f2(fCoul, 2, 2, statsCode) * f2(fCoul, 2, 2, statsCode)
                    + statFactor * (f2(fCoul, 2, 1, 1) * f2(fCoul, 2, 1, 1)
                    + f2(fCoul, 2, 2, 1) * f2(fCoul, 2, 2, 1)
                    - f2(fCoul, 2, 1, 2) * f2(fCoul, 2, 1, 2)
                    - f2(fCoul, 2, 2, 2) * f2(fCoul, 2, 2, 2))) * 10.;
                sigma = (fTempVector[2*(statsCode-1)] * fTempVector[2*(statsCode-1)]
                    + fTempVector[2*(statsCode-1)+1] * fTempVector[2*(statsCode-1)+1]
                    + statFactor * (fTempVector[0]*fTempVector[0] + fTempVector[1]*fTempVector[1]
                    - fTempVector[2]*fTempVector[2] - fTempVector[3]*fTempVector[3])) * 10.;
                if (printSwitch) {
                sigLow = (fLowVector[2*(statsCode-1)] * fLowVector[2*(statsCode-1)]
                    + fLowVector[2*(statsCode-1)+1] * fLowVector[2*(statsCode-1)+1]
                    + statFactor * (fLowVector[0]*fLowVector[0] + fLowVector[1]*fLowVector[1]
                    - fLowVector[2]*fLowVector[2] - fLowVector[3]*fLowVector[3])) * 10.;
                sigHigh = (fHighVector[2*(statsCode-1)] * fHighVector[2*(statsCode-1)]
                    + fHighVector[2*(statsCode-1)+1] * fHighVector[2*(statsCode-1)+1]
                    + statFactor * (fHighVector[0]*fHighVector[0] + fHighVector[1]*fHighVector[1]
                    - fHighVector[2]*fHighVector[2] - fHighVector[3]*fHighVector[3])) * 10.;
                } // end if (printSwitch)
                } // end else isIdenticalParticles (L250 path)

                //
                // SAVE F IF NECESSARY. (L260 merge)
                //
                if (keepFAmplitude) {
                {
                    //
                    // fSize/2 = nSpline EXCEPT =3 FOR IDENTICAL PARTICLES
                    //
                    int fPairCount = fSize / 2;
                    for (int K = 1; K <= fPairCount; K++) {
                        fOutBase[2 * angIndex + 2 * angleCount * K] = fTempVector[2*(K-1)];
                        fOutBase[2 * angIndex + 2 * angleCount * K + 1] = fTempVector[2*(K-1)+1];
                    }
                }
                } // end if (keepFAmplitude)

                //
                // COMPUTE RELATIVE ERRORS FROM LOWER AND UPPER L-TRUNCATIONS
                //
                {
                    double errLow = 50 * (1 - sigLow / sigma);
                    double errHigh = 50 * (1 - sigHigh / sigma);
                    double sigmaReal = 0;
                    if (ruth != 0) sigmaReal = sigma / ruth;

                    if (printSwitch) {
                    { double sigmaLab = aJacob * sigma;

                    //
                    // PRINT EVERYTHING IN MILLIBARNS
                    //
                    lineCount = lineCount + 1;
                    if (!(lineCount <= lineMax || angIndex == angleCount)) {
                    // Reprint header
                    printElasticDcsHeader(reaction_, eLab);
                    lineCount = 1;
                    } // end "if need header"
                    std::printf("%7.2f%7.2f%14.4G%14.4G%12.4G%15.4G%13.2f%8.2f\n",
                        angle, angleLab, sigmaReal, sigma, sigmaLab, ruth, errLow, errHigh);

                    if (std::fmod(angle + angleStep + 1.e-7, angleBlock) <= 1.e-5) {
                        std::printf(" \n");
                        lineCount = lineCount + 1;
                        if (lineCount + numBlocks > lineMax) lineCount = lineMax + 1;
                    }
                    } // end scope for sigmaLab
                    } // end if (printSwitch)
                    torutBase[angIndex - 1] = sigmaReal;
                    ruthBase[angIndex - 1]  = ruth;
                    if (outputInLab != 0 && !printSwitch)
                        sigma = aJacob * sigma;
                    crossBase[angIndex - 1] = sigma;
                }
            }
        }

        if (printSwitch) {
            std::printf("\n0TOTAL REACTION CROSS SECTION =%13.3f MB\n"
                        "  NUCLEAR TOTAL CROSS SECTION =%13.3f MB\n", sigReaction, sigNuc);
        }

        //
        // FREE ALLOCATED SPACE
        //
        // plmVector, betasVector, fTempVector, fLowVector, fHighVector, contRIVector, fErrVector, fEpsLowVector auto-destruct here.
        // Caller-owned torutOut/ruthOut/crossOut remain populated for downstream readers
        // (xSection elastic-channel sigma reads / TORUTH reads; Elastic.cpp passes throwaway).
    }
    return;
}
