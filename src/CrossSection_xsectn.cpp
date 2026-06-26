// CrossSection_xsectn.cpp — XSECTN: transfer/inelastic cross sections from S-matrix elements.

#include "CrossSectionCalc.h"
#include "ptolemy_types.h"
#include "MathTables.h"
#include "Timing.h"
#include "Reaction.h"
#include "Constants.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

// Fortran Gw.d format emulation.
// F-mode (0.1 <= |x| < 10^d): F(w-4).(d-n) right-justified + 4 trailing spaces.
//   n = floor(log10(|x|))+1 for |x|>=1, 0 for 0.1<=|x|<1, 1 for x=0 (gfortran behavior).
// E-mode otherwise: [sign]0.dddddE±ee right-justified in field w.
static std::string fG(double x, int w, int d) {
    double absVal = std::fabs(x);
    double limit = std::pow(10.0, (double)d);
    bool fMode;
    int n;
    if (x == 0.0) {
        fMode = true; n = 1;
    } else if (absVal >= 0.1 && absVal < limit) {
        fMode = true;
        n = (int)std::floor(std::log10(absVal)) + 1;
        if (n < 0) n = 0;
    } else {
        fMode = false; n = 0;
    }
    char buf[128];
    if (fMode) {
        int fw = w - 4, fd = d - n;
        if (fd < 0) fd = 0;
        if (fd == 0) {
            // Fortran Fw.0 prints a trailing decimal point; C's %.0f does not
            char num[64];
            std::snprintf(num, sizeof(num), "%.0f.", x);
            int numLen = (int)std::strlen(num);
            int pad = fw - numLen; if (pad < 0) pad = 0;
            return std::string(pad, ' ') + num + "    ";
        }
        std::snprintf(buf, sizeof(buf), "%*.*f", fw, fd, x);
        return std::string(buf) + "    ";
    }
    // E-mode: Fortran notation 0.ddddd x 10^expDecimal
    int expDecimal = (absVal > 0.0) ? (int)std::floor(std::log10(absVal)) + 1 : 0;
    double mantissa = (absVal > 0.0) ? absVal / std::pow(10.0, (double)expDecimal) : 0.0;
    if (mantissa >= 1.0)      { expDecimal++; mantissa /= 10.0; }
    if (mantissa > 0.0 && mantissa < 0.1) { expDecimal--; mantissa *= 10.0; }
    char mantStr[64];
    std::snprintf(mantStr, sizeof(mantStr), "%.*f", d, mantissa);
    if (mantStr[0] == '1') { expDecimal++; mantissa /= 10.0; std::snprintf(mantStr, sizeof(mantStr), "%.*f", d, mantissa); }
    char signChar = (x < 0.0) ? '-' : ' ';
    char expSignChar = (expDecimal >= 0)  ? '+' : '-';
    char content[64];
    std::snprintf(content, sizeof(content), "%c%sE%c%02d",
                  signChar, mantStr, expSignChar, std::abs(expDecimal));
    int cLen = (int)std::strlen(content);
    int pad = w - cLen; if (pad < 0) pad = 0;
    return std::string(pad, ' ') + content;
}

// start and again on every page break; extracted verbatim. The caller still owns
// the post-header `lineCount = 4; headerPrintSwitch = true;` bookkeeping.
static void printXsectnHeader(const Reaction& reaction_, double eLab)
{
    std::printf("1%54sP T O L E M Y\n"
                "%10sCOMPUTATION OF CROSS SECTIONS"
                "%38sSOPHISTICATED PROGRAMS CONTAIN SUBTLE ERRORS\n"
                "0%.45s%7.2f MEV     %.65s\n\n",
                "", "", "",
                &reaction_.reactStr[1], eLab, &reaction_.header[1]);
}

void CrossSectionCalc::xSection()
{
    const auto RADIAN = Constants::RADIAN;
    const auto smlNum = Constants::smlNum;
    auto& lMin   = reaction_.angMom.lMin;
    auto& lMax   = reaction_.angMom.lMax;
    constexpr int leBack = 15;
    auto& angleMin = reaction_.rxn.angleMin;
    auto& angleMax = reaction_.rxn.angleMax;
    auto& angleStep = reaction_.rxn.angleStep;
    auto& eLab   = reaction_.energies.eLab;
    auto& eCm    = reaction_.energies.eCm;
    auto& mass_a    = reaction_.masses.massesArr[0];
    auto& mass_b    = reaction_.masses.massesArr[1];
    auto& mass_A = reaction_.masses.massesArr[2];
    auto& mass_B = reaction_.masses.massesArr[3];
    auto& Q      = reaction_.energies.Q;
    auto& outputInLab = reaction_.flags.outputInLab;
    auto& lInMax = reaction_.internalState.lInMax;
    auto& akIn    = reaction_.kin.akIn;
    auto& akOut    = reaction_.kin.akOut;
    auto& lOutMax = reaction_.kin.lOutMax;
    auto& lxMax  = reaction_.inelastic.lxMax;
    auto& lxStep = reaction_.inelastic.lxStep;
    auto& nSpl  = reaction_.inelastic.nSpl;
    // CHARACTER*8 data
    // FENAME dropped — F-amplitude buffers (FIN/FOUT) now caller-owned std::vector
    // live entries were INHEAD[1]="INCOMING" and OUTHED[1]="OUTGOING",
    char8 labWord("LAB.    ");
    char8 cmWord("C.M.    ");
    char8 BLANK("        ");
    char8 RUTWRD("/RUTH   ");

    // Local variables
    double tStart, t1, t2, t3, t4, t5, t6, t7, t8;
    double tBeta, tCross, flopTotal;
    bool isTransferReaction, debugSwitch, convergencePrintSwitch, headerPrintSwitch, keepFAmplitude;
    int verbosity, angleCount, lineCount, lineMax = 57;
    double ruthFactor, tau;
    auto& aBar = reaction_.kin.aBar;
    // TORUT/RUTH/sigma/F per-channel storage from elasticDcs (caller-owned vectors).
    // crossSectionInVector/crossSectionOutVector replace the pool slots formerly NAMLOC'd as "CROSSSEC"
    // and renamed via NAMCOM.NAMES[...] = "SIGMAIN "/"SIGMAOUT" (the rename was dead).
    std::vector<double> torutInVector, torutOutVector, ruthInVector;
    std::vector<double> crossSectionInVector, crossSectionOutVector;
    std::vector<double> fInVector, fOutVector, fXferVector;
    double reactionIn, reactionOut, sigTotal;
    int identicalParticles;
    char8 labCm;
    int lxScanMin, lxScanMax, lDeltaMax, jpMax, lxRange, mlxCount, tempsCount;
    int plmSize;
    // LTORT2_l/LSIGIN no longer needed — readers go through 1-based vector aliases.
    double angleBlock;
    int numBlocks, linesPerBlock, mlxHeaderCount;

    tStart = (float)second();
    // discarded iret after the call, and xSection only ever wrote returnCode=1 at the end.
    flopTotal = 0;
    isTransferReaction = (reaction_.internalState.stripPickup != 0);
    verbosity = reaction_.flags.printLevel % 10;
    debugSwitch = (verbosity >= 4);
    convergencePrintSwitch = (verbosity >= 3);

    // Print phase shifts (allocates sMag/sPhase as side-effect for XSECTN)
    phasePrint(verbosity >= 1);
    t1 = (float)second();
    tBeta = 0;
    tCross = 0;

    // Print header
    printXsectnHeader(reaction_, eLab);
    lineCount = 4;

    headerPrintSwitch = true;

    labCm = (outputInLab != 0) ? labWord : cmWord;



    // Lab angle range limits
    aBar = 1.0e+5;
    if (outputInLab != 0) {

    {
        double sinThetaLab = mass_A / mass_a;
        double temp = sqrt(sinThetaLab * mass_B * (mass_a + mass_A) / (mass_b * (mass_b + mass_B))
                         * (eCm + Q) / eCm);
        if (debugSwitch) std::printf("\nLIMITS ON SIN(THETALAB):%15.5G%15.5G\n", sinThetaLab, temp);
        sinThetaLab = std::min(sinThetaLab, temp);
        if (sinThetaLab <= 1) {
        double angBar = RADIAN * asin(sinThetaLab);
        aBar = angleStep * ((int)((angBar - angleMin) / angleStep + 1.0e-10) + 0.5) + angleMin;
        double twoABar = 2 * aBar;
        angleMax = std::min(angleMax, twoABar);
        }
    }
    } // end outputInLab != 0
    angleCount = (int)((angleMax - angleMin) / angleStep + 1.5);

    // Compute elastic cross sections
    keepFAmplitude = reaction_.distortedWave.channel[1].hasSpinorbit;
    elasticDcs(akIn, reaction_.kin.etaCh[1], angleMin, angleMax, angleStep,
          lMin, lInMax, reaction_.distortedWave.channel[1].lSkips, reaction_.distortedWave.channel[1].statsCode, reaction_.distortedWave.channel[1].twoSpin,
          reaction_.distortedWave.channel[1].smatArr.data(), reaction_.distortedWave.channel[1].toceArr.data(),
          reaction_.distortedWave.channel[1].nJStates, reaction_.distortedWave.channel[1].sigmaArr.data(),
          eLab, reaction_.kin.tauRatio[1], 0 /*false*/,
          torutInVector, keepFAmplitude ? 1 : 0, fInVector,
          reactionIn, ruthInVector, crossSectionInVector);

    // 2nd elastic channel — IPROB collapsed to literal 2 after the CC IPROB=1
    // override went (problemType>=23 unreachable). AKS_p was a 1-vs-2 dispatch that
    // now always picks akOut.
    keepFAmplitude = reaction_.distortedWave.channel[2].hasSpinorbit;
    {
        // ruthOutDiscardVector — elasticDcs always writes ruthOut, but no caller
        // reads the second-channel result. Discarded after the call.
        std::vector<double> ruthOutDiscardVector;
        elasticDcs(akOut, reaction_.kin.etaCh[2], angleMin, angleMax, angleStep,
              lMin, lInMax, reaction_.distortedWave.channel[2].lSkips, reaction_.distortedWave.channel[2].statsCode,
              reaction_.distortedWave.channel[2].twoSpin, reaction_.distortedWave.channel[2].smatArr.data(),
              reaction_.distortedWave.channel[2].toceArr.data(),
              reaction_.distortedWave.channel[2].nJStates, reaction_.distortedWave.channel[2].sigmaArr.data(),
              eLab, reaction_.kin.tauRatio[2],
              0 /*false*/,
              torutOutVector, keepFAmplitude ? 1 : 0, fOutVector,
              reactionOut, ruthOutDiscardVector, crossSectionOutVector);
    }

    t2 = (float)second();

    // Setup for calculations. NCHNDF/MCHN/LTOCOF/LSOFF/LFOFF dropped

    // Rutherford conversion factor
    ruthFactor = 1;
    if (reaction_.kin.etaCh[1] != 0)
        ruthFactor = reaction_.kin.etaCh[2] * akIn / (reaction_.kin.etaCh[1] * akOut);


    identicalParticles = reaction_.inelastic.densitySwitch ? 1 : 0;

    {

        // Scan TOCS for extrema
        lxScanMin = 1000;
        lxScanMax = -1000;
        lDeltaMax = -1000;
        jpMax = -1000;
        {
        int* tocsPointer = reaction_.inelastic.tocsPointer;  // 1-based + per-channel offset
        for (int kOffset = 1; kOffset <= nSpl; kOffset++) {
            if (tocsPointer[4 * kOffset] < 0) continue;
            int lx = tocsPointer[4 * kOffset - 2];
            lxScanMin = std::min(lxScanMin, lx);
            lxScanMax = std::max(lxScanMax, lx);
            int lDelta = tocsPointer[4 * kOffset - 3];
            lDeltaMax = std::max(lDeltaMax, lDelta);
            jpMax = std::max(jpMax, tocsPointer[4 * kOffset - 1]);
        }
        }
        lxRange = lxScanMax - lxScanMin + 1;
        mlxCount = (lxScanMax - lxScanMin + lxStep) / lxStep;
        if (debugSwitch) std::printf(" LXMN, LXMX, LDELMX, JPMX, NMLX:%8d%8d%8d%8d%8d\n",
                                lxScanMin, lxScanMax, lDeltaMax, jpMax, mlxCount);

        tempsCount = (lxScanMax + 1) * lxRange * (lDeltaMax + 1);
        std::vector<double> tempsVector(tempsCount + 1, 0.0);
        std::vector<double> lxTotalVector(lxScanMax + 2, 0.0);
        std::vector<double> mxTotalVector(nSpl + 1, 0.0);

        // Compute BETA's
        t3 = (float)second();
        setLog(2 * (lOutMax + lxMax));

        // sMag/sPhase live in reaction_.inelastic.smagArr / sphaseArr (0-based).
        float* smagBasePointer   = reaction_.inelastic.smagArr.data();
        float* sphaseBasePointer = reaction_.inelastic.sphaseArr.data();
        // non-CC reads reaction_.inelastic.sMatrixArr (0-based).
        float* sBasePointer = reaction_.inelastic.sMatrixArr.data();
        int* tocsPointer = reaction_.inelastic.tocsPointer;  // 1-based + per-channel offset


        for (int qq = 0; qq < nSpl * (lInMax + 1); qq++) {
            if (std::isnan((double)smagBasePointer[qq])) smagBasePointer[qq] = 0.0f;
            if (std::isnan((double)sphaseBasePointer[qq])) sphaseBasePointer[qq] = 0.0f;
        }
        double dummyD = 0;
        double dummy8[6] = {};  // aLowFc(2,2)=4 elems, fCoul(2,3)=6 elems
        betCalc(0 /*false*/, akIn, reaction_.distortedWave.channel[1].twoSpin, nSpl, lMin, lInMax,
               reaction_.inelastic.lSkip,
               lxScanMin, lxScanMax, lDeltaMax, tempsCount, reaction_.distortedWave.channel[1].statsCode, debugSwitch ? 1 : 0,
               sigTotal, dummyD,
               tocsPointer + 1, dummy8, smagBasePointer, sphaseBasePointer,
               reaction_.distortedWave.channel[1].sigmaArr.data(), reaction_.distortedWave.channel[2].sigmaArr.data(),
               sBasePointer, lxTotalVector.data(), mxTotalVector.data(),
               dummy8, tempsVector.data());

        t4 = (float)second();
        tBeta = tBeta + t4 - t3;

        // Compute tau for lab/CM conversion
        tau = sqrt(mass_a * mass_b * (mass_b + mass_B) * eCm /
                    (mass_A * mass_B * (mass_a + mass_A) * (eCm + Q)));

        plmSize = lInMax + 1 + ((2 * lInMax + 1 - lxScanMax) * lxScanMax) / 2;
        std::vector<double> plmVector(plmSize + 1, 0.0);
        std::vector<double> contRIVector(2 * (lInMax + 1) + 1, 0.0);

        // V_LXCRO/P_LXCRO (formerly LXCROSSS/MXCROSSS pool slots) and V_MXCON
        keepFAmplitude = std::min({(int)reaction_.angMom.js[1], jpMax, lxScanMax}) > 0;
        if (keepFAmplitude) fXferVector.assign(2 * angleCount * nSpl, 0.0);
        std::vector<double> lxContribVector(lxScanMax + 2, 0.0);
        std::vector<double> fTempVector(2 * nSpl + 1, 0.0);
        std::vector<double> fLowVector(2 * nSpl + 1, 0.0);   // init 0 (replaces INIT8)
        std::vector<double> fHighVector(2 * nSpl + 1, 0.0);
        std::vector<double> fErrVector(nSpl + 1, 0.0);
        std::vector<double> fEpsLowVector(2 * leBack + 2, 0.0);

        setLog(2 * (lOutMax + lxMax));

        // Get locations
        // non-CC reads reaction_.inelastic.sMatrixArr (0-based).
        float* sBasePointer2 = reaction_.inelastic.sMatrixArr.data();
        // fXferPointer aliases fXferVector with the same 1-based offset as the
        double* fXferPointer = fXferVector.empty() ? nullptr : fXferVector.data() - 2 * angleCount - 2;
        // TORT/RUTH/sigma(elastic) read via 1-based vector aliases owned above.
        // When zProj*zTarget==0 (neutral projectile/target) the TORT2 fallback uses crossSectionOutVector
        double* torutInPointer = torutInVector.data();     // torutInPointer[angIndex-1] = torutInVector[angIndex-1]
        double* ruthInPointer = ruthInVector.data();
        double* crossSectionInPointer = crossSectionInVector.data();
        double* torutOutPointer = (reaction_.charges.zProj * reaction_.charges.zTarget == 0)
            ? crossSectionOutVector.data()
            : torutOutVector.data();
        int* tocsPointer2 = reaction_.inelastic.tocsPointer;  // 1-based + per-channel offset

        // Setup print spacing
        angleBlock = 300;
        if (angleCount > 40) {
            if (angleStep <= 0.201)     angleBlock = 1;
            else if (angleStep < 0.99)  angleBlock = 5;
            else                        angleBlock = 10;
        }

        mlxHeaderCount = std::min(mlxCount, 3);
        const char* ruthLabel = (reaction_.charges.zProj * reaction_.charges.zTarget == 0) ? "     " : "/RUTH";

        // linesPerBlock: number of subheader lines per page (depends on mlxCount)
        linesPerBlock = (mlxCount + 2) / 3;

        // Loop over angles
        numBlocks = (int)(angleBlock / angleStep + 0.5);
        numBlocks = numBlocks * linesPerBlock;

        for (int angIndex = 1; angIndex <= angleCount; angIndex++) {
            if (lineCount >= lineMax) {
                printXsectnHeader(reaction_, eLab);
                lineCount = 4;
                headerPrintSwitch = true;
            }
            if (headerPrintSwitch) {
                std::printf("0  %.4s  REACTION     REACTION   LOW L  HIGH L   %% FROM\n", labCm.data);
                std::printf("+%58s%-8.8s ELASTIC%6s%-8.8s ",
                            "", "INCOMING", "", "OUTGOING");
                for (int i = 0; i < mlxHeaderCount; i++) {
                    if (i < mlxHeaderCount - 1)
                        std::printf("   REACTION   ");
                    else
                        std::printf("   REACTION");
                }
                std::printf("\n");
                std::printf("  ANGLE  %.4s MB        /RUTH     %%/L   %% ERROR  L>LMAX\n", labCm.data);
                std::printf("+%58s%.4s MB      /RUTH      %.5s\n",
                            "", labCm.data, ruthLabel);
                if (mlxCount > 1) {
                    bool fl = true;
                    int count = 0;
                    for (int lx = lxScanMin; lx <= lxScanMax; lx += lxStep) {
                        if (count % 3 == 0) {
                            if (fl) { std::printf("+%89s", ""); fl = false; }
                            else      std::printf("\n %89s", "");
                        }
                        std::printf("    LX =%2d", lx);
                        if (count % 3 < 2 && (lx + lxStep) <= lxScanMax)
                            std::printf("    ");
                        count++;
                    }
                    std::printf("\n");
                }
                // FORMAT 548: blank line
                std::printf(" \n");
                lineCount += linesPerBlock + 1;
                headerPrintSwitch = false;
            }
            // process angle
            double angle = angleMin + (angIndex - 1) * angleStep;
            double angleCm = angle;
            double aJacob = 1.0;

            if (outputInLab != 0) {
                // Lab to CM conversion
                double angleOrig = angle;
                if (angle > aBar) angle = 2 * aBar - angle;
                double angleRadians = Constants::DEGREE * angle;
                double temp = 1 - (tau * sin(angleRadians)) * (tau * sin(angleRadians));
                if (temp < 0) temp = 0;
                temp = cos(angleRadians) * sqrt(temp);
                if (angleOrig > aBar) temp = -temp;
                angleCm = acos(temp - tau * sin(angleRadians) * sin(angleRadians));
                temp = tau * (tau + 2 * cos(angleCm)) + 1;
                if (temp < 0) temp = 0;
                aJacob = temp * sqrt(temp) / fabs(tau * cos(angleCm) + 1);
                angleCm = RADIAN * angleCm;
            }

            // Calculate transition amplitudes
            int lHigh = lMax;
            int flopCount;
            ampCalc(angleCm, 0 /*false*/, nSpl, lMin, lInMax, reaction_.inelastic.lSkip,
                   lxScanMax, identicalParticles, 0.0, 0.0, 0.0, 1 /*true*/, lHigh,
                   tocsPointer2 + 1, sBasePointer2, dummy8,
                   fTempVector.data(), fLowVector.data(), fHighVector.data(),
                   fErrVector.data(), dummy8,
                   plmVector.data(), contRIVector.data(), contRIVector.data() + lInMax + 1,
                   fEpsLowVector.data(), flopCount);
            flopTotal += flopCount;

            double sigma = 0;
            double sumLow = 0, sumHigh = 0, sumErr = 0;

            for (int lx = 0; lx <= lxScanMax; lx++)
                lxContribVector[lx] = 0;

            // Sum over kOffset = jT, jProj, lx, mX
            for (int kOffset = 1; kOffset <= nSpl; kOffset++) {
                if (tocsPointer2[4 * kOffset] < 0) continue;
                int lx = tocsPointer2[4 * kOffset - 2];
                int mX = (tocsPointer2[4 * kOffset - 3] + lx + 1) / 2;

                double mPairMult = (mX != 0) ? 2.0 : 1.0;

                double fReal = fTempVector[2 * (kOffset - 1)];
                double fImag = fTempVector[2 * (kOffset - 1) + 1];
                if (keepFAmplitude && fXferPointer != nullptr) {
                    fXferPointer[2 * angIndex + 2 * angleCount * kOffset]     = fReal;
                    fXferPointer[2 * angIndex + 2 * angleCount * kOffset + 1] = fImag;
                }

                double contrib = (10 * aJacob) * (fReal * fReal + fImag * fImag);
                contrib = contrib * mPairMult;
                lxContribVector[lx] = contrib + lxContribVector[lx];
                sigma += contrib;

                sumLow += (10 * mPairMult * aJacob) *
                    (fLowVector[2 * (kOffset - 1)] * fLowVector[2 * (kOffset - 1)]
                   + fLowVector[2 * (kOffset - 1) + 1] * fLowVector[2 * (kOffset - 1) + 1]);
                sumHigh += (10 * mPairMult * aJacob) *
                    (fHighVector[2 * (kOffset - 1)] * fHighVector[2 * (kOffset - 1)]
                   + fHighVector[2 * (kOffset - 1) + 1] * fHighVector[2 * (kOffset - 1) + 1]);
                sumErr += fErrVector[kOffset - 1] * contrib;

                // P_LXCRO[angIndex + angleCount*lx] = lxContribVector[lx] write dropped
            }


            // Compute errors
            double safe = std::max(sigma, smlNum);
            double errLow = 100 * (sigma - sumLow) / (safe + sumLow);
            double hiCont = 100 * (sigma - sumHigh) / safe;
            sumErr = 200 * sumErr / safe;

            // Print cross sections
            double ruth = ruthInPointer[angIndex - 1];
            double elas = crossSectionInPointer[angIndex - 1];
            double sigmaReal = 0;
            if (ruth != 0) sigmaReal = sigma / (ruthFactor * aJacob * ruth);

            std::printf("%7.2f%s%10.6f%8.2f%8.2f%8.1f%s%10.6f%12.6f",
                        angle, fG(sigma, 13, 5).c_str(),
                        sigmaReal, errLow, sumErr, hiCont,
                        fG(elas, 14, 5).c_str(),
                        torutInPointer[angIndex - 1], torutOutPointer[angIndex - 1]);
            if (mlxCount > 1) {
                // (T91, 3G14.5): first group at col 91 (already there), overflow → new line
                int lxPrintIndex = 0;
                for (int lx = lxScanMin; lx <= lxScanMax; lx += lxStep) {
                    if (lxPrintIndex > 0 && lxPrintIndex % 3 == 0)
                        std::printf("\n %89s", "");  // ' ' CC + 89 spaces → T91
                    std::printf("%s", fG(lxContribVector[lx], 14, 5).c_str());
                    lxPrintIndex++;
                }
            }
            std::printf("\n");

            lineCount += linesPerBlock;
            if (convergencePrintSwitch) lineCount += nSpl;
            double origAngle = angleMin + (angIndex - 1) * angleStep;
            if (fmod(origAngle + angleStep + 1e-7, angleBlock) <= 1e-5) {
                if (lineCount < lineMax) std::printf(" \n");
                lineCount++;
                if (lineCount + numBlocks > lineMax) lineCount = 1000;
            }
        } // end angle loop

        t5 = (float)second();
        tCross += t5 - t4;
        flopTotal += 26 * angleCount;

        std::printf("\n0TOTAL:%s\n", fG(sigTotal, 15, 5).c_str());
        double safe = std::max(sigTotal, smlNum);
        std::printf("+%53s%s%s\n", "", fG(reactionIn, 14, 5).c_str(), fG(reactionOut, 24, 5).c_str());
        if (mlxCount > 1) {
            int count = 0;
            for (int lx = lxScanMin; lx <= lxScanMax; lx += lxStep) {
                if (count % 3 == 0) {
                    if (count == 0) std::printf("+%89s", "");  // '+' + T91
                    else            std::printf("\n %89s", ""); // ' ' CC + T91
                }
                std::printf("%s", fG(lxTotalVector[lx], 14, 5).c_str());
                count++;
            }
            std::printf("\n");
        }

        if (mlxCount > 1) {
            for (int lx = lxScanMin; lx <= lxScanMax; lx += lxStep)
                lxTotalVector[lx] = 100.0 * lxTotalVector[lx] / safe;
            int count = 0;
            for (int lx = lxScanMin; lx <= lxScanMax; lx += lxStep) {
                if (count % 3 == 0) {
                    if (count > 0) std::printf("\n");
                    std::printf(" %85s", "");
                }
                std::printf("%13.2f%%", lxTotalVector[lx]);
                count++;
            }
            std::printf("\n");
        }

        if (nSpl > 1) {
            std::printf("\n0 JP  JT  LX MX      TOTAL      PERCENT"
                        "%5s(VALUES FOR M > 0 NOT DOUBLED.)\n\n", "");
            for (int kOffset = 1; kOffset <= nSpl; kOffset++) {
                int jTt = tocsPointer2[4 * kOffset];
                if (jTt < 0) continue;
                int lx = tocsPointer2[4 * kOffset - 2];
                int mX = (tocsPointer2[4 * kOffset - 3] + lx + 1) / 2;
                double percent = 100 * mxTotalVector[kOffset - 1] / safe;
                std::printf("%3d/2%3d/2%3d%3d%s%8.2f\n",
                            tocsPointer2[4 * kOffset - 1], jTt,
                            lx, mX, fG(mxTotalVector[kOffset - 1], 15, 5).c_str(), percent);
            }
        }

    } // end main channel block

    // In non-CC mode IRDINT is the 0 sentinel (betas data lives in reaction_.inelastic.sMatrixArr) and
    // the rename would have written into NAMES[0], an unused 1-based-array slot.
    // plmVector, contRIVector, V_MXCON, lxContribVector, fTempVector, fLowVector, fHighVector, fErrVector, fEpsLowVector,

    // Analyzing powers
    t6 = (float)second();
    t7 = t6;

    if (reaction_.angMom.js[1] != 0) {
    // Reaction analyzing powers
    if (lxScanMax == 0 || jpMax == 0) {
        std::printf("\n---- ALL REACTION ANALYZING POWERS ARE ZERO ----\n");
    } else {
        int lxParity = reaction_.angMom.parities[1] * reaction_.angMom.parities[2] * reaction_.angMom.parities[3] * reaction_.angMom.parities[4];
        if (lxParity == 0 && isTransferReaction)
            lxParity = 1 - 2 * ((reaction_.boundState.vertex[1].lBound + reaction_.boundState.vertex[2].lBound) % 2);
        if (lxParity == 0) {
            std::printf("\n**** REACTION ANALYZING POWERS CANNOT BE CALCULATED"
                        " BECAUSE THE PARITY CHANGE IS UNKNOWN ****\n");
        } else {
            analyzingPower(angleMin, angleMax, angleStep, (int)reaction_.angMom.js[1], (int)reaction_.angMom.js[2],
                   (int)reaction_.angMom.js[3], (int)reaction_.angMom.js[4], nSpl, lxParity,
                   debugSwitch ? 1 : 0, eLab, "REACTION",
                   fXferVector,
                   reaction_.inelastic.tocsPointer);  // 1-based class-owned
            t7 = (float)second();
        }
    }

    // Incoming elastic analyzing powers
    if (!reaction_.distortedWave.channel[1].hasSpinorbit) {
        std::printf("\n---- ALL INCOMING ELASTIC ANALYZING POWERS ARE ZERO ----\n");
    } else {
        analyzingPower(angleMin, angleMax, angleStep, (int)reaction_.angMom.js[1], (int)reaction_.angMom.js[1],
               (int)reaction_.angMom.js[3], (int)reaction_.angMom.js[3], reaction_.distortedWave.channel[1].nJStates, 1,
               debugSwitch ? 1 : 0, eLab, "INCOMING",
               fInVector, reaction_.distortedWave.channel[1].toceArr.data() - 1);
    }
    } // end js[1] != 0

    // Outgoing elastic analyzing powers
    if (reaction_.angMom.js[2] != 0) {
        if (reaction_.distortedWave.channel[2].hasSpinorbit) {
            double eLabOut = (1 + reaction_.kin.tauRatio[2]) * (reaction_.distortedWave.channel[2].Ecm + Q);
            analyzingPower(angleMin, angleMax, angleStep, (int)reaction_.angMom.js[2], (int)reaction_.angMom.js[2],
                   (int)reaction_.angMom.js[4], (int)reaction_.angMom.js[4], reaction_.distortedWave.channel[2].nJStates, 1,
                   debugSwitch ? 1 : 0, eLabOut, "OUTGOING",
                   fOutVector, reaction_.distortedWave.channel[2].toceArr.data() - 1);
        } else {
            std::printf("\n---- ALL OUTGOING ELASTIC ANALYZING POWERS ARE ZERO ----\n");
        }
    }
    t8 = (float)second();

    // Print timing
    {
        double totalTime = second() - tStart;
        double tPhase = t1 - tStart;
        double tElastic = t2 - t1;
        t6 = t7 - t6;
        t7 = t8 - t7;
        double t9 = totalTime - tPhase - tElastic - tBeta - tCross - t6 - t7;
        flopTotal = 1.0e-6 * flopTotal;
        double flopsPerSec = (tCross > 0) ? flopTotal / tCross : 0;
        std::printf("0\n0\n0CPU TIMES FOR CROSS SECTIONS (SECONDS):\n"
                    "0PRINTING PHASE SHIFTS%7s%10.3f\n"
                    " ELASTIC CROSS SECTIONS%6s%10.3f\n"
                    " BETA'S (SUM ON LI)%10s%10.3f\n"
                    " REACTION CROSS SECTIONS%5s%10.3f%10s(%7.2f MFLOPS =%7.3f MFLOP/SEC )\n"
                    " REACTION ANALYZING POWERS%3s%10.3f\n"
                    " ELASTIC ANALYZING POWERS%4s%10.3f\n"
                    " ALL OTHER TIME%14s%10.3f\n"
                    " TOTAL TIME%18s%17.3f\n\n0\n",
                    "", tPhase, "", tElastic, "", tBeta,
                    "", tCross, "", flopTotal, flopsPerSec,
                    "", t6, "", t7, "", t9, "", totalTime);
    }

}
