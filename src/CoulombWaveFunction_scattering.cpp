// CoulombWaveFunction_scattering.cpp — COULST: Coulomb integrals (FF,FG,GF,GG)
// for inelastic DWBA / coupled-channels; integrates pure-Coulomb tails for the
// extrapolation correction.

#include "CoulombWaveFunction.h"
#include "Timing.h"
#include "Reaction.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>

// COULWORK (R8) scratch-buffer length — identical formula in the FF1 and FF2
// setup paths; pure integer buffer sizing (not a physics formula).
static int coulombWorkSize(int nCoulombPoints, int maxCoulomb, int lMaxMax1, int maxDl1) {
    return std::max({2 * nCoulombPoints + 4 * maxCoulomb,
                     2 * (lMaxMax1 + maxDl1) + 1,
                     2 * nCoulombPoints + lMaxMax1 + maxDl1}) + 1;
}

// ============================================================================
// SECTION 4: CoulombWaveFunction::scattering — COULST
//            CoulombWaveFunction::solveRTXLNX — RTXLNX (4-arg)
//            CoulombWaveFunction::generateBasisIndex — GENBNX_full (stub)
//            CoulombWaveFunction::setBasisFactors    — SETBFC_full (stub)
// ============================================================================
//   1. Determines lInMax (the maximum L needed for Coulomb extrapolation)
//      using asymptotic formulae (Alder et al., Rev.Mod.Phys. 28, 432, 1956).
//   2. For non-CC calculations, allocates ICL1FF/FG/GF/GG and integrates
//      INTEGRAL(sumMax -> inf) FF, FG, GF, GG by calling COULIN for each lx.
//   3. Allocates ICL2FF and computes pure-Coulomb (FF-only) integrals out to
//      lOutMax for the Coulomb extrapolation correction.
//   4. For coupled-channels (problemType==24), delegates to GENBNX, SETBFC, SETFG
//      and allocates CC-specific arrays (HOMO, INHR, PADE, SMATITER …).
//   Returns localRc=1 on success, 0 on error.



void CoulombWaveFunction::scattering(int& returnCode, Reaction& reaction) {
    //
    // COULOMB INTEGRALS FOR INELASTIC DWBA & COUPLED CHANNELS
    //

    int debugSwitch;
    int isInfoPrint;
    int hasCoulombCoupling;

    // Local scalars
    int    verbosity;
    int    lx, li, lo, i, k, deltaIndex, liIndex;
    // MBASCP/LBASDF/NBASDF/MBASDF/LCHNDF/MCHNDF/IBINDX_loc dropped earlier
    // same day — only CC dead-block readers.
    int    maxDl1, lDlDimension, lMinMin1, lMaxMax, lMaxMax1;
    std::vector<double> startTableVector;   // COULSTRT: 16*lDlDimension
    std::vector<double> ffWorkVector;   // FFWORK  : 4*dim1 (non-CC) or dim1 (FF2)
    std::vector<double> fiWorkVector;   // FIWORK  : 4*dim2 (non-CC) or dim2 (FF2)
    std::vector<double> foWorkVector;   // FOWORK  : dim2
    std::vector<double> giWorkVector;   // GIWORK  : dim2
    std::vector<double> goWorkVector;   // GOWORK  : dim2
    std::vector<double> coulombWorkVector;   // COULWORK : max(...)
    int    kBase, lDeltaMin, deltaCount, deltaMax;
    int    deltaIndexStart, lDelta, deltaRowIndex, lColIndex, ffBase;
    int    dim1, dim2;
    int    lMinMin;
    int    lExtStart;
    int    negLx;
    int    l;
    int    localRc;

    // double temporaries
    double clTime;
    double lCritDouble;
    double aTerm, bTerm, cTerm;
    double coulombPhaseIn, coulombPhaseOut;
    double* sigma1Pointer = nullptr;
    double* sigma2Pointer = nullptr;

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    int&    printLevel  = reaction.flags.printLevel;
    int&    lMin    = reaction.angMom.lMin;
    int&    lMax    = reaction.angMom.lMax;
    int&    maxLExtrap  = reaction.integrationGrid.maxLExtrap;
    int&    nCoulombPoints  = reaction.integrationGrid.nCoulombPoints;
    constexpr int maxCoulomb = 80;


    int&    lInMax  = reaction.internalState.lInMax;
    int&    iExcit  = reaction.internalState.iExcit;

    double& sumMax  = reaction.integrationGrid.sumMax;
    double& accuracyInel   = reaction.integrationGrid.accuracyInel;

    auto& sigmaArr1 = reaction.distortedWave.channel[1].sigmaArr;
    auto& sigmaArr2 = reaction.distortedWave.channel[2].sigmaArr;
    int&    lOutMax  = reaction.kin.lOutMax;
    double* etaCh     = reaction.kin.etaCh;     // 1-based

    // GRIDCM — using the inelastic-code C++ member names that correspond to
    // the COULST Fortran /gridcmf/ view.  Mapping derived from comparing the
    //
    // COULST Fortran name → C++ GridcmCommon field (same physical slot)
    //   IBINDX → IMSVAL,  MBINDX → IINTS,   NBINDX → INTOFF
    //   NPTMIN  → IWIO,    NUMPT  → nRiRoInterp
    //
    // IBINDX/MBINDX/NBINDX/NMFFAC/MXLXGS/NMBFAC/NPTMIN/IPADE GridData aliases

    int&    lxMin   = reaction.inelastic.lxMin;
    int&    lxMax   = reaction.inelastic.lxMax;
    auto&   cl2ffArr = reaction.inelastic.cl2ffArr;
    int&    nSmatPerL   = reaction.inelastic.nSmatPerL;
    int&    nLx     = reaction.inelastic.nLx;
    int&    nSpl    = reaction.inelastic.nSpl;
    // poolBetas now holds std::vector<double>* into reaction.named.
    // Only [2]=BETARATS is read in this TU (for the betarPointer base pointer).
    std::vector<double>** poolBetas = reaction.inelastic.poolBetas;


    // r2s[1]=real radius, r2s[2]=imag radius, r2s[3]=Coulomb radius, r2s[4]=VC.

    float* times    = reaction.timing.times;    // 1-based [1..8]

    // Decode the (lx, iExcit) channel into the indxs table: derive k, then read
    // kBase/lDeltaMin/deltaCount. Identical preamble open-coded at the FF1 and FF2
    // integral-store sites. Pure integer index arithmetic (no math touched).
    auto decodeChannelK = [&]() {
        int jBaseExcit = (&reaction.inelastic.jpMin)[4 * (iExcit - 1) + 2];
        k = lx + 1 + nLx * (lx - jBaseExcit / 2);
        int* idx = reaction.inelastic.indxsPointer;
        kBase      = idx[3*k - 3];
        lDeltaMin  = idx[3*k - 2];
        deltaCount = idx[3*k - 1];
    };

    // FF/FG/GF/GG work-array index decode (integer addressing only), identical
    auto computeFfBase = [&]() {
        lDelta = lDeltaMin - 2 + 2 * deltaIndex;
        lo   = li + lDelta;
        deltaRowIndex   = ((deltaMax - lDelta) >> 1);
        lColIndex   = lo + li - 2 * lMinMin;
        ffBase   = deltaRowIndex + lColIndex * lDlDimension;
    };

    // Identical COULIN error-bailout guard at both the inelastic and Coulomb
    // call sites; returns true (caller then returns) on a nonzero return code.
    auto bailIfCoulinError = [&]() {
        if (localRc != 0) { std::printf("0**** ERROR RETURN FROM COULIN:%12d\n", localRc); return true; }
        return false;
    };


    // =========================================================================
    // Begin executable code
    // =========================================================================
    localRc = 0;

    verbosity  = printLevel % 10;
    debugSwitch = (verbosity >= 3);
    isInfoPrint = (verbosity >= 2);

    //
    //

    //
    // THESE INTEGRALS WILL BE JUST ZERO IF THE BETA(COULOMB) IS ZERO.
    // IN THIS CASE THE LARGE COULOMB lMax IS ALSO NOT NEEDED.
    //
    lInMax        = lMax;
    times[5]      = 0.0f;
    times[6]      = 0.0f;
    times[7]      = 0.0f;
    times[8]      = 0.0f;
    clTime        = 0.0;
    hasCoulombCoupling        = 0;   // .FALSE.

    // Was: IBINDX/LBINDX/MBINDX setup + LBASCP/LBASDF/LCHNDF reads +
    // generateBasisIndex(1,...) pass-1 call + NOBFAC=1 if no B factors needed.
    //
    // NON-CC CASE: CHECK WHETHER ANY COULOMB BETA IS NONZERO
    // (NCHNDF was set to 2 here; dropped along with NCHN loop unroll below.)
    //
    // betarPointer base = BETARATS vector start, offset by -lxMin/2 so
    // betarPointer[lx/2] == BETARATS[(lx-lxMin)/2]. Was: int LBETAR =
    { double* betarPointer = poolBetas[2]->data() - lxMin / 2;
    for (lx = lxMin; lx <= lxMax; lx += 2) {
        if (betarPointer[lx / 2] != 0.0) {
            hasCoulombCoupling = 1;   // .TRUE.
            break;
        }
    } }


    // -----------------------------------------------------------------------
    // COULOMB COUPLING EXISTS: compute lInMax using asymptotic formula
    // (non-CC: run iff hasCoulombCoupling; CC: run iff hasCoulombCoupling && maxLExtrap != 0)
    // -----------------------------------------------------------------------

    //
    // FIRST WE FIND THE MAXIMUM li TO BE USED FOR THE COULOMB INTEGRALS
    // WE USE AN ASYMPTOTIC FORM (Alder et al., Rev.Mod.Phys. 28, 432, 1956,
    // formula ii E.83).
    // FOR CC THIS IS DONE SEPARATELY FOR EACH CHANNEL AND THE RESULT
    // STORED IN THE CHANNEL DEFINITION ARRAY.
    //

    if (hasCoulombCoupling) {
    if (isInfoPrint) std::printf("-CHN LX  LO-LI  MAXIMUM L NEEDED\n");

    lCritDouble = (double)reaction.kin.lCrit;
    lx = lxMin;

    {
        li = lMax;


        negLx = -lx;
        for (int lDelta = negLx; lDelta <= lx; lDelta += 2) {
            // aTerm = std::fabs( (eta_ch(2)-eta_ch(1))/eta_ch(1) )
            aTerm = std::fabs((etaCh[2] - etaCh[1]) / etaCh[1]);
            bTerm = (lx + lDelta + 1.0) / 2.0;
            cTerm = std::log(reaction.integrationGrid.dwCutoff) - aTerm * (lCritDouble + 0.5) - bTerm * std::log(lCritDouble + 0.5);
            // L = RTXLNX( aTerm, BTERM, CTERM, .1D0 ) + .5
            l = (int)(CoulombWaveFunction::solveRTXLNX(aTerm, bTerm, cTerm, 0.1) + 0.5);
            li = std::max(li, l);

            if (isInfoPrint) std::printf("%4d%3d%6d%13d\n", 2, lx, lDelta, l);
            if (debugSwitch) std::printf("+%49s%15.5G%15.5G%15.5G\n", "", aTerm, bTerm, cTerm);
        }

        //
        // LIMIT BY maxLExtrap
        //
        li = std::min(li, lMax + maxLExtrap);

        lInMax = std::max(lInMax, li);
    }  // end channel block
    }  // end L200 work (Coulomb lInMax)

    // -----------------------------------------------------------------------

    lOutMax = lInMax + lxMax;

    //

    //
    // EXPAND THE COULOMB PHASE SHIFT ARRAYS
    // WE ARE SLOPPY AND FIND SOME MORE THAN ARE REALLY NECESSARY
    //
    {
        int lMaxMaxSigma = lOutMax + 2 * lxMax;
        sigmaArr1.resize(lMaxMaxSigma + 1, 0.0);

        //
        // NON-CC: expand both sigma arrays and fill in new entries
        //
        sigmaArr2.resize(lMaxMaxSigma + 1, 0.0);

        sigma1Pointer = sigmaArr1.data();
        sigma2Pointer = sigmaArr2.data();

        coulombPhaseIn = sigma1Pointer[lMax];
        coulombPhaseOut = sigma2Pointer[lMax];

        lExtStart = lMax + 1;
        for (l = lExtStart; l <= lMaxMaxSigma; l++) {
            coulombPhaseIn = coulombPhaseIn + std::atan(etaCh[1] / (double)l);
            coulombPhaseOut = coulombPhaseOut + std::atan(etaCh[2] / (double)l);
            sigma1Pointer[l] = coulombPhaseIn;
            sigma2Pointer[l] = coulombPhaseOut;
        }

        //
        // WE SETUP INTEGRAL(R TO INF) FF, FG, GF AND GG FOR lMin TO lMax,
        // AND JUST FF FOR (0 TO INF) AND lMin TO lInMax.
        // THIS IS FOR NON-CC CALCULATIONS — FIRST THE (R, INF) STUFF
        //

        reaction.inelastic.cl1ffArr.assign(nSmatPerL, 0.0);
        reaction.inelastic.cl1fgArr.assign(nSmatPerL, 0.0);
        reaction.inelastic.cl1gfArr.assign(nSmatPerL, 0.0);
        reaction.inelastic.cl1ggArr.assign(nSmatPerL, 0.0);
    }

    // -----------------------------------------------------------------------
    // COMMON PATH FOR BOTH CC AND NON-CC: set up work array dimensions
    // -----------------------------------------------------------------------

    {
        // Compute work-array dimensions (valid for both CC and non-CC paths)
        maxDl1  = std::max(2, lxMax);
        lDlDimension  = maxDl1 + 1;
        lMinMin1  = std::max(lMin - maxDl1, 0);
        startTableVector.assign(16 * lDlDimension, 0.0);  // COULSTRT work array

        if (hasCoulombCoupling) {

        // Non-CC branch: work arrays for FF/FG/GF/GG integrals
        lMaxMax   = lMax + lxMax;
        lMaxMax1  = lMax + maxDl1;
        dim1   = lDlDimension * (2 * (lMaxMax1 - lMinMin1) + 1);
        ffWorkVector.assign(4 * dim1, 0.0);  // FFWORK (non-CC, 4*dim1)
        dim2   = std::max(lMaxMax1 + maxDl1, 4 * maxCoulomb) + 1;
        fiWorkVector.assign(4 * dim2, 0.0);  // FIWORK (non-CC, 4*dim2)
        coulombWorkVector.assign(coulombWorkSize(nCoulombPoints, maxCoulomb, lMaxMax1, maxDl1), 0.0);  // COULWORK

        sigma1Pointer = sigmaArr1.data();
        sigma2Pointer = sigmaArr2.data();


        // NON-CC: set up FG, GF, GG, and FO/GI/GO work arrays
        // ffWorkVector[0..dim1): FF; [dim1..2*dim1): FG; [2*dim1..3*dim1): GF; [3*dim1..4*dim1): GG
        // fiWorkVector[0..dim2): FI; [dim2..2*dim2): FO; [2*dim2..3*dim2): GI; [3*dim2..4*dim2): GO
        {

        //
        for (lx = lxMin; lx <= lxMax; lx += 2) {
            lMinMin = std::max(lMin - (lx + 1) / 2, lx / 2);
            CoulombWaveFunction::coulombIntegral(
                lx + 1, lx, lMinMin, lMaxMax,
                etaCh[2], reaction.kin.akOut, sigma2Pointer,
                etaCh[1], reaction.kin.akIn, sigma1Pointer,
                sumMax, 1 /*TRUE*/,
                ffWorkVector.data(),            // FF
                ffWorkVector.data() + dim1,    // FG
                ffWorkVector.data() + 2*dim1,  // GF
                ffWorkVector.data() + 3*dim1,  // GG
                lDlDimension,
                accuracyInel, nCoulombPoints,     // NTERMS dropped
                coulombWorkVector.data(),            // WK
                fiWorkVector.data(),            // FI
                fiWorkVector.data() + dim2,    // FO
                fiWorkVector.data() + 2*dim2,  // GI
                fiWorkVector.data() + 3*dim2,  // GO
                startTableVector.data(),            // ST
                printLevel / 100 % 10, localRc, clTime, reaction);

            if (bailIfCoulinError()) return;

            deltaMax = std::max(2, lx);

            //
            // STORE THE INTEGRALS IN THE SAME ORDER AS THE i(lx,li,lo) ARRAYS.
            // MULTIPLY IN THE CHARGE AND RADIUS FACTORS (BUT NOT THE BETA'S)
            //
            {
                decodeChannelK();  // k/kBase/lDeltaMin/deltaCount from indxsPointer

                for (liIndex = 1; liIndex <= reaction.inelastic.nLValues; liIndex++) {
                    li = reaction.inelastic.lisArr[liIndex];  // 1-based
                    deltaIndexStart = (li < lx) ? lx + 1 - li : 1;

                    for (deltaIndex = deltaIndexStart; deltaIndex <= deltaCount; deltaIndex++) {
                        computeFfBase();  // lDelta/lo/deltaRowIndex/lColIndex/ffBase
                        i    = nSpl * (liIndex - 1) + kBase - 2 + deltaIndex;

                        double r2s4 = reaction.inelastic.r2s[4];  // Coulomb VC
                        // non-CC FF/FG/GF/GG live in reaction.inelastic.CL1*_arr (0-based, idx=i).
                        reaction.inelastic.cl1ffArr[i] = -r2s4 * ffWorkVector[ffBase];
                        reaction.inelastic.cl1fgArr[i] = -r2s4 * ffWorkVector[ffBase + dim1];
                        reaction.inelastic.cl1gfArr[i] = -r2s4 * ffWorkVector[ffBase + 2*dim1];
                        reaction.inelastic.cl1ggArr[i] = -r2s4 * ffWorkVector[ffBase + 3*dim1];
                    }
                }
            }

        }  // end DO 439 lx loop

        }  // end non-CC inner block

        //
        //  490 times(5) = clTime
        //      times(7) = second() - TT
        //
        times[5] = (float)clTime;
        times[7] = 0.0f;  // second() - TT stub

        }  // end if (hasCoulombCoupling)

        {

        // -----------------------------------------------------------------------
        // NOW THE PURE COULOMB (FF) INTEGRALS OUT TO lOutMax
        // -----------------------------------------------------------------------

        // ICL2FF allocation: n_spl*(lInMax+1) elements for li=lMin..lInMax × n_spl channels.
        // Old formula NMLOLX*(lInMax+1) was too small when n_spl > NMLOLX, causing OOB writes
        // that corrupted adjacent pool allocations (like IRDINT) producing NaN.
        dim1   = nSpl * (lInMax - lMin + 1);  // correct size: same as i = nSpl*(li-lMin)+kOffset
        cl2ffArr.assign(dim1, 0.0);

        lMaxMax   = lOutMax;
        lMaxMax1  = lMaxMax + maxDl1 - lxMax;

        if (!hasCoulombCoupling) { localRc = 1; return; }

        dim1   = lDlDimension * (2 * (lMaxMax1 - lMinMin1) + 1);
        ffWorkVector.assign(dim1, 0.0);   // FFWORK (FF2, dim1)
        dim2   = std::max(lMaxMax1 + maxDl1, 4 * maxCoulomb) + 1;
        fiWorkVector.assign(dim2, 0.0);   // FIWORK (FF2, dim2)
        foWorkVector.assign(dim2, 0.0);   // FOWORK (FF2, dim2)
        giWorkVector.assign(dim2, 0.0);   // GIWORK (FF2, dim2)
        goWorkVector.assign(dim2, 0.0);   // GOWORK (FF2, dim2)
        coulombWorkVector.assign(coulombWorkSize(nCoulombPoints, maxCoulomb, lMaxMax1, maxDl1), 0.0);  // COULWORK (FF2)

        {

        //
        for (lx = lxMin; lx <= lxMax; lx += 2) {
            lMinMin = std::max(lMin - (lx + 1) / 2, lx / 2);

            static double dummy8[2] = {0.0, 0.0};

            CoulombWaveFunction::coulombIntegral(
                lx + 1, lx, lMinMin, lMaxMax,
                etaCh[2], reaction.kin.akOut, sigmaArr2.data(),
                etaCh[1], reaction.kin.akIn, sigmaArr1.data(),
                0.0 /*sumMax=0*/, 0 /*FALSE*/,
                ffWorkVector.data(), dummy8, dummy8, dummy8,
                lDlDimension,
                accuracyInel, nCoulombPoints,     // NTERMS dropped
                coulombWorkVector.data(),
                fiWorkVector.data(),
                foWorkVector.data(),
                giWorkVector.data(),
                goWorkVector.data(),
                startTableVector.data(),
                printLevel / 100 % 10, localRc, clTime, reaction);

            if (bailIfCoulinError()) return;

            deltaMax = std::max(2, lx);

            {
                decodeChannelK();  // k/kBase/lDeltaMin/deltaCount from indxsPointer

                for (li = lMin; li <= lInMax; li++) {
                    deltaIndexStart = (li < lx) ? lx + 1 - li : 1;

                    for (deltaIndex = deltaIndexStart; deltaIndex <= deltaCount; deltaIndex++) {
                        computeFfBase();  // lDelta/lo/deltaRowIndex/lColIndex/ffBase
                        i    = nSpl * (li - lMin) + kBase - 2 + deltaIndex;

                        double r2s4 = reaction.inelastic.r2s[4];  // Coulomb VC
                        cl2ffArr[i] = -r2s4 * ffWorkVector[ffBase];
                    }
                }
            }

        }  // end DO 599 lx loop
        }  // end FF2 block

        times[6] = (float)clTime;
        times[8] = 0.0f;  // second() - TT stub

        // Work arrays ffWorkVector/FI/FO/GI/GO/WK/ST freed automatically (vectors)

        localRc = 1;
        return;

        }  // end FF2 path

        // Was ~150 lines: SETBFC call + DSGMAL elastic phase reestablish +
        // SETFG full / FANDG packing + the entire CC allocations-for-COUPLN
        // block (CCHOMO/CCINHR/CCINHI/CCRHSR/CCRHSI/CCSMSQ/PADE/IHOMOA/B/
        // IINHR8/IINHI8). All consumers (CC numerical kernels) were never
        // ported; vectors lived as write-only stubs.
    }  // end L400 block

    localRc = 1;
    return;
}
