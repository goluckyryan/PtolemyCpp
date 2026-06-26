// inelastic_dwba.cpp — INELDC: transfer-reaction radial integrals (DWBA transfer
// amplitudes). Computes H(li,lo,lx) = integral of A12, V, and the bound-state
// wavefunctions, integrates H over the distorted waves, then renormalizes to
// S(li,lo,lx,jProj,jT).

#include "MathTables.h"
#include "Reaction_potentials.h"
#include "ScatteringSolver.h"
#include "Timing.h"
#include "angular_coupling.h"
#include "math/spline.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "Reaction.h"
#include "Constants.h"

// A12 (angularCoupling12) — declared in angular_coupling.h.

// ineld2 — the INELDC epilogue: frees work arrays and prints timing statistics.
// Sole caller is DWBAGrid::inelDc below, so it lives here as a file-static.
static void ineld2(float tStart, int hInterpCount, int secondCallCount, int innerLoopPassCount, int cosineIterCount,
            int m1m2LoopPassCount, int phiLoopPassCount, int riRoProcessedCount, int hComputedCount,
            Reaction& reaction)
{
    auto& nPhiPoints   = reaction.gridData.nPhiPoints;

    // print timing info
    //   1 = WAVELJ initialization, 2 = WAVELJ loop, 3 = WAVELJ interp/norm
    //   4 = A12, 5 = phi/A12/first-U loops, 6 = H-interpolations
    //   7 = complete V-loop (5+6), 8 = complete lo-loop (1+2+3+4+7)
    float* times = reaction.timing.times;  // 1-based: times[1..8]

    double totalTime = (double)second() - (double)tStart;
    double otherTime = totalTime - (double)times[8];
    times[2] = times[2] + times[1];
    times[8] = times[8] - times[7] - times[4] - times[3] - times[2];
    times[7] = times[7] - times[6] - times[5];

    std::printf("0%-33s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7d CALLS\n",
                "TIME IN A12",           (double)times[4],
                " IN PHI & A12 LOOP",   (double)times[5],
                " INTERPOLATING H'S",   (double)times[6],
                " COMPUTING SCATTERING WAVES",     (double)times[2],
                " INTERPOLATING SCATTERING WAVES", (double)times[3],
                " REST OF RI & RO LOOP",           (double)times[7],
                " REST OF LO LOOP",                (double)times[8],
                " ALL OTHER TIME",                 otherTime,
                " TOTAL TIME",                     totalTime,
                " NUMBER OF CALLS TO second",      secondCallCount);

    innerLoopPassCount = innerLoopPassCount * nPhiPoints;
    cosineIterCount = cosineIterCount * nPhiPoints;
    // times[7] = average H time (milliseconds per H computed)
    times[7] = (hComputedCount > 0) ? 1.e3f * (times[4] + times[5] + times[7]) / hComputedCount : 0.f;
    // times[4] = average term time (microseconds)
    times[4] = 1.e6f * times[5] / (float)(innerLoopPassCount + 1);
    // times[5] = average cosine time (microseconds)
    times[5] = 1.e6f * times[5] / (float)(cosineIterCount + 1);
    // times[6] = average H-interpolation time (microseconds)
    times[6] = (hInterpCount > 0) ? 1.e6f * times[6] / (float)hInterpCount : 0.f;
    m1m2LoopPassCount = nPhiPoints * ((m1m2LoopPassCount + 4 * phiLoopPassCount) / 5);
    phiLoopPassCount = phiLoopPassCount * nPhiPoints;
    double skipPercent = (double)reaction.inelastic.nLValues * (double)reaction.gridData.nPhiSum * (double)reaction.gridData.nPhiDifference;
    skipPercent = (skipPercent > 0) ? 100.0 * (skipPercent - riRoProcessedCount) / skipPercent : 0.0;

    // FORMAT 2587: blank line
    std::printf(" \n");

    std::printf("%12d PASSES WERE MADE THROUGH THE INNERMOST LOOP;       "
                " AVERAGE term TIME =%9.3f MICROSECONDS\n",
                innerLoopPassCount, (double)times[4]);
    // Col 63: I12(12) + label(41) = 53; T63 fills 9 spaces to col 63; ' AVERAGE' space at 63, A at 64
    std::printf("%12d ITERATIONS OF COS(ARG-MU*PHI) WERE MADE;         "
                " AVERAGE COSINE TIME =%9.3f MICROSECONDS\n",
                cosineIterCount, (double)times[5]);
    std::printf("%12d PASSES WERE MADE THROUGH THE (M1, M2) LOOP\n"
                "%12d PASSES WERE MADE THE PHI LOOP\n", m1m2LoopPassCount, phiLoopPassCount);

    std::printf("%12d H'S WERE COMPUTED OR READ                              "
                "AVERAGE H TIME =%9.3f MILLISECONDS\n",
                hComputedCount, (double)times[7]);
    // Col 55: I12(12) + label(32) = 44; T55 fills 10 spaces; "AVERAGE" (no lead space) at col 55
    std::printf("%12d H'S WERE FOUND BY INTERPOLATION          "
                "AVERAGE H-INTERPOLATION TIME =%9.3f MICROSECONDS\n",
                hInterpCount, (double)times[6]);
    std::printf("%12.2f%% OF THE (RI,RO) POINTS WERE SKIPPED DUE TO DWCUTOFF.\n\n\n",
                skipPercent);
}

void DWBAGrid::inelDc(Reaction& reaction) {
    //
    // COMPUTES THE TRANSFER REACTION RADIAL INTEGRALS.
    //
    // Local variables
    int debugSwitch;
    int isInfoPrint;
    int isAfterFirstLi;

    float tStart, tt, tt2, tt5, tt7;

    double halfNeg = -0.50;

    int&    printLevel  = reaction.flags.printLevel;
    int&    lMin    = reaction.angMom.lMin;
    int&    lMax    = reaction.angMom.lMax;
    int&    nPhiSum   = reaction.gridData.nPhiSum;
    int&    nPhiDifference   = reaction.gridData.nPhiDifference;
    int&    nPhiPoints   = reaction.gridData.nPhiPoints;
    constexpr int nCosin = 256;


    int&    lxMax   = reaction.inelastic.lxMax;
    int&    numLx   = reaction.inelastic.numLx;
    int&    nLValues    = reaction.inelastic.nLValues;

    int&    nWfo    = reaction.gridData.nWfo;
    int&    nRiRoInterp  = reaction.gridData.nRiRoInterp;
    int&    nInterpPoints  = reaction.gridData.nInterpPoints;
    int&    cosinQuarter  = reaction.gridData.cosinQuarter;

    // Local scalars
    int verbosity;
    // INIA12 dropped — A12's INIT parameter went; sole caller (here) was
    // setting it to 5 to enable the setup block, which always runs now.
    int liMin, liL, liParity;
    int liIndex, li;
    int trkWriteOffset, kWO, kWI;
    int lo, lasO, jPo, jPi;
    // LASI dropped — was only used as dispatchPartialWave's LAS arg
    int i, ii, srcIndex, vIndex, uIndex, hIndex;
    int hBlockTag;
    int riRoIndex, smivlIndex;
    int phiPointIndex;
    int hBlockSize;
    int jA12M, jA12N, jA12An;
    int hCount, loMinMin, loMaxMax;
    // lA12Of/lhStart init to 0 silences -Wmaybe-uninitialized: both are
    // assigned inside the jA12M / kA12m inner loop and read at the loop-end
    // sanity printf (`i = lhStart + lA12Of - la12Vl - 1`). gcc can't
    // prove the inner loop iterates at least once.
    int lA12Of = 0, lhStart = 0, hStride, muCount;
    int kA12m, mu, lh;
    int kDw, dwCount, numIi;
    int secondCallCount, cosineIterCount, innerLoopPassCount;
    int m1m2LoopPassCount, phiLoopPassCount, riRoProcessedCount, hComputedCount;
    int hInterpCount;
    int cosIndex;
    int jA12Of, kA12mEnd;

    double factor;
    double dwAmpMin, dwAmpMax, dwRioC, dwAmpCount, dwAmpPeak, dwAmpThis;
    double phi, phiP, phiT, trapWeight;
    double arg, deltaC, cosIndexReal;
    double sin2Phi, cos2Phi, sin2PhiSq;
    double cosTheta, sinSin, cosNew;
    double fiR, fiI, foR, foI;   // foR to avoid keyword clash
    double dwReal, dwImag;
    double riRoWeight;
    double pi256;

    // Allocator pointer locals (computed once at start).
    int la12Vl;
    int lhEnd;

    // =========================================================================
    // Begin executable code
    // =========================================================================

    tStart = (float)second();
    // iret after the call, and inelDc only ever wrote returnCode=1 at the end.

    for (int timeIndex = 1; timeIndex <= 8; timeIndex++) {
        reaction.timing.times[timeIndex] = 0.0f;
    }

    secondCallCount = 4;
    cosineIterCount = 0;
    innerLoopPassCount = 0;
    m1m2LoopPassCount = 0;
    phiLoopPassCount = 0;
    riRoProcessedCount = 0;
    hComputedCount = 0;
    hInterpCount = 0;

    verbosity = printLevel % 10;
    debugSwitch = (verbosity >= 4);
    isInfoPrint = (verbosity >= 2);

    //
    // WRITE header
    //   10X, 'COMPUTATION OF TRANSFER S-MATRIX ELEMENTS', T95, 'WHEELS WITHIN WHEELS' /
    //   '0', 45A1, 'eLab =', F7.2, ' MEV', 5X, 65A1 / )
    //
    std::printf("1%58sP T O L E M Y\n", "");
    std::printf("%10sCOMPUTATION OF TRANSFER S-MATRIX ELEMENTS%43sWHEELS WITHIN WHEELS\n", "", "");
    std::printf("0%.45sELAB =%7.2f MEV     %.65s\n\n",
                &reaction.reactStr[1], reaction.energies.eLab, &reaction.header[1]);

    setLog(2 * (reaction.kin.lOutMax + lxMax));

    //
    // WE WILL HAVE NO MORE IALLOC CALLS SO GET ADDRESSES
    //
    // walks them directly via msvalPointer/ja12sPointer 1-based pointers; jA12Of is
    // the literal 1, kA12mEnd is just jA12M, and the loop starts at 0
    double* msvalPointer = reaction.gridData.msvalArr.data();  // 0-based: [kA12m]=DMSVAL(kA12m+1)
    int*    ja12sPointer  = reaction.gridData.ja12sArr.data();  // 0-based: accessed [jA12Of+kA12m-1 ..]
    // pass &intsArr[intsOffset] as the 1-based "INTS+1" pointer.
    double* intsPointer = reaction.gridData.intsArr.data() + reaction.gridData.intsOffset;
    // for lh in 1..hCount and hsA12Pointer[lA12Of+lh] (= A12[jA12On+(lh-LXTMP)])
    // for the A12 fold.
    la12Vl  = reaction.inelastic.nMloLx;
    double* hsA12Pointer = reaction.gridData.hsA12Arr.data();  // 0-based
    double* a12T1Pointer = reaction.dwbaGrid.a12tm_.data();
    double* a12T2Pointer = reaction.dwbaGrid.a12tm_.data() + numLx;
    int*    a12T3p1Pointer = reinterpret_cast<int*>(reaction.dwbaGrid.a12tm_.data() + numLx + 2 * lxMax + 2);
    // DW class-owned by DWBAGrid::dw_; reaction.gridData.dwPointer = dw_.data() - 1 (1-based).
    double* dwPointer = reaction.gridData.dwPointer;
    int* lisPointer = reaction.inelastic.lisArr.data();  // 1-based (lisPointer[liIndex] valid for 1..nLValues)
    // rlv2Pointer/imv2Pointer/cent2Pointer are 0-based pointers into channel[2]'s vectors.
    double* rlv2Pointer = reaction.distortedWave.channel[2].rlvsArr.data();
    double* imv2Pointer = reaction.distortedWave.channel[2].imvsArr.data();
    double* cent2Pointer = reaction.distortedWave.channel[2].centrArr.data();
    // ATERM now class-owned by InelasticData; sFromI gets (data()-1) for 1-based access.
    // int* for WFI/WFO/DW TOC arrays
    int* wfIoPointer = reaction.gridData.iDwfoPointer;
    int* wfiIPointer = reaction.gridData.iDwfiPointer;

    // direct pool pointers — valid until arrays freed at function exit
    // ptr[i] == pool[lBase + i] for all i in range
    // LILOR/LILOI class-owned; use reaction.inelastic.liloRPointer/liloIPointer (0-based)
    double* liloRPointer  = reaction.inelastic.liloRPointer;
    double* liloIPointer  = reaction.inelastic.liloIPointer;
    double* abs1Pointer   = reaction.gridData.abs1Pointer;  // class-owned
    double* hsMlPointer   = hsA12Pointer;  // H values at hsA12Arr[0..NMLOLX-1], 0-based ptr[0]=H[1]
    // IHINT/IHABS class-owned; use reaction.gridData.hintPointer/habsPointer (1-based: ptr[hIndex]=H[hIndex])
    double* hintPointer   = reaction.gridData.hintPointer;
    double* habsPointer   = reaction.gridData.habsPointer;
    double* cosinPointer  = reaction.gridData.cosinPointer;  // class-owned 0-based cos table
    double* smhvlPointer  = reaction.gridData.smhvlPointer;
    double* rioExPointer  = reaction.gridData.rioExPointer;   // 0-based (accessed [riRoIndex - 1])
    // LIR/LII/LOR/LOI class-owned; 0-based float pointers
    float*  lirPointer    = reaction.gridData.lirPointer;
    float*  liiPointer    = reaction.gridData.liiPointer;
    float*  lorPointer    = reaction.gridData.lorPointer;
    float*  loiPointer    = reaction.gridData.loiPointer;
    float*  wioPointer    = reaction.gridData.wioPointer;  // class-owned
    // use GRIDCM pointer fields (set by gridSet)
    float*  phiPointer    = reaction.gridData.phiPointer;
    float*  phiPPointer   = reaction.gridData.phiPPointer;
    float*  phiTPointer   = reaction.gridData.phiTPointer;
    float*  trapWeightPointer    = reaction.gridData.trapWeightPointer;
    // pool pointer aliases
    int*    dwiPointer  = reaction.gridData.dwiPointer;
    int*    iiindxPointer = reaction.gridData.iiindxPointer;

    //
    // WE ASSUME THAT L1+L2+lx <= 256  AND  NUMCOSIN <= 16384
    //
    pi256 = 256.0 * Constants::PI;

    // DUMMY1() Fortran overlay-hint (force ThreeJ COMMON load) deleted

    factor = 2.0 * sqrt(reaction.kin.akIn * reaction.kin.akOut / (reaction.distortedWave.channel[1].Ecm * reaction.distortedWave.channel[2].Ecm));
    if (reaction.inelastic.densitySwitch) factor = sqrt(2.0) * factor;

    //
    // THE LOOPS BEGIN NOW
    //
    tt7 = (float)second();

    //
    // LOOP OVER li
    //
    // TWO PASSES OVER li LOOP NOW; EVEN li FIRST, THEN ODD
    //
    liMin = lMin + 1;
    liL   = lMin;
    if (lMin % 2 != 1) {
        liMin = lMin;
        liL   = lMin + 1;
    }

    // Fill one out-channel wavefunction slot. Identical body in the first-li
    // and subsequent-li branches below; captures lo/lasO/jPo/trkWriteOffset by
    // ref so it reproduces the inline loop body exactly.
    auto fillOutChannelWave = [&](int kSlot) {
        lo   = wfIoPointer[4*kSlot - 3] + li;
        lasO = wfIoPointer[4*kSlot - 2] + li;
        jPo  = wfIoPointer[4*kSlot - 1] + 2*li;
        wfIoPointer[4*kSlot] = trkWriteOffset;
        dispatchPartialWave(reaction, lo, jPo, 2, nRiRoInterp, reaction.gridData.roPointer,
              reaction.gridData.lorPointer + trkWriteOffset, reaction.gridData.loiPointer + trkWriteOffset,
              reaction.distortedWave.scatteringSolver.wavRPointer, reaction.distortedWave.scatteringSolver.wavIPointer,
              rlv2Pointer - 1, imv2Pointer - 1, cent2Pointer - 1);
        trkWriteOffset = trkWriteOffset + nRiRoInterp;
    };

    //
    // liParity = 1 FOR li EVEN, = 2 FOR li ODD
    //
    for (liParity = 1; liParity <= 2; liParity++) {

        if (liMin > lMax) { liMin = liL; continue; }  // skip to next liParity
        isAfterFirstLi = false;
        liIndex = 1;

        for (li = liMin; li <= lMax; li += 2) {

            //
            // GET THE REQUIRED SCATTERING WAVEFUNCTIONS.
            //
            if (!isAfterFirstLi) {
                //
                // FOR FIRST li, FILL UP THE OUT-CHANNEL WAVEFUNCTION ARRAY
                //
                trkWriteOffset = 0;
                for (kWO = 1; kWO <= nWfo; kWO++) {
                    fillOutChannelWave(kWO);
                }
                trkWriteOffset = 0;
                isAfterFirstLi = true;
            } else {
                //
                // ON SUBSEQUENT li'S WE NEED OUT-CHANNEL WAVE FUNCTIONS FOR ONE
                // NEW lo. DETERMINE WHERE THE NEW lo STARTS.
                //
                lo   = wfIoPointer[1] + 2;
                lasO = wfIoPointer[2] + 2;
                jPo  = wfIoPointer[3] + 4;
                i = 0;  // default: no match found
                for (ii = 1; ii <= nWfo; ii++) {
                    if (lo   == wfIoPointer[4*ii - 3]
                     && lasO == wfIoPointer[4*ii - 2]
                     && jPo  == wfIoPointer[4*ii - 1]) {
                        i = nWfo - ii + 1;
                        break;
                    }
                }
                // Move pointers down
                for (kWO = 1; kWO <= i; kWO++) {
                    srcIndex = kWO + (nWfo - i);
                    wfIoPointer[4*kWO] = wfIoPointer[4*srcIndex];
                }
                // Put new wave functions into circular buffer
                i = i + 1;
                for (kWO = i; kWO <= nWfo; kWO++) {
                    fillOutChannelWave(kWO);
                    if (trkWriteOffset > reaction.gridData.nCrit) trkWriteOffset = 0;
                }
            }

            //
            // FOR EACH li, GET THE IN-CHANNEL WAVE FUNCTIONS
            //
            for (kWI = 1; kWI <= reaction.gridData.nWfi; kWI++) {
                jPi  = wfiIPointer[3*kWI - 1] + 2*li;
                i    = wfiIPointer[3*kWI] + 1;
                dispatchPartialWave(reaction, li, jPi, 1, nRiRoInterp, reaction.gridData.riPointer,  // class-owned
                      reaction.gridData.lirPointer + i - 1, reaction.gridData.liiPointer + i - 1,
                      reaction.distortedWave.scatteringSolver.wavRPointer, reaction.distortedWave.scatteringSolver.wavIPointer,
                      reaction.distortedWave.channel[1].rlvsArr.data() - 1,
                      reaction.distortedWave.channel[1].imvsArr.data() - 1,
                      reaction.distortedWave.channel[1].centrArr.data() - 1);  // 1-based ptrs into class-owned vectors
            }

            //
            // WAVEFUNCTIONS ALL FOUND, WILL WE DO THIS VALUE OF li
            //
            // Search lisPointer for current li
            {
                bool found = false;
                while (true) {
                    int diff = li - lisPointer[liIndex];
                    if (diff < 0) break;
                    if (diff == 0) { found = true; break; }
                    if (liIndex >= nLValues) break;
                    liIndex++;
                }
                if (!found) continue;  // skip to next li
            }

            //
            // THIS IS AN li TO DO, FIND THE A12 COEFFICIENTS
            //
            // TT1 / times[4] dropped — times[4] was incremented here but never
            // read by any input_reader/ineld2 printout (times[1..3,5..8] only).

            angularCoupling12(li, reaction.inelastic.lxMin, lxMax, lMin, lMax,
                reaction.dwbaGrid.xlam.data(),
                reaction.dwbaGrid.nlam.data(),
                hsA12Pointer + la12Vl,                    // hsA12Arr base + A12 region start; A12 does a12Vl=arg-1
                reaction.gridData.msvalArr.data(),           // class-owned; A12 does msval=arg-1
                reaction.gridData.ja12sArr.data(),           // class-owned; A12 does JA12S=arg-1
                jA12M, jA12N, jA12An,
                reaction.gridData.iiindxPointer + 1,
                hCount, loMinMin, loMaxMax,
                intsPointer,
                a12T1Pointer,                           // reaction.dwbaGrid.a12tm_.data()
                a12T2Pointer,                           // reaction.dwbaGrid.a12tm_.data() + numLx
                a12T3p1Pointer,                         // int region into reaction.dwbaGrid.a12tm_
                la12Vl, printLevel,
                reaction.gridData.dwiPointer + 1,
                dwCount, numIi,
                reaction.gridData.iDwfiPointer + 1,
                &wfIoPointer[1],                          // iDwfo: 0-based macro
                reaction);
            hBlockSize = hCount * nPhiSum;

            if (jA12M > reaction.inelastic.a12mSize || jA12N > reaction.gridData.a12nSize) {
                std::printf(" **** ERROR *** A12 HAS TOO MANY TERMS: %8d%8d%8d%8d\n",
                            reaction.inelastic.a12mSize, jA12M, reaction.gridData.a12nSize, jA12N);
            }
            secondCallCount = secondCallCount + 2;

            if (loMinMin > loMaxMax) continue;  // skip this li

            //
            // NOW PROCEED TO THE (RI, rO) INTEGRATION LOOP
            //
            lhEnd = hCount;  // was LHSM1 + hCount with LHSM1 = 0
            phiPointIndex = 1;
            dwAmpMin = 10.0;
            dwAmpMax = 0.0;
            dwRioC = 0.0;
            dwAmpCount = 0.0;

            //
            // INITIALIZE "I1REAL", "I1IMAG" AND "RIROABS" FOR THIS li.
            //
            // use liloRPointer/liloIPointer/abs1Pointer
            for (ii = 1; ii <= numIi; ii++) {
                liloRPointer[ii - 1] = 0.0;
                liloIPointer[ii - 1] = 0.0;
                abs1Pointer[ii - 1] = 0.0;
            }

            tt5 = (float)second();

            //
            // LOOP OVER THE V-GRID AND FOR EACH ONE DO A COMPLETE U-INTEGRAL.
            //
            for (vIndex = 1; vIndex <= nPhiDifference; vIndex++) {

                //
                // GET THE H'S FOR INTERPOLATION INPUT
                //
                // INDEX TO THE H'S FILE
                //
                hBlockTag = 1000 * (1000 * liParity + li) + vIndex;

                //
                // LOOP OVER U-VALUES AND COMPUTE THE H'S
                //
                riRoIndex = nPhiSum * (vIndex - 1);
                tt2 = (float)second();

                for (uIndex = 1; uIndex <= nPhiSum; uIndex++) {

                    riRoIndex = riRoIndex + 1;

                    for (hIndex = 1; hIndex <= hCount; hIndex++) {
                        hsMlPointer[hIndex - 1] = 0.0;
                        hintPointer[hIndex - 1] = 0.0;
                        habsPointer[hIndex - 1] = 0.0;
                    }

                    //
                    // THE DX INTEGRAL (AT LAST)
                    //
                    for (ii = 1; ii <= nPhiPoints; ii++) {

                        // use phiPointer/phiPPointer/phiTPointer/trapWeightPointer
                        phi  = (double)phiPointer[phiPointIndex];
                        kA12mEnd = jA12M;  // LMSVAL=0 was inlined
                        phiP = (double)phiPPointer[phiPointIndex];
                        arg  = phi + phi;
                        phiT = (double)phiTPointer[phiPointIndex];
                        trapWeight = (double)trapWeightPointer[phiPointIndex];
                        jA12Of = 1;  // LJA12S=LMSVAL=0 inlined

                        //
                        // HIGH SPEED INLINE COSINE AND SINE(arg)
                        // arg = 2*phi, finding COS(2PHI) and SIN(2PHI)
                        //
                        cosIndexReal = (double)(int)(arg * reaction.gridData.stepInverse + 0.5);
                        cosIndex = (int)cosIndexReal;
                        cosIndex = cosIndex % nCosin;
                        deltaC = arg - cosIndexReal * reaction.gridData.cosStep;
                        // use cosinPointer
                        sin2Phi = ((1.0 + halfNeg * deltaC * deltaC)
                                * cosinPointer[std::abs(cosIndex - cosinQuarter)]
                                + deltaC * cosinPointer[cosIndex]);
                        cos2Phi = (1.0 + halfNeg * deltaC * deltaC)
                                * cosinPointer[cosIndex]
                                - deltaC * cosinPointer[std::abs(cosIndex - cosinQuarter)];

                        // THE COSINE AND SINE(arg) HAS BEEN FOUND
                        sin2PhiSq = sin2Phi * sin2Phi;

                        //
                        // LOOP OVER TERMS IN A12
                        //
                        for (kA12m = 0; kA12m <= kA12mEnd; kA12m += 5) {  // LMSVAL=0 inlined

                            //
                            // THE MAGIC angle FOR THE FIRST VALUE OF MU-2
                            // WE ADD 256*PI TO FORCE IT POSITIVE
                            //
                            arg = pi256 + msvalPointer[kA12m] * phiT  // msvalArr 0-based
                                        + msvalPointer[kA12m + 1] * phiP
                                        - msvalPointer[kA12m + 2] * phi;

                            lA12Of = ja12sPointer[jA12Of + kA12m - 1];  // ja12sArr 0-based

                            //
                            // HIGH SPEED INLINE COSINE AND SINE(arg)
                            //
                            muCount = ja12sPointer[jA12Of + kA12m];  // ja12sArr 0-based
                            cosIndexReal = (double)(int)(arg * reaction.gridData.stepInverse + 0.5);
                            cosIndex = (int)cosIndexReal;
                            cosIndex = cosIndex % nCosin;
                            deltaC = arg - cosIndexReal * reaction.gridData.cosStep;

                            //
                            // cosTheta = COS(MT*phiT+MP*phiP-MU*phi)
                            // sinSin = SIN(MT*phiT+MP*phiP-MU*phi) * SIN(2*phi)
                            //
                            // use cosinPointer
                            cosTheta = (1.0 + halfNeg * deltaC * deltaC)
                                    * cosinPointer[cosIndex]
                                    - deltaC * cosinPointer[std::abs(cosIndex - cosinQuarter)];
                            sinSin = ((1.0 + halfNeg * deltaC * deltaC)
                                    * cosinPointer[std::abs(cosIndex - cosinQuarter)]
                                    + deltaC * cosinPointer[cosIndex]) * sin2Phi;

                            // THE COSINE AND SINE(arg) HAS BEEN FOUND
                            lhStart = ja12sPointer[jA12Of + kA12m + 1];  // ja12sArr 0-based
                            cosTheta = trapWeight * cosTheta;
                            hStride = ja12sPointer[jA12Of + kA12m + 2];    // ja12sArr 0-based
                            sinSin = trapWeight * sinSin;

                            //
                            // INNERMOST LOOP OVER MU'S.
                            // EACH MU IS 2 GREATER THAN THE PREVIOUS, SO THETA IS DECREMENTED
                            // BY 2PHI AND COS(THETA) IS FOUND ITERATIVELY.
                            //
                            for (mu = 1; mu <= muCount; mu++) {
                                cosNew = cosTheta * cos2Phi + sinSin;
                                sinSin = sinSin * cos2Phi - cosTheta * sin2PhiSq;
                                cosTheta = cosNew;

                                // LOOP OVER ALL (lx, lo) PAIRS FOR THIS li
                                for (lh = lhStart; lh <= lhEnd; lh++) {
                                    // pool[lh] = pool[lh] + cosNew*pool4[lA12Of+lh]
                                    hsA12Pointer[lh - 1] += cosNew * hsA12Pointer[lA12Of + lh - 1];  // H at vec[0..NMLOLX-1], A12 at vec[NMLOLX..]
                                }
                                lA12Of = lA12Of + hStride;
                            }

                        }  // end of kA12m loop (A12 terms)

                        phiPointIndex = phiPointIndex + 1;

                        //
                        // ACCUMULATE THE DESIRED INTEGRAL IN LHINT AND THE INTEGRAL
                        // OF |INTEGRAND| IN LHABS.
                        //
                        // use hintPointer/habsPointer/hsMlPointer
                        for (hIndex = 1; hIndex <= hCount; hIndex++) {
                            hintPointer[hIndex - 1] += hsMlPointer[hIndex - 1];
                            habsPointer[hIndex - 1] += fabs(hsMlPointer[hIndex - 1]);
                            hsMlPointer[hIndex - 1] = 0.0;
                        }

                    }  // end of ii=1,nPhiPoints (DX integral terms)

                    //
                    // END OF DX INTEGRAL TERMS
                    //
                    innerLoopPassCount = innerLoopPassCount + jA12N;
                    cosineIterCount = cosineIterCount + jA12An;
                    m1m2LoopPassCount = m1m2LoopPassCount + jA12M;
                    phiLoopPassCount = phiLoopPassCount + 1;

                    //
                    // STORE THE H'S WITH THE factor EXP(-ALPHAP*rP - ALPHAT*rT)
                    // REMOVED IN THE INTERPOLATION INPUT ARRAY
                    //
                    // use smhvlPointer/hintPointer/rioExPointer
                    { double rioEx = rioExPointer[riRoIndex - 1];
                      for (hIndex = 1; hIndex <= hCount; hIndex++)
                          smhvlPointer[uIndex + nPhiSum * (hIndex - 1)] = hintPointer[hIndex - 1] * rioEx; }
                    if (verbosity >= 9) {
                        std::printf(" CMPTED: %9d%5d", hBlockTag, riRoIndex);
                        for (hIndex = 1; hIndex <= hCount; hIndex++) {
                            if ((hIndex - 1) % 7 == 0 && hIndex > 1) std::printf("\n%22s", "");
                            std::printf(" %14.5G", hintPointer[hIndex - 1]);
                        }
                        std::printf("\n");
                    }

                }  // end of uIndex=1,nPhiSum (U loop)

                //
                // END OF LOOP ON U
                //
                reaction.timing.times[5] = reaction.timing.times[5] + (float)second() - tt2;
                secondCallCount = secondCallCount + 2;


                riRoProcessedCount = riRoProcessedCount + nPhiSum;
                hComputedCount = hComputedCount + hBlockSize;

                //
                // NOW INTERPOLATE
                //
                tt2 = (float)second();
                for (hIndex = 1; hIndex <= hCount; hIndex++) {
                    // use GRIDCM pointer fields instead of lookups
                    naturalCubicSpline(nPhiSum, reaction.gridData.smhptsPointer, reaction.gridData.smhvlPointer + (hIndex-1)*nPhiSum,
                           reaction.gridData.smhwkPointer, reaction.gridData.smhwkPointer + nPhiSum,
                           reaction.gridData.smhwkPointer + 2 * nPhiSum);
                    cubicSplineInterp(nPhiSum, reaction.gridData.smhptsPointer, reaction.gridData.smhvlPointer + (hIndex-1)*nPhiSum,
                           reaction.gridData.smhwkPointer, reaction.gridData.smhwkPointer + nPhiSum,
                           reaction.gridData.smhwkPointer + 2 * nPhiSum,
                           nInterpPoints, reaction.gridData.smiptsPointer,
                           reaction.gridData.smivlPointer + (hIndex-1)*nInterpPoints);
                    hInterpCount = hInterpCount + nInterpPoints;
                }
                reaction.timing.times[6] = reaction.timing.times[6] + (float)second() - tt2;
                secondCallCount = secondCallCount + 2;

                //
                // NOW DO INTEGRAL DU FOR FIXED V.
                //
                riRoIndex = nInterpPoints * (vIndex - 1);
                for (uIndex = 1; uIndex <= nInterpPoints; uIndex++) {
                    riRoIndex = riRoIndex + 1;

                    //
                    // COMPUTE IMAG AND REAL PARTS OF PRODUCT OF THE DISTORTED
                    // SCATTERING WAVES.
                    //
                    dwAmpPeak = 0.0;
                    for (kDw = 1; kDw <= dwCount; kDw++) {
                        i = dwiPointer[4*kDw - 1] + riRoIndex;
                        // use lirPointer/liiPointer/lorPointer/loiPointer
                        fiR = (double)lirPointer[i];
                        fiI = (double)liiPointer[i];
                        i = dwiPointer[4*kDw] + riRoIndex;
                        foR = (double)lorPointer[i];
                        foI  = (double)loiPointer[i];
                        dwReal = foR * fiR - foI * fiI;
                        dwImag = foR * fiI + foI * fiR;
                        // Guard against NaN from WAVELJ overflow (outgoing L=0)
                        if (std::isnan(dwReal) || std::isinf(dwReal)) dwReal = 0.0;
                        if (std::isnan(dwImag) || std::isinf(dwImag)) dwImag = 0.0;
                        dwAmpThis = std::max(fabs(dwReal), fabs(dwImag));
                        if (dwAmpThis > dwAmpMax) dwAmpMax = dwAmpThis;
                        if (dwAmpThis < dwAmpMin && dwAmpThis != 0.0) dwAmpMin = dwAmpThis;
                        dwRioC = dwRioC + dwAmpThis;
                        dwAmpPeak = std::max(dwAmpPeak, dwAmpThis);
                        dwPointer[2*kDw - 2] = dwReal;
                        dwPointer[2*kDw - 1] = dwImag;
                        dwAmpCount = dwAmpCount + 1.0;
                    }

                    //
                    // SKIP THE hasNextBlock PART IF ALL DISTORTED-WAVE PRODUCTS ARE SMALL.
                    //
                    if (dwAmpPeak < 1.0e-30) {
                        // skip — DW product too small
                    } else {

                    riRoWeight = (double)wioPointer[riRoIndex];  // use wioPointer

                    //
                    // NOW ADD THESE H'S INTO THE APPROPRIATE i(JO,JI,lo,lx)
                    // ALONG WITH THE PRODUCT OF THE DISTORTED WAVES
                    //
                    for (ii = 1; ii <= numIi; ii++) {
                        smivlIndex = iiindxPointer[4*ii - 3] + uIndex;
                        // iIndex(2,ii) now holds the 1-based dw_ subscript for the
                        // imaginary half of pair (writer: angular_coupling.cpp:249, kDw = 2*KDW_orig);
                        // dwPointer[kDw-2] = dwReal, dwPointer[kDw-1] = dwImag (dwPointer = dw_.data(), 0-based).
                        kDw   = iiindxPointer[4*ii - 2];
                        double dwReal = dwPointer[kDw - 2], dwImag = dwPointer[kDw - 1];
                        if (std::isnan(dwReal) || std::isinf(dwReal)) dwReal = 0.0;
                        if (std::isnan(dwImag) || std::isinf(dwImag)) dwImag = 0.0;
                        {
                        // smivlArr is class-owned; smivlIndex is now a 0-based offset (writer at angular_coupling:248).
                        double smivlValue = reaction.gridData.smivlArr.data()[smivlIndex];
                        double addR = riRoWeight * smivlValue * dwReal;
                        double addI = riRoWeight * smivlValue * dwImag;
                        liloRPointer[ii - 1] += addR;
                        liloIPointer[ii - 1] += addI;
                        }

                        //
                        // COMPUTE AN INDICATION OF LOSS OF SIGNIFICANCE:
                        //   "RIROABS" = INTEGRAL(RI,rO) |INTEGRAL(PHI) INTEGRAND|
                        //
                        if (isInfoPrint) {
                            abs1Pointer[ii - 1] += fabs(riRoWeight * reaction.gridData.smivlArr.data()[smivlIndex])
                                * (fabs(dwPointer[kDw - 2]) + fabs(dwPointer[kDw - 1]));
                        }
                    }

                    }
                }  // end of uIndex=1,nInterpPoints

            }  // end of vIndex=1,nPhiDifference

            reaction.timing.times[7] = reaction.timing.times[7] + (float)second() - tt5;
            secondCallCount = secondCallCount + 2;

            //
            // END OF (RI,rO) INTEGRATION
            //
            dwRioC = dwRioC / dwAmpCount;
            if (debugSwitch) {
                std::printf(" SMALLEST, LARGEST AND AVERAGE std::fabs(DWRIRO)= %15.5G%15.5G%15.5G\n\n",
                            dwAmpMin, dwAmpMax, dwRioC);
            }

            //
            //
            sFromI(li, liIndex,
                   reaction.inelastic.smatRPointer,  // 0-based class-owned
                   reaction.inelastic.smatIPointer,  // 0-based class-owned
                   reaction.inelastic.indxsArr.data(),  // 0-based INDXS
                   reaction.inelastic.liloRPointer,          // class-owned
                   reaction.inelastic.liloIPointer,          // class-owned
                   reaction.gridData.iiindxPointer + 1,
                   numIi,
                   reaction.gridData.dwiPointer + 1,
                   reaction.gridData.iDwfiPointer + 1,
                   &wfIoPointer[1],
                   reaction.gridData.abs1Pointer,           // class-owned
                   reaction.inelastic.atermArr.data(),  // 0-based ATERM[LXP] via class-owned vector
                   factor, isInfoPrint, reaction);


            if (phiPointIndex - 1 > reaction.gridData.maxCount) {
                std::printf(" **** ERROR ****  PHI COUNTS MESSED UP: %10d%10d\n",
                            phiPointIndex, reaction.gridData.maxCount);
            }
            i = lhStart + lA12Of - la12Vl - 1;
            if (i != jA12N) {
                std::printf(" ******* A12 COUNTS MESSED UP: %10d%10d%10d%10d%10d%10d%10d%10d\n",
                            i, jA12N, li, jA12M, hCount, lhStart, lA12Of, la12Vl);
            }

        }  // end of li=liMin,lMax,2

        //
        // END OF THE li LOOP FOR A GIVEN parity OF li
        //
        liMin = liL;

    }  // end of liParity=1,2

    tt = (float)second();
    reaction.timing.times[8] = tt - tt7;

    //
    // END OF THE 2 li LOOPS (EVEN AND ODD parity OF li)
    //
    // SAVEHS/USEHS post-li cleanup (savehs_mode==1 dummy-record write and
    //
    // A L L    D O N E
    //
    ineld2(tStart, hInterpCount, secondCallCount, innerLoopPassCount, cosineIterCount,
           m1m2LoopPassCount, phiLoopPassCount, riRoProcessedCount, hComputedCount, reaction);
}
