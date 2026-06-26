// input_reader.cpp — INRDIN: inelastic excitation radial integrals i(li,lo,lx).
// Integrates U(R)*H(R)*U(R) (H = nuclear + Coulomb coupling) over the distorted
// waves, adds the asymptotic correction, and subtracts the pure-Coulomb piece.

#include "MathTables.h"
#include "Timing.h"
#include "math/angular_momentum_coeff.h"
#include <cstdio>
#include <cmath>
#include <complex>
#include "Reaction.h"
#include "Constants.h"

// INDXJS(4,2) <=> jpMin (INELCM)
// Fortran: INDXJS(i,J) = jpMin + 4*(J-1) + (i-1)
// C++: indxJs(reaction, i, J) = (&reaction.inelastic.jpMin)[ 4*(J-1) + (i-1) ]
inline int indxJs(Reaction& reaction, int i, int J) {
    return (&reaction.inelastic.jpMin)[ 4*(J-1) + (i-1) ];
}

// Print the INRDIN page header (P T O L E M Y / COMPUTATION OF INELASTIC
// S-MATRIX ELEMENTS / reaction string + ELAB + header). Emitted at routine
// start and again on each page break; extracted verbatim. reactStr/header are
// 1-based char pointers (entries [1]..[45] / [1]..[65]).
static void printInelasticHeader(const char* reactStr, const char* header, double eLab)
{
    char react45[46]; react45[45] = '\0';
    char header65[66]; header65[65] = '\0';
    for (int ii = 0; ii < 45; ii++) react45[ii]  = reactStr[ii+1];
    for (int ii = 0; ii < 65; ii++) header65[ii] = header[ii+1];
    std::printf("1%58sP T O L E M Y\n", "");
    std::printf("%49sCOMPUTATION OF INELASTIC S-MATRIX ELEMENTS\n", "");
    std::printf("0%.45sELAB =%7.2f MEV     %.65s\n\n",
                react45, eLab, header65);
}

void DWBAGrid::inelasticRadialIntegrals(Reaction& reaction)
{
    //
    // COMPUTES THE INELASTIC EXCITATION RADIAL INTEGRALLS  i(li,lo,lx)
    //
    // COMPUTES THE INELASTIC EXCITATION INTEGRALS
    //   INTEGRAL(0 TO sumMax) DR U(R) H U(R)
    // WHERE H IS THE FULL H
    // AND THEN COMPUTES
    //   i(NUCLEAR) = i(0 TO sumMax) + i(sumMax TO INF) - i(COULOMB)
    //

    // WARNING:  VAX VERSION MAY BE INACCURATE (NO COMPLEX*16).
    // imagUnit = (0, 1) = i
    std::complex<double> imagUnit(0.0, 1.0);

    // FOLLOWING ARE  -(i)**lx
    // PHASES(4):  -(i)^1 = -i = (0,-1)
    //             -(i)^2 = -(-1) = 1... wait — Fortran says:
    // These are indexed 1..4 (1-based), so phases[0]=(0,-1), [1]=(-1,0), [2]=(0,1), [3]=(1,0)
    static const std::complex<double> phases[4] = {
        {0.0, -1.0},   // PHASES(1)
        {-1.0, 0.0},   // PHASES(2)
        {0.0,  1.0},   // PHASES(3)
        {+1.0, 0.0}    // PHASES(4)
    };


    // /FLOAT/
    double& eLab   = reaction.energies.eLab;

    int& lMin    = reaction.angMom.lMin;
    int& lMax    = reaction.angMom.lMax;
    int& printLevel  = reaction.flags.printLevel;

    // /INTRNL/
    int& iExcit  = reaction.internalState.iExcit;

    // /INELCM/
    int& lxMin   = reaction.inelastic.lxMin;
    int& lxMax   = reaction.inelastic.lxMax;
    int& nMloLx  = reaction.inelastic.nMloLx;
    int& nLValues    = reaction.inelastic.nLValues;
    int& nLx     = reaction.inelastic.nLx;
    int& nSpl    = reaction.inelastic.nSpl;
    // poolBetas now holds std::vector<double>* into reaction.named (set by probe_print).
    // 0-based: [0]=BETA, [1]=BETACOUL (sentinel nullptr after probe_print), [2]=BETARATS.
    std::vector<double>** poolBetas = reaction.inelastic.poolBetas;
    int& densitySwitch  = reaction.inelastic.densitySwitch;

    double& akIn    = reaction.kin.akIn;
    double& akOut    = reaction.kin.akOut;
    int&    lOutMax = reaction.kin.lOutMax;

    float*  times  = reaction.timing.times;   // 1-based: times[1]..times[8]

    // /HEDCOM/
    char*  reactStr  = reaction.reactStr;   // 1-based: reactStr[1]..reactStr[45]
    char*  header = reaction.header;  // 1-based: header[1]..header[65]

    // INRDIN-specific aliases for inelastic GRIDCM fields. The Fortran
    // 1-based Fortran offset arithmetic that the overlay used to encode.
    double& rcuEff = reaction.gridData.jacobian;        // rcuEff
    // permanent 0 sentinels in this path; only the pool-free guards used them.
    int&   numPoint   = reaction.gridData.nRiRoInterp;          // NUMPT
    int&   nCrit = reaction.gridData.nCrit;

    // --- Local variables ---
    int debugSwitch, isInfoPrint;   // LOGICAL (0/1)
    int isAfterFirstLi;

    float tStart;

    int verbosity;
    int jProj;

    double factor;

    int lRange;
    int liMin, liL, liParity;
    int liIndex, li;
    int loMin, loMax;
    // loBase init to 0 silences -Wmaybe-uninitialized for the empty-li-loop
    // edge — loBase is set on the first li iteration (line ~325 in the
    // !isAfterFirstLi branch) and updated conditionally on subsequent iterations.
    // Line ~420's read assumes the for-li loop ran at least once.
    int trkWriteOffset;
    int loBase = 0;
    int lo, lx;
    int hIndex, uIndex;
    int lineCount;
    int i, i1, i2;
    int k, kBase, lDeltaMin, kOffset;
    int trkReadOffset;
    int loMinMin, loMaxMax;
    int loMn, loMx;
    int dwOffset;

    double rValue, rat;
    double fiR, fiI, foR, foI;
    double dwReal, dwImag, dwAmpThis;
    double dwAmpMin, dwAmpMax, dwRioC, dwAmpCount;
    double hNucReal, hNucImag, hCoulomb;
    double hReal, hImag;
    double termReal, termImag;

    double betaRatio;
    double coefficient, ffI;
    double ampCoulomb, ampNuclear, ampTotal;
    double phaseCoulomb, phaseNuclear, phaseTotal, absRatio;

    std::complex<double> nuclearAmp, compAmp, rToInAmp, totalAmp;
    std::complex<double> sIn, sOut;   // S-matrix scalar elements; camelCase capital I disambigs sIn from C math sin
    std::complex<double> phase;

    double totalTime, waveTime, otherTime, recurTime1, recurTime2;

    // =========================================================================
    // Begin executable code
    // =========================================================================

    tStart = (float)second();
    // discarded iret after the call, and inelasticRadialIntegrals only ever wrote returnCode=1
    // at the end.

    for (int iii = 1; iii <= 3; iii++) times[iii] = 0.0f;

    verbosity = printLevel % 10;
    debugSwitch = (verbosity >= 4);
    isInfoPrint = (verbosity >= 2);

    //
    //   49X, 'COMPUTATION OF INELASTIC S-MATRIX ELEMENTS' /
    //   '0', 45A1, 'eLab =', F7.2, ' MEV', 5X, 65A1 / )
    //
    printInelasticHeader(reactStr, header, eLab);

    lineCount = 1000;

    //
    // NOTE -- NO SPIN ORBIT FORCE
    //
    jProj = 1;

    //
    // CALL setLog ( 2*(lOutMax+lxMax) )
    //
    setLog(2 * (lOutMax + lxMax));

    //
    // HERE WE START THE REAL GUTS OF THE CALCULATION.
    // WE WILL HAVE NO MORE IALLOC CALLS SO GET ADDRESSES
    //
    // (smiptsPointer = 1-based RPTS, rpts4Pointer = 1-based RPTS4).
    // LIR/LII/LOR/LOI replaced by reaction.gridData.lirPointer etc (class-owned float* vectors)
    // LDW no longer used — dwPointer below derives directly from reaction.gridData.dwPointer.
    // cl2ffArr is 0-based class-owned vector.
    // bratPointer base = BETARATS vector start, offset by -lxMin/2 so bratPointer[lx/2]
    double* bratBase = poolBetas[2]->data() - lxMin/2;
    // rlvsArr/imvsArr/centrArr are per-channel class-owned vectors (0-based).
    double* rlv1Pointer = reaction.distortedWave.channel[1].rlvsArr.data();
    double* rlv2Pointer = reaction.distortedWave.channel[2].rlvsArr.data();
    double* imv1Pointer = reaction.distortedWave.channel[1].imvsArr.data();
    double* imv2Pointer = reaction.distortedWave.channel[2].imvsArr.data();
    double* cent1Pointer = reaction.distortedWave.channel[1].centrArr.data();
    double* cent2Pointer = reaction.distortedWave.channel[2].centrArr.data();

    int*    lisPointer = reaction.inelastic.lisArr.data();  // 1-based, lisPointer[K] valid for 1..nLValues
    int*    indxsPointer = reaction.inelastic.indxsPointer;  // 0-based class-owned
    // direct pointers — all valid until the backing arrays are freed at function exit
    double* rptsPointer   = reaction.gridData.smiptsPointer;    // 1-based
    double* nuclearHPointer   = reaction.gridData.nuclearHArr.data();  // 0-based, nuclearHPointer[2*uIndex-2] = nuclearHArr[2*uIndex-2]
    double* coulombHPointer  = reaction.gridData.coulombHArr.data();  // 0-based, coulombHPointer[uIndex-1] = coulombHArr[uIndex-1]
    double* smatRPointer  = reaction.inelastic.smatRPointer;   // 0-based class-owned
    double* smatIPointer  = reaction.inelastic.smatIPointer;   // 0-based class-owned
    double* abs1Pointer   = reaction.gridData.abs1Pointer;  // class-owned
    double* bnratPointer  = reaction.inelastic.betnrArr.data() - lxMin/2;  // class-owned BETANRAT
    double* bratPointer    = bratBase;            // bratPointer[lx/2] == BETARATS[(lx-lxMin)/2]
    float*  lirPointer    = reaction.gridData.lirPointer;  // class-owned
    float*  liiPointer    = reaction.gridData.liiPointer;  // class-owned
    float*  lorPointer    = reaction.gridData.lorPointer;  // class-owned
    float*  loiPointer    = reaction.gridData.loiPointer;  // class-owned
    // dwPointer now 0-based view of DWBAGrid::dw_.
    // reaction.gridData.dwPointer is now 0-based (dw_.data()), used directly as the base.
    double* dwPointer      = reaction.gridData.dwPointer;   // dwPointer[lo-LOMNMN] == dw_[lo-LOMNMN]

    // Decode the channel index K for the current (lx, iExcit) and read its
    //   K = lx+1 + nLx*INDXJS(4,3-IEXCIT)*(lx-INDXJS(3,IEXCIT)/2)
    auto computeKBase = [&]() {
        k = lx + 1 + nLx * indxJs(reaction, 4, 3-iExcit) * (lx - indxJs(reaction, 3, iExcit)/2);
        kBase     = indxsPointer[3*k - 3];
        lDeltaMin = indxsPointer[3*k - 2];
    };

    // factor = std::sqrt( akIn*akOut/(Ecm_ch(1)*Ecm_ch(2)*PI) )
    factor = std::sqrt(akIn * akOut / (reaction.distortedWave.channel[1].Ecm * reaction.distortedWave.channel[2].Ecm * Constants::PI));
    if (densitySwitch) factor = std::sqrt(2.0) * factor;

    //
    lRange = reaction.gridData.noFlo - 1;

    //
    // THE LOOPS BEGIN NOW
    //
    //
    // LOOP OVER li
    //
    // TWO PASSES OVER li LOOP NOW; EVEN li FIRST, THEN ODD
    //
    // liMin = lMin + 1
    liMin = lMin + 1;
    // liL   = lMin
    liL   = lMin;
    if (((lMin) % (2)) != 1) {
        // liMin = lMin (only when lMin is even)
        liMin = lMin;
        // liL   = lMin+1
        liL   = lMin + 1;
    }

    // OUT-CHANNEL wavefunction solve — the identical 13-arg solvePartialWave
    // call used both to fill the array for the first li and to add each
    // subsequent li's new out-channel wavefunc. trkWriteOffset is captured by
    // ref so each call reads its current value, exactly as the inline code did.
    auto solveOutChannel = [&](int loArg) {
        reaction.distortedWave.scatteringSolver.solvePartialWave(loArg, jProj, 2, numPoint, reaction.gridData.rpts4Pointer,
               lorPointer + trkWriteOffset, loiPointer + trkWriteOffset,
               reaction.distortedWave.scatteringSolver.wavRPointer, reaction.distortedWave.scatteringSolver.wavIPointer,
               (rlv2Pointer - 1), (imv2Pointer - 1), (cent2Pointer - 1), reaction);
    };

    // liParity = 1 FOR li EVEN,  = 2 FOR li ODD
    for (liParity = 1; liParity <= 2; liParity++) {

        if (liMin > lMax) { liMin = liL; continue; }  // skip to next liParity
        isAfterFirstLi = 0;
        // liIndex = 1
        liIndex = 1;

        for (li = liMin; li <= lMax; li += 2) {

            //
            // GET THE REQUIRED SCATTERING WAVEFUNCTIONS.
            // WE ALWAYS DO THIS EVEN IF li WON'T BE USED SO AS TO GET THE
            // ELASTIC S-MATRIX ELEMENTS.
            //
            loMin = li - lRange;
            loMax = li + lRange;

            if (!isAfterFirstLi) {
                //
                // FOR FIRST li, FILL UP THE OUT-CHANNEL WAVEFUNCTION ARRAY
                //
                trkWriteOffset = 0;

                for (lo = loMin; lo <= loMax; lo += 2) {
                    solveOutChannel(lo);
                    trkWriteOffset = trkWriteOffset + numPoint;
                }

                trkWriteOffset   = 0;
                loBase = loMin;
                isAfterFirstLi = 1;
            } else {
                //
                // ON SUBSEQUENT li'S WE NEED ONE NEW OUT-CHANNEL WAVEFUNC.
                //
                lo = li + lRange;
                solveOutChannel(lo);
                if (trkWriteOffset == 0) loBase = lo;
                trkWriteOffset = trkWriteOffset + numPoint;
                if (trkWriteOffset > nCrit) trkWriteOffset = 0;
            }

            //
            // FOR EACH li, WE MUST GET THE INCOMING CHANNEL WAVE FUNCTION
            //
            reaction.distortedWave.scatteringSolver.solvePartialWave(li, jProj, 1, numPoint, reaction.gridData.rpts4Pointer,
                   lirPointer, liiPointer,
                   reaction.distortedWave.scatteringSolver.wavRPointer, reaction.distortedWave.scatteringSolver.wavIPointer,
                   (rlv1Pointer - 1), (imv1Pointer - 1), (cent1Pointer - 1), reaction);

            //
            // WAVEFUNCTIONS ALL FOUND, WILL WE DO THIS VALUE OF li
            //
            // Search lisPointer for current li
            {
                bool found = false;
                while (true) {
                    int diff = li - lisPointer[liIndex];
                    if (diff < 0) break;       // li not in list — skip
                    if (diff == 0) { found = true; break; }  // found
                    if (liIndex >= nLValues) break;  // past end — skip
                    liIndex++;
                }
                if (!found) continue;  // skip to next li in DO 959 loop
            }

            //
            // THIS IS AN li TO DO
            //
            loMinMin = std::max(((loMax) % (2)), loMin);
            loMaxMax = loMax;
            dwAmpMin = 10.0;
            dwAmpMax = 0.0;
            dwRioC = 0.0;
            dwAmpCount = 0.0;

            //
            // IN ABS1 WE ACCUMULATE INTEGRAL |INTEGRAND| FOR THIS li.
            //
            for (hIndex = 1; hIndex <= nMloLx; hIndex++) {
                abs1Pointer[hIndex - 1] = 0.0;
            }

            //
            // LOOP OVER R-VALUES AND COMPUTE THE H'S
            //
            for (uIndex = 1; uIndex <= numPoint; uIndex++) {
                // uses rptsPointer
                rValue = rptsPointer[uIndex];
                // rat = 1/rValue
                rat = 1.0 / rValue;
                if (rValue < rcuEff) rat = rValue / (rcuEff * rcuEff);
                // rat = rat**2
                rat = rat * rat;

                //
                // COMPUTE IMAG AND REAL PARTS OF PRODUCT OF THE DISTORTED
                // SCATTERING WAVES.
                //
                // uses lirPointer
                fiR = (double)lirPointer[uIndex];
                // uses liiPointer
                fiI = (double)liiPointer[uIndex];
                // Debug: print first li, first uIndex
                if (std::fabs(fiR) + std::fabs(fiI) < 1.0e-20) {
                    // Skip H calculation — distorted wave product too small
                } else {

                //
                // UNDO THE WRAP AROUND STORAGE AND COMPUTE ALL PSI(A)*PSI(B)'S
                // FOR THIS (RI, rO).
                //
                for (lo = loMinMin; lo <= loMaxMax; lo += 2) {
                    // trkReadOffset = NUMPT*( (lo-loBase)/2 )
                    trkReadOffset = numPoint * ((lo - loBase) / 2);
                    if (trkReadOffset < 0) trkReadOffset = trkReadOffset + numPoint + nCrit;
                    foR = (double)lorPointer[trkReadOffset + uIndex];
                    // uses loiPointer
                    foI  = (double)loiPointer[trkReadOffset + uIndex];
                    if (std::fabs(foR) < 1.0e-20) foR = 0.0;
                    if (std::fabs(foI)  < 1.0e-20) foI  = 0.0;
                    // dwReal=FOR*fiR-foI*fiI
                    dwReal = foR * fiR - foI * fiI;
                    // dwImag=FOR*fiI+foI*fiR
                    dwImag  = foR * fiI + foI * fiR;
                    // dwAmpThis=std::max(std::fabs(dwReal),std::fabs(dwImag))
                    dwAmpThis = std::max(std::fabs(dwReal), std::fabs(dwImag));
                    if (dwAmpThis > dwAmpMax) dwAmpMax = dwAmpThis;
                    if ((dwAmpThis < dwAmpMin) && (dwAmpThis != 0.0)) dwAmpMin = dwAmpThis;
                    // dwRioC=dwRioC+dwAmpThis
                    dwRioC = dwRioC + dwAmpThis;
                    { int off = lo - loMinMin; dwPointer[off] = dwReal; dwPointer[off+1] = dwImag; }
                    dwAmpCount = dwAmpCount + 1.0;
                }

                // uses nuclearHPointer
                hNucReal = nuclearHPointer[2*uIndex - 2];
                // uses nuclearHPointer
                hNucImag = nuclearHPointer[2*uIndex - 1];
                // uses coulombHPointer
                hCoulomb = coulombHPointer[uIndex - 1];

                //
                // NOW ADD THESE H'S INTO THE APPROPRIATE i(li,lo,lx)
                // ALONG WITH THE PRODUCT OF THE DISTORTED WAVES
                // NOTE THAT THE BETA ARRAYS ALSO CONTAIN
                //      CLEBSCH( jIn lx K 0 : JOUT K )
                //
                hIndex = 0;
                for (lx = lxMin; lx <= lxMax; lx += 2) {
                    { int lxHalf = lx/2; hReal = hNucReal * bnratPointer[lxHalf] + hCoulomb * bratPointer[lxHalf];
                    hImag = hNucImag * bnratPointer[lxHalf]; }
                    // hCoulomb = rat * hCoulomb
                    hCoulomb = rat * hCoulomb;
                    // loMn = std::abs(li-lx)
                    loMn = std::abs(li - lx);
                    // loMx = li+lx
                    loMx = li + lx;
                    // K / KBASE / LDELMN decode (see computeKBase lambda above)
                    computeKBase();
                    for (lo = loMn; lo <= loMx; lo += 2) {
                        // kOffset = KBASE + (lo-li-LDELMN)/2
                        kOffset = kBase + (lo - li - lDeltaMin) / 2;
                        // i = (liIndex-1)*nSpl + kOffset
                        i = (liIndex - 1) * nSpl + kOffset;
                        dwOffset = lo - loMinMin;  // now offset into dwPointer
                        // hIndex = hIndex + 1
                        hIndex = hIndex + 1;
                        termReal = hReal * dwPointer[dwOffset] - hImag * dwPointer[dwOffset + 1];
                        smatRPointer[i - 1] += termReal;
                        termImag = hReal * dwPointer[dwOffset + 1] + hImag * dwPointer[dwOffset];
                        smatIPointer[i - 1] += termImag;


                        //
                        // COMPUTE INDICATION OF LOSS OF SIGNIFICANCE:
                        abs1Pointer[hIndex - 1] += std::fabs(termReal) + std::fabs(termImag);
                    }
                }

                }  // end else (fiR/fiI large enough)
            }

            //
            // END OF INTEGRATION
            //
            // dwRioC=dwRioC/DWCONT
            if (dwAmpCount > 0.0) dwRioC = dwRioC / dwAmpCount;
            if (debugSwitch) {
                std::printf(" SMALLEST, LARGEST AND AVERAGE std::fabs(DWRIRO)= %15.5G%15.5G%15.5G\n\n",
                            dwAmpMin, dwAmpMax, dwRioC);
            }

            if (isInfoPrint) {
            lineCount = lineCount + nMloLx + 1;
            if (lineCount > 58) {
            if (lineCount < 1000) {
                printInelasticHeader(reactStr, header, eLab);
            }
            //   T55, 'COULOMB', T70, 'NUCLEAR PART', T94, 'TOTAL AMPLITUDE' /
            //   '  IN OUT', 11X, 'MAG.', 6X, 'PHASE', 4X, 'CANCELLATION',
            //   T54, 'AMPLITUDE',
            //   T68, 2( 'MAG.', 6X, 'PHASE', 9X ) / )
            std::printf("  L   L  LX%14sINTEGRAL(0, SUMMAX)%29sCOULOMB%12sNUCLEAR PART%21sTOTAL AMPLITUDE\n",
                        "", "", "", "");
            std::printf("  IN OUT%11sMAG.      PHASE    CANCELLATION%13sAMPLITUDE%7sMAG.      PHASE         MAG.      PHASE\n",
                        "", "", "");
            lineCount = 7 + nMloLx;
            }  // end if (lineCount > 58)
            }  // end if (isInfoPrint)

            //
            // ADD ON PIECE FROM sumMax TO INFINITY,
            // MULTIPLY IN clebschGordan GORDEN factor,
            // COMPUTE CANCELLATIONS,
            // AND SEPARATE INTO NUCLEAR AND COULOMB PIECES
            // ALSO PUT IN THE i**(-lx-1)
            //
            { double* smatInPointer = reaction.distortedWave.channel[1].smatArr.data();  // 0-based class-owned
              sIn = std::complex<double>(smatInPointer[2*li], smatInPointer[2*li + 1]); }
            hIndex = 0;
            for (lx = lxMin; lx <= lxMax; lx += 2) {
                // loMn = std::abs(li-lx)
                loMn = std::abs(li - lx);
                // loMx = li+lx
                loMx = li + lx;
                // uses bratPointer
                betaRatio = bratPointer[lx/2];
                // Fortran: iand(lx,3) gives lx mod 4 (0-based), +1 gives 1-based index
                phase = phases[lx & 3];  // C++ 0-based: phases[0..3]
                // K / KBASE / LDELMN decode (see computeKBase lambda above)
                computeKBase();
                for (lo = loMn; lo <= loMx; lo += 2) {
                    // kOffset = KBASE + (lo-li-LDELMN)/2
                    kOffset = kBase + (lo - li - lDeltaMin) / 2;
                    // hIndex = hIndex + 1
                    hIndex = hIndex + 1;
                    // i1 = (liIndex-1)*nSpl + kOffset
                    i1 = (liIndex - 1) * nSpl + kOffset;
                    // i2 = (li-lMin)*nSpl + kOffset
                    i2 = (li - lMin) * nSpl + kOffset;
                    // C = factor *  std::fabs( clebschGordan( 2*li, 2*lx, 0, 0, 2*lo, 0 ) )
                    coefficient = factor * std::fabs(clebschGordan(2*li, 2*lx, 0, 0, 2*lo, 0));
                    compAmp = coefficient * std::complex<double>(smatRPointer[i1 - 1], smatIPointer[i1 - 1]);
                    { double* smatOutPointer = reaction.distortedWave.channel[2].smatArr.data();  // 0-based class-owned
                      sOut = std::complex<double>(smatOutPointer[2*lo], smatOutPointer[2*lo + 1]); }
                    // 1-based pointers into class-owned CL1*_arr.
                    { double* cl1ffPointer = reaction.inelastic.cl1ffArr.data(); double* cl1ggPointer = reaction.inelastic.cl1ggArr.data();  // 0-based (accessed [i1-1])
                      double* cl1fgPointer = reaction.inelastic.cl1fgArr.data(); double* cl1gfPointer = reaction.inelastic.cl1gfArr.data();
                    rToInAmp = 0.25 * coefficient * betaRatio * (
                        (1.0 + sOut) * (1.0 + sIn) * cl1ffPointer[i1 - 1]
                        - (1.0 - sOut) * (1.0 - sIn) * cl1ggPointer[i1 - 1]
                        + imagUnit * ((1.0 + sOut) * (1.0 - sIn) * cl1fgPointer[i1 - 1]
                                   + (1.0 - sOut) * (1.0 + sIn) * cl1gfPointer[i1 - 1])
                    );
                    // cl2ffArr is 0-based: [I2] == cl2ffArr[I2-1].
                    ffI = coefficient * betaRatio * reaction.inelastic.cl2ffArr[i2 - 1]; }
                    // totalAmp = compAmp + rToInAmp
                    totalAmp = compAmp + rToInAmp;
                    // nuclearAmp = totalAmp - ffI
                    nuclearAmp = totalAmp - ffI;
                    // nuclearAmp = PHASE*nuclearAmp
                    nuclearAmp = phase * nuclearAmp;
                    smatRPointer[i1 - 1] = nuclearAmp.real();
                    smatIPointer[i1 - 1] = nuclearAmp.imag();
                    if (isInfoPrint) {
                    // compAmp = PHASE*compAmp
                    compAmp  = phase * compAmp;
                    // totalAmp = PHASE*totalAmp
                    totalAmp = phase * totalAmp;
                    ampCoulomb = std::abs(compAmp);
                    absRatio  = (ampCoulomb > 0.0) ? coefficient * abs1Pointer[hIndex - 1] / ampCoulomb : 0.0;
                    // ampNuclear = CDABS(nuclearAmp)
                    ampNuclear = std::abs(nuclearAmp);
                    // ampTotal = CDABS(totalAmp)
                    ampTotal = std::abs(totalAmp);
                    // phaseCoulomb = std::atan2( compAmp.imag(), compAmp.real() )
                    phaseCoulomb = std::atan2(compAmp.imag(), compAmp.real());
                    // phaseNuclear = std::atan2( nuclearAmp.imag(), nuclearAmp.real() )
                    phaseNuclear = std::atan2(nuclearAmp.imag(), nuclearAmp.real());
                    // phaseTotal = std::atan2( totalAmp.imag(), totalAmp.real() )
                    phaseTotal = std::atan2(totalAmp.imag(), totalAmp.real());

                    std::printf("%4d%4d%3d%16.4G%11.4G%9.2f%16.4G%13.4G%11.4G%13.4G%11.4G\n",
                                li, lo, lx, ampCoulomb, phaseCoulomb,
                                absRatio, ffI, ampNuclear, phaseNuclear, ampTotal, phaseTotal);

                    }  // end if (isInfoPrint)
                }
            }

            if (isInfoPrint) std::printf(" \n");

            // (end of li loop body)
        }  // end DO 959 li=liMin,lMax,2

        //
        // END OF THE li LOOP FOR A GIVEN parity OF li
        // IF li IS EVEN HAVE ANOTHER BASH AT li LOOP WITH ODD li
        //
        // liMin=liL (end of liParity loop body)
        liMin = liL;

    }  // end DO 989 liParity=1,2

    //
    // END OF THE 2 li LOOPS (EVEN AND ODD parity OF li)
    //
    // FREE WORK AREAS
    //
    if (((printLevel / 100) % (10)) < 1) {

    // are permanently 0 sentinels (nuclearHArr/coulombHArr class-owned vectors

    }  // end if (printLevel < 100) — pool freeing

    //
    // TIMING INFO
    //
    totalTime = second() - (double)tStart;
    waveTime = (double)times[1] + (double)times[2] + (double)times[3];
    recurTime1 = (double)times[7] - (double)times[5];
    recurTime2 = (double)times[8] - (double)times[6];
    otherTime = totalTime - waveTime;
    totalTime = totalTime + (double)times[7] + (double)times[8];

    //   '0SCATTERING WAVE INITIALIZATIONS:', T35, F9.3 /
    //   ' SCATTERING WAVE LOOP:', T35, F9.3 /
    //   ' SCATTERING WAVE INTERPOLATIONS:', T35, F9.3 /
    //   ' TOTAL SCATTERING WAVE TIME', T45, F9.3 /
    //   ' INTEGRATIONS AND CLEBSCH''S:', T45, F9.3 /
    //   ' BELLINGS (SUMMAX, INFINITY):', T35, F9.3 /
    //   ' RECURSION (sumMax, INFINITY):', T35, F9.3 /
    //   ' TOTAL COULOMB (sumMax, INFINITY):', T45, F9.3 /
    //   ' BELLINGS (0, INFINITY):', T35, F9.3 /
    //   ' RECURSION (0, INFINITY):', T35, F9.3 /
    //   ' TOTAL PURE COULOMB (0, INFINITY):', T45, F9.3 /
    //   ' TOTAL:', T55, F9.3, ' SEC.' / )
    std::printf("\n0RADIAL INTEGRAL COMPUTATION TIMES (SECONDS):\n");
    std::printf("0SCATTERING WAVE INITIALIZATIONS:%26s%9.3f\n", "", (double)times[1]);
    std::printf(" SCATTERING WAVE LOOP:%33s%9.3f\n",            "", (double)times[2]);
    std::printf(" SCATTERING WAVE INTERPOLATIONS:%28s%9.3f\n",  "", (double)times[3]);
    std::printf(" TOTAL SCATTERING WAVE TIME%38s%9.3f\n",       "", waveTime);
    std::printf(" INTEGRATIONS AND CLEBSCH'S:%38s%9.3f\n",      "", otherTime);
    std::printf(" BELLINGS (SUMMAX, INFINITY):%26s%9.3f\n",     "", (double)times[5]);
    std::printf(" RECURSION (sumMax, INFINITY):%25s%9.3f\n",    "", recurTime1);
    std::printf(" TOTAL COULOMB (sumMax, INFINITY):%21s%9.3f\n","", (double)times[7]);
    std::printf(" BELLINGS (0, INFINITY):%31s%9.3f\n",          "", (double)times[6]);
    std::printf(" RECURSION (0, INFINITY):%30s%9.3f\n",         "", recurTime2);
    std::printf(" TOTAL PURE COULOMB (0, INFINITY):%21s%9.3f\n","", (double)times[8]);
    std::printf(" TOTAL:%49s%9.3f SEC.\n\n",                    "", totalTime);

    //
    // A L L    D O N E
    //
}

// ---------------------------------------------------------------------------
// Backward-compat free function wrapper — calls the class method
// ---------------------------------------------------------------------------
