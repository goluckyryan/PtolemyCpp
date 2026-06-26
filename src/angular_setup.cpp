// angular_setup.cpp — ANGSET: sets up arrays for the angular transforms in the
// radial integrals. Called after WAVPOT and before INELDC; allocates arrays used
// by A12. returnCode = 0 on error, 1 otherwise.

#include "ptolemy_types.h"
#include "CrossSectionCalc.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "Reaction.h"
#include "Constants.h"

void CrossSectionCalc::angleSet()
{
    // =====================================================================
    // Local variables
    // =====================================================================
    int verbosity;
    int liExtraCount, i;
    int liIndex;
    int xlamSize;
    int jpMin, jpMax, jPi, jPo;
    int loMax, loMin;
    int lo;

    logical debugSwitch;

    // CHARACTER*8 DATA initializations
    // NAMER/NAMEI for WAVSAV* pool slots deleted. The per-channel



    // =====================================================================
    // =====================================================================
    constexpr int facFr4 = 2;

    auto& lMax      = reaction_.angMom.lMax;
    auto& lMin      = reaction_.angMom.lMin;
    auto& lStep     = reaction_.integrationGrid.lStep;


    auto& stripPickup    = reaction_.internalState.stripPickup;

    auto& lBoundProj       = reaction_.boundState.vertex[1].lBound;
    auto& lBoundTarg       = reaction_.boundState.vertex[2].lBound;

    auto& nSmatPerL     = reaction_.inelastic.nSmatPerL;
    auto& nAspli    = reaction_.inelastic.nAspli;
    auto& lxMin     = reaction_.inelastic.lxMin;
    auto& lxMax     = reaction_.inelastic.lxMax;
    auto& mStop     = reaction_.inelastic.mStop;
    auto& nMloLx    = reaction_.inelastic.nMloLx;
    auto& nLValues    = reaction_.inelastic.nLValues;
    auto& lSkip     = reaction_.inelastic.lSkip;
    auto& liFit     = reaction_.inelastic.liFit;
    auto& liFitIndex    = reaction_.inelastic.liFitIndex;
    auto& a12mSize     = reaction_.inelastic.a12mSize;

    // indxsSize was an alias of reaction_.iSize — actually local-in-disguise:
    // angular_setup writes indxsSize=0, then SETSPT fills it with the indxs
    // count, then it sizes inelastic.indxsArr. Repurposing the pool-size
    // field as a local counter; now a plain local int. (Field dropped
    int indxsSize;

    auto& a12nSize     = reaction_.gridData.a12nSize;
    auto& noFlo     = reaction_.gridData.noFlo;
    auto& nRiRoInterp    = reaction_.gridData.nRiRoInterp;
    auto& nCrit     = reaction_.gridData.nCrit;
    auto& nWfi      = reaction_.gridData.nWfi;
    auto& nWfo      = reaction_.gridData.nWfo;



    // =====================================================================
    // BEGIN EXECUTABLE CODE
    // =====================================================================
    // ever wrote returnCode=1 at the end.
//
    verbosity = ((reaction_.flags.printLevel) % (10));
    debugSwitch = verbosity >= 3;
//
//
//
//
//     WE DO NOT NECESSARILY COMPUTE EVERY i( li, lo, lx ).  RATHER
//     WE FIND ONLY THOSE FOR
//        li = lMin, lMin+lStep, lMin+2LSTEP, ..., lMax
//     AND
//        li = lMin, lMin+1, ..., lxMax
//     WHERE THE SERIES TERMINATES AT OR BEFORE lMax.  ALL lo AND lx
//     FOR THE ABOVE li'S ARE FOUND.  THE ARRAY "LIS" CONTAINS THOSE
//     li'S FORWHICH THE CALCULATION IS TO BE DONE
//     NOTE THAT IF lMin < lxMax, WE EFFECTIVELY START OUT WITH
//     lStep=1 TO POVIDE A SMOOTH STARTING POINT FOR INTERPOLATION
//
//**************************************************
//
//     NOTE:  FOR COUPLED CHANNELS CALCULATIONS THESE ARE
//            REALLY THE TOTAL J VALUES.  HOWEVER WE DO
//            NOT YET DOUBLE THEM.
//
//**************************************************
//
//
    // warning were never reached on any test input.
//
//     COMPUTE THE NUMBER OF li'S INCLUDING THE li<lxMax PART
//
    nLValues = (lMax - lMin + lStep) / lStep;
    liExtraCount = 0;
    if (lxMax > lMin) liExtraCount = (lxMax - lMin) / lSkip
        - (lxMax - lMin) / lStep;
    nLValues = nLValues + liExtraCount;
    // lisArr: 1-based, indices [1..nLValues] valid.
    reaction_.inelastic.lisArr.assign(nLValues + 1, 0);
    int* lisPointer = reaction_.inelastic.lisArr.data();
    for (i = 1; i <= nLValues; i++) {
        int l = lMin + (i - 1) * lSkip;
        if (l > lxMax) l = (i - liExtraCount - 1) * lStep + lMin;
        lisPointer[i] = l;
    } // 309
//
//
//     DETERMINE THE RANGE OF li'S TO BE USED FOR EXTRAPOLATION.
//
    // lisPointer already points at reaction_.inelastic.lisArr.data() (1-based)
    // LBACK was a Reaction struct field bulk-init'd to NOTDEF_INT by DEFALT; the
    // `if (LBACK == NOTDEF_INT)` branch was always taken in this port (no user-input
    int lBackoff = lMax - lisPointer[nLValues - 3];
    liFit = lMax - lBackoff;
    for (liIndex = 1; liIndex <= nLValues; liIndex++) {
        if (liFit < lisPointer[liIndex]) break;
    }
    liFitIndex = liIndex - 1;
    liFit = lisPointer[liFitIndex];
    reaction_.inelastic.nolFit = nLValues - liFitIndex + 1;
//
//
//     "INDXS" WILL CONTAIN ALL INFORMATION NECESSARY TO FIND A
//     GIVEN S-MATRIX ELEMENT.  IT HAS 3 VALUES FOR EACH
//     (lx, jProj, jT) TRIPLE, INDEXED BY
//        K = 1 + lx + nLx*( jProj-jpBase + nJp*(jT-jtBase) )/2
//
//        (3*K-2) = kOffset = OFFSET INTO S-MATRIX (=0 IF NONEXISTENT).
//
//
//     "TOCS" WILL CONTAIN A TABLE OF CONTENTS FOR THE S MATRIX.
//     IT HAS 4 VALUES IDENTIFYING EACH S-MATRIX ELEMENT USED, INDEXED
//     BY ITS POSITION (kOffset + (lo-li+LDELMN)/2 ) IN THE S MATRIX:
//
//        (4*kOffset-3) = lo - li
//        (4*kOffset-2) = lx
//        (4*kOffset-1) = jProj
//        (4*kOffset  ) = jT
//
//
//     "SMATR", "SMATI" WILL HOLD THE CALCULATED (NOT INTERPOLATED)
//     REACTION S MATRIX.  FOR A GIVEN (lo, lx, jProj, jT, li), THE
//     S MATRIX IS INDEXED AS
//
//        (li-lMin)*n_spl/lStep + kOffset + ( lo - li - LDELMN )/2
//
//     WHERE kOffset AND LDELMN ARE TAKED FROM "INDXS" AS DESCRIBED ABOVE.
//
//
//     FIRST PASS, ALLOCATE "INDXS" AND "TOCS" AND INITIALIZE POINTERS.
//     FIRST FIND INDXS size
//
    indxsSize = 0;
    reaction_.setupInelasticAngMomTable(indxsSize, nullptr, nullptr, 1);
//
    reaction_.inelastic.indxsArr.assign(indxsSize + 1, 0);          // +1 tail pad; 0-based via indxsPointer[0..indxsSize-1]
    reaction_.inelastic.indxsPointer = reaction_.inelastic.indxsArr.data();
//     nAspli IS n_spl SUMMED OVER ALL CHANNELS
    nAspli = 0;
//
//     NOW FILL IN "INDXS" FOR ALL S-MATRIX ELEMENTS NEEDED ON THIS
//     AND SUBSEQUENT PASSES.
//
    reaction_.setupInelasticAngMomTable(nAspli, reaction_.inelastic.indxsPointer, nullptr, 2);
//
//     ALLOCATE AND INITIALIZE THE S-MATRIX AND ITS TABLE OF CONTENTS.
//     DELTASR, DELTASI  ARE FOR S - S(BORN).  WE START OUT WITH
//     THEM DEFINED AS THE FULL S-MATRIX ARRAYS.
//
    nSmatPerL = nAspli * nLValues;
    // smatrArr/smatiArr: +1 element for a 1-based safety pad.
    reaction_.inelastic.smatrArr.assign(nSmatPerL + 1, 0.0);
    reaction_.inelastic.smatiArr.assign(nSmatPerL + 1, 0.0);
    reaction_.inelastic.smatRPointer = reaction_.inelastic.smatrArr.data();  // 0-based
    reaction_.inelastic.smatIPointer = reaction_.inelastic.smatiArr.data();  // 0-based
    reaction_.inelastic.tocsArr.assign(4 * nAspli + 1, -1);     // +1 safety pad; -1 sentinel marks "skip"
    reaction_.inelastic.tocsPointer = reaction_.inelastic.tocsArr.data() - 1;
//
//
//     ALL PASSES:  FILL IN "TOCS" FOR THOSE S-MATRIX ELEMENTS WHICH
//     ARE USED ON THIS PASS.
//
    {
        int dummyCounter = 0;
        reaction_.setupInelasticAngMomTable(dummyCounter, reaction_.inelastic.indxsPointer, reaction_.inelastic.tocsPointer, 3);
    }
//
//     THIS IS THE NUMBER OF (lo', lx') PAIRS PER li IN THE INNER
//     INTEGRATION LOOPS.  lx'S DUE ONLY TO SPIN-ORBIT DISTORTIONS
//     ARE NOT INCLUDED.
//
    nMloLx = ((lxMax + lxMin + 1) * (lxMax - lxMin + 1)
              + ((lxMin + lBoundProj + lBoundTarg + 1) % (2))) / 2;
//
//
//     STORAGE FOR THE H'S (TRANSFER ONLY)
//
    if (stripPickup != 0) {
//
    // The full vector is sized below once a12nSize is known.
    reaction_.dwbaGrid.allocateHint(nMloLx, reaction_);
    reaction_.dwbaGrid.allocateHabs(nMloLx, reaction_);
//
//     ALLOCATE ARRAYS FOR  A12
//
    mStop = std::max({lxMax + lBoundProj, lxMax + lBoundTarg, lBoundProj + lBoundTarg + 1, 2 * lxMax});
    mStop = 2 * (mStop / 2);
    xlamSize = (mStop / 2 + 1) * (reaction_.kin.lOutMax + 1);
    // Size lOutMax+2 covers A12's nLam[1..lOutMax+1] access (A12 does nLam=arg-1).
    reaction_.dwbaGrid.allocateNlam(reaction_.kin.lOutMax + 2);  // INLAM = 0
    reaction_.dwbaGrid.allocateXlam(xlamSize + 1);  // IXLAM = 0; +1 covers A12's 1-based [1..xlamSize]
//
//     FOLLOWING IS AN UPPER LIMIT ON THE NUMBER OF TERMS IN THE
//     A12 SUMS.
//     FIRST THE MAXIMUM NUMBER OF M1 AND M2
//     COMBINATIONS.  WE USE
//        M1  <= (lBoundProj, lx+lBoundTarg)    M2  <= (lBoundTarg, lx+lBoundProj)
//        mX  <= (lx,lBoundProj+lBoundTarg)   M1+lBoundProj AND M2+lBoundTarg  EVEN
//
    a12mSize = (std::min(lxMax, lBoundProj + lBoundTarg) + 1) * (std::min(lBoundProj, lBoundTarg) + 1);
//
//     THEN FOR EACH OF THE ABOVE MU TAKES ON VALUES SUCH THAT
//       0 <= MU <= lo <= lOutMax,    MU+lo EVEN
//     AND WE HAVE  NMLOLX  SETS OF SUCH A12 VALUES.
//
    a12nSize = a12mSize * (reaction_.kin.lOutMax / 2 + 1) * nMloLx;
//
//
//     FOR EACH (M1, M2) PAIR WE HAVE UP TO FIVE ITEMS TO STORE
//
    a12mSize = 5 * a12mSize;
    // hsA12Arr layout: H at [0..NMLOLX-1], A12 at [NMLOLX..NMLOLX+a12nSize-1].
    reaction_.gridData.hsA12Arr.assign(nMloLx + a12nSize, 0.0);
    reaction_.gridData.msvalArr.assign(a12mSize, 0.0);
    reaction_.gridData.ja12sArr.assign(((a12mSize + 1) / facFr4) * facFr4, 0);
//
//     FOLLOWING IS 3 ARRAYS; ONLY THE THIRD IS USED FOR USEHS RUNS.
//     REAL*8   numLx
//     REAL*8   2*lxMax+2
//     INTEGER  3*numLx  ( WE TREAT AS REAL*8 TO ALLOW FOR CDC )
//
    int a12tmSize = 4 * reaction_.inelastic.numLx + 2 * lxMax + 2;
    reaction_.dwbaGrid.a12tm_.assign(a12tmSize, 0.0);
//
//     "iDwfi" WILL CONTAIN INFORMATION ABOUT THE INCIDENT
//     SCATTERING WAVE FUNCTIONS, FOR A GIVEN li:
//
//        (3*K-2) = LASI - li
//        (3*K-1) = jPi - 2*li  (MAX. VALUE IF NO S. O.)
//        (3*K  ) = INDEX TO BASE OF WF IN ALLOCATOR
//
    } // end if (stripPickup != 0)
    reaction_.dwbaGrid.allocateIwfii((3 * reaction_.distortedWave.channel[1].nJStates + 1) / facFr4, reaction_);
    // channel[1].TCSWS permanently 0 so LASMIN=LASMAX=0; the asymptotic-L
    // axis degenerates to a single point and the body's MIN0/std::max(2*LAS,0)
    // reduce to 0 with LAS=0.
    nWfi = 0;
//     THE INNER LOOP IS OVER jProj = TOTAL PROJECTILE J
    jpMax = reaction_.distortedWave.channel[1].twoSpin;
    jpMin = reaction_.distortedWave.channel[1].hasSpinorbit ? -reaction_.distortedWave.channel[1].twoSpin : jpMax;
    for (jPi = jpMin; jPi <= jpMax; jPi += 2) { // 1809
        nWfi = nWfi + 1;
        reaction_.gridData.iDwfiPointer[3 * nWfi - 2] = 0;  // LAS, permanently 0
        reaction_.gridData.iDwfiPointer[3 * nWfi - 1] = jPi;
        reaction_.gridData.iDwfiPointer[3 * nWfi]     = (nWfi - 1) * nRiRoInterp;
    } // 1809
    i = 3 * nWfi;
    if (debugSwitch) {
        for (int k = 1; k <= i; k++)
            std::printf(" iDwfi%10d\n", reaction_.gridData.iDwfiPointer[k]);
    }
//
//
//     "iDwfo" WILL CONTAIN INFORMATION ABOUT THE OUTGOING
//     SCATTERING WAVE FUNCTIONS, FOR A GIVEN li:
//
//        (4*K-3) = lo - li
//        (4*K-2) = LASO - li
//        (4*K-1) = JPO - 2*li
//        (4*K  ) = INDEX TO BASE OF WF ARRAY (SET IN INELDC).
//
    i = (stripPickup != 0) ? ((lBoundProj + lBoundTarg + lxMax) % (2)) : 0;
    noFlo = lxMax + 1 - i;
    reaction_.dwbaGrid.allocateIwfio(4 * noFlo * reaction_.distortedWave.channel[2].nJStates / facFr4, reaction_);
    loMax = noFlo - 1;
    loMin = -loMax;
    nWfo = 0;
    for (lo = loMin; lo <= loMax; lo += 2) { // 1869
        jpMax = 2 * lo + reaction_.distortedWave.channel[2].twoSpin;
        jpMin = reaction_.distortedWave.channel[2].hasSpinorbit ? 2 * lo - reaction_.distortedWave.channel[2].twoSpin : jpMax;
        for (jPo = jpMin; jPo <= jpMax; jPo += 2) { // 1849
//     ELIMINATE INITIAL JPO'S TOO LOW TO COUPLE TO li.
            if (reaction_.distortedWave.channel[2].hasSpinorbit && nWfo == 0 && jPo + reaction_.distortedWave.channel[1].twoSpin + reaction_.boundState.vertex[2].jB < 0) continue; // goto 1849
            nWfo = nWfo + 1;
            reaction_.gridData.iDwfoPointer[4 * nWfo - 3] = lo;
            reaction_.gridData.iDwfoPointer[4 * nWfo - 2] = lo;   // was LAS = lo
            reaction_.gridData.iDwfoPointer[4 * nWfo - 1] = jPo;
        } // 1849
    } // 1869
//     ELIMINATE FINAL jProj'S TOO LARGE TO COUPLE TO li.
    if (reaction_.distortedWave.channel[2].hasSpinorbit) {
    i = nWfo;
    for (int k = 1; k <= i; k++) { // 1879
        if (reaction_.gridData.iDwfoPointer[4 * nWfo - 1] <= reaction_.distortedWave.channel[1].twoSpin + reaction_.boundState.vertex[2].jB) break;
        nWfo = nWfo - 1;
    } // 1879
    if (nWfo < noFlo * reaction_.distortedWave.channel[2].nJStates) {
        reaction_.dwbaGrid.resizeIwfio(4 * nWfo / facFr4, reaction_);
        // reaction_.gridData.iDwfoPointer re-synced by resizeIwfio
    }
    } // end if (reaction_.distortedWave.channel[2].hasSpinorbit)
    i = 4 * nWfo;
    if (debugSwitch) {
        for (int k = 1; k <= i; k++)
            std::printf(" iDwfo%10d\n", reaction_.gridData.iDwfoPointer[k]);
    }
//
//     NOW ALLOCATE SPACE FOR THE WAVE FUNCTIONS (AT GAUSS POINTS)
//
    // ILIR/ILII use DWBAGrid class member vectors (WAVEAR/WAVEAI → reaction_.dwbaGrid.lir/lii)
    i = (nWfi * nRiRoInterp + 1) / facFr4;
    reaction_.dwbaGrid.allocateLir(i, reaction_);  // sets lirPointer=lir.data()
    reaction_.dwbaGrid.allocateLii(i, reaction_);  // sets liiPointer=lii.data()
    // ILOR/ILOI use DWBAGrid class member vectors (WAVEBR/WAVEBI → reaction_.dwbaGrid.lor_/loi)
    // Bug fix: use MAX(nWfo,noFlo)*nRiRoInterp since inelasticRadialIntegrals fills noFlo L-channels in its first-li loop
    // but nWfo may be trimmed below noFlo, making the buffer too small.
    i = (std::max(nWfo, noFlo) * nRiRoInterp + 1) / facFr4;

    reaction_.dwbaGrid.allocateLor(i, reaction_);  // sets lorPointer=lor_.data()
    reaction_.dwbaGrid.allocateLoi(i, reaction_);  // sets loiPointer=loi.data()
    // lirPointer/liiPointer/lorPointer/loiPointer set by class allocators (0-based float*).
    // ILIR/ILII/ILOR/ILOI int fields removed (always 0 sentinels).
//     THIS POINTER CONTROLS WF STORAGE IN THE CIRCULAR BUFFER.
//     When noFlo > nWfo (spin-orbit trimming), the circular buffer must be sized
//     to noFlo slots (to hold the first-li fill), so nCrit must match.
    nCrit = std::max(nWfo, noFlo) * nRiRoInterp - nRiRoInterp;
//
//     "DW" WILL CONTAIN PRODUCTS OF SCATTERING WAVE FUNCTIONS.
//     "indxDw" WILL CONTAIN POINTERS USED TO CALCULATE AND USE "DW".
//
    // DW → DWBAGrid::dw_ class member vector (dwPointer=dw_.data(), 0-based);
    reaction_.dwbaGrid.allocateDw(nWfi * nWfo * 2, reaction_);
    reaction_.dwbaGrid.allocateDwi(nWfi * nWfo * 4 / facFr4, reaction_);
    // COUPSW CC LCHNDF pool tagging (no outgoing wave arrays) dropped

    // --- ANGSET lines 1545-1627 ---

    // Compute max number of H's (local to angleSet — was gridData.NUMHS)
    const int hCount = reaction_.inelastic.nMloLx;

    i = hCount * reaction_.distortedWave.channel[1].nJStates * reaction_.distortedWave.channel[2].nJStates;

    // "RIROABS" for round-off error checking
    reaction_.dwbaGrid.allocateAbs1(i, reaction_);  // sets abs1Pointer=abs1.data() (0-based)

    if (stripPickup != 0) {

    // "SUMHVALS" and "SUMIVALS" — H's before/after interpolation
    reaction_.dwbaGrid.allocateSmhvl(reaction_.gridData.nPhiSum * hCount, reaction_);
    // stores a 0-based offset into iIndex(1, numIi); reader (inelastic_dwba) indexes smivlArr.
    reaction_.gridData.smivlArr.assign(reaction_.gridData.nInterpPoints * hCount, 0.0);
    reaction_.gridData.smivlPointer = reaction_.gridData.smivlArr.data() - 1;  // 1-based

    // "I1REAL", "I1IMAG" — final 3D integrals
    // "iIndex" — pointers into H and DW
    reaction_.dwbaGrid.allocateLiloR(i, reaction_);
    reaction_.dwbaGrid.allocateLiloI(i, reaction_);
    // liloRPointer/liloIPointer set by allocate methods
    reaction_.dwbaGrid.allocateIiindx(4 * i / facFr4, reaction_);

    // Setup for high speed cosine
    {
        constexpr int nCosin = 256;
        constexpr double diCoss = (double)nCosin;
        reaction_.gridData.cosStep = 2.0 * Constants::PI / diCoss;
        reaction_.gridData.cosinQuarter = nCosin / 4;
        reaction_.gridData.stepInverse = diCoss / (2.0 * Constants::PI);
        int cosinSize = nCosin + 1;
        // 0-based cos table
        reaction_.dwbaGrid.allocateCosin(cosinSize, reaction_);
        double* cosinPointer = reaction_.dwbaGrid.cosin.data();  // 0-based for fill (accessed [ii-1])
        double cosArg = -reaction_.gridData.cosStep;
        for (int ii = 1; ii <= cosinSize; ii++) {
            cosArg = cosArg + reaction_.gridData.cosStep;
            cosinPointer[ii - 1] = std::cos(cosArg);  // fills cosin[0..cosinSize-1]
        }
    }

    // Table of integers in double precision form
    {
        i = reaction_.angMom.lMax + reaction_.inelastic.lxMax + reaction_.boundState.vertex[1].lBound + reaction_.boundState.vertex[2].lBound;
        int intsMin = std::max(i, 2 * (reaction_.boundState.vertex[1].lBound + reaction_.boundState.vertex[2].lBound) + 20);
        reaction_.gridData.intsOffset = intsMin + 1;
        int intsMax = std::max(intsMin, 2 * i + 2);
        reaction_.gridData.intsArr.assign(intsMax + intsMin + 1, 0.0);
        double* intsPointer = reaction_.gridData.intsArr.data() + reaction_.gridData.intsOffset;  // 0-based (accessed [i-1])
        intsMin = -intsMin;
        double iDouble = intsMin;
        for (i = intsMin; i <= intsMax; i++) {
            intsPointer[i - 1] = iDouble;
            iDouble = iDouble + 1.0;
        }
    }

    } // end if (stripPickup != 0)
}

// ---------------------------------------------------------------------------
// Backward-compat free function wrapper — calls the class method
// ---------------------------------------------------------------------------
