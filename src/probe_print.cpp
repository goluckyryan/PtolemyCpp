// probe_print.cpp — PRBPRT: consistency checks + prints the problem summary.

#include "math/angular_momentum_coeff.h"
#include "Reaction.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "Constants.h"

// aMs(r) — 0-based pointer into the massesArr std::array<double, 5>
// One of three TUs with positional aMs(r)[i] callers
// (others: set_channels.cpp, parameters.cpp); each anon-namespace copy is
// independent, so this one is 0-based while the others stay 1-based.
namespace {
double* aMs(Reaction& r) { return r.masses.massesArr.data(); }  // 0-based
// Print character c, n times — replaces repeated `for(_i...) printf` rule/box lines.
void printRepeatedChar(char c, int n) { for (int i = 0; i < n; i++) std::printf("%c", c); }
// Standard "invalid bound-state angular momenta" error header (2 printf lines,
// byte-identical at the lxMax<lxMin and jxMax<jxMin checks). The differing
// "RESULTS IN ..." follow-up line stays at each call site.
void printInvalidBoundStateMomenta(int jBp, int jBt, int lBoundProj, int lBoundTarg) {
    std::printf("\n**** INVALID BOUND STATE ANGULAR MOMENTA:\n"
                " JP, JT, LP, LT =%5d/2%5d/2%5d%5d\n",
                jBp, jBt, lBoundProj, lBoundTarg);
}
}

void Reaction::probePrint(int& returnCode)
{

    // =========================================================================
    // Local variable declarations
    // =========================================================================

    // Static DATA arrays (initialized once)
    static const int ind[6] = { 0, 3, 1, 2, 4, 5 };  // 1-based

    static const char particleNames[6][9] = {
        "",
        "A       ", "B       ", "BIGA    ", "BIGB    ", "X       "
    };

    static const char xWords[3][9] = {
        "",
        "        ", ";   X   "
    };

    static const char excitationWords[5][9] = {
        "",
        "PROJECTI", "LE      ", "TARGET  ", "        "
    };

    static const char betaNames[5][9] = {
        "",
        "BETA    ", "BETACOUL", "BETARATS", "BELX    "
    };

    static const char parityWord[4][9] = {
        "",
        "     -1 ", "UNKNOWN ", "     +1 "
    };

    int printSwitch;
    int hasBeta[5] = {0,0,0,0,0}; // LOGICAL, 1-based
    int isTransferReaction;

    // Local scalars
    int i, ii, verbosity;
    int jBpMin, jBpMax, jBtMin, jBtMax;
    int lxParity;
    int jIn, jOut;
    // lMinInput is set only in the else (estimate-lMin-from-lCrit) branch
    // at line 657 (lMinInput = lMin). In the LVALUES branch (NAMLOC return
    // != 0), it stays uninitialized — but lMin there comes from the
    // LVALUES pool array so the `if (lMinInput == NOTDEF_INT)` rounding at
    int lMinInput = 0;
    int hasLvalues = 0, lPrevious, len;
    int lParityBase = 0;
    int numLxI;
    int jxMin, jxMax;
    double beta, charge, aTerm, temp;
    double r2Mass, coulombCoupling;
    // r2s[1..4] now lives in this->inelastic.r2s.
    double* r2s = this->inelastic.r2s;  // 1-based, [0] unused
    int mProj;
    double* betaPointer  = nullptr;
    double* betaCPointer = nullptr;
    double* betarPointer = nullptr;
    double* belxPointer  = nullptr;
    double* betnrPointer = nullptr;
    int phaseParity, numParticles;

    // =========================================================================
    // =========================================================================

    // aMs(i): 0-based pointer into masses.massesArr (mass_a, mass_b, mass_A, mass_B, mass_x)
    // aMsPointer[0]=mass_a, aMsPointer[1]=mass_b, etc. (callers shift their 1-based index by -1).
    double* aMsPointer = aMs(*this);

    // js index map (1-based): js[1]=JA, js[2]=JB, js[3]=jResidual, js[4]=jBigB,
    // js[5]=JX. Stored as angMom.js[6] (double); for integer comparisons cast to
    // int: (int)this->angMom.js[i] != NOTDEF_INT.

    // zArray: this->charges.zArray[i], 1-based (zArray[1]..zArray[5])
    int* zArray = this->charges.zArray;

    // parities: this->angMom.parities[i], 1-based
    int* parities = this->angMom.parities;



    // r2s[1]=real radius, r2s[2]=imag radius, r2s[3]=Coulomb radius, r2s[4]=VC.

    // lCrits: this->kin.lCrits[1], this->kin.lCrits[2]
    // lnkAd2: this->linkuleData.lnkAd2[i][J][K], 1-based

    // Convenient references to frequently used COMMON members
    double& undefValue     = this->internalState.undefValue;
    int&    stripPickup    = this->internalState.stripPickup;
    int&    iExcit    = this->internalState.iExcit;
    double& r0Mass    = this->internalState.r0Mass;
    double& eLab      = this->energies.eLab;
    double& mass_a    = this->masses.massesArr[0];
    double& mass_b    = this->masses.massesArr[1];
    double& mass_A    = this->masses.massesArr[2];
    double& mass_B    = this->masses.massesArr[3];
    double& massX     = this->masses.massesArr[4];
    double& R         = this->integrationGrid.R;
    double& rI        = this->opticalPotentialParams.rI;
    double& R0        = this->opticalPotentialParams.R0;
    double& rI0       = this->opticalPotentialParams.rI0;
    double& rC        = this->opticalPotentialParams.rC;
    double& rC0       = this->opticalPotentialParams.rC0;
    double& alMnMt    = this->opticalPotentialParams.alMnMt;
    double& alMxMt    = this->opticalPotentialParams.alMxMt;
    // tvReal, tvImag, taReal, taImag from this->distortedWave.scatteringSolver.potentialWork
    double& tvReal       = this->distortedWave.scatteringSolver.potentialWork.tvReal;
    double& tvImag       = this->distortedWave.scatteringSolver.potentialWork.tvImag;
    double& taReal       = this->distortedWave.scatteringSolver.potentialWork.taReal;
    double& taImag       = this->distortedWave.scatteringSolver.potentialWork.taImag;
    int&    nuConL    = this->flags.nuConL;
    int&    printLevel    = this->flags.printLevel;
    int&    lMin      = this->angMom.lMin;
    int&    lMax      = this->angMom.lMax;
    int&    lMinSub    = this->integrationGrid.lMinSub;
    int&    lMaxAdditional    = this->angMom.lMaxAdditional;
    int&    lStep     = this->integrationGrid.lStep;
    int&    L         = this->angMom.L;
    int&    lBoundProj       = this->boundState.vertex[1].lBound;
    int&    lBoundTarg       = this->boundState.vertex[2].lBound;
    int&    jBp       = this->boundState.vertex[1].jB;
    int&    jBt       = this->boundState.vertex[2].jB;
    int&    lCrit     = this->kin.lCrit;
    int&    lxMin     = this->inelastic.lxMin;
    int&    lxMax     = this->inelastic.lxMax;
    int&    lxStep    = this->inelastic.lxStep;
    int&    nLValues    = this->inelastic.nLValues;
    int&    jpMin     = this->inelastic.jpMin;
    int&    jpMax     = this->inelastic.jpMax;
    int&    jtMin     = this->inelastic.jtMin;
    int&    jtMax     = this->inelastic.jtMax;
    int&    lSkip     = this->inelastic.lSkip;
    // the `if (lx != notDefSentinel)` branch (this->lx has zero writers).
    int&    densitySwitch    = this->inelastic.densitySwitch;

    // js references (Fortran js(1)=JA, js(2)=JB, js(3)=jResidual, js(4)=jBigB, js(5)=JX)
    // In C++ this->angMom.js[1..5]
    double& ja    = this->angMom.js[1];  // js(1) = JA
    double& jb    = this->angMom.js[2];  // js(2) = JB
    double& jBigA = this->angMom.js[3];  // js(3) = jResidual
    double& jBigB = this->angMom.js[4];  // js(4) = jBigB
    double& jx    = this->angMom.js[5];  // js(5) = JX

    // Cast macros for half-integer spin values stored as doubles
    #define jAInt    ((int)ja)
    #define jBInt    ((int)jb)
    #define jBigA_INT ((int)jBigA)
    #define jBigB_INT ((int)jBigB)
    #define jXInt    ((int)jx)

    // exs: this->energies.exs[1..5], 1-based
    // lSkips, statsCode now per-channel: this->distortedWave.channel[ch].lSkips/statsCode

    // lnkAd2(i,J,K): this->linkuleData.lnkAd2[i][J][K], 1-based
    // Access pattern: this->linkuleData.lnkAd2[i][j][k]

    // =========================================================================
    //      returnCode = 1
    // =========================================================================
    returnCode = 1;

    //      verbosity = (( printLevel) % (10 ))
    verbosity = ((printLevel) % (10));
    printSwitch = (verbosity >= 1);

    //    1  '   CONSTRUCTION OF THE INTEGRATION GRIDS', T70,
    //    2  'MERELY TO CONCEIVE OF SUCH THINGS MAKES THEM APPEAR RIDICULOUS'
    //    3   /
    //    4  '0', 45A1, 'eLab =', F7.2, ' MEV', 5X, 65A1 / )
    //
    // FORMAT '1' = page break (just print content); '0' = blank line before
    if (printSwitch) {
        std::printf("1%47sP T O L E M Y\n", "");
        std::printf("   CONSTRUCTION OF THE INTEGRATION GRIDS%29sMERELY TO CONCEIVE OF SUCH THINGS MAKES THEM APPEAR RIDICULOUS\n", "");
        std::printf("0");
        for (int _i = 1; _i <= 45; _i++) std::printf("%c", this->reactStr[_i]);
        std::printf("ELAB =%7.2f MEV     ", eLab);
        for (int _i = 1; _i <= 65; _i++) std::printf("%c", this->header[_i]);
        std::printf("\n\n");
    }

    // =========================================================================
    // =========================================================================
    isTransferReaction = (stripPickup != 0);

    //      lxStep = 1
    lxStep = 1;
    if (!isTransferReaction) {
        //      lBoundProj = std::abs(JB-JA)/2
        //      lBoundTarg = std::abs(jBigB-jResidual)/2
        lxStep = 2;
        lBoundProj = std::abs(jBInt - jAInt) / 2;
        lBoundTarg = std::abs(jBigB_INT - jBigA_INT) / 2;
    }

    // =========================================================================
    // =========================================================================
    for (i = 1; i <= 4; i++) {
        if (aMsPointer[i - 1] == undefValue) {
            std::printf("\n***** ERROR:  M%.4s IS NOT DEFINED.\n", particleNames[i]);
            returnCode = 0;
        }
        if ((int)this->angMom.js[i] == NOTDEF_INT) {  // js stored as double, compared as int (Fortran quirk)
            std::printf("\n***** ERROR:  J%.4s IS NOT DEFINED.\n", particleNames[i]);
            returnCode = 0;
        }
        if (zArray[i] == NOTDEF_INT) {
            std::printf("\n***** ERROR:  Z%.4s IS NOT DEFINED.\n", particleNames[i]);
            returnCode = 0;
        }
    }

    if (massX == undefValue) massX = std::abs(mass_a - mass_b);

    if (std::abs(mass_a + mass_A - mass_b - mass_B) >= 0.3) {
        //    1  ' MBIGA(MA, MB)MBIGB =', 4G15.5 )
        std::printf("\n**** THE REACTION DOES NOT CONSERVE MASS:\n"
                    " MBIGA(MA, MB)MBIGB =%15.5G%15.5G%15.5G%15.5G\n",
                    mass_A, mass_a, mass_b, mass_B);
        returnCode = 0;
    }

    if (zArray[1] + zArray[3] != zArray[2] + zArray[4]) {
        //    1  ' ZBIGA(ZA, ZB)ZBIGB =', 4I8 )
        std::printf("\n**** THE REACTION DOES NOT CONSERVE CHARGE:\n"
                    " ZBIGA(ZA, ZB)ZBIGB =%8d%8d%8d%8d\n",
                    zArray[3], zArray[1], zArray[2], zArray[4]);
        returnCode = 0;
    }

    if (std::abs(massX - std::abs(mass_a - mass_b)) >= 0.3) {
        //    1  'OTHER FOUR MASSES.' )
        std::printf("\n**** THE EXCHANGED MASS IS INCONSISTANT WITH THE OTHER FOUR MASSES.\n");
        returnCode = 0;
    }

    if (massX <= 0.7) {
    if (isTransferReaction) {
        std::printf("\n**** THE EXCHANGED MASS MUST NOT BE 0 OR NEGATIVE:%15.5G\n", massX);
        returnCode = 0;
    } else if (zArray[5] != 0) {
        std::printf("\n**** CANNOT TRANSFER CHARGE WITHOUT MASS TRANSFER: mass_x = 0, ZX =%4d\n",
                    zArray[5]);
        returnCode = 0;
    }
    }

    {

    // =========================================================================
    //      FOR TRANSFER CANNOT HAVE NUCLEAR CORE CORRECTIONS FOR LINKULES
    // =========================================================================
    if (isTransferReaction &&
        (this->linkuleData.lnkAd2[3][1][1] != 0 ||
         this->linkuleData.lnkAd2[3][2][1] != 0) &&
        nuConL >= 3) {
        //    1  ' BE COMPUTED FOR LINKULE-COMPUTED POTENTIALS.' /
        //    2  6X, 'ONLY COULOMB CORE-CORE CORRECTIONS WILL BE USED.' )
        std::printf("\n**** WARNING - NUCLEAR CORE-CORE CORRECTIONS CANNOT"
                    " BE COMPUTED FOR LINKULE-COMPUTED POTENTIALS.\n"
                    "      ONLY COULOMB CORE-CORE CORRECTIONS WILL BE USED.\n");
        nuConL = 2;
    }

    // =========================================================================
    //      CHECK THAT JBP AND JBT ARE INSIDE ALLOWED BOUNDS
    // =========================================================================
    // 200  jBpMin=std::abs(JA-JB)
    //      jBpMax=JA+JB
    //      jBtMin=std::abs(jResidual-jBigB)
    //      jBtMax=jResidual+jBigB
    jBpMin = std::abs(jAInt - jBInt);
    jBpMax = jAInt + jBInt;
    jBtMin = std::abs(jBigA_INT - jBigB_INT);
    jBtMax = jBigA_INT + jBigB_INT;

    if (isTransferReaction) {
    // =========================================================================
    //      FOR TRANSFER, CHECK (JBP, JBT) PAIR
    // =========================================================================
    if (!(jBp >= jBpMin  &&  jBp <= jBpMax  &&
          jBt >= jBtMin  &&  jBt <= jBtMax)) {
        std::printf("\n***** ERROR IN INPUT:\n"
                    " JA, JB, JBIGA, JBIGB =%5d/2%5d/2%5d/2%5d/2\n",
                    jAInt, jBInt, jBigA_INT, jBigB_INT);
        std::printf(" RESULTS IN\n"
                    " jBpMin, jBpMax, jBtMin, jBtMax =%5d/2%5d/2%5d/2%5d/2\n",
                    jBpMin, jBpMax, jBtMin, jBtMax);
        returnCode = 0;

        if (jBp == NOTDEF_INT  ||  jBt == NOTDEF_INT) {
            std::printf(" THUS THERE IS MORE THAN ONE POSSIBILITY FOR"
                        " JP FOR THE PROJECTILE OR TARGET BOUND STATE.\n");
            std::printf(" IN SUCH AMBIGUOUS CASES JP MUST BE EXPLICITLY INDICATED.\n"
                        " WE PRESENTLY HAVE  JP(PROJ) = %3d/2,    JP(TARG) = %3d/2\n",
                        jBp, jBt);
            std::printf(" WHERE THE *** INDICATES THE AMBIGUOUS CASE.\n");
        } else {
            std::printf(" BUT  JBP, JBT =%5d/2%5d/2  AND ARE OUT OF BOUNDS.\n",
                        jBp, jBt);
        }
    }

    // =========================================================================
    //      Find lx RANGE FOR TRANSFER AND MAKE NECESSARY CHECKS
    // =========================================================================
    // 260  lxMin=std::abs(JBP-JBT)
    //      lxMax=JBP+JBT
    lxMin = std::abs(jBp - jBt);
    lxMax = jBp + jBt;
    if (lxMin < std::abs(2*(lBoundProj-lBoundTarg))) lxMin = std::abs(2*(lBoundProj-lBoundTarg));
    if (lxMax > 2*(lBoundProj+lBoundTarg)) lxMax = 2*(lBoundProj+lBoundTarg);

    //      lxMax=lxMax/2
    //      lxMin=lxMin/2
    lxMax = lxMax / 2;
    lxMin = lxMin / 2;
    lxParity = 1 - 2*((lBoundProj+lBoundTarg) % (2));
    jpMin = jBp;
    jpMax = jBp;
    jtMin = jBt;
    jtMax = jBt;

    if (lxMax < lxMin) {
        printInvalidBoundStateMomenta(jBp, jBt, lBoundProj, lBoundTarg);
        std::printf(" RESULTS IN  LXMIN, LXMAX =%8d%8d\n", lxMin, lxMax);
        returnCode = 0;
    }
    } else {
    // =========================================================================
    //      SETUP lx RANGE FOR INELASTIC SCATTERING
    // =========================================================================
    if (this->energies.exs[4] == this->energies.exs[3]) {
        // PROJECTILE EXCITATION
        iExcit = 1;
        lxMax = jBpMax;
        lxMin = jBpMin;
        jpMin = jBpMin;
        jpMax = jBpMax;
        jtMin = 0;
        jtMax = 0;
        jIn  = jAInt;
        jOut = jBInt;
        lxParity = parities[1] * parities[2];

        if (this->energies.exs[2] == this->energies.exs[1]) {
            std::printf("\n**** E*B OR E*BIGB MUST BE DEFINED\n");
            returnCode = 0;
        }
    } else {
        // TARGET EXCITATION
        iExcit = 2;
        lxMax = jBtMax;
        lxMin = jBtMin;
        jpMin = 0;
        jpMax = 0;
        jtMin = jBtMin;
        jtMax = jBtMax;
        jIn  = jBigA_INT;
        jOut = jBigB_INT;
        lxParity = parities[3] * parities[4];

        if (this->energies.exs[2] != this->energies.exs[1]) {
            std::printf("\n**** CANNOT SIMULTANEOUSLY EXCITE BOTH PROJECTILE AND TARGET.\n");
            returnCode = 0;
        }
    }
    if (((lxMin) % (2)) != 0) {
        //    1  ' INTEGER CHANGE IN NUCLEAR SPINS', 3I8 )
        std::printf("\n**** INELASTIC EXCITATION MUST INVOLVE AN"
                    " INTEGER CHANGE IN NUCLEAR SPINS%8d%8d%8d\n",
                    lxMin, lxMax, iExcit);
        returnCode = 0;
    }

    // 285  lxMin = lxMin/2
    //      lxMax = lxMax/2
    lxMin = lxMin / 2;
    lxMax = lxMax / 2;

    if (lxParity != 0) {
        //      i = (IPARIT+3)/2
        //      lxMin = lxMin + ((i+lxMin) % (2))
        //      lxMax = lxMax - ((i+lxMax) % (2))
        int parityShift = (lxParity + 3) / 2;
        lxMin = lxMin + ((parityShift + lxMin) % (2));
        lxMax = lxMax - ((parityShift + lxMax) % (2));
        this->boundState.vertex[iExcit].lBound = lxMin;
    }

    if (lxMax != lxMin) {
        if (lxParity == 0) {
            // warning: parities unknown, assume (-1)**|J(IN)-J(OUT)|
            std::printf("\n");
            printRepeatedChar('*', 120);
            std::printf("\n");
            std::printf(" **** WARNING, THE PARITIES OF THE STATES MUST BE KNOWN"
                        " TO UNIQUELY DETERMIN THE MULTIPOLARITY OF THE EXCITATION\n");
            std::printf(" **** PTOLEMY ASSUMES THAT THE CHANGE IN THE INTRINSIC PARITY"
                        " IS JUST  (-1)**|J(IN)-J(OUT)|.\n");
            std::printf(" ");
            printRepeatedChar('*', 120);
            std::printf("\n");
            lxMax = lxMax - ((lxMin + lxMax) % (2));
        } else if (lxMax <= lxMin) {
            // 295: parities and J-values incompatible
            std::printf("\n**** ERROR:  PARITIES AND J-VALUES RESULT IN"
                        " IMPOSSIBLE EXCITATION:%4d%4d/2%4d/2%4d%4d\n",
                        lxParity, jIn, jOut, lxMin, lxMax);
            returnCode = 0;
        }
    }

    } // end if (isTransferReaction) / else (inelastic)

    // ONELSW=(lx!=notDefSentinel) "user specified a single lx" branch dropped

    this->inelastic.numLx = lxMax - lxMin + 1;
    numLxI = (lxMax - lxMin) / 2 + 1;
    if (isTransferReaction) {
        // =====================================================================
        //      sum OVER JX
        // =====================================================================
        jxMax = 2*lBoundProj + jBp;
        if (jxMax > 2*lBoundTarg + jBt) jxMax = 2*lBoundTarg + jBt;
        jxMin = std::abs(2*lBoundProj - jBp);
        if (jxMin < std::abs(2*lBoundTarg - jBt)) jxMin = std::abs(2*lBoundTarg - jBt);

        if (jxMax < jxMin) {
            printInvalidBoundStateMomenta(jBp, jBt, lBoundProj, lBoundTarg);
            std::printf(" RESULTS IN  jxMin, jxMax =%7d/2%7d/2\n", jxMin, jxMax);
            returnCode = 0;
        }

        if (jXInt == NOTDEF_INT) {
            // 342: JX not set
            if (jxMin == jxMax) {
                jx = (double)jxMin;
            } else {
                std::printf("\n***** ERROR: MORE THAN ONE JX IS POSSIBLE BUT JX WAS"
                            " NOT SET.  jxMin, jxMax =%5d/2%5d/2\n\n",
                            jxMin, jxMax);
                returnCode = 0;
            }
        } else if (!(jxMin <= jXInt && jXInt <= jxMax)) {
            // 345: validate user-specified JX is in [jxMin, jxMax]
            std::printf("\n**** ERROR, jxMin, jxMax, JX =%5d/2%5d/2%5d/2\n\n",
                        jxMin, jxMax, jXInt);
            returnCode = 0;
        }
    }

    }  // end non-CC block

    // =========================================================================
    //      COMPUTE L CRITICAL
    // =========================================================================
    // 400  lCrit = (lCrits(1)+lCrits(2))/2
    lCrit = (this->kin.lCrits[1] + this->kin.lCrits[2]) / 2;

    // =========================================================================
    //      IF AN ARRAY NAMED LVALUES IS DEFINED, WE USE IT
    // =========================================================================
    const std::vector<double>* lvalues = this->named.find("LVALUES  ");
    hasLvalues = (lvalues != nullptr) ? 1 : 0;
    if (lvalues != nullptr) {
        // lisArr: 1-based, indices [1..nLValues] valid.
        nLValues = (int)lvalues->size();
        this->inelastic.lisArr.assign(nLValues + 1, 0);

        // The reference implementation re-read the first L slot (the *same*
        // element) inside the i=1..n loop, so the entire LVALUES array was
        // effectively the first value. We preserve that behaviour bit-identically
        // by sourcing both lMin and the loop's L from lvalues->front().
        lMin = (int)(*lvalues)[0];

        lPrevious = -1;
        lStep = 1;

        for (i = 1; i <= nLValues; i++) {
            L = (int)(*lvalues)[0];
            if (!(L > lPrevious)) {
                std::printf("\n*** LVALUES ARRAY MUST BE INCREASING:  "
                            "ELEMENT%4d IS %5d, PREVIOUS IS%4d\n",
                            i, L, lPrevious);
                return;
            }
            this->inelastic.lisArr[i] = L;
            lStep = std::max(lStep, L - lPrevious);
            lPrevious = L;
        }

        lMax = L;
    } else {
        // =========================================================================
        //      NOW ESTIMATE lMin, lMax FROM L CRITICAL IF THEY WERE NOT GIVEN
        // =========================================================================
        lMinInput = lMin;
        if (lMin == NOTDEF_INT) {
            // NOTDEF_INT is therefore unreachable here.
            //      lMin = lCrit*alMnMt
            lMin = (int)(lCrit * alMnMt);
            //      lMin = std::max( 0, std::min( lMin, lCrit-lMinSub ) )
            lMin = std::max(0, std::min(lMin, lCrit - lMinSub));
        }

        if (lMax == NOTDEF_INT) {
            //      lMax = alMxMt*lCrit
            lMax = (int)(alMxMt * lCrit);
            //      lMax = std::max( lMax, lCrit + std::max(4*lStep, lMaxAdditional) )
            lMax = std::max(lMax, lCrit + std::max(4*lStep, lMaxAdditional));
        }
    }

    // =========================================================================
    //      IF THERE ARE IDENTICAL PARTICLES IN EXIT OR ENTRANCE CHANNEL,
    //      MUST FIX UP CHOICE OF L'S NOW.
    // =========================================================================
    //      lSkip = std::max( lSkips(1), lSkips(2) )
    densitySwitch = false;
    lSkip = std::max(this->distortedWave.channel[1].lSkips, this->distortedWave.channel[2].lSkips);
    bool skipLStepSetup = false;
    if (this->distortedWave.channel[1].statsCode != 3 && jAInt == 0) {
        // 470 branch via jAInt==0 fast path
        lParityBase = 0;
        densitySwitch = true;
        lStep = lStep + ((lStep) % (2));
    } else if (this->distortedWave.channel[2].statsCode == 3) {
        // 460/480 fall-through: skip L470 body
    } else if (jBInt > 0) {
        std::printf("\n**** ONLY SPIN-0 IDENTICAL PARTICLES MAY BE USED FOR NOW.\n");
        returnCode = 0;
        skipLStepSetup = true;
    } else {
        // 460 → fall through to L470
        lParityBase = ((lBoundProj + lBoundTarg) % (2));
        densitySwitch = true;
        lStep = lStep + ((lStep) % (2));
    }

    // =========================================================================
    //      MAKE lMin AND lMax lStep MULTIPLES
    // =========================================================================
    if (!skipLStepSetup) {
    if (lMinInput == NOTDEF_INT) lMin = lStep * (lMin / lStep);
    if (lSkip != 1) lMin = std::abs(lMin - lParityBase);
    lMax = lStep * ((lMax - lMin + lStep - 1) / lStep) + lMin;
    }

    // =========================================================================
    // Part 2: Process spectroscopic or deformation factors (statement 500+)
    // =========================================================================
    if (!isTransferReaction) {

    // Setup excitation potential parameters
        if (R0 != undefValue) R = r0Mass * R0;
        if (rI0 != undefValue) rI = r0Mass * rI0;
        if (this->opticalPotentialParams.rSi0 != undefValue) this->opticalPotentialParams.rSi = r0Mass * this->opticalPotentialParams.rSi0;
        taReal = this->opticalPotentialParams.A;
        taImag = this->opticalPotentialParams.aI;
        tvReal = this->opticalPotentialParams.V;
        tvImag = this->opticalPotentialParams.vI;
        if (rC0 != undefValue) rC = r0Mass * rC0;

        // Always use point and sphere Coulomb for effective excitation
        this->masses.rcProj = 0.0;
        this->masses.rcTarget = rC;

        if (rI == undefValue) rI = R;

        if (R == undefValue || rC == undefValue) {
            std::printf("\n**** BOTH R AND RC (OR R0 AND RC0) MUST BE"
                        " DEFINED TO COMPUTE DEFORMATION LENGTHS:\n"
                        "      R, RC = %8.3f %8.3f\n", R, rC);
            returnCode = 0;
            return;
        }

        R0 = R / r0Mass;
        rI0 = rI / r0Mass;
        rC0 = rC / r0Mass;
        coulombCoupling = -3.0 * zArray[1] * zArray[3] * Constants::hbar_c / Constants::fine_structure_inv;
        r2s[4] = coulombCoupling;

        {

        // Deformation lengths use only target or projectile radius
        r2Mass = std::pow(aMsPointer[2*iExcit - 2], 1.0/3.0);
        r2s[1] = R0 * r2Mass;
        r2s[2] = rI0 * r2Mass;
        r2s[3] = rC0 * r2Mass;

        // Walk the four BETA-family arrays (BETA / BETACOUL / BETARATS / BELX):
        //   - look up in reaction.named (filled by DEFINE for user-supplied
        //     values, or by alloc() below for the implicit broadcast case);
        //   - if missing OR a single scalar, alloc numLxI doubles and
        //     broadcast the scalar (default 1.0) across the lx grid;
        //   - if present with the wrong length, abort with the legacy error.
        //
        // The four resulting vector pointers are then offset by -lxMin/2 to
        // build the 1-based-style betaPointer/betaCPointer/betarPointer/belxPointer aliases the
        // loops below use as P_X[lx/2].
        std::vector<double>* betaVecs[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
        beta = 1.0;
        for (ii = 0; ii <= 3; ii++) {
            std::vector<double>* v = this->named.find(betaNames[ii+1]);
            hasBeta[ii+1] = (v != nullptr);
            len = 0;
            bool needAlloc = (v == nullptr);

            if (!needAlloc) {
                len = (int)v->size();
                if (len >= 1) beta = (*v)[0];
                if (len != numLxI) {
                    if (len <= 1) {
                        needAlloc = true;
                    } else {
                        std::printf("\n**** %.8s WAS DEFINED TO HAVE%4d"
                                    " VALUES BUT THERE ARE%4d VALUES OF LX\n",
                                    betaNames[ii+1], len, numLxI);
                        returnCode = 0;
                    }
                }
            }
            if (needAlloc) {
                v = &this->named.alloc(betaNames[ii+1], numLxI);
                for (i = 0; i < numLxI; i++) (*v)[i] = beta;
            }

            betaVecs[ii+1] = v;
            if (ii < 3) this->inelastic.poolBetas[ii] = v;
        }

        charge = (double)zArray[2*iExcit - 1];
        this->inelastic.betnrArr.assign(numLxI, 0.0);
        { int off = -lxMin/2;
        betaPointer  = betaVecs[1]->data() + off;   // BETA
        betaCPointer = betaVecs[2]->data() + off;   // BETACOUL
        betarPointer = betaVecs[3]->data() + off;   // BETARATS
        belxPointer  = betaVecs[4]->data() + off;   // BELX
        betnrPointer = this->inelastic.betnrArr.data() + off;
        }

        mProj = jIn % 2;

        // Make sure BETA(Coulomb) are defined
        if (!hasBeta[2]) {
            if (hasBeta[4]) {
                for (int lx = lxMin; lx <= lxMax; lx += 2) {
                    betaCPointer[lx/2] =
                        (4.0*Constants::PI / (3.0*charge))
                        * std::sqrt(belxPointer[lx/2])
                        / (std::pow(r2s[3]/10.0, (double)lx)
                           * clebschGordan(jIn, 2*lx, mProj, 0, jOut, mProj));
                }
            } else {
                std::printf("\n**** WARNING:  BETA COULOMB IS BEING CHOOSEN SUCH THAT:\n"
                            "                RC*BETACOULOMB = R*BETA"
                            "    WHERE RC AND R REFER ONLY TO THE EXCITED NUCLEUS\n");
                for (int lx = lxMin; lx <= lxMax; lx += 2)
                    betaCPointer[lx/2] = (r2s[1] / r2s[3]) * betaPointer[lx/2];
            }
        }

        if (!hasBeta[1]) {
            if (hasBeta[2] || hasBeta[4]) {
                for (int lx = lxMin; lx <= lxMax; lx += 2)
                    betaPointer[lx/2] = (r2s[3] / r2s[1]) * betaCPointer[lx/2];
            } else {
                std::printf("\n");
                printRepeatedChar('*', 120);
                std::printf("\n **");
                printRepeatedChar(' ', 115);
                std::printf("**\n");
                std::printf(" **    WARNING:  BETA ASSUMED TO BE 1 --"
                            " MUST MULTIPLY CROSS SECTIONS BY SQUARE OF ACTUAL"
                            " BETA TO GET TRUE VALUES.");
                printRepeatedChar(' ', 19);
                std::printf("**\n ** ");
                printRepeatedChar(' ', 115);
                std::printf("**\n ");
                printRepeatedChar('*', 120);
                std::printf("\n");
            }
        }

        // Define auxiliary quantities
        for (int lx = lxMin; lx <= lxMax; lx += 2) {
            temp = clebschGordan(jIn, 2*lx, mProj, 0, jOut, mProj);
            betarPointer[lx/2] = betaCPointer[lx/2] * temp
                           * std::pow(r2s[3], (double)lx) / (2*lx + 1.0);
            betnrPointer[lx/2] = betaPointer[lx/2] * temp;
            belxPointer[lx/2]  = std::pow(
                3.0 * charge * betaCPointer[lx/2] * temp
                * std::pow(r2s[3]/10.0, (double)lx) / (4.0*Constants::PI),
                2.0);
        }

        // poolBetas[1] sentinel reset to nullptr.
        this->inelastic.poolBetas[1] = nullptr;
        }  // end beta setup block
    } else {

    // =========================================================================
    // Process spectroscopic info for transfer
    // =========================================================================
    if (this->spec.specAmpProj == 1.0) this->spec.specAmpProj = std::sqrt(this->spec.specFactorProj);
    this->spec.specFactorProj = this->spec.specAmpProj * this->spec.specAmpProj;
    if (this->spec.specAmpTgt == 1.0) this->spec.specAmpTgt = std::sqrt(this->spec.specFactorTgt);
    this->spec.specFactorTgt = this->spec.specAmpTgt * this->spec.specAmpTgt;

    this->inelastic.atermArr.assign(lxMax + 1, 0.0);

    temp = (jBigB + 1.0) / (jBigA + 1.0);
    if (stripPickup == -1) temp = (jb + 1.0) / (ja + 1.0);
    temp = std::sqrt(temp);

    for (int lx = lxMin; lx <= lxMax; lx++) {
        double racahCoefficient = racah(2*lBoundTarg, jBt, 2*lBoundProj, jBp, jXInt, 2*lx);
        aTerm = temp * std::sqrt(2*lx + 1.0) * this->spec.specAmpProj * this->spec.specAmpTgt
              * racahCoefficient;
        phaseParity = jXInt - jBp + 2*(lBoundProj + lBoundTarg);
        if (stripPickup == -1)
            phaseParity = jXInt - jBt + jAInt + jBigA_INT - jBInt - jBigB_INT;
        phaseParity = phaseParity/2 + 1;
        if (phaseParity % 2 != 0) aTerm = -aTerm;
        this->inelastic.atermArr[lx] = aTerm;
        if (verbosity >= 4)
            std::printf(" LX =%3d     SIGN, aTerm =%4d%15.8G\n",
                        lx, phaseParity, aTerm);
    }

    }  // end if(isTransferReaction)/else — spec/inelastic dispatch

    // =========================================================================
    // Produce a nice summary of the reaction
    // =========================================================================
    if (printSwitch) {

    {
        numParticles = 4 + std::abs(stripPickup);

        std::printf("0%20sSUMMARY OF THE REACTION\n", "");
        std::printf("0%23s", "");
        for (int ch = 1; ch <= 45; ch++) std::printf("%c", this->reactStr[ch]);
        std::printf("\n");
        std::printf("0%17sBIGA   (   A   ,    B   )  BIGB   %s\n",
                    "", xWords[numParticles - 3]);
        std::printf("0M (AMU)      ");
        for (int j = 1; j <= numParticles; j++)
            std::printf("%9.2f", aMsPointer[ind[j] - 1]);
        std::printf("\n");

        std::printf(" Z         ");
        for (int j = 1; j <= numParticles; j++)
            std::printf("%9d", zArray[ind[j]]);
        std::printf("\n");

        std::printf(" E* (MEV)       ");
        for (int j = 1; j <= numParticles; j++)
            std::printf("%9.4f", this->energies.exs[ind[j]]);
        std::printf("\n");

        std::printf(" J           ");
        for (int j = 1; j <= numParticles; j++)
            std::printf("%7d/2", (int)this->angMom.js[ind[j]]);
        std::printf("\n");

        std::printf("PARITY        ");
        for (int j = 1; j <= numParticles; j++) {
            std::printf("%7.7s", parityWord[parities[ind[j]] + 2]);
            if (j < numParticles) std::printf("  ");
        }
        std::printf("\n");
    }

    if (isTransferReaction) {
        // Bound state properties (transfer only)
        std::printf("\n0     BOUND STATE PROPERTIES\n"
                    "0           PROJECTILE   TARGET\n");
        std::printf(" E         %10.4f%10.4f          Q =%9.4f\n",
                    this->internalState.eBnds[1], this->internalState.eBnds[2], this->energies.Q);
        std::printf(" JP      %8d/2%8d/2\n", jBp, jBt);
        std::printf(" L      %10d%10d\n", lBoundProj, lBoundTarg);
        std::printf(" NODES  %10d%10d\n", this->boundState.vertex[1].nodeCount, this->boundState.vertex[2].nodeCount);
        std::printf(" SPEC. AMP.%10.4f%10.4f\n",
                    this->spec.specAmpProj, this->spec.specAmpTgt);
        std::printf(" SPEC. factor%8.4f%10.4f\n",
                    this->spec.specFactorProj, this->spec.specFactorTgt);

        if (stripPickup ==  1) std::printf("0THIS IS A STRIPPING REACTION\n");
        if (stripPickup == -1) std::printf("0THIS IS A PICKUP REACTION\n");

        std::printf("0THE EXCHANGE POTENTIAL CONSISTS OF\n");
        // IVRTEX permanently 1: always print PROJECTILE VERTEX.
        std::printf("    THE NUCLEAR POTENTIAL AT THE PROJECTILE VERTEX\n");
        // nuConL∈{2,3} (defaults seeds 3, the linkule-clash branch above
        std::printf("    AND THE COULOMB POTENTIAL INCLUDING CORE-CORE"
                    " CORRECTION TERMS\n");
        if (nuConL == 3)
            std::printf("    AND THE CORE-CORE CORRECTIONS FROM THE REAL"
                        " PART OF THE NUCLEAR OPTICAL POTENTIAL\n");
    } else {
        // Statement 680: Printout for inelastic scattering
        std::printf("\nTHE REACTION IS INELASTIC EXCITATION OF THE %s%s\n",
                    excitationWords[2*iExcit-1], excitationWords[2*iExcit]);
        std::printf("\nLX         DEFORMATION PARAMETERS%21sB(E(LX))\n"
                    "            NUCLEAR         COULOMB        (E**2 BARN**LX)\n", "");
        for (int lx = lxMin; lx <= lxMax; lx += 2) {
            std::printf("%3d     %15.5G%15.5G%15.5G\n",
                        lx,
                        betaPointer[lx/2],
                        betaCPointer[lx/2],
                        belxPointer[lx/2]);
        }
        std::printf("\nDEFORMATION RADII:\n"
                    "\n  REAL%14.4f\n"
                    "   IMAGINARY%9.4f\n"
                    "   COULOMB%11.4f\n", r2s[1], r2s[2], r2s[3]);
    }

    // =========================================================================
    // Statement 750: Critical-L summary and L-range printout
    // =========================================================================
    std::printf("\n0ESTIMATED CRITICAL L'S:    INCOMING =%4d"
                "     OUTGOING =%4d     AVERAGE =%4d\n",
                this->kin.lCrits[1], this->kin.lCrits[2], lCrit);

    std::printf("0%5d  =< LX =<%3d\n"
                "+%23s%5d  =<  lIn, lOut  =<%4d\n",
                lxMin, lxMax, "", lMin, lMax);

    if (hasLvalues == 0)
        std::printf("+%58sLSTEP =%3d\n", "", lStep);
    if (hasLvalues != 0)
        std::printf("+%58sUSING SPECIFIED VALUES WITH A MAXIMUM LSTEP OF%4d\n",
                    "", lStep);

    }  // end if (printSwitch) — summary printout

    // =========================================================================
    // Q-value consistency warnings
    // =========================================================================
    if (isTransferReaction &&
        std::abs(this->internalState.eBnds[1] - this->internalState.eBnds[2] - (double)stripPickup * this->energies.Q) > 0.001)
        std::printf("\n**** WARNING:  Q IS NOT EQUAL TO THE DIFFERENCE"
                    " OF THE BOUND STATE ENERGIES.\n");

    if (std::abs(this->energies.eCm + this->energies.Q - this->distortedWave.channel[2].Ecm) > 0.001)
        std::printf("\n**** WARNING:  THE DIFFERENCE OF THE C.M. SCATTERING"
                    " ENERGIES IS NOT EQUAL TO Q.\n");

    this->kin.lOutMax = lMax + lxMax;
    this->internalState.lInMax = lMax;

    if (lCrit == 0) lCrit = lMax / 2;

    return;

    #undef jAInt
    #undef jBInt
    #undef jBigA_INT
    #undef jBigB_INT
    #undef jXInt
} // end Reaction::probePrint

// ---------------------------------------------------------------------------
// Backward-compat free function wrapper — calls the class method
// ---------------------------------------------------------------------------
