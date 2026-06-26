// channel_setup.cpp — CHANEL: sets up per-channel parameters.

#include "ptolemy_types.h"
#include "masstable.h"
#include <cstdio>
#include "Reaction.h"
#include "InputParser.h"

// File-scope alias to the InputParser-owned scanner buffer. Was an
namespace { auto& inputBuffer = InputParser::buffer(); }

// ============================================================================
//
//     SCANS A CHANNEL SPECIFICATION
//
//     FORMS ARE
//        N15 + 209BI(9/2- 1.3)     ( FOR SCATTERING )
//        P + 208PB = BI209(7/2 .9)    ( FOR BOUND STATES )
//     NOTE THAT PROJECTILE ALWAYS COMES FIRST.
//
//
//     AND STORES  Z'S,  M'S,  J'S,  E*'S   OF
//     THE 2 OR 3 PARTICLES INTO THE REACTION STRUCTURES.
//     FOR BOUNDSTATES THE TOTAL J AND BINDING ENERGY ARE ALSO STORED.
//
//     9/22/75 - BASED ON REACTN - S. PIEPER
//     4/24/77 - STORE E* IN exsPt
//     9/24/77 - RECOGNIZE E*=0; ALLOW LATTER DEF OF MASS EXCESS
//     12/10/77 - FIX ERROR IN ABOVE TO DEFINE J FOR B.S.
//     5/6/79 - G.S. MASS EXCESS, PARITIES, L & nNodes FROM LEVEL - S.P.
//     12/31/79 - STANDARD MORTRAN - RPG
//     1/25/80 - CHARACTER SYMBOL, ALPHA - RPG
//
// ============================================================================

bool Reaction::setupChannel()
{
    Reaction& reaction = *this;
    double* amPts = &reaction.masses.massProj;   // amPts[1]=massProj, amPts[2]=massTgt (1-based)

    // izPts(2) <=> zProj, zTarget  (contiguous ints in INTGER)
    int* izPts = &reaction.charges.zProj;            // izPts[1]=zProj, izPts[2]=zTarget (1-based)

    // spinProj and spinTarget are contiguous doubles; jSpts is a 1-based
    // pointer view over them.
    double* jSpts = &reaction.angMom.spinProj;         // jSpts[1]=spinProj, jSpts[2]=spinTarget (1-based)

    // Local variables
    char8  guy[4];         // 1-based: guy[1]..guy[3]
    double eStars[4];      // 1-based: eStars[1]..eStars[3]
    double excesses[4];      // 1-based: excesses[1]..excesses[3]
    int    nodeVals[4];      // 1-based
    int    lVals[4];       // 1-based
    int    jVals[4];       // 1-based
    int    iParities[4];      // 1-based

    // 1-based: names[1]..names[3] (sole reader is the %-10s warning printf
    // for ground-state spin mismatch); [0] unused.
    static const char* names[4] = { "", "PROJECTILE", "TARGET", "BOUND STATE" };

    int    iz = 0, ia = 0, returnCode, numSymbols, i, N;
    int    inChStart;
    int    aBound = 0, zBound = 0;
    double atomicMassExcess;
    int    jj;
    double eStar;
    int    jVal, gsParity;

    //
    //     GET PAST ANY INITIAL DELIMITORS
    //
    while (true) {
        char myChar = inputBuffer.iBuf[inputBuffer.inCh];
        if ((myChar >= 'A' && myChar <= 'Z') ||
            (myChar >= 'a' && myChar <= 'z') ||
            (myChar >= '0' && myChar <= '9')) break;
        inputBuffer.inCh = inputBuffer.inCh + 1;
        if (inputBuffer.inCh > inputBuffer.nOch) {
            std::printf("\n **** CHANNEL KEYWORD AND COMPLETE SPECIFICATION"
                        " MUST BE ON ONE LINE.\n");
            return false;
        }
    }
    //
    //     START OF SPECIFICATION FOUND.  BREAK IT INTO 2 OR 3 PIECES.
    //
    inChStart = inputBuffer.inCh;
    channelScan(guy, eStars, nodeVals, lVals, jVals, iParities, returnCode, reaction);
    if (returnCode < 0)  return false;
    //
    //     NOW GET THE A AND Z VALUES
    //
    numSymbols = returnCode;
    for (i = 1; i <= numSymbols; i++) {
        azCode(guy[i].data, iz, ia, returnCode);
        if (returnCode != 0) {
            if (returnCode != -2)  std::printf("\n **** A SYMBOL HAS INCORRECT SYNTAX:  "
                                         "SYMBOL = %.8s\n", guy[i].data);
            if (returnCode == -2)  std::printf("\n **** THE SYMBOL %.6s IS NOT A KNOWN ELEMENT.\n",
                                         guy[i].data);
            return false;
        }
        // returnCode == 0 (label L230):
        if (i == 3) {
            // bound state
            reaction.angMom.J = jVals[3];
            reaction.angMom.parity = iParities[3];
            reaction.angMom.nNodes = nodeVals[3];
            reaction.angMom.L = lVals[3];
            aBound = ia;
            zBound = iz;
            if (!(aBound == (int)(amPts[0] + amPts[1]) && iz == izPts[0] + izPts[1])) {
                std::printf("\n **** BOUNDSTATE DOES NOT CONSERVE NUCLEON NUMBER"
                            " OR CHARGE:\n"
                            " ZP, ZT, Z(BOUND) =%8d%8d%8d\n"
                            " MP, MT, M(BOUND) =%8.0f%8.0f%8d\n",
                            izPts[0], izPts[1], iz,
                            amPts[0], amPts[1], ia);
                return false;
            }
        } else {
            // i == 1 or 2
            izPts[i-1] = iz;
            amPts[i-1] = ia;
            // Only i==2 slot is read downstream (set_channels). i==1 write dropped.
            if (i == 2) {
                reaction.internalState.nodePt2 = nodeVals[i];
                reaction.internalState.lSpcPt2 = lVals[i];
            }
            jSpts[i-1] = jVals[i];
            reaction.energies.exsPt[i] = eStars[i];
            if (reaction.energies.exsPt[i] == reaction.internalState.undefValue)  reaction.energies.exsPt[i] = 0;
            reaction.angMom.parityPt[i] = iParities[i];
        }
    }

    //
    //     NOW GET GROUND STATE MASS EXCESS AND DETERMIN IF THE NUCLEII
    //     EXIST.  (L300 block)
    //
    for (i = 1; i <= numSymbols; i++) {
        if (i == 3) {
            ia = aBound;
            iz = zBound;
        } else {
            ia = (int)amPts[i-1];
            iz = izPts[i-1];
        }
        N = ia - iz;
        atomicMassExcess = excess(iz, ia, returnCode);
        if (returnCode != 0) {
            std::printf("\n **** THE NUCLEUS WITH A =%4d,  Z =%4d"
                        ",  N =%4d  IS NOT BOUND "
                        "ACCORDING TO PTOLEMY\"S MASS TABLE.\n", ia, iz, N);
            // fall through to L359 (continue loop)
        } else {
            excesses[i] = atomicMassExcess;
            if (i <= 2) {
                reaction.masses.amxgPt[i] = atomicMassExcess;
                excesses[i] = excesses[i] + reaction.energies.exsPt[i];
            }
            if (i == 3 && eStars[3] != reaction.internalState.undefValue)
                excesses[3] = excesses[3] + eStars[3];
        }
        // continue
    }
    //
    //     NOW GET J FOR EACH OF THE GROUND STATES
    //
    for (i = 1; i <= numSymbols; i++) {
        //
        if (i == 3) {
            ia = aBound;
            iz = zBound;
            jj = (int)reaction.angMom.J;
        } else {
            ia = (int)amPts[i-1];
            iz = izPts[i-1];
            jj = (int)jSpts[i-1];
        }
        eStar = eStars[i];
        //
        //     IF E* HAS BEEN ENTERED THEN G.S. J IS MEANINGLESS
        //
        if (eStar > 1.0e-5 && eStar != reaction.internalState.undefValue) continue;  // goto L459

        N = ia - iz;
        groundStateInfo(iz, ia, jVal, gsParity, returnCode);

        if (returnCode == 0 && jVal != NOTDEF_INT) {
            // J is known from mass table
            if (jj != NOTDEF_INT) {
                // jj already defined
                if (jj == jVal || eStar != reaction.internalState.undefValue) {
                    // goto L459 (continue)
                } else {
                    std::printf("\n **** WARNING:  GROUND STATE SPIN FOR %-10s"
                                " IS%3d/2.  YOU HAVE SPECIFIED J"
                                " =%3d/2, BUT HAVE NOT SPECIFIED E*\n",
                                names[i], jVal, jj);
                }
                continue;  // goto L459
            }
            // jj == NOTDEF_INT: fall to L450
        } else {
            // returnCode != 0 or jVal == NOTDEF_INT:
            //     No error if J is already defined or explicitly given.
            if (jj != NOTDEF_INT) continue;  // goto L459

            //     DINEUTRON OR DIPROTON IS J = 0
            if (ia == 2 && iz != 1) {
                jVal = 0;
                // fall to L450
            } else {
                std::printf("\n **** WARNING:  THE GROUND STATE SPIN OF THE"
                            " NUCLEUS WITH A =%4d,  Z =%4d,  N =%4d"
                            "  IS NOT KNOWN TO PTOLEMY.\n", ia, iz, N);
                continue;  // goto L459
            }
        }

        // store J and parity
        if (i == 3) {
            reaction.angMom.J = jVal;
            if (gsParity != 0)  reaction.angMom.parity = gsParity;
        } else {
            jSpts[i-1] = jVal;
            if (gsParity != 0)  reaction.angMom.parityPt[i] = gsParity;
        }
        // continue
    }
    //
    //     DEFINE BINDING ENERGY
    //
    if (numSymbols == 3)
        reaction.energies.E = excesses[3] - excesses[1] - excesses[2];
    //
    //     IT WAS SUCCESFUL, COPY REACTION FOR header
    //
    for (i = 1; i <= 45; i++) {
        reaction.reactStr[i] = ' ';
    }
    N = std::min(45, inputBuffer.inCh - inChStart);
    for (i = 1; i <= N; i++) {
        reaction.reactStr[i] = inputBuffer.iBuf[inChStart - 1 + i];
    }
    return true;
}

// ============================================================================
// channelScan — Scans an input line for channel specification
// (folded in from source_channels.cpp)
void channelScan(char8* guy, double* eStars, int* nodeVals, int* lVals, int* jVals,
            int* iParities, int& returnCode, Reaction& reaction)
{
    // Character constants — only plus has 2+ uses, the rest were single-use
    const char plus  = '+';

    int i, chIndex;
    int inChStart;
    bool seenPlus, seenEqual;
    int tokenCode, intValue, messageLength;
    char stopChar;
    char8 charValue;
    double value, dummyValue;
    char dummy[4] = {};   // unused message parameter for mScan

    returnCode    = 0;
    inChStart  = inputBuffer.inCh;
    seenPlus  = false;
    seenEqual  = false;

    // Shared error epilogue: print "WAS PROCESSING FIELD i ..."
    // and the captured inChStart..inCh range, then set returnCode=-1.
    auto reportErrorAndReturn = [&]() {
        int inChM1 = inputBuffer.inCh - 1;
        std::printf("      WAS PROCESSING FIELD%2d OF THE SPECIFICATION\n"
                    "      SPECIFICATION UP TO THE ERROR IS: ", i);
        for (chIndex = inChStart; chIndex <= inChM1 && chIndex <= iBufSize; chIndex++)
            std::printf("%c", inputBuffer.iBuf[chIndex]);
        std::printf("\n");
        returnCode = -1;
    };

    // tokenCode == -8 means the channel spec ran off the end of the input line.
    // Same message + error epilogue appeared at two scan points (symbol scan and
    // the /2 re-scan); shared here. Caller still issues the `return`.
    auto reportIncompleteSpec = [&]() {
        std::printf("0**** COMPLETE SPECIFICATION MUST BE ON SAME LINE AS"
                    " CHANNEL KEYWORD.\n");
        reportErrorAndReturn();
    };

    // The "CHANNEL SPECIFICATION REQUIRES TWO NUCLEAR SYMBOLS" message + error
    // epilogue appears at two end-of-scan points (mid-loop blank-skip exhaustion
    // and the final delimiter check); shared here. Caller still issues any return.
    auto reportTwoSymbolsRequired = [&]() {
        std::printf("0**** CHANNEL SPECIFICATION REQUIRES TWO"
                    " NUCLEAR SYMBOLS SEPARATED BY A \"+\".\n");
        reportErrorAndReturn();
    };

    // The mScan(18) field read plus its -9/-8 bail guards appear at two scan
    // points in the excited-state loop (loop top and the /2 bypass re-scan);
    // shared here. Returns true if the caller should bail (caller issues return).
    auto scanExcitedFieldOrBail = [&]() -> bool {
        mScan(18, tokenCode, charValue, value, intValue, dummy, messageLength, stopChar);
        if (tokenCode == -9) { reportErrorAndReturn(); return true; }
        if (tokenCode == -8) { reportIncompleteSpec(); return true; }
        return false;
    };

    // Initialise all 3 slots (1-based)
    char8 blank8;   // default ctor fills .data with 8 spaces
    for (i = 1; i <= 3; i++) {
        guy[i]    = blank8;
        nodeVals[i] = NOTDEF_INT;
        lVals[i]  = NOTDEF_INT;
        iParities[i] = 0;
        jVals[i]  = NOTDEF_INT;
        eStars[i] = reaction.internalState.undefValue;
    }

    i = 1;

// -----------------------------------------------------------------------
//  LOOP ON REACTION PARTICIPANTS
// -----------------------------------------------------------------------
    while (true) {
    mScan(-20, tokenCode, guy[i], dummyValue, intValue, dummy, messageLength, stopChar);
    if (messageLength > 5) {
        std::printf("0**** TOO MANY CHARACTERS IN SYMBOL AND ATOMIC MASS.\n");
        { reportErrorAndReturn(); return; }
    }
    if (messageLength == 0) {
        std::printf("0**** UNEXPECTED DELIMITOR OR END OF LINE ENCOUNTERED"
                    " WHEN LOOKING FOR ELEMENT SYMBOL.\n");
        { reportErrorAndReturn(); return; }
    }
    inputBuffer.inCh = inputBuffer.inCh + 1;   // advance past stop char

    // If stopped by '(' → process excited-state / J specification
    if (stopChar == '(') {

// -----------------------------------------------------------------------
//  PROCESS EXCITED STATE INFO  (E*, J, parity, level descriptor)
// -----------------------------------------------------------------------
    while (true) {
    if (scanExcitedFieldOrBail()) return;
    if (tokenCode == -2) {
        // It is E*
        if (eStars[i] != reaction.internalState.undefValue) {
            std::printf("0**** ATTEMPTING TO ENTER 2 VALUES FOR E*\n");
            { reportErrorAndReturn(); return; }
        }
        eStars[i] = value;
        continue;
    }

    if (tokenCode == -5) {
        // It is J (half-integer token)
        if (jVals[i] != NOTDEF_INT) {
            std::printf("0**** ATTEMPTING TO ENTER 2 VALUES FOR J\n");
            { reportErrorAndReturn(); return; }
        }
        jVals[i] = 2 * intValue;
        if (stopChar == '/') {
            // Bypass /2: advance past '/', then verify the '2'
            inputBuffer.inCh = inputBuffer.inCh + 1;
            jVals[i] = intValue;
            if (scanExcitedFieldOrBail()) return;
            if (!(tokenCode == -5 && intValue == 2)) {
                std::printf("0**** INVALID J FIELD\n");
                { reportErrorAndReturn(); return; }
            }
        }

        // Parity may be given by terminal + or -
        if (stopChar == plus)  { iParities[i] = +1; inputBuffer.inCh++; }
        else if (stopChar == '-') { iParities[i] = -1; inputBuffer.inCh++; }
        continue;
    }

    if (tokenCode == -3) {
        // It is a nodes-and-L level descriptor (e.g. 2P3/2-)
        if (lVals[i] != NOTDEF_INT) {
            std::printf("0**** ATTEMPTING TO ENTER TWO LEVEL DESCRIPTORS\n");
            { reportErrorAndReturn(); return; }
        }
        lVals[i]  = intValue;
        nodeVals[i] = (int)value;
        continue;
    }

    // Is it the closing ')'?
    if (tokenCode != 11) {
        std::printf("0**** UNRECOGNIZABLE EXCITED STATE DESCRIPTOR (TOKEN =%4d)\n", tokenCode);
        { reportErrorAndReturn(); return; }
    }
    // Done: bypass the ')' already consumed by mScan; set stopChar to char after it
    stopChar         = inputBuffer.iBuf[inputBuffer.inCh];
    inputBuffer.inCh    = inputBuffer.inCh + 1;
    break;
    }  // end while (excited-state scan)
    }

// -----------------------------------------------------------------------
//  READY TO PROCESS THE hasNextBlock FIELD (or end of channel spec)
// -----------------------------------------------------------------------
    inputBuffer.inCh = inputBuffer.inCh - 1;
    if (seenPlus && seenEqual) { returnCode = 3; return; }  // seenEqual is guaranteed true inside the && branch.

    // Skip blanks
    while (true) {
        inputBuffer.inCh = inputBuffer.inCh + 1;
        if (stopChar != ' ') break;
        if (inputBuffer.inCh >= inputBuffer.nOch) {
            if (seenPlus) { returnCode = seenEqual ? 3 : 2; return; }
            { reportTwoSymbolsRequired(); return; }
        }
        stopChar = inputBuffer.iBuf[inputBuffer.inCh];
    }

    // Only '+' and '=' are valid delimiters inside the specification
    if (stopChar == '=') {
        // '=' SIGN — bound-state result
        if (seenEqual) {
            std::printf("0**** TWO \"=\" SIGNS ARE NOT ALLOWED\n");
            { reportErrorAndReturn(); return; }
        }
        seenEqual = true;
        // Always store result of '=' in 3rd slot
        if (i == 2) {
            i = 3;
        } else {
            guy[3]    = guy[1];
            jVals[3]  = jVals[1];
            eStars[3] = eStars[1];
            // i stays the same (1), re-scan guy[1]
        }
        continue;
    }
    if (stopChar == plus) {
        // '+' SIGN
        if (seenPlus) {
            std::printf("0**** TWO \"+\" SIGNS ARE NOT ALLOWED.\n");
            { reportErrorAndReturn(); return; }
        }
        seenPlus = true;
        i = i + 1;
        continue;
    }
    break;
    }  // end while (reaction-participants)

    // End of channel spec: stopChar will be reconsidered by CONTRL
    inputBuffer.inCh = inputBuffer.inCh - 1;
    if (seenPlus) { returnCode = 2; if (seenEqual) returnCode = 3; return; }
    reportTwoSymbolsRequired();
}
