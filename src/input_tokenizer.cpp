// input_tokenizer.cpp — input scanner / tokenizer: card buffer handling and the
// token scanners (mScan, nxWord, newCard, qvScan, InputParser::defineArray).

#include "ptolemy_types.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include "Reaction.h"
#include "InputParser.h"

// File-scope alias to the InputParser-owned scanner buffer.
namespace { auto& inputBuffer = InputParser::buffer(); }

// ============================================================================
// mScan — Field scanner for input
// ============================================================================
void mScan(int keyEx, int& tokenCode, char8& charValue, double& value, int& intValue,
           char* message, int& messageLength, char& stop)
{
    // Local character constants — only the multi-use letter/digit range
    // boundaries survive; DOLLAR/DOT/STAR were single-use and inlined
    const char blank = ' ';
    const char num0 = '0';
    const char num9 = '9';
    static const char lChars[10] = { ' ', 'S', 'P', 'D', 'F', 'G', 'H', 'i', 'J', 'K' };
    int numLs = 9;

    char8 word, blank8(" ");
    char word1[9]; // 1-based
    // wordLength is normally reset to 0 in processAlpha's inPhase==1 branch
    // before any read; init to 0 so the (theoretical) "alpha_done called
    // before any alpha char processed" path returns messageLength=0 instead of UB.
    // Silences the pre-existing -Wmaybe-uninitialized warning.
    int wordLength = 0;

    static const char chars[22] = { ' ', '%', '#', '&', ':', ';',
        '-', '+', '/', '*', '(',
        ')', '=', '@', '<', '>',
        '|', '~', '?', '"', '.',
        ',' };

    static const char delimiters[5] = { ' ', '+', '-', '\'', '$' };

    char curChar;
    double val, power, sign, power10, value1;
    int digitCount, decimalSeen, inPhase, inComment, expSign, expValue;
    int inChStart, i;

    // Presets
    inChStart = inputBuffer.inCh;
    sign = +1.0;
    val = 0.0;
    digitCount = 0;
    inputBuffer.inCh = inputBuffer.inCh - 1;
    messageLength = 0;
    inComment = 0;
    decimalSeen = 0;

    // Phases: 1=nothing, 2=alpha, 3=numeric
    inPhase = 1;

    // Shared error epilogue: print "FIELD UP TO THE ERROR IS: ..."
    // then set tokenCode = -9.
    auto reportFieldErrorAndReturn = [&]() {
        std::printf(" FIELD UP TO THE ERROR IS: ");
        for (int ii = inChStart; ii <= inputBuffer.inCh && ii <= iBufSize; ii++)
            std::printf("%c", inputBuffer.iBuf[ii]);
        std::printf("\n");
        tokenCode = -9;
    };

    // "INVALID DECIMAL POINT" message + field-error epilogue recurs at two
    // period-handling points (inPhase==2, and inPhase==3 with a decimal already
    // seen); shared here. Caller still issues the return.
    auto reportInvalidDecimal = [&]() {
        std::printf("\n**** INVALID DECIMAL POINT.\n");
        reportFieldErrorAndReturn();
    };

    // "INVALID ALPHABETIC IN NUMERIC FIELD" message + field-error epilogue
    // recurs at two points in the alpha-character handler (level-char lookup
    // miss, and a non-E/non-D alpha); shared here. Caller still issues return.
    auto reportInvalidAlpha = [&]() {
        std::printf("\n**** INVALID ALPHABETIC IN NUMERIC FIELD\n");
        reportFieldErrorAndReturn();
    };

    // Shared field-found epilogue: record the next-char stop.
    auto fieldFound = [&]() {
        stop = (inputBuffer.inCh >= 1 && inputBuffer.inCh <= iBufSize)
                   ? inputBuffer.iBuf[inputBuffer.inCh]
                   : ' ';
    };

    // Shared alpha-done epilogue: emit collected word as charValue.
    auto alphaDoneAndExit = [&]() {
        messageLength = wordLength;
        std::memset(word.data, ' ', 8);
        for (int k = 1; k <= wordLength && k <= 8; k++) word.data[k-1] = word1[k];
        charValue = word;
        tokenCode = -1;
        fieldFound();
    };

    // Shared numeric-done epilogue: compute value/intValue and exit.
    auto numericDoneAndExit = [&]() {
        tokenCode = -2;
        value1 = val / power;
        value = value1 * sign;
        intValue = (int)value;
        if (decimalSeen == 0) tokenCode = -5;
        if (digitCount == 0) {
            std::printf("\n**** MISPLACED DECIMAL POINT OR PLUS OR MINUS SIGN.\n");
            reportFieldErrorAndReturn();
            return;
        }
        fieldFound();
    };

    // Alpha character in numeric field: either a level descriptor
    // (when decimalSeen == 0) or the E/D of E-format exponent. Always exits mScan
    // (via fieldFound, numericDoneAndExit, or reportFieldErrorAndReturn).
    auto alphaInNumericField = [&]() {
        if (curChar != 'E') {
            if (decimalSeen == 0) {
                // Alphabetic found in integer numeric field - level specification
                value = sign * val - 1.0;
                for (i = 1; i <= numLs; i++) {
                    if (curChar == lChars[i]) break;
                }
                if (i > numLs) { reportInvalidAlpha(); return; }
                intValue = i - 1;
                tokenCode = -3;
                inputBuffer.inCh = inputBuffer.inCh + 1;
                fieldFound();
                return;
            }
            if (curChar != 'D') { reportInvalidAlpha(); return; }
        }

        // E or D format
        expSign = 1;
        decimalSeen = 1;
        expValue = 0;
        inputBuffer.inCh = inputBuffer.inCh + 1;
        if (inputBuffer.inCh > inputBuffer.nOch) { std::printf("\n**** E OR D MUST BE FOLLOWED BY NUMBER\n"); reportFieldErrorAndReturn(); return; }
        curChar = inputBuffer.iBuf[inputBuffer.inCh];
        {
        bool needAdvance = true;  // L334 entry advances first; L335 entry skips
        if (curChar == delimiters[1]) {
            // explicit '+' — consume it; L334 will advance
        } else if (curChar == delimiters[2]) {
            expSign = -1;  // consume '-'; L334 will advance
        } else {
            needAdvance = false;  // no sign — check current curChar as digit (L335 entry)
        }
        while (true) {
            if (needAdvance) {
                inputBuffer.inCh = inputBuffer.inCh + 1;
                if (inputBuffer.inCh > inputBuffer.nOch) break;  // was goto L336
                curChar = inputBuffer.iBuf[inputBuffer.inCh];
            }
            needAdvance = true;
            if (curChar < num0 || curChar > num9) break;  // was goto L336
            expValue = expValue * 10 + (curChar - num0);
            if (expValue > 300) { std::printf("\n**** INVALID OR TOO LARGE EXPONENT FIELD.\n"); reportFieldErrorAndReturn(); return; }
        }
        }
        if (expValue != 0) {
            power10 = std::pow(10.0, (double)expValue);
            if (expSign < 0) {
                power = power * power10;
            } else {
                power = power / power10;
            }
        }
        numericDoneAndExit();
    };

    // Process current alphabetic char into the active word.
    // Returns true if caller should `continue` the scanner loop, false if
    // it should `return` (alphaInNumericField exited mScan).
    auto processAlpha = [&]() -> bool {
        if (inPhase == 3) { alphaInNumericField(); return false; }
        if (inPhase == 1) {
            word = blank8;
            std::memcpy(word1 + 1, word.data, 8);
            wordLength = 0;
        }
        // append char to current word
        if (wordLength < 8) { word1[wordLength + 1] = curChar; wordLength = wordLength + 1; }
        inPhase = 2;
        return true;
    };

    // Shared special-character dispatch: either skip the char
    // (returns true → caller continues scanner) or emit it as a found field
    // (returns false → caller returns from mScan).
    auto handleSpecialChar = [&]() -> bool {
        if (tokenCode > std::abs(keyEx)) return true;
        inputBuffer.inCh = inputBuffer.inCh + 1;
        fieldFound();
        return false;
    };

    // Main scanner loop. `continue;` jumps back to the top of the body
    // (which immediately advances inCh and reads the next char).
    while (true) {
    inputBuffer.inCh = inputBuffer.inCh + 1;
    if (inputBuffer.inCh > inputBuffer.nOch) {
        // end of line — dispatch on inPhase, always returns
        switch (inPhase) {
            case 1: tokenCode = -8; return;
            case 2: { alphaDoneAndExit(); return; }
            case 3: { numericDoneAndExit(); return; }
        }
    }
    curChar = inputBuffer.iBuf[inputBuffer.inCh];

    // Is a literal coming in
    if (messageLength > 0) {
        message[messageLength] = curChar;
        messageLength = messageLength + 1;
        if (messageLength > 80) { std::printf("\n**** FIELD IS TOO LONG.\n"); { reportFieldErrorAndReturn(); return; } }
        if (curChar != delimiters[3]) continue;
        messageLength = messageLength - 2;
        message[messageLength + 1] = blank;
        tokenCode = -4;
        inputBuffer.inCh = inputBuffer.inCh + 1;
        { fieldFound(); return; }
    }

    // Is a comment coming in
    if (inComment != 0) {
        if (curChar == '$') inComment = 0;
        continue;
    }

    // Blank?
    if (curChar == blank) {
        switch (inPhase) {
            case 1: continue;
            case 2: { alphaDoneAndExit(); return; }
            case 3: { numericDoneAndExit(); return; }
        }
    }

    // Number?
    if (curChar >= num0 && curChar <= num9) {
        if (keyEx < 0) { if (processAlpha()) continue; return; }
        if (inPhase == 1) {
            // Start of numeric field
            decimalSeen = 0;
            val = (double)(curChar - num0);
            digitCount = 1;
            power = 1.0; inPhase = 3;
            continue;
        }
        if (inPhase == 2) {
            // digit char extends current alpha word
            if (wordLength < 8) { word1[wordLength + 1] = curChar; wordLength = wordLength + 1; }
            inPhase = 2;
            continue;
        }
        // inPhase == 3: continuing numeric field
        val = 10.0 * val + (double)(curChar - num0);
        if (decimalSeen != 0) power = power * 10.0;
        digitCount = digitCount + 1;
        continue;
    }

    // Alphabetic?
    if ((curChar >= 'A' && curChar <= 'Z') || curChar == '*' || curChar == '_') {
        if (processAlpha()) continue;
        return;
    }

    // Convert lowercase to uppercase
    if (curChar >= 'a' && curChar <= 'z') {
        curChar = (char)(curChar + 'A' - 'a');
        if (processAlpha()) continue;
        return;
    }

    // Is it a period
    if (curChar == '.') {
        tokenCode = 20;
        if (keyEx < 0) { if (handleSpecialChar()) continue; return; }
        switch (inPhase) {
            case 1: decimalSeen = 1; val = 0.0; power = 1.0; inPhase = 3; continue;
            case 2: { reportInvalidDecimal(); return; }
            case 3:
                // Decimal pt encountered in numeric field
                if (decimalSeen != 0) { reportInvalidDecimal(); return; }
                decimalSeen = 1;
                continue;
        }
    }

    // All other special characters stop a field in progress
    switch (inPhase) {
        case 1: {
            // Check for + - ' $
            bool nextTokenSignal = false;
            for (i = 1; i <= 4; i++) {
                if (curChar == delimiters[i]) {
                    // dispatch on the matched delimiter index;
                    // each case set nextTokenSignal to break for-loop and continue scanner.
                    switch (i) {
                        case 1: sign = +1.0; tokenCode = 7; if (keyEx < 0) { if (handleSpecialChar()) { nextTokenSignal = true; break; } return; } power = 1.0; inPhase = 3; nextTokenSignal = true; break;
                        case 2: sign = -1.0; tokenCode = 6; if (keyEx < 0) { if (handleSpecialChar()) { nextTokenSignal = true; break; } return; } power = 1.0; inPhase = 3; nextTokenSignal = true; break;
                        case 3: messageLength = 1; nextTokenSignal = true; break;
                        case 4: inComment = 1; nextTokenSignal = true; break;
                    }
                    if (nextTokenSignal) break;
                }
            }
            if (nextTokenSignal) continue;
            // Look for tokens
            bool tokenMatch = false;
            for (tokenCode = 1; tokenCode <= 21; tokenCode++) {
                if (curChar == chars[tokenCode]) { tokenMatch = true; break; }
            }
            if (tokenMatch) {
                if (handleSpecialChar()) continue;
                return;
            }
            { std::printf("\n**** UNRECOGNIZABLE CHARACTER '%c' (HEX IS %02X).\n", curChar, (unsigned char)curChar); { reportFieldErrorAndReturn(); return; } }
        }
        case 2: { alphaDoneAndExit(); return; }
        case 3: { numericDoneAndExit(); return; }
    }

    }  // end while (scanner loop)
}


// ============================================================================
// NXINT + entries (nxValue, nxWord, newCard) — Fortran-era multi-entry parser.
// Uses enum dispatch pattern. NXINT/NXVALF/NXHINT/FITKEY entries deleted
// ============================================================================

// Helper: read next card on tokenCode==-8, echo if requested. Returns true to
// continue scanning, false to return EOF (1) from caller.
static bool nxintReadCard() {
    char line[201];
    if (std::fgets(line, sizeof(line), stdin) == nullptr) return false;
    int len = (int)std::strlen(line);
    if (len > 0 && line[len-1] == '\n') { line[--len] = '\0'; }
    if (len > iBufSize) len = iBufSize;
    for (int i = 1; i <= len; i++) inputBuffer.iBuf[i] = line[i-1];
    for (int i = len + 1; i <= iBufSize; i++) inputBuffer.iBuf[i] = ' ';
    for (inputBuffer.nOch = iBufSize; inputBuffer.nOch >= 1; inputBuffer.nOch--) {
        if (inputBuffer.iBuf[inputBuffer.nOch] != ' ') break;
    }
    if (inputBuffer.nOch < 1) inputBuffer.nOch = 1;
    inputBuffer.inCh = 1;
    std::printf("0INPUT... ");
    for (int i = 1; i <= inputBuffer.nOch; i++) std::printf("%c", inputBuffer.iBuf[i]);
    std::printf("\n");
    return true;
}

// nxValue — read next numeric value from input buffer.
// Returns: 0=normal, 1/2/3=alternate (EOF / end-of-record / rescan)
static int nxValue(double& valArg) {
    int tokenCode, intValue, messageLength, inChStart;
    double val;
    char8 charValue;
    char stop;
    char messageBuffer[2];
    while (true) {
        inChStart = inputBuffer.inCh;
        mScan(11, tokenCode, charValue, val, intValue, messageBuffer, messageLength, stop);
        if (tokenCode == -8) {
            if (!nxintReadCard()) return 1;
            continue;
        }
        if (tokenCode == 0) continue;
        // KEY=3 (nxValue): only tokenCode ∈ {-2 (number), -5 (decimal-missing)} hand back a value;
        // every other token (incl. tokenCode==5 end-of-list) folds into a non-zero return code.
        if (tokenCode != -2 && tokenCode != -5) {
            if (tokenCode == 5) return 2;
            inputBuffer.inCh = inChStart;
            return tokenCode < 0 ? 1 : 3;
        }
        // tokenCode is -2 or -5 — stop='/' means slash-terminator → also rewind + return 1.
        if (stop == '/') {
            inputBuffer.inCh = inChStart;
            return 1;
        }
        valArg = val;
        return 0;
    }
}

// nxWord — read next word from input buffer.
int nxWord(char* cvArg) {
    int tokenCode, intValue, messageLength, inChStart;
    double val;
    char8 charValue;
    char stop;
    char messageBuffer[2];
    while (true) {
        inChStart = inputBuffer.inCh;
        mScan(11, tokenCode, charValue, val, intValue, messageBuffer, messageLength, stop);
        if (tokenCode == -8) {
            if (!nxintReadCard()) return 1;
            continue;
        }
        if (tokenCode == 0) continue;
        // KEY=8 (nxWord)
        if (tokenCode != -1) {
            if (tokenCode == 5) return 2;
            inputBuffer.inCh = inChStart;
            return tokenCode < 0 ? 1 : 3;
        }
        std::memcpy(cvArg, charValue.data, 8);
        return 0;
    }
}

// newCard — reset input buffer so the next read fetches a fresh card.
void newCard() {
    inputBuffer.nOch = 0;
    inputBuffer.inCh = 1;
}





// ============================================================================
// InputParser::defineArray — read input values from the input buffer into
// reaction.named under the 8-char key name (formerly DEFINE free function;
// ============================================================================
void InputParser::defineArray(char8 name, Reaction& reaction)
{
    std::vector<double> values;
    double value;

    while (true) {
        int r = nxValue(value);
        if (r == 2) {
            // End of list -- backup for semicolon
            inputBuffer.inCh = inputBuffer.inCh - 1;
        }
        if (r != 0) break;
        values.push_back(value);
    }

    reaction.named.define(std::string_view(name.data, 8), std::move(values));
}


// ============================================================================
// qvScan — Q-value / reaction scan
// ============================================================================
void qvScan(char8* guy, double* eStars, int* nodeVals, int* lVals,
            int* jVals, int* iParities, int& returnCode, Reaction& reaction)
{
    static const char stopS[5] = { ' ', '(', ',', ')', ' ' }; // 1-based
    const char blank = ' ';   // only blank survives (2 uses); the rest were single-use
    char8 blank8(" ");
    char8 charValue;
    double value, dummy;
    int tokenCode, intValue, messageLength;
    char stopChar;
    int stopType;  // set on every for-loop iteration before the post-loop read

    returnCode = 0;
    int inChStart = inputBuffer.inCh;

    // Shared error epilogue: print "WAS PROCESSING FIELD ..." +
    // inChStart..inCh range, then set returnCode=-1.
    auto reportErrorAndReturn = [&]() {
        int inChM1 = inputBuffer.inCh - 1;
        std::printf("      WAS PROCESSING FIELD%2d OF THE REACTION\n", 0);
        std::printf("      REACTION UP TO THE ERROR IS: ");
        for (int ii = inChStart; ii <= inChM1 && ii <= iBufSize; ii++)
            std::printf("%c", inputBuffer.iBuf[ii]);
        std::printf("\n");
        returnCode = -1;
    };

    // tokenCode == -8 (spec ran off the input line) hits the same
    // "COMPLETE REACTION MUST BE ON SAME LINE" message + error epilogue at two
    // scan points (the E*/J inner loop and its /2 bypass re-scan); shared here.
    // Caller still issues the return.
    auto reportIncompleteReaction = [&]() {
        std::printf("\n**** COMPLETE REACTION MUST BE ON SAME LINE AS REACTION KEYWORD.\n");
        reportErrorAndReturn();
    };

    // The E*/J inner loop and its /2-bypass re-scan both do an mScan(18) followed
    // by the same -9 (error) / -8 (incomplete) bail-outs; shared here as a
    // returning guard — returns true if the caller should issue its return.
    auto scanJFieldOrBail = [&]() -> bool {
        mScan(18, tokenCode, charValue, value, intValue, inputBuffer.iBuf, messageLength, stopChar);
        if (tokenCode == -9) { reportErrorAndReturn(); return true; }
        if (tokenCode == -8) { reportIncompleteReaction(); return true; }
        return false;
    };

    // Loop on the reaction participants
    for (int i = 1; i <= 4; i++) {
        guy[i] = blank8;
        nodeVals[i] = NOTDEF_INT;
        lVals[i] = NOTDEF_INT;
        jVals[i] = NOTDEF_INT;
        eStars[i] = reaction.internalState.undefValue;
        iParities[i] = 0;

        mScan(-20, tokenCode, guy[i], value, intValue, inputBuffer.iBuf, messageLength, stopChar);
        if (messageLength > 5) {
            std::printf("\n**** TOO MANY CHARACTERS IN SYMBOL AND ATOMIC MASS.\n");
            { reportErrorAndReturn(); return; }
        }
        if (messageLength == 0) { std::printf("\n**** UNEXPECTED DELIMITOR OR END OF LINE ENCOUNTERED WHEN LOOKING FOR ELEMENT SYMBOL.\n"); { reportErrorAndReturn(); return; } }
        inputBuffer.inCh = inputBuffer.inCh + 1;

        // Only process E*/J when stopped by '(' and not first field
        if (stopChar == stopS[1] && i != 1) {
            // Get E* and J — inner loop, exits when tokenCode==11
            while (true) {
                if (scanJFieldOrBail()) return;
                if (tokenCode == -2) {
                    // It is E*
                    if (eStars[i] != reaction.internalState.undefValue) {
                        std::printf("\n**** ATTEMPTING TO ENTER 2 VALUES FOR E*\n");
                        { reportErrorAndReturn(); return; }
                    }
                    eStars[i] = value;
                    continue;
                }

                if (tokenCode == -5) {
                    // It is J
                    if (jVals[i] != NOTDEF_INT) {
                        std::printf("\n**** ATTEMPTING TO ENTER 2 VALUES FOR J\n");
                        { reportErrorAndReturn(); return; }
                    }
                    jVals[i] = 2 * intValue;
                    if (stopChar == '/') {
                        // Bypass /2
                        inputBuffer.inCh = inputBuffer.inCh + 1;
                        jVals[i] = intValue;
                        if (scanJFieldOrBail()) return;
                        if (!(tokenCode == -5 && intValue == 2)) {
                            std::printf("\n**** INVALID J FIELD\n");
                            { reportErrorAndReturn(); return; }
                        }
                    }

                    // Bypass terminal + or -
                    if (stopChar == '+') {
                        iParities[i] = +1;
                        inputBuffer.inCh = inputBuffer.inCh + 1;
                    } else if (stopChar == '-') {
                        iParities[i] = -1;
                        inputBuffer.inCh = inputBuffer.inCh + 1;
                    }
                    continue;
                }

                // non-J excited-state token
                if (tokenCode == -3) {
                    // Level descriptor
                    if (lVals[i] != NOTDEF_INT) {
                        std::printf("\n**** ATTEMPTING TO ENTER TWO LEVEL DESCRIPTORS\n");
                        { reportErrorAndReturn(); return; }
                    }
                    lVals[i] = intValue;
                    nodeVals[i] = (int)value;
                    continue;
                }

                if (tokenCode == 11) {
                    stopChar = inputBuffer.iBuf[inputBuffer.inCh];
                    inputBuffer.inCh = inputBuffer.inCh + 1;
                    break;
                }

                std::printf("\n**** UNRECOGNIZABLE EXCITED STATE DESCRIPTOR (TOKEN =%4d)\n", tokenCode);
                { reportErrorAndReturn(); return; }
            }
        }

        // Bypass blanks
        if (stopChar == blank && i != 2 && i != 4) {
            mScan(-20, tokenCode, charValue, dummy, intValue, inputBuffer.iBuf, messageLength, stopChar);
            stopChar = inputBuffer.iBuf[inputBuffer.inCh - 1];
            if (messageLength != 0) {
                std::printf("\n**** UNEXPECTED FIELD ENCOUNTERED WHEN LOOKING FOR DELIMITOR\n");
                { reportErrorAndReturn(); return; }
            }
        }

        if ((i == 1 || i == 3) && stopChar != stopS[i]) {
            std::printf("\n**** MISSING PARENTHESIS.\n");
            { reportErrorAndReturn(); return; }
        }
        stopType = (stopChar == blank || stopChar == ',') ? 1 : 2;
        if (i == 2 && stopType != 1) {
            std::printf("\n**** INVALID DELIMITOR AFTER INCOMING PROJECTILE.\n");
            { reportErrorAndReturn(); return; }
        }
    }

    if (stopType != 1 && inputBuffer.inCh <= inputBuffer.nOch) {
        std::printf("\n**** ONLY 4 PARTICIPANTS IN REACTION ALLOWED.\n");
        reportErrorAndReturn();
    }
}


// keyword-list-with-values dump. Never called from any C++ caller.
// Was the only consumer of positional FLOAT_arr(i)/INTGER_arr_f(i)/
// SWITCH_arr indexing for keyword display (real users of these fields
// access them by name).
