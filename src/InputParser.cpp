// InputParser.cpp — InputParser implementation: parses Ptolemy .in files and
// applies them to the reaction's COMMON-derived structs.

#include "InputParser.h"
#include "ptolemy_types.h"
#include "Reaction.h"
#include "Elastic.h"
#include "Inelastic.h"
#include "Transfer.h"
#include "DwbaInputExpander.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>

// File-scope alias to the InputParser-owned scanner buffer. Was an
namespace { auto& inputBuffer = InputParser::buffer(); }

// ============================================================================
// Static utility helpers
// ============================================================================

std::string InputParser::toUpper(const std::string& s) {
    std::string r = s;
    for (char& c : r) c = (char)std::toupper((unsigned char)c);
    return r;
}

std::string InputParser::trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Split on whitespace and commas; strip comments after '!' or '$'.
std::vector<std::string> InputParser::tokenize(const std::string& s) {
    std::vector<std::string> rawToks;
    std::string line = s;
    // strip '!' comment
    size_t ci = line.find('!');
    if (ci != std::string::npos) line = line.substr(0, ci);
    // strip '$' comment (Ptolemy uses $ for inline comments)
    ci = line.find('$');
    if (ci != std::string::npos) line = line.substr(0, ci);
    std::istringstream ss(line);
    std::string t;
    while (ss >> t) {
        // split on commas: replace comma with space and re-split
        bool hasComma = false;
        for (char& c : t) if (c == ',') { c = ' '; hasComma = true; }
        if (hasComma) {
            std::istringstream ss2(t);
            std::string t2;
            while (ss2 >> t2) rawToks.push_back(t2);
        } else {
            rawToks.push_back(t);
        }
    }

    // Pattern: token[i] ends without '=', token[i+1] == "=", token[i+2] is a value.
    std::vector<std::string> toks;
    for (size_t i = 0; i < rawToks.size(); i++) {
        // Check if next token is bare "="
        if (i + 2 < rawToks.size() && rawToks[i+1] == "=") {
            // Merge: "KEY" + "=" + "VALUE" → "KEY=VALUE"
            toks.push_back(rawToks[i] + "=" + rawToks[i+2]);
            i += 2;
        } else if (i + 1 < rawToks.size() && rawToks[i+1] == "=") {
            // "KEY" + "=" with no following value → bare "KEY="
            toks.push_back(rawToks[i] + "=");
            i += 1;
        } else {
            toks.push_back(rawToks[i]);
        }
    }
    return toks;
}

// Join toks[start..] into a single string separated by single spaces.
std::string InputParser::joinTokensFrom(const std::vector<std::string>& toks, size_t start) {
    std::string rest;
    for (size_t i = start; i < toks.size(); i++) {
        if (i > start) rest += " ";
        rest += toks[i];
    }
    return rest;
}

// Try to match "KEY=VALUE" in a single token (case-insensitive).
// Also accepts "KEY=" with empty value (trailing '=') — caller handles next-token fallback.
bool InputParser::matchKeyValue(const std::string& tok, const std::string& key, std::string& valStr) {
    std::string u = toUpper(tok);
    std::string k = toUpper(key) + "=";
    if (u.rfind(k, 0) != 0) return false;
    valStr = tok.substr(k.size());
    return !valStr.empty();  // bare "KEY=" needs next-token
}

bool InputParser::getDouble(const std::string& tok, const std::string& key, double& out) {
    std::string valStr;
    if (!matchKeyValue(tok, key, valStr)) return false;
    try { out = std::stod(valStr); return true; }
    catch (...) { return false; }
}

bool InputParser::getInt(const std::string& tok, const std::string& key, int& out) {
    std::string valStr;
    if (!matchKeyValue(tok, key, valStr)) return false;
    try { out = std::stoi(valStr); return true; }
    catch (...) { return false; }
}

// Check if tok is exactly "KEY=" (bare equals, value follows in next token)
bool InputParser::isBareKey(const std::string& tok, const std::string& key) {
    std::string u = toUpper(tok);
    std::string k = toUpper(key) + "=";
    return (u == k);
}

// Parse half-integer spin, e.g. "3/2-" → twoJ=3 (parity discarded).
// Also handles "0+" "3/2" "7/2-" "1/2+" "0" etc.
bool InputParser::parseHalfInt(const std::string& tok, int& twoJ) {
    std::string s = tok;
    // strip trailing parity (consumed but discarded — no caller reads it)
    if (!s.empty() && (s.back() == '+' || s.back() == '-')) s.pop_back();
    // check for fraction
    size_t slash = s.find('/');
    if (slash != std::string::npos) {
        int num = std::stoi(s.substr(0, slash));
        int den = std::stoi(s.substr(slash + 1));
        if (den == 2) { twoJ = num; return true; }
        return false;
    }
    // integer: J = N, twoJ = 2*N
    try { twoJ = 2 * std::stoi(s); return true; }
    catch (...) { return false; }
}

// ============================================================================
// parse() — main entry point
// ============================================================================
bool InputParser::parse(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::fprintf(stderr, "InputParser: cannot open '%s'\n", path.c_str());
        return false;
    }
    return parse(f);
}

bool InputParser::parse(std::istream& f) {
    section_   = Section::TOPLEVEL;
    inBlock_  = false;

    std::string line;
    int lineNo = 0;
    while (std::getline(f, line)) {
        ++lineNo;
        std::string tl = trim(line);
        if (tl.empty()) continue;
        // comment lines ($ and ! are Ptolemy comment chars; # also)
        if (tl[0] == '!' || tl[0] == '#' || tl[0] == '$') continue;

        if (!parseLine(tl)) {
            std::fprintf(stderr, "InputParser: error at line %d: %s\n", lineNo, tl.c_str());
            // non-fatal — keep going
        }
    }
    return true;
}

// ============================================================================
// hasDwbaKeyword() — true if the first non-blank, non-comment line is "DWBA".
// ============================================================================
bool InputParser::hasDwbaKeyword(const std::string& content) {
    std::istringstream probe(content);
    std::string line;
    while (std::getline(probe, line)) {
        std::string t = trim(line);
        if (t.empty()) continue;
        if (t[0] == '!' || t[0] == '#' || t[0] == '$') continue;
        return toUpper(t) == "DWBA";
    }
    return false;
}

// ============================================================================
// parseDwba() — expand a human-readable DWBA reaction description to a Ptolemy
// deck (via DwbaInputExpander) and parse it. A leading "DWBA" keyword line, if
// present, is stripped before expansion.
// ============================================================================
bool InputParser::parseDwba(const std::string& content) {
    std::fprintf(stderr, "ptolemy: DWBA input detected, expanding...\n");
    std::string deck = content;
    std::istringstream probe(content);
    std::string line;
    while (std::getline(probe, line)) {
        std::string t = trim(line);
        if (t.empty() || t[0] == '!' || t[0] == '#' || t[0] == '$') continue;
        if (toUpper(t) == "DWBA") {
            size_t p = content.find(line);
            if (p != std::string::npos)
                deck = content.substr(0, p) + content.substr(p + line.size());
        }
        break;
    }
    std::string expanded = DwbaExpander::expand(deck);
    std::istringstream in(expanded);
    return parse(in);
}

// ============================================================================
// parseFromArgs() — CLI entry.
// Recognizes optional flags before the input file path:
//   --fixedLS  use physics-standard <L*S> spin-orbit coupling.
//   --dwba     force DWBA input mode (human-readable reaction description).
//   --help     print usage.
// If no input file given (only flags or nothing), reads from stdin.
//
// DWBA auto-detection:
//   * file argument — slurp and use content-based detection (handles leading
//     comments / "DWBA" keyword); native decks are parsed by parse(path).
//   * stdin — detection MUST NOT extract from std::cin: the legacy linkule
//     memory pool used by the elastic path is sensitive to heap/stream layout
//     and any pre-read shifts it into a failing state. std::cin.peek() inspects
//     the next byte without extracting. A DWBA reaction line starts with a
//     mass-number digit; native Ptolemy decks start with a keyword letter or
//     comment char. Leading digit => DWBA, else native (parse std::cin). For a
//     DWBA stdin stream with leading comments or the "DWBA" keyword, pass --dwba
//     or use a file argument.
// ============================================================================
bool InputParser::parseFromArgs(int argc, char** argv) {
    bool forceDwba = false;
    int firstNonFlag = 1;
    while (firstNonFlag < argc) {
        const char* a = argv[firstNonFlag];
        if (std::strcmp(a, "--fixedLS") == 0 || std::strcmp(a, "--fixed-ls") == 0) {
            d_.cliFixedLS = true;
            ++firstNonFlag;
            continue;
        }
        if (std::strcmp(a, "--dwba") == 0 || std::strcmp(a, "--DWBA") == 0) {
            forceDwba = true;
            ++firstNonFlag;
            continue;
        }
        if (std::strcmp(a, "--help") == 0 || std::strcmp(a, "-h") == 0) {
            std::fprintf(stderr,
                "Usage: %s [--fixedLS] [--dwba] [input_file]\n"
                "  --fixedLS    use physics-standard <L*S> spin-orbit coupling\n"
                "               (default: Cleopatra-faithful sigma*L convention)\n"
                "  --dwba       force DWBA input mode (human-readable reaction\n"
                "               description); otherwise auto-detected.\n"
                "  input_file   path to input deck; if omitted, reads stdin.\n"
                "               Accepts a native Ptolemy deck or a DWBA reaction\n"
                "               description.\n",
                argv[0]);
            return false;
        }
        break;  // first non-flag token → input file path
    }

    // ---- file argument ----
    if (firstNonFlag < argc) {
        const char* srcName = argv[firstNonFlag];
        std::ifstream f(srcName);
        if (!f.is_open()) {
            std::fprintf(stderr, "ptolemy: cannot open '%s'\n", srcName);
            return false;
        }
        std::stringstream ss; ss << f.rdbuf();
        std::string content = ss.str();
        bool ok;
        if (forceDwba || DwbaExpander::looksLikeDwba(content) || hasDwbaKeyword(content)) {
            ok = parseDwba(content);
        } else {
            ok = parse(srcName);  // unchanged native path
        }
        if (!ok) { std::fprintf(stderr, "ptolemy: parse failed on '%s'\n", srcName); return false; }
        return true;
    }

    // ---- stdin ----
    bool isDwba = forceDwba;
    if (!forceDwba) {
        int c = std::cin.peek();   // inspect only; does not extract
        if (c != EOF && std::isdigit((unsigned char)c)) isDwba = true;
    }
    if (isDwba) {
        std::string content((std::istreambuf_iterator<char>(std::cin)),
                             std::istreambuf_iterator<char>());
        if (!parseDwba(content)) {
            std::fprintf(stderr, "ptolemy: parse failed on stdin\n");
            return false;
        }
        return true;
    }
    if (!parse(std::cin)) {        // unchanged native path
        std::fprintf(stderr, "ptolemy: parse failed on stdin\n");
        return false;
    }
    return true;
}

// ============================================================================
// parseLine() — dispatch based on current section and keyword
// ============================================================================
bool InputParser::parseLine(const std::string& raw) {
    // Strip $ inline comments before processing
    std::string clean = raw;
    size_t dc = clean.find('$');
    if (dc != std::string::npos) clean = clean.substr(0, dc);
    clean = trim(clean);
    if (clean.empty()) return true;

    std::string u = toUpper(clean);
    auto toks = tokenize(clean);
    if (toks.empty()) return true;

    std::string kw = toUpper(toks[0]);

    // ---- Semicolon or bare ';': end-of-block trigger ----
    if (kw == ";" || clean == ";") {
        inBlock_ = false;
        section_  = Section::TOPLEVEL;
        return true;
    }

    // ---- RESET ----
    if (kw == "RESET") {
        d_ = ParsedInput{};
        section_ = Section::TOPLEVEL;
        inBlock_ = false;
        return true;
    }

    // ---- FIXEDLS (input-deck switch) ----
    // Bare keyword that selects the physics-standard <L*S> spin-orbit coupling
    // for this run, same effect as the --fixedLS CLI flag. See README
    // "Spin-orbit convention" for the relationship between input Vso under
    // the two modes. May appear anywhere at top level.
    if (kw == "FIXEDLS" || kw == "FIXED-LS") {
        d_.cliFixedLS = true;
        return true;
    }

    // ---- END / RETURN / QUIT / STOP ----
    if (kw == "END" || kw == "RETURN" || kw == "QUIT" || kw == "STOP") {
        return true;
    }

    // ---- PRINT level ----
    if (kw == "PRINT") {
        if (toks.size() > 1) {
            try { d_.printLevel = std::stoi(toks[1]); } catch (...) {}
        } else {
            d_.printLevel = 1;
        }
        return true;
    }

    // ---- CHANNEL: keyword (elastic scattering format) ----
    // Format: "CHANNEL: A + 208PB" or "CHANNEL: P + 208PB"
    if (kw == "CHANNEL:" || kw == "CHANNEL") {
        // Parse "projectile + target" from the rest
        // toks after CHANNEL: keyword are like: ["A", "+", "208PB"]
        // or after stripping colon from first tok
        std::string rest;
        size_t cPos = raw.find(':');
        if (cPos != std::string::npos)
            rest = trim(raw.substr(cPos + 1));
        else
            rest = joinTokensFrom(toks, 1);
        // Strip $ from rest too
        size_t dPos = rest.find('$');
        if (dPos != std::string::npos) rest = rest.substr(0, dPos);
        rest = trim(rest);
        // Parse "PROJ + TARGET" (+ is separator)
        // Use whitespace tokenization
        auto rtoks = tokenize(rest);
        // Remove "+" tokens
        std::vector<std::string> parts;
        for (auto& t : rtoks) if (t != "+") parts.push_back(t);
        if (parts.size() >= 1) d_.channelProjectile = toUpper(parts[0]);
        if (parts.size() >= 2) d_.channelTarget       = toUpper(parts[1]);
        // After CHANNEL keyword, top-level OM params belong to INCOMING
        section_  = Section::INCOMING;
        inBlock_ = true;
        d_.incoming.set = true;
        d_.isElastic = true;
        return true;
    }

    // ---- ELASTIC SCATTERING keyword ----
    if (u == "ELASTIC SCATTERING" || u == "ELASTIC") {
        d_.isElastic = true;
        // End of OM params block for elastic
        inBlock_ = false;
        section_  = Section::TOPLEVEL;
        return true;
    }

    // ---- Block headers ----
    if (kw == "PROJECTILE" || kw == "PROJECTI") {
        section_ = Section::PROJECTILE;
        inBlock_ = true;
        d_.projectileBs.set = true;
        return true;
    }
    if (kw == "TARGET" || kw == "ZIELSCHE") {
        section_ = Section::TARGET;
        inBlock_ = true;
        d_.targetBs.set = true;
        return true;
    }
    if (kw == "INCOMING" || kw == "EINGANG") {
        section_ = Section::INCOMING;
        inBlock_ = true;
        d_.incoming.set = true;
        return true;
    }
    if (kw == "OUTGOING" || kw == "AUSGANG") {
        section_ = Section::OUTGOING;
        inBlock_ = true;
        d_.outgoing.set = true;
        return true;
    }

    // ---- REACTION: line ----
    if (kw == "REACTION:" || kw == "REACTION" || kw == "REAKTION") {
        // rest of line after the keyword
        size_t pos = raw.find(':');
        std::string rest;
        if (pos != std::string::npos)
            rest = trim(raw.substr(pos + 1));
        else if (toks.size() > 1)
            // rebuild rest from toks[1..]
            rest = joinTokensFrom(toks, 1);
        // Strip $ comments from reaction line
        size_t dPos = rest.find('$');
        if (dPos != std::string::npos) rest = rest.substr(0, dPos);
        return parseReactionLine(trim(rest));
    }

    // ---- PARAMETERSET ----
    if (kw == "PARAMETERSET" || kw == "PARAMETERE" || kw == "PARAMETER") {
        // rest of line
        std::string rest = joinTokensFrom(toks, 1);
        return parseParametersetLine(rest);
    }

    // ---- header: / TITLE: ignored ----
    if (kw == "HEADER:" || kw == "HEADER" || kw == "TITLE") return true;

    // ---- ASYMPTOPIA as standalone keyword (not just ASYMPTOPIA=X) ----
    // Format: "ASYMPTOPIA=35" (handled in global keyvals)
    // Also handled below in parseGlobalKeyvals via getDouble

    // ---- Dispatch to section parsers ----
    // Always check for global-only keywords first (eLab, ANGLEMIN, etc.).
    // These can appear anywhere (e.g. between CHANNEL: and ELASTIC SCATTERING).
    {
        // Extract the keyword part of toks[0] (before '=' or full token)
        std::string uFirst = toUpper(toks[0]);
        size_t eqPos = uFirst.find('=');
        if (eqPos != std::string::npos) uFirst = uFirst.substr(0, eqPos);
        // Known global-only keywords (not valid OM params)
        static const char* globalOnly[] = {
            "ELAB", "ECM", "ANGLEMIN", "ANGLEMAX", "ANGLESTEP", "ANGLESTE",
            "LMAX", "LMIN", "LSTEP", "MAXLEXTRAP", "LMAXADD",
            "ASYMPTOPIA", "STEPSPER", "JBIGA", "BELX", "BETACOUL",
            "SPFACT", "SPFACP", "PRINT", nullptr
        };
        for (int gi = 0; globalOnly[gi]; gi++) {
            if (uFirst == globalOnly[gi]) {
                return parseGlobalKeyvals(clean);
            }
        }
    }

    switch (section_) {
        case Section::PROJECTILE:
            return parseKeyvals(clean, d_.projectileBs);
        case Section::TARGET:
            return parseKeyvals(clean, d_.targetBs);
        case Section::INCOMING:
            return parseKeyvals(clean, d_.incoming);
        case Section::OUTGOING:
            return parseKeyvals(clean, d_.outgoing);
        case Section::TOPLEVEL:
        default:
            return parseGlobalKeyvals(clean);
    }
}

// ============================================================================
// parseReactionLine() — parse "48CA(D,P)49CA(3/2- 0.000) ELAB=30"
// Also handles "ELAB= 14.780" (space after '=')
// ============================================================================
bool InputParser::parseReactionLine(const std::string& rest) {
    // Store the full reaction string for REACTN
    d_.reactionLine = rest;

    // Parse eLab=X from the token list
    // Handle both "ELAB=30" and "ELAB= 14.780" (split token) forms
    auto toks = tokenize(rest);
    for (size_t i = 0; i < toks.size(); i++) {
        double v;
        if (getDouble(toks[i], "ELAB", v)) { d_.eLab = v; continue; }
        if (getDouble(toks[i], "ECM",  v)) { d_.eCm  = v; continue; }
        // Handle bare "ELAB=" followed by value in next token
        if (isBareKey(toks[i], "ELAB") && i + 1 < toks.size()) {
            try { d_.eLab = std::stod(toks[i+1]); i++; } catch (...) {}
            continue;
        }
        if (isBareKey(toks[i], "ECM") && i + 1 < toks.size()) {
            try { d_.eCm  = std::stod(toks[i+1]); i++; } catch (...) {}
            continue;
        }
    }

    // Parse the reaction notation: A(a,b)B(jProj Ex)
    // Find the main parentheses structure
    std::string s = trim(rest);
    // Strip eLab=... from the string for nucleus parsing
    size_t ePos = toUpper(s).find("ELAB=");
    if (ePos != std::string::npos) s = trim(s.substr(0, ePos));

    // validate the structural format so callers get a true/false parse result.
    size_t p1 = s.find('(');
    if (p1 == std::string::npos) return false;
    size_t p2 = s.find(')', p1);
    if (p2 == std::string::npos) return false;
    std::string ab = s.substr(p1 + 1, p2 - p1 - 1);
    return ab.find(',') != std::string::npos;
}

// ============================================================================
// parseParametersetLine() — "dpsb r0Target [key=val ...]"
// Handles: r0Target, lMax=N, lMin=N, lStep=N, maxLExtrap=N, asymptopia=X, lMaxAdditional=N
// ============================================================================
bool InputParser::parseParametersetLine(const std::string& rest) {
    auto toks = tokenize(rest);
    if (toks.empty()) return false;
    d_.parameterSet = toUpper(toks[0]);
    for (size_t i = 1; i < toks.size(); i++) {
        std::string u = toUpper(toks[i]);
        if (u == "R0TARGET" || u == "R0TARG") {
            d_.r0Target = true; continue;
        }
        double dv; int iv;
        if (getDouble(toks[i], "ASYMPTOPIA", dv)) { d_.asymptopia = dv; continue; }
        if (getInt(toks[i], "LMAX",       iv)) { d_.lMax       = iv; continue; }
        if (getInt(toks[i], "LMIN",       iv)) { d_.lMin       = iv; continue; }
        if (getInt(toks[i], "LSTEP",      iv)) { d_.lStep      = iv; continue; }
        if (getInt(toks[i], "MAXLEXTRAP", iv)) { d_.maxLExtrap = iv; continue; }
        if (getInt(toks[i], "LMAXADD",   iv)) { d_.lMaxAdditional    = iv; continue; }
        // Bare key= followed by next token
        if (isBareKey(toks[i], "ASYMPTOPIA") && i+1 < toks.size()) {
            try { d_.asymptopia = std::stod(toks[i+1]); i++; } catch (...) {} continue;
        }
        if (isBareKey(toks[i], "LMAX") && i+1 < toks.size()) {
            try { d_.lMax = std::stoi(toks[i+1]); i++; } catch (...) {} continue;
        }
        if (isBareKey(toks[i], "LMIN") && i+1 < toks.size()) {
            try { d_.lMin = std::stoi(toks[i+1]); i++; } catch (...) {} continue;
        }
        if (isBareKey(toks[i], "LSTEP") && i+1 < toks.size()) {
            try { d_.lStep = std::stoi(toks[i+1]); i++; } catch (...) {} continue;
        }
        if (isBareKey(toks[i], "MAXLEXTRAP") && i+1 < toks.size()) {
            try { d_.maxLExtrap = std::stoi(toks[i+1]); i++; } catch (...) {} continue;
        }
        if (isBareKey(toks[i], "LMAXADD") && i+1 < toks.size()) {
            try { d_.lMaxAdditional = std::stoi(toks[i+1]); i++; } catch (...) {} continue;
        }
    }
    return true;
}

// ============================================================================
// parseKeyvals() — bound state block (PROJECTILE / TARGET)
// ============================================================================
bool InputParser::parseKeyvals(const std::string& raw, ParsedInput::BSParams& bs) {
    auto toks = tokenize(raw);
    for (size_t i = 0; i < toks.size(); i++) {
        const std::string& tok = toks[i];
        std::string u = toUpper(tok);
        double dv; int iv;

        // wavefunction <type>
        if (u == "WAVEFUNCTION" || u == "WAVEFUNC" || u == "WAVEF") {
            // next token is the wf type
            if (i + 1 < toks.size()) { bs.wavefunction = toUpper(toks[++i]); }
            continue;
        }
        if (u == "AV18" || u == "RCWFN") {
            bs.wavefunction = u; continue;
        }

        if (getDouble(tok, "V",    dv)) { bs.V    = dv; bs.hasV    = true; continue; }
        if (getDouble(tok, "R0",   dv)) { bs.r0   = dv; bs.hasR0   = true; continue; }
        if (getDouble(tok, "A",    dv)) { bs.a    = dv; bs.hasA    = true; continue; }
        if (getDouble(tok, "VSO",  dv)) { bs.vSo  = dv; bs.hasVSO  = true; continue; }
        if (getDouble(tok, "RSO0", dv)) { bs.rSo0 = dv; bs.hasRSO0 = true; continue; }
        if (getDouble(tok, "ASO",  dv)) { bs.aSo  = dv; bs.hasASO  = true; continue; }
        if (getDouble(tok, "RC0",  dv)) { bs.rC0  = dv; bs.hasRC0  = true; continue; }

        if (getInt(tok, "L",     iv)) { bs.l     = iv; continue; }
        if (getInt(tok, "NODES", iv)) { bs.nodes = iv; continue; }

        // Handle bare "KEY=" followed by value in next token
        auto tryBareDouble = [&](const std::string& key, double& destination, bool& flag) -> bool {
            if (isBareKey(tok, key) && i + 1 < toks.size()) {
                try { destination = std::stod(toks[i+1]); flag = true; i++; return true; } catch (...) {}
            }
            return false;
        };
        auto tryBareInteger = [&](const std::string& key, int& destination) -> bool {
            if (isBareKey(tok, key) && i + 1 < toks.size()) {
                try { destination = std::stoi(toks[i+1]); i++; return true; } catch (...) {}
            }
            return false;
        };
        if (tryBareDouble("V",    bs.V,    bs.hasV))    continue;
        if (tryBareDouble("R0",   bs.r0,   bs.hasR0))   continue;
        if (tryBareDouble("A",    bs.a,    bs.hasA))    continue;
        if (tryBareDouble("VSO",  bs.vSo,  bs.hasVSO))  continue;
        if (tryBareDouble("RSO0", bs.rSo0, bs.hasRSO0)) continue;
        if (tryBareDouble("ASO",  bs.aSo,  bs.hasASO))  continue;
        if (tryBareDouble("RC0",  bs.rC0,  bs.hasRC0))  continue;
        if (tryBareInteger("L",     bs.l))      continue;
        if (tryBareInteger("NODES", bs.nodes))  continue;

        // jp=3/2 or jp=7/2-  (half-integer J); JP= and J= share the same body
        auto parseBsJ = [&](const std::string& jStr) {
            int twoJ = 0;
            if (!jStr.empty()) {
                if (parseHalfInt(jStr, twoJ)) bs.j = twoJ / 2.0;
            } else if (i + 1 < toks.size()) {
                if (parseHalfInt(toUpper(toks[++i]), twoJ)) bs.j = twoJ / 2.0;
            }
        };
        if (u.rfind("JP=", 0) == 0) { parseBsJ(u.substr(3)); continue; }
        if (u.rfind("J=",  0) == 0) { parseBsJ(u.substr(2)); continue; }

        // jResidual=0+ (for TARGET)
        if (u.rfind("JBIGA=", 0) == 0) {
            std::string jStr = u.substr(6);
            int twoJ = 0;
            if (!jStr.empty()) {
                if (parseHalfInt(jStr, twoJ)) {
                    d_.jResidual = twoJ;
                }
            } else if (i + 1 < toks.size()) {
                if (parseHalfInt(toUpper(toks[++i]), twoJ)) {
                    d_.jResidual = twoJ;
                }
            }
            continue;
        }
    }

    return true;
}

// ============================================================================
// parseKeyvals() — optical model block (INCOMING / OUTGOING)
// Handles both "KEY=VALUE" and "KEY= VALUE" (bare-key + next-token) forms.
// ============================================================================
bool InputParser::parseKeyvals(const std::string& raw, ParsedInput::OMParams& om) {
    auto toks = tokenize(raw);
    for (size_t i = 0; i < toks.size(); i++) {
        const std::string& tok = toks[i];
        double dv;
        if (getDouble(tok, "V",    dv)) { om.V    = dv; om.hasV    = true; continue; }
        if (getDouble(tok, "R0",   dv)) { om.r0   = dv; om.hasR0   = true; continue; }
        if (getDouble(tok, "A",    dv)) { om.a    = dv; om.hasA    = true; continue; }
        if (getDouble(tok, "VI",   dv)) { om.vI   = dv; om.hasVI   = true; continue; }
        if (getDouble(tok, "RI0",  dv)) { om.rI0  = dv; om.hasRI0  = true; continue; }
        if (getDouble(tok, "AI",   dv)) { om.aI   = dv; om.hasAI   = true; continue; }
        if (getDouble(tok, "VSI",  dv)) { om.vSi  = dv; om.hasVSI  = true; continue; }
        if (getDouble(tok, "rSi0", dv)) { om.rSi0 = dv; om.hasRSI0 = true; continue; }
        if (getDouble(tok, "ASI",  dv)) { om.aSi  = dv; om.hasASI  = true; continue; }
        if (getDouble(tok, "VSO",  dv)) { om.vSo  = dv; om.hasVSO  = true; continue; }
        if (getDouble(tok, "RSO0", dv)) { om.rSo0 = dv; om.hasRSO0 = true; continue; }
        if (getDouble(tok, "ASO",  dv)) { om.aSo  = dv; om.hasASO  = true; continue; }
        if (getDouble(tok, "VSOI", dv)) { om.vSoi = dv; om.hasVSOI = true; continue; }
        if (getDouble(tok, "RSOI0",dv)) { om.rSoi0= dv; om.hasRSOI0= true; continue; }
        if (getDouble(tok, "ASOI", dv)) { om.aSoi = dv; om.hasASOI = true; continue; }
        if (getDouble(tok, "RC0",  dv)) { om.rC0  = dv; om.hasRC0  = true; continue; }
        // Also accept "W" as alias for vI
        if (getDouble(tok, "W",    dv)) { om.vI   = dv; om.hasVI   = true; continue; }

        // Handle bare "KEY=" followed by value in next token
        auto tryBare = [&](const std::string& key, double& destination, bool& flag) -> bool {
            if (isBareKey(tok, key) && i + 1 < toks.size()) {
                try { destination = std::stod(toks[i+1]); flag = true; i++; return true; } catch (...) {}
            }
            return false;
        };
        if (tryBare("V",    om.V,    om.hasV))    continue;
        if (tryBare("R0",   om.r0,   om.hasR0))   continue;
        if (tryBare("A",    om.a,    om.hasA))    continue;
        if (tryBare("VI",   om.vI,   om.hasVI))   continue;
        if (tryBare("RI0",  om.rI0,  om.hasRI0))  continue;
        if (tryBare("AI",   om.aI,   om.hasAI))   continue;
        if (tryBare("VSI",  om.vSi,  om.hasVSI))  continue;
        if (tryBare("rSi0", om.rSi0, om.hasRSI0)) continue;
        if (tryBare("ASI",  om.aSi,  om.hasASI))  continue;
        if (tryBare("VSO",  om.vSo,  om.hasVSO))  continue;
        if (tryBare("RSO0", om.rSo0, om.hasRSO0)) continue;
        if (tryBare("ASO",  om.aSo,  om.hasASO))  continue;
        if (tryBare("VSOI", om.vSoi, om.hasVSOI)) continue;
        if (tryBare("RSOI0",om.rSoi0,om.hasRSOI0))continue;
        if (tryBare("ASOI", om.aSoi, om.hasASOI)) continue;
        if (tryBare("RC0",  om.rC0,  om.hasRC0))  continue;
        if (tryBare("W",    om.vI,   om.hasVI))   continue;
    }
    return true;
}

// ============================================================================
// parseGlobalKeyvals() — top-level key=value pairs (lMax, elab, angles, etc.)
// Handles both "KEY=VALUE" and "KEY= VALUE" (space-split) forms.
// ============================================================================
bool InputParser::parseGlobalKeyvals(const std::string& raw) {
    auto toks = tokenize(raw);
    for (size_t i = 0; i < toks.size(); i++) {
        const std::string& tok = toks[i];
        std::string u = toUpper(tok);
        double dv; int iv;

        // --- Handle "KEY=VALUE" single-token forms ---
        if (getInt(tok,    "LMAX",       iv)) { d_.lMax       = iv; continue; }
        if (getInt(tok,    "LMIN",       iv)) { d_.lMin       = iv; continue; }
        if (getInt(tok,    "LSTEP",      iv)) { d_.lStep      = iv; continue; }
        if (getInt(tok,    "MAXLEXTRAP", iv)) { d_.maxLExtrap = iv; continue; }
        if (getInt(tok,    "LMAXADD",   iv)) { d_.lMaxAdditional    = iv; continue; }
        if (getInt(tok,    "PRINT",      iv)) { d_.printLevel= iv; continue; }
        if (getDouble(tok, "ASYMPTOPIA", dv)) { d_.asymptopia = dv; continue; }
        if (getDouble(tok, "STEPSPER",   dv)) { d_.stepsPer   = dv; continue; }
        if (getDouble(tok, "ELAB",       dv)) { d_.eLab       = dv; continue; }
        if (getDouble(tok, "ECM",        dv)) { d_.eCm        = dv; continue; }
        if (getDouble(tok, "ANGLEMIN",   dv)) { d_.angleMin  = dv; d_.hasAngleMin   = true; continue; }
        if (getDouble(tok, "ANGLEMAX",   dv)) { d_.angleMax  = dv; d_.hasAngleMax   = true; continue; }
        if (getDouble(tok, "ANGLESTEP",  dv)) { d_.angleStep = dv; d_.hasAngleStep = true; continue; }
        if (getDouble(tok, "ANGLESTE",   dv)) { d_.angleStep = dv; d_.hasAngleStep = true; continue; }
        if (getDouble(tok, "BELX",       dv)) { d_.belx       = dv; continue; }
        if (getDouble(tok, "BETACOUL",   dv)) { d_.betaCoul   = dv; continue; }
        if (getDouble(tok, "SPFACT",     dv)) { d_.spFactor     = dv; continue; }
        if (getDouble(tok, "SPFACP",     dv)) { d_.spFactorProj     = dv; continue; }

        // --- Handle bare "KEY=" followed by value in next token ---
        auto tryBareInteger = [&](const std::string& key, int& destination) -> bool {
            if (isBareKey(tok, key) && i + 1 < toks.size()) {
                try { destination = std::stoi(toks[i+1]); i++; return true; } catch (...) {}
            }
            return false;
        };
        auto tryBareDouble = [&](const std::string& key, double& destination) -> bool {
            if (isBareKey(tok, key) && i + 1 < toks.size()) {
                try { destination = std::stod(toks[i+1]); i++; return true; } catch (...) {}
            }
            return false;
        };

        if (tryBareDouble("ELAB",       d_.eLab))       continue;
        if (tryBareDouble("ECM",        d_.eCm))        continue;
        if (tryBareDouble("ASYMPTOPIA", d_.asymptopia)) continue;
        if (tryBareDouble("STEPSPER",   d_.stepsPer))   continue;
        if (tryBareDouble("ANGLEMIN",   d_.angleMin))  { d_.hasAngleMin   = true; continue; }
        if (tryBareDouble("ANGLEMAX",   d_.angleMax))  { d_.hasAngleMax   = true; continue; }
        if (tryBareDouble("ANGLESTEP",  d_.angleStep)) { d_.hasAngleStep = true; continue; }
        if (tryBareDouble("BELX",       d_.belx))       continue;
        if (tryBareInteger("LMAX",       d_.lMax))       continue;
        if (tryBareInteger("LMIN",       d_.lMin))       continue;
        if (tryBareInteger("LSTEP",      d_.lStep))      continue;
        if (tryBareInteger("MAXLEXTRAP", d_.maxLExtrap)) continue;
        if (tryBareInteger("LMAXADD",   d_.lMaxAdditional))    continue;
        if (tryBareInteger("PRINT",      d_.printLevel)) continue;

        // --- ASYMPTOPIA as standalone keyword (no '='), value is next token ---
        if (u == "ASYMPTOPIA" && i + 1 < toks.size()) {
            try { d_.asymptopia = std::stod(toks[i+1]); i++; } catch (...) {}
            continue;
        }

        // LABANGLES / LABANGLE — switch keyword: output angles in lab frame
        if (u == "LABANGLES" || u == "LABANGLE") {
            d_.labAngles = true;
            continue;
        }

        // jResidual=0+ at top level
        if (u.rfind("JBIGA=", 0) == 0) {
            std::string jStr = u.substr(6);
            int twoJ = 0;
            if (parseHalfInt(jStr, twoJ)) {
                d_.jResidual = twoJ;
            }
            continue;
        }
        // Bare jResidual= followed by value
        if (u == "JBIGA=" && i + 1 < toks.size()) {
            int twoJ = 0;
            if (parseHalfInt(toks[i+1], twoJ)) {
                d_.jResidual = twoJ; i++;
            }
            continue;
        }
    }
    return true;
}


// ============================================================================
// Helper: stuff a string into INPBUF so nxWord/nxValue/REACTN can read it
// ============================================================================
static void stuffInpbuf(const std::string& line) {
    int len = (int)std::min(line.size(), (size_t)iBufSize);
    for (int i = 1; i <= len; i++) inputBuffer.iBuf[i] = line[i-1];
    for (int i = len + 1; i <= iBufSize; i++) inputBuffer.iBuf[i] = ' ';
    // Find nOch (last non-blank)
    inputBuffer.nOch = len;
    while (inputBuffer.nOch > 1 && inputBuffer.iBuf[inputBuffer.nOch] == ' ') inputBuffer.nOch--;
    inputBuffer.inCh = 1;
}

// Helper: set OM params in FLOAT_common from a ParsedInput::OMParams block.
// Uses the `has<X>` flags (set true by parseKeyvals when the user explicitly
// wrote the keyword) so that `VSI=0.0` correctly clears the prior-channel
// value, distinct from the user simply omitting VSI.
static void setOMparams(const ParsedInput::OMParams& om, Reaction& reaction) {
    if (om.hasV)    reaction.opticalPotentialParams.V    = om.V;
    if (om.hasR0)   reaction.opticalPotentialParams.R0   = om.r0;
    if (om.hasA)    reaction.opticalPotentialParams.A    = om.a;
    if (om.hasVI)   reaction.opticalPotentialParams.vI   = om.vI;
    if (om.hasRI0)  reaction.opticalPotentialParams.rI0  = om.rI0;
    if (om.hasAI)   reaction.opticalPotentialParams.aI   = om.aI;
    if (om.hasVSI)  reaction.opticalPotentialParams.vSi  = om.vSi;
    if (om.hasRSI0) reaction.opticalPotentialParams.rSi0 = om.rSi0;
    if (om.hasASI)  reaction.opticalPotentialParams.aSi  = om.aSi;
    if (om.hasVSO)  reaction.opticalPotentialParams.vSo  = om.vSo;
    if (om.hasRSO0) reaction.opticalPotentialParams.rSo0 = om.rSo0;
    if (om.hasASO)  reaction.opticalPotentialParams.aSo  = om.aSo;
    if (om.hasVSOI) reaction.opticalPotentialParams.vSoi = om.vSoi;
    if (om.hasRSOI0)reaction.opticalPotentialParams.rSoi0= om.rSoi0;
    if (om.hasASOI) reaction.opticalPotentialParams.aSoi = om.aSoi;
    if (om.hasRC0)  reaction.opticalPotentialParams.rC0  = om.rC0;
}

// Helper: the 3-statement head shared by the CONTRL incoming/outgoing channel-setup
// paths — flag the next block, push OM params (when present), then trigger setChannel.
static void applyOMandSetChannel(int channel, const ParsedInput::OMParams& om,
                                 Reaction& reaction, int& returnCode) {
    reaction.flags.hasNextBlock = channel;
    if (om.set) setOMparams(om, reaction);
    reaction.setChannel(channel, returnCode);
}

// Helper: apply OM params + set scattering channel, run SETPOT and scattering
// waves, then mark iDone. Returns false (caller should return) on setChannel or
// SETPOT failure. tag is " [inelastic]" for the CONTRL inelastic path, "" otherwise.
// Caller keeps its own clearChannel(channel) (provenance comments differ per site).
static bool setupScatteringChannel(int channel, const ParsedInput::OMParams& om,
                                   Reaction& reaction, int& returnCode, const char* tag) {
    applyOMandSetChannel(channel, om, reaction, returnCode);
    if (returnCode == 0) {
        std::fprintf(stderr, "applyToCommons: reaction.setChannel(%d)%s failed\n", channel, tag);
        return false;
    }
    if (!reaction.setupOpticalPotential()) {
        std::fprintf(stderr, "applyToCommons: SETPOT(%d)%s failed\n", channel, tag);
        return false;
    }
    reaction.internalState.waveChannel = channel - 2;   // ch3->1, ch4->2
    reaction.distortedWave.scatteringSolver.setupScatteringWaves(returnCode, 0, reaction);
    reaction.internalState.iDone |= (1 << (channel - 1));
    return true;
}

// Helper: set BS (bound state) potential params in FLOAT_common.
// Uses `has<X>` flags so explicit `=0.0` overrides distinguish from the user
// omitting the keyword (same pattern as setOMparams).
static void setBSparams(const ParsedInput::BSParams& bs, Reaction& reaction) {
    if (bs.hasV)    reaction.opticalPotentialParams.V   = bs.V;
    if (bs.hasR0)   reaction.opticalPotentialParams.R0  = bs.r0;
    if (bs.hasA)    reaction.opticalPotentialParams.A   = bs.a;
    if (bs.hasVSO)  reaction.opticalPotentialParams.vSo = bs.vSo;
    if (bs.hasRSO0) reaction.opticalPotentialParams.rSo0= bs.rSo0;
    if (bs.hasASO)  reaction.opticalPotentialParams.aSo = bs.aSo;
    if (bs.hasRC0)  reaction.opticalPotentialParams.rC0 = bs.rC0;
}

// ============================================================================
// applyInputOverrides() — trivial direct-write field overrides only.
//
// These are the parsed values that map 1:1 onto a Reaction-class field with no
// transformation needed. Extracted from applyToCommons() so that future
// refactors can call this standalone (e.g. once setChannel/SETPOT/BOUND no longer
// require the rest of applyToCommons' DEFALT/REACTN/PARAM chain).
//
// Anything that needs masstable lookups, pool allocations, or potential
// construction does NOT belong here — keep it in applyToCommons().
// ============================================================================
void InputParser::applyInputOverrides(Reaction& reaction) {
    const auto& d = d_;

    if (d.r0Target)          reaction.flags.r0Type           = 1;
    if (d.labAngles)        reaction.flags.outputInLab           = 1;
    if (d.spFactor >= 0.0)     reaction.spec.specFactorTgt  = d.spFactor;
    if (d.spFactorProj >= 0.0)     reaction.spec.specFactorProj = d.spFactorProj;
    if (d.lMax >= 0)         reaction.angMom.lMax             = d.lMax;
    if (d.lMin >= 0)         reaction.angMom.lMin             = d.lMin;
    if (d.lStep > 0)         reaction.integrationGrid.lStep            = d.lStep;
    if (d.maxLExtrap >= 0)   reaction.integrationGrid.maxLExtrap           = d.maxLExtrap;
    if (d.lMaxAdditional > 0)       reaction.angMom.lMaxAdditional    = d.lMaxAdditional;
    if (d.asymptopia > 0)    reaction.integrationGrid.asymptopia   = d.asymptopia;
    if (d.stepsPer > 0)      reaction.integrationGrid.stepsPerUnit = d.stepsPer;
    if (d.printLevel > 0)   reaction.flags.printLevel           = d.printLevel;

    // Only override angle params if explicitly specified in the input
    if (d.hasAngleMin)   reaction.rxn.angleMin  = d.angleMin;
    if (d.hasAngleMax)   reaction.rxn.angleMax  = d.angleMax;
    if (d.hasAngleStep) reaction.rxn.angleStep = d.angleStep;

    // jResidual (= js[3] = J of target nucleus = particle 3 in PTOLEMY indexing)
    // jResidual=0+ means target ground state is J=0.
    if (d.jResidual >= 0 && !d.isElastic)
        reaction.angMom.js[3] = d.jResidual;
}

// ============================================================================
// applyToCommons() — bypass CONTRL entirely
//
// Supports: transfer (d,p / p,d / 3He,a), elastic (CHANNEL: format),
//           inelastic (p,p' / d,d' / a,a').
//
// Flow:
//   DEFALT → (CHANEL or REACTN) → PARAM → global params → ISTART
//   → per-channel SETCHN+SETPOT+BOUND/WAVSET+CLRCHN
// ============================================================================
void InputParser::applyToCommons(Reaction& reaction) {
    const auto& d = d_;

    // -----------------------------------------------------------------------
    // 1. DEFALT
    // -----------------------------------------------------------------------
    // KFRMTP/IBRNSB/ISAVSM batch + IBSPAS + IASYMP/ICHECK/MASTYP/ITSO + IFIT
    // NUMRAN/NBLKSZ/NFIROF) + the 5 write-only IPARM1..IPARM5 trailers
    // → 31 → 30 (NBACK + MAXITR + NAITKN + LOOKST + NVPOLY + MAPSUM +
    // (FITMUL/FITACC/DERIVS/FITRAT) removed + STEP1R/STEP1I folded inline
    // VTR/VTRI/VTL/VTLI/VTP/VTPI (tensor central depths, all 0, no input)
    // RTR/RTRI/RTL/RTLI/RTP/RTPI (tensor radii, dead after central depths
    // RTR0/RTRI0/RTL0/RTLI0/RTP0/RTPI0 (tensor radii-at-zero, same cascade)
    // ATR/ATRI/ATL/ATLI/ATP/ATPI (tensor diffuseness, same cascade) dropped
    reaction.applyDefaults();
    newCard();

    // CLI flag propagation: --fixedLS switches spin-orbit coupling to
    // physics-standard <L*S> form (default is Cleopatra-faithful sigma*L).
    reaction.flags.fixedLS = d.cliFixedLS;

    // -----------------------------------------------------------------------
    // 2a. Elastic: use CHANEL to parse "projectile + target"
    // 2b. Transfer/inelastic: use REACTN
    // -----------------------------------------------------------------------
    if (d.isElastic) {
        // Build "PROJ + TARGET" string for CHANEL
        std::string chstr = d.channelProjectile + " + " + d.channelTarget;
        stuffInpbuf(chstr);
        if (!reaction.setupChannel()) {
            std::fprintf(stderr, "applyToCommons: CHANEL failed on '%s'\n", chstr.c_str());
        }
    } else {
        std::string rLine = d.reactionLine;
        std::string ru = InputParser::toUpper(rLine);
        size_t ePos = ru.find("ELAB=");
        if (ePos == std::string::npos) ePos = ru.find("ECM=");
        if (ePos != std::string::npos) rLine = InputParser::trim(rLine.substr(0, ePos));
        stuffInpbuf(rLine);
        if (!reaction.parseReactionString()) {
            std::fprintf(stderr, "applyToCommons: REACTN failed on '%s'\n", rLine.c_str());
        }
    }

    // -----------------------------------------------------------------------
    // 3. eLab or eCm. set_channels (formerly setChannel) at L450 reads whichever
    // is != undefValue; if both are set, eCm wins. So write only the field the
    // input file specified — leave the other as undefValue (from DEFALT).
    // -----------------------------------------------------------------------
    if (d.eCm >= 0.0) {
        reaction.energies.eCm = d.eCm;
    } else {
        reaction.energies.eLab = d.eLab;
    }

    // -----------------------------------------------------------------------
    // 4. PARAMETERSET
    // -----------------------------------------------------------------------
    if (!d.parameterSet.empty()) {
        stuffInpbuf(d.parameterSet);
        if (!reaction.applyParameterSet()) {
            std::fprintf(stderr, "applyToCommons: PARAM('%s') failed\n",
                         d.parameterSet.c_str());
        }
    }

    // -----------------------------------------------------------------------
    // 5. Global parameters (override PARAM defaults)
    // -----------------------------------------------------------------------
    applyInputOverrides(reaction);

    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // 7. isElastic flag for elastic
    // -----------------------------------------------------------------------
    if (d.isElastic) {
        reaction.flags.isElastic = 1;
    }

    int returnCode = 0;

    // =======================================================================
    // ELASTIC path: one channel (reaction.setChannel(6)), then WAVSET
    // =======================================================================
    if (d.isElastic) {
        reaction.flags.hasNextBlock   = 6;
        reaction.flags.problemType = 6;   // standalone elastic
        reaction.internalState.iDone  = 0;

        if (d.incoming.set) setOMparams(d.incoming, reaction);

        reaction.setChannel(6, returnCode);
        if (returnCode == 0) {
            std::fprintf(stderr, "applyToCommons: reaction.setChannel(6) failed\n");
            return;
        }
        if (!reaction.setupOpticalPotential()) {
            std::fprintf(stderr, "applyToCommons: SETPOT(6) failed\n");
            return;
        }
        reaction.internalState.waveChannel = 1;
        reaction.distortedWave.scatteringSolver.setupScatteringWaves(returnCode, 0, reaction);
        // NOTE: Do NOT call CLRCHN for elastic. In ptolemy, CONTRL sets IGOTO=13
        // and returns directly to main which calls WAVEF. CLRCHN is not called
        // until after WAVEF completes. WAVEF needs R, E, etc. to still be valid.
        reaction.flags.hasNextBlock = 0;
        return;  // caller uses WAVEF for DCS
    }

    // =======================================================================
    // INELASTIC path: two channels (3=INCOMING, 4=OUTGOING), no BS channels
    // CONTRL flow: REACTN → PARAM → DEFINE(BELX) → global → reaction.setChannel(3) → SETPOT → WAVSET
    //              → CLRCHN(3) → reaction.setChannel(4) → SETPOT → WAVSET → CLRCHN(4)
    //              → then computes INGRST→ANGSET→COULST→INRDIN→LINTRP→XSECTN
    // =======================================================================
    if (d.belx > 0.0) {
        reaction.flags.problemType = 22;   // inelastic
        reaction.internalState.iDone  = 0;

        // Store BELX deformation parameter in allocator pool (CONTRL does DEFINE).
        // DEFINE reads values from INPBUF via nxValue.
        {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%.10g", d.belx);
            stuffInpbuf(buf);
            char8 belxName("BELX    ");
            InputParser::defineArray(belxName, reaction);
        }

        // BETACOUL (Coulomb deformation) — DEFINE-style allocation, like BELX.
        // probe_print.cpp looks it up via NAMLOC("BETACOUL").
        if (d.betaCoul > 0.0) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%.10g", d.betaCoul);
            stuffInpbuf(buf);
            char8 bcName("BETACOUL");
            InputParser::defineArray(bcName, reaction);
        }

        // ---- INCOMING ----
        // In CONTRL, reaction.setChannel(3) is triggered at the ";" that closes the INCOMING
        // block, by which time all OM params have been read into reaction.rxn.
        // Set OM params from d.incoming BEFORE reaction.setChannel(3).
        if (!setupScatteringChannel(3, d.incoming, reaction, returnCode, " [inelastic]")) return;
        reaction.clearChannel(3);
        // After CLRCHN(3) for inelastic, V, R0, A etc. survive
        // (CLRCHN goes via L400 path when iDone==(1<<(3-1))).

        // ---- OUTGOING ----
        // "OUTGOING" triggers reaction.setChannel(4) with OM params still in FLOAT_common
        // from the INCOMING block (they survive CLRCHN(3)).
        // If OUTGOING block has explicit params, they override.
        if (!setupScatteringChannel(4, d.outgoing, reaction, returnCode, " [inelastic]")) return;

        // Match legacy CONTRL: nothing saved/restored around CLRCHN(4).
        // CLRCHN clears E/aM/massProj/massTgt to undefValue (and eLab if
        // wasSet[16] is set), but its L800 block (stripPickup==0 &&
        // iDone==12) restores all the per-channel-1 values the prologue
        // needs. Saving the OUTGOING-channel values across CLRCHN(4) and
        // re-applying them after pollutes commons relative to legacy.
        reaction.clearChannel(4);

        reaction.flags.hasNextBlock = 0;
        return;  // caller uses elastic→PRBPRT→GETSCT→INGRST→ANGSET→COULST→INRDIN→LINTRP→XSECTN
    }

    // =======================================================================
    // TRANSFER path (d,p / p,d / 3He,a / etc.)
    // Channels: 1=PROJECTILE BS, 2=TARGET BS, 3=INCOMING, 4=OUTGOING
    // =======================================================================
    reaction.flags.problemType = 20;   // assume DWBA transfer
    reaction.internalState.iDone  = 0;

    // ---- Channel 1: PROJECTILE ----
    reaction.flags.hasNextBlock = 1;
    if (d.projectileBs.set) {
        reaction.internalState.lSpecs[1] = d.projectileBs.l;
        reaction.internalState.nodesP[1] = d.projectileBs.nodes;
        reaction.angMom.L         = d.projectileBs.l;
        reaction.angMom.nNodes     = d.projectileBs.nodes;
        if (d.projectileBs.j >= 0)
            reaction.angMom.jProj = (int)(2.0 * d.projectileBs.j);
        if (!d.projectileBs.wavefunction.empty()) {
            stuffInpbuf(d.projectileBs.wavefunction);
            char8 wfKey("WAVEFUNC");
            reaction.loadLinkule(wfKey);
        }
        setBSparams(d.projectileBs, reaction);
    }
    reaction.setChannel(1, returnCode);
    if (returnCode == 0) {
        std::fprintf(stderr, "applyToCommons: reaction.setChannel(1) failed\n");
        return;
    }
    // After reaction.setChannel(1) we know stripPickup
    reaction.flags.problemType = (reaction.internalState.stripPickup == 0) ? 22 : 20;

    if (!reaction.setupOpticalPotential()) {
        std::fprintf(stderr, "applyToCommons: SETPOT(1) failed\n");
        return;
    }
    reaction.internalState.boundChannel = 1;
    reaction.boundState.solve(returnCode, reaction);
    reaction.internalState.iDone |= (1 << 0);
    reaction.clearChannel(1);

    // ---- Channel 2: TARGET ----
    // For reaction.setChannel(2) the boundIndex index depends on stripPickup:
    //   stripping (>0): boundIndex=4 → lSpecs[4], nodesP[4]; reaction.angMom.J = jp of transferred nucleon
    //   pickup    (<0): boundIndex=3 → lSpecs[3], nodesP[3]; reaction.angMom.J left as notDefSentinel (SETCHN reads js[3]=J(target))
    //                             reaction.angMom.jProj = jp of transferred nucleon
    reaction.flags.hasNextBlock = 2;
    if (d.targetBs.set) {
        reaction.angMom.L     = d.targetBs.l;
        reaction.angMom.nNodes = d.targetBs.nodes;

        if (reaction.internalState.stripPickup > 0) {
            // Stripping: boundIndex=4, IT=3
            reaction.internalState.lSpecs[4] = d.targetBs.l;
            reaction.internalState.nodesP[4] = d.targetBs.nodes;
            // reaction.angMom.J = jp of the transferred nucleon (= J of the residual's state)
            if (d.targetBs.j >= 0)
                reaction.angMom.J = (int)(2.0 * d.targetBs.j);
        } else {
            // Pickup: boundIndex=3, IT=4
            // reaction.angMom.J is left notDefSentinel → SETCHN uses js[3] = J(target nucleus = 48Ca) = 0
            // jProj = j of the transferred nucleon (the angular momentum being picked up)
            reaction.internalState.lSpecs[3] = d.targetBs.l;
            reaction.internalState.nodesP[3] = d.targetBs.nodes;
            if (d.targetBs.j >= 0)
                reaction.angMom.jProj = (int)(2.0 * d.targetBs.j);
            // Leave reaction.angMom.J at notDefSentinel so SETCHN uses js[boundIndex=3] = J(composite/target)
        }
        if (!d.targetBs.wavefunction.empty()) {
            stuffInpbuf(d.targetBs.wavefunction);
            char8 wfKey("WAVEFUNC");
            reaction.loadLinkule(wfKey);
        }
        setBSparams(d.targetBs, reaction);
    }
    reaction.setChannel(2, returnCode);
    if (returnCode == 0) {
        std::fprintf(stderr, "applyToCommons: reaction.setChannel(2) failed\n");
        return;
    }
    if (!reaction.setupOpticalPotential()) {
        std::fprintf(stderr, "applyToCommons: SETPOT(2) failed\n");
        return;
    }
    reaction.internalState.boundChannel = 2;
    reaction.boundState.solve(returnCode, reaction);
    reaction.internalState.iDone |= (1 << 1);
    reaction.clearChannel(2);

    // ---- Channel 3: INCOMING ----
    if (!setupScatteringChannel(3, d.incoming, reaction, returnCode, "")) return;
    reaction.clearChannel(3);

    // ---- Channel 4: OUTGOING ----
    if (!setupScatteringChannel(4, d.outgoing, reaction, returnCode, "")) return;
    reaction.clearChannel(4);

    reaction.flags.hasNextBlock = 0;
    if ((reaction.internalState.iDone & 15) != 15) {
        std::fprintf(stderr, "applyToCommons: WARNING: IDONE=%d, expected 15\n",
                     reaction.internalState.iDone);
    }
}

// ============================================================================
// dispatchCalculation() — pick the right Reaction class and run it.
//
// Selection rules (matched the legacy ptolemy_main dispatch):
//   isElastic    → Elastic            (single-channel ELASTIC SCATTERING)
//   BELX > 0.0   → InelasticReaction  (inelastic DWBA)
//   default      → Transfer           (DWBA pickup/transfer)
// ============================================================================
bool InputParser::dispatchCalculation(Reaction& reaction) const {
    if (d_.isElastic) return Elastic(reaction).calculate();
    if (d_.belx > 0.0) return InelasticReaction(reaction).calculate();
    return Transfer(reaction).calculate();
}


