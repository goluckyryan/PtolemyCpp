#pragma once
// InputParser.h — InputParser class for Ptolemy .in files.
//
// The parser reads the standard Ptolemy .in format:
//   REACTION: A(d,p)B(jProj Ex) ELAB=E
//   PARAMETERSET dpsb [r0Target]
//   [lMax=N lMin=N lStep=N maxLExtrap=N asymptopia=X]
//   PROJECTILE
//     wavefunction av18
//     r0=X a=X l=N rc0=X
//   ;
//   TARGET
//     nodes=N l=N jp=N/2
//     v=X r0=X a=X vso=X rso0=X aso=X rc0=X
//   ;
//   INCOMING
//     v=X r0=X a=X vi=X ri0=X ai=X vsi=X rsi0=X asi=X vso=X rso0=X aso=X rc0=X
//   ;
//   OUTGOING
//     v=X r0=X a=X vi=X ri0=X ai=X vsi=X rsi0=X asi=X vso=X rso0=X aso=X rc0=X
//   ;
//   anglemin=X anglemax=X anglestep=X
//   ;
//   end

#include <string>
#include <vector>
#include <iosfwd>

#include "ptolemy_types.h"  // char8

class Reaction;

// ============================================================================
// InputBuffer — single card image used by the Fortran-port input scanner
// (mScan/NXINT/nxWord/qvScan/channelScan/...). Owned by InputParser as a private
// static member; accessed by free-function scanners via InputParser::buffer().
// ============================================================================
constexpr int iBufSize = 200;

struct InputBuffer {
    int  nOch;                  // number of characters in the current word
    int  inCh;                  // 1-based scan position in iBuf
    char iBuf[iBufSize + 1];    // 1-based [1..200]: current input line
};

// ============================================================================
// ParsedInput — all data extracted from a Ptolemy .in file
// ============================================================================
struct ParsedInput {
    // Reaction string (raw, for REACTN-style processing)
    std::string reactionLine;    // e.g. "48CA(D,P)49CA(3/2- 0.000)"

    double eLab     = 0.0;       // lab energy (MeV) — used when input gave eLab
    double eCm      = -1.0;      // CM energy (MeV) — set when input gave eCm, else -1

    // Parameterset
    std::string parameterSet;    // e.g. "DPSB"
    bool r0Target   = false;     // R0TARGET modifier

    // Global grid/L params (overrides from input)
    int    lMax        = -1;     // -1 = not set
    int    lMin        = -1;     // -1 = not set (matches lMax sentinel)
    int    lStep       = -1;     // -1 = not set
    int    maxLExtrap  = -1;     // -1 = not set (PARAMETERSET wins); 0 is a valid explicit override
    double asymptopia  = -1.0;   // -1 = not set
    double stepsPer    = -1.0;   // STEPSPER (steps per wavelength); -1 = not set
    bool   labAngles  = false;  // LABANGLES switch — output in lab frame
    double betaCoul    = -1.0;   // BETACOUL=N (Coulomb deformation); -1 = not set
    double spFactor      = -1.0;   // SPFACT=N (target spectroscopic factor); -1 = not set
    double spFactorProj      = -1.0;   // SPFACP=N (projectile spectroscopic factor); -1 = not set

    // Bound-state block (projectile side). Each potential-parameter field has a
    // paired `has<X>` bool marking whether the user explicitly wrote the keyword
    // — same pattern as OMParams to distinguish `VSO=0.0` (explicit zero) from
    // "not present in input" (leave the prior value alone).
    struct BSParams {
        int    nodes   = 0;
        int    l       = 0;
        double j       = -1.0;   // half-integer allowed; -1 = not set
        double V       = 0.0;    // real depth (MeV); 0 = use eigenvalue search
        double r0      = 1.0;    // radius parameter (fm)
        double a       = 0.5;    // diffuseness (fm)
        double vSo     = 0.0;
        double rSo0    = 0.0;
        double aSo     = 0.0;
        double rC0     = 0.0;    // Coulomb radius param
        double E       = 0.0;    // binding-energy override (MeV); < 0 for bound
        bool hasV    = false;  bool hasR0   = false;  bool hasA    = false;
        bool hasVSO  = false;  bool hasRSO0 = false;  bool hasASO  = false;
        bool hasRC0  = false;  bool hasE    = false;
        std::string wavefunction; // "av18", "rcwfn", "phiffer", or ""
        bool   set     = false;  // true if the block was present in input
    };
    BSParams projectileBs;      // PROJECTILE block
    BSParams targetBs;          // TARGET block
    int      jResidual   = 0;        // 2*J for residual (read as "jResidual=0+")

    // Optical potential block. Each field has a paired `has<X>` bool that
    // marks whether the user explicitly wrote the keyword in the input deck.
    // The bool is true even when the explicit value is 0.0 — this distinguishes
    // `VSI=0.0` (user asked for zero surface absorption) from the user leaving
    // VSI unspecified (the prior-channel value should persist).
    struct OMParams {
        double V    = 0.0;   double r0   = 0.0;   double a    = 0.0;
        double vI   = 0.0;   double rI0  = 0.0;   double aI   = 0.0;
        double vSi  = 0.0;   double rSi0 = 0.0;   double aSi  = 0.0;
        double vSo  = 0.0;   double rSo0 = 0.0;   double aSo  = 0.0;
        double vSoi = 0.0;   double rSoi0= 0.0;   double aSoi = 0.0;
        double rC0  = 0.0;
        bool hasV    = false;  bool hasR0   = false;  bool hasA    = false;
        bool hasVI   = false;  bool hasRI0  = false;  bool hasAI   = false;
        bool hasVSI  = false;  bool hasRSI0 = false;  bool hasASI  = false;
        bool hasVSO  = false;  bool hasRSO0 = false;  bool hasASO  = false;
        bool hasVSOI = false;  bool hasRSOI0= false;  bool hasASOI = false;
        bool hasRC0  = false;
        bool   set   = false; // true if the block was present in input
    };
    OMParams incoming;
    OMParams outgoing;

    // Angle output (use hasAngle* to know if explicitly set)
    double angleMin  =   0.0;
    double angleMax  = 180.0;
    double angleStep =   2.0;
    bool hasAngleMin  = false;
    bool hasAngleMax  = false;
    bool hasAngleStep = false;

    // Inelastic (not used for d,p but keep for completeness)
    double belx  = 0.0;

    // Reaction type flags. Inelastic-vs-transfer dispatch happens in
    // dispatchCalculation by inspecting parameterSet names, not via parsed flags.
    bool isElastic     = false;  // ELASTIC SCATTERING

    // CHANNEL: keyword target/projectile (elastic format)
    std::string channelTarget;       // e.g. "208PB"
    std::string channelProjectile;   // e.g. "P"

    // PRINT level
    int printLevel = 0;

    // LMAXADD
    int lMaxAdditional = 0;

    // CLI flag (from --fixedLS) — use physics-standard <L*S> coupling.
    // Propagated to reaction.flags.fixedLS in applyToCommons().
    bool cliFixedLS = false;
};

// ============================================================================
// InputParser — parses a Ptolemy .in file and applies it to the Reaction state
// ============================================================================
class InputParser {
public:
    // Parse the file at `path`. Returns true on success.
    bool parse(const std::string& path);

    // Parse from an arbitrary input stream (eg std::cin). Returns true
    // on success. Used by the main binary's stdin path so test_sakura.sh
    // can keep its `./ptolemy < input.in` invocation shape.
    bool parse(std::istream& in);

    // Parse from argc/argv exactly like the ptolemy CLI: argv[1] = file path,
    // or read from std::cin when no positional arg is given. Emits a stderr
    // diagnostic on parse failure. Returns true on success.
    bool parseFromArgs(int argc, char** argv);

    // After applyToCommons(): dispatch to the right Reaction calculation
    // (Elastic / Inelastic / Transfer) based on parsed flags. Returns the
    // bool from the chosen calculate() method.
    bool dispatchCalculation(Reaction& reaction) const;

    // Apply parsed values to the reaction's COMMON-derived structs.
    // Must be called after parse() succeeds.
    void applyToCommons(Reaction& reaction);

    // Read a sequence of doubles from the input buffer and store them in
    // reaction.named under the 8-char key 'name'. Used during DEFINE BELX / BC
    // input cards.
    static void defineArray(char8 name, Reaction& reaction);

    // Access the parsed data
    const ParsedInput& data() const { return d_; }

    // Access the shared input scanner buffer. The storage is a private static
    // member of InputParser; this accessor is the only way free-function
    // scanners (mScan/nxWord/qvScan/channelScan/...) reach it.
    static InputBuffer& buffer() { return buffer_; }

private:
    ParsedInput d_;

    // Private static storage for the input scanner buffer. Accessed only via
    // the public static buffer() accessor above. inline-static means single
    // definition across all TUs and reliable static-init.
    inline static InputBuffer buffer_ = {};

    // Apply ONLY the trivial direct-write overrides (lMax, lMin, angles, ...).
    // Called from applyToCommons() as the "step 5" overrides block.
    void applyInputOverrides(Reaction& reaction);

    // DWBA input support (Phase C). hasDwbaKeyword() is true when the first
    // content line is the literal "DWBA". parseDwba() expands a human-readable
    // DWBA reaction description to a Ptolemy deck (DwbaInputExpander) and parses
    // it. Detection/dispatch live in parseFromArgs.
    static bool hasDwbaKeyword(const std::string& content);
    bool parseDwba(const std::string& content);

    // Parsing helpers
    bool parseLine(const std::string& raw);
    bool parseReactionLine(const std::string& rest);
    bool parseParametersetLine(const std::string& rest);
    bool parseKeyvals(const std::string& line, ParsedInput::BSParams& bs);
    bool parseKeyvals(const std::string& line, ParsedInput::OMParams& om);
    bool parseGlobalKeyvals(const std::string& line);

    // State machine for section tracking
    enum class Section {
        TOPLEVEL, PROJECTILE, TARGET, INCOMING, OUTGOING
    };
    Section section_ = Section::TOPLEVEL;
    bool    inBlock_ = false; // true when inside a block (waiting for ;)

    // Token utilities
    static std::string toUpper(const std::string& s);
    static std::string trim(const std::string& s);
    static std::vector<std::string> tokenize(const std::string& s);
    // Join toks[start..] into one string with single-space separators.
    static std::string joinTokensFrom(const std::vector<std::string>& toks, size_t start);
    // Parse "3/2" → 3, "7/2-" → 7 (denominator=2, returns 2*J as integer).
    // Trailing +/- is consumed but discarded.
    static bool parseHalfInt(const std::string& tok, int& twoJ);
    // Match a "KEY=value" token (case-insensitive); on success set valStr to the
    // value substring and return true. Returns false on no-match or empty value.
    static bool matchKeyValue(const std::string& tok, const std::string& key, std::string& valStr);
    // Parse a single key=value token; return true if matched
    static bool getDouble(const std::string& tok, const std::string& key, double& out);
    static bool getInt   (const std::string& tok, const std::string& key, int& out);
    // Check if token is exactly "KEY=" (value in next token)
    static bool isBareKey(const std::string& tok, const std::string& key);

};

// Input-scanner free functions (operate on InputParser::buffer()).
void mScan(int keyEx, int& tokenCode, char8& charValue, double& value, int& intValue,
           char* message, int& messageLength, char& stop);
int  nxWord(char* cvArg);
void newCard();
void qvScan(char8* guy, double* eStars, int* nodeVals, int* lVals,
            int* jVals, int* iParities, int& returnCode, Reaction& reaction);
void channelScan(char8* guy, double* eStars, int* nodeVals, int* lVals, int* jVals,
            int* iParities, int& returnCode, Reaction& reaction);
