// ptolemy_main.cpp — main entry for the `ptolemy` binary.
//
// Reads from argv[1] or from stdin (legacy `./ptolemy < input.in`), parses,
// applies to globals, then dispatches to Elastic / InelasticReaction /
// Transfer via InputParser::dispatchCalculation(). No CONTRL state machine.

#include "InputParser.h"
#include "Reaction.h"

int main(int argc, char** argv) {
    InputParser parser;
    if (!parser.parseFromArgs(argc, argv)) return 1;
    Reaction reaction;
    parser.applyToCommons(reaction);
    return parser.dispatchCalculation(reaction) ? 0 : 1;
}
