#pragma once
// print_utils.h — formatted-output helpers (parameter_print.cpp).
class Reaction;
void print_G(int w, int d, double val);
void parameterPrint(int i, const char* name, double vPar, double rPar, double aPar,
            double rPar0, Reaction& reaction);
