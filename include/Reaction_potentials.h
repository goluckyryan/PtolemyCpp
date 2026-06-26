#pragma once
// Reaction_potentials.h — sFromI S-matrix renormalization helper (Reaction_potentials.cpp).
class Reaction;
void sFromI(int li, int liIndex, double* sMatR, double* sMatI, int* indxs,
            double* xiReal, double* xiImag, int* iIndex, int numIi,
            int* indxDw, int* iDwfi, int* iDwfo, double* abs1,
            double* aTerm, double factor, int isInfoPrint,
            Reaction& reaction);
