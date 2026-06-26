#pragma once
// angular_coupling.h — A12 angular-momentum coupling (angular_coupling.cpp).
class Reaction;
void angularCoupling12(int li, int lxMin, int lxMax, int lMin, int lMax, double* xLam,
         int* nLam, double* a12Vl, double* msval,
         int* jA12S, int& jA12M, int& jA12N, int& jA12An, int* iIndex,
         int& hCount, int& loMinMin, int& loMaxMax, double* dInts, double* outTemp,
         double* xLoTemp, int* lxTemp, int la12Vl, int printLevel,
         int* indxDw, int& dwCount, int& numIi, int* iDwfi, int* iDwfo,
         Reaction& reaction);
