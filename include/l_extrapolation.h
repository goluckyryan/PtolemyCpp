#pragma once
// l_extrapolation.h — L-extrapolation helpers for LINTRP (l_extrapolation.cpp).
double lxTrpM(int excitationType, double barA, double b, double barL,
              double lMaxDouble, double weeBoy);
void lxTrp1(int excitationType, int N, int& convergenceCode, int printLevel, double* xVals, double* sVals,
            double& flCrit, double& aVal, double& width,
            double& barL, double& barA, double& b, double& barC, double lMaxDouble,
            double& chiSq, int lx, int lDelta, int jProj, int jT);
void lxTrp2(int excitationType, double barA, double b, double barC, double barL,
            double lMaxDouble, int li, double& size);
