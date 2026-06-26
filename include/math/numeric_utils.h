#pragma once
// numeric_utils.h — numeric helper routines (MathFunctions.cpp).
void cubMap(int mapType, double xLow, double xMidIn, double xHigh, double gamma,
            double* args, double* weights, int nPts);
void epsLon(double* xIn, int pointCount, double* fRet, double& relativeError);
void linLsq(int fitForm, int pointCount, double* xVals, double* sVals, double& cVal,
             double& aVal, double& b, double& chiSq, int debugSwitch);
