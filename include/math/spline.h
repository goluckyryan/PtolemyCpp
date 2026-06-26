#pragma once

// Cubic spline interpolation

void naturalCubicSpline(int pointCount, double* xPts, double* yIn, double* bCoef, double* cCoef, double* dCoef);
void cubicSplineInterp(int cubicCount, double* xCubes, double* aCoef, double* bCoef, double* cCoef,
            double* dCoef, int nPts, double* xPts, double* yOut);
