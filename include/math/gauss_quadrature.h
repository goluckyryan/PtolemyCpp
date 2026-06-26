#pragma once

// Gauss-Legendre and Gauss-Laguerre quadrature

void gaussL(int n, double* x, double* w);

// CSX/CSW/TSX/TSW out-params dropped — sole caller never read them.
void laguerre(int nn, double* x, double* w, double alpha);
// LAGBC/LGRECR/LGROOT — internal helpers; defined static in MathFunctions.cpp.
