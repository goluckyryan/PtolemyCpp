#pragma once

// Matrix inversion, linear equation solvers, least-squares polynomial fit

// (the whole sub-tree had no external callers).
// sum out-param dropped — sole caller (grid_setup poly-fit stage) never read it.
void lsqPol(double* X, double* Y, double* W, double* residual, int nSub,
            int lSub, double* A, double* B, int mSub);
