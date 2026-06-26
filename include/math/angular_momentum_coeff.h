#pragma once

// Angular momentum coupling coefficients
// Clebsch-Gordan, 3J, 6J, Racah, 9J symbols

double clebschGordan(int a, int b, int x, int y, int cIn, int zIn);
double threeJ(int a, int b, int cIn, int x, int y, int zIn);
double sixJ(int a, int b, int c, int x, int y, int z);
double racah(int a, int b, int y, int x, int c, int z);
double wig9J(int j1, int j2, int j3, int j4, int j5, int j6,
             int j7, int j8, int j9);
