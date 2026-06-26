#pragma once

#include "ptolemy_types.h"

// Complex continued fraction interpolation

void cfracInit(int nPts, complex16* xs, complex16* ys);
void cfracEval(complex16* xs, complex16* ys, complex16 x, complex16& y);
