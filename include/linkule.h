#pragma once
// linkule.h — linkule (linked-potential) dispatcher (linkule.cpp).
#include "ptolemy_types.h"  // char8
class Reaction;
void linkule(int linkuleIndex, char8 alias, int* linkuleInts, int potType, int requestCode,
             int& callStatus, int L, double J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, char* iD,
             Reaction& reaction);
