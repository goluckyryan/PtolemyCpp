#pragma once
// masstable.h — nuclide mass-excess and ground-state lookups
// (masstable.cpp, parameters.cpp).
double excess(int kZ, int kA, int& notTabulated);
void   azCode(char* symIn, int& iz, int& ia, int& returnCode);
void   groundStateInfo(int kZ, int kA, int& jVal, int& parity, int& notTabulated);
