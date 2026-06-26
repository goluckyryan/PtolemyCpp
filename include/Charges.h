#pragma once
// Charges.h — nuclear charges (Z) for the reaction particles.

struct Charges {
    // zProj and zTarget MUST stay adjacent and in this order: channel_setup.cpp
    // takes `int* izPts = &charges.zProj` and reads izPts[0]/izPts[1] as
    // zProj/zTarget.
    int zProj = 0;       // projectile charge
    int zTarget = 0;     // target charge
    int zArray[6] = {};  // 1-based per-particle charges [1..5]
};
