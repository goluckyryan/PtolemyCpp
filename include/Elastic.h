#pragma once
// Elastic.h — Elastic scattering class
//
// Full elastic scattering calculation:
//   Loop over L → solvePartialWave → S-matrix → ELDCS → dσ/dΩ(θ)
//
// calculate() runs the full L-loop + S-matrix + DCS pipeline.

class Reaction;

class Elastic {
public:
    explicit Elastic(Reaction& reaction);
    bool calculate();

private:
    Reaction& reaction_;
};
