#pragma once
// transfer.h — Transfer reaction class (DWBA)
//
// Handles: (d,p), (p,d), (³He,α), (t,p), (⁶Li,d), (⁶Li,α), etc.
//
// Takes: Reaction + BoundState(×2) + DistortedWave(×2)
// Does:  Integration grid setup, DWBA radial integrals,
//        angular coupling, L-interpolation
// Output: dσ/dΩ(θ), analyzing powers, spin observables
//
// Absorbs: GRDSET, INELDC, LINTRP, XSECTN (for transfer),
//          ANGSET, BETCAL, CALANG, AMPCAL, ANAPOW

class Reaction;

class Transfer {
public:
    explicit Transfer(Reaction& reaction);
    // Full transfer calculation:
    //   PRBPRT + GETSCT + GRDSET + INELDC (radial integrals T_L) +
    //   LINTRP + XSECTN + ANAPOW.
    bool calculate();

private:
    Reaction& reaction_;
};
