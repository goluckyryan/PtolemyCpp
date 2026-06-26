// DWBAGrid.cpp — DWBAGrid: integration-grid array storage for gridSet / angleSet.

#include "DWBAGrid.h"
#include "Reaction.h"


// ---------------------------------------------------------------------------
// allocateRioEx(size)
// Resize rioEx to size+1 elements (0-based; trailing element spare).
// Updates reaction.gridData.rioExPointer (class-owned).
// ---------------------------------------------------------------------------
void DWBAGrid::allocateRioEx(int size, Reaction& reaction)
{
    if (size <= 0) return;
    rioEx.assign(size + 1, 0.0);
    reaction.gridData.rioExPointer = rioEx.data();  // 0-based
}

// ---------------------------------------------------------------------------
// allocateSmhpts(nPhiSum)
// Sets reaction.gridData.smhptsPointer = smhpts.data()-1.
// ---------------------------------------------------------------------------
void DWBAGrid::allocateSmhpts(int nPhiSum, Reaction& reaction)
{
    smhpts.assign(nPhiSum + 1, 0.0);
    reaction.gridData.smhptsPointer = smhpts.data() - 1;
}

// ---------------------------------------------------------------------------
// allocateSmhwk(nPhiSum)
// Sets reaction.gridData.smhwkPointer = smhwk.data()-1.
// ---------------------------------------------------------------------------
void DWBAGrid::allocateSmhwk(int nPhiSum, Reaction& reaction)
{
    smhwk.assign(3 * nPhiSum + 1, 0.0);
    reaction.gridData.smhwkPointer = smhwk.data() - 1;
}

// ---------------------------------------------------------------------------
// allocateSmipts(nPhiSumI)
// Sets reaction.gridData.smiptsPointer = smipts.data()-1.
// ---------------------------------------------------------------------------
void DWBAGrid::allocateSmipts(int nPhiSumI, Reaction& reaction)
{
    smipts.assign(nPhiSumI + 1, 0.0);
    reaction.gridData.smiptsPointer = smipts.data() - 1;
}

// ---------------------------------------------------------------------------
// allocateSmivl(nPhiSumI)
// Sets reaction.gridData.smivlPointer = smivl.data()-1.
// ---------------------------------------------------------------------------
void DWBAGrid::allocateSmivl(int nPhiSumI, Reaction& reaction)
{
    smivl.assign(nPhiSumI + 1, 0.0);
    reaction.gridData.smivlPointer = smivl.data() - 1;
}

// ---------------------------------------------------------------------------
// allocateSmhvl(nPhiSumHCount)
// Sets reaction.gridData.smhvlPointer = smhvl.data()-1.
// ---------------------------------------------------------------------------
void DWBAGrid::allocateSmhvl(int nPhiSumHCount, Reaction& reaction)
{
    smhvl.assign(nPhiSumHCount + 1, 0.0);
    reaction.gridData.smhvlPointer = smhvl.data() - 1;
}

// ---------------------------------------------------------------------------
// allocateHint(nMloLx)
// hintPointer = hint.data(), 0-based (accessed [hIndex-1]); the +1-element
// pad keeps a spare slot at the tail.
// ---------------------------------------------------------------------------
void DWBAGrid::allocateHint(int nMloLx, Reaction& reaction)
{
    hint.assign(nMloLx + 1, 0.0);
    reaction.gridData.hintPointer = hint.data();
}

// ---------------------------------------------------------------------------
// allocateHabs(nMloLx)
// habsPointer = habs.data(), 0-based (accessed [hIndex-1]); the +1-element
// pad keeps a spare slot at the tail.
// ---------------------------------------------------------------------------
void DWBAGrid::allocateHabs(int nMloLx, Reaction& reaction)
{
    habs.assign(nMloLx + 1, 0.0);
    reaction.gridData.habsPointer = habs.data();
}

// ---------------------------------------------------------------------------
// allocateAbs1(size)
// Sets abs1Pointer = abs1.data() (0-based).
// ---------------------------------------------------------------------------
void DWBAGrid::allocateAbs1(int size, Reaction& reaction)
{
    abs1.assign(size + 1, 0.0);
    reaction.gridData.abs1Pointer = abs1.data();
}

// ---------------------------------------------------------------------------
// allocateIiindx(n)
// Sets iiindxPointer=iiindx.data() (0-based int*).
// ---------------------------------------------------------------------------
void DWBAGrid::allocateIiindx(int n, Reaction& reaction)
{
    iiindx.assign(2 * n + 2, 0);   // FACFR4=2: 2*n+2 ints
    reaction.gridData.iiindxPointer = iiindx.data();
}

// ---------------------------------------------------------------------------
// allocateDw(n)
// Sets dwPointer=dw_.data() (0-based double*).
// ---------------------------------------------------------------------------
void DWBAGrid::allocateDw(int n, Reaction& reaction)
{
    dw_.assign(n + 1, 0.0);  // 0-based: [0..n-1] valid (+1 spare element kept)
    reaction.gridData.dwPointer = dw_.data();  // 0-based double*
}

void DWBAGrid::allocateDwi(int n, Reaction& reaction)
{
    dwi.assign(2 * n + 2, 0);
    reaction.gridData.dwiPointer = dwi.data();
}
void DWBAGrid::allocateIwfii(int n, Reaction& reaction)
{
    iwfii.assign(2 * n + 2, 0);
    reaction.gridData.iDwfiPointer = iwfii.data();
}
void DWBAGrid::allocateIwfio(int n, Reaction& reaction)
{
    iwfio.assign(2 * n + 2, 0);
    reaction.gridData.iDwfoPointer = iwfio.data();
}
void DWBAGrid::allocateNlam(int size)
{
    nlam.assign(size, 0);
}
void DWBAGrid::allocateXlam(int size)
{
    xlam.assign(size, 0.0);
}
void DWBAGrid::resizeIwfio(int n, Reaction& reaction)
{
    // Shrink (replaces IREDEF(4*NWFO/FACFR4, IWFIO) in angleSet).
    // n = new logical size in doubles; keep at least 2*n+2 ints.
    int needed = 2 * n + 2;
    if ((int)iwfio.size() > needed) iwfio.resize(needed);
    reaction.gridData.iDwfoPointer = iwfio.data();  // re-sync pointer after potential resize
}

// ---------------------------------------------------------------------------
// allocateLir/Lii/Lor/Loi(nRiRoH)
// Sets lirPointer/liiPointer/lorPointer/loiPointer = vector.data().
// ---------------------------------------------------------------------------
void DWBAGrid::allocateLir(int nRiRoH, Reaction& reaction)
{
    lir.assign(2 * nRiRoH + 2, 0.0f);
    reaction.gridData.lirPointer = lir.data();
}
void DWBAGrid::allocateLii(int nRiRoH, Reaction& reaction)
{
    lii.assign(2 * nRiRoH + 2, 0.0f);
    reaction.gridData.liiPointer = lii.data();
}
void DWBAGrid::allocateLor(int nRiRoH, Reaction& reaction)
{
    lor_.assign(2 * nRiRoH + 2, 0.0f);
    reaction.gridData.lorPointer = lor_.data();
}
void DWBAGrid::allocateLoi(int nRiRoH, Reaction& reaction)
{
    loi.assign(2 * nRiRoH + 2, 0.0f);
    reaction.gridData.loiPointer = loi.data();
}

// ---------------------------------------------------------------------------
// allocateCosin(cosinSize)
// Sets reaction.gridData.cosinPointer = cosin.data() (0-based).
// ---------------------------------------------------------------------------
void DWBAGrid::allocateCosin(int cosinSize, Reaction& reaction)
{
    cosin.assign(cosinSize + 1, 0.0);
    reaction.gridData.cosinPointer = cosin.data();
}

// ---------------------------------------------------------------------------
// allocateLiloR(size)
// Sets reaction.inelastic.liloRPointer = liloR.data() (0-based).
// ---------------------------------------------------------------------------
void DWBAGrid::allocateLiloR(int size, Reaction& reaction)
{
    liloR.assign(size + 1, 0.0);
    reaction.inelastic.liloRPointer = liloR.data();
}

// ---------------------------------------------------------------------------
// allocateLiloI(size)
// Sets reaction.inelastic.liloIPointer = liloI.data() (0-based).
// ---------------------------------------------------------------------------
void DWBAGrid::allocateLiloI(int size, Reaction& reaction)
{
    liloI.assign(size + 1, 0.0);
    reaction.inelastic.liloIPointer = liloI.data();
}


// ---------------------------------------------------------------------------
// allocatePhiArrays(jSize)
// jSize = size in doubles.
// Each array gets 2*jSize+2 floats (0-based; valid indices 0..2*jSize).
// Sets phiTPointer = phiT.data() (0-based float*).
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// allocateRiRoWio(n4rio)
// Each array: 2*n4rio+2 floats (0-based).
// ---------------------------------------------------------------------------
void DWBAGrid::allocateRiRoWio(int n4rio, Reaction& reaction)
{
    int nFloat = 2 * n4rio + 2;
    ri_.assign(nFloat, 0.0f);
    ro_.assign(nFloat, 0.0f);
    wio_.assign(nFloat, 0.0f);
    reaction.gridData.riPointer  = ri_.data();
    reaction.gridData.roPointer  = ro_.data();
    reaction.gridData.wioPointer = wio_.data();
}

void DWBAGrid::allocatePhiArrays(int jSize, Reaction& reaction)
{
    int nFloat = 2 * jSize + 2;
    phiT.assign(nFloat, 0.0f);
    phiP.assign(nFloat, 0.0f);
    phi_.assign(nFloat, 0.0f);
    trapWeight_.assign(nFloat, 0.0f);
    reaction.gridData.phiTPointer  = phiT.data();
    reaction.gridData.phiPPointer  = phiP.data();
    reaction.gridData.phiPointer   = phi_.data();
    reaction.gridData.trapWeightPointer = trapWeight_.data();
}
