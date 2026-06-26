// channel_clear.cpp — CLRCHN: clears channel-specific parameters between channels.

#include "Reaction.h"

void Reaction::clearChannel(int channelIndex)
{
    Reaction& reaction = *this;

    //
    // CLEAR THE DEFINITIONS OF CHANNEL VARIABLES ( P, T, AND POTS )
    //
    // THIS ROUTINE IS CALLED AFTER EACH BOUNDSTATE AND OPTICAL
    // MODEL SCATTERING STATE.  IT SETS CHANNEL VARIABLES THAT ARE
    // NOT LIKELY TO BE SIGNIFICANT FOR LATER CALCULATIONS TO
    // UNDEFINED STATUS.
    //
    // FOR STAND ALONES, THINGS ARE SET TO UNDEFINED ONLY IF SOMETHING
    // ELSE CAN STILL DEFINE THEM (E.G., R IS SET UNDEFINED IF R0,
    // MP AND MT ARE ALL DEFINED).
    //

    // Local references for convenience
    double& A      = reaction.opticalPotentialParams.A;
    double& aI     = reaction.opticalPotentialParams.aI;
    double& aSo    = reaction.opticalPotentialParams.aSo;
    double& aSoi   = reaction.opticalPotentialParams.aSoi;
    double& asymptopia = reaction.integrationGrid.asymptopia;
    double& E      = reaction.energies.E;
    double& eLab   = reaction.energies.eLab;
    double& aM     = reaction.masses.aM;
    double& massProj    = reaction.masses.massProj;
    double& massTgt    = reaction.masses.massTgt;
    double& R      = reaction.integrationGrid.R;
    double& R0     = reaction.opticalPotentialParams.R0;
    double& rC     = reaction.opticalPotentialParams.rC;
    double& rC0    = reaction.opticalPotentialParams.rC0;
    double& rcProj    = reaction.masses.rcProj;
    double& rcTarget    = reaction.masses.rcTarget;
    double& rc0Proj   = reaction.masses.rc0Proj;
    double& rc0Target   = reaction.masses.rc0Target;
    double& rI     = reaction.opticalPotentialParams.rI;
    double& rI0    = reaction.opticalPotentialParams.rI0;
    double& rSo    = reaction.opticalPotentialParams.rSo;
    double& rSo0   = reaction.opticalPotentialParams.rSo0;
    double& rSoi   = reaction.opticalPotentialParams.rSoi;
    double& rSoi0  = reaction.opticalPotentialParams.rSoi0;
    double& rSi    = reaction.opticalPotentialParams.rSi;
    double& rSi0   = reaction.opticalPotentialParams.rSi0;
    double& V      = reaction.opticalPotentialParams.V;
    double& vI     = reaction.opticalPotentialParams.vI;
    double& vSi    = reaction.opticalPotentialParams.vSi;
    double& aSi    = reaction.opticalPotentialParams.aSi;

    int& L       = reaction.angMom.L;
    int& parity  = reaction.angMom.parity;

    double& J    = reaction.angMom.J;
    double& jProj   = reaction.angMom.jProj;

    int& iDone   = reaction.internalState.iDone;
    int& stripPickup  = reaction.internalState.stripPickup;
    double undefValue = reaction.internalState.undefValue;


    // Helper: clear A-type potential parameters (L310 block)
    auto clearA = [&]() {
        A     = undefValue;
        aI    = undefValue;
        aSo   = undefValue;
        aSoi  = undefValue;
        aSi   = undefValue;
        R0    = undefValue;
        rI0   = undefValue;
        rSo0  = undefValue;
        rSoi0 = undefValue;
        rSi0  = undefValue;
        rC0   = undefValue;
        rc0Proj  = undefValue;
        rc0Target  = undefValue;
        asymptopia = undefValue;
    };

    // Helper: check iDone skip conditions (L300 block)
    // Returns true if we should skip to L400 (skip clearA)
    auto idoneSkip = [&]() -> bool {
        const int bit = 1 << (channelIndex - 1);
        return iDone == bit || iDone == 3 + bit || iDone == 12 + bit;
    };

    if (channelIndex >= 5) {
        //
        // CLEAR ONLY THINGS THAT CAN BE REDEFINED LATER
        //
        if (massProj != undefValue && massTgt != undefValue) {
            aM = undefValue;
            if (R0  != undefValue) R    = undefValue;
            if (rI0 != undefValue) rI   = undefValue;
            if (rSo0  != undefValue) rSo  = undefValue;
            if (rSoi0 != undefValue) rSoi = undefValue;
            if (rSi0  != undefValue) rSi  = undefValue;
            if (rC0  != undefValue) rC  = undefValue;
            if (rc0Proj != undefValue) rcProj = undefValue;
            if (rc0Target != undefValue) rcTarget = undefValue;
            // after VTR..VTPI removal). The loop was a 1-based off-by-one
            // (touched RTRI..ATR-1) and the targets are now gone anyway.
            if (eLab != undefValue) E = undefValue;
        }
        if (reaction.energies.eCm != undefValue) E = undefValue;
        // fall through to L400 block below
    } else {
        //
        // WHEN DOING DWBA; FREE LOTS  (L200 block)
        //
        E    = undefValue;
        aM   = undefValue;
        massProj  = undefValue;
        massTgt  = undefValue;
        R    = undefValue;
        rC   = undefValue;
        rcProj  = undefValue;
        rcTarget  = undefValue;
        rI   = undefValue;
        rSo  = undefValue;
        rSoi = undefValue;
        rSi  = undefValue;
        reaction.spec.spam = undefValue;
        L     = NOTDEF_INT;
        reaction.angMom.nNodes = NOTDEF_INT;
        reaction.charges.zProj   = NOTDEF_INT;
        reaction.charges.zTarget   = NOTDEF_INT;
        J     = reaction.internalState.notDefSentinel;
        jProj    = reaction.internalState.notDefSentinel;
        reaction.angMom.spinProj   = reaction.internalState.notDefSentinel;
        reaction.angMom.spinTarget   = reaction.internalState.notDefSentinel;
        parity = 0;
        for (int i = 1; i <= 2; i++) {
            reaction.angMom.parityPt[i] = 0;
            reaction.masses.amxgPt[i] = undefValue;
            reaction.energies.exsPt[i] = undefValue;
        }
        reaction.internalState.lSpcPt2 = NOTDEF_INT;
        reaction.internalState.nodePt2 = NOTDEF_INT;

        //
        // DO NOT SET POTENTIALS TO 0 WHEN GOING FROM ONE SCAT. STATE TO
        // THE hasNextBlock.
        //
        if (!(channelIndex >= 3 && iDone / 4 != 3)) {
            reaction.opticalPotentialParams.vSo  = 0;
            reaction.opticalPotentialParams.vSoi = 0;
            vSi  = 0;
            // (and the original loop was off-by-one anyway: 1-based [1..6] hit
            // VTRI..PARAM1, not VTR..VTPI; masked by the fields all being 0).
            V  = 0;
            vI = 0;

            //
            // CANCEL LINKULE REFERENCES EXCEPT FROM ONE SCAT STATE TO ANOTHER
            // ALSO DO NOT CANCEL PRIOR TO EFFECTIVE EXCITATION POTENTIAL
            //
            // Original: if stripPickup==0 goto L300 (skip linkule loop)
            //           else: do linkule loop, then if channelIndex==0 goto L310
            //                 else fall to L300 (check iDone)
            // if iDone matches → goto L400 (skip clearA)
            // clearA()
            //
            if (stripPickup != 0) {
                for (int i = 1; i <= numLinkules; i++) {
                    reaction.linkuleData.linkuleAddr[i][3] = 0;
                }

                //
                // WE LEAVE THE A'S DEFINED FROM ONE BOUND STATE TO ANOTHER AND
                // FROM ONE SCATTERING TO ANOTHER BUT NOT FROM BOUND TO SCATTER
                //
                if (channelIndex == 0) {
                    // clearA unconditionally
                    clearA();
                } else if (!idoneSkip()) {
                    // fall through: check iDone
                    clearA();
                }
            } else if (!idoneSkip()) {
                // stripPickup == 0: goto L300 (no linkule loop; clear A when iDone allows)
                clearA();
            }
        }
        // fall through to L400 block below
    }

    //
    // WE ALWAYS UNDEFINE PARAMETERS THAT WERE SET TO OTHER PARAMETERS
    // BY SETPOT  (L400 block)
    //
    if (reaction.internalState.wasSet[1]) {
        rI0   = undefValue;
    }
    if (reaction.internalState.wasSet[2]) {
        aI    = undefValue;
    }
    if (reaction.internalState.wasSet[3])  rSo    = undefValue;
    if (reaction.internalState.wasSet[4])  aSo    = undefValue;
    if (reaction.internalState.wasSet[5])  rSoi   = undefValue;
    if (reaction.internalState.wasSet[6])  aSoi   = undefValue;
    if (reaction.internalState.wasSet[7])  rSi    = undefValue;
    if (reaction.internalState.wasSet[8])  aSi    = undefValue;
    if (reaction.internalState.wasSet[9]) {
        rC0    = undefValue;
        rC     = undefValue;
    }
    // 12..16,29; source_potentials only writes 1..9). Both guards always false.
    if (reaction.internalState.wasSet[12]) jProj     = reaction.internalState.notDefSentinel;
    if (reaction.internalState.wasSet[13]) J      = reaction.internalState.notDefSentinel;
    if (reaction.internalState.wasSet[14]) L      = NOTDEF_INT;
    if (reaction.internalState.wasSet[15]) asymptopia = undefValue;
    if (reaction.internalState.wasSet[16]) eLab   = undefValue;
    // RTPI and ATR..ATPI fields all gone. The corresponding wasSet[17..29]
    // odd/even slots are no longer set by anything (source_potentials SETPOT
    // tensor parameter loop was deleted with the VTR..VTPI cascade).
    if (reaction.internalState.wasSet[29]) parity = 0;

    //
    // IF THIS IS INELASTIC SCATTERING, SET BACK TO FIRST CHANNEL
    // AFTER BOTH CHANNELS ARE DONE.  (L800 block)
    //
    if (stripPickup == 0 && iDone == 12) {
        for (int i = 1; i <= numLinkules; i++) {
            for (int ii = 1; ii <= 6; ii++) {
                reaction.linkuleData.linkuleAddr[i][ii] = reaction.linkuleData.lnkAd2[i][1][ii];
            }
        }
        R     = reaction.kin.rScts[1];
        rI    = reaction.distortedWave.channel[1].rI;
        rC    = reaction.kin.rcScts[1];
        A     = reaction.kin.aScts[1];
        aI    = reaction.distortedWave.channel[1].aI;
        V     = reaction.distortedWave.channel[1].v0R;
        vI    = reaction.distortedWave.channel[1].v0I;
        rSi   = reaction.distortedWave.channel[1].rSi;
        aSi   = reaction.distortedWave.channel[1].aSi;
        vSi   = reaction.distortedWave.channel[1].v0Si;
    }

    return;
}
