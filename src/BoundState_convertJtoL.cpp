// BoundState_convertJtoL.cpp — JPTOLX: translates the elastic S-matrix from the
// (L',L,jProj) system to (L',L,lx); spin-orbit case multiplies by Racah coeffs.

#include "math/angular_momentum_coeff.h"
#include "BoundState.h"
#include "Reaction.h"

void BoundState::convertJtoL(int L, int jProj, int channelIndex, double sJr, double sJi,
            int* indxePointer, double* sLx, Reaction& reaction)
{
    // IF NO SPIN-ORBIT, JUST TRANSFER IT.

    if ( !reaction.distortedWave.channel[channelIndex].hasSpinorbit ) {
        int i = 2*L;  // 0-based: sLx[2*L]/sLx[2*L+1] (matches smatArr[2*L] elsewhere)
        sLx[i] = sJr;
        sLx[i+1] = sJi;
        return;
    }

    // SPIN-ORBIT EXISTS.  MULTIPLY BY RACAH'S ETC., AND ADD
    // TO EXISTING VALUES FOR EACH lx.
    int lxMax = std::min( reaction.distortedWave.channel[channelIndex].twoSpin, 2*L );
    double coefficient = (reaction.distortedWave.channel[channelIndex].twoSpin + 1) * (2*L + 1);
    coefficient = (jProj + 1) / std::sqrt(coefficient);
    if ( ((( reaction.distortedWave.channel[channelIndex].twoSpin - jProj + 2*L ) / 2) % (2)) != 0 ) coefficient = -coefficient;
    for (int lx = 0; lx <= lxMax; lx++) {
        int lx2 = 2*lx;
        double t = coefficient * std::sqrt((double)(lx2 + 1)) *
            racah( 2*L, 2*L, reaction.distortedWave.channel[channelIndex].twoSpin, reaction.distortedWave.channel[channelIndex].twoSpin,
                   lx2, jProj );
        int k = 1 + lx*( reaction.distortedWave.channel[channelIndex].twoSpin + 2 );
        // indxePointer is 0-based pointer (Fortran INDX(3,*) → C++ indxePointer[(k-1)*3 + (i-1)])
        if ( indxePointer[(k-1)*3 + 0] <= 0 )  continue;
        int i = L*reaction.distortedWave.channel[channelIndex].nJStates + indxePointer[(k-1)*3 + 0] + (-indxePointer[(k-1)*3 + 1]) / 2;
        i = 2*i - 2;  // 0-based: sLx is now smatArr.data()
        sLx[i] = sLx[i] + t*sJr;
        sLx[i+1] = sLx[i+1] + t*sJi;
    }
}
