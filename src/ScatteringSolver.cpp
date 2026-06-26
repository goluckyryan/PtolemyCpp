// scattering.cpp — ScatteringSolver class + scattering matching functions
// Merged from: scattering_solver.cpp, scattering_matching.cpp

#include "Reaction.h"
#include "ptolemy_types.h"
#include "ScatteringSolver.h"
#include "Timing.h"
#include "math/gauss_quadrature.h"
#include "math/spline.h"
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <vector>

// ============================================================================
// Part 1: ScatteringSolver class methods
// ============================================================================


// ScatteringSolver now lives as a member of reaction.distortedWave (composition).
// See include/DistortedWave.h.

// ---------------------------------------------------------------------------
// reallocate(size)
//
// Resize the wavefunction work arrays to `size` elements (0-based vector).
// Element [0] is unused (1-based indexing convention); valid range [1..size-1].
//
// Updates reaction.distortedWave.scatteringSolver.wavRPointer / wavIPointer to point here.
// ---------------------------------------------------------------------------
void ScatteringSolver::reallocate(int size)
{
    if (size <= 0) return;

    // Only grow, never shrink (matches PTOLEMY's behavior: pool never shrinks mid-run)
    if (size > static_cast<int>(wavR.size())) {
        wavR.assign(size, 0.0);
        wavI.assign(size, 0.0);
    }

    // Update WAVCOM pointers (1-based: data()-1 so ptr[1] == wavR[1])
    wavRPointer = wavR.data() - 1;
    wavIPointer = wavI.data() - 1;
}

// ---------------------------------------------------------------------------
// allocateVWork/Sors/Sois/Tens — class-owned scratch vectors
// (allocateRlvs/Imvs/Centr were no-op shims for an earlier migration step
// ---------------------------------------------------------------------------
void ScatteringSolver::allocateVWork(int n)
{
    vWork.assign(n + 1, 0.0);
    vWorkPointer = vWork.data() - 1;
}
void ScatteringSolver::allocateSors(int k, int n, Reaction& reaction)
{
    // soRPointer uses 0-based convention.
    // source_potentials: soR0Pointer[i-1] writes element i-1 (0-based)
    // solvePartialWave: p[i] reads element i (0-based)
    // Grow-only: MAKPOT may write to NSTEPP which can exceed n (gridPointCount).
    // TODO (Option 3): Track actual max index used and size correctly.
    if ((int)sors[k].size() < n + 2) sors[k].resize(std::max(n + 2, 10000), 0.0);
    else std::fill(sors[k].begin(), sors[k].begin() + n + 1, 0.0);
    reaction.distortedWave.channel[k].soRPointer = sors[k].data();  // 0-based: p[0]=sors[0], != nullptr
}
void ScatteringSolver::allocateSois(int k, int n, Reaction& reaction)
{
    if ((int)sois[k].size() < n + 2) sois[k].resize(std::max(n + 2, 10000), 0.0);
    else std::fill(sois[k].begin(), sois[k].begin() + n + 1, 0.0);
    reaction.distortedWave.channel[k].soIPointer = sois[k].data();  // 0-based
}

// leftover as BoundState::reset(); zero callers since main_control died

// ============================================================================
// Part 2: lCritL, sWkb. TMATCH deleted.
// ============================================================================

// sWkb is file-static (only lCritL calls it) — forward decl for use-before-def.
static void sWkb(double fK, int l, double eps, int printLevel,
          double eta, double rC,
          int nSteps, double roStart, double h,
          double* vArray, double* vCoefsFlat, double* wArray, double* wCoefsFlat,
          int pointCount, double* points, double* weights,
          double& deltaReal, double& deltaImag, double& theta, double& thetaC,
          double& rhoT, double& roTc);

// ============================================================================
//
// USE CLASSICAL DEFLECTION AND WKB PHASES TO ESTIMATE CRITICAL L
// LC IS TAKEN TO BE GREATER OF TWO ESTIMATES lC1 AND lC2
// lC1 IS THE FIRST L FOR WHICH |S(L)| > (|S|MIN + 1)/2.
// lC2 IS L-VALUE OF ALGEBRAICALLY MINIMUM THETA(L)
// WHERE THETA(L) IS THE CLASSICAL DEFLECTION FUNCTION.
// ============================================================================
// rStart parameter dropped — sole caller passes 0.0; roStart folded to 0.
// NGAUSS parameter dropped — sole caller passes 48; now a function-local
// constexpr.
void lCritL(double fK, double eta, double rC, int nSteps,
            double stepSize, const double* vInPointer, const double* wInPointer,
            int printLevel, int& lC1, int& lC2, int& lC)
{

    // Sole caller passed 48; baked in here so sWkb's pointCount arg stays a
    // local constant.
    constexpr int nGauss = 48;

    // Local variables
    int isPastDeflectionMin;
    double eps = 1.0e-12;

    // SET MINIMUM AND MAXIMUM L FOR SEARCH
    // row IS THE L WHERE COULOMB CLASSICAL TURNING POINT IS MAX(rC,1.2)

    isPastDeflectionMin = FALSE_F;
    double rcl = std::max(rC, 1.20e0);
    double row = fK * rcl;
    row = row * row - 2 * eta * row + .25;
    if (row < 0) row = 0;
    if (row > 0) row = 0.5 + std::sqrt(row);
    int lMin = (int)(.2 * row - 4.);
    lMin = std::max(lMin, 0);
    int lMax = (int)(1.5 * row);
    lMax = std::max(lMax, lMin + 50);
    lC1 = 0;
    lC2 = 0;
    lC = 0;
    double dlSMatrix = 0;
    double dlDflC = 0;
    int lCount = lMax - lMin + 1;
    double ampMin = 1.0e20;


    // gaussPointsVector/gaussWeightsVector: Gauss points and weights
    // gaussL uses 1-based indexing: X[1..N], so allocate N+1 elements
    std::vector<double> gaussPointsVector(nGauss + 1, 0.0);
    std::vector<double> gaussWeightsVector(nGauss + 1, 0.0);

    // semiclassicalVector: 5 sub-arrays of lCount elements each — amplitude, nuclear phase,
    //          deflection function, classical TP (nuclear), classical TP (Coulomb)
    int semiCount = 5 * lCount;
    std::vector<double> semiclassicalVector(semiCount + 1, 0.0);  // 0-based via semiPointer
    double* semiPointer = semiclassicalVector.data();  // semiPointer[0..semiCount-1] valid
    // Sub-array base offsets (0-based section starts; paired with 0-based loop index i):
    int ampStart  = 1;                // 1-based section start (amp); lAmp below is 0-based
    int lAmp  = ampStart - 1;        // = 0: amp[i] = semiPointer[lAmp + i] (0-based i)
    int nuclearPhaseStart = ampStart + lCount;     // nuclear-phase sub-array
    int ldNuc = nuclearPhaseStart - 1;
    int deflectionStart = nuclearPhaseStart + lCount;    // deflection function
    int lDeflection = deflectionStart - 1;
    int nuclearCtpStart  = deflectionStart + lCount;    // nuclear classical TP
    int lCtp  = nuclearCtpStart - 1;
    int coulombCtpStart = nuclearCtpStart + lCount;     // Coulomb classical TP
    int lCtpc = coulombCtpStart - 1;

    // wkbSplineVector: rho-grid + cubic spline coefficients for real and imaginary
    //          potentials (7 sections of nSteps elements each)
    std::vector<double> wkbSplineVector(7 * nSteps + 1, 0.0);  // 0-based via cubPointer
    double* cubPointer = wkbSplineVector.data();  // cubPointer[0..7*nSteps-1] valid
    int lRval = 0;               // rho-grid section starts at index 0
    int lvCubic = lRval + nSteps;  // real-V cubic coeff sections
    int lwCubic = lvCubic + 3 * nSteps;  // imag-V cubic coeff sections
    // vPointer/wPointer passed in by caller (0-based class-owned vectors).
    // const_cast: legacy sWkb signature is non-const; data is treated as read-only here.
    double* vPointer = const_cast<double*>(vInPointer);
    double* wPointer = const_cast<double*>(wInPointer);

    // DETERMINE UNMAPPED GAUSS POINTS FOR (-1,1)

    gaussL(nGauss, gaussPointsVector.data(), gaussWeightsVector.data());

    // GET CUBIC SPLINES FOR THE POTENTIALS

    // rStart parameter dropped (sole caller passes 0.0) — roStart folded to 0.
    constexpr double roStart = 0.0;
    double rho = roStart;
    double h = fK * stepSize;
    for (int i = 1; i <= nSteps; i++) {
        cubPointer[lRval - 1 + i] = rho;  // = cubPointer[i-1] = wkbSplineVector[i-1]
        rho = rho + h;
    }
    double roEnd = rho - h - eps;

    // naturalCubicSpline expects: nSteps, rho[], V[], cubic_a[], cubic_b[], cubic_c[]
    //   all as ALLOC_base-style pointers (ptr[0]=first element)
    naturalCubicSpline(nSteps, &cubPointer[lRval], vPointer, &cubPointer[lvCubic],
           &cubPointer[lvCubic + nSteps], &cubPointer[lvCubic + 2 * nSteps]);

    naturalCubicSpline(nSteps, &cubPointer[lRval], wPointer, &cubPointer[lwCubic],
           &cubPointer[lwCubic + nSteps], &cubPointer[lwCubic + 2 * nSteps]);

    // CALCULATE AND STORE WKB S-MATRIX ELEMENTS (AMP AND ARG)
    // CLASSICAL DEFLECTION FUNCTIONS(DEFL) AND CLASSICAL TPS(CTP CTPC)

    // rhoT init silences -Wmaybe-uninitialized at sWkb's `x = rhoT - rho;`
    // (~line 454). sWkb's outer while always exits with rhoT set — either
    // the uEff>1 break (then if(rho>0) sets rhoT=rho) or the rho<=0 break
    // (which sets rhoT=0 or h directly) — but GCC analyzes through the
    // by-ref param and can't see the cases are exhaustive.
    double deltaReal, deltaImag, theta, thetaC, roTc;
    double rhoT = 0;
    double amplitude;
    int l;

    for (int i = 0; i < lCount; i++) {
        l = lMin + i;
        // sWkb uses 0-based array macros (vArr[N-1], points[i-1])
        sWkb(fK, l, eps, printLevel, eta, rC,
             nSteps, roStart, h, vPointer, &cubPointer[lvCubic], wPointer,
             &cubPointer[lwCubic], nGauss, gaussPointsVector.data(), gaussWeightsVector.data(),
             deltaReal, deltaImag, theta, thetaC, rhoT, roTc);
        if (!((deltaReal != 5.0e0) || (deltaImag != 5.0e0))) {
            std::printf("0**** ERROR REturn IN LCRITL FROM SWKB \n");
            break;  // goto L640 — exit loop, dlSMatrix/dlDflC stay 0
        }
        semiPointer[lCtp + i] = rhoT;
        semiPointer[lCtpc + i] = roTc;
        amplitude = 0;
        if (deltaImag < 20) amplitude = std::exp(-2.0e0 * deltaImag);
        semiPointer[lAmp + i] = amplitude;
        ampMin = std::min(ampMin, amplitude);
        if (theta > thetaC) theta = thetaC;
        semiPointer[lDeflection + i] = theta;
        semiPointer[ldNuc + i] = 2 * deltaReal;
        bool hitRoEnd = (rhoT >= roEnd);
        if (!hitRoEnd && !(amplitude < .05 || i == 0)) {
            if (!isPastDeflectionMin) {
                isPastDeflectionMin = (theta > semiPointer[lDeflection + i - 1] + eps);
            } else if (!(theta > semiPointer[lDeflection + i - 1] - eps)) {
                // Deflection function is falling from max after min; |S| close to 1: stop
                if (!(amplitude < .95 + .05 * ampMin))
                    hitRoEnd = true;  // lMax = L; break
            }
        }
        if (hitRoEnd) {
            lMax = l;
            break;  // goto L300
        }
    }

    lCount = lMax - lMin + 1;

    // PRINT WKB PHASES AND CLASSICAL DEFLECTION FUNCTION

    if (printLevel != 0) {
        std::printf("\n     ---  WKB S-MATRIX ELEMENTS AND CLASSICAL TPS ---\n"
                    "  L      AMP(S(L))        ARG(S(L))        DEFL.ANGLE "
                    "       CLASSICAL TP      COULOMB TP\n");
        for (int i = 0; i < lCount; i++) {
            l = lMin + i;
            double ctr  = semiPointer[lCtp  + i] / fK;
            double ctrC = semiPointer[lCtpc + i] / fK;
            double amp  = semiPointer[lAmp  + i];
            double nuclearPhase = semiPointer[ldNuc + i];
            double deflection = semiPointer[lDeflection + i];
            std::printf("  %3d%16.5G%16.5G%16.5G%16.5G%16.5G\n", l, amp, nuclearPhase, deflection, ctr, ctrC);
        }
    }

    // ESTIMATE CRITICAL L FROM PHASE-SHIFTS AND DEFLECTION FN

    isPastDeflectionMin = FALSE_F;
    ampMin = 0.5 * (ampMin + 1.);

    for (int i = 1; i <= lCount; i++) {
        int i1 = lCount - i;   // 0-based index counting down from the top L
        l = lMax + 1 - i;
        // Find dlSMatrix: half-amplitude point
        if (dlSMatrix == 0 && !(semiPointer[lAmp + i1] > ampMin) && !(ampMin > .99)) {
            if (i <= 1) {
                std::printf("\n**** LCRITL PICKED LMAX TOO SMALL TO FIND"
                            " |S| > .5*(MIN(|S|)+1):%5d%15.5G"
                            "  lCrit SET TO LMAX\n", lMax, semiPointer[lAmp + i1]);
                dlSMatrix = lMax;
            } else {
                dlSMatrix = l + (ampMin - semiPointer[lAmp + i1]) / (semiPointer[lAmp + i1 + 1]
                         - semiPointer[lAmp + i1]);
            }
        }
        // Find dlDflC: outermost minimum of deflection function
        if (dlDflC == 0 && l != lMax) {
            if (!isPastDeflectionMin) {
                isPastDeflectionMin = (semiPointer[lDeflection + i1] < semiPointer[lDeflection + i1 + 1] - eps);
            } else if (!(semiPointer[lDeflection + i1] < semiPointer[lDeflection + i1 + 1] + eps)) {
                double c = semiPointer[lDeflection + i1] - 2 * semiPointer[lDeflection + i1 + 1]
                         + semiPointer[lDeflection + i1 + 2];
                dlDflC = l + 1 + (semiPointer[lDeflection + i1] - semiPointer[lDeflection + i1 + 2]) / (2 * c);
            }
        }
    }

    lC1 = (int)(dlSMatrix + .5);
    lC2 = (int)(dlDflC + .5);
    lC = std::max(lC1, lC2);

    if (printLevel != 0)
        std::printf("\n   --- 1/2 POINT OF |S(L)|  =  %15.5G"
                    "     MINIMUM OF DEFLECTION FUNC =%15.5G"
                    "     L CRIT =%5d\n", dlSMatrix, dlDflC, lC);

    // if no critical L found, estimate from classical Coulomb TP
    if (dlSMatrix == 0 && dlDflC == 0) {
        lC = (int)(row + .5);
        std::printf("\n**** NO CRITICAL L FOUND; SET LC =%5d"
                    " = L FOR WHICH CLASSICAL COULOMB TURNING POINT =%10.3f\n", lC, rcl);
    }
    // No ALLOC frees needed — all work arrays are local std::vectors
    // (gaussPointsVector, gaussWeightsVector, semiclassicalVector, wkbSplineVector auto-destruct on return)
    return;
}


// ============================================================================
//
// EVALUATES WKB PHASE-SHIFT FOR PARTIAL WAVE L
// AND THE CORRESPONDING CLASSICAL DEFLECTION FUNCTION
// ============================================================================
static void sWkb(double fK, int l, double eps, int printLevel,
          double eta, double rC,
          int nSteps, double roStart, double h,
          double* vArray, double* vCoefsFlat, double* wArray, double* wCoefsFlat,
          int pointCount, double* points, double* weights,
          double& deltaReal, double& deltaImag, double& theta, double& thetaC,
          double& rhoT, double& roTc)
{
    // duEff out-param dropped — never read by lCritL (sole caller); used
    // only inside sWkb for the L11/CTP debug printf. Now a plain local.
    double duEff;
    // VCOEFS and WCOEFS are (nSteps,3) arrays — access as VCOEFS[col*nSteps + row]
    // But in Fortran column-major: VCOEFS(N,J) -> vCoefsFlat[(J-1)*nSteps + (N-1)]
    // Using 1-based pointer arithmetic: element (N,J) is at offset (J-1)*nSteps + (N-1)
    // We'll define a macro-like lambda for access
    // VCOEFS(N,J) -> vCoefsFlat[ (J-1)*nSteps + N - 1 ]  (1-based N, 1-based J)
    // which are 1-based. So the pointer already points to element [1].
    // In the Fortran, VARRAY(N) means the N-th element starting from 1.
    // The C++ pointer starts at element [0] of the passed array.
    // So vArray[N-1] = Fortran VARRAY(N).

    // Helper macros (0-based pointer, 1-based indices)
    #define vArr(n) vArray[(n)-1]
    #define wArr(n) wArray[(n)-1]
    #define vCoef(n,j) vCoefsFlat[((j)-1)*nSteps + (n)-1]
    #define wCoef(n,j) wCoefsFlat[((j)-1)*nSteps + (n)-1]

    double pi = 3.14159265358979300e0;
    double tp = 1.0 / 1.41421356237309500e0;   // = 1/sqrt(2)
    int iterMax = 50;
    double rhoC = fK * rC;
    double uC0 = (2 * eta) / rhoC;

    // NOTE QUANTUM-MECHANICAL VALUE.  THIS AVOIDS SOME TROUBLES

    double fLl = l * (l + 1.0e0);
    double ecc = std::sqrt(eta * eta + fLl);
    deltaReal = 5.0e0;
    deltaImag = 5.0e0;

    // OUTER LIMIT OF THE WKB INTEGRALS WILL BE ASYMPTOPIA

    double roMax = roStart + h * (nSteps - 1);

    // CALCULATE CTP IN PURE COULOMB POTENTIAL

    roTc = eta + ecc;
    if (roTc < rhoC && eta != 0) {
        double f1 = (3 * eta - rhoC) / (2 * eta);
        roTc = std::sqrt(f1 * f1 + fLl / (eta * rhoC)) + f1;
        if (roTc < 0) roTc = 0;
        roTc = std::sqrt(roTc) * rhoC;
    }

    // CALCULATE CTP IN FULL(NUCLEAR+COULOMB) POTENTIAL
    // FIND INTERVAL (RHO,RHO-H) CONTAINING CTP

    { double uLast = 0;
      double rho = roMax;
      int n = nSteps;
      int iterCount = 0;
      double uEff;
      double x;

      // Scan inward to find interval containing CTP
      while (true) {
          uEff = ((12.0 / (h * h)) + 1.0) - (12.0 / (h * h)) * vArr(n) + fLl / (rho * rho);
          if (uEff > 1) break;  // found bracket
          uLast = uEff;
          rho = rho - h;
          n = n - 1;
          if (rho <= 0) {
              rhoT = 0;
              rho = 0;
              if (l != 0) {
                  std::printf("\n   --- RHOT < H FOR L = %3d ---\n", l);
                  rhoT = h;
              }
              break;
          }
      }

      if (rho > 0) {
          // uEff > 1 at RHO: CTP is between RHO and RHO+H
          rhoT = rho;
          if (rho < roMax && std::fabs(uEff - 1.0e0) > eps) {
              // USE REGULA FALSI TO LOCATE CTP WITHIN (RHO,RHO+H)
              double rO1 = rho;
              double rO2 = rho + h;
              double u1 = uEff;
              double u2 = uLast;
              double rO3, rO, acc;
              while (true) {
                  rO3 = (rO2 * (u1 - 1) + rO1 * (1 - u2)) / (u1 - u2);
                  rO = rO3;
                  x = (rO - roStart - (n - 1) * h);
                  uEff = ((12.0 / (h * h)) + 1.0) - (12.0 / (h * h)) * (vArr(n)
                      + x * (vCoef(n, 1) + x * (vCoef(n, 2) + x * vCoef(n, 3))))
                      + fLl / (rO * rO);
                  if (std::fabs(uEff - 1.0e0) <= eps) {
                      rhoT = rO3;
                      break;
                  }
                  if (uEff >= 1) { rO1 = rO3; u1 = uEff; }
                  else            { rO2 = rO3; u2 = uEff; }
                  iterCount = iterCount + 1;
                  if (iterCount > iterMax) {
                      acc = std::fabs(uEff - 1.0e0);
                      std::printf("\n**** SWKB: FOR L = %3d AFTER %4d ITERATIONS\n"
                                  "   REQUESTED ACCURACY OF %12.4E CANNOT BE ACHIEVED "
                                  "   ACCURACY ACHIEVED = %12.4E\n", l, iterMax, eps, acc);
                      rhoT = rO3;
                      break;
                  }
              }
          }
      }

      // EVALUATE D(uEff)/D(RHO) AT CTP

      // L11: evaluate duEff at CTP
      x = (rhoT - rho);
      { double duN = -(12.0 / (h * h)) * (vCoef(n, 1) + x * (2 * vCoef(n, 2) +
                3 * x * vCoef(n, 3)));
        double duL = 0;
        if (l != 0) duL = -(2 * fLl) / (rhoT * rhoT * rhoT);
        duEff = duN + duL;
        double rT = rhoT / fK;
        if (printLevel > 1)
            std::printf("\n   --- CTP FOR L = %3d FOUND IN %4d"
                        " ITERATIONS   RHOT = %13.5G RT = %13.5G\n"
                        " DUEFF =%14.5G    DUN+DUC = %13.5G DUL = %13.5G\n",
                        l, iterCount, rhoT, rT, duEff, duN, duL);
        if (duEff >= 0 && rhoT != 0)
            std::printf("\n****** DUEFF=%15.5G ERROR IN SWKB ***\n", duEff);
      }

      // INTEGRATE FROM rhoT TO roMax TO EVALUATE WKB PHASE-SHIFT
      // AND CLASSICAL DEFLECTION FUNCTION

      // WKB phase-shift integration
      { double phaseC, temp, phi;
        if (eta == 0) {
            phaseC = 0;
            temp = std::sqrt(fLl) / roMax;
        } else {
            temp = eta / ecc;
            phaseC = std::acos(temp);
            temp = (eta + fLl / roMax) / ecc;
        }
        if (temp > 1) temp = 1;
        phi = std::acos(temp);

                double turn;
        if (roMax <= rhoT) {
            deltaReal = 0;
            deltaImag = 0;
            turn = phaseC;
        } else {
        { double aP = (roMax - rhoT) / 2.0;
          double bP = (roMax + rhoT) / 2.0;
          double aD = roMax / 2.0;
          double bD = 0;
          if (rhoT != 0) {
              aD = 1.0 / rhoT - 1.0 / roMax;
              aD = .5 * std::sqrt(aD);
              bD = 1.0 / rhoT;
          }
          double tempReal = 0;
          double tempImag = 0;
          double temp2 = 0;

          for (int i = 1; i <= pointCount; i++) {

              double wtP = aP * weights[i - 1];
              double wtD = aD * weights[i - 1];
              double rO = aP * points[i - 1] + bP;
              double xVal = aD * (1 + points[i - 1]);
              n = (int)((rO - roStart) / h) + 1;
              x = (rO - roStart - (n - 1) * h);
              uEff = ((12.0 / (h * h)) + 1.0) - (12.0 / (h * h)) * (vArr(n)
                  + x * (vCoef(n, 1) + x * (vCoef(n, 2) + x * vCoef(n, 3))))
                  + fLl / (rO * rO);
              // --- return to label 185 ---
              double uEffC = fLl / (rO * rO);
              double y = rO / rhoC;
              if (rO < rhoC) uEffC = uEffC + uC0 * (1.5 - .5 * (y * y));
              if (rO >= rhoC) uEffC = uEffC + uC0 / y;
              double w = -(12.0 / (h * h)) * (wArr(n) +
                  x * (wCoef(n, 1) + x * (wCoef(n, 2) + x * wCoef(n, 3))));
              double uE1 = 1 - uEff;
              if (uE1 < 0) uE1 = 0;
              double uE2 = std::sqrt(uE1 * uE1 + w * w);
              if ((uE2 - uE1) < 0) uE2 = uE1;
              double fIntI = tp * std::sqrt(uE2 - uE1);
              double fIntR = tp * std::sqrt(uE2 + uE1);
              if (rO > roTc) {
                  double uEc = 1 - uEffC;
                  if (uEc < 0) uEc = 0;
                  fIntR = fIntR - std::sqrt(uEc);
              }
              if (rhoT != 0) {
                  rO = bD - xVal * xVal;
                  rO = 1.0 / rO;
              } else {
                  rO = .5 * aD * (1 + points[i - 1]) * (1 + points[i - 1]);
              }
              n = (int)((rO - roStart) / h) + 1;
              x = (rO - roStart - (n - 1) * h);
              uEff = ((12.0 / (h * h)) + 1.0) - (12.0 / (h * h)) * (vArr(n)
                  + x * (vCoef(n, 1) + x * (vCoef(n, 2) + x * vCoef(n, 3))))
                  + fLl / (rO * rO);
              // --- return to label 200 ---
              { double uEr = 1 - uEff;
                if (uEr <= 0) {
                    std::printf("\n**** SWKB ERROR AT STMNT 205:%5d%5d\n"
                                "      %15.5G%15.5G%15.5G%15.5G%15.5G%15.5G\n",
                                l, i, rhoT, rhoT/fK, xVal, rO, uEff, uEr);
                    return;
                }
                double fIntD = xVal / std::sqrt(uEr);
                tempImag = tempImag + wtP * fIntI;
                tempReal = tempReal + wtP * fIntR;
                temp2 = temp2 + wtD * fIntD;
              }
              continue;

          }

          deltaReal = tempReal;
          deltaImag = tempImag;
          turn = 2 * (l + .5) * temp2 + phaseC - phi;
        } // end Gauss integration
        } // end if (roMax > rhoT)

        // L230_inline: merge point
        theta = pi - 2 * turn;
        thetaC = pi - 2 * phaseC;
      }
    } // end uLast scope

    #undef vArr
    #undef wArr
    #undef vCoef
    #undef wCoef

    return;
}


// TMATCH deleted — tensor-coupling S-matrix renormalizer (used to set up
// the 4-by-4 → 2-by-2 reduction for the waveReal/waveImag matched waves and the
// per-(jProj, L) S-matrix entries via convertJtoL). Sole call site was the
// other callers after that purge.

