// Contains: aitken, BKGPTElp, DEFORMEd, FIXEDWOo, GAUSSIAn,
//           JDEPEN, JDEPENWS, LAGRANGE, LTSTELp, PARITWOO,
//           SHAPE, SPLINE, TWOSHApe
//
// NOTE: These routines access shared state via the Reaction struct
//       (FLOAT_common, INTGER, CNSTNT, etc.) rather than through
//       the passed array parameters, to avoid struct padding issues.
//

#include "linkule.h"
#include "math/spline.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Constants.h"
#include "Reaction.h"
#include "OpticalPotential.h"
#include "LinkulePlugin.h"

// Factories at end of file; the linkule() switch constructs via them.
struct GaussianLinkulePlugin : LinkulePlugin {
    void run(char8 alias, int* linkuleInts, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, Reaction& reaction) override;
};
struct FixedWoodsSaxonLinkulePlugin : LinkulePlugin {
    void run(char8 alias, int* linkuleInts, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, Reaction& reaction) override;
};
struct JDependentWoodsSaxonLinkulePlugin : LinkulePlugin {
    void run(char8 alias, int* linkuleInts, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, Reaction& reaction) override;
};
struct JDependentWoodsSaxonFermiLinkulePlugin : LinkulePlugin {
    void run(char8 alias, int* linkuleInts, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, Reaction& reaction) override;
};
struct ParityWoodsSaxonLinkulePlugin : LinkulePlugin {
    void run(char8 alias, int* linkuleInts, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, Reaction& reaction) override;
};
struct ShapeLinkulePlugin : LinkulePlugin {
    void run(char8 alias, int* linkuleInts, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, Reaction& reaction) override;
};
struct SplineLinkulePlugin : LinkulePlugin {
    void run(char8 alias, int* linkuleInts, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, Reaction& reaction) override;
};
struct LagrangeLinkulePlugin : LinkulePlugin {
    void run(char8 alias, int* linkuleInts, int potType, int requestCode,
             int& callStatus, int L, double& J, double rStart, double stepSize,
             int nPts, double* array1, double* array2, Reaction& reaction) override;
};

// params(reaction) — 0-based pointer to the PARAM_arr[0..20] block in
// callers; every other reader uses the named linkuleParams.PARAM_arr[N-1]
// single std::array<double, 21> instead of PARAM1..PARAM5 + PAR620[16].
namespace {
double* params(Reaction& r) { return r.linkuleParams.PARAM_arr.data(); }

// Shared validation-failure report, emitted byte-identically by
// fixedWoodsSaxon / jDependentWoodsSaxon / parityWoodsSaxon when R, V, A
// or PARAM1 is undefined. The `return;` stays at each call site.
void reportRvaParam1Undefined(int& callStatus) {
    printf("\n**** R, V, A, AND PARAM1 MUST BE DEFINED.\n");
    callStatus = -1;
}

// Shared two-point validation-failure report, emitted byte-identically by
// splineLinkule ("SPLINE") and lagrange ("LAGRANGE") at both their pair-count
// guards. The `return;` stays at each call site.
void reportTwoPointsRequired(const char* method, int& callStatus) {
    printf("\n**** POTENTIAL MUST BE DEFINED ON AT LEAST TWO POINTS FOR %s.\n", method);
    callStatus = -1;
}

// Shared positive-radius/diffuseness guard, emitted byte-identically by
// jDependentWoodsSaxon / jDependentWoodsSaxonFermi / parityWoodsSaxon right
// after the undefined-parameter check. Returns true (the caller's `return;`
// stays at the call site) when r or a is non-positive.
bool bailIfBadRadiusOrDiffuse(double r, double a, int& callStatus) {
    if (r <= 0.0 || a <= 0.0) { callStatus = -1; return true; }
    return false;
}

// Shared point-collection preamble of splineLinkule ("SPLINE") and lagrange
// ("LAGRANGE"). param is a 2x10 array (param[2*k]=R, param[2*k+1]=V); reads up
// to 10 pairs, drops any with an undefined R or V, and insertion-sorts the
// survivors into rPoints[1..]/vPoints[1..] ascending in R. Returns the pair
// count, or -1 after reporting (the caller's `return;` stays at the call site)
// when no first pair is found or fewer than 3 valid pairs survive.
int collectOrderedPoints(double* param, double* rPoints, double* vPoints,
                         double undefValue, const char* method, int& callStatus) {
    int pairCount = 0;
    int pointIndex = 0;

    // Find first valid pair
    while (pointIndex < 10) {
        double rCandidate = param[2 * pointIndex];
        double vCandidate = param[2 * pointIndex + 1];
        if (rCandidate != undefValue && vCandidate != undefValue) {
            pairCount = 1;
            rPoints[1] = rCandidate;
            vPoints[1] = vCandidate;
            pointIndex++;
            break;
        }
        pointIndex++;
    }
    if (pairCount == 0) {
        reportTwoPointsRequired(method, callStatus);
        return -1;
    }

    // Find remaining valid pairs and insert in order
    while (pointIndex < 10) {
        double rCandidate = param[2 * pointIndex];
        double vCandidate = param[2 * pointIndex + 1];
        pointIndex++;
        if (rCandidate == undefValue || vCandidate == undefValue) continue;

        // Find insertion point
        int ins = pairCount + 1;
        for (int i = 1; i <= pairCount; i++) {
            if (rCandidate < rPoints[i]) {
                ins = i;
                break;
            }
        }
        // Shift right
        for (int ii = pairCount; ii >= ins; ii--) {
            rPoints[ii + 1] = rPoints[ii];
            vPoints[ii + 1] = vPoints[ii];
        }
        pairCount++;
        rPoints[ins] = rCandidate;
        vPoints[ins] = vCandidate;
    }

    if (pairCount <= 2) {
        reportTwoPointsRequired(method, callStatus);
        return -1;
    }
    return pairCount;
}

// requestCode==3 "set up the two pieces": generate a fresh Woods-Saxon shape
// into a scratch workarrs slot, then copy it into array1[1..nPts]. Byte-identical
// across the three J/L-dependent fitters (jDependentWoodsSaxon,
// jDependentWoodsSaxonFermi, parityWoodsSaxon); the caller's `return;` stays inline.
void setupWoodsSaxonPiece(Reaction& reaction, const int* linkuleInts, int nPts,
                          double rStart, double stepSize, double* array1,
                          double v, double r, double a) {
    std::vector<double>& workSlot = reaction.linkuleData.workarrs[linkuleInts[0] - 1];  // workarrs (1-based index)
    double* workArrayPointer = workSlot.data() - 1;
    // resize() zeroes the buffer so ADD-semantics fillWoodsSaxon == old SET; the
    // filled slot is reused by requestCode==4 (J-dependent add), so write it.
    OpticalPotential pot;
    pot.resize(nPts, rStart, stepSize);
    pot.fillWoodsSaxon(v, r, a);
    std::copy(pot.values.begin(), pot.values.end(), workSlot.begin());
    for (int i = 1; i <= nPts; i++) {
        array1[i] = workArrayPointer[i];
    }
}

// Read the Woods-Saxon depth (v), diffuseness (a), and radius (r) for this
// potType from the scattering-solver scratch — TEMPVS(potType)/TEMPVS(4+potType)
// and FLOAT(43)/FLOAT(47). Pure field reads, shared by the WS fitter variants
// (gaussian, jDependentWoodsSaxon, jDependentWoodsSaxonFermi, parityWoodsSaxon).
void readWoodsSaxonVAR(Reaction& reaction, int potType, double& v, double& a, double& r) {
    v = (&reaction.distortedWave.scatteringSolver.potentialWork.tvReal)[potType - 1];  // TEMPVS(potType)
    a = (&reaction.distortedWave.scatteringSolver.potentialWork.tvReal)[3 + potType];   // TEMPVS(4+potType)
    r = (potType == 2) ? reaction.opticalPotentialParams.rI : reaction.integrationGrid.R;     // FLOAT(43)/FLOAT(47)
}

// Append a zero-filled nPts-long work array to LinkuleData.workarrs and return
// its 1-based slot index (the value the requestCode==3/4 readers index by).
int allocateWorkSlot(Reaction& reaction, int nPts) {
    reaction.linkuleData.workarrs.emplace_back(nPts, 0.0);
    return (int)reaction.linkuleData.workarrs.size();
}

// Append a fresh nPts-long work array, fill it with the radial coordinate grid
// (rStart, rStart+stepSize, ...), and store its 1-based slot index in
// linkuleInts[0] (the value the requestCode==3 readers index by); linkuleInts[1]=0.
// Shared by splineLinkule()/lagrange() requestCode==1.
void buildRadialGridSlot(Reaction& reaction, int* linkuleInts, int nPts,
                         double rStart, double stepSize) {
    reaction.linkuleData.workarrs.emplace_back(nPts, 0.0);
    double* arr = reaction.linkuleData.workarrs.back().data();
    double rr = rStart;
    for (int ii = 0; ii < nPts; ii++) {
        arr[ii] = rr;
        rr = rr + stepSize;
    }
    linkuleInts[0] = (int)reaction.linkuleData.workarrs.size();
    linkuleInts[1] = 0;
}

// requestCode==3 "Setup output X's": clamp array1 below rPoints[1] to -vPoints[1],
// zero it above rPoints[pairCount], and record the [startIndex, lastIndex] window
// that the interpolation will fill. Pure boundary/index bookkeeping over the
// radial grid (no potential math). Shared by splineLinkule()/lagrange()
// requestCode==3 (byte-identical). Returns outputCount; startIndex/lastIndex
// returned by reference.
int setupOutputWindow(double* array1, const double* rPoints, const double* vPoints,
                      int pairCount, int nPts, double rStart, double stepSize,
                      int& startIndex, int& lastIndex) {
    double rr = rStart;
    startIndex = 1;
    lastIndex = 1;
    for (int n = 1; n <= nPts; n++) {
        if (rr < rPoints[1]) {
            array1[n] = -vPoints[1];
            startIndex = n + 1;
        } else if (rr > rPoints[pairCount]) {
            array1[n] = 0.0;
        } else {
            lastIndex = n;
        }
        rr = rr + stepSize;
    }
    return lastIndex + 1 - startIndex;
}

// requestCode==2 point-table dump: print each (R, V) pair five-per-row. Pure
// output (no math). Shared by splineLinkule()/lagrange() case 2 (byte-identical
// loop; only the preceding header line differs, kept at the call site).
void printPointTable(const double* rPoints, const double* vPoints, int pairCount) {
    for (int i = 1; i <= pairCount; i++) {
        printf("%10.3f%13.5g", rPoints[i], vPoints[i]);
        if (i % 5 == 0 || i == pairCount) printf("\n");
    }
}
}

// The DEFORMED linkule stub still returns "not yet operational".


// ============================================================================
// Nearest surrounding points are used, biased toward center.
//
// May 15, 1976 - first version - S. Pieper
// Nov 5, 1976 - always use surrounding points
// ============================================================================

static void aitken(int maxDegree, int inputCount,
            double* xIn, double* fIn,
            int outputCount, double* xOut, double* fOut, double* work)
{
    // work is dimensioned (2, maxDegree+1) minimum — treated as work[2*(maxDegree+1)]
    // work(1,j) = work[2*(j-1)], work(2,j) = work[2*(j-1)+1]
    // Using 1-based: work1(j) = work[2*(j-1)], work2(j) = work[2*(j-1)+1]
    #define work1(j) work[2*((j)-1)]
    #define work2(j) work[2*((j)-1)+1]

    bool isMovedUp;

    int nearestIndex = 1;
    double xMid = xIn[(inputCount + 1) / 2];
    int degree = std::min(maxDegree, inputCount - 1);

    for (int outputIndex = 1; outputIndex <= outputCount; outputIndex++) {
        double xTarget = xOut[outputIndex];
        double fInterp;

        // Find nearest point
        isMovedUp = false;
        double bestDistance = std::fabs(xTarget - xIn[nearestIndex]);

        // Scan forward as long as next xIn is closer
        while (nearestIndex != inputCount) {
            double probeDistance = std::fabs(xTarget - xIn[nearestIndex + 1]);
            if (probeDistance >= bestDistance) break;
            isMovedUp = true;
            nearestIndex = nearestIndex + 1;
            bestDistance = probeDistance;
        }

        // If we never moved forward, scan backward
        if (!isMovedUp) {
            while (nearestIndex != 1) {
                double probeDistance = std::fabs(xTarget - xIn[nearestIndex - 1]);
                if (probeDistance >= bestDistance) break;
                nearestIndex = nearestIndex - 1;
                bestDistance = probeDistance;
            }
        }

        // Have found nearest xIn
        {
            int lowerIndex = nearestIndex - 1;
            int upperIndex = nearestIndex + 1;
            fInterp = fIn[nearestIndex];
            work1(1) = xTarget - xIn[nearestIndex];
            isMovedUp = work1(1) < 0.0;
            work2(1) = fInterp;

            // Do up to `degree` iterations.
            if (degree > 0) {
            for (int degreeIndex = 1; degreeIndex <= degree; degreeIndex++) {
                int i;

                // Bias choice of 3rd point towards center
                if (degreeIndex == 2) isMovedUp = xTarget < xMid;

                // Pick next point: alternate upper/lower based on isMovedUp,
                // fall back to the other side when one runs out of bounds.
                while (true) {
                    if (!isMovedUp) {
                        i = upperIndex;
                        upperIndex = upperIndex + 1;
                        isMovedUp = true;
                        if (i <= inputCount) break;
                    } else {
                        i = lowerIndex;
                        lowerIndex = lowerIndex - 1;
                        isMovedUp = false;
                        if (i >= 1) break;
                    }
                }

                {
                    double del = xTarget - xIn[i];
                    fInterp = fIn[i];
                    for (int priorSlotIndex = 1; priorSlotIndex <= degreeIndex; priorSlotIndex++) {
                        fInterp = (work1(priorSlotIndex) * fInterp - del * work2(priorSlotIndex)) / (work1(priorSlotIndex) - del);
                    }

                    work1(degreeIndex + 1) = del;
                    work2(degreeIndex + 1) = fInterp;
                }
            } // end degreeIndex loop
            }  // end if (degree > 0)
        }

        fOut[outputIndex] = fInterp;
    } // end outputIndex loop

    #undef work1
    #undef work2
    return;
}


// ============================================================================
//   POT(X) = -V * EXP( (R^(2P) - X^(2P)) / A^(2P) )
//
// 5/30/78 - made from FIXEDWOOD - S.P.
// ============================================================================

void GaussianLinkulePlugin::run(char8 /*alias*/, int* linkuleInts, int potType, int requestCode,
              int& callStatus, int /*L*/, double& /*J*/, double /*rStart*/,
              double stepSize, int nPts,
              double* array1, double* /*array2*/,
              Reaction& reaction)
{
    double undefValue = reaction.internalState.undefValue;

    callStatus = 0;

    // Get parameters based on potential type
    double v, a, r;
    readWoodsSaxonVAR(reaction, potType, v, a, r);

    if (r == undefValue || a == undefValue || v == undefValue) {
        printf("\n**** R, V, AND A MUST BE DEFINED.\n");
        callStatus = -1;
        return;
    }
    if (r <= 0.0 || a <= 0.0) {
        callStatus = -1;
        return;
    }

    switch (requestCode) {
    case 1:
        // Generate radius array for cubicSplineInterp
        linkuleInts[0] = 0;
        linkuleInts[1] = 0;
        return;

    case 2:
        // Print it out
        printf(" THE POTENTIAL IS A GENERALIZED GAUSSIAN:\n"
               " V =%13.5g AT R =%6.2f FM;   A =%8.3f FM\n",
               v, r, a);
        return;

    case 3: {
        int endIndex = (int)(std::sqrt(r * r + a * a * (46.0 + std::log(std::fabs(v)))) / stepSize);
        endIndex = std::min(nPts, endIndex);
        int startIndex = 1;
        double temp = r * r - a * a * (23.0 - std::log(std::fabs(v)));
        if (temp > 0.0) startIndex = (int)(std::sqrt(temp) / stepSize);
        startIndex = std::max(startIndex, 1);

        if (startIndex != 1) {
            for (int i = 1; i <= startIndex; i++) {
                array1[i] = -1.0e+10;
            }
        }

        double rValue = (startIndex - 1) * stepSize / a;
        for (int i = startIndex; i <= endIndex; i++) {
            double rOverA = r / a;
            array1[i] = -v * std::exp(rOverA * rOverA - rValue * rValue);
            rValue = rValue + (stepSize / a);
        }

        if (endIndex >= nPts) return;
        endIndex = endIndex + 1;
        for (int i = endIndex; i <= nPts; i++) {
            array1[i] = 0.0;
        }
        return;
    }

    default:
        return;
    }
}


// ============================================================================
// V(R=PARAM1) = VFIX, with shape from R, A, POWER.
//
// 4/6/78 - based on SHAPE - S.P.
// ============================================================================

void FixedWoodsSaxonLinkulePlugin::run(char8 /*alias*/, int* linkuleInts, int potType, int requestCode,
              int& callStatus, int /*L*/, double& /*J*/,
              double rStart, double stepSize, int nPts,
              double* array1, double* /*array2*/,
              Reaction& reaction)
{
    double undefValue = reaction.internalState.undefValue;

    callStatus = 0;

    double rFix = (potType == 1) ? reaction.linkuleParams.PARAM_arr[0]
                                 : reaction.linkuleParams.PARAM_arr[1];
    double vFix = (&reaction.distortedWave.scatteringSolver.potentialWork.tvReal)[potType - 1];
    double a = (&reaction.distortedWave.scatteringSolver.potentialWork.tvReal)[3 + potType];
    double r = (potType == 2) ? reaction.opticalPotentialParams.rI : reaction.integrationGrid.R;

    if (r == undefValue || a == undefValue || vFix == undefValue || rFix == undefValue) {
        reportRvaParam1Undefined(callStatus);
        return;
    }
    if (r <= 0.0 || a <= 0.0 || rFix <= 0.0) {
        callStatus = -1;
        return;
    }

    double v = vFix * (1.0 + std::exp((rFix - r) / a));   // pow(.., 1.0) folded out

    if (requestCode == 1) {
        linkuleInts[0] = 0;
        linkuleInts[1] = 0;
        return;
    }

    if (requestCode == 2) {
        printf(" THE POTENTIAL IS A WOODS-SAXON:%15.5g%15.5g%15.5g\n"
               " WITH V =%15.5g     AT R =%15.5g\n",
               v, r, a, vFix, rFix);
        return;
    }

    if (requestCode == 3) {
        // resize() zeroes the buffer so ADD-semantics fillWoodsSaxon == old SET;
        // array1 isn't pre-zeroed, so fill a fresh grid then copy into array1[1..nPts].
        OpticalPotential pot;
        pot.resize(nPts, rStart, stepSize);
        pot.fillWoodsSaxon(v, r, a);
        std::copy(pot.values.begin(), pot.values.end(), array1 + 1);
        return;
    }
}


// ============================================================================
//   V = (1 + PARAM1*J + PARAM2*J^2) * WS(V, R, A)
//
// 11/16/81 - JDEPEN based on PARITWOOD - S.P.
// ============================================================================

void JDependentWoodsSaxonLinkulePlugin::run(char8 /*alias*/, int* linkuleInts, int potType, int requestCode,
            int& callStatus, int /*L*/, double& /*J*/,
            double rStart, double stepSize, int nPts,
            double* array1, double* array2,
            Reaction& reaction)
{
    double undefValue = reaction.internalState.undefValue;

    callStatus = 0;

    double param1 = reaction.linkuleParams.PARAM_arr[0];
    double param2 = reaction.linkuleParams.PARAM_arr[1];
    if (param2 == undefValue) param2 = 0.0;
    double v, a, r;
    readWoodsSaxonVAR(reaction, potType, v, a, r);

    if (r == undefValue || a == undefValue || v == undefValue || param1 == undefValue) {
        reportRvaParam1Undefined(callStatus);
        return;
    }
    if (bailIfBadRadiusOrDiffuse(r, a, callStatus)) return;

    if (requestCode != 1) {
        if (requestCode == 2) {
            printf(" THE POTENTIAL IS A J-DEPENDANT WOODS-SAXON:%15.5g%15.5g%15.5g\n"
                   " WITH THE J-DEPENDANT DEPTH factor  1 + J *%13.5g + J**2 *%13.5g\n",
                   v, r, a, param1, param2);
            return;
        }
        if (requestCode == 3) {
            // Set up the two pieces
            setupWoodsSaxonPiece(reaction, linkuleInts, nPts, rStart, stepSize, array1, v, r, a);
            return;
        }
        if (requestCode == 4) {
            // Add in the J-dependent part
            double* workArrayPointer = reaction.linkuleData.workarrs[linkuleInts[0] - 1].data();  // 0-based
            double halfJ = 0.5 * reaction.angMom.J;
            double fac = array2[1] * halfJ * (param1 + halfJ * param2);
            int n1 = (int)(rStart / stepSize + 0.5);
            double* shiftedWorkPointer = workArrayPointer + n1;
            for (int i = 1; i <= nPts; i++) {
                array1[i] = array1[i] + fac * shiftedWorkPointer[i - 1];
            }
            return;
        }
        return;
    }

    // Generate work array to save the W.S.
    // linkuleInts[0] is a 1-based index
    // into LinkuleData.workarrs; the NAMLOC-tagged name was dead (no test
    // looks up the fitter iD via REDEF).
    {
        linkuleInts[0] = allocateWorkSlot(reaction, nPts);
        linkuleInts[1] = 0;
    }
    // This is an L-dependent potential
    callStatus = +1;
    return;
}


// ============================================================================
//   V = F(J) * WS(V, R, A)
//   F(J) = 1 / (1 + EXP((J - PARM100)/PARM101))
//
// 5/11/82 - JDEPENWS made from JDEPEN - S.P.
// ============================================================================

void JDependentWoodsSaxonFermiLinkulePlugin::run(char8 /*alias*/, int* linkuleInts, int potType, int requestCode,
              int& callStatus, int /*L*/, double& /*J*/,
              double rStart, double stepSize, int nPts,
              double* array1, double* array2,
              Reaction& reaction)
{
    double undefValue = reaction.internalState.undefValue;

    callStatus = 0;

    // param(10) and param(11) — these are params(reaction)[9] and params(reaction)[10]
    double param10 = params(reaction)[9];
    double param11 = params(reaction)[10];
    double v, a, r;
    readWoodsSaxonVAR(reaction, potType, v, a, r);

    if (r == undefValue || a == undefValue || v == undefValue ||
        param10 == undefValue || param11 == undefValue) {
        printf("\n**** R, V, A, PARAM10 & PARAM11 MUST BE DEFINED.\n");
        callStatus = -1;
        return;
    }
    if (bailIfBadRadiusOrDiffuse(r, a, callStatus)) return;

    if (requestCode == 1) {
        // linkuleInts[0] holds the 1-based index for the requestCode==3/4 readers.
        linkuleInts[0] = allocateWorkSlot(reaction, nPts);
        linkuleInts[1] = 0;
        callStatus = +1;
        return;
    }
    if (requestCode == 2) {
        printf(" THE POTENTIAL IS A J-DEPENDANT WOODS-SAXON:%15.5g%15.5g%15.5g\n"
               " WITH THE J-DEPENDANT DEPTH factor  1 / ( 1 + EXP( ( J - %12.5g) / %12.5g) )\n",
               v, r, a, param10, param11);
        return;
    }
    if (requestCode == 3) {
        setupWoodsSaxonPiece(reaction, linkuleInts, nPts, rStart, stepSize, array1, v, r, a);
        return;
    }
    if (requestCode == 4) {
        // Add in the J-dependent part
        double* workArrayPointer = reaction.linkuleData.workarrs[linkuleInts[0] - 1].data();  // 0-based
        double halfJ = 0.5 * reaction.angMom.J;
        double fac = 1.0;
        double x = (halfJ - param10) / param11;
        if (x >= -30.0) {
            fac = 0.0;
            if (x <= 30.0) {
                fac = 1.0 / (1.0 + std::exp(x));
            }
        }
        fac = array2[1] * (fac - 1.0);
        int n1 = (int)(rStart / stepSize + 0.5);
        double* shiftedWorkPointer = workArrayPointer + n1;
        for (int i = 1; i <= nPts; i++) {
            array1[i] = array1[i] + fac * shiftedWorkPointer[i - 1];
        }
        return;
    }
}


// ============================================================================
//   V = (1 + PARAM1*(-1)^L) * WS(V, R, A)
//
// 6/13/78 - based on FIXEDWOOD - S.P.
// ============================================================================

void ParityWoodsSaxonLinkulePlugin::run(char8 /*alias*/, int* linkuleInts, int potType, int requestCode,
              int& callStatus, int L, double& /*J*/, double rStart, double stepSize, int nPts,
              double* array1, double* array2,
              Reaction& reaction)
{
    double undefValue = reaction.internalState.undefValue;

    callStatus = 0;

    double parityCoefficient = (potType == 1) ? reaction.linkuleParams.PARAM_arr[0]
                                              : reaction.linkuleParams.PARAM_arr[1];
    double v, a, r;
    readWoodsSaxonVAR(reaction, potType, v, a, r);

    if (r == undefValue || a == undefValue || v == undefValue || parityCoefficient == undefValue) {
        reportRvaParam1Undefined(callStatus);
        return;
    }
    if (bailIfBadRadiusOrDiffuse(r, a, callStatus)) return;

    if (requestCode != 1) {
        if (requestCode == 2) {
            printf(" THE POTENTIAL IS A PARITY WOODS-SAXON:%15.5g%15.5g%15.5g\n"
                   " WITH THE PARITY-DEPENDANT DEPTH factor  1 + (-PARAM1)**L:  PARAM1 =%15.5g\n",
                   v, r, a, parityCoefficient);
            return;
        }
        if (requestCode == 3) {
            setupWoodsSaxonPiece(reaction, linkuleInts, nPts, rStart, stepSize, array1, v, r, a);
            return;
        }
        if (requestCode == 4) {
            // Add in the L-dependent part
            double* workArrayPointer = reaction.linkuleData.workarrs[linkuleInts[0] - 1].data();  // 0-based
            double fac = array2[1] * parityCoefficient;
            if ((((L) >> (0)) & 1)) fac = -fac;
            for (int i = 1; i <= nPts; i++) {
                array1[i] = array1[i] + fac * workArrayPointer[i - 1];
            }
            return;
        }
        return;
    }

        // linkuleInts[0] holds the 1-based index for the requestCode==3/4 readers.
    {
        linkuleInts[0] = allocateWorkSlot(reaction, nPts);
        linkuleInts[1] = 0;
    }
    callStatus = +1;
    return;
}


// ============================================================================
// Shape is read from named objects (REALSHAP, IMAGSHAP, etc.)
//
// 12/8/77 - first version - S.P.
// ============================================================================

void ShapeLinkulePlugin::run(char8 /*alias*/, int* linkuleInts, int potType, int requestCode,
           int& callStatus, int /*L*/, double& /*J*/,
           double rStart, double stepSize, int nPts,
           double* array1, double* /*array2*/,
           Reaction& reaction)
{
    static const char names[6][9] = {
        "", "REALSHAP", "IMAGSHAP", "REALSOSH", "IMAGSOSH", "SHAPE   "
    };
    static const char scalNames[6][9] = {
        "", "REALSCAL", "IMAGSCAL", "REALSOSC", "IMAGSOSC", "SHAPESCA"
    };

    callStatus = 0;

    double cons = (&reaction.distortedWave.scatteringSolver.potentialWork.tvReal)[potType - 1];
    if (cons == reaction.internalState.undefValue) {
        printf("\n**** THE WELL DEPTH PARAMETER (V, vI, vSo, VSOI) MUST BE DEFINED FOR THE SHAPE LINKULE.\n");
        callStatus = -1;
        return;
    }
    cons = -cons;

    if (requestCode == 2) {
        printf(" THE POTENTIAL SHAPE STORED IN OBJECT %.8s IS BEING MULTIPLIED BY%15.5g\n",
               names[linkuleInts[1]], cons);
        return;
    }
    if (requestCode == 3) {
        // Read the stored shape from the workarrs slot allocated in requestCode==1.
        const std::vector<double>& shape =
            reaction.linkuleData.workarrs[linkuleInts[0] - 1];
        for (int i = 1; i <= nPts; i++) {
            array1[i] = cons * shape[i - 1];
        }
        return;
    }
    if (requestCode != 1) return;

    // requestCode == 1: locate the user-supplied shape array. Try names[potType]
    // first; fall back to names[5]="SHAPE   ".
    int ii = potType;
    const std::vector<double>* src = reaction.named.find(names[ii]);
    if (src == nullptr) {
        ii = 5;
        src = reaction.named.find(names[ii]);
        if (src == nullptr) {
            printf("\n**** AN OBJECT WITH THE NAME %.8s OR SHAPE MUST BE DEFINED TO USE THE SHAPE LINKULE.\n",
                   names[potType]);
            callStatus = -1;
            return;
        }
    }

    int inputCount = (int)src->size();
    double rMax = rStart + (nPts - 1) * stepSize;
    const std::vector<double>* scal = reaction.named.find(scalNames[ii]);

    bool skipInterp = false;
    bool isFullRArray = false;
    double r1 = 0.0, r2 = 0.0;

    if (scal == nullptr) {
        if (nPts == inputCount) {
            printf(" %.8s IS BEING ASSUMED TO BE DEFINED FOR%15.5g =< R =<%15.5g WITH STEPSIZE =%15.5g\n",
                   names[ii], rStart, rMax, stepSize);
            skipInterp = true;
        } else {
            printf("\n**** ASYMPTOPIA AND STEPSIZE =%15.5g%15.5g AND REQUIRE%6d POINTS.\n",
                   rMax, stepSize, nPts);
            printf("      HOWEVER, %.8s HAS%6d POINTS AND %.8s IS NOT DEFINED.\n",
                   names[ii], inputCount, scalNames[ii]);
            callStatus = -1;
            return;
        }
    } else {
        int scalSize = (int)scal->size();
        if (scalSize == 2) {
            r1 = (*scal)[0];
            r2 = (*scal)[1];
            if (r1 == rStart && inputCount == nPts &&
                std::fabs(r2 - rMax) / stepSize < 1.0e-5) {
                skipInterp = true;
            }
        } else if (scalSize == inputCount) {
            isFullRArray = true;
            r1 = (*scal)[0];
            r2 = (*scal)[inputCount - 1];
        } else {
            printf("\n**** SCALLING ARRAY %.8s SHOULD HAVE 2 OR%4d ELEMENTS, BUT IT HAS%6d ELEMENTS.\n",
                   scalNames[ii], inputCount, scalSize);
            callStatus = -1;
            return;
        }
    }

    // Validation passed — allocate the output workarrs slot (nPts long).
    int slotIndex = allocateWorkSlot(reaction, nPts);
    std::vector<double>& dst = reaction.linkuleData.workarrs[slotIndex - 1];

    if (skipInterp) {
        for (int i = 0; i < nPts; i++) dst[i] = (*src)[i];
    } else {
        // aitken interpolation onto the (rStart, stepSize, nPts) grid.
        int workSize = inputCount + nPts;
        std::vector<double> xWorkVector(workSize + 1, 0.0);

        double* yInPointer  = const_cast<double*>(src->data()) - 1;  // 1-based
        double* yOutPointer = dst.data() - 1;                        // 1-based

        double* xInPointer  = xWorkVector.data() - 1;
        double* xOutPointer = xWorkVector.data() + inputCount - 1;
        if (isFullRArray) xInPointer = const_cast<double*>(scal->data()) - 1;
        if (!isFullRArray) {
            double rr = r1;
            double stepIn = (r2 - r1) / (inputCount - 1);
            for (int n = 1; n <= inputCount; n++) {
                xInPointer[n] = rr;
                rr = rr + stepIn;
            }
        }

        // Get extrapolating function for points past r2.
        double vStep = 0.0;
        double aVal = 0.0;
        double ynM1 = yInPointer[inputCount - 1];
        double yn   = yInPointer[inputCount];
        if (yn * ynM1 > 0.0) {
            aVal = yn;
            double bVal = std::log(ynM1 / aVal) / (r2 - xInPointer[inputCount - 1]);
            vStep = std::exp(-bVal * stepSize);
        }

        // Output X's: clamp below r1 to val0, extrapolate above r2.
        double rr = rStart;
        int startIndex = 1;
        int lastIndex = 1;
        double val0 = yInPointer[1];
        for (int n = 1; n <= nPts; n++) {
            if (rr >= r1) {
                if (rr <= r2) {
                    xOutPointer[n] = rr;
                    lastIndex = n;
                } else {
                    aVal = aVal * vStep;
                    if (std::fabs(aVal) < Constants::smlNum) aVal = 0.0;
                    yOutPointer[n] = aVal;
                }
            } else {
                yOutPointer[n] = val0;
                startIndex = n + 1;
            }
            rr = rr + stepSize;
        }
        int outputCount = lastIndex + 1 - startIndex;

        double work[20];
        aitken(5, inputCount, xInPointer, (yInPointer + 1),
               outputCount, (xOutPointer + startIndex - 1), (yOutPointer + startIndex),
               work);
    }

    linkuleInts[0] = slotIndex;
    linkuleInts[1] = ii;
}


// ============================================================================
// Uses param pairs (R1,V1)...(R10,V10).
//
// 4/6/78 - based on SHAPE - S.P.
// ============================================================================

void SplineLinkulePlugin::run(char8 /*alias*/, int* linkuleInts, int /*potType*/, int requestCode,
                    int& callStatus, int /*L*/, double& /*J*/,
                    double rStart, double stepSize, int nPts,
                    double* array1, double* /*array2*/,
                    Reaction& reaction)
{
    double undefValue = reaction.internalState.undefValue;

    // Local arrays for sorting points
    double rPoints[21], vPoints[21], bCoef[21], cCoef[21], dCoef[21];

    callStatus = 0;

    // Collect & order the defined R/V point pairs (shared with lagrange()).
    int pairCount = collectOrderedPoints(params(reaction), rPoints, vPoints,
                                         undefValue, "SPLINE", callStatus);
    if (pairCount < 0) return;

    switch (requestCode) {
    case 1:
        buildRadialGridSlot(reaction, linkuleInts, nPts, rStart, stepSize);
        return;

    case 2:
        printf(" THE POTENTIAL IS DEFINED BY THE POINTS:\n");
        printPointTable(rPoints, vPoints, pairCount);
        return;

    case 3: {
        double* workArrayPointer = reaction.linkuleData.workarrs[linkuleInts[0] - 1].data() - 2;  // workarrs (1-based index)

        // Setup output X's
        int startIndex, lastIndex;
        int outputCount = setupOutputWindow(array1, rPoints, vPoints, pairCount,
                                            nPts, rStart, stepSize, startIndex, lastIndex);

        // Interpolate the logs
        for (int i = 1; i <= pairCount; i++) {
            vPoints[i] = std::log(vPoints[i]);
        }

        naturalCubicSpline(pairCount, rPoints, vPoints, bCoef, cCoef, dCoef);
        cubicSplineInterp(pairCount, rPoints, vPoints, bCoef, cCoef, dCoef,
               outputCount, (workArrayPointer + 1 + startIndex - 1), &array1[startIndex]);

        for (int i = startIndex; i <= lastIndex; i++) {
            array1[i] = -std::exp(array1[i]);
        }

        // Fill in the end with exponential decay
        double term = array1[lastIndex];
        lastIndex = lastIndex + 1;
        if (lastIndex > nPts) return;
        double step = std::exp(bCoef[pairCount] * stepSize);
        term = term * step;
        for (int i = lastIndex; i <= nPts; i++) {
            array1[i] = term;
            if (term > -1.0e-30) return;
            term = term * step;
        }
        return;
    }

    default:
        return;
    }
}


// ============================================================================
// Uses param pairs and aitken interpolation.
//
// 4/10/78 - based on SPLINE - S.P.
// ============================================================================

void LagrangeLinkulePlugin::run(char8 /*alias*/, int* linkuleInts, int /*potType*/, int requestCode,
              int& callStatus, int /*L*/, double& /*J*/,
              double rStart, double stepSize, int nPts,
              double* array1, double* /*array2*/,
              Reaction& reaction)
{
    double undefValue = reaction.internalState.undefValue;
    double rPoints[21], vPoints[21], aitkenWork[31];

    callStatus = 0;

    // Interpolation order — was reaction.MAXFUN, set by DEFALT to 50 with
    int maxOrder = 50;

    // Collect & order the defined R/V point pairs (shared with splineLinkule()).
    int pairCount = collectOrderedPoints(params(reaction), rPoints, vPoints,
                                         undefValue, "LAGRANGE", callStatus);
    if (pairCount < 0) return;

    if (requestCode == 1) {
        buildRadialGridSlot(reaction, linkuleInts, nPts, rStart, stepSize);
        return;
    }
    if (requestCode == 2) {
        printf(" THE POTENTIAL IS DEFINED BY LAGRANGE INTERPOLATION OF ORDER%2d BASED ON THE POINTS:\n",
               maxOrder);
        printPointTable(rPoints, vPoints, pairCount);
        return;
    }
    if (requestCode != 3) return;

    {
        double* workArrayPointer = reaction.linkuleData.workarrs[linkuleInts[0] - 1].data() - 2;  // workarrs (1-based index)

        int startIndex, lastIndex;
        int outputCount = setupOutputWindow(array1, rPoints, vPoints, pairCount,
                                            nPts, rStart, stepSize, startIndex, lastIndex);

        // Does the potential change sign
        int useOrder = maxOrder;
        double vSign = -std::copysign(1.0, vPoints[pairCount]);
        int lStartIndex = 1;

        int signChangeIndex = 0;
        for (int ii = 2; ii <= pairCount; ii++) {
            int i = pairCount + 2 - ii;
            if (vSign * vPoints[i - 1] < 0.0) {
                signChangeIndex = i;
                break;
            }
        }

        if (signChangeIndex > 0) {
            // Yes, use straight interpolation up to the last change
            lStartIndex = signChangeIndex;
            int startIndex2 = (int)((rPoints[lStartIndex] - rStart) / stepSize + 2);
            aitken(std::min(maxOrder, lStartIndex - 1), lStartIndex, rPoints, vPoints,
                   startIndex2 - startIndex, (workArrayPointer + 1 + startIndex - 1), &array1[startIndex],
                   aitkenWork);
            for (int i = startIndex; i <= startIndex2; i++) {
                array1[i] = -array1[i];
            }

            useOrder = std::min(maxOrder, pairCount - lStartIndex);
            if (useOrder <= 0) return;
            startIndex = startIndex2;
            outputCount = lastIndex + 1 - startIndex2;
        }

        // Interpolate the logs
        for (int i = lStartIndex; i <= pairCount; i++) {
            vPoints[i] = std::log(std::max(std::fabs(vPoints[i]), 1.0e-15));
        }

        aitken(useOrder, pairCount - lStartIndex + 1,
               &rPoints[lStartIndex], &vPoints[lStartIndex],
               outputCount, (workArrayPointer + 1 + startIndex - 1), &array1[startIndex],
               aitkenWork);

        for (int i = startIndex; i <= lastIndex; i++) {
            double x = array1[i];
            x = std::max(x, -69.0);
            x = std::min(x, 69.0);
            array1[i] = vSign * std::exp(x);
        }

        // Fill in the end with exponential decay
        double term = array1[lastIndex];
        lastIndex = lastIndex + 1;
        if (lastIndex > nPts) return;
        double step = term / array1[lastIndex - 2];
        term = term * step;
        if (1.0 - step < 1.0e-5) return;
        for (int i = lastIndex; i <= nPts; i++) {
            array1[i] = term;
            if (term > -1.0e-30) return;
            term = term * step;
        }
    }
    return;
}


// They printed "NOT YET OPERATIONAL" and returned callStatus=-1; no test
// exercises them. LINKUL's default branch now handles their linkuleIndex values
// (1/5/10/11) with "LINKULE NOT AVAILABLE" + FSTOP.

std::unique_ptr<LinkulePlugin> makeFixedWoodsSaxonPlugin()        { return std::make_unique<FixedWoodsSaxonLinkulePlugin>(); }          // 2
std::unique_ptr<LinkulePlugin> makeGaussianPlugin()              { return std::make_unique<GaussianLinkulePlugin>(); }                // 3
std::unique_ptr<LinkulePlugin> makeLagrangePlugin()             { return std::make_unique<LagrangeLinkulePlugin>(); }               // 4
std::unique_ptr<LinkulePlugin> makeShapePlugin()               { return std::make_unique<ShapeLinkulePlugin>(); }                  // 8
std::unique_ptr<LinkulePlugin> makeSplinePlugin()             { return std::make_unique<SplineLinkulePlugin>(); }                 // 9
std::unique_ptr<LinkulePlugin> makeJDependentWoodsSaxonPlugin()      { return std::make_unique<JDependentWoodsSaxonLinkulePlugin>(); }     // 12
std::unique_ptr<LinkulePlugin> makeJDependentWoodsSaxonFermiPlugin() { return std::make_unique<JDependentWoodsSaxonFermiLinkulePlugin>(); }// 13
std::unique_ptr<LinkulePlugin> makeParityWoodsSaxonPlugin()         { return std::make_unique<ParityWoodsSaxonLinkulePlugin>(); }        // 15
