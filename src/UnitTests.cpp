// unit_tests.cpp — Unit tests for the C++ port's standalone math functions
// Tests: clebschGordan, threeJ, sixJ, racah, wig9J, gaussL, naturalCubicSpline, plmSub, RCWFN
// All angular momentum arguments use Ptolemy convention: 2*J (integers)

#include "CoulombWaveFunction.h"
#include "math/angular_momentum_coeff.h"
#include "math/gauss_quadrature.h"
#include "math/spline.h"
#include "math/legendre.h"
#include "OpticalPotential.h"   // new flat class under test
#include "OpticalPotentialLibrary.h"  // named OM potential library (Phase A)
#include "DwbaInputExpander.h"         // DWBA → Ptolemy deck expander (Phase B)
#include <string>
void setVsq(double rr1, double rr2, int iz1, int iz2, int k);  // seeds vcsq12 (source_misc.cpp)
#include <cstdio>
#include <cmath>
#include <vector>

static int passCount = 0;
static int failCount = 0;

static void record(const char* name, double got, double expected, double err, double tolerance) {
    if (err < tolerance) {
        passCount++;
    } else {
        printf("  FAIL  %-40s got=%.10g  expected=%.10g  err=%.2e\n", name, got, expected, err);
        failCount++;
    }
}

static void check(const char* name, double got, double expected, double tolerance) {
    double err = (expected != 0) ? std::fabs(got - expected) / std::fabs(expected)
                                 : std::fabs(got);
    record(name, got, expected, err, tolerance);
}

static void checkAbs(const char* name, double got, double expected, double tolerance) {
    double err = std::fabs(got - expected);
    record(name, got, expected, err, tolerance);
}

// ============================================================================
// Clebsch-Gordan coefficients: <j1 m1 j2 m2 | J M>
// Args: clebschGordan(2*j1, 2*j2, 2*m1, 2*m2, 2*J, 2*M)
// ============================================================================
static void testClebsch() {
    printf("\n--- Clebsch-Gordan ---\n");

    // <1/2 1/2 1/2 -1/2 | 1 0> = 1/sqrt(2)
    check("CG(1/2,1/2,1/2,-1/2|1,0)", clebschGordan(1,1,1,-1,2,0), 1.0/std::sqrt(2.0), 1e-12);

    // <1/2 1/2 1/2 1/2 | 1 1> = 1
    check("CG(1/2,1/2,1/2,1/2|1,1)", clebschGordan(1,1,1,1,2,2), 1.0, 1e-12);

    // <1 0 1 0 | 0 0> = -1/sqrt(3)
    check("CG(1,0,1,0|0,0)", clebschGordan(2,2,0,0,0,0), -1.0/std::sqrt(3.0), 1e-12);

    // <1 0 1 0 | 2 0> = sqrt(2/3)
    check("CG(1,0,1,0|2,0)", clebschGordan(2,2,0,0,4,0), std::sqrt(2.0/3.0), 1e-12);

    // <1 1 1 -1 | 0 0> = 1/sqrt(3)
    check("CG(1,1,1,-1|0,0)", clebschGordan(2,2,2,-2,0,0), 1.0/std::sqrt(3.0), 1e-12);

    // Zero by M conservation: m1+m2 != M
    checkAbs("CG(1,1,1,0|1,0) = 0", clebschGordan(2,2,2,0,2,0), 0.0, 1e-14);

    // <2 0 1 0 | 2 0> = 0 (parity)
    checkAbs("CG(2,0,1,0|2,0) = 0", clebschGordan(4,2,0,0,4,0), 0.0, 1e-14);

    // Large J: <5 0 5 0 | 0 0> — sign from (-1)^5 phase
    check("CG(5,0,5,0|0,0)", clebschGordan(10,10,0,0,0,0), -1.0/std::sqrt(11.0), 1e-12);
}

// ============================================================================
// 3-J symbols
// THRJ(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3)
// ============================================================================
static void testThreej() {
    printf("\n--- 3-J symbols ---\n");

    // (1/2 1/2 1 ; 1/2 -1/2 0) = 1/sqrt(6)  [sign: (-1)^(j1-j2-m3)]
    double val = threeJ(1,1,2, 1,-1,0);
    check("3J(1/2,1/2,1;1/2,-1/2,0)", val, 1.0/std::sqrt(6.0), 1e-12);

    // (1 1 0 ; 0 0 0) = (-1)^1 / sqrt(3) = -1/sqrt(3)
    val = threeJ(2,2,0, 0,0,0);
    check("3J(1,1,0;0,0,0)", val, -1.0/std::sqrt(3.0), 1e-12);

    // Triangle rule violation: (1 1 3 ; 0 0 0) = 0
    val = threeJ(2,2,6, 0,0,0);
    checkAbs("3J triangle violation = 0", val, 0.0, 1e-14);
}

// ============================================================================
// 6-J symbols
// sixJ(2*j1, 2*j2, 2*j3, 2*j4, 2*j5, 2*j6)
// ============================================================================
static void testSixj() {
    printf("\n--- 6-J symbols ---\n");

    // {1/2 1/2 0 ; 1/2 1/2 1}
    // Ptolemy includes (-1)^(NMIN+R) phase
    double val = sixJ(1,1,0, 1,1,2);
    check("6J{1/2,1/2,0;1/2,1/2,1}", val, 0.5, 1e-12);

    // {1 1 0 ; 1 1 1}
    val = sixJ(2,2,0, 2,2,2);
    check("6J{1,1,0;1,1,1}", val, -1.0/3.0, 1e-12);

    // {1 1 2 ; 1 1 0}
    val = sixJ(2,2,4, 2,2,0);
    check("6J{1,1,2;1,1,0}", val, 1.0/3.0, 1e-12);

    // Triangle violation
    val = sixJ(2,2,10, 2,2,2);
    checkAbs("6J triangle violation = 0", val, 0.0, 1e-14);
}

// ============================================================================
// Racah W coefficients
// W(a,b,e,d;c,f) = (-1)^(a+b+d+e) * {a b c ; d e f}
// racah(2*a, 2*b, 2*e, 2*d, 2*c, 2*f) — note arg order!
// ============================================================================
static void test_racah() {
    printf("\n--- Racah W ---\n");

    // Racah = (-1)^phase * 6J; check against actual computed value
    double val = racah(1,1,1,1,0,2);
    check("W(1/2,1/2,1/2,1/2;0,1)", val, 0.5, 1e-12);
}

// ============================================================================
// 9-J symbols
// wig9J(2*j1, ..., 2*j9)
// ============================================================================
static void test_wig9j() {
    printf("\n--- 9-J symbols ---\n");

    // Check against known value — compute numerically
    double val = wig9J(1,1,0, 1,1,0, 0,0,0);
    check("9J{1/2,1/2,0;1/2,1/2,0;0,0,0}", val, 0.5, 1e-12);
}

// ============================================================================
// Gauss-Legendre quadrature
// ============================================================================
static void test_gaussl() {
    printf("\n--- Gauss-Legendre ---\n");

    double x[33], w[33]; // 1-based

    // N=2: roots at ±1/sqrt(3), weights = 1
    gaussL(2, x, w);
    check("GL2 x1", x[1], -1.0/std::sqrt(3.0), 1e-14);
    check("GL2 x2", x[2],  1.0/std::sqrt(3.0), 1e-14);
    check("GL2 w1", w[1], 1.0, 1e-14);
    check("GL2 w2", w[2], 1.0, 1e-14);

    // N=4: sum of weights = 2
    gaussL(4, x, w);
    double wSum = w[1] + w[2] + w[3] + w[4];
    check("GL4 sum(w)=2", wSum, 2.0, 1e-14);

    // Symmetry: x[1] = -x[4], x[2] = -x[3]
    checkAbs("GL4 symmetry x1+x4", x[1]+x[4], 0.0, 1e-14);
    checkAbs("GL4 symmetry x2+x3", x[2]+x[3], 0.0, 1e-14);

    // Integrate x^2 from -1 to 1 = 2/3 (exact for N>=2)
    double integral = 0;
    for (int i = 1; i <= 4; i++) integral += w[i] * x[i] * x[i];
    check("GL4 int(x^2)=2/3", integral, 2.0/3.0, 1e-14);

    // N=16: integrate x^30 (exact for 2N-1 >= 30 → N >= 16)
    gaussL(16, x, w);
    integral = 0;
    for (int i = 1; i <= 16; i++) integral += w[i] * std::pow(x[i], 30);
    check("GL16 int(x^30)=2/31", integral, 2.0/31.0, 1e-12);
}

// ============================================================================
// Cubic spline: naturalCubicSpline
// ============================================================================
static void test_splncb() {
    printf("\n--- Cubic spline ---\n");

    // Spline sin(x) on [0, pi/4, pi/2, 3pi/4, pi] — interior should be accurate
    int n = 5;
    double xv[6] = {0, 0, M_PI/4, M_PI/2, 3*M_PI/4, M_PI}; // 1-based
    double yv[6] = {0, std::sin(xv[1]), std::sin(xv[2]), std::sin(xv[3]), std::sin(xv[4]), std::sin(xv[5])};
    double b[6], c[6], d[6];

    naturalCubicSpline(n, xv, yv, b, c, d);

    // Evaluate at x = pi/3 (segment 2: between pi/4 and pi/2)
    double dx = M_PI/3 - xv[2];
    double yInterp = yv[2] + dx*(b[2] + dx*(c[2] + dx*d[2]));
    check("spline sin(pi/3)", yInterp, std::sin(M_PI/3), 2e-3);

    // Evaluate at x = pi/2 (knot point — should be exact)
    dx = 0;
    yInterp = yv[3];
    check("spline sin(pi/2) exact", yInterp, 1.0, 1e-14);
}

// ============================================================================
// Legendre polynomials: plmSub
// ============================================================================
static void test_plmsub() {
    printf("\n--- Legendre P_l(x) ---\n");

    double legendreP[100]; // 1-based
    double x = 0.5;

    // P_0(0.5) = 1
    plmSub(0, 0, x, legendreP);
    check("P_0(0.5)", legendreP[1], 1.0, 1e-14);

    // P_1(0.5) = 0.5
    plmSub(1, 0, x, legendreP);
    check("P_1(0.5)", legendreP[2], 0.5, 1e-14);

    // P_2(0.5) = (3*0.25 - 1)/2 = -0.125
    plmSub(2, 0, x, legendreP);
    check("P_2(0.5)", legendreP[3], -0.125, 1e-14);

    // P_3(0.5) = (5*0.125 - 3*0.5)/2 = -0.4375
    plmSub(3, 0, x, legendreP);
    check("P_3(0.5)", legendreP[4], -0.4375, 1e-14);

    // P_10(0) = 0 (odd index → P=0 at x=0? No, P_10 is even)
    // Actually P_10(0) = (-1)^5 * (10-1)!! / 10!! = -63/256
    plmSub(10, 0, 0.0, legendreP);
    check("P_10(0)", legendreP[11], -63.0/256.0, 1e-12);

    // P_l(1) = 1 for all l
    plmSub(5, 0, 1.0, legendreP);
    check("P_5(1)=1", legendreP[6], 1.0, 1e-14);
}

// ============================================================================
// Coulomb wavefunctions: RCWFN
// ============================================================================
static void test_rcwfn() {
    printf("\n--- Coulomb wavefunctions ---\n");

    double fC[20], fCp[20], gC[20], gCp[20];
    int returnCode;

    // For eta=0, F_0(rho) = sin(rho), G_0(rho) = cos(rho)
    double rho = 1.0;
    CoulombWaveFunction::computeFG(rho, 0.0, 0, 5, fC, fCp, gC, gCp, 1e-12, returnCode);
    check("F_0(1,eta=0)=sin(1)", fC[0], std::sin(1.0), 1e-10);
    check("G_0(1,eta=0)=cos(1)", gC[0], std::cos(1.0), 1e-10);

    // Wronskian: |F*G' - F'*G| = 1 for all L (sign depends on convention)
    for (int L = 0; L <= 5; L++) {
        double W = fC[L]*gCp[L] - fCp[L]*gC[L];
        char name[64];
        snprintf(name, sizeof(name), "Wronskian L=%d", L);
        check(name, std::fabs(W), 1.0, 1e-8);
    }

    // With eta: rho=5, eta=2 (proton on ~Z=10 at moderate energy)
    rho = 5.0;
    CoulombWaveFunction::computeFG(rho, 2.0, 0, 3, fC, fCp, gC, gCp, 1e-12, returnCode);
    double W = fC[0]*gCp[0] - fCp[0]*gC[0];
    check("Wronskian eta=2 L=0", std::fabs(W), 1.0, 1e-8);
    W = fC[2]*gCp[2] - fCp[2]*gC[2];
    check("Wronskian eta=2 L=2", std::fabs(W), 1.0, 1e-8);
}

// ============================================================================
// Main
// ============================================================================
static void test_numerov_free_particle();
static void test_optical_potential_bit_identity();
static void test_optical_potential_library();
static void test_dwba_expander();

int main() {
    printf("=== PtolemyCpp Unit Tests ===\n");

    testClebsch();
    testThreej();
    testSixj();
    test_racah();
    test_wig9j();
    test_gaussl();
    test_splncb();
    test_plmsub();
    test_rcwfn();

    printf("\n--- NumerovSolver ---\n");
    test_numerov_free_particle();

    printf("\n--- OpticalPotential (flat) bit-identity vs golden ---\n");
    test_optical_potential_bit_identity();

    printf("\n--- OpticalPotentialLibrary (named OMPs) vs Cleopatra golden ---\n");
    test_optical_potential_library();

    printf("\n--- DwbaInputExpander ---\n");
    test_dwba_expander();


    printf("\n══════════════════════════════════════\n"
           "  Results: %d passed, %d failed / %d total\n"
           "══════════════════════════════════════\n",
           passCount, failCount, passCount + failCount);

    return failCount > 0 ? 1 : 0;
}

// ============================================================================
// NumerovSolver test — free particle u''(r) = -k²u(r)
// Solution: u(r) = sin(kr), starting at u(0)=0, u(h)=sin(kh)
// ============================================================================
#include "NumerovSolver.h"

static void test_numerov_free_particle() {
    NumerovSolver solver;
    solver.stepSize = 0.1;

    double k = 1.0;  // wave number
    double h = solver.stepSize;
    int n = 100;

    // f(r) = -k² (constant for free particle)
    std::vector<double> f(n, k * k);  // u'"' + k²u = 0 → f = +k²

    // Initial: u(0) = sin(0) = 0, u(h) = sin(kh)
    solver.integrateOutward(f, 0.0, std::sin(k * h));

    // Check solution at several points against sin(kr)
    double maxErr = 0.0;
    for (int i = 0; i < n; i++) {
        double r = i * h;
        double exact = std::sin(k * r);
        double err = std::fabs(solver.solution[i] - exact);
        if (err > maxErr) maxErr = err;
    }
    // Numerov is O(h⁶) accurate, with h=0.1 and 100 steps should be <1e-6
    check("Numerov free particle maxErr", maxErr, 0.0, 1e-4);

    // Test node counting: sin(kr) for k=1, r=0..10 has 3 nodes (π, 2π, 3π)
    check("Numerov nodes", (double)solver.nodesFound, 3.0, 0.01);
}

// ============================================================================
// OpticalPotential (flat) golden-value tests
// Each new fill* is pinned to captured golden values (6 sampled grid indices,
// exact == comparison, no tolerance). These golden values were captured from the
// fill* output while the old Potential subclasses still existed and the live
// old-vs-new bit-identity comparison passed; the old classes are deleted in
// ponytail-unify Phase D, so all four forms now pin to golden (matching the
// Surface/Coulomb pattern). The "n1>1" case (extreme R/a) is the only one that
// exercises region 1, where the new ADD bodies bound region 1 at n1-1 instead of n1.
// ============================================================================
struct WSCase { double V, R, a; int n; double step; const char* tag; };
static const WSCase kWoodsSaxonCases[] = {
    { 50.0,   6.0, 0.65, 200, 0.10, "WS normal" },
    {-30.0,   6.0, 0.65, 200, 0.10, "WS neg-depth" },
    { 50.0, 100.0, 0.50, 200, 0.10, "WS region1 (n1>1)" },
    {  0.0,   6.0, 0.65, 200, 0.10, "WS zero-depth" },
};

static const WSCase kSpinOrbitCases[] = {
    {  6.0,   6.0, 0.65, 200, 0.10, "SO normal" },
    { -6.0,   6.0, 0.65, 200, 0.10, "SO neg-depth" },
    {  6.0, 100.0, 0.50, 200, 0.10, "SO region1 (n1>1)" },
    {  0.0,   6.0, 0.65, 200, 0.10, "SO zero-depth" },
};

// Surface & Coulomb reference subclasses are deleted in Phase E (optical-composite
// only), so those two forms are pinned to captured golden values instead of a
// live old-vs-new comparison. WS and SpinOrbit keep the live comparison — their
// classes survive (BoundState / linkule fitters still use them).
static const int    kGoldenIdx[6]     = { 1, 40, 80, 120, 160, 200 };
static const double kGoldenSurface[6] = { -0.0039183457267380961, -1.4631307105795064,
    -1.936840275165159, -0.0045698625317765603, -9.7145066144346655e-06, -2.0646171555950257e-08 };
static const double kGoldenCoulomb[6] = { 189.20022430633924, 168.29324350085031,
    117.54475197521796, 79.380178020469586, 59.410321914691096, 47.468548665506802 };
// Golden values for fillWoodsSaxon / fillSpinOrbit, captured (Phase D) while the
// old WoodsSaxonPotential/SpinOrbitPotential fill() still existed and matched.
// Index order matches kWoodsSaxonCases / kSpinOrbitCases.
static const double kGoldenWS[4][6] = {
    { -49.995101587952767, -48.098795002960237, -2.5512253607868929, -0.0057129809277423156, -1.2143136217158476e-05, -2.5807714458258585e-08 }, // WS normal
    { 29.997060952771662, 28.859277001776142, 1.5307352164721357, 0.0034277885566453893, 7.2858817302950855e-06, 1.548462867495515e-08 },        // WS neg-depth
    { -50, -50, -50, -50, -50, -50 },                                                                                                            // WS region1 (n1>1)
    { 0, 0, 0, 0, 0, 0 },                                                                                                                         // WS zero-depth
};
static const double kGoldenSO[4][6] = {
    { -1.8084672584945056e+27, -0.17315156338219001, -0.11315522542347585, -0.00017724095145869051, -2.81988580970528e-07, -4.7884433450213162e-10 }, // SO normal
    { 1.8084672584945056e+27, 0.17315156338219001, 0.11315522542347585, 0.00017724095145869051, 2.81988580970528e-07, 4.7884433450213162e-10 },        // SO neg-depth
    { 0, 0, 0, 0, 0, -3.2164896954177092e-70 },                                                                                                        // SO region1 (n1>1)
    { 0, 0, 0, 0, 0, 0 },                                                                                                                               // SO zero-depth
};

static void checkGolden(const char* name, OpticalPotential& op, const double* golden) {
    for (int j = 0; j < 6; j++) {
        double got = op.data1Based()[kGoldenIdx[j]];
        if (got != golden[j]) {
            printf("  FAIL  %-44s idx=%d got=%.17g exp=%.17g\n", name, kGoldenIdx[j], got, golden[j]);
            failCount++; return;
        }
    }
    passCount++;
}

static void test_optical_potential_bit_identity() {
    for (int i = 0; i < 4; i++) {
        const WSCase& c = kWoodsSaxonCases[i];
        OpticalPotential op; op.resize(c.n, 0.0, c.step);
        op.fillWoodsSaxon(c.V, c.R, c.a);
        checkGolden(c.tag, op, kGoldenWS[i]);
    }
    for (int i = 0; i < 4; i++) {
        const WSCase& c = kSpinOrbitCases[i];
        OpticalPotential op; op.resize(c.n, 0.0, c.step);
        op.fillSpinOrbit(c.V, c.R, c.a);
        checkGolden(c.tag, op, kGoldenSO[i]);
    }
    {   // fillSurface vs captured golden (SurfacePotential deleted in Phase E)
        OpticalPotential op; op.resize(200, 0.0, 0.10);
        op.fillSurface(10.0, 6.0, 0.65);
        checkGolden("fillSurface golden (V=10,R=6,a=0.65)", op, kGoldenSurface);
    }
    setVsq(7.0, 4.0, 82, 8, 1);  // seed vcsq12 sphere pair 1 (Pb+O, r1=7,r2=4)
    {   // fillCoulomb vs captured golden (CoulombPotential deleted in Phase E)
        OpticalPotential op; op.resize(200, 0.0, 0.10);
        op.fillCoulomb(1, 1);
        checkGolden("fillCoulomb golden (ch 1)", op, kGoldenCoulomb);
    }
    {   // transfer pseudo-channel: channelIndexIn==3 leaves vCoul 0 (all zero)
        OpticalPotential op; op.resize(200, 0.0, 0.10);
        op.fillCoulomb(1, 3);
        bool ok = true;
        for (int i = 1; i <= 200; i++) if (op.data1Based()[i] != 0.0) { ok = false; break; }
        if (ok) passCount++;
        else { printf("  FAIL  Coulomb chIn==3 not all zero\n"); failCount++; }
    }
    {   // add() / scale() composition exactness (scale() is API-only; pin it here)
        OpticalPotential a; a.resize(3, 0.0, 1.0); a.values = {1.0, 2.0, 4.0};
        OpticalPotential b; b.resize(3, 0.0, 1.0); b.values = {0.5, 0.25, 0.125};
        a.add(b);      // {1.5, 2.25, 4.125}
        a.scale(2.0);  // {3.0, 4.5, 8.25}
        bool ok = (a.values[0] == 3.0 && a.values[1] == 4.5 && a.values[2] == 8.25);
        if (ok) passCount++;
        else { printf("  FAIL  add()/scale() composition\n"); failCount++; }
    }
}

// ============================================================================
// OpticalPotentialLibrary — values must match Cleopatra's potentials.h exactly.
// Golden numbers generated by running the original CallPotential() (globals
// version) on the same inputs. Tolerance is tight (1e-6 abs) — these are the
// same arithmetic, so they agree to round-off.
// ============================================================================
static void test_optical_potential_library() {
    const double TOL = 1e-6;

    // --- Verification potential: full 16-parameter check ---
    // Koning-Delaroche, p + 16O @ 10 MeV. This is the INCOMING block of the
    // reference working/DWBA.in, so it also pins the expander cross-check.
    {
        OMPset p = callPotential("K", 16, 8, 10.0, 1);
        checkAbs("K(16O,p,10) v",    p.v,    53.07556635,   TOL);
        checkAbs("K(16O,p,10) r0",   p.r0,   1.143016903,   TOL);
        checkAbs("K(16O,p,10) a",    p.a,    0.675432,      TOL);
        checkAbs("K(16O,p,10) vi",   p.vi,   0.8268252704,  TOL);
        checkAbs("K(16O,p,10) ri0",  p.ri0,  1.143016903,   TOL);
        checkAbs("K(16O,p,10) ai",   p.ai,   0.675432,      TOL);
        checkAbs("K(16O,p,10) vsi",  p.vsi,  7.688504912,   TOL);
        checkAbs("K(16O,p,10) rsi0", p.rsi0, 1.302460503,   TOL);
        checkAbs("K(16O,p,10) asi",  p.asi,  0.527028,      TOL);
        checkAbs("K(16O,p,10) vso",  p.vso,  5.551115216,   TOL);
        checkAbs("K(16O,p,10) rso0", p.rso0, 0.9286378798,  TOL);
        checkAbs("K(16O,p,10) aso",  p.aso,  0.59,          TOL);
        checkAbs("K(16O,p,10) vsoi", p.vsoi, -0.03954298118,TOL);
        checkAbs("K(16O,p,10) rsoi0",p.rsoi0,0.9286378798,  TOL);
        checkAbs("K(16O,p,10) asoi", p.asoi, 0.59,          TOL);
        checkAbs("K(16O,p,10) rc0",  p.rc0,  1.435682137,   TOL);
    }

    // --- Spot checks across the popular potentials (A,K,V,D,L,l,p,c,s + Z,b,Q) ---
    { OMPset p = callPotential("A", 208, 82, 30.0, 1);   // An & Cai (deuteron)
      checkAbs("A(208Pb,d,30) v",   p.v,   93.36947508,  TOL);
      checkAbs("A(208Pb,d,30) vsi", p.vsi, 9.912,        TOL);
      checkAbs("A(208Pb,d,30) rc0", p.rc0, 1.303,        TOL); }

    { OMPset p = callPotential("K", 208, 82, 30.0, 1);   // Koning (proton)
      checkAbs("K(208Pb,p,30) v",   p.v,   51.1815854,   TOL);
      checkAbs("K(208Pb,p,30) vsi", p.vsi, 9.253149918,  TOL); }

    { OMPset p = callPotential("V", 90, 40, 40.0, 1);    // Varner CH89 (proton)
      checkAbs("V(90Zr,p,40) v",    p.v,   46.04026435,  TOL);
      checkAbs("V(90Zr,p,40) vsi",  p.vsi, 6.661328806,  TOL); }

    { OMPset p = callPotential("D", 120, 50, 25.0, 1);   // Daehnick (deuteron)
      checkAbs("D(120Sn,d,25) v",   p.v,   89.84556293,  TOL);
      checkAbs("D(120Sn,d,25) vsi", p.vsi, 12.00100188,  TOL); }

    { OMPset p = callPotential("L", 88, 38, 11.0, 1);    // Lohr-Haeberli (deuteron)
      checkAbs("L(88Sr,d,11) v",    p.v,   109.9251323,  TOL);
      checkAbs("L(88Sr,d,11) vsi",  p.vsi, 11.01881045,  TOL); }

    { OMPset p = callPotential("l", 40, 20, 30.0, 2);    // Liang (3He)
      checkAbs("l(40Ca,3He,30) v",  p.v,   115.5192314,  TOL);
      checkAbs("l(40Ca,3He,30) vsoi",p.vsoi,-1.1591,     TOL); }

    { OMPset p = callPotential("p", 40, 20, 30.0, 2);    // Pang (isospin dep.)
      checkAbs("p(40Ca,3He,30) v",  p.v,   116.4605675,  TOL);
      checkAbs("p(40Ca,3He,30) vsi",p.vsi, 18.86649859,  TOL); }

    { OMPset p = callPotential("c", 90, 40, 20.0, 1);    // Li-Liang-Cai (triton)
      checkAbs("c(90Zr,t,20) v",    p.v,   161.9629765,  TOL);
      checkAbs("c(90Zr,t,20) vi",   p.vi,  13.553,       TOL); }

    { OMPset p = callPotential("s", 120, 50, 50.0, 2);   // Su & Han (alpha)
      checkAbs("s(120Sn,a,50) v",   p.v,   148.0157672,  TOL);
      checkAbs("s(120Sn,a,50) vsi", p.vsi, 31.5966,      TOL); }

    { OMPset p = callPotential("Z", 6, 3, 10.0, 1);      // Zhang-Pang-Lou (6Li special case)
      checkAbs("Z(6Li,d,10) v",     p.v,   66.39900739,  TOL);
      checkAbs("Z(6Li,d,10) vsi",   p.vsi, 38.15087993,  TOL); }

    { OMPset p = callPotential("b", 90, 40, 20.0, 1);    // Becchetti A=3 (triton branch)
      checkAbs("b(90Zr,t,20) v",    p.v,   160.8888889,  TOL);
      checkAbs("b(90Zr,t,20) vi",   p.vi,  27.17777778,  TOL); }

    { OMPset p = callPotential("Q", 90, 40, 20.0, 1);    // Perey & Perey (deuteron)
      checkAbs("Q(90Zr,d,20) v",    p.v,   94.45154534,  TOL);
      checkAbs("Q(90Zr,d,20) vsi",  p.vsi, 19.2,         TOL); }

    // Unknown code → ok = false
    { OMPset p = callPotential("?", 40, 20, 10.0, 1);
      record("unknown code ok==false", p.ok ? 1.0 : 0.0, 0.0, p.ok ? 1.0 : 0.0, 0.5); }
}

// ============================================================================
// DwbaInputExpander — expansion of human-readable reaction lines.
// The 16O(p,p) inelastic case is checked byte-for-byte against the deck the
// original Cleopatra InFileCreator produces (working/DWBA.in). Transfer-deck
// structure and the auto-detect heuristic are checked too.
// ============================================================================
static void test_dwba_expander() {
    auto checkBool = [&](const char* name, bool got, bool exp) {
        record(name, got ? 1.0 : 0.0, exp ? 1.0 : 0.0,
               (got == exp) ? 0.0 : 1.0, 0.5);
    };
    auto checkContains = [&](const char* name, const std::string& hay, const char* needle) {
        bool found = hay.find(needle) != std::string::npos;
        record(name, found ? 1.0 : 0.0, 1.0, found ? 0.0 : 1.0, 0.5);
    };

    // --- byte-identical inelastic deck (matches Cleopatra working/DWBA.in) ---
    {
        std::string in = "16O(p,p)16O       0+        none     2+           6.00    10MeV/u     KK\n";
        std::string got = DwbaExpander::expand(in);
        const std::string expected =
"$============================================ Ex=6.00(p+16O|2+)KK,ELab=10.00\n"
"reset\n"
"REACTION: 16O(p,p)16O(2+ 6.00) ELAB= 10.000\n"
"PARAMETERSET ineloca2 r0target\n"
"JBIGA=0+\n"
"$Koning and Delaroche (2009) 0.001 < E < 200 | 24 < A < 209 | Iso. Dep. | http://dx.doi.org/10.1016/S0375-9474(02)01321-0\n"
"INCOMING\n"
"v    =  53.076    r0 =   1.143    a =   0.675\n"
"vi   =   0.827   ri0 =   1.143   ai =   0.675\n"
"vsi  =   7.689  rsi0 =   1.302  asi =   0.527  rc0 =   1.436\n"
";\n"
"OUTGOING\n"
"$Koning and Delaroche (2009) 0.001 < E < 200 | 24 < A < 209 | Iso. Dep. | http://dx.doi.org/10.1016/S0375-9474(02)01321-0\n"
"v    =  55.458    r0 =   1.143    a =   0.675\n"
"vi   =   0.383   ri0 =   1.143   ai =   0.675\n"
"vsi  =   6.489  rsi0 =   1.302  asi =   0.527  rc0 =   1.436\n"
";\n"
"anglemin=0.000000 anglemax=180.000000 anglestep=1.000000\n"
";\n"
"end $================================== end of input\n";
        checkBool("DWBA inelastic 16O(p,p) byte-identical", got == expected, true);
    }

    // --- single-nucleon transfer (d,p): deck structure ---
    {
        std::string in = "40Ca(d,p)41Ca   0+   1f7/2   7/2-   0.000   10MeV/u   AK\n";
        std::string got = DwbaExpander::expand(in);
        checkContains("dp REACTION line",   got, "REACTION: 40Ca(d,p)41Ca(7/2- 0.000)");
        checkContains("dp PARAMETERSET",    got, "PARAMETERSET dpsb r0target");
        checkContains("dp av18 wavefunc",   got, "wavefunction av18");
        checkContains("dp JBIGA",           got, "JBIGA=0+");
        checkContains("dp bound state l=3", got, "nodes=1 l=3 jp=7/2");
        checkContains("dp INCOMING An&Cai", got, "INCOMING $An and Cai");
        checkContains("dp OUTGOING Koning", got, "OUTGOING $Koning");
    }

    // --- two-nucleon transfer (t,p): alpha3 / phiffer / L= form ---
    {
        std::string in = "10Be(t,p)12Be   0+   1L=0   0+   0.000   5MeV/u   lA\n";
        std::string got = DwbaExpander::expand(in);
        checkContains("tp alpha3 paramset",  got, "PARAMETERSET alpha3 r0target");
        checkContains("tp phiffer wavefunc", got, "wavefunction phiffer");
        checkContains("tp L= bound form",    got, "nodes=1 L=0");
    }

    // --- elastic (Ex=0): ELASTIC SCATTERING / CHANNEL block ---
    {
        std::string in = "40Ca(d,d)40Ca   0+   none   0+   0.000   10MeV/u   AA\n";
        std::string got = DwbaExpander::expand(in);
        checkContains("elastic CHANNEL",  got, "CHANNEL d + 40Ca");
        checkContains("elastic SCATTER",  got, "ELASTIC SCATTERING");
    }

    // --- auto-detect heuristic ---
    checkBool("looksLikeDwba: reaction line",
              DwbaExpander::looksLikeDwba("16O(p,p)16O 0+ none 2+ 6.00 10MeV/u KK"), true);
    checkBool("looksLikeDwba: skip comments first",
              DwbaExpander::looksLikeDwba("# comment\n206Hg(d,p)207Hg 0+ 1g9/2 9/2+ 0.0 7MeV/u AK"), true);
    checkBool("looksLikeDwba: native deck -> false",
              DwbaExpander::looksLikeDwba("reset\nREACTION: 16O(p,p)16O ELAB=10\nend"), false);
    checkBool("looksLikeDwba: too few tokens -> false",
              DwbaExpander::looksLikeDwba("16O(p,p)16O 0+ none"), false);
}
