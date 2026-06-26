// OpticalPotential.cpp — the flat computed-grid optical potential.
// Phase B: grid management + composition. The fill* builders (one per V(r) form)
// are added in Phase C, each copied verbatim from the old Potential subclass.

#include "OpticalPotential.h"
#include "Constants.h"

#include <algorithm>
#include <cmath>

// vcsq12 — Coulomb sphere/point split kernel (source_misc.cpp). Per-point:
// writes x = Coulomb potential at rValue for sphere pair channelIndex.
void vcsq12(double rValue, double& x, int channelIndex);

void OpticalPotential::resize(int n, double r0, double step)
{
    nPts     = n;
    rStart   = r0;
    stepSize = step;
    values.assign(n, 0.0);
}

void OpticalPotential::clear()
{
    std::fill(values.begin(), values.end(), 0.0);
}

void OpticalPotential::add(const OpticalPotential& other)
{
    for (std::size_t i = 0; i < values.size(); ++i)
        values[i] += other.values[i];
}

void OpticalPotential::scale(double k)
{
    for (double& v : values)
        v *= k;
}


// fillWoodsSaxon — straight Woods-Saxon (potForm 1), copied from
// WoodsSaxonPotential::fill (source_potentials.cpp): region-bracket via n1/n2,
// then the multiplicative exp recurrence -V/(1+y), y*=yStep. ADDS to values
// (data1Based, 1..nPts); region-1 bounded at n1-1 (region-2 owns n1) so the
// composed result is byte-identical to the former overwrite fill. The grid spec
// (rStart, stepSize) is read from the members set by resize().
void OpticalPotential::fillWoodsSaxon(double V, double R, double a)
{
    const double A = a;
    double* vRay = data1Based();

    int n1 = 0;
    int n2 = 0;
    if (V != 0.0)
    {
        double xStep = stepSize / A;
        double yStep = std::exp(xStep);
        double x = (rStart - stepSize - R) / A;

        // Break region into three regions
        n1 = (int)(-(100.0 + x) / xStep);
        n1 = std::max(n1, 1);
        n1 = std::min(n1, nPts);
        double yVal = 46.0 + std::log(std::fabs(V));
        yVal = std::min(yVal, Constants::BIGLOG);
        n2 = (int)((yVal - x) / xStep);
        n2 = std::min(n2, nPts);
        n2 = std::max(n2, n1);

        // Region 1 — straight W.S. saturates at -V (n1-1: region 2 owns n1)
        if (n1 != 1) {
            for (int i = 1; i <= n1 - 1; i++) vRay[i] += -V;
        }

        // Region 2 — straight W.S.
        x = x + n1 * xStep;
        double y = std::exp(x);
        for (int i = n1; i <= n2; i++) {
            vRay[i] += -V / (1.0 + y);
            y = y * yStep;
        }
    }

    // Fill out rest with 0 (ADD: leaves any pre-existing contribution untouched)
    if (n2 == nPts) return;
    for (int i = n2 + 1; i <= nPts; i++) vRay[i] += 0.0;
}


// fillCoulomb — built-in Coulomb via the vcsq12 sphere/point split, copied from
// CoulombPotential::fill. Per-point: rValue runs 0, step, 2*step, ... (step is
// the member stepSize, i.e. the channel rStart). channelIndexIn == 3 (the transfer
// pseudo-channel) leaves vCoul at 0 → adds 0. ADD semantics; this is the first
// fill of the real-central buffer (onto zeroed values, so += vCoul == vCoul).
void OpticalPotential::fillCoulomb(int channelIndex, int channelIndexIn)
{
    double* vRay = data1Based();
    double rValue = -stepSize;
    double vCoul = 0.0;
    for (int i = 1; i <= nPts; i++) {
        rValue = rValue + stepSize;
        if (channelIndexIn != 3) vcsq12(rValue, vCoul, channelIndex);
        vRay[i] += vCoul;
    }
}


// fillSurface — derivative-Woods-Saxon surface-imaginary (potForm 3), copied
// from SurfacePotential::fill: same deriv-WS recurrence as fillSpinOrbit's
// region 2 but WITHOUT the (0.5/A)/r spin-orbit conversion. ADD semantics,
// region-1 bounded at n1-1; built standalone (imv surface term).
void OpticalPotential::fillSurface(double V, double R, double a)
{
    const double A = a;
    double* vRay = data1Based();

    int n1 = 0;
    int n2 = 0;
    if (V != 0.0)
    {
        double xStep = stepSize / A;
        double yStep = std::exp(xStep);
        double x = (rStart - stepSize - R) / A;

        // Break region into three regions
        n1 = (int)(-(100.0 + x) / xStep);
        n1 = std::max(n1, 1);
        n1 = std::min(n1, nPts);
        double yVal = 46.0 + std::log(std::fabs(V));
        yVal = std::min(yVal, Constants::BIGLOG);
        n2 = (int)((yVal - x) / xStep);
        n2 = std::min(n2, nPts);
        n2 = std::max(n2, n1);

        // Region 1 — derivative shape is zero in the saturated interior
        if (n1 != 1) {
            for (int i = 1; i <= n1 - 1; i++) vRay[i] += 0.0;
        }

        // Region 2 — derivative of W.S.
        x = x + n1 * xStep;
        double y = std::exp(x);
        for (int i = n1; i <= n2; i++) {
            vRay[i] += (-4.0 * V) * (y / (1.0 + y)) / (1.0 + y);
            y = y * yStep;
        }
    }

    // Fill out rest with 0 (ADD: leaves any pre-existing contribution untouched)
    if (n2 == nPts) return;
    for (int i = n2 + 1; i <= nPts; i++) vRay[i] += 0.0;
}


// fillSpinOrbit — derivative-Woods-Saxon converted to spin-orbit (potForm 2),
// copied from SpinOrbitPotential::fill. Region-1 zero (n1-1 bound, as fillWoodsSaxon).
// The old fill ran two region-2 loops — first storing the deriv-WS recurrence
// -4V*(y/(1+y))/(1+y), then a separate sweep applying the (0.5/A)/rT spin-orbit
// conversion in place. Here the two are FUSED into one ADD: storing then
// reloading a double is exact, so deriv reused directly gives the identical
// final value while keeping the contribution additive (no in-place overwrite).
// Spin-orbit is always built standalone (soR/soI), so the zeroed-buffer no-op
// region-1/tail adds match the old overwrite.
void OpticalPotential::fillSpinOrbit(double V, double R, double a)
{
    const double A = a;
    double* vRay = data1Based();

    int n1 = 0;
    int n2 = 0;
    if (V != 0.0)
    {
        double xStep = stepSize / A;
        double yStep = std::exp(xStep);
        double x = (rStart - stepSize - R) / A;

        // Break region into three regions
        n1 = (int)(-(100.0 + x) / xStep);
        n1 = std::max(n1, 1);
        n1 = std::min(n1, nPts);
        double yVal = 46.0 + std::log(std::fabs(V));
        yVal = std::min(yVal, Constants::BIGLOG);
        n2 = (int)((yVal - x) / xStep);
        n2 = std::min(n2, nPts);
        n2 = std::max(n2, n1);

        // Region 1 — derivative shape is zero in the saturated interior
        if (n1 != 1) {
            for (int i = 1; i <= n1 - 1; i++) vRay[i] += 0.0;
        }

        // Region 2 — derivative of W.S., fused with the spin-orbit (0.5/A)/r
        // conversion.
        x = x + n1 * xStep;
        double y = std::exp(x);
        double rT = rStart + (n1 - 1) * stepSize + 1.0e-30;
        for (int i = n1; i <= n2; i++) {
            double deriv = (-4.0 * V) * (y / (1.0 + y)) / (1.0 + y);
            vRay[i] += (0.5 / A) * deriv / rT;
            y = y * yStep;
            rT = rT + stepSize;
        }
    }

    // Fill out rest with 0 (ADD: leaves any pre-existing contribution untouched)
    if (n2 == nPts) return;
    for (int i = n2 + 1; i <= nPts; i++) vRay[i] += 0.0;
}
