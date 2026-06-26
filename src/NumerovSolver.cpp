// NumerovSolver.cpp — Numerov method implementation
//
// The Numerov recursion for u''(r) = f(r) * u(r):
//   (1 + h²f_{n+1}/12) * u_{n+1} = 2*(1 - 5h²f_n/12) * u_n
//                                    - (1 + h²f_{n-1}/12) * u_{n-1}

#include "NumerovSolver.h"

// ============================================================================
// Outward Numerov integration
// ============================================================================
void NumerovSolver::integrateOutward(const std::vector<double>& f,
                                      double u0, double u1)
{
    int pointCount = (int)f.size();
    solution.resize(pointCount);
    solution[0] = u0;
    if (pointCount < 2) return;
    solution[1] = u1;

    double h2 = stepSize * stepSize;

    for (int n = 1; n < pointCount - 1; n++) {
        double fnM1  = f[n - 1];
        double fn    = f[n];
        double fnP1  = f[n + 1];

        double numer = 2.0 * (1.0 - 5.0 * h2 * fn / 12.0) * solution[n]
                     - (1.0 + h2 * fnM1 / 12.0) * solution[n - 1];
        double denom = 1.0 + h2 * fnP1 / 12.0;

        solution[n + 1] = numer / denom;
    }

    nodesFound = countNodes();
}

// ============================================================================
// Count nodes (zero crossings)
// ============================================================================
int NumerovSolver::countNodes() const
{
    int nodes = 0;
    for (size_t i = 1; i < solution.size(); i++) {
        if (solution[i] * solution[i - 1] < 0.0) {
            nodes++;
        }
    }
    return nodes;
}
