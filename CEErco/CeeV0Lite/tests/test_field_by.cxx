#include "CeeV0LiteTransport.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

namespace {
using namespace CeeV0Lite;

// Independent reference: integrate dr/ds=p/|p| and dp/ds=k*q*(p/|p| x B).
// No helix formula or production transport function is used by this oracle.
void Derivative(const double state[6], int charge, double by, double d[6])
{
  const double p = std::sqrt(state[3]*state[3] + state[4]*state[4] +
                             state[5]*state[5]);
  const double b[3] = {0., by, 0.};
  for (int i = 0; i < 3; ++i) d[i] = state[3+i] / p;
  const double kq = 0.00299792458 * charge;
  d[3] = kq * (d[1]*b[2] - d[2]*b[1]);
  d[4] = kq * (d[2]*b[0] - d[0]*b[2]);
  d[5] = kq * (d[0]*b[1] - d[1]*b[0]);
}

void Integrate(const TrackState& input, double path, double by, double out[6])
{
  for (int i = 0; i < 3; ++i) {
    out[i] = input.positionCm[i];
    out[3+i] = input.momentumGeV[i];
  }
  const int steps = std::max(1, int(std::ceil(std::fabs(path) / 0.2)));
  const double h = path / steps;
  for (int step = 0; step < steps; ++step) {
    double k1[6], k2[6], k3[6], k4[6], tmp[6];
    Derivative(out, input.charge, by, k1);
    for (int i = 0; i < 6; ++i) tmp[i] = out[i] + 0.5*h*k1[i];
    Derivative(tmp, input.charge, by, k2);
    for (int i = 0; i < 6; ++i) tmp[i] = out[i] + 0.5*h*k2[i];
    Derivative(tmp, input.charge, by, k3);
    for (int i = 0; i < 6; ++i) tmp[i] = out[i] + h*k3[i];
    Derivative(tmp, input.charge, by, k4);
    for (int i = 0; i < 6; ++i)
      out[i] += h*(k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i])/6.;
  }
}

TrackState MakeState(const Vec3& momentum, int charge)
{
  TrackState state;
  state.positionCm = Vec3(0.3, -0.7, 1.2);
  state.momentumGeV = momentum;
  state.charge = charge;
  // Non-diagonal positive-definite covariance, including position-momentum
  // correlations: this detects accidental covariance-axis permutations.
  double lower[6][6] = {{0.}};
  for (int i = 0; i < 6; ++i) {
    lower[i][i] = i < 3 ? 0.02 : 0.003;
    for (int j = 0; j < i; ++j)
      lower[i][j] = lower[i][i] * 0.04 * (i-j);
  }
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      for (int k = 0; k < 6; ++k)
        state.covariance.m[i][j] += lower[i][k] * lower[j][k];
  return state;
}
}

int main()
{
  using namespace CeeV0Lite;
  TransportConfig config;
  config.maximumAbsPathCm = 2500.;

  // A positive particle moving along the CEE beam (+z) bends toward -x
  // in +By. Reversing either the charge or field reverses the deflection.
  for (int charge = -1; charge <= 1; charge += 2) {
    for (int polarity = -1; polarity <= 1; polarity += 2) {
      TrackState start = MakeState(Vec3(0., 0., 1.), charge), out;
      assert(PropagateState(start, 100., 0.5*polarity, config, out));
      assert(out.momentumGeV[0] * charge * polarity < 0.);
      assert((out.positionCm[0]-start.positionCm[0]) * charge * polarity < 0.);
      assert(std::fabs(out.momentumGeV[1]) < 1.e-14);
      assert(std::fabs(out.positionCm[1]-start.positionCm[1]) < 1.e-14);
    }
  }

  // A particle parallel to the field does not curve, even for nonzero By.
  TrackState parallel = MakeState(Vec3(0., -0.6, 0.), +1), out;
  assert(PropagateState(parallel, 40., 0.8, config, out));
  assert(Norm(out.positionCm - parallel.positionCm - Vec3(0., -40., 0.)) < 1.e-12);
  assert(Norm(out.momentumGeV-parallel.momentumGeV) < 1.e-12);
  double nearestPath = 0.;
  ClosestApproachResult closest;
  assert(TransportToClosestApproach(parallel, parallel.positionCm + Vec3(3.,-17.,4.),
    0.8, config, out, nearestPath, &closest));
  assert(std::fabs(nearestPath-17.) < 1.e-6);

  // Check all components, path reversal, |p| and py conservation against
  // the Lorentz-force ODE for both charges, field polarities and path signs.
  for (int charge = -1; charge <= 1; ++charge) {
    for (int polarity = -1; polarity <= 1; polarity += 2) {
      for (int sign = -1; sign <= 1; sign += 2) {
        TrackState start = MakeState(Vec3(0.2, 0.3, 1.1), charge), back;
        const double by = 0.8*polarity, path = 270.*sign;
        double ref[6];
        Integrate(start, path, by, ref);
        assert(PropagateState(start, path, by, config, out));
        for (int i = 0; i < 3; ++i) {
          assert(std::fabs(out.positionCm[i]-ref[i]) < 2.e-8);
          assert(std::fabs(out.momentumGeV[i]-ref[3+i]) < 1.e-10);
        }
        assert(std::fabs(Norm(out.momentumGeV)-Norm(start.momentumGeV)) < 1.e-12);
        assert(out.momentumGeV[1] == start.momentumGeV[1]);
        assert(PropagateState(out, -path, by, config, back));
        assert(Norm(back.positionCm-start.positionCm) < 1.e-8);
        assert(Norm(back.momentumGeV-start.momentumGeV) < 1.e-12);
      }
    }
  }

  // The numerical covariance must follow derivatives of the independent ODE
  // in the original (x,y,z,px,py,pz) order, including off-diagonal entries.
  TrackState start = MakeState(Vec3(0.2, 0.3, 1.1), +1);
  const double path = 65., by = -0.8;
  double jacobian[6][6] = {{0.}};
  for (int j = 0; j < 6; ++j) {
    TrackState plus = start, minus = start;
    const double h = j < 3 ? 1.e-4 : 1.e-6;
    if (j < 3) { plus.positionCm[j] += h; minus.positionCm[j] -= h; }
    else { plus.momentumGeV[j-3] += h; minus.momentumGeV[j-3] -= h; }
    double yp[6], ym[6];
    Integrate(plus, path, by, yp); Integrate(minus, path, by, ym);
    for (int i = 0; i < 6; ++i) jacobian[i][j] = (yp[i]-ym[i])/(2.*h);
  }
  assert(PropagateState(start, path, by, config, out));
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double reference = 0.;
      for (int a = 0; a < 6; ++a)
        for (int b = 0; b < 6; ++b)
          reference += jacobian[i][a]*start.covariance.m[a][b]*jacobian[j][b];
      assert(std::fabs(out.covariance.m[i][j]-reference) <
        1.e-8*std::max(1., std::fabs(reference)));
    }
  }
  assert(IsPositiveSemidefinite(out.covariance));

  // Independent target beyond one turn with nonzero pitch along y.
  start = MakeState(Vec3(0.16, 0.12, 0.), +1);
  const double targetPath = 1.25 * 2.*std::acos(-1.)*0.2/(0.00299792458*0.5);
  double target[6];
  Integrate(start, targetPath, 0.5, target);
  assert(TransportToClosestApproach(start, Vec3(target[0],target[1],target[2]),
    0.5, config, out, nearestPath, &closest));
  assert(closest.valid && closest.converged);
  assert(std::fabs(nearestPath-targetPath) < 1.e-4);
  std::cout << "CEE By Lorentz-force, polarity, covariance and multi-turn POCA tests passed\n";
}
