#include "CeeV0LiteTransport.h"
#include "CeeV0LiteVertexFitter.h"

#include <cassert>
#include <cmath>
#include <iostream>

namespace
{
bool Close(double a, double b, double tolerance)
{
  return std::fabs(a - b) < tolerance;
}

CeeV0Lite::TrackState MakeState(const CeeV0Lite::Vec3& position,
                                const CeeV0Lite::Vec3& momentum,
                                int charge)
{
  CeeV0Lite::TrackState state;
  state.positionCm = position;
  state.momentumGeV = momentum;
  state.charge = charge;
  for (int i = 0; i < 3; ++i) {
    state.covariance.m[i][i] = 0.01 * 0.01;
    state.covariance.m[i + 3][i + 3] = 0.001 * 0.001;
  }
  return state;
}

class FixedRandom
{
public:
  explicit FixedRandom(unsigned int seed) : fState(seed) {}

  double Uniform(double minimum, double maximum)
  {
    fState = 1664525u * fState + 1013904223u;
    const double unit = double(fState) / 4294967295.;
    return minimum + (maximum - minimum) * unit;
  }

private:
  unsigned int fState;
};
}

int main()
{
  CeeV0Lite::TransportConfig transportConfig;
  const CeeV0Lite::TrackState start = MakeState(
    CeeV0Lite::Vec3(0., 0., 0.), CeeV0Lite::Vec3(1., 0., 0.), +1);

  CeeV0Lite::TrackState straight;
  assert(CeeV0Lite::PropagateState(start, 10., 0., transportConfig, straight));
  assert(Close(straight.positionCm[0], 10., 1.e-10));
  assert(Close(straight.positionCm[1], 0., 1.e-10));
  assert(Close(CeeV0Lite::Norm(straight.momentumGeV), 1., 1.e-12));

  const double path = 100.;
  const double angle = 0.00299792458 * path;
  CeeV0Lite::TrackState helix;
  assert(CeeV0Lite::PropagateState(start, path, 1., transportConfig, helix));
  assert(Close(helix.momentumGeV[0], std::cos(angle), 1.e-10));
  assert(Close(helix.momentumGeV[1], 0., 1.e-10));
  assert(Close(helix.momentumGeV[2], std::sin(angle), 1.e-10));
  assert(Close(CeeV0Lite::Norm(helix.momentumGeV), 1., 1.e-12));

  CeeV0Lite::TrackState returned;
  assert(CeeV0Lite::PropagateState(helix, -path, 1., transportConfig, returned));
  assert(CeeV0Lite::Norm(returned.positionCm - start.positionCm) < 1.e-8);
  assert(CeeV0Lite::Norm(returned.momentumGeV - start.momentumGeV) < 1.e-10);

  // Closest approach now reports convergence explicitly.
  CeeV0Lite::TrackState atPoint;
  double closestPath = 0.;
  CeeV0Lite::ClosestApproachResult closest;
  assert(CeeV0Lite::TransportToClosestApproach(
    start, CeeV0Lite::Vec3(3., 4., 0.), 0., transportConfig,
    atPoint, closestPath, &closest));
  assert(closest.valid && closest.converged && !closest.hitPathLimit);
  assert(Close(closestPath, 3., 1.e-10));
  assert(Close(closest.stationarityResidualCm, 0., 1.e-10));

  // A root outside the allowed path is rejected by default but remains
  // inspectable for diagnostics.
  CeeV0Lite::TransportConfig bounded = transportConfig;
  bounded.maximumAbsPathCm = 1.;
  assert(!CeeV0Lite::TransportToClosestApproach(
    start, CeeV0Lite::Vec3(5., 1., 0.), 0., bounded,
    atPoint, closestPath, &closest));
  assert(closest.valid && !closest.converged && closest.hitPathLimit);
  assert(Close(closestPath, 1., 1.e-12));
  bounded.requireClosestApproachConvergence = false;
  assert(CeeV0Lite::TransportToClosestApproach(
    start, CeeV0Lite::Vec3(5., 1., 0.), 0., bounded,
    atPoint, closestPath, &closest));
  assert(closest.valid && !closest.converged && closest.hitPathLimit);

  // With no Newton iteration the tangent seed is insufficient for a curved
  // track, so non-convergence is reported and rejected.
  CeeV0Lite::TransportConfig noClosestIterations = transportConfig;
  noClosestIterations.maximumClosestApproachIterations = 0;
  assert(!CeeV0Lite::TransportToClosestApproach(
    start, CeeV0Lite::Vec3(80., 20., 0.), 1., noClosestIterations,
    atPoint, closestPath, &closest));
  assert(closest.valid && !closest.converged && closest.iterations == 0);
  noClosestIterations.requireClosestApproachConvergence = false;
  assert(CeeV0Lite::TransportToClosestApproach(
    start, CeeV0Lite::Vec3(80., 20., 0.), 1., noClosestIterations,
    atPoint, closestPath, &closest));
  assert(closest.valid && !closest.converged && closest.iterations == 0);

  // A one-seed Newton solve used to accept s=0 here, although it is the
  // farthest point on the circle.  The bounded global POCA must inspect the
  // other helix roots and choose either of the equivalent half-turn minima.
  const double lowMomentum = 0.2;
  const double testField = 0.5;
  const double helixRadius = lowMomentum /
    (0.00299792458 * testField);
  const CeeV0Lite::TrackState lowMomentumTrack = MakeState(
    CeeV0Lite::Vec3(0., 0., 0.),
    CeeV0Lite::Vec3(lowMomentum, 0., 0.), +1);
  CeeV0Lite::TransportConfig multiTurn = transportConfig;
  multiTurn.maximumAbsPathCm = 1000.;
  const CeeV0Lite::Vec3 oppositeCirclePoint(0., 0., 2. * helixRadius);
  assert(CeeV0Lite::TransportToClosestApproach(
    lowMomentumTrack, oppositeCirclePoint, testField, multiTurn,
    atPoint, closestPath, &closest));
  assert(closest.valid && closest.converged && !closest.hitPathLimit);
  assert(CeeV0Lite::Norm(atPoint.positionCm - oppositeCirclePoint) < 1.e-7);
  assert(Close(std::fabs(closestPath), std::acos(-1.) * helixRadius, 1.e-6));
  assert(std::fabs(closestPath) > 100.);

  // The same global search must work for negative curvature.  Reversing the
  // charge reflects the circle and its opposite point across the x axis.
  const CeeV0Lite::TrackState negativeTrack = MakeState(
    CeeV0Lite::Vec3(0., 0., 0.),
    CeeV0Lite::Vec3(lowMomentum, 0., 0.), -1);
  const CeeV0Lite::Vec3 oppositeNegativeCirclePoint(
    0., 0., -2. * helixRadius);
  assert(CeeV0Lite::TransportToClosestApproach(
    negativeTrack, oppositeNegativeCirclePoint, testField, multiTurn,
    atPoint, closestPath, &closest));
  assert(closest.valid && closest.converged && !closest.hitPathLimit);
  assert(CeeV0Lite::Norm(atPoint.positionCm -
                         oppositeNegativeCirclePoint) < 1.e-7);
  assert(Close(std::fabs(closestPath), std::acos(-1.) * helixRadius, 1.e-6));

  // A helix with longitudinal motion has geometrically different turns.
  // A target placed after more than one full turn must select that unique
  // turn, not the nearest transverse copy.
  CeeV0Lite::TrackState pitchedTrack = MakeState(
    CeeV0Lite::Vec3(0., 0., 0.),
    CeeV0Lite::Vec3(0.16, 0.12, 0.), +1);
  CeeV0Lite::TransportConfig longPath = transportConfig;
  longPath.maximumAbsPathCm = 2500.;
  const double pitchedPeriod = 2. * std::acos(-1.) * helixRadius;
  const double targetPath = 1.25 * pitchedPeriod;
  CeeV0Lite::TrackState pitchedTargetState;
  assert(CeeV0Lite::PropagateState(pitchedTrack, targetPath, testField,
                                   longPath, pitchedTargetState, false));
  assert(CeeV0Lite::TransportToClosestApproach(
    pitchedTrack, pitchedTargetState.positionCm, testField, longPath,
    atPoint, closestPath, &closest));
  assert(closest.valid && closest.converged && !closest.hitPathLimit);
  assert(CeeV0Lite::Norm(atPoint.positionCm -
                         pitchedTargetState.positionCm) < 1.e-7);
  assert(Close(closestPath, targetPath, 1.e-6));

  // If the true half-turn minimum is outside the allowed interval, s=0 is
  // still a maximum and must not be returned.  The global constrained result
  // is the path boundary and remains a diagnostic-only state by default.
  CeeV0Lite::TransportConfig shortArc = multiTurn;
  shortArc.maximumAbsPathCm = 0.25 * std::acos(-1.) * helixRadius;
  assert(!CeeV0Lite::TransportToClosestApproach(
    lowMomentumTrack, oppositeCirclePoint, testField, shortArc,
    atPoint, closestPath, &closest));
  assert(closest.valid && !closest.converged && closest.hitPathLimit);
  assert(Close(std::fabs(closestPath), shortArc.maximumAbsPathCm, 1.e-8));
  shortArc.requireClosestApproachConvergence = false;
  assert(CeeV0Lite::TransportToClosestApproach(
    lowMomentumTrack, oppositeCirclePoint, testField, shortArc,
    atPoint, closestPath, &closest));
  assert(closest.valid && !closest.converged && closest.hitPathLimit);

  // Fixed-seed global cross-check: positive/negative curvature, transverse
  // (uy=0) and pitched helices are compared with every point of a fine grid
  // over the complete allowed path.  The analytic/root-refined result must
  // never be farther away than any sampled path or either boundary.
  FixedRandom random(0x5eed1234u);
  for (int trial = 0; trial < 48; ++trial) {
    const double azimuth = random.Uniform(-std::acos(-1.), std::acos(-1.));
    const double longitudinal = (trial % 3 == 0)
      ? 0. : random.Uniform(-0.8, 0.8);
    const double transverse = std::sqrt(1. - longitudinal * longitudinal);
    const double magnitude = random.Uniform(0.02, 0.20);
    const CeeV0Lite::Vec3 momentum(
      magnitude * transverse * std::cos(azimuth),
      magnitude * longitudinal,
      magnitude * transverse * std::sin(azimuth));
    const CeeV0Lite::TrackState randomTrack = MakeState(
      CeeV0Lite::Vec3(random.Uniform(-20., 20.),
                      random.Uniform(-20., 20.),
                      random.Uniform(-20., 20.)),
      momentum, (trial % 2 == 0) ? +1 : -1);
    const CeeV0Lite::Vec3 randomPoint(random.Uniform(-200., 200.),
                                      random.Uniform(-200., 200.),
                                      random.Uniform(-200., 200.));
    CeeV0Lite::TransportConfig globalCheck = transportConfig;
    globalCheck.maximumAbsPathCm = 600.;
    globalCheck.requireClosestApproachConvergence = false;
    const double randomField = (trial % 4 < 2) ? 0.5 : -0.5;
    assert(CeeV0Lite::TransportToClosestApproach(
      randomTrack, randomPoint, randomField, globalCheck,
      atPoint, closestPath, &closest));
    assert(closest.valid);
    const double solvedDistanceSquared =
      CeeV0Lite::Norm2(atPoint.positionCm - randomPoint);
    double sampledDistanceSquared = 1.e300;
    const int samples = 20000;
    for (int sample = 0; sample <= samples; ++sample) {
      const double sampledPath = -globalCheck.maximumAbsPathCm +
        2. * globalCheck.maximumAbsPathCm * sample / samples;
      CeeV0Lite::TrackState sampledState;
      assert(CeeV0Lite::PropagateState(
        randomTrack, sampledPath, randomField, globalCheck,
        sampledState, false));
      sampledDistanceSquared = std::min(sampledDistanceSquared,
        CeeV0Lite::Norm2(sampledState.positionCm - randomPoint));
    }
    const double gridTolerance = 1.e-7 *
      std::max(1., sampledDistanceSquared);
    assert(solvedDistanceSquared <= sampledDistanceSquared + gridTolerance);
  }

  // TOF beta boundaries must fail cleanly rather than produce infinities.
  assert(CeeV0Lite::EstimateTofMomentum(0.938, 100., 10., 0.1, 0.02).valid);
  assert(!CeeV0Lite::EstimateTofMomentum(0.938, 100., 0., 0.1, 0.02).valid);
  assert(!CeeV0Lite::EstimateTofMomentum(0.938, 100., 3., 0.1, 0.02).valid);

  const CeeV0Lite::Vec3 trueVertex(1., 2., 3.);
  CeeV0Lite::Vec3 firstDirection(1., 0.2, 0.5);
  CeeV0Lite::Vec3 secondDirection(-0.3, 1., 0.4);
  CeeV0Lite::Normalize(firstDirection, firstDirection);
  CeeV0Lite::Normalize(secondDirection, secondDirection);
  const CeeV0Lite::TrackState first = MakeState(
    trueVertex + 20. * firstDirection, firstDirection, +1);
  const CeeV0Lite::TrackState second = MakeState(
    trueVertex + 15. * secondDirection, secondDirection, -1);

  CeeV0Lite::VertexFitter fitter;
  CeeV0Lite::VertexFitResult fit;
  assert(fitter.Fit(first, second, 0., fit));
  assert(fit.valid);
  assert(fit.converged);
  assert(fit.status == CeeV0Lite::kVertexFitConverged);
  assert(CeeV0Lite::Norm(fit.vertexCm - trueVertex) < 1.e-5);
  assert(fit.daughterDcaCm < 1.e-5);

  // Parallel and almost-parallel tracks have an ill-conditioned vertex and
  // are rejected before the 1/sin(opening angle)^2 seed can explode.
  const CeeV0Lite::TrackState almostParallel = MakeState(
    CeeV0Lite::Vec3(0., 1., 0.), CeeV0Lite::Vec3(1., 1.e-5, 0.), -1);
  assert(!fitter.Fit(start, almostParallel, 0., fit));
  assert(!fit.valid && !fit.converged);
  assert(fit.status == CeeV0Lite::kVertexFitDegenerateGeometry);

  // Maximum-iteration exit is a hard failure by default and is available as
  // a diagnostic-only fit when explicitly requested.
  CeeV0Lite::VertexFitterConfig oneIteration;
  oneIteration.maximumIterations = 1;
  oneIteration.convergenceCm = 0.;
  CeeV0Lite::VertexFitter strictOneIteration(oneIteration, transportConfig);
  assert(!strictOneIteration.Fit(first, second, 0., fit));
  assert(!fit.converged);
  assert(fit.status == CeeV0Lite::kVertexFitMaximumIterations);
  oneIteration.requireConvergence = false;
  CeeV0Lite::VertexFitter diagnosticOneIteration(oneIteration,
                                                  transportConfig);
  assert(diagnosticOneIteration.Fit(first, second, 0., fit));
  assert(fit.valid && !fit.converged);
  assert(fit.status == CeeV0Lite::kVertexFitMaximumIterations);

  std::cout << "math, transport boundaries, convergence and vertex tests passed"
            << std::endl;
  return 0;
}
