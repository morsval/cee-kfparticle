#include "CeeV0LiteTransport.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace
{
// For cm, GeV/c and Tesla with dp/dt = q v x B:
// CEE beam axis is z; the uniform field is (0, By, 0).
// d(px)/ds = -q*k*By*pz/p, d(py)/ds = 0, d(pz)/ds = q*k*By*px/p.
const double kMagneticRotationPerCm = 0.00299792458;

double Clamp(double value, double minimum, double maximum)
{
  return std::max(minimum, std::min(maximum, value));
}

void StateToArray(const CeeV0Lite::TrackState& state, double values[6])
{
  for (int i = 0; i < 3; ++i) {
    values[i] = state.positionCm[i];
    values[i + 3] = state.momentumGeV[i];
  }
}

void ArrayToState(const double values[6], CeeV0Lite::TrackState& state)
{
  for (int i = 0; i < 3; ++i) {
    state.positionCm[i] = values[i];
    state.momentumGeV[i] = values[i + 3];
  }
}

bool PropagateParameters(const CeeV0Lite::TrackState& input,
                         double signedPathCm,
                         double magneticFieldByTesla,
                         CeeV0Lite::TrackState& output)
{
  output = input;
  const double momentum = CeeV0Lite::Norm(input.momentumGeV);
  if (!(momentum > 1.e-12) || !std::isfinite(momentum)) return false;

  const double omega = kMagneticRotationPerCm * input.charge *
                       magneticFieldByTesla / momentum;
  if (std::fabs(omega) < 1.e-12) {
    const CeeV0Lite::Vec3 direction = input.momentumGeV / momentum;
    output.positionCm = input.positionCm + signedPathCm * direction;
    return CeeV0Lite::IsFinite(output.positionCm);
  }

  const double angle = omega * signedPathCm;
  const double sine = std::sin(angle);
  const double cosine = std::cos(angle);
  const double px = input.momentumGeV[0];
  const double py = input.momentumGeV[1];
  const double pz = input.momentumGeV[2];
  const double denominator = momentum * omega;

  output.positionCm[0] += (px * sine - pz * (1. - cosine)) / denominator;
  output.positionCm[1] += signedPathCm * py / momentum;
  output.positionCm[2] += (pz * sine + px * (1. - cosine)) / denominator;

  output.momentumGeV[0] = px * cosine - pz * sine;
  output.momentumGeV[1] = py;
  output.momentumGeV[2] = pz * cosine + px * sine;
  return CeeV0Lite::IsFinite(output.positionCm) &&
         CeeV0Lite::IsFinite(output.momentumGeV);
}

struct PathEvaluation
{
  PathEvaluation()
    : valid(false), distanceSquaredCm2(std::numeric_limits<double>::infinity()),
      stationarityResidualCm(0.), stationarityDerivative(0.)
  {
  }

  bool valid;
  double distanceSquaredCm2;
  double stationarityResidualCm;
  double stationarityDerivative;
  CeeV0Lite::TrackState state;
};

bool EvaluatePath(const CeeV0Lite::TrackState& input,
                  const CeeV0Lite::Vec3& pointCm,
                  double magneticFieldByTesla,
                  double signedPathCm,
                  PathEvaluation& evaluation)
{
  evaluation = PathEvaluation();
  if (!PropagateParameters(input, signedPathCm, magneticFieldByTesla,
                           evaluation.state)) return false;
  const double momentum = CeeV0Lite::Norm(evaluation.state.momentumGeV);
  if (!(momentum > 1.e-12) || !std::isfinite(momentum)) return false;

  const CeeV0Lite::Vec3 direction =
    evaluation.state.momentumGeV / momentum;
  const CeeV0Lite::Vec3 residual =
    evaluation.state.positionCm - pointCm;
  const double omega = kMagneticRotationPerCm * evaluation.state.charge *
                       magneticFieldByTesla / momentum;
  const CeeV0Lite::Vec3 directionDerivative(-omega * direction[2],
                                            0., omega * direction[0]);
  evaluation.distanceSquaredCm2 = CeeV0Lite::Norm2(residual);
  evaluation.stationarityResidualCm = CeeV0Lite::Dot(residual, direction);
  evaluation.stationarityDerivative = 1. +
    CeeV0Lite::Dot(residual, directionDerivative);
  evaluation.valid = std::isfinite(evaluation.distanceSquaredCm2) &&
                     std::isfinite(evaluation.stationarityResidualCm) &&
                     std::isfinite(evaluation.stationarityDerivative);
  return evaluation.valid;
}

struct ClosestCandidate
{
  ClosestCandidate()
    : pathCm(0.), converged(false), atPathBoundary(false),
      singularDerivative(false), iterations(0)
  {
  }

  double pathCm;
  bool converged;
  bool atPathBoundary;
  bool singularDerivative;
  int iterations;
  PathEvaluation evaluation;
};

// Add all copies of a periodic path which lie in [minimumPath, maximumPath].
// The caller proves that this interval is no wider than two periods, including
// the exactly transverse case.
void AppendPeriodicPaths(double basePath,
                         double period,
                         double minimumPath,
                         double maximumPath,
                         std::vector<double>& paths)
{
  if (!(period > 0.) || !std::isfinite(period) ||
      !(minimumPath <= maximumPath)) return;
  double offset = std::fmod(basePath - minimumPath, period);
  if (offset < 0.) offset += period;
  double path = minimumPath + offset;
  const double pathScale = std::max(period,
    std::max(std::fabs(minimumPath), std::fabs(maximumPath)));
  const double pathTolerance = 64. *
    std::numeric_limits<double>::epsilon() * pathScale;
  if (path > maximumPath + pathTolerance) return;

  // The caller searches an interval no wider than two periods, so this loop
  // adds at most three representable copies and needs no truncation cap.
  for (;;) {
    if (path >= minimumPath - pathTolerance &&
        path <= maximumPath + pathTolerance)
      paths.push_back(Clamp(path, minimumPath, maximumPath));
    const double next = path + period;
    if (!(next > path) || next > maximumPath + pathTolerance) break;
    path = next;
  }
}

double HelixStationarity(double pathCm,
                         double omega,
                         double coefficientCos,
                         double coefficientSin,
                         double coefficientConstant,
                         double longitudinalDirectionSquared)
{
  const double angle = omega * pathCm;
  return coefficientCos * std::cos(angle) +
         coefficientSin * std::sin(angle) + coefficientConstant +
         longitudinalDirectionSquared * pathCm;
}

double HelixStationarityDerivative(double pathCm,
                                   double omega,
                                   double coefficientCos,
                                   double coefficientSin,
                                   double longitudinalDirectionSquared)
{
  const double angle = omega * pathCm;
  return omega * (-coefficientCos * std::sin(angle) +
                   coefficientSin * std::cos(angle)) +
         longitudinalDirectionSquared;
}

// Refine one sign-changing root of d(distance^2)/ds / 2.  Newton is used
// only inside its bracket; otherwise the step falls back to bisection.
ClosestCandidate RefineStationaryPath(
    const CeeV0Lite::TrackState& input,
    const CeeV0Lite::Vec3& pointCm,
    double magneticFieldByTesla,
    double lowerPath,
    double upperPath,
    double lowerValue,
    double upperValue,
    double omega,
    double coefficientCos,
    double coefficientSin,
    double coefficientConstant,
    double longitudinalDirectionSquared,
    const CeeV0Lite::TransportConfig& config)
{
  ClosestCandidate candidate;
  double path = (std::fabs(lowerValue) < std::fabs(upperValue))
    ? lowerPath : upperPath;
  double value = (path == lowerPath) ? lowerValue : upperValue;

  if (std::fabs(value) <= config.closestApproachToleranceCm) {
    candidate.converged = true;
  } else {
    for (int iteration = 0;
         iteration < config.maximumClosestApproachIterations;
         ++iteration) {
      candidate.iterations = iteration + 1;
      const double derivative = HelixStationarityDerivative(
        path, omega, coefficientCos, coefficientSin,
        longitudinalDirectionSquared);
      double proposed = std::numeric_limits<double>::quiet_NaN();
      if (std::isfinite(derivative) && std::fabs(derivative) > 1.e-12)
        proposed = path - value / derivative;
      else
        candidate.singularDerivative = true;
      if (!std::isfinite(proposed) || !(proposed > lowerPath) ||
          !(proposed < upperPath))
        proposed = 0.5 * (lowerPath + upperPath);

      const double proposedValue = HelixStationarity(
        proposed, omega, coefficientCos, coefficientSin,
        coefficientConstant, longitudinalDirectionSquared);
      if (!std::isfinite(proposedValue)) break;
      path = proposed;
      value = proposedValue;
      if (std::fabs(value) <= config.closestApproachToleranceCm) {
        candidate.converged = true;
        break;
      }
      if ((lowerValue <= 0. && value >= 0.) ||
          (lowerValue >= 0. && value <= 0.)) {
        upperPath = path;
        upperValue = value;
      } else {
        lowerPath = path;
        lowerValue = value;
      }
    }
  }

  candidate.pathCm = path;
  EvaluatePath(input, pointCm, magneticFieldByTesla, path,
               candidate.evaluation);
  if (candidate.evaluation.valid &&
      std::fabs(candidate.evaluation.stationarityResidualCm) <=
        config.closestApproachToleranceCm)
    candidate.converged = true;
  return candidate;
}

bool IsLocalMinimum(const ClosestCandidate& candidate,
                    const CeeV0Lite::TrackState& input,
                    const CeeV0Lite::Vec3& pointCm,
                    double magneticFieldByTesla,
                    double minimumPath,
                    double maximumPath,
                    double period,
                    double stationarityTolerance)
{
  if (!candidate.evaluation.valid) return false;
  if (candidate.evaluation.stationarityDerivative > 1.e-10) return true;
  if (candidate.evaluation.stationarityDerivative < -1.e-10) return false;

  // A degenerate stationary point is accepted only if the distance does not
  // decrease to either side.  This retains flat minima and rejects maxima or
  // tangent roots of the stationarity equation.
  const double probe = std::max(1.e-5,
    std::min(0.01 * period, 10. * stationarityTolerance));
  const double leftPath = std::max(minimumPath, candidate.pathCm - probe);
  const double rightPath = std::min(maximumPath, candidate.pathCm + probe);
  PathEvaluation left;
  PathEvaluation right;
  if (!EvaluatePath(input, pointCm, magneticFieldByTesla, leftPath, left) ||
      !EvaluatePath(input, pointCm, magneticFieldByTesla, rightPath, right))
    return false;
  const double tolerance = 1.e-12 * std::max(1.,
    candidate.evaluation.distanceSquaredCm2);
  return candidate.evaluation.distanceSquaredCm2 <=
           left.distanceSquaredCm2 + tolerance &&
         candidate.evaluation.distanceSquaredCm2 <=
           right.distanceSquaredCm2 + tolerance;
}

bool IsBetterCandidate(const ClosestCandidate& candidate,
                       const ClosestCandidate& best,
                       double preferredPath)
{
  if (!candidate.evaluation.valid) return false;
  if (!best.evaluation.valid) return true;
  const double tolerance = 1.e-12 * std::max(1., std::max(
    candidate.evaluation.distanceSquaredCm2,
    best.evaluation.distanceSquaredCm2));
  if (candidate.evaluation.distanceSquaredCm2 <
      best.evaluation.distanceSquaredCm2 - tolerance) return true;
  if (candidate.evaluation.distanceSquaredCm2 >
      best.evaluation.distanceSquaredCm2 + tolerance) return false;
  if (candidate.converged != best.converged) return candidate.converged;
  if (candidate.atPathBoundary != best.atPathBoundary)
    return !candidate.atPathBoundary;
  const double candidateOffset = std::fabs(candidate.pathCm - preferredPath);
  const double bestOffset = std::fabs(best.pathCm - preferredPath);
  if (candidateOffset < bestOffset - 1.e-12) return true;
  if (candidateOffset > bestOffset + 1.e-12) return false;
  return candidate.pathCm > best.pathCm;
}

bool FindClosestApproachParameters(
    const CeeV0Lite::TrackState& input,
    const CeeV0Lite::Vec3& pointCm,
    double magneticFieldByTesla,
    const CeeV0Lite::TransportConfig& config,
    CeeV0Lite::TrackState& output,
    CeeV0Lite::ClosestApproachResult& result)
{
  result = CeeV0Lite::ClosestApproachResult();
  CeeV0Lite::Vec3 initialDirection;
  if (!CeeV0Lite::Normalize(input.momentumGeV, initialDirection) ||
      !CeeV0Lite::IsFinite(pointCm) ||
      !(config.maximumAbsPathCm > 0.) ||
      !std::isfinite(config.maximumAbsPathCm) ||
      !(config.closestApproachToleranceCm >= 0.) ||
      !std::isfinite(config.closestApproachToleranceCm) ||
      config.maximumClosestApproachIterations < 0) return false;

  const double unconstrainedPath = CeeV0Lite::Dot(
    pointCm - input.positionCm, initialDirection);
  if (!std::isfinite(unconstrainedPath)) return false;

  const double pathLimit = config.maximumAbsPathCm;
  const double preferredPath = Clamp(unconstrainedPath, -pathLimit, pathLimit);
  const double inputMomentum = CeeV0Lite::Norm(input.momentumGeV);
  const double omega = kMagneticRotationPerCm * input.charge *
                       magneticFieldByTesla / inputMomentum;

  // Preserve the exact straight-line projection, including boundary and
  // diagnostic behavior.  This is the path used for the current zero-field
  // CEE data and for neutral tracks.
  if (std::fabs(omega) < 1.e-12) {
    PathEvaluation straight;
    if (!EvaluatePath(input, pointCm, magneticFieldByTesla, preferredPath,
                      straight)) {
      result.signedPathCm = preferredPath;
      return false;
    }
    output = straight.state;
    result.valid = true;
    result.signedPathCm = preferredPath;
    result.stationarityResidualCm = straight.stationarityResidualCm;
    result.converged = std::fabs(straight.stationarityResidualCm) <=
      config.closestApproachToleranceCm;
    result.hitPathLimit = !result.converged &&
      std::fabs(preferredPath) >= pathLimit * (1. - 1.e-12);
    return result.converged || !config.requireClosestApproachConvergence;
  }

  const double ux = initialDirection[0];
  const double uy = initialDirection[1];
  const double uz = initialDirection[2];
  const double longitudinalDirectionSquared = uy * uy;
  // "Longitudinal" here means parallel to the field (y), not the beam (z).
  // The x-z helix projection is a circle with this centre. Expressing the
  // stationarity equation through the centre removes a cancelling radius
  // term and gives
  //   g(s) = A cos(omega s) + B sin(omega s) + C + uy^2 s.
  const double centreX = input.positionCm[0] - uz / omega;
  const double centreZ = input.positionCm[2] + ux / omega;
  const double centreDx = centreX - pointCm[0];
  const double centreDz = centreZ - pointCm[2];
  const double coefficientCos = centreDx * ux + centreDz * uz;
  const double coefficientSin = -centreDx * uz + centreDz * ux;
  const double coefficientConstant =
    uy * (input.positionCm[1] - pointCm[1]);
  const double transverseAmplitude = std::hypot(coefficientCos,
                                                 coefficientSin);
  const double period = 2. * std::acos(-1.) / std::fabs(omega);

  std::vector<ClosestCandidate> candidates;
  candidates.reserve(10);
  const double boundaryPaths[2] = {-pathLimit, pathLimit};
  for (int i = 0; i < 2; ++i) {
    ClosestCandidate boundary;
    boundary.pathCm = boundaryPaths[i];
    boundary.atPathBoundary = true;
    if (EvaluatePath(input, pointCm, magneticFieldByTesla, boundary.pathCm,
                     boundary.evaluation)) {
      boundary.converged = std::fabs(
        boundary.evaluation.stationarityResidualCm) <=
          config.closestApproachToleranceCm;
      candidates.push_back(boundary);
    }
  }

  // The transverse distance is exactly periodic.  For any transverse phase,
  // adding or subtracting one period leaves that part unchanged and moves y
  // closer to its unconstrained optimum until the nearest representative (or
  // a path boundary) is reached.  Therefore the bounded global minimum must
  // lie within one period of the longitudinal optimum; remote turns cannot
  // win and never need to be enumerated.
  double representativePath = preferredPath;
  if (uy != 0.)
    representativePath = Clamp(
      (pointCm[1] - input.positionCm[1]) / uy,
      -pathLimit, pathLimit);
  const double searchMinimum = std::max(-pathLimit,
                                        representativePath - period);
  const double searchMaximum = std::min(pathLimit,
                                        representativePath + period);

  std::vector<double> breakpoints;
  breakpoints.reserve(8);
  if (searchMinimum <= searchMaximum) {
    breakpoints.push_back(searchMinimum);
    breakpoints.push_back(searchMaximum);
    if (transverseAmplitude > 0.) {
      const double ratio = -longitudinalDirectionSquared /
        (omega * transverseAmplitude);
      if (std::fabs(ratio) <= 1.) {
        const double phase = std::atan2(-coefficientCos, coefficientSin);
        const double angle = std::acos(Clamp(ratio, -1., 1.));
        AppendPeriodicPaths((phase + angle) / omega, period,
                            searchMinimum, searchMaximum, breakpoints);
        AppendPeriodicPaths((phase - angle) / omega, period,
                            searchMinimum, searchMaximum, breakpoints);
      }
    }
  }
  std::sort(breakpoints.begin(), breakpoints.end());
  // Only bit-identical values are merged.  In particular, a very small helix
  // period or a shallow minimum/maximum pair must not disappear behind an
  // absolute path tolerance.
  breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end()),
                    breakpoints.end());

  std::vector<ClosestCandidate> stationaryCandidates;
  stationaryCandidates.reserve(8);
  for (std::size_t i = 0; i < breakpoints.size(); ++i) {
    const double value = HelixStationarity(
      breakpoints[i], omega, coefficientCos, coefficientSin,
      coefficientConstant, longitudinalDirectionSquared);
    if (std::isfinite(value) &&
        std::fabs(value) <= config.closestApproachToleranceCm) {
      ClosestCandidate root;
      root.pathCm = breakpoints[i];
      root.converged = true;
      EvaluatePath(input, pointCm, magneticFieldByTesla, root.pathCm,
                   root.evaluation);
      stationaryCandidates.push_back(root);
    }
  }
  for (std::size_t i = 0; i + 1 < breakpoints.size(); ++i) {
    const double lowerPath = breakpoints[i];
    const double upperPath = breakpoints[i + 1];
    const double lowerValue = HelixStationarity(
      lowerPath, omega, coefficientCos, coefficientSin,
      coefficientConstant, longitudinalDirectionSquared);
    const double upperValue = HelixStationarity(
      upperPath, omega, coefficientCos, coefficientSin,
      coefficientConstant, longitudinalDirectionSquared);
    if (!std::isfinite(lowerValue) || !std::isfinite(upperValue) ||
        lowerValue == 0. || upperValue == 0. ||
        ((lowerValue < 0.) == (upperValue < 0.))) continue;
    stationaryCandidates.push_back(RefineStationaryPath(
      input, pointCm, magneticFieldByTesla, lowerPath, upperPath,
      lowerValue, upperValue, omega, coefficientCos, coefficientSin,
      coefficientConstant, longitudinalDirectionSquared, config));
  }

  for (std::size_t i = 0; i < stationaryCandidates.size(); ++i) {
    const ClosestCandidate& candidate = stationaryCandidates[i];
    if (!IsLocalMinimum(candidate, input, pointCm, magneticFieldByTesla,
                        -pathLimit, pathLimit, period,
                        config.closestApproachToleranceCm)) continue;
    candidates.push_back(candidate);
  }

  ClosestCandidate best;
  for (std::size_t i = 0; i < candidates.size(); ++i)
    if (IsBetterCandidate(candidates[i], best, preferredPath))
      best = candidates[i];
  if (!best.evaluation.valid) return false;

  output = best.evaluation.state;
  result.valid = true;
  result.converged = best.converged;
  result.hitPathLimit = best.atPathBoundary && !best.converged;
  result.singularDerivative = best.singularDerivative;
  result.iterations = best.iterations;
  result.signedPathCm = best.pathCm;
  result.stationarityResidualCm =
    best.evaluation.stationarityResidualCm;
  return result.converged || !config.requireClosestApproachConvergence;
}
// Derivative of exactly the same fixed-path By map used for central values.
// The path is a signed 3D length in cm (not STAR's s/p parameter).
bool FixedPathJacobian(const CeeV0Lite::TrackState& input,
                       double path, double field,
                       const CeeV0Lite::TrackState& output,
                       double jacobian[6][6])
{
  const double p = CeeV0Lite::Norm(input.momentumGeV);
  if (!(p > 1.e-8) || !std::isfinite(p)) return false;
  const double k = kMagneticRotationPerCm * input.charge * field;
  const double omega = k / p;
  if (!std::isfinite(omega)) return false;
  // Avoid differentiating across the legacy straight/helix branch threshold.
  if (std::fabs(omega) > 0.5e-12 && std::fabs(omega) < 2.e-12)
    return false;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) jacobian[i][j] = 0.;
  for (int i = 0; i < 3; ++i) jacobian[i][i] = 1.;
  if (std::fabs(omega) < 1.e-12) {
    const CeeV0Lite::Vec3 u = input.momentumGeV / p;
    for (int i = 0; i < 3; ++i) {
      jacobian[i+3][i+3] = 1.;
      for (int j = 0; j < 3; ++j)
        jacobian[i][j+3] = path / p *
          (((i == j) ? 1. : 0.) - u[i]*u[j]);
    }
  } else {
    const double angle = omega * path;
    const double sn = std::sin(angle), cs = std::cos(angle);
    const double rotation[3][3] = {{cs,0.,-sn},{0.,1.,0.},{sn,0.,cs}};
    const double direct[3][3] = {{sn/k,0.,-(1.-cs)/k},
                                  {0.,path/p,0.},
                                  {(1.-cs)/k,0.,sn/k}};
    for (int j = 0; j < 3; ++j) {
      const double dangle = -angle * (input.momentumGeV[j]/p) / p;
      jacobian[0][j+3] = direct[0][j] + output.momentumGeV[0]/k*dangle;
      jacobian[1][j+3] = direct[1][j] - path/p *
        (input.momentumGeV[1]/p)*(input.momentumGeV[j]/p);
      jacobian[2][j+3] = direct[2][j] + output.momentumGeV[2]/k*dangle;
      jacobian[3][j+3] = rotation[0][j] - output.momentumGeV[2]*dangle;
      jacobian[4][j+3] = rotation[1][j];
      jacobian[5][j+3] = rotation[2][j] + output.momentumGeV[0]*dangle;
    }
  }
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      if (!std::isfinite(jacobian[i][j])) return false;
  return true;
}

// Two small matrix products instead of a four-index sum. Same full C,
// including position-momentum correlations; noise is added by the caller.
CeeV0Lite::Mat6 PropagateCovariance(const double jacobian[6][6],
                                    const CeeV0Lite::Mat6& covariance)
{
  double jc[6][6] = {{0.}};
  CeeV0Lite::Mat6 result;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      for (int k = 0; k < 6; ++k)
        jc[i][j] += jacobian[i][k]*covariance.m[k][j];
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      for (int k = 0; k < 6; ++k)
        result.m[i][j] += jc[i][k]*jacobian[j][k];
  return result;
}
}

namespace CeeV0Lite
{

TransportConfig::TransportConfig()
  : maximumClosestApproachIterations(12),
    closestApproachToleranceCm(1.e-5),
    maximumAbsPathCm(1000.),
    requireClosestApproachConvergence(true),
    useAnalyticFixedPathJacobian(true),
    positionProcessNoiseCm2PerCm(0.),
    angularProcessNoiseRad2PerCm(0.)
{
}

ClosestApproachResult::ClosestApproachResult()
  : valid(false),
    converged(false),
    hitPathLimit(false),
    singularDerivative(false),
    iterations(0),
    signedPathCm(0.),
    stationarityResidualCm(0.)
{
}

bool PropagateState(const TrackState& input,
                    double signedPathCm,
                    double magneticFieldByTesla,
                    const TransportConfig& config,
                    TrackState& output,
                    bool propagateCovariance)
{
  if (!std::isfinite(signedPathCm) ||
      std::fabs(signedPathCm) > config.maximumAbsPathCm) return false;
  if (!PropagateParameters(input, signedPathCm, magneticFieldByTesla, output))
    return false;
  if (!propagateCovariance) return true;

  double jacobian[6][6] = {{0.}};
  double central[6];
  StateToArray(input, central);

  const bool analytic = config.useAnalyticFixedPathJacobian &&
    FixedPathJacobian(input, signedPathCm, magneticFieldByTesla, output, jacobian);
  for (int column = 0; !analytic && column < 6; ++column) {
    double step = (column < 3)
      ? 1.e-4
      : std::max(1.e-7, std::fabs(central[column]) * 1.e-5);
    TrackState plus = input;
    TrackState minus = input;
    double plusValues[6];
    double minusValues[6];
    StateToArray(input, plusValues);
    StateToArray(input, minusValues);
    plusValues[column] += step;
    minusValues[column] -= step;
    ArrayToState(plusValues, plus);
    ArrayToState(minusValues, minus);

    TrackState propagatedPlus;
    TrackState propagatedMinus;
    if (!PropagateParameters(plus, signedPathCm, magneticFieldByTesla,
                             propagatedPlus) ||
        !PropagateParameters(minus, signedPathCm, magneticFieldByTesla,
                             propagatedMinus)) return false;
    StateToArray(propagatedPlus, plusValues);
    StateToArray(propagatedMinus, minusValues);
    for (int row = 0; row < 6; ++row)
      jacobian[row][column] = (plusValues[row] - minusValues[row]) / (2. * step);
  }

  Mat6 propagatedCovariance = PropagateCovariance(jacobian, input.covariance);

  const double absPath = std::fabs(signedPathCm);
  const double positionNoise = config.positionProcessNoiseCm2PerCm * absPath;
  for (int i = 0; i < 3; ++i)
    propagatedCovariance.m[i][i] += positionNoise;

  const double momentum = Norm(output.momentumGeV);
  Vec3 direction;
  if (Normalize(output.momentumGeV, direction)) {
    const double angularNoise = momentum * momentum *
      config.angularProcessNoiseRad2PerCm * absPath;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        propagatedCovariance.m[i + 3][j + 3] += angularNoise *
          (((i == j) ? 1. : 0.) - direction[i] * direction[j]);
      }
    }
  }

  Symmetrize(propagatedCovariance);
  output.covariance = propagatedCovariance;
  return true;
}

bool ReanchorAtClosestApproach(const TrackState& input,
                              const Vec3& pointCm,
                              double magneticFieldByTesla,
                              const TransportConfig& config,
                              TrackState& output,
                              ClosestApproachResult* diagnostics)
{
  ClosestApproachResult result;
  TrackState nominal, propagated;
  TransportConfig reference = config;
  reference.requireClosestApproachConvergence = true;
  reference.positionProcessNoiseCm2PerCm = 0.;
  reference.angularProcessNoiseRad2PerCm = 0.;
  bool ok = IsFinite(input.positionCm) && IsFinite(input.momentumGeV) &&
    std::isfinite(magneticFieldByTesla) &&
    IsPositiveSemidefinite(input.covariance) &&
    FindClosestApproachParameters(input, pointCm, magneticFieldByTesla,
                                   reference, nominal, result) &&
    result.converged && !result.hitPathLimit;
  if (ok) {
    ok = PropagateState(input, result.signedPathCm, magneticFieldByTesla,
                         reference, propagated) &&
         IsPositiveSemidefinite(propagated.covariance);
  }
  if (diagnostics) *diagnostics = result;
  if (!ok) return false;
  output = propagated;
  return true;
}

bool TransportToClosestApproach(const TrackState& input, const Vec3& pointCm,
                                double magneticFieldByTesla, const TransportConfig& config,
                                TrackState& output, double& signedPathCm,
                                ClosestApproachResult* diagnostics)
{
  return TransportToClosestApproach(input, pointCm, magneticFieldByTesla, config,
                                    output, signedPathCm, diagnostics, 0);
}

bool TransportToClosestApproach(const TrackState& input,
                                const Vec3& pointCm,
                                double magneticFieldByTesla,
                                const TransportConfig& config,
                                TrackState& output,
                                double& signedPathCm,
                                ClosestApproachResult* diagnostics,
                                double (*stateJacobian)[6])
{
  ClosestApproachResult result;
  if (!FindClosestApproachParameters(input, pointCm, magneticFieldByTesla,
                                     config, output, result)) {
    signedPathCm = result.signedPathCm;
    if (diagnostics) *diagnostics = result;
    return false;
  }
  signedPathCm = result.signedPathCm;

  // Recompute the closest path for every input perturbation.  This includes
  // ds_closest/d(state) and makes a reference-point slide along the same
  // trajectory invisible to the transported covariance.
  double jacobian[6][6] = {{0.}};
  double central[6];
  StateToArray(input, central);
  for (int column = 0; column < 6; ++column) {
    const double step = column < 3
      ? 1.e-4 : std::max(1.e-7, std::fabs(central[column]) * 1.e-5);
    double plusValues[6];
    double minusValues[6];
    StateToArray(input, plusValues);
    StateToArray(input, minusValues);
    plusValues[column] += step;
    minusValues[column] -= step;
    TrackState plus = input;
    TrackState minus = input;
    ArrayToState(plusValues, plus);
    ArrayToState(minusValues, minus);
    TrackState propagatedPlus;
    TrackState propagatedMinus;
    ClosestApproachResult plusResult;
    ClosestApproachResult minusResult;
    if (!FindClosestApproachParameters(plus, pointCm, magneticFieldByTesla,
                                       config, propagatedPlus, plusResult) ||
        !FindClosestApproachParameters(minus, pointCm, magneticFieldByTesla,
                                       config, propagatedMinus, minusResult)) {
      if (diagnostics) *diagnostics = result;
      return false;
    }
    StateToArray(propagatedPlus, plusValues);
    StateToArray(propagatedMinus, minusValues);
    for (int row = 0; row < 6; ++row)
      jacobian[row][column] =
        (plusValues[row] - minusValues[row]) / (2. * step);
  }

  Mat6 propagatedCovariance = PropagateCovariance(jacobian, input.covariance);

  const double absPath = std::fabs(signedPathCm);
  const double positionNoise = config.positionProcessNoiseCm2PerCm * absPath;
  for (int i = 0; i < 3; ++i)
    propagatedCovariance.m[i][i] += positionNoise;
  const double momentum = Norm(output.momentumGeV);
  Vec3 direction;
  if (Normalize(output.momentumGeV, direction)) {
    const double angularNoise = momentum * momentum *
      config.angularProcessNoiseRad2PerCm * absPath;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        propagatedCovariance.m[i + 3][j + 3] += angularNoise *
          (((i == j) ? 1. : 0.) - direction[i] * direction[j]);
  }
  Symmetrize(propagatedCovariance);
  output.covariance = propagatedCovariance;
  if (stateJacobian)
    for (int i=0;i<6;++i) for (int j=0;j<6;++j) stateJacobian[i][j]=jacobian[i][j];
  if (diagnostics) *diagnostics = result;
  return true;
}

} // namespace CeeV0Lite
