#include "CeeV0LiteVertexFitter.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
void ResetMotherJacobianDiagnostics(CeeV0Lite::MotherState& mother)
{
  mother.jacobianFitsAttempted = 0;
  mother.jacobianFitsSucceeded = 0;
  mother.jacobianFailureCount = 0;
  mother.jacobianPlusFailureMask = 0u;
  mother.jacobianMinusFailureMask = 0u;
  mother.jacobianFailedColumn = -1;
  mother.jacobianFailedSide = 0;
  mother.jacobianFirstFailureFitStatus = 0;
}

void InvalidateMotherCovariance(
    CeeV0Lite::MotherState& mother,
    CeeV0Lite::MotherCovarianceStatus status)
{
  const double invalid = std::numeric_limits<double>::quiet_NaN();
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      mother.covariance.m[i][j] = invalid;
  mother.covarianceValid = false;
  mother.covarianceStatus = status;
}

bool ClosestLineSeed(const CeeV0Lite::TrackState& first,
                     const CeeV0Lite::TrackState& second,
                     double minimumOpeningAngleRad,
                     CeeV0Lite::Vec3& seed,
                     bool& degenerate)
{
  degenerate = false;
  CeeV0Lite::Vec3 u;
  CeeV0Lite::Vec3 v;
  if (!CeeV0Lite::Normalize(first.momentumGeV, u) ||
      !CeeV0Lite::Normalize(second.momentumGeV, v)) return false;

  const CeeV0Lite::Vec3 w = first.positionCm - second.positionCm;
  const double uv = CeeV0Lite::Dot(u, v);
  const double uw = CeeV0Lite::Dot(u, w);
  const double vw = CeeV0Lite::Dot(v, w);
  const double denominator = 1. - uv * uv;
  const double minimumSine = std::sin(minimumOpeningAngleRad);
  if (!std::isfinite(minimumSine) || minimumOpeningAngleRad < 0. ||
      minimumOpeningAngleRad > 1.5707963267948966 ||
      denominator < minimumSine * minimumSine) {
    degenerate = true;
    return false;
  }
  const double firstPath = (uv * vw - uw) / denominator;
  const double secondPath = (vw - uv * uw) / denominator;

  seed = 0.5 * ((first.positionCm + firstPath * u) +
                (second.positionCm + secondPath * v));
  return CeeV0Lite::IsFinite(seed);
}

void FillProjection(const CeeV0Lite::Vec3& firstBasis,
                    const CeeV0Lite::Vec3& secondBasis,
                    double projection[2][3])
{
  for (int j = 0; j < 3; ++j) {
    projection[0][j] = firstBasis[j];
    projection[1][j] = secondBasis[j];
  }
}

void ProjectCovariance(const CeeV0Lite::Mat3& covariance,
                       const double projection[2][3],
                       double output[2][2])
{
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      output[i][j] = 0.;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          output[i][j] += projection[i][a] * covariance.m[a][b] *
                          projection[j][b];
    }
  }
}

void ProjectPoint(const CeeV0Lite::Vec3& point,
                  const double projection[2][3],
                  double output[2])
{
  for (int i = 0; i < 2; ++i) {
    output[i] = 0.;
    for (int j = 0; j < 3; ++j) output[i] += projection[i][j] * point[j];
  }
}

bool InvertPositiveDefinite2(const double input[2][2], double inverse[2][2])
{
  const double determinant = input[0][0] * input[1][1] -
                             input[0][1] * input[1][0];
  const double scale = std::max(std::fabs(input[0][0]),
                                std::fabs(input[1][1]));
  if (!(input[0][0] > 0.) || !(input[1][1] > 0.) || !(scale > 0.) ||
      !(determinant > 1.e-14 * scale * scale) ||
      !std::isfinite(determinant)) return false;
  return CeeV0Lite::InvertSymmetric2(input, inverse);
}

bool BatchUpdateVertex(const CeeV0Lite::TrackState trackAtPoint[2],
                       double varianceFloor,
                       CeeV0Lite::Vec3& vertex,
                       CeeV0Lite::Mat3& vertexCovariance)
{
  CeeV0Lite::Mat3 information;
  CeeV0Lite::Vec3 rightHandSide;
  for (int track = 0; track < 2; ++track) {
    CeeV0Lite::Vec3 direction;
    if (!CeeV0Lite::Normalize(trackAtPoint[track].momentumGeV, direction))
      return false;
    CeeV0Lite::Vec3 e1;
    CeeV0Lite::Vec3 e2;
    CeeV0Lite::MakeTransverseBasis(direction, e1, e2);
    double h[2][3];
    FillProjection(e1, e2, h);

    double measurementCovariance[2][2];
    ProjectCovariance(CeeV0Lite::PositionBlock(
      trackAtPoint[track].covariance), h, measurementCovariance);
    measurementCovariance[0][0] += varianceFloor;
    measurementCovariance[1][1] += varianceFloor;
    double inverseMeasurement[2][2];
    if (!InvertPositiveDefinite2(measurementCovariance,
                                 inverseMeasurement)) return false;

    double measured[2];
    ProjectPoint(trackAtPoint[track].positionCm, h, measured);
    for (int i = 0; i < 3; ++i) {
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
          rightHandSide[i] += h[a][i] * inverseMeasurement[a][b] *
                              measured[b];
      for (int j = 0; j < 3; ++j)
        for (int a = 0; a < 2; ++a)
          for (int b = 0; b < 2; ++b)
            information.m[i][j] += h[a][i] * inverseMeasurement[a][b] *
                                    h[b][j];
    }
  }
  if (!CeeV0Lite::InvertPositiveDefinite3(information,
                                           vertexCovariance)) return false;
  vertex = CeeV0Lite::Multiply(vertexCovariance, rightHandSide);
  return CeeV0Lite::IsFinite(vertex);
}

bool ConstructMotherAtCommonVertex(
    const CeeV0Lite::TrackState propagated[2],
    double varianceFloor,
    CeeV0Lite::TrackState daughterAtVertex[2],
    CeeV0Lite::MotherState& mother,
    double& chi2)
{
  // Joint input x=(r1,p1,r2,p2).  At this stage the two transported tracks
  // contain independent detector errors; common event-time uncertainty is
  // added once at the mother level by the Finder.
  double x[12] = {0.};
  double covariance[12][12] = {{0.}};
  for (int daughter = 0; daughter < 2; ++daughter) {
    const int offset = 6 * daughter;
    for (int i = 0; i < 3; ++i) {
      x[offset + i] = propagated[daughter].positionCm[i];
      x[offset + 3 + i] = propagated[daughter].momentumGeV[i];
    }
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j)
        covariance[offset + i][offset + j] =
          propagated[daughter].covariance.m[i][j];
  }

  // Each track has one free path parameter.  After eliminating those two
  // nuisances, a two-prong vertex has one effective residual: the separation
  // along the normal to both local directions.  Constraining all three raw
  // position components would spuriously constrain motion along the tracks.
  CeeV0Lite::Vec3 firstDirection;
  CeeV0Lite::Vec3 secondDirection;
  CeeV0Lite::Vec3 normal;
  if (!CeeV0Lite::Normalize(propagated[0].momentumGeV, firstDirection) ||
      !CeeV0Lite::Normalize(propagated[1].momentumGeV, secondDirection) ||
      !CeeV0Lite::Normalize(CeeV0Lite::Cross(firstDirection,
                                             secondDirection), normal))
    return false;
  double h[12] = {0.};
  for (int i = 0; i < 3; ++i) {
    h[i] = -normal[i];
    h[6 + i] = normal[i];
  }
  double residual = 0.;
  for (int i = 0; i < 12; ++i) residual += h[i] * x[i];
  double innovation = 2. * varianceFloor;
  for (int i = 0; i < 12; ++i)
    for (int j = 0; j < 12; ++j)
      innovation += h[i] * covariance[i][j] * h[j];
  if (!(innovation > 0.) || !std::isfinite(innovation)) return false;

  double covarianceXH[12] = {0.};
  for (int i = 0; i < 12; ++i)
    for (int j = 0; j < 12; ++j)
      covarianceXH[i] += covariance[i][j] * h[j];

  double fitted[12] = {0.};
  for (int i = 0; i < 12; ++i)
    fitted[i] = x[i] - covarianceXH[i] * residual / innovation;

  double fittedCovariance[12][12] = {{0.}};
  for (int i = 0; i < 12; ++i) {
    for (int j = 0; j < 12; ++j) {
      fittedCovariance[i][j] = covariance[i][j] -
        covarianceXH[i] * covarianceXH[j] / innovation;
    }
  }
  for (int i = 0; i < 12; ++i) {
    for (int j = 0; j < i; ++j) {
      const double value = 0.5 * (fittedCovariance[i][j] +
                                  fittedCovariance[j][i]);
      fittedCovariance[i][j] = value;
      fittedCovariance[j][i] = value;
    }
  }

  for (int daughter = 0; daughter < 2; ++daughter) {
    const int offset = 6 * daughter;
    daughterAtVertex[daughter] = propagated[daughter];
    daughterAtVertex[daughter].positionCm = CeeV0Lite::Vec3(
      fitted[offset], fitted[offset + 1], fitted[offset + 2]);
    daughterAtVertex[daughter].momentumGeV = CeeV0Lite::Vec3(
      fitted[offset + 3], fitted[offset + 4], fitted[offset + 5]);
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j)
        daughterAtVertex[daughter].covariance.m[i][j] =
          fittedCovariance[offset + i][offset + j];
    CeeV0Lite::Symmetrize(daughterAtVertex[daughter].covariance);
  }

  // Project the correlated 12D fitted state into the compact mother state
  // y=(SV, p1+p2, E1+E2).  Averaging the two fitted positions keeps the
  // construction invariant under daughter order.
  mother = CeeV0Lite::MotherState();
  mother.positionCm = 0.5 * (
    daughterAtVertex[0].positionCm + daughterAtVertex[1].positionCm);
  mother.momentumGeV = daughterAtVertex[0].momentumGeV +
                       daughterAtVertex[1].momentumGeV;
  const double daughterMomentum[2] = {
    CeeV0Lite::Norm(daughterAtVertex[0].momentumGeV),
    CeeV0Lite::Norm(daughterAtVertex[1].momentumGeV)
  };
  const double daughterEnergy[2] = {
    std::sqrt(daughterMomentum[0] * daughterMomentum[0] +
      propagated[0].massGeV * propagated[0].massGeV),
    std::sqrt(daughterMomentum[1] * daughterMomentum[1] +
      propagated[1].massGeV * propagated[1].massGeV)
  };
  if (!(daughterEnergy[0] > 0.) || !(daughterEnergy[1] > 0.)) return false;
  mother.energyGeV = daughterEnergy[0] + daughterEnergy[1];

  double projection[7][12] = {{0.}};
  for (int i = 0; i < 3; ++i) {
    projection[i][i] = 0.5;
    projection[i][6 + i] = 0.5;
    projection[3 + i][3 + i] = 1.;
    projection[3 + i][9 + i] = 1.;
    projection[6][3 + i] =
      daughterAtVertex[0].momentumGeV[i] / daughterEnergy[0];
    projection[6][9 + i] =
      daughterAtVertex[1].momentumGeV[i] / daughterEnergy[1];
  }
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      for (int a = 0; a < 12; ++a)
        for (int b = 0; b < 12; ++b)
          mother.covariance.m[i][j] += projection[i][a] *
            fittedCovariance[a][b] * projection[j][b];
  CeeV0Lite::Symmetrize(mother.covariance);
  mother.covarianceValid = CeeV0Lite::IsFinite(mother.covariance) &&
                           CeeV0Lite::IsPositiveSemidefinite(
                             mother.covariance) &&
                           CeeV0Lite::IsFinite(mother.positionCm) &&
                           CeeV0Lite::IsFinite(mother.momentumGeV) &&
                           std::isfinite(mother.energyGeV);
  mother.covarianceStatus = mother.covarianceValid
    ? CeeV0Lite::kMotherCovarianceValid
    : CeeV0Lite::kMotherCovarianceNotPositiveSemidefinite;
  if (!mother.covarianceValid) return false;

  chi2 = residual * residual / innovation;
  return std::isfinite(chi2) && chi2 >= 0.;
}
}

namespace CeeV0Lite
{

VertexFitterConfig::VertexFitterConfig()
  : maximumIterations(12),
    convergenceCm(1.e-4),
    minimumOpeningAngleRad(1.e-4),
    requireConvergence(true),
    initialVertexVarianceCm2(1.e6),
    measurementVarianceFloorCm2(1.e-6),
    motherJacobianPositionStepCm(1.e-3),
    motherJacobianMomentumRelativeStep(1.e-4),
    motherJacobianMinimumMomentumStepGeV(1.e-6)
{
}

VertexFitResult::VertexFitResult()
  : valid(false),
    converged(false),
    status(kVertexFitNotRun),
    iterations(0),
    daughterDcaCm(0.),
    chi2(0.),
    ndf(1)
{
}

VertexFitter::VertexFitter(const VertexFitterConfig& fitterConfig,
                           const TransportConfig& transportConfig)
  : fFitterConfig(fitterConfig),
    fTransportConfig(transportConfig)
{
}

bool VertexFitter::Fit(const TrackState& first,
                       const TrackState& second,
                       double magneticFieldByTesla,
                       VertexFitResult& result) const
{
  result = VertexFitResult();
  Vec3 vertex;
  bool degenerate = false;
  if (!ClosestLineSeed(first, second,
                       fFitterConfig.minimumOpeningAngleRad,
                       vertex, degenerate)) {
    result.status = degenerate ? kVertexFitDegenerateGeometry
                               : kVertexFitNumericalFailure;
    return false;
  }

  Mat3 vertexCovariance;
  TrackState propagated[2];
  const TrackState inputs[2] = {first, second};

  for (int iteration = 0; iteration < fFitterConfig.maximumIterations; ++iteration) {
    const Vec3 previous = vertex;
    // Both daughters are linearized at the same frozen vertex.  The batch
    // transverse update is therefore invariant under daughter order.
    for (int i = 0; i < 2; ++i) {
      double path = 0.;
      ClosestApproachResult closestApproach;
      if (!TransportToClosestApproach(inputs[i], previous,
                                      magneticFieldByTesla,
                                      fTransportConfig,
                                      propagated[i], path,
                                      &closestApproach)) {
        result.status = kVertexFitTransportFailure;
        return false;
      }
    }
    if (!BatchUpdateVertex(propagated,
                           fFitterConfig.measurementVarianceFloorCm2,
                           vertex, vertexCovariance)) {
      result.status = kVertexFitNumericalFailure;
      return false;
    }

    result.iterations = iteration + 1;
    if (Norm(vertex - previous) < fFitterConfig.convergenceCm) {
      result.converged = true;
      break;
    }
  }

  if (!result.converged && fFitterConfig.requireConvergence) {
    result.vertexCm = vertex;
    result.covarianceCm2 = vertexCovariance;
    result.status = kVertexFitMaximumIterations;
    return false;
  }

  for (int i = 0; i < 2; ++i) {
    double path = 0.;
    ClosestApproachResult closestApproach;
    if (!TransportToClosestApproach(inputs[i], vertex,
                                    magneticFieldByTesla,
                                    fTransportConfig,
                                    propagated[i], path,
                                    &closestApproach)) {
      result.status = kVertexFitTransportFailure;
      return false;
    }
  }

  result.daughterDcaCm = Norm(propagated[0].positionCm -
                              propagated[1].positionCm);
  if (!ConstructMotherAtCommonVertex(
        propagated, fFitterConfig.measurementVarianceFloorCm2,
        result.daughterAtVertex, result.mother, result.chi2)) {
    result.mother.covarianceStatus = kMotherCovarianceVertexUpdateFailed;
    result.status = kVertexFitNumericalFailure;
    return false;
  }
  result.vertexCm = result.mother.positionCm;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      result.covarianceCm2.m[i][j] = result.mother.covariance.m[i][j];
  result.ndf = 1; // two tracks x two transverse measurements - three vertex parameters
  result.valid = IsFinite(vertex) && std::isfinite(result.chi2);
  if (!result.valid) {
    result.status = kVertexFitNumericalFailure;
    return false;
  }
  result.status = result.converged ? kVertexFitConverged
                                   : kVertexFitMaximumIterations;
  return result.converged || !fFitterConfig.requireConvergence;
}

bool VertexFitter::RecalculateMotherCovariance(
    const TrackState& first,
    const TrackState& second,
    double magneticFieldByTesla,
    MotherState& mother) const
{
  ResetMotherJacobianDiagnostics(mother);
  // J C J^T below starts from the two input track covariances.  Non-zero
  // stochastic transport noise would need extra independent noise modes and
  // must not be silently omitted.
  if (fTransportConfig.positionProcessNoiseCm2PerCm != 0. ||
      fTransportConfig.angularProcessNoiseRad2PerCm != 0.) {
    InvalidateMotherCovariance(
      mother, kMotherCovarianceProcessNoiseUnsupported);
    return false;
  }
  if (!IsPositiveSemidefinite(first.covariance) ||
      !IsPositiveSemidefinite(second.covariance)) {
    InvalidateMotherCovariance(
      mother, kMotherCovarianceNotPositiveSemidefinite);
    return false;
  }
  if (!(fFitterConfig.motherJacobianPositionStepCm > 0.) ||
      !(fFitterConfig.motherJacobianMomentumRelativeStep > 0.) ||
      !(fFitterConfig.motherJacobianMinimumMomentumStepGeV > 0.) ||
      !std::isfinite(fFitterConfig.motherJacobianPositionStepCm) ||
      !std::isfinite(fFitterConfig.motherJacobianMomentumRelativeStep) ||
      !std::isfinite(fFitterConfig.motherJacobianMinimumMomentumStepGeV)) {
    InvalidateMotherCovariance(
      mother, kMotherCovarianceNumericalPropagationFailed);
    return false;
  }

  const TrackState nominal[2] = {first, second};
  double jacobian[7][12] = {{0.}};
  for (int column = 0; column < 12; ++column) {
    const int daughter = column / 6;
    const int parameter = column % 6;
    TrackState plus[2] = {nominal[0], nominal[1]};
    TrackState minus[2] = {nominal[0], nominal[1]};
    const double central = parameter < 3
      ? nominal[daughter].positionCm[parameter]
      : nominal[daughter].momentumGeV[parameter - 3];
    const double step = parameter < 3
      ? fFitterConfig.motherJacobianPositionStepCm
      : std::max(fFitterConfig.motherJacobianMinimumMomentumStepGeV,
          std::fabs(central) *
            fFitterConfig.motherJacobianMomentumRelativeStep);
    if (parameter < 3) {
      plus[daughter].positionCm[parameter] += step;
      minus[daughter].positionCm[parameter] -= step;
    }
    else {
      plus[daughter].momentumGeV[parameter - 3] += step;
      minus[daughter].momentumGeV[parameter - 3] -= step;
    }

    VertexFitResult plusFit;
    VertexFitResult minusFit;
    ++mother.jacobianFitsAttempted;
    const bool plusOk = Fit(plus[0], plus[1], magneticFieldByTesla, plusFit);
    if (plusOk) {
      ++mother.jacobianFitsSucceeded;
    }
    else {
      ++mother.jacobianFailureCount;
      mother.jacobianPlusFailureMask |= (1u << column);
      if (mother.jacobianFailedColumn < 0) {
        mother.jacobianFailedColumn = column;
        mother.jacobianFailedSide = +1;
        mother.jacobianFirstFailureFitStatus =
          static_cast<int>(plusFit.status);
      }
    }
    ++mother.jacobianFitsAttempted;
    const bool minusOk = Fit(minus[0], minus[1], magneticFieldByTesla,
                             minusFit);
    if (minusOk) {
      ++mother.jacobianFitsSucceeded;
    }
    else {
      ++mother.jacobianFailureCount;
      mother.jacobianMinusFailureMask |= (1u << column);
      if (mother.jacobianFailedColumn < 0) {
        mother.jacobianFailedColumn = column;
        mother.jacobianFailedSide = -1;
        mother.jacobianFirstFailureFitStatus =
          static_cast<int>(minusFit.status);
      }
    }
    if (!plusOk || !minusOk) continue;
    double plusValues[7] = {
      plusFit.mother.positionCm[0], plusFit.mother.positionCm[1],
      plusFit.mother.positionCm[2], plusFit.mother.momentumGeV[0],
      plusFit.mother.momentumGeV[1], plusFit.mother.momentumGeV[2],
      plusFit.mother.energyGeV
    };
    double minusValues[7] = {
      minusFit.mother.positionCm[0], minusFit.mother.positionCm[1],
      minusFit.mother.positionCm[2], minusFit.mother.momentumGeV[0],
      minusFit.mother.momentumGeV[1], minusFit.mother.momentumGeV[2],
      minusFit.mother.energyGeV
    };
    for (int row = 0; row < 7; ++row)
      jacobian[row][column] =
        (plusValues[row] - minusValues[row]) / (2. * step);
  }

  if (mother.jacobianFailureCount != 0) {
    InvalidateMotherCovariance(
      mother, kMotherCovarianceJacobianPerturbationFailed);
    return false;
  }

  double inputCovariance[12][12] = {{0.}};
  for (int daughter = 0; daughter < 2; ++daughter) {
    const int offset = 6 * daughter;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        const double value = nominal[daughter].covariance.m[i][j];
        if (!std::isfinite(value)) {
          InvalidateMotherCovariance(
            mother, kMotherCovarianceNumericalPropagationFailed);
          return false;
        }
        inputCovariance[offset + i][offset + j] = value;
      }
    }
    // Treat the configured transverse measurement floor as an independent
    // per-track position floor.  The complete POCA rerun eliminates its
    // unobservable along-track component in the output Jacobian.
    for (int i = 0; i < 3; ++i)
      inputCovariance[offset + i][offset + i] +=
        fFitterConfig.measurementVarianceFloorCm2;
  }

  Mat7 propagatedCovariance;
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      for (int a = 0; a < 12; ++a)
        for (int b = 0; b < 12; ++b)
          propagatedCovariance.m[i][j] += jacobian[i][a] *
            inputCovariance[a][b] * jacobian[j][b];
  Symmetrize(propagatedCovariance);
  if (!IsFinite(propagatedCovariance)) {
    InvalidateMotherCovariance(mother, kMotherCovarianceNonFinite);
    return false;
  }
  if (!IsPositiveSemidefinite(propagatedCovariance)) {
    InvalidateMotherCovariance(
      mother, kMotherCovarianceNotPositiveSemidefinite);
    return false;
  }
  mother.covariance = propagatedCovariance;
  mother.covarianceValid = true;
  mother.covarianceStatus = kMotherCovarianceValid;
  return true;
}

} // namespace CeeV0Lite
