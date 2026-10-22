#include "CeeV0LiteTypes.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
const double kSpeedOfLightCmPerNs = 29.9792458;
}

namespace CeeV0Lite
{

TrackInput::TrackInput()
  : id(-1),
    sigmaDirectionRad(0.003),
    lambdaRoleDcaBoundaryCm(1.2),
    hasTrackingMomentum(false),
    hasTrackingCovariance(false),
    chargeKnown(false),
    charge(0),
    hasTof(false),
    tofTimeNs(0.),
    sigmaTofTimeNs(0.10),
    sigmaPathLengthCm(0.5),
    hypothesisMask(kPionHypothesis | kProtonHypothesis),
    qualityChi2(0.),
    qualityNdf(1)
{
  positionCovarianceCm2.m[0][0] = 0.01;
  positionCovarianceCm2.m[1][1] = 0.01;
  positionCovarianceCm2.m[2][2] = 0.0225;
}

PrimaryVertex::PrimaryVertex()
  : chi2(0.), nContributors(1)
{
  covarianceCm2.m[0][0] = 0.01;
  covarianceCm2.m[1][1] = 0.01;
  covarianceCm2.m[2][2] = 0.01;
}

EventInput::EventInput()
  : startTimeNs(0.), sigmaStartTimeNs(0.)
{
}

TrackState::TrackState()
  : id(-1), pdg(0), charge(0), massGeV(0.)
{
}

MotherState::MotherState()
  : energyGeV(0.),
    covarianceValid(false),
    covarianceStatus(kMotherCovarianceNotCalculated),
    jacobianFitsAttempted(0),
    jacobianFitsSucceeded(0),
    jacobianFailureCount(0),
    jacobianPlusFailureMask(0u),
    jacobianMinusFailureMask(0u),
    jacobianFailedColumn(-1),
    jacobianFailedSide(0),
    jacobianFirstFailureFitStatus(0)
{
}

TofMomentumEstimate::TofMomentumEstimate()
  : valid(false),
    beta(0.),
    momentumGeV(0.),
    sigmaMomentumGeV(0.),
    pathLengthCm(0.),
    flightTimeNs(0.)
{
}

Candidate::Candidate()
  : type(kK0Short),
    massGeV(0.),
    sigmaMassGeV(std::numeric_limits<double>::quiet_NaN()),
    daughterDcaCm(0.),
    geometricChi2(0.),
    geometricNdf(1),
    geometricChi2Ndf(0.),
    topologicalChi2Ndf(0.),
    topologicalChi2(0.),
    topologicalNdf(3),
    primaryVertexConstraintChi2(0.),
    primaryVertexConstraintNdf(2),
    primaryVertexConstraintChi2Ndf(0.),
    decayLengthCm(0.),
    sigmaDecayLengthCm(0.),
    decayLengthSignificance(0.),
    pointingCosine(-1.),
    decayTimeNs(0.),
    tofIterations(0),
    tofConverged(false),
    motherJacobianComplete(false),
    covarianceDependentCutsEvaluated(false)
{
  daughterId[0] = daughterId[1] = -1;
  daughterPdg[0] = daughterPdg[1] = 0;
  daughterDcaToPvCm[0] = daughterDcaToPvCm[1] = 0.;
  daughterChi2FromPv[0] = daughterChi2FromPv[1] = 0.;
  daughterPvSignificance[0] = daughterPvSignificance[1] = 0.;
  daughterDeviationFromPv[0] = daughterDeviationFromPv[1] = 0.;
}

TofMomentumEstimate EstimateTofMomentum(double massGeV,
                                        double pathLengthCm,
                                        double flightTimeNs,
                                        double sigmaPathLengthCm,
                                        double sigmaFlightTimeNs,
                                        double maximumBeta)
{
  TofMomentumEstimate out;
  out.pathLengthCm = pathLengthCm;
  out.flightTimeNs = flightTimeNs;
  if (!(massGeV > 0.) || !(pathLengthCm > 0.) || !(flightTimeNs > 0.)) return out;

  const double beta = pathLengthCm / (kSpeedOfLightCmPerNs * flightTimeNs);
  if (!(beta > 0.) || !(beta < maximumBeta) || !(maximumBeta < 1.)) return out;

  const double gamma = 1. / std::sqrt(1. - beta * beta);
  const double momentum = massGeV * beta * gamma;
  const double dBetaDLength = 1. / (kSpeedOfLightCmPerNs * flightTimeNs);
  const double dBetaDTime = -pathLengthCm /
    (kSpeedOfLightCmPerNs * flightTimeNs * flightTimeNs);
  const double sigmaBeta2 = dBetaDLength * dBetaDLength *
                              sigmaPathLengthCm * sigmaPathLengthCm +
                            dBetaDTime * dBetaDTime *
                              sigmaFlightTimeNs * sigmaFlightTimeNs;
  const double sigmaMomentum = massGeV * gamma * gamma * gamma *
                               std::sqrt(std::max(0., sigmaBeta2));

  out.valid = std::isfinite(momentum) && std::isfinite(sigmaMomentum) &&
              momentum > 0.;
  out.beta = beta;
  out.momentumGeV = momentum;
  out.sigmaMomentumGeV = sigmaMomentum;
  return out;
}

Mat6 BuildTofStateCovariance(const TrackInput& track,
                             const Vec3& unitDirection,
                             double momentumGeV,
                             double sigmaMomentumGeV)
{
  Mat6 covariance;
  SetPositionBlock(covariance, track.positionCovarianceCm2);

  Mat3 momentumCovariance;
  const double sigmaP2 = sigmaMomentumGeV * sigmaMomentumGeV;
  const double angularVariance = momentumGeV * momentumGeV *
                                 track.sigmaDirectionRad *
                                 track.sigmaDirectionRad;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double uu = unitDirection[i] * unitDirection[j];
      momentumCovariance.m[i][j] = sigmaP2 * uu +
        angularVariance * (((i == j) ? 1. : 0.) - uu);
    }
  }
  SetMomentumBlock(covariance, momentumCovariance);
  return covariance;
}

bool MotherMassAndError(const MotherState& mother,
                        double& massGeV,
                        double& sigmaMassGeV,
                        double roundoffToleranceMultiplier)
{
  const double mass2 = mother.energyGeV * mother.energyGeV -
                       Norm2(mother.momentumGeV);
  if (!(mass2 > 0.) || !std::isfinite(mass2)) return false;
  massGeV = std::sqrt(mass2);
  if (!mother.covarianceValid || !IsFinite(mother.covariance) ||
      !IsPositiveSemidefinite(mother.covariance) ||
      !(roundoffToleranceMultiplier >= 0.) ||
      !std::isfinite(roundoffToleranceMultiplier)) return false;

  // Gradient of M with respect to (x,y,z,px,py,pz,E).
  double gradient[7] = {0., 0., 0.,
    -mother.momentumGeV[0] / massGeV,
    -mother.momentumGeV[1] / massGeV,
    -mother.momentumGeV[2] / massGeV,
     mother.energyGeV / massGeV};
  double variance = QuadraticForm(gradient, mother.covariance);
  if (!std::isfinite(variance)) return false;
  double absoluteSum = 0.;
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      absoluteSum += std::fabs(gradient[i] * mother.covariance.m[i][j] *
                               gradient[j]);
  const double tolerance = roundoffToleranceMultiplier *
    std::numeric_limits<double>::epsilon() * absoluteSum;
  if (variance < -tolerance) return false;
  if (variance < 0.) variance = 0.;
  sigmaMassGeV = std::sqrt(variance);
  return std::isfinite(massGeV) && std::isfinite(sigmaMassGeV);
}

bool AddCommonNuisanceCovariance(MotherState& mother,
                                 const double derivative[7],
                                 double variance)
{
  if (!mother.covarianceValid || !(variance >= 0.) ||
      !std::isfinite(variance)) return false;
  for (int i = 0; i < 7; ++i) {
    if (!std::isfinite(derivative[i])) return false;
    for (int j = 0; j < 7; ++j)
      mother.covariance.m[i][j] += derivative[i] * variance * derivative[j];
  }
  Symmetrize(mother.covariance);
  if (!IsFinite(mother.covariance)) {
    mother.covarianceValid = false;
    mother.covarianceStatus = kMotherCovarianceNonFinite;
    return false;
  }
  if (!IsPositiveSemidefinite(mother.covariance)) {
    mother.covarianceValid = false;
    mother.covarianceStatus = kMotherCovarianceNotPositiveSemidefinite;
    return false;
  }
  return true;
}

} // namespace CeeV0Lite
