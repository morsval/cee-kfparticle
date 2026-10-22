#include "CeeV0LiteFinder.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
const double kPionMassGeV = 0.13957039;
const double kProtonMassGeV = 0.9382720813;
const double kSpeedOfLightCmPerNs = 29.9792458;

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

bool HasHypothesis(const CeeV0Lite::TrackInput& track, unsigned int hypothesis)
{
  return (track.hypothesisMask & hypothesis) != 0u;
}

bool OrientedDirection(const CeeV0Lite::TrackInput& track,
                       CeeV0Lite::Vec3& direction)
{
  if (!CeeV0Lite::Normalize(track.direction, direction)) return false;
  if (track.hasTof &&
      CeeV0Lite::Dot(track.tofPositionCm - track.positionCm, direction) < 0.)
    direction = -1. * direction;
  return true;
}

bool RawLineSeed(const CeeV0Lite::TrackInput& first,
                 const CeeV0Lite::TrackInput& second,
                 double minimumOpeningAngleRad,
                 CeeV0Lite::Vec3& seed)
{
  CeeV0Lite::Vec3 u;
  CeeV0Lite::Vec3 v;
  if (!OrientedDirection(first, u) || !OrientedDirection(second, v)) return false;
  const CeeV0Lite::Vec3 w = first.positionCm - second.positionCm;
  const double uv = CeeV0Lite::Dot(u, v);
  const double uw = CeeV0Lite::Dot(u, w);
  const double vw = CeeV0Lite::Dot(v, w);
  const double denominator = 1. - uv * uv;
  const double minimumSine = std::sin(minimumOpeningAngleRad);
  if (minimumOpeningAngleRad < 0. || !std::isfinite(minimumSine) ||
      minimumOpeningAngleRad > 1.5707963267948966 ||
      denominator < minimumSine * minimumSine) return false;
  const double firstPath = (uv * vw - uw) / denominator;
  const double secondPath = (vw - uv * uw) / denominator;
  seed = 0.5 * ((first.positionCm + firstPath * u) +
                (second.positionCm + secondPath * v));
  return CeeV0Lite::IsFinite(seed);
}

bool UseTrackingMomentum(const CeeV0Lite::TrackInput& track,
                         const CeeV0Lite::FinderConfig& config)
{
  if (config.momentumMode == CeeV0Lite::kTofMomentum) return false;
  return track.hasTrackingMomentum;
}

bool ChargesAllow(const CeeV0Lite::TrackInput& first,
                  int expectedFirst,
                  const CeeV0Lite::TrackInput& second,
                  int expectedSecond,
                  const CeeV0Lite::FinderConfig& config)
{
  const bool fieldOn = std::fabs(config.magneticFieldByTesla) >
                       config.zeroFieldThresholdTesla;
  if (fieldOn && config.requireMeasuredChargeInMagneticField &&
      (!first.chargeKnown || !second.chargeKnown)) return false;
  if (first.chargeKnown && first.charge != expectedFirst) return false;
  if (second.chargeKnown && second.charge != expectedSecond) return false;
  return true;
}

bool BuildTrackState(const CeeV0Lite::TrackInput& track,
                     int pdg,
                     int expectedCharge,
                     double massGeV,
                     const CeeV0Lite::Vec3& decayVertex,
                     double decayTimeNs,
                     const CeeV0Lite::FinderConfig& config,
                     CeeV0Lite::TrackState& state)
{
  state = CeeV0Lite::TrackState();
  state.id = track.id;
  state.pdg = pdg;
  state.massGeV = massGeV;
  state.positionCm = track.positionCm;
  state.charge = track.chargeKnown ? track.charge : expectedCharge;

  const bool fieldOn = std::fabs(config.magneticFieldByTesla) >
                       config.zeroFieldThresholdTesla;
  const bool useTracking = UseTrackingMomentum(track, config);
  if (config.momentumMode == CeeV0Lite::kTrackingMomentum &&
      !track.hasTrackingMomentum) return false;
  if (fieldOn && !useTracking && !config.allowTofMomentumInMagneticField)
    return false;

  if (useTracking) {
    const double momentum = CeeV0Lite::Norm(track.trackingMomentumGeV);
    if (!(momentum > 1.e-10)) return false;
    state.momentumGeV = track.trackingMomentumGeV;
    if (track.hasTrackingCovariance) {
      state.covariance = track.trackingCovariance;
    }
    else {
      CeeV0Lite::Vec3 direction;
      CeeV0Lite::Normalize(track.trackingMomentumGeV, direction);
      const double sigmaP = std::max(1.e-6,
        config.defaultRelativeTrackingMomentumError * momentum);
      state.covariance = CeeV0Lite::BuildTofStateCovariance(
        track, direction, momentum, sigmaP);
    }
    return true;
  }

  if (!track.hasTof) return false;
  CeeV0Lite::Vec3 direction;
  if (!OrientedDirection(track, direction)) return false;
  const double pathLength = CeeV0Lite::Dot(track.tofPositionCm - decayVertex,
                                           direction);
  const double flightTime = track.tofTimeNs - decayTimeNs;
  // Only per-track timing error belongs in this marginal covariance.  T0 and
  // the common decay-time model error are shared by both daughters and are
  // propagated once into the mother covariance after the fixed-point solve.
  const double sigmaFlightTime = track.sigmaTofTimeNs;
  const CeeV0Lite::TofMomentumEstimate estimate =
    CeeV0Lite::EstimateTofMomentum(massGeV, pathLength, flightTime,
                                   track.sigmaPathLengthCm,
                                   sigmaFlightTime,
                                   config.maximumBeta);
  if (!estimate.valid) return false;
  state.momentumGeV = estimate.momentumGeV * direction;
  state.covariance = CeeV0Lite::BuildTofStateCovariance(
    track, direction, estimate.momentumGeV, estimate.sigmaMomentumGeV);
  return true;
}

double UpdatedDecayTime(const CeeV0Lite::EventInput& event,
                        const CeeV0Lite::VertexFitResult& fit)
{
  const CeeV0Lite::Vec3 parentMomentum = fit.mother.momentumGeV;
  const double momentum = CeeV0Lite::Norm(parentMomentum);
  if (!(momentum > 1.e-12)) return std::numeric_limits<double>::quiet_NaN();
  const double beta = momentum / fit.mother.energyGeV;
  if (!(beta > 0.) || !(beta < 1.))
    return std::numeric_limits<double>::quiet_NaN();
  const CeeV0Lite::Vec3 direction = parentMomentum / momentum;
  const double signedLength = CeeV0Lite::Dot(
    fit.vertexCm - event.primaryVertex.positionCm, direction);
  if (!(signedLength >= 0.)) return std::numeric_limits<double>::quiet_NaN();
  return event.startTimeNs + signedLength / (beta * kSpeedOfLightCmPerNs);
}

double PrimaryVertexConstraintChi2(const CeeV0Lite::Candidate& candidate,
                                   const CeeV0Lite::EventInput& event,
                                   double varianceFloor)
{
  CeeV0Lite::Vec3 direction;
  if (!CeeV0Lite::Normalize(candidate.momentumGeV, direction))
    return std::numeric_limits<double>::infinity();
  CeeV0Lite::Vec3 e1;
  CeeV0Lite::Vec3 e2;
  CeeV0Lite::MakeTransverseBasis(direction, e1, e2);
  const CeeV0Lite::Vec3 displacement = candidate.decayVertexCm -
                                       event.primaryVertex.positionCm;
  const double residual[2] = {
    CeeV0Lite::Dot(e1, displacement), CeeV0Lite::Dot(e2, displacement)
  };
  const CeeV0Lite::Vec3 basis[2] = {e1, e2};
  const double momentum = CeeV0Lite::Norm(candidate.momentumGeV);
  if (!(momentum > 0.) || !candidate.mother.covarianceValid)
    return std::numeric_limits<double>::infinity();
  const double signedLength = CeeV0Lite::Dot(direction, displacement);
  double jacobian[2][7] = {{0.}};
  for (int a = 0; a < 2; ++a) {
    for (int i = 0; i < 3; ++i) {
      jacobian[a][i] = basis[a][i];
      // f=e_perp dot (SV-PV) with the parent direction allowed to vary:
      // df/dP = -(longitudinal length / |P|) e_perp.
      jacobian[a][3 + i] = -signedLength * basis[a][i] / momentum;
    }
  }
  double projected[2][2] = {{0.}};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int a = 0; a < 7; ++a)
        for (int b = 0; b < 7; ++b)
          projected[i][j] += jacobian[i][a] *
            candidate.mother.covariance.m[a][b] * jacobian[j][b];
      projected[i][j] += CeeV0Lite::Dot(basis[i],
        CeeV0Lite::Multiply(event.primaryVertex.covarianceCm2, basis[j]));
    }
  }
  projected[0][0] += varianceFloor;
  projected[1][1] += varianceFloor;
  double inverse[2][2];
  if (!InvertPositiveDefinite2(projected, inverse))
    return std::numeric_limits<double>::infinity();
  const double chi2 = residual[0] *
    (inverse[0][0] * residual[0] + inverse[0][1] * residual[1]) +
    residual[1] *
    (inverse[1][0] * residual[0] + inverse[1][1] * residual[1]);
  return std::isfinite(chi2) && chi2 >= 0.
    ? chi2 : std::numeric_limits<double>::infinity();
}

struct DaughterPvMetrics
{
  DaughterPvMetrics()
    : valid(false), dcaCm(0.), chi2(0.), significance(0.) {}
  bool valid;
  double dcaCm;
  double chi2;
  double significance;
};

DaughterPvMetrics MeasureDaughterAtPv(
    const CeeV0Lite::TrackState& state,
    const CeeV0Lite::EventInput& event,
    const CeeV0Lite::FinderConfig& config)
{
  DaughterPvMetrics metrics;
  CeeV0Lite::TrackState atPv;
  double path = 0.;
  CeeV0Lite::ClosestApproachResult closestApproach;
  if (!CeeV0Lite::TransportToClosestApproach(
        state, event.primaryVertex.positionCm,
        config.magneticFieldByTesla, config.transportConfig,
        atPv, path, &closestApproach)) return metrics;
  CeeV0Lite::Vec3 direction;
  if (!CeeV0Lite::Normalize(atPv.momentumGeV, direction)) return metrics;
  CeeV0Lite::Vec3 e1;
  CeeV0Lite::Vec3 e2;
  CeeV0Lite::MakeTransverseBasis(direction, e1, e2);
  const CeeV0Lite::Mat3 covariance =
    CeeV0Lite::PositionBlock(atPv.covariance) +
    event.primaryVertex.covarianceCm2;
  const CeeV0Lite::Vec3 residual3 = atPv.positionCm -
                                    event.primaryVertex.positionCm;
  metrics.dcaCm = CeeV0Lite::Norm(residual3);
  const CeeV0Lite::Vec3 basis[2] = {e1, e2};
  const double residual[2] = {
    CeeV0Lite::Dot(e1, residual3), CeeV0Lite::Dot(e2, residual3)
  };
  double projected[2][2] = {{0.}};
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      projected[i][j] = CeeV0Lite::Dot(basis[i],
        CeeV0Lite::Multiply(covariance, basis[j]));
  projected[0][0] += config.vertexFitterConfig.measurementVarianceFloorCm2;
  projected[1][1] += config.vertexFitterConfig.measurementVarianceFloorCm2;
  double inverse[2][2];
  if (!InvertPositiveDefinite2(projected, inverse)) return metrics;
  const double chi2 = residual[0] *
    (inverse[0][0] * residual[0] + inverse[0][1] * residual[1]) +
    residual[1] *
    (inverse[1][0] * residual[0] + inverse[1][1] * residual[1]);
  if (!std::isfinite(chi2) || chi2 < 0.) return metrics;
  // Match KFParticle::GetDeviationFromVertex(): the historical STAR method
  // name says "deviation", but the returned and cut quantity is raw chi2.
  metrics.chi2 = chi2;
  metrics.significance = std::sqrt(chi2);
  metrics.valid = std::isfinite(metrics.dcaCm) &&
                  std::isfinite(metrics.chi2) &&
                  std::isfinite(metrics.significance);
  return metrics;
}

bool PassCuts(const CeeV0Lite::Candidate& candidate,
              const CeeV0Lite::ChannelCuts& cuts)
{
  if (!cuts.enabled) return false;
  const bool passesCentralCuts =
    candidate.massGeV >= cuts.minimumMassGeV &&
    candidate.massGeV <= cuts.maximumMassGeV &&
    candidate.daughterDcaCm <= cuts.maximumDaughterDcaCm &&
    candidate.geometricChi2Ndf <= cuts.maximumGeometricChi2Ndf &&
    candidate.decayLengthCm >= cuts.minimumDecayLengthCm &&
    candidate.pointingCosine >= cuts.minimumPointingCosine &&
    candidate.daughterChi2FromPv[0] >= cuts.minimumDaughterChi2FromPv &&
    candidate.daughterChi2FromPv[1] >= cuts.minimumDaughterChi2FromPv &&
    candidate.daughterPvSignificance[0] >=
      cuts.minimumDaughterDeviationFromPv &&
    candidate.daughterPvSignificance[1] >=
      cuts.minimumDaughterDeviationFromPv;
  if (!passesCentralCuts) return false;

  if (!candidate.covarianceDependentCutsEvaluated) {
    return !cuts.requireCompleteMotherJacobian &&
      !candidate.motherJacobianComplete &&
      candidate.mother.covarianceStatus ==
        CeeV0Lite::kMotherCovarianceJacobianPerturbationFailed;
  }
  return candidate.motherJacobianComplete &&
         candidate.topologicalChi2Ndf <= cuts.maximumTopologicalChi2Ndf &&
         candidate.decayLengthSignificance >=
           cuts.minimumDecayLengthSignificance;
}

void InvalidateCovarianceDependentObservables(
    CeeV0Lite::Candidate& candidate)
{
  const double invalid = std::numeric_limits<double>::quiet_NaN();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      candidate.decayVertexCovarianceCm2.m[i][j] = invalid;
  candidate.sigmaMassGeV = invalid;
  candidate.sigmaDecayLengthCm = invalid;
  candidate.decayLengthSignificance = invalid;
  candidate.primaryVertexConstraintChi2 = invalid;
  candidate.primaryVertexConstraintChi2Ndf = invalid;
  candidate.topologicalChi2 = invalid;
  candidate.topologicalChi2Ndf = invalid;
  candidate.covarianceDependentCutsEvaluated = false;
}

struct PairSolution
{
  CeeV0Lite::VertexFitResult fit;
  CeeV0Lite::TrackState inputStates[2];
  double decayTimeNs;
  int iterations;
  bool converged;
  bool anyTof;
};

bool SolvePair(const CeeV0Lite::TrackInput& firstTrack,
               int firstPdg,
               int firstCharge,
               double firstMass,
               const CeeV0Lite::TrackInput& secondTrack,
               int secondPdg,
               int secondCharge,
               double secondMass,
               const CeeV0Lite::EventInput& event,
               const CeeV0Lite::FinderConfig& config,
               PairSolution& solution)
{
  CeeV0Lite::Vec3 decayVertex;
  if (!RawLineSeed(firstTrack, secondTrack,
                   config.vertexFitterConfig.minimumOpeningAngleRad,
                   decayVertex)) return false;
  const double seedBeta = std::max(0.05, std::min(0.999,
                                                  config.seedParentBeta));
  double decayTime = event.startTimeNs +
    CeeV0Lite::Norm(decayVertex - event.primaryVertex.positionCm) /
    (seedBeta * kSpeedOfLightCmPerNs);

  solution.anyTof = !UseTrackingMomentum(firstTrack, config) ||
                    !UseTrackingMomentum(secondTrack, config);
  const int maximumIterations = solution.anyTof
    ? std::max(1, config.maximumTofIterations) : 1;
  CeeV0Lite::VertexFitter fitter(config.vertexFitterConfig,
                                 config.transportConfig);
  solution.converged = !solution.anyTof;
  solution.iterations = 0;
  solution.decayTimeNs = decayTime;

  for (int iteration = 0; iteration < maximumIterations; ++iteration) {
    if (!BuildTrackState(firstTrack, firstPdg, firstCharge, firstMass,
                         decayVertex, decayTime, config,
                         solution.inputStates[0]) ||
        !BuildTrackState(secondTrack, secondPdg, secondCharge, secondMass,
                         decayVertex, decayTime, config,
                         solution.inputStates[1])) return false;
    if (!fitter.Fit(solution.inputStates[0], solution.inputStates[1],
                    config.magneticFieldByTesla, solution.fit)) return false;
    ++solution.iterations;
    if (!solution.anyTof) break;

    const double updatedTime = UpdatedDecayTime(event, solution.fit);
    if (!std::isfinite(updatedTime)) return false;
    const double vertexShift = CeeV0Lite::Norm(solution.fit.vertexCm -
                                               decayVertex);
    const double timeShift = std::fabs(updatedTime - decayTime);
    decayVertex = solution.fit.vertexCm;
    decayTime = updatedTime;
    solution.decayTimeNs = updatedTime;
    if (vertexShift < config.tofVertexConvergenceCm &&
        timeShift < config.tofTimeConvergenceNs) {
      solution.converged = true;
      break;
    }
  }
  if (solution.anyTof && config.requireTofConvergence &&
      !solution.converged) return false;
  if (solution.anyTof) {
    // The convergence test above compares the old fit to an updated
    // (vertex,time).  Refit once at that final point so the returned mother,
    // daughter states and reported decayTime all belong to the same
    // fixed-point iterate.
    if (!BuildTrackState(firstTrack, firstPdg, firstCharge, firstMass,
                         decayVertex, decayTime, config,
                         solution.inputStates[0]) ||
        !BuildTrackState(secondTrack, secondPdg, secondCharge, secondMass,
                         decayVertex, decayTime, config,
                         solution.inputStates[1]) ||
        !fitter.Fit(solution.inputStates[0], solution.inputStates[1],
                    config.magneticFieldByTesla, solution.fit)) return false;
    solution.decayTimeNs = decayTime;
    if (solution.converged) {
      const double checkTime = UpdatedDecayTime(event, solution.fit);
      if (!std::isfinite(checkTime) ||
          CeeV0Lite::Norm(solution.fit.vertexCm - decayVertex) >=
            config.tofVertexConvergenceCm ||
          std::fabs(checkTime - decayTime) >=
            config.tofTimeConvergenceNs) return false;
    }
  }
  else {
    const double fittedDecayTime = UpdatedDecayTime(event, solution.fit);
    if (std::isfinite(fittedDecayTime))
      solution.decayTimeNs = fittedDecayTime;
  }
  return solution.fit.valid && solution.fit.mother.covarianceValid;
}

void MotherParameters(const CeeV0Lite::MotherState& mother, double values[7])
{
  for (int i = 0; i < 3; ++i) {
    values[i] = mother.positionCm[i];
    values[i + 3] = mother.momentumGeV[i];
  }
  values[6] = mother.energyGeV;
}

bool AddSharedTimeUncertainty(
    const CeeV0Lite::TrackInput& firstTrack,
    int firstPdg,
    int firstCharge,
    double firstMass,
    const CeeV0Lite::TrackInput& secondTrack,
    int secondPdg,
    int secondCharge,
    double secondMass,
    const CeeV0Lite::EventInput& event,
    const CeeV0Lite::FinderConfig& config,
    const PairSolution& nominal,
    CeeV0Lite::MotherState& mother)
{
  if (!config.propagateSharedTimeUncertainty || !nominal.anyTof) return true;
  if (!(event.sigmaStartTimeNs >= 0.) ||
      !(config.additionalDecayTimeSigmaNs >= 0.)) return false;
  const double variance = event.sigmaStartTimeNs * event.sigmaStartTimeNs +
    config.additionalDecayTimeSigmaNs * config.additionalDecayTimeSigmaNs;
  if (variance == 0.) return true;
  if (!(variance > 0.) || !std::isfinite(variance) ||
      !(config.sharedTimeDerivativeStepNs > 0.) ||
      !std::isfinite(config.sharedTimeDerivativeStepNs) ||
      !(config.sharedTimeDerivativeStabilityTolerance >= 0.) ||
      !std::isfinite(config.sharedTimeDerivativeStabilityTolerance))
    return false;

  CeeV0Lite::EventInput plusEvent = event;
  CeeV0Lite::EventInput minusEvent = event;
  plusEvent.startTimeNs += config.sharedTimeDerivativeStepNs;
  minusEvent.startTimeNs -= config.sharedTimeDerivativeStepNs;
  PairSolution plus;
  PairSolution minus;
  if (!SolvePair(firstTrack, firstPdg, firstCharge, firstMass,
                 secondTrack, secondPdg, secondCharge, secondMass,
                 plusEvent, config, plus) ||
      !SolvePair(firstTrack, firstPdg, firstCharge, firstMass,
                 secondTrack, secondPdg, secondCharge, secondMass,
                 minusEvent, config, minus)) {
    mother.covarianceValid = false;
    mother.covarianceStatus =
      CeeV0Lite::kMotherCovarianceSharedTimeDerivativeFailed;
    return false;
  }

  CeeV0Lite::EventInput halfPlusEvent = event;
  CeeV0Lite::EventInput halfMinusEvent = event;
  halfPlusEvent.startTimeNs += 0.5 * config.sharedTimeDerivativeStepNs;
  halfMinusEvent.startTimeNs -= 0.5 * config.sharedTimeDerivativeStepNs;
  PairSolution halfPlus;
  PairSolution halfMinus;
  if (!SolvePair(firstTrack, firstPdg, firstCharge, firstMass,
                 secondTrack, secondPdg, secondCharge, secondMass,
                 halfPlusEvent, config, halfPlus) ||
      !SolvePair(firstTrack, firstPdg, firstCharge, firstMass,
                 secondTrack, secondPdg, secondCharge, secondMass,
                 halfMinusEvent, config, halfMinus)) {
    mother.covarianceValid = false;
    mother.covarianceStatus =
      CeeV0Lite::kMotherCovarianceSharedTimeDerivativeFailed;
    return false;
  }

  double plusValues[7];
  double minusValues[7];
  double halfPlusValues[7];
  double halfMinusValues[7];
  MotherParameters(plus.fit.mother, plusValues);
  MotherParameters(minus.fit.mother, minusValues);
  MotherParameters(halfPlus.fit.mother, halfPlusValues);
  MotherParameters(halfMinus.fit.mother, halfMinusValues);
  double derivative[7] = {0.};
  const double denominator = 2. * config.sharedTimeDerivativeStepNs;
  for (int i = 0; i < 7; ++i) {
    const double fullStep = (plusValues[i] - minusValues[i]) / denominator;
    const double halfStep = (halfPlusValues[i] - halfMinusValues[i]) /
                            config.sharedTimeDerivativeStepNs;
    const double scale = std::max(1.e-6,
      std::max(std::fabs(fullStep), std::fabs(halfStep)));
    if (!std::isfinite(fullStep) || !std::isfinite(halfStep) ||
        std::fabs(fullStep - halfStep) >
          config.sharedTimeDerivativeStabilityTolerance * scale) {
      mother.covarianceValid = false;
      mother.covarianceStatus =
        CeeV0Lite::kMotherCovarianceSharedTimeDerivativeFailed;
      return false;
    }
    // The half-step estimate has the smaller O(h^2) truncation error.
    derivative[i] = halfStep;
  }
  if (!CeeV0Lite::AddCommonNuisanceCovariance(mother, derivative, variance)) {
    mother.covarianceValid = false;
    mother.covarianceStatus =
      CeeV0Lite::kMotherCovarianceSharedTimeDerivativeFailed;
    return false;
  }
  mother.covarianceStatus = CeeV0Lite::kMotherCovarianceValid;
  return true;
}

void BuildProductionCopy(
    const CeeV0Lite::TrackInput& firstTrack, int firstPdg, int firstCharge, double firstMass,
    const CeeV0Lite::TrackInput& secondTrack, int secondPdg, int secondCharge, double secondMass,
    const CeeV0Lite::EventInput& event, const CeeV0Lite::FinderConfig& config,
    const PairSolution& nominal, CeeV0Lite::Candidate& candidate)
{
  using namespace CeeV0Lite;
  if (!config.buildProductionState) {
    candidate.motherTopo.status=kProductionDisabled; return;
  }
  if (!candidate.mother.covarianceValid) {
    candidate.motherTopo.status=kProductionInvalidMother; return;
  }
  ProductionVertexInput input;
  input.mother=candidate.mother;
  input.vertex=event.primaryVertex;
  const ProductionCorrelationModel model=nominal.anyTof
    ? kProductionTofPvResponse : kProductionTrackingIndependentPv;
  candidate.motherTopo.correlationModel=model;
  if (nominal.anyTof) {
    const double step=config.productionVertexConfig.pvDerivativeStepCm;
    const double tolerance=config.productionVertexConfig.pvDerivativeStabilityTolerance;
    if (!nominal.converged || !(step>0.) || !std::isfinite(step) ||
        !(tolerance>=0.) || !std::isfinite(tolerance)) {
      candidate.motherTopo.status=kProductionPvResponseFailed; return;
    }
    double response[7][3]={{0.}};
    for (int column=0; column<3; ++column) {
      // A fixed (exact) PV coordinate has no uncertainty to propagate.
      if (event.primaryVertex.covarianceCm2.m[column][column]==0.) continue;
      double derivative[2][7]={{0.}};
      for (int scale=0; scale<2; ++scale) {
        const double h=step/(scale+1.);
        EventInput plusEvent=event, minusEvent=event;
        plusEvent.primaryVertex.positionCm[column]+=h;
        minusEvent.primaryVertex.positionCm[column]-=h;
        PairSolution plus, minus;
        if (!SolvePair(firstTrack,firstPdg,firstCharge,firstMass,
                       secondTrack,secondPdg,secondCharge,secondMass,
                       plusEvent,config,plus) || !plus.converged ||
            !SolvePair(firstTrack,firstPdg,firstCharge,firstMass,
                       secondTrack,secondPdg,secondCharge,secondMass,
                       minusEvent,config,minus) || !minus.converged) {
          candidate.motherTopo.status=kProductionPvResponseFailed; return;
        }
        double yp[7], ym[7];
        MotherParameters(plus.fit.mother,yp); MotherParameters(minus.fit.mother,ym);
        for (int i=0; i<7; ++i) derivative[scale][i]=(yp[i]-ym[i])/(2.*h);
      }
      for (int i=0; i<7; ++i) {
        const double size=std::max(1.e-6,std::max(std::fabs(derivative[0][i]),std::fabs(derivative[1][i])));
        if (!std::isfinite(derivative[0][i]) || !std::isfinite(derivative[1][i]) ||
            std::fabs(derivative[0][i]-derivative[1][i])>tolerance*size) {
          candidate.motherTopo.status=kProductionPvResponseFailed; return;
        }
        response[i][column]=derivative[1][i];
      }
    }
    // Joint prior = [[C_cond + D V D^T, D V], [V D^T, V]]. The existing
    // seven-dimensional mother is left unchanged for legacy outputs/cuts.
    for (int i=0; i<7; ++i) for (int a=0; a<3; ++a)
      for (int b=0; b<3; ++b)
        input.motherVertexCovariance[i][a]+=response[i][b]*input.vertex.covarianceCm2.m[b][a];
    for (int i=0; i<7; ++i) for (int j=0; j<7; ++j)
      for (int a=0; a<3; ++a)
        input.mother.covariance.m[i][j]+=input.motherVertexCovariance[i][a]*response[j][a];
    Symmetrize(input.mother.covariance);
  }
  const bool ok=FitProductionVertex(input,config.productionVertexConfig,
    candidate.geometricChi2,candidate.geometricNdf,candidate.motherTopo,
    firstCharge+secondCharge,config.magneticFieldByTesla);
  candidate.motherTopo.correlationModel=model;
  if (ok) GetFlightObservables(candidate.motherTopo,candidate.topoFlight);
}

bool BuildCandidate(const CeeV0Lite::TrackInput& firstTrack,
                    int firstPdg,
                    int firstCharge,
                    double firstMass,
                    const CeeV0Lite::TrackInput& secondTrack,
                    int secondPdg,
                    int secondCharge,
                    double secondMass,
                    CeeV0Lite::CandidateType type,
                    const CeeV0Lite::ChannelCuts& cuts,
                    const CeeV0Lite::EventInput& event,
                    const CeeV0Lite::FinderConfig& config,
                    CeeV0Lite::Candidate& candidate)
{
  if (!ChargesAllow(firstTrack, firstCharge, secondTrack, secondCharge, config))
    return false;

  PairSolution solution;
  if (!SolvePair(firstTrack, firstPdg, firstCharge, firstMass,
                 secondTrack, secondPdg, secondCharge, secondMass,
                 event, config, solution)) return false;

  candidate = CeeV0Lite::Candidate();
  candidate.type = type;
  candidate.daughterId[0] = firstTrack.id;
  candidate.daughterId[1] = secondTrack.id;
  candidate.daughterPdg[0] = firstPdg;
  candidate.daughterPdg[1] = secondPdg;
  candidate.daughterAtVertex[0] = solution.fit.daughterAtVertex[0];
  candidate.daughterAtVertex[1] = solution.fit.daughterAtVertex[1];
  candidate.mother = solution.fit.mother;
  candidate.decayVertexCm = candidate.mother.positionCm;
  candidate.momentumGeV = candidate.mother.momentumGeV;
  candidate.daughterDcaCm = solution.fit.daughterDcaCm;
  candidate.geometricChi2 = solution.fit.chi2;
  candidate.geometricNdf = std::max(1, solution.fit.ndf);
  candidate.geometricChi2Ndf = candidate.geometricChi2 /
                               candidate.geometricNdf;
  candidate.decayTimeNs = solution.decayTimeNs;
  candidate.tofIterations = solution.iterations;
  candidate.tofConverged = solution.converged;

  const double mass2 = candidate.mother.energyGeV *
                       candidate.mother.energyGeV -
                       CeeV0Lite::Norm2(candidate.momentumGeV);
  if (!(mass2 > 0.) || !std::isfinite(mass2)) return false;
  candidate.massGeV = std::sqrt(mass2);

  const CeeV0Lite::Vec3 displacement = candidate.decayVertexCm -
                                       event.primaryVertex.positionCm;
  candidate.decayLengthCm = CeeV0Lite::Norm(displacement);
  if (!(candidate.decayLengthCm > 0.)) return false;
  const double parentMomentum = CeeV0Lite::Norm(candidate.momentumGeV);
  if (!(parentMomentum > 0.)) return false;
  candidate.pointingCosine = CeeV0Lite::Dot(displacement,
                                             candidate.momentumGeV) /
    (candidate.decayLengthCm * parentMomentum);

  // Cheap central-value cuts are evaluated before the two extra shared-T0
  // solves.  None of these depends on the final mother covariance.
  if (!cuts.enabled || candidate.massGeV < cuts.minimumMassGeV ||
      candidate.massGeV > cuts.maximumMassGeV ||
      candidate.daughterDcaCm > cuts.maximumDaughterDcaCm ||
      candidate.geometricChi2Ndf > cuts.maximumGeometricChi2Ndf ||
      candidate.decayLengthCm < cuts.minimumDecayLengthCm ||
      candidate.pointingCosine < cuts.minimumPointingCosine) return false;

  const DaughterPvMetrics firstPv = MeasureDaughterAtPv(
    solution.inputStates[0], event, config);
  const DaughterPvMetrics secondPv = MeasureDaughterAtPv(
    solution.inputStates[1], event, config);
  if (!firstPv.valid || !secondPv.valid) return false;
  candidate.daughterDcaToPvCm[0] = firstPv.dcaCm;
  candidate.daughterDcaToPvCm[1] = secondPv.dcaCm;
  candidate.daughterChi2FromPv[0] = firstPv.chi2;
  candidate.daughterChi2FromPv[1] = secondPv.chi2;
  candidate.daughterPvSignificance[0] = firstPv.significance;
  candidate.daughterPvSignificance[1] = secondPv.significance;
  candidate.daughterDeviationFromPv[0] = firstPv.significance;
  candidate.daughterDeviationFromPv[1] = secondPv.significance;
  if (candidate.daughterChi2FromPv[0] <
        cuts.minimumDaughterChi2FromPv ||
      candidate.daughterChi2FromPv[1] <
        cuts.minimumDaughterChi2FromPv ||
      candidate.daughterPvSignificance[0] <
        cuts.minimumDaughterDeviationFromPv ||
      candidate.daughterPvSignificance[1] <
        cuts.minimumDaughterDeviationFromPv) return false;

  // Mirror the existing CEE zero-field Lambda role assignment. Daughter 0 is
  // the (anti)proton hypothesis and daughter 1 is the pion hypothesis.
  if ((type == CeeV0Lite::kLambda || type == CeeV0Lite::kAntiLambda) &&
      config.useLambdaDcaRoleAssignment) {
    if (!(candidate.daughterDcaToPvCm[0] <=
          firstTrack.lambdaRoleDcaBoundaryCm)) return false;
    if (!(candidate.daughterDcaToPvCm[1] >
          secondTrack.lambdaRoleDcaBoundaryCm)) return false;
  }

  // Re-run the complete POCA/vertex estimator for each input-state
  // perturbation.  This outer Jacobian contains ds_closest/d(track state),
  // the key path correlation that a fixed-path 6x6 transport misses.
  CeeV0Lite::VertexFitter covarianceFitter(config.vertexFitterConfig,
                                            config.transportConfig);
  const bool covarianceComplete =
    covarianceFitter.RecalculateMotherCovariance(
      solution.inputStates[0], solution.inputStates[1],
      config.magneticFieldByTesla, candidate.mother);
  candidate.motherJacobianComplete = covarianceComplete;

  if (!covarianceComplete) {
    // Only this exact outer-Jacobian failure may be retained, and only when
    // explicitly requested.  The nominal POCA, common vertex, TOF fixed point
    // and central mother state have already succeeded above.
    if (candidate.mother.covarianceStatus !=
          CeeV0Lite::kMotherCovarianceJacobianPerturbationFailed ||
        cuts.requireCompleteMotherJacobian) return false;
    InvalidateCovarianceDependentObservables(candidate);
  }
  else {
    if (!AddSharedTimeUncertainty(
          firstTrack, firstPdg, firstCharge, firstMass,
          secondTrack, secondPdg, secondCharge, secondMass,
          event, config, solution, candidate.mother)) return false;
    candidate.decayVertexCm = candidate.mother.positionCm;
    candidate.momentumGeV = candidate.mother.momentumGeV;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        candidate.decayVertexCovarianceCm2.m[i][j] =
          candidate.mother.covariance.m[i][j];
    if (!CeeV0Lite::MotherMassAndError(candidate.mother,
                                       candidate.massGeV,
                                       candidate.sigmaMassGeV)) return false;

    const CeeV0Lite::Vec3 lengthDirection = displacement /
                                           candidate.decayLengthCm;
    const CeeV0Lite::Mat3 lengthCovariance =
      candidate.decayVertexCovarianceCm2 +
      event.primaryVertex.covarianceCm2;
    const double lengthVariance = CeeV0Lite::QuadraticForm(
      lengthDirection, lengthCovariance);
    candidate.sigmaDecayLengthCm = std::sqrt(std::max(0., lengthVariance));
    if (!(candidate.sigmaDecayLengthCm > 0.) ||
        !std::isfinite(candidate.sigmaDecayLengthCm)) return false;
    candidate.decayLengthSignificance = candidate.decayLengthCm /
                                        candidate.sigmaDecayLengthCm;
    candidate.primaryVertexConstraintChi2 = PrimaryVertexConstraintChi2(
      candidate, event,
      config.vertexFitterConfig.measurementVarianceFloorCm2);
    candidate.primaryVertexConstraintNdf = 2;
    candidate.primaryVertexConstraintChi2Ndf =
      candidate.primaryVertexConstraintChi2 /
      candidate.primaryVertexConstraintNdf;
    // A two-track decay vertex has one geometric degree of freedom.  Adding a
    // production vertex removes the unobservable mother path and contributes
    // two effective constraints, matching KFParticle::SetProductionVertex's
    // NDF accounting. Keep the pure PV residual above for diagnostics.
    candidate.topologicalChi2 = candidate.geometricChi2 +
                                candidate.primaryVertexConstraintChi2;
    candidate.topologicalNdf = candidate.geometricNdf +
                               candidate.primaryVertexConstraintNdf;
    candidate.topologicalChi2Ndf = candidate.topologicalChi2 /
                                   candidate.topologicalNdf;
    candidate.covarianceDependentCutsEvaluated = true;
  }

  if (!std::isfinite(candidate.massGeV) ||
      !std::isfinite(candidate.geometricChi2Ndf)) return false;
  if (candidate.covarianceDependentCutsEvaluated &&
      (!std::isfinite(candidate.sigmaMassGeV) ||
       !std::isfinite(candidate.topologicalChi2Ndf))) return false;
  if (!PassCuts(candidate, cuts)) return false;
  BuildProductionCopy(firstTrack,firstPdg,firstCharge,firstMass,
    secondTrack,secondPdg,secondCharge,secondMass,event,config,solution,candidate);
  return true;
}
}

namespace CeeV0Lite
{

ChannelCuts::ChannelCuts()
  : enabled(true),
    minimumMassGeV(0.),
    maximumMassGeV(10.),
    maximumDaughterDcaCm(3.),
    maximumGeometricChi2Ndf(25.),
    maximumTopologicalChi2Ndf(25.),
    minimumDecayLengthCm(0.2),
    minimumDecayLengthSignificance(0.),
    minimumPointingCosine(0.90),
    minimumDaughterChi2FromPv(0.),
    minimumDaughterDeviationFromPv(0.),
    requireCompleteMotherJacobian(true)
{
}

FinderConfig::FinderConfig()
  : momentumMode(kAutoMomentum),
    magneticFieldByTesla(0.),
    reconstructK0Short(false),
    reconstructLambda(true),
    reconstructAntiLambda(false),
    useLambdaDcaRoleAssignment(true),
    requireMeasuredChargeInMagneticField(true),
    allowTofMomentumInMagneticField(false),
    zeroFieldThresholdTesla(1.e-6),
    usePrimaryVertexDcaReference(true),
    maximumTofIterations(10),
    tofVertexConvergenceCm(1.e-3),
    tofTimeConvergenceNs(1.e-4),
    maximumBeta(0.9995),
    additionalDecayTimeSigmaNs(0.),
    propagateSharedTimeUncertainty(false),
    sharedTimeDerivativeStepNs(1.e-3),
    sharedTimeDerivativeStabilityTolerance(0.05),
    seedParentBeta(0.7),
    requireTofConvergence(true),
    defaultRelativeTrackingMomentumError(0.05),
    maximumCandidates(100000),
    buildProductionState(true)
{
  referenceTransportConfig.maximumClosestApproachIterations = 64;
  referenceTransportConfig.closestApproachToleranceCm = 1.e-8;
  referenceTransportConfig.maximumAbsPathCm = 1000.;
  k0Cuts.minimumMassGeV = 0.35;
  k0Cuts.maximumMassGeV = 0.65;
  lambdaCuts.minimumMassGeV = 1.05;
  lambdaCuts.maximumMassGeV = 1.25;
}

TrackReferenceDiagnostic::TrackReferenceDiagnostic()
  : trackId(-1), status(kReferenceUnchanged),
    dcaCm(std::numeric_limits<double>::quiet_NaN())
{
}

bool PrepareTrackReference(const TrackInput& input, const EventInput& event,
                           const FinderConfig& config, TrackInput& output,
                           TrackReferenceDiagnostic* diagnostic)
{
  TrackReferenceDiagnostic info;
  info.trackId = input.id;
  const bool enabled = config.usePrimaryVertexDcaReference &&
    std::fabs(config.magneticFieldByTesla) > config.zeroFieldThresholdTesla;
  if (!enabled) {
    output = input;
    if (diagnostic) *diagnostic = info;
    return true;
  }
  // Never guess a magnetic trajectory or construct a new covariance here.
  if (!UseTrackingMomentum(input, config) || !input.hasTrackingCovariance ||
      !input.chargeKnown || input.charge == 0 ||
      !IsFinite(input.positionCm) || !IsFinite(input.trackingMomentumGeV) ||
      !(Norm(input.trackingMomentumGeV) > 1.e-10) ||
      !IsPositiveSemidefinite(input.trackingCovariance) ||
      !IsFinite(event.primaryVertex.positionCm)) {
    info.status = kReferenceInvalidInput;
    if (diagnostic) *diagnostic = info;
    return false;
  }
  TrackState state, anchored;
  state.id = input.id;
  state.charge = input.charge;
  state.positionCm = input.positionCm;
  state.momentumGeV = input.trackingMomentumGeV;
  state.covariance = input.trackingCovariance;
  if (!ReanchorAtClosestApproach(state, event.primaryVertex.positionCm,
        config.magneticFieldByTesla, config.referenceTransportConfig,
        anchored, &info.closestApproach)) {
    info.status = kReferenceTransportFailed;
    if (diagnostic) *diagnostic = info;
    return false;
  }
  output = input;
  output.positionCm = anchored.positionCm;
  output.trackingMomentumGeV = anchored.momentumGeV;
  Normalize(anchored.momentumGeV, output.direction);
  output.trackingCovariance = anchored.covariance;
  output.positionCovarianceCm2 = PositionBlock(anchored.covariance);
  info.status = kReferenceAtPvDca;
  info.dcaCm = Norm(anchored.positionCm - event.primaryVertex.positionCm);
  if (diagnostic) *diagnostic = info;
  return true;
}

Finder::Finder(const FinderConfig& config)
  : fConfig(config)
{
}

void Finder::Reconstruct(const std::vector<TrackInput>& tracks,
                         const EventInput& event,
                         std::vector<Candidate>& candidates,
                         std::vector<TrackReferenceDiagnostic>* diagnostics) const
{
  candidates.clear();
  if (diagnostics) diagnostics->clear();
  // Cache once before pairing (not per mass hypothesis or numerical refit).
  std::vector<TrackInput> prepared;
  prepared.reserve(tracks.size());
  for (std::size_t i = 0; i < tracks.size(); ++i) {
    TrackInput track;
    TrackReferenceDiagnostic info;
    if (PrepareTrackReference(tracks[i], event, fConfig, track, &info))
      prepared.push_back(track);
    if (diagnostics) diagnostics->push_back(info);
  }
  for (std::size_t i = 0; i < prepared.size(); ++i) {
    for (std::size_t j = i + 1; j < prepared.size(); ++j) {
      if (candidates.size() >= fConfig.maximumCandidates) return;
      const TrackInput& first = prepared[i];
      const TrackInput& second = prepared[j];
      Candidate candidate;

      if (fConfig.reconstructK0Short &&
          HasHypothesis(first, kPionHypothesis) &&
          HasHypothesis(second, kPionHypothesis)) {
        const TrackInput* positive = &first;
        const TrackInput* negative = &second;
        if (first.chargeKnown && second.chargeKnown &&
            first.charge < 0 && second.charge > 0) {
          positive = &second;
          negative = &first;
        }
        if (BuildCandidate(*positive, 211, +1, kPionMassGeV,
                           *negative, -211, -1, kPionMassGeV,
                           kK0Short, fConfig.k0Cuts, event, fConfig,
                           candidate)) {
          candidates.push_back(candidate);
          if (candidates.size() >= fConfig.maximumCandidates) return;
        }
      }

      if (fConfig.reconstructLambda) {
        if (HasHypothesis(first, kProtonHypothesis) &&
            HasHypothesis(second, kPionHypothesis) &&
            BuildCandidate(first, 2212, +1, kProtonMassGeV,
                           second, -211, -1, kPionMassGeV,
                           kLambda, fConfig.lambdaCuts, event, fConfig,
                           candidate)) {
          candidates.push_back(candidate);
          if (candidates.size() >= fConfig.maximumCandidates) return;
        }
        if (HasHypothesis(second, kProtonHypothesis) &&
            HasHypothesis(first, kPionHypothesis) &&
            BuildCandidate(second, 2212, +1, kProtonMassGeV,
                           first, -211, -1, kPionMassGeV,
                           kLambda, fConfig.lambdaCuts, event, fConfig,
                           candidate)) {
          candidates.push_back(candidate);
          if (candidates.size() >= fConfig.maximumCandidates) return;
        }
      }

      if (fConfig.reconstructAntiLambda) {
        if (HasHypothesis(first, kProtonHypothesis) &&
            HasHypothesis(second, kPionHypothesis) &&
            BuildCandidate(first, -2212, -1, kProtonMassGeV,
                           second, 211, +1, kPionMassGeV,
                           kAntiLambda, fConfig.lambdaCuts, event, fConfig,
                           candidate)) {
          candidates.push_back(candidate);
          if (candidates.size() >= fConfig.maximumCandidates) return;
        }
        if (HasHypothesis(second, kProtonHypothesis) &&
            HasHypothesis(first, kPionHypothesis) &&
            BuildCandidate(second, -2212, -1, kProtonMassGeV,
                           first, 211, +1, kPionMassGeV,
                           kAntiLambda, fConfig.lambdaCuts, event, fConfig,
                           candidate)) {
          candidates.push_back(candidate);
          if (candidates.size() >= fConfig.maximumCandidates) return;
        }
      }
    }
  }
}

} // namespace CeeV0Lite
