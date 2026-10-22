#include "CeeV0LiteFinder.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace
{
const double kC = 29.9792458;
const double kLambdaMass = 1.115683;
const double kK0ShortMass = 0.497611;
const double kProtonMass = 0.9382720813;
const double kPionMass = 0.13957039;

double TwoBodyMomentum(double mother, double first, double second)
{
  const double a = mother * mother - (first + second) * (first + second);
  const double b = mother * mother - (first - second) * (first - second);
  return std::sqrt(a * b) / (2. * mother);
}

struct SyntheticDecay
{
  CeeV0Lite::Vec3 vertex;
  CeeV0Lite::Vec3 protonMomentum;
  CeeV0Lite::Vec3 pionMomentum;
  double parentBeta;
};

SyntheticDecay MakeBoostedLambda()
{
  SyntheticDecay decay;
  decay.vertex = CeeV0Lite::Vec3(0., 0., 5.);
  decay.parentBeta = 0.6;
  const double gamma = 1. / std::sqrt(1. - decay.parentBeta * decay.parentBeta);
  const double pStar = TwoBodyMomentum(kLambdaMass, kProtonMass, kPionMass);
  const double protonEnergyStar = std::sqrt(pStar * pStar +
                                             kProtonMass * kProtonMass);
  const double pionEnergyStar = std::sqrt(pStar * pStar +
                                           kPionMass * kPionMass);
  decay.protonMomentum = CeeV0Lite::Vec3(
    pStar, 0., gamma * decay.parentBeta * protonEnergyStar);
  decay.pionMomentum = CeeV0Lite::Vec3(
    -pStar, 0., gamma * decay.parentBeta * pionEnergyStar);
  return decay;
}

SyntheticDecay MakeBoostedK0Short()
{
  SyntheticDecay decay;
  decay.vertex = CeeV0Lite::Vec3(0., 0., 4.);
  decay.parentBeta = 0.7;
  const double gamma = 1. / std::sqrt(1. - decay.parentBeta * decay.parentBeta);
  const double pStar = TwoBodyMomentum(kK0ShortMass, kPionMass, kPionMass);
  const double pionEnergyStar = std::sqrt(pStar * pStar +
                                           kPionMass * kPionMass);
  decay.protonMomentum = CeeV0Lite::Vec3(
    pStar, 0., gamma * decay.parentBeta * pionEnergyStar);
  decay.pionMomentum = CeeV0Lite::Vec3(
    -pStar, 0., gamma * decay.parentBeta * pionEnergyStar);
  return decay;
}

CeeV0Lite::TrackInput MakeDaughter(const CeeV0Lite::Vec3& vertex,
                                   const CeeV0Lite::Vec3& momentum,
                                   double mass,
                                   int id,
                                   int charge,
                                   unsigned int mask,
                                   bool trackingMode,
                                   double parentFlightTimeNs)
{
  CeeV0Lite::Vec3 direction;
  CeeV0Lite::Normalize(momentum, direction);
  CeeV0Lite::TrackInput track;
  track.id = id;
  track.positionCm = vertex + 20. * direction;
  track.direction = direction;
  track.positionCovarianceCm2 = CeeV0Lite::Mat3::Identity(0.02 * 0.02);
  track.sigmaDirectionRad = 0.0005;
  track.hypothesisMask = mask;
  track.chargeKnown = trackingMode;
  track.charge = charge;
  track.hasTrackingMomentum = trackingMode;
  track.trackingMomentumGeV = momentum;
  if (trackingMode) {
    track.hasTrackingCovariance = true;
    for (int i = 0; i < 3; ++i) {
      track.trackingCovariance.m[i][i] = 0.02 * 0.02;
      track.trackingCovariance.m[i + 3][i + 3] = 0.001 * 0.001;
    }
  }

  track.hasTof = true;
  track.tofPositionCm = vertex + 100. * direction;
  const double energy = std::sqrt(CeeV0Lite::Norm2(momentum) + mass * mass);
  const double beta = CeeV0Lite::Norm(momentum) / energy;
  track.tofTimeNs = parentFlightTimeNs + 100. / (beta * kC);
  track.sigmaTofTimeNs = 0.02;
  track.sigmaPathLengthCm = 0.05;
  return track;
}

CeeV0Lite::TrackInput MakeCurvedDaughter(
    const CeeV0Lite::Vec3& vertex,
    const CeeV0Lite::Vec3& momentumAtVertex,
    double mass,
    int id,
    int charge,
    unsigned int mask,
    double magneticFieldByTesla,
    double referencePathCm)
{
  CeeV0Lite::TrackInput track = MakeDaughter(
    vertex, momentumAtVertex, mass, id, charge, mask, true, 0.);
  CeeV0Lite::TrackState atVertex;
  atVertex.id = id;
  atVertex.charge = charge;
  atVertex.massGeV = mass;
  atVertex.positionCm = vertex;
  atVertex.momentumGeV = momentumAtVertex;
  atVertex.covariance = track.trackingCovariance;
  CeeV0Lite::TrackState atReference;
  CeeV0Lite::TransportConfig transport;
  assert(CeeV0Lite::PropagateState(atVertex, referencePathCm,
                                   magneticFieldByTesla, transport,
                                   atReference));
  track.positionCm = atReference.positionCm;
  assert(CeeV0Lite::Normalize(atReference.momentumGeV, track.direction));
  track.trackingMomentumGeV = atReference.momentumGeV;
  track.trackingCovariance = atReference.covariance;
  // The tracking-momentum path must not use a straight TOF direction to orient
  // a curved reference state.
  track.hasTof = false;
  return track;
}

void CheckProductionCopy(const CeeV0Lite::Candidate& c)
{
  if (!c.motherTopo.valid) std::cerr << "production status " << c.motherTopo.status << "\n";
  assert(c.motherTopo.valid && c.motherTopo.atProductionVertex);
  assert(c.topoFlight.valid && c.topoFlight.sigmaCtauCm>0.);
  assert(CeeV0Lite::IsPositiveSemidefinite(c.motherTopo.covariance));
  assert(c.motherTopo.totalNdf==c.geometricNdf+2);
  assert(c.motherTopo.flightParameterS>0.);
  CeeV0Lite::ProductionState atDecay=c.motherTopo;
  assert(CeeV0Lite::TransportToDecayVertex(atDecay));
  assert(CeeV0Lite::Norm(atDecay.positionCm-c.decayVertexCm)<3.e-3);
  assert(std::fabs(c.topoFlight.ctauCm-c.motherTopo.flightParameterS*c.topoFlight.massGeV)<1.e-12);
}

void CheckCandidate(const CeeV0Lite::Candidate& candidate,
                    const CeeV0Lite::Vec3& trueVertex,
                    double massTolerance,
                    double vertexTolerance)
{
  CheckProductionCopy(candidate);
  assert(candidate.type == CeeV0Lite::kLambda);
  assert(candidate.mother.covarianceValid);
  assert(candidate.mother.covarianceStatus ==
         CeeV0Lite::kMotherCovarianceValid);
  assert(candidate.motherJacobianComplete);
  assert(candidate.covarianceDependentCutsEvaluated);
  assert(candidate.mother.jacobianFitsAttempted == 24);
  assert(candidate.mother.jacobianFitsSucceeded == 24);
  assert(candidate.mother.jacobianFailureCount == 0);
  assert(CeeV0Lite::IsFinite(candidate.mother.covariance));
  assert(candidate.mother.energyGeV > CeeV0Lite::Norm(candidate.momentumGeV));
  assert(std::isfinite(candidate.sigmaMassGeV));
  assert(candidate.sigmaMassGeV > 0.);
  assert(std::fabs(candidate.massGeV - kLambdaMass) < massTolerance);
  assert(CeeV0Lite::Norm(candidate.decayVertexCm - trueVertex) < vertexTolerance);
  assert(candidate.pointingCosine > 0.999);
  assert(candidate.daughterDcaToPvCm[0] <= 1.2);
  assert(candidate.daughterDcaToPvCm[1] > 1.2);
  assert(candidate.geometricNdf == 1);
  assert(std::fabs(candidate.geometricChi2Ndf -
                   candidate.geometricChi2 / candidate.geometricNdf) < 1.e-12);
  assert(candidate.primaryVertexConstraintNdf == 2);
  assert(std::fabs(candidate.primaryVertexConstraintChi2Ndf -
                   candidate.primaryVertexConstraintChi2 /
                     candidate.primaryVertexConstraintNdf) < 1.e-12);
  assert(candidate.topologicalNdf == candidate.geometricNdf +
                                      candidate.primaryVertexConstraintNdf);
  assert(std::fabs(candidate.topologicalChi2 -
                   (candidate.geometricChi2 +
                    candidate.primaryVertexConstraintChi2)) < 1.e-12);
  assert(std::fabs(candidate.topologicalChi2Ndf -
                   candidate.topologicalChi2 / candidate.topologicalNdf) <
         1.e-12);
  for (int i = 0; i < 2; ++i) {
    assert(std::fabs(candidate.daughterChi2FromPv[i] -
                     candidate.daughterPvSignificance[i] *
                       candidate.daughterPvSignificance[i]) < 1.e-10);
    assert(candidate.daughterDeviationFromPv[i] ==
           candidate.daughterPvSignificance[i]);
  }
}

void CheckK0ShortCandidate(const CeeV0Lite::Candidate& candidate,
                           const CeeV0Lite::Vec3& trueVertex,
                           double massTolerance,
                           double vertexTolerance)
{
  CheckProductionCopy(candidate);
  assert(candidate.type == CeeV0Lite::kK0Short);
  assert(candidate.mother.covarianceValid);
  assert(std::isfinite(candidate.sigmaMassGeV));
  assert(candidate.sigmaMassGeV > 0.);
  assert(std::fabs(candidate.massGeV - kK0ShortMass) < massTolerance);
  assert(CeeV0Lite::Norm(candidate.decayVertexCm - trueVertex) < vertexTolerance);
  assert(candidate.pointingCosine > 0.999);
}

void MotherParameters(const CeeV0Lite::MotherState& mother, double values[7])
{
  for (int i = 0; i < 3; ++i) {
    values[i] = mother.positionCm[i];
    values[3 + i] = mother.momentumGeV[i];
  }
  values[6] = mother.energyGeV;
}
}

int main()
{
  // Library defaults are the ideal simulation policy: generator time origin
  // t=0, no measured/shared T0 uncertainty.  The independent stop-detector
  // timing resolution remains a TrackInput property.
  const CeeV0Lite::EventInput defaultEvent;
  const CeeV0Lite::FinderConfig defaultFinderConfig;
  assert(defaultEvent.startTimeNs == 0.);
  assert(defaultEvent.sigmaStartTimeNs == 0.);
  assert(defaultFinderConfig.additionalDecayTimeSigmaNs == 0.);
  assert(!defaultFinderConfig.propagateSharedTimeUncertainty);

  const SyntheticDecay decay = MakeBoostedLambda();
  const double parentFlightTime = CeeV0Lite::Norm(decay.vertex) /
                                  (decay.parentBeta * kC);

  CeeV0Lite::EventInput event;
  event.primaryVertex.positionCm = CeeV0Lite::Vec3(0., 0., 0.);
  event.primaryVertex.covarianceCm2 = CeeV0Lite::Mat3::Identity(0.01 * 0.01);
  event.startTimeNs = 0.;
  event.sigmaStartTimeNs = 0.01;

  // Future/curvature-momentum path.
  std::vector<CeeV0Lite::TrackInput> trackingTracks;
  trackingTracks.push_back(MakeDaughter(
    decay.vertex, decay.protonMomentum, kProtonMass, 0, +1,
    CeeV0Lite::kProtonHypothesis, true, parentFlightTime));
  trackingTracks.push_back(MakeDaughter(
    decay.vertex, decay.pionMomentum, kPionMass, 1, -1,
    CeeV0Lite::kPionHypothesis, true, parentFlightTime));
  CeeV0Lite::FinderConfig trackingConfig;
  trackingConfig.momentumMode = CeeV0Lite::kTrackingMomentum;
  trackingConfig.magneticFieldByTesla = 0.;
  CeeV0Lite::Finder trackingFinder(trackingConfig);
  std::vector<CeeV0Lite::Candidate> candidates;
  trackingFinder.Reconstruct(trackingTracks, event, candidates);
  assert(candidates.size() == 1u);
  CheckCandidate(candidates[0], decay.vertex, 2.e-4, 1.e-4);
  assert(std::fabs(candidates[0].decayTimeNs - parentFlightTime) < 1.e-3);
  const CeeV0Lite::Candidate trackingBaseline = candidates[0];

  // The strict default rejects a candidate whose nominal fit succeeded but
  // whose outer numerical-Jacobian refits are incomplete.  The explicit QA
  // mode retains the same central candidate, marks the failure, and writes
  // every covariance-dependent observable as NaN so it cannot be mistaken
  // for an exceptionally precise result.
  const CeeV0Lite::Vec3 qaVertex(5., 5., 0.);
  std::vector<CeeV0Lite::TrackInput> jacobianFailureTracks;
  jacobianFailureTracks.push_back(MakeDaughter(
    qaVertex, CeeV0Lite::Vec3(0.116, 0., 0.), kProtonMass, 20, +1,
    CeeV0Lite::kProtonHypothesis, true, 0.));
  jacobianFailureTracks.push_back(MakeDaughter(
    qaVertex, CeeV0Lite::Vec3(0., 0.116, 0.), kPionMass, 21, -1,
    CeeV0Lite::kPionHypothesis, true, 0.));
  CeeV0Lite::FinderConfig jacobianFailureConfig = trackingConfig;
  jacobianFailureConfig.useLambdaDcaRoleAssignment = false;
  jacobianFailureConfig.lambdaCuts.minimumDecayLengthCm = 0.;
  jacobianFailureConfig.lambdaCuts.minimumPointingCosine = -1.;
  jacobianFailureConfig.vertexFitterConfig.
    motherJacobianMinimumMomentumStepGeV = 0.116;
  assert(jacobianFailureConfig.lambdaCuts.requireCompleteMotherJacobian);
  CeeV0Lite::Finder jacobianFailureFinder(jacobianFailureConfig);
  jacobianFailureFinder.Reconstruct(
    jacobianFailureTracks, event, candidates);
  assert(candidates.empty());
  jacobianFailureConfig.lambdaCuts.requireCompleteMotherJacobian = false;
  jacobianFailureFinder.SetConfig(jacobianFailureConfig);
  jacobianFailureFinder.Reconstruct(
    jacobianFailureTracks, event, candidates);
  assert(candidates.size() == 1u);
  const CeeV0Lite::Candidate& retained = candidates[0];
  assert(!retained.motherJacobianComplete);
  assert(!retained.covarianceDependentCutsEvaluated);
  assert(!retained.mother.covarianceValid);
  assert(retained.mother.covarianceStatus ==
         CeeV0Lite::kMotherCovarianceJacobianPerturbationFailed);
  assert(retained.mother.jacobianFitsAttempted == 24);
  assert(retained.mother.jacobianFailureCount == 2);
  assert(std::isfinite(retained.massGeV));
  assert(std::isfinite(retained.decayVertexCm[0]));
  assert(std::isfinite(retained.daughterDcaCm));
  assert(std::isfinite(retained.geometricChi2Ndf));
  assert(std::isfinite(retained.decayLengthCm));
  assert(std::isfinite(retained.pointingCosine));
  assert(std::isnan(retained.sigmaMassGeV));
  assert(std::isnan(retained.decayVertexCovarianceCm2.m[0][0]));
  assert(std::isnan(retained.sigmaDecayLengthCm));
  assert(std::isnan(retained.decayLengthSignificance));
  assert(std::isnan(retained.primaryVertexConstraintChi2));
  assert(std::isnan(retained.primaryVertexConstraintChi2Ndf));
  assert(std::isnan(retained.topologicalChi2));
  assert(std::isnan(retained.topologicalChi2Ndf));
  assert(!retained.motherTopo.valid && retained.motherTopo.status==CeeV0Lite::kProductionInvalidMother);
  assert(std::isnan(retained.motherTopo.flightParameterS));

  // The opt-in is deliberately narrow: it must not retain other covariance
  // failures such as unsupported stochastic process noise.
  jacobianFailureConfig.transportConfig.positionProcessNoiseCm2PerCm = 1.e-6;
  jacobianFailureFinder.SetConfig(jacobianFailureConfig);
  jacobianFailureFinder.Reconstruct(
    jacobianFailureTracks, event, candidates);
  assert(candidates.empty());

  // The canonical daughter-PV cut now uses raw chi2, while the deprecated
  // compatibility cut retains the old sqrt(chi2) convention.
  const double smallerDaughterPvChi2 = std::min(
    trackingBaseline.daughterChi2FromPv[0],
    trackingBaseline.daughterChi2FromPv[1]);
  CeeV0Lite::FinderConfig metricCutConfig = trackingConfig;
  metricCutConfig.lambdaCuts.minimumDaughterChi2FromPv =
    smallerDaughterPvChi2 * (1. - 1.e-6);
  CeeV0Lite::Finder metricCutFinder(metricCutConfig);
  metricCutFinder.Reconstruct(trackingTracks, event, candidates);
  assert(candidates.size() == 1u);
  metricCutConfig.lambdaCuts.minimumDaughterChi2FromPv =
    smallerDaughterPvChi2 * (1. + 1.e-6);
  metricCutFinder.SetConfig(metricCutConfig);
  metricCutFinder.Reconstruct(trackingTracks, event, candidates);
  assert(candidates.empty());

  metricCutConfig = trackingConfig;
  metricCutConfig.lambdaCuts.minimumDaughterDeviationFromPv =
    std::sqrt(smallerDaughterPvChi2) * (1. + 1.e-6);
  metricCutFinder.SetConfig(metricCutConfig);
  metricCutFinder.Reconstruct(trackingTracks, event, candidates);
  assert(candidates.empty());

  CeeV0Lite::EventInput topologyEvent = event;
  topologyEvent.primaryVertex.positionCm[0] += 0.02;
  metricCutConfig = trackingConfig;
  metricCutConfig.lambdaCuts.maximumTopologicalChi2Ndf = 1.e30;
  metricCutFinder.SetConfig(metricCutConfig);
  metricCutFinder.Reconstruct(trackingTracks, topologyEvent, candidates);
  assert(candidates.size() == 1u);
  const double nonzeroTopologicalChi2Ndf = candidates[0].topologicalChi2Ndf;
  assert(nonzeroTopologicalChi2Ndf > 0.);
  metricCutConfig.lambdaCuts.maximumTopologicalChi2Ndf =
    nonzeroTopologicalChi2Ndf * (1. + 1.e-6);
  metricCutFinder.SetConfig(metricCutConfig);
  metricCutFinder.Reconstruct(trackingTracks, topologyEvent, candidates);
  assert(candidates.size() == 1u);
  metricCutConfig.lambdaCuts.maximumTopologicalChi2Ndf =
    nonzeroTopologicalChi2Ndf * (1. - 1.e-6);
  metricCutFinder.SetConfig(metricCutConfig);
  metricCutFinder.Reconstruct(trackingTracks, topologyEvent, candidates);
  assert(candidates.empty());

  // Current zero-field TOF-momentum path.
  std::vector<CeeV0Lite::TrackInput> tofTracks;
  tofTracks.push_back(MakeDaughter(
    decay.vertex, decay.protonMomentum, kProtonMass, 0, +1,
    CeeV0Lite::kProtonHypothesis | CeeV0Lite::kPionHypothesis,
    false, parentFlightTime));
  tofTracks.push_back(MakeDaughter(
    decay.vertex, decay.pionMomentum, kPionMass, 1, -1,
    CeeV0Lite::kProtonHypothesis | CeeV0Lite::kPionHypothesis,
    false, parentFlightTime));
  CeeV0Lite::FinderConfig tofConfig;
  tofConfig.momentumMode = CeeV0Lite::kTofMomentum;
  tofConfig.magneticFieldByTesla = 0.;
  tofConfig.requireTofConvergence = true;
  // Enable the real-data-style common-time term explicitly for the dedicated
  // shared-T0 regression below.  Simulation defaults leave this disabled.
  tofConfig.propagateSharedTimeUncertainty = true;
  CeeV0Lite::Finder tofFinder(tofConfig);
  tofFinder.Reconstruct(tofTracks, event, candidates);
  assert(candidates.size() == 1u);
  CheckCandidate(candidates[0], decay.vertex, 2.e-3, 1.e-3);
  assert(candidates[0].tofConverged);
  assert(std::fabs(candidates[0].decayTimeNs - parentFlightTime) < 1.e-3);

  const CeeV0Lite::Candidate topoEnabled=candidates[0];
  assert(topoEnabled.motherTopo.correlationModel==CeeV0Lite::kProductionTofPvResponse);
  CeeV0Lite::FinderConfig withoutTopoConfig=tofConfig;
  withoutTopoConfig.buildProductionState=false;
  CeeV0Lite::Finder withoutTopo(withoutTopoConfig);
  withoutTopo.Reconstruct(tofTracks,event,candidates);
  assert(candidates.size()==1u);
  assert(candidates[0].motherTopo.status==CeeV0Lite::kProductionDisabled);
  assert(!candidates[0].motherTopo.valid && std::isnan(candidates[0].topoFlight.ctauCm));
  assert(candidates[0].massGeV==topoEnabled.massGeV);
  assert(candidates[0].topologicalChi2==topoEnabled.topologicalChi2);
  for(int i=0;i<7;++i)for(int j=0;j<7;++j)
    assert(candidates[0].mother.covariance.m[i][j]==topoEnabled.mother.covariance.m[i][j]);
  // A failed auxiliary derivative must be explicit and preserve the old row.
  withoutTopoConfig.buildProductionState=true;
  withoutTopoConfig.productionVertexConfig.pvDerivativeStepCm=0.;
  withoutTopo.SetConfig(withoutTopoConfig);
  withoutTopo.Reconstruct(tofTracks,event,candidates);
  assert(candidates.size()==1u && candidates[0].massGeV==topoEnabled.massGeV);
  assert(candidates[0].motherTopo.status==CeeV0Lite::kProductionPvResponseFailed);
  assert(!candidates[0].motherTopo.valid && std::isnan(candidates[0].motherTopo.flightParameterS));

  // A shared T0 is propagated once as a common nuisance.  Doubling sigma(T0)
  // must multiply its covariance contribution by four while leaving the
  // central mass unchanged.
  CeeV0Lite::FinderConfig sharedTimeConfig = tofConfig;
  sharedTimeConfig.additionalDecayTimeSigmaNs = 0.;
  CeeV0Lite::Finder sharedTimeFinder(sharedTimeConfig);
  CeeV0Lite::EventInput zeroT0Event = event;
  zeroT0Event.sigmaStartTimeNs = 0.;
  sharedTimeFinder.Reconstruct(tofTracks, zeroT0Event, candidates);
  assert(candidates.size() == 1u);
  const CeeV0Lite::Candidate zeroT0 = candidates[0];
  CeeV0Lite::FinderConfig disabledSharedTimeConfig = sharedTimeConfig;
  disabledSharedTimeConfig.propagateSharedTimeUncertainty = false;
  CeeV0Lite::Finder disabledSharedTimeFinder(disabledSharedTimeConfig);
  disabledSharedTimeFinder.Reconstruct(tofTracks, zeroT0Event, candidates);
  assert(candidates.size() == 1u);
  const CeeV0Lite::Candidate disabledSharedTime = candidates[0];
  assert(CeeV0Lite::Norm(zeroT0.mother.positionCm -
                         disabledSharedTime.mother.positionCm) < 1.e-14);
  assert(CeeV0Lite::Norm(zeroT0.mother.momentumGeV -
                         disabledSharedTime.mother.momentumGeV) < 1.e-14);
  assert(std::fabs(zeroT0.mother.energyGeV -
                   disabledSharedTime.mother.energyGeV) < 1.e-14);
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      assert(std::fabs(zeroT0.mother.covariance.m[i][j] -
                       disabledSharedTime.mother.covariance.m[i][j]) < 1.e-14);
  CeeV0Lite::EventInput oneT0Event = event;
  oneT0Event.sigmaStartTimeNs = 0.01;
  sharedTimeFinder.Reconstruct(tofTracks, oneT0Event, candidates);
  assert(candidates.size() == 1u);
  const CeeV0Lite::Candidate oneT0 = candidates[0];
  CeeV0Lite::EventInput twoT0Event = event;
  twoT0Event.sigmaStartTimeNs = 0.02;
  sharedTimeFinder.Reconstruct(tofTracks, twoT0Event, candidates);
  assert(candidates.size() == 1u);
  const CeeV0Lite::Candidate twoT0 = candidates[0];
  assert(std::fabs(zeroT0.massGeV - oneT0.massGeV) < 1.e-12);
  assert(std::fabs(zeroT0.massGeV - twoT0.massGeV) < 1.e-12);

  // Prove that the common T0 contribution is exactly one rank-one term, not
  // merely that it happens to obey the expected sigma(T0)^2 scaling.  The
  // external half-step derivative below is the same complete fixed-point
  // derivative used by Finder, evaluated with sigma(T0)=0.
  CeeV0Lite::EventInput halfPlusT0Event = zeroT0Event;
  CeeV0Lite::EventInput halfMinusT0Event = zeroT0Event;
  halfPlusT0Event.startTimeNs +=
    0.5 * sharedTimeConfig.sharedTimeDerivativeStepNs;
  halfMinusT0Event.startTimeNs -=
    0.5 * sharedTimeConfig.sharedTimeDerivativeStepNs;
  sharedTimeFinder.Reconstruct(tofTracks, halfPlusT0Event, candidates);
  assert(candidates.size() == 1u);
  const CeeV0Lite::Candidate halfPlusT0 = candidates[0];
  sharedTimeFinder.Reconstruct(tofTracks, halfMinusT0Event, candidates);
  assert(candidates.size() == 1u);
  const CeeV0Lite::Candidate halfMinusT0 = candidates[0];
  double halfPlusValues[7];
  double halfMinusValues[7];
  MotherParameters(halfPlusT0.mother, halfPlusValues);
  MotherParameters(halfMinusT0.mother, halfMinusValues);
  double t0Derivative[7];
  for (int i = 0; i < 7; ++i)
    t0Derivative[i] = (halfPlusValues[i] - halfMinusValues[i]) /
                      sharedTimeConfig.sharedTimeDerivativeStepNs;
  const double oneT0Variance = oneT0Event.sigmaStartTimeNs *
                               oneT0Event.sigmaStartTimeNs;
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      const double observed = oneT0.mother.covariance.m[i][j] -
                              zeroT0.mother.covariance.m[i][j];
      const double expected = t0Derivative[i] * oneT0Variance *
                              t0Derivative[j];
      assert(std::fabs(observed - expected) <
             1.e-13 + 1.e-8 * std::fabs(expected));
    }
  }
  double deltaOne = 0.;
  double deltaTwo = 0.;
  for (int i = 3; i < 7; ++i) {
    deltaOne += oneT0.mother.covariance.m[i][i] -
                zeroT0.mother.covariance.m[i][i];
    deltaTwo += twoT0.mother.covariance.m[i][i] -
                zeroT0.mother.covariance.m[i][i];
  }
  assert(deltaOne > 0.);
  assert(std::fabs(deltaTwo / deltaOne - 4.) < 1.e-8);

  // TOF convergence is required by default.  An intentionally impossible
  // one-iteration tolerance rejects the candidate; the diagnostic override
  // preserves it while marking tofConverged=false.
  CeeV0Lite::FinderConfig nonConvergedTofConfig;
  assert(nonConvergedTofConfig.requireTofConvergence);
  nonConvergedTofConfig.momentumMode = CeeV0Lite::kTofMomentum;
  nonConvergedTofConfig.maximumTofIterations = 1;
  nonConvergedTofConfig.tofVertexConvergenceCm = 0.;
  nonConvergedTofConfig.tofTimeConvergenceNs = 0.;
  CeeV0Lite::Finder strictTofFinder(nonConvergedTofConfig);
  strictTofFinder.Reconstruct(tofTracks, event, candidates);
  assert(candidates.empty());
  nonConvergedTofConfig.requireTofConvergence = false;
  strictTofFinder.SetConfig(nonConvergedTofConfig);
  strictTofFinder.Reconstruct(tofTracks, event, candidates);
  assert(candidates.size() == 1u);
  assert(!candidates[0].tofConverged);
  assert(candidates[0].motherTopo.status==CeeV0Lite::kProductionPvResponseFailed);

  // End-to-end future-field Lambda: build the two measured reference states
  // by helix transport from a common decay vertex, then reconstruct it.
  const double futureFieldTesla = 0.5;
  std::vector<CeeV0Lite::TrackInput> curvedTracks;
  curvedTracks.push_back(MakeCurvedDaughter(
    decay.vertex, decay.protonMomentum, kProtonMass, 10, +1,
    CeeV0Lite::kProtonHypothesis, futureFieldTesla, 20.));
  curvedTracks.push_back(MakeCurvedDaughter(
    decay.vertex, decay.pionMomentum, kPionMass, 11, -1,
    CeeV0Lite::kPionHypothesis, futureFieldTesla, 20.));
  CeeV0Lite::FinderConfig curvedConfig;
  curvedConfig.momentumMode = CeeV0Lite::kTrackingMomentum;
  curvedConfig.magneticFieldByTesla = futureFieldTesla;
  CeeV0Lite::Finder curvedFinder(curvedConfig);
  curvedFinder.Reconstruct(curvedTracks, event, candidates);
  assert(candidates.size() == 1u);
  CheckCandidate(candidates[0], decay.vertex, 5.e-4, 2.e-3);

  // K0S -> pi+ pi- is covered by both momentum sources as well.
  const SyntheticDecay k0Decay = MakeBoostedK0Short();
  const double k0ParentFlightTime = CeeV0Lite::Norm(k0Decay.vertex) /
                                    (k0Decay.parentBeta * kC);
  trackingTracks.clear();
  trackingTracks.push_back(MakeDaughter(
    k0Decay.vertex, k0Decay.protonMomentum, kPionMass, 2, +1,
    CeeV0Lite::kPionHypothesis, true, k0ParentFlightTime));
  trackingTracks.push_back(MakeDaughter(
    k0Decay.vertex, k0Decay.pionMomentum, kPionMass, 3, -1,
    CeeV0Lite::kPionHypothesis, true, k0ParentFlightTime));
  trackingConfig.reconstructLambda = false;
  trackingConfig.reconstructK0Short = true;
  trackingFinder.SetConfig(trackingConfig);
  trackingFinder.Reconstruct(trackingTracks, event, candidates);
  assert(candidates.size() == 1u);
  CheckK0ShortCandidate(candidates[0], k0Decay.vertex, 2.e-4, 1.e-4);

  tofTracks.clear();
  tofTracks.push_back(MakeDaughter(
    k0Decay.vertex, k0Decay.protonMomentum, kPionMass, 2, +1,
    CeeV0Lite::kPionHypothesis, false, k0ParentFlightTime));
  tofTracks.push_back(MakeDaughter(
    k0Decay.vertex, k0Decay.pionMomentum, kPionMass, 3, -1,
    CeeV0Lite::kPionHypothesis, false, k0ParentFlightTime));
  tofConfig.reconstructLambda = false;
  tofConfig.reconstructK0Short = true;
  tofFinder.SetConfig(tofConfig);
  tofFinder.Reconstruct(tofTracks, event, candidates);
  assert(candidates.size() == 1u);
  CheckK0ShortCandidate(candidates[0], k0Decay.vertex, 2.e-3, 1.e-3);
  assert(candidates[0].tofConverged);
  assert(std::fabs(candidates[0].decayTimeNs - k0ParentFlightTime) < 1.e-3);

  std::cout << "tracking, zero-field TOF and finite-field Lambda/K0S tests passed"
            << std::endl;
  return 0;
}
