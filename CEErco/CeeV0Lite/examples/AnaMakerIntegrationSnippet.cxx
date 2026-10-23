// Integration snippet for AnaMaker.cxx. It documents the intended boundary;
// it is not built by the standalone Makefile because the CEE classes live in
// the experiment software environment.

#include "CeeV0LiteCeeAdapter.h"
#include "CeeV0LiteFinder.h"

#include "CeeTpcTrack.h"
#include "CeeTpcVertex.h"

#include <vector>

// Call this before TOF matching and build the projection line exclusively from
// the returned state.  The same adapterConfig must then be passed to
// AddV0Track; this prevents matching a legacy GetX/P line and fitting a
// different persisted vertex-seed state.
static bool SelectV0MatchingLine(
    const CeeTpcTrack& tpcTrack,
    const CeeV0Lite::CeeAdapterConfig& adapterConfig,
    CeeV0Lite::Vec3& positionCm,
    CeeV0Lite::Vec3& direction,
    bool& usedVertexSeed)
{
  CeeV0Lite::CeeTrackReferenceState reference;
  CeeV0Lite::CeeAdapterStatus status;
  if (!CeeV0Lite::SelectTrackReferenceState(
        tpcTrack, CeeV0Lite::kCeeZeroFieldTofMode, reference,
        adapterConfig, &status)) return false;
  positionCm = reference.positionCm;
  direction = reference.direction;
  usedVertexSeed = reference.usedVertexSeed;
  return true;
}

// Call this after matching the line selected above to one unique iTOF/eTOF
// hit.  The hit classes share the same (position,time) boundary here.
static bool AddV0Track(const CeeTpcTrack& tpcTrack,
                       int reconstructedTrackIndex,
                       const CeeV0Lite::Vec3& matchedTofPositionCm,
                       double matchedTofTimeNs,
                       unsigned int pidMask,
                       const CeeV0Lite::CeeAdapterConfig& adapter,
                       std::vector<CeeV0Lite::TrackInput>& output)
{
  CeeV0Lite::CeeTofMeasurement tof;
  tof.available = true;
  tof.positionCm = matchedTofPositionCm;
  tof.timeNs = matchedTofTimeNs;

  CeeV0Lite::TrackInput converted;
  CeeV0Lite::CeeAdapterStatus status;
  if (!CeeV0Lite::MakeTrackInput(
        tpcTrack, reconstructedTrackIndex,
        CeeV0Lite::kCeeZeroFieldTofMode, tof, pidMask, converted,
        adapter, &status)) return false;
  output.push_back(converted);
  return true;
}

static bool ReconstructV0s(const CeeTpcVertex& tpcVertex,
                           bool simulationData,
                           double calibratedT0Ns,
                           double calibratedT0SigmaNs,
                           const std::vector<CeeV0Lite::TrackInput>& tracks,
                           const CeeV0Lite::CeeAdapterConfig& adapter,
                           std::vector<CeeV0Lite::Candidate>& candidates)
{
  CeeV0Lite::EventInput event;
  CeeV0Lite::CeeAdapterStatus status;
  if (!CeeV0Lite::MakePrimaryVertex(tpcVertex, event.primaryVertex, 0,
                                    adapter, &status)) return false;
  CeeV0Lite::FinderConfig config;
  config.momentumMode = CeeV0Lite::kTofMomentum;
  config.magneticFieldByTesla = 0.;
  config.reconstructK0Short = false;
  config.reconstructLambda = true;
  config.reconstructAntiLambda = false;
  config.useLambdaDcaRoleAssignment = true;
  if (simulationData) {
    // MC TOF times are referenced to the exact generator collision time.  This
    // is a coordinate convention, not information from a CeeT0Hit detector.
    event.startTimeNs = 0.;
    event.sigmaStartTimeNs = 0.;
    config.additionalDecayTimeSigmaNs = 0.;
    config.propagateSharedTimeUncertainty = false;
  }
  else {
    event.startTimeNs = calibratedT0Ns;
    event.sigmaStartTimeNs = calibratedT0SigmaNs;
    config.additionalDecayTimeSigmaNs = 0.; // set only from a calibrated model
    config.propagateSharedTimeUncertainty = calibratedT0SigmaNs > 0.;
  }

  CeeV0Lite::Finder finder(config);
  finder.Reconstruct(tracks, event, candidates);
  return true;
}

// Configure this once per job and share it between SelectV0MatchingLine,
// AddV0Track and ReconstructV0s. Replace the fallbacks with CEE residual and
// timing calibrations; a valid seed+6x6 covariance is consumed atomically.
static CeeV0Lite::CeeAdapterConfig MakeCeeAdapterConfig()
{
  CeeV0Lite::CeeAdapterConfig adapter;
  adapter.fallbackTrackSigmaXcm = 0.10;
  adapter.fallbackTrackSigmaYcm = 0.10;
  adapter.fallbackTrackSigmaZcm = 0.15;
  adapter.fallbackDirectionSigmaRad = 0.003;
  adapter.fallbackTofTimeSigmaNs = 0.10;
  adapter.fallbackPathLengthSigmaCm = 0.5;
  adapter.lambdaRoleDcaBoundaryCm = 1.2;
  return adapter;
}

// Example PID masks when building the event track list:
//   pion-like   -> CeeV0Lite::kPionHypothesis
//   proton-like -> CeeV0Lite::kProtonHypothesis
//   unresolved  -> CeeV0Lite::kPionHypothesis |
//                  CeeV0Lite::kProtonHypothesis
//
// The returned Candidate already contains mass, mass error, decay vertex,
// momentum, daughter DCA, geometric chi2/NDF, total topological chi2/NDF,
// the separate PV-constraint chi2/NDF, decay length/significance, pointing
// cosine, raw daughter-to-PV chi2, its sqrt(chi2) significance, TOF iteration
// count and daughter ids.  In 0.6.x topologicalChi2Ndf is the total
// (geometricChi2 + primaryVertexConstraintChi2)/(geometricNdf + 2).
// Use minimumDaughterChi2FromPv for a raw-chi2 cut (the STAR-like preset uses
// 18.6).  minimumDaughterDeviationFromPv keeps the old sqrt(chi2) meaning only
// for backward compatibility.
