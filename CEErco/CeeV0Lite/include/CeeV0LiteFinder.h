#ifndef CEE_V0_LITE_FINDER_H
#define CEE_V0_LITE_FINDER_H

#include "CeeV0LiteVertexFitter.h"
#include "CeeV0LiteProductionVertex.h"

#include <cstddef>
#include <vector>

namespace CeeV0Lite
{

struct ChannelCuts
{
  ChannelCuts();

  bool enabled;
  double minimumMassGeV;
  double maximumMassGeV;
  double maximumDaughterDcaCm;
  double maximumGeometricChi2Ndf;
  double maximumTopologicalChi2Ndf;
  double minimumDecayLengthCm;
  double minimumDecayLengthSignificance;
  double minimumPointingCosine;
  // Raw KFParticle-style daughter chi2 with respect to the primary vertex.
  double minimumDaughterChi2FromPv;
  // Deprecated compatibility cut on sqrt(chi2).  It is applied in addition to
  // the raw-chi2 cut and defaults to zero.
  double minimumDaughterDeviationFromPv;
  // By default a candidate must complete all +/- refits used by the numerical
  // mother Jacobian.  Setting this false retains that one exact failure class
  // for QA; its covariance-dependent observables remain invalid/NaN.  All
  // other covariance failures are still rejected.
  bool requireCompleteMotherJacobian;
};

struct FinderConfig
{
  FinderConfig();

  MomentumMode momentumMode;
  // Signed uniform field along CEE y; beam remains along z. Positive is +y.
  double magneticFieldByTesla;
  bool reconstructK0Short;
  bool reconstructLambda;
  bool reconstructAntiLambda;
  bool useLambdaDcaRoleAssignment;

  bool requireMeasuredChargeInMagneticField;
  bool allowTofMomentumInMagneticField;
  double zeroFieldThresholdTesla;

  // v0.6.5: once per event/track, use measured track--PV DCA as the magnetic
  // reference. Requires measured charge, momentum and full 6D covariance.
  // Disable only to replay the stored-reference convention. Zero B bypasses.
  bool usePrimaryVertexDcaReference;
  TransportConfig referenceTransportConfig;

  int maximumTofIterations;
  double tofVertexConvergenceCm;
  double tofTimeConvergenceNs;
  double maximumBeta;
  // Both quantities below are common to the two daughters.  Defaults are zero
  // and disabled, corresponding to an exact t=0 generator time in simulation.
  // Real-data analysis must enable and calibrate them explicitly.  They are not
  // inserted into two independent track covariances; one shared derivative
  // is propagated through the complete TOF fixed-point solution instead.
  double additionalDecayTimeSigmaNs;
  bool propagateSharedTimeUncertainty;
  double sharedTimeDerivativeStepNs;
  double sharedTimeDerivativeStabilityTolerance;
  double seedParentBeta;
  bool requireTofConvergence;
  double defaultRelativeTrackingMomentumError;

  std::size_t maximumCandidates;
  ChannelCuts k0Cuts;
  ChannelCuts lambdaCuts;
  TransportConfig transportConfig;
  VertexFitterConfig vertexFitterConfig;
  // Build an additional 8D PV-constrained state after legacy cuts.
  // Failed extra fits are flagged and do not discard the original candidate.
  bool buildProductionState;
  ProductionVertexConfig productionVertexConfig;
};

enum TrackReferenceStatus
{
  kReferenceUnchanged = 0,
  kReferenceAtPvDca = 1,
  kReferenceInvalidInput = 2,
  kReferenceTransportFailed = 3
};

struct TrackReferenceDiagnostic
{
  TrackReferenceDiagnostic();
  int trackId;
  TrackReferenceStatus status;
  ClosestApproachResult closestApproach;
  double dcaCm;
};

// Same preparation used by Finder::Reconstruct; no mass hypothesis needed.
// Returns false for rejected tracks; output is unchanged on failure.
bool PrepareTrackReference(const TrackInput& input, const EventInput& event,
                           const FinderConfig& config, TrackInput& output,
                           TrackReferenceDiagnostic* diagnostic = 0);

class Finder
{
 public:
  explicit Finder(const FinderConfig& config = FinderConfig());

  void SetConfig(const FinderConfig& config) { fConfig = config; }
  const FinderConfig& GetConfig() const { return fConfig; }

  void Reconstruct(const std::vector<TrackInput>& tracks,
                   const EventInput& event,
                   std::vector<Candidate>& candidates,
                   std::vector<TrackReferenceDiagnostic>* diagnostics = 0) const;

 private:
  FinderConfig fConfig;
};

} // namespace CeeV0Lite

#endif
