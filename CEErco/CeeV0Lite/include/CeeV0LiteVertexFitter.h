#ifndef CEE_V0_LITE_VERTEX_FITTER_H
#define CEE_V0_LITE_VERTEX_FITTER_H

#include "CeeV0LiteTransport.h"

namespace CeeV0Lite
{

struct VertexFitterConfig
{
  VertexFitterConfig();

  int maximumIterations;
  double convergenceCm;
  // Reject ill-conditioned two-track vertices.  For unit directions this is
  // the minimum opening angle modulo pi (parallel and anti-parallel are both
  // degenerate).
  double minimumOpeningAngleRad;
  bool requireConvergence;
  // Retained for source compatibility with 0.4.x.  The 0.5 batch fit has no
  // order-dependent prior and does not use this value.
  double initialVertexVarianceCm2;
  double measurementVarianceFloorCm2;
  // Central-difference steps used by RecalculateMotherCovariance().  Each
  // perturbation reruns the complete POCA/vertex solve, thereby retaining the
  // path derivative that a fixed-path transport covariance does not contain.
  double motherJacobianPositionStepCm;
  double motherJacobianMomentumRelativeStep;
  double motherJacobianMinimumMomentumStepGeV;
};

enum VertexFitStatus
{
  kVertexFitNotRun = 0,
  kVertexFitConverged,
  kVertexFitDegenerateGeometry,
  kVertexFitTransportFailure,
  kVertexFitNumericalFailure,
  kVertexFitMaximumIterations
};

struct VertexFitResult
{
  VertexFitResult();

  bool valid;
  bool converged;
  VertexFitStatus status;
  int iterations;
  Vec3 vertexCm;
  Mat3 covarianceCm2;
  TrackState daughterAtVertex[2];
  MotherState mother;
  double daughterDcaCm;
  double chi2;
  int ndf;
};

class VertexFitter
{
 public:
  VertexFitter(const VertexFitterConfig& fitterConfig = VertexFitterConfig(),
               const TransportConfig& transportConfig = TransportConfig());

  bool Fit(const TrackState& first,
           const TrackState& second,
           double magneticFieldByTesla,
           VertexFitResult& result) const;

  bool RecalculateMotherCovariance(const TrackState& first,
                                   const TrackState& second,
                                   double magneticFieldByTesla,
                                   MotherState& mother) const;

 private:
  VertexFitterConfig fFitterConfig;
  TransportConfig fTransportConfig;
};

} // namespace CeeV0Lite

#endif
