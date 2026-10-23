#ifndef CEE_V0_LITE_TRANSPORT_H
#define CEE_V0_LITE_TRANSPORT_H

#include "CeeV0LiteTypes.h"

namespace CeeV0Lite
{

struct TransportConfig
{
  TransportConfig();

  int maximumClosestApproachIterations;
  double closestApproachToleranceCm;
  double maximumAbsPathCm;
  // A closest-approach state obtained at the iteration/path boundary is not
  // silently accepted by default.  Set this to false only for diagnostics.
  bool requireClosestApproachConvergence;
  // v0.6.4: exact fixed-path map derivative; false retains the FD reference.
  // This switch does not change the POCA derivative or root solver.
  bool useAnalyticFixedPathJacobian;
  double positionProcessNoiseCm2PerCm;
  double angularProcessNoiseRad2PerCm;
};

struct ClosestApproachResult
{
  ClosestApproachResult();

  bool valid;
  bool converged;
  bool hitPathLimit;
  bool singularDerivative;
  int iterations;
  double signedPathCm;
  // At a true closest approach this stationarity residual, (x-x0).u, is zero.
  double stationarityResidualCm;
};

// Propagate Cartesian (x,y,z,px,py,pz) by signed 3D path length s.
// CEE uses B=(0,By,0), with the beam along z; no coordinate reordering.
// Units: cm, GeV/c and Tesla. Fixed-path covariance uses an analytic Jacobian
// by default, with the original numerical derivative retained as a reference.
bool PropagateState(const TrackState& input,
                    double signedPathCm,
                    double magneticFieldByTesla,
                    const TransportConfig& config,
                    TrackState& output,
                    bool propagateCovariance = true);

// Choose a reference point on the measured helix nearest to pointCm.
// Solve the nominal path once, then propagate the COMPLETE 6D covariance at
// that fixed path. This is a reversible reparameterization, not a point/PV
// measurement update. No point covariance or process noise is added.
// Unlike TransportToClosestApproach, it does not differentiate the POCA path.
// Convergence is mandatory; output is unchanged on failure.
bool ReanchorAtClosestApproach(const TrackState& input,
                              const Vec3& pointCm,
                              double magneticFieldByTesla,
                              const TransportConfig& config,
                              TrackState& output,
                              ClosestApproachResult* diagnostics = 0);

bool TransportToClosestApproach(const TrackState& input,
                                const Vec3& pointCm,
                                double magneticFieldByTesla,
                                const TransportConfig& config,
                                TrackState& output,
                                double& signedPathCm,
                                ClosestApproachResult* diagnostics = 0);

// Reuse the already computed POCA derivative for 7D energy cross-covariance.
bool TransportToClosestApproach(const TrackState& input,
                                const Vec3& pointCm,
                                double magneticFieldByTesla,
                                const TransportConfig& config,
                                TrackState& output,
                                double& signedPathCm,
                                ClosestApproachResult* diagnostics,
                                double (*stateJacobian)[6]);

} // namespace CeeV0Lite

#endif
