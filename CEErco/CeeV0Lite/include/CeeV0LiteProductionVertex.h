#ifndef CEE_V0_LITE_PRODUCTION_VERTEX_H
#define CEE_V0_LITE_PRODUCTION_VERTEX_H

#include "CeeV0LiteTypes.h"

namespace CeeV0Lite {

struct ProductionVertexConfig
{
  ProductionVertexConfig();
  int maximumIterations;
  double relativeStateTolerance;
  double constraintToleranceCm;
  // Finder's full-TOF-solve derivative with respect to measured PV position.
  double pvDerivativeStepCm;
  double pvDerivativeStabilityTolerance;
  // Charged fit is a local decay-flight solution, not a winding-number fit.
  double maximumAbsFlightPathCm;
  double maximumAbsBendingAngleRad;
};

// The mother's covariance here is its MARGINAL covariance, already including
// any response to the uncertain PV. Cross block is Cov(mother, PV).
// Zero cross block is valid only under an explicit independence assumption.
struct ProductionVertexInput
{
  ProductionVertexInput();
  MotherState mother;
  PrimaryVertex vertex;
  double motherVertexCovariance[7][3];
};

// Fit two pointing constraints on a joint 10D Gaussian prior. Re-linearize
// against the SAME prior, never re-counting the PV at each iteration.
// Compatibility overload: omitting the field is allowed only for neutral
// mothers. A nonzero charge here fails with kProductionMissingMagneticField.
bool FitProductionVertex(const ProductionVertexInput& input,
                         const ProductionVertexConfig& config,
                         double geometricChi2, int geometricNdf,
                         ProductionState& output, int motherCharge = 0);

// Explicit uniform CEE field B=(0,By,0). Charged mothers use a helix; By=0 is
// permitted only when explicitly supplied. Inputs are never modified.
// S is signed arc length / |p|, not the straight-line endpoint projection.
bool FitProductionVertex(const ProductionVertexInput& input,
                         const ProductionVertexConfig& config,
                         double geometricChi2, int geometricNdf,
                         ProductionState& output, int motherCharge,
                         double magneticFieldByTesla);

// Deterministic endpoint transport in the field stored in the state, including
// momentum rotation and all position/momentum derivatives with respect to S.
// Repeated calls to the same endpoint are idempotent.
bool TransportToDecayVertex(ProductionState& state);
bool TransportToProductionVertex(ProductionState& state);
bool GetFlightObservables(const ProductionState& state,
                           FlightObservables& output);

} // namespace CeeV0Lite
#endif
