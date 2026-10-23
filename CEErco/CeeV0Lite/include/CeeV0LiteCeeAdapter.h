#ifndef CEE_V0_LITE_CEE_ADAPTER_H
#define CEE_V0_LITE_CEE_ADAPTER_H

#include "CeeV0LiteTypes.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

namespace CeeV0Lite
{

// In the present zero-field data the TPC momentum vector is used only to
// define the track direction; its magnitude must not enter the V0 fit.
// In magnetic-field mode the fitted momentum and its full covariance are
// used, so a persisted vertex-seed state is mandatory.
enum CeeTrackInputMode
{
  kCeeZeroFieldTofMode = 0,
  kCeeMagneticFieldTrackingMode = 1
};

enum CeeAdapterStatus
{
  kCeeAdapterSuccess = 0,
  kCeeAdapterSuccessVertexSeed,
  kCeeAdapterSuccessLegacyFallback,
  kCeeAdapterNoValidVertex,
  kCeeAdapterInvalidVertexUncertainty,
  kCeeAdapterInvalidTrackState,
  kCeeAdapterTrackFitNotConverged,
  kCeeAdapterMissingVertexSeed,
  kCeeAdapterInvalidVertexSeed,
  kCeeAdapterInvalidTrackingCovariance,
  kCeeAdapterMissingCharge,
  kCeeAdapterMissingTof,
  kCeeAdapterInvalidTof,
  kCeeAdapterInvalidHypothesisMask
};

inline const char* CeeAdapterStatusMessage(CeeAdapterStatus status)
{
  switch (status) {
    case kCeeAdapterSuccess: return "success";
    case kCeeAdapterSuccessVertexSeed:
      return "success: persisted vertex-seed reference state";
    case kCeeAdapterSuccessLegacyFallback:
      return "success: legacy track-state fallback";
    case kCeeAdapterNoValidVertex: return "no finite primary vertex";
    case kCeeAdapterInvalidVertexUncertainty:
      return "invalid primary-vertex uncertainty";
    case kCeeAdapterInvalidTrackState: return "invalid track state";
    case kCeeAdapterTrackFitNotConverged: return "track fit did not converge";
    case kCeeAdapterMissingVertexSeed: return "missing persisted vertex seed";
    case kCeeAdapterInvalidVertexSeed: return "invalid persisted vertex seed";
    case kCeeAdapterInvalidTrackingCovariance:
      return "invalid vertex-seed covariance";
    case kCeeAdapterMissingCharge: return "missing fitted charge";
    case kCeeAdapterMissingTof: return "missing TOF measurement";
    case kCeeAdapterInvalidTof: return "invalid TOF measurement";
    case kCeeAdapterInvalidHypothesisMask: return "empty mass-hypothesis mask";
  }
  return "unknown adapter status";
}

inline bool CeeAdapterStatusIsSuccess(CeeAdapterStatus status)
{
  return status == kCeeAdapterSuccess ||
         status == kCeeAdapterSuccessVertexSeed ||
         status == kCeeAdapterSuccessLegacyFallback;
}

struct CeeTofMeasurement
{
  CeeTofMeasurement()
    : available(false), timeNs(0.), sigmaTimeNs(0.),
      sigmaPathLengthCm(0.) {}

  bool available;
  Vec3 positionCm;
  // Absolute detector time in ns. EventInput::startTimeNs is subtracted by
  // the finder; this adapter intentionally does not guess the T0 convention.
  double timeNs;
  // Non-positive values request the configured fallback resolutions.
  double sigmaTimeNs;
  double sigmaPathLengthCm;
};

struct CeeAdapterConfig
{
  CeeAdapterConfig()
    : fallbackTrackSigmaXcm(0.10), fallbackTrackSigmaYcm(0.10),
      fallbackTrackSigmaZcm(0.15), fallbackVertexSigmaXcm(0.10),
      fallbackVertexSigmaYcm(0.10), fallbackVertexSigmaZcm(0.15),
      fallbackDirectionSigmaRad(0.003), fallbackTofTimeSigmaNs(0.10),
      fallbackPathLengthSigmaCm(0.5), lambdaRoleDcaBoundaryCm(1.2),
      minimumMomentumGeV(1.e-9), minimumCovarianceDiagonal(1.e-16),
      fittedChargeThreshold(0.25), allowMissingVertexRank(true),
      allowVertexErrorFallback(true), allowLegacyStateInZeroField(true),
      requireConvergedTrackFit(false), requirePositiveDefiniteCovariance(true)
  {}

  double fallbackTrackSigmaXcm;
  double fallbackTrackSigmaYcm;
  double fallbackTrackSigmaZcm;
  double fallbackVertexSigmaXcm;
  double fallbackVertexSigmaYcm;
  double fallbackVertexSigmaZcm;
  double fallbackDirectionSigmaRad;
  double fallbackTofTimeSigmaNs;
  double fallbackPathLengthSigmaCm;
  double lambdaRoleDcaBoundaryCm;
  double minimumMomentumGeV;
  double minimumCovarianceDiagonal;
  double fittedChargeThreshold;
  bool allowMissingVertexRank;
  bool allowVertexErrorFallback;
  bool allowLegacyStateInZeroField;
  bool requireConvergedTrackFit;
  bool requirePositiveDefiniteCovariance;
};

// The one authoritative reference-state decision shared by TOF matching and
// MakeTrackInput.  A caller must use positionCm/direction from this object for
// hit matching instead of independently checking HasVertexSeed().
struct CeeTrackReferenceState
{
  CeeTrackReferenceState()
    : sigmaDirectionRad(0.), usedVertexSeed(false),
      hasTrackingCovariance(false) {}

  Vec3 positionCm;
  Vec3 direction;
  Mat3 positionCovarianceCm2;
  double sigmaDirectionRad;
  bool usedVertexSeed;
  bool hasTrackingCovariance;
  Mat6 trackingCovariance;
};

namespace CeeAdapterDetail
{

inline void SetStatus(CeeAdapterStatus value, CeeAdapterStatus* output)
{
  if (output) *output = value;
}

inline bool IsPositiveFinite(double value)
{
  return std::isfinite(value) && value > 0.;
}

inline bool IsFiniteCovariance(const Mat6& covariance,
                               double minimumDiagonal,
                               bool requirePositiveDefinite)
{
  double maximumDiagonal = 0.;
  for (int row = 0; row < 6; ++row) {
    for (int column = 0; column < 6; ++column) {
      if (!std::isfinite(covariance.m[row][column])) return false;
    }
    if (!(covariance.m[row][row] > minimumDiagonal)) return false;
    maximumDiagonal = std::max(maximumDiagonal, covariance.m[row][row]);
  }
  if (!requirePositiveDefinite) return true;

  // Cholesky test with a small scale-dependent tolerance. Small negative
  // pivots caused solely by roundoff are floored; a genuinely indefinite
  // matrix is rejected before it reaches the fitter.
  double lower[6][6] = {{0.}};
  const double tolerance = std::max(minimumDiagonal,
                                    maximumDiagonal * 1.e-12);
  for (int row = 0; row < 6; ++row) {
    for (int column = 0; column <= row; ++column) {
      double value = covariance.m[row][column];
      for (int k = 0; k < column; ++k)
        value -= lower[row][k] * lower[column][k];
      if (row == column) {
        if (value < -tolerance) return false;
        lower[row][column] = std::sqrt(std::max(value, tolerance));
      }
      else {
        lower[row][column] = value / lower[column][column];
      }
    }
  }
  return true;
}

inline bool FillFallbackPositionCovariance(const CeeAdapterConfig& config,
                                           Mat3& covariance)
{
  if (!IsPositiveFinite(config.fallbackTrackSigmaXcm) ||
      !IsPositiveFinite(config.fallbackTrackSigmaYcm) ||
      !IsPositiveFinite(config.fallbackTrackSigmaZcm)) return false;
  covariance = Mat3();
  covariance.m[0][0] = config.fallbackTrackSigmaXcm *
                       config.fallbackTrackSigmaXcm;
  covariance.m[1][1] = config.fallbackTrackSigmaYcm *
                       config.fallbackTrackSigmaYcm;
  covariance.m[2][2] = config.fallbackTrackSigmaZcm *
                       config.fallbackTrackSigmaZcm;
  return true;
}

template <typename MemberFunction> struct FirstArgument;

template <typename ClassT, typename ReturnT, typename ArgumentT>
struct FirstArgument<ReturnT (ClassT::*)(ArgumentT) const>
{
  typedef typename std::remove_cv<
    typename std::remove_reference<ArgumentT>::type>::type type;
};

template <class CeeTrackT>
CeeAdapterStatus CopyVertexSeed(const CeeTrackT& track,
                                Vec3& position,
                                Vec3& momentum,
                                Mat6& covariance,
                                const CeeAdapterConfig& config)
{
  if (!track.HasVertexSeed()) return kCeeAdapterMissingVertexSeed;

  const auto& sourcePosition = track.GetVertexSeedPos();
  const auto& sourceMomentum = track.GetVertexSeedMom();
  position = Vec3(sourcePosition.X(), sourcePosition.Y(), sourcePosition.Z());
  momentum = Vec3(sourceMomentum.X(), sourceMomentum.Y(), sourceMomentum.Z());
  if (!IsFinite(position) || !IsFinite(momentum) ||
      !(Norm(momentum) > config.minimumMomentumGeV))
    return kCeeAdapterInvalidVertexSeed;

  typedef typename FirstArgument<
    decltype(&CeeTrackT::GetVertexSeedCov)>::type SourceCovariance;
  SourceCovariance sourceCovariance;
  track.GetVertexSeedCov(sourceCovariance);
  if (sourceCovariance.GetNrows() < 6)
    return kCeeAdapterInvalidTrackingCovariance;

  covariance = Mat6();
  for (int row = 0; row < 6; ++row) {
    for (int column = 0; column < 6; ++column) {
      const double a = sourceCovariance(row, column);
      const double b = sourceCovariance(column, row);
      covariance.m[row][column] = 0.5 * (a + b);
    }
  }
  return IsFiniteCovariance(covariance, config.minimumCovarianceDiagonal,
                            config.requirePositiveDefiniteCovariance)
    ? kCeeAdapterSuccess : kCeeAdapterInvalidTrackingCovariance;
}

inline double DirectionSigmaFromMomentumCovariance(const Vec3& momentum,
                                                   const Mat6& covariance,
                                                   double fallback)
{
  const double momentum2 = Norm2(momentum);
  if (!(momentum2 > 0.)) return fallback;
  const Vec3 unit = momentum / std::sqrt(momentum2);
  double trace = 0.;
  double longitudinal = 0.;
  for (int row = 0; row < 3; ++row) {
    trace += covariance.m[row + 3][row + 3];
    for (int column = 0; column < 3; ++column)
      longitudinal += unit[row] * covariance.m[row + 3][column + 3] *
                      unit[column];
  }
  const double transverse = trace - longitudinal;
  if (!(transverse > 0.) || !std::isfinite(transverse)) return fallback;
  const double sigma = std::sqrt(transverse / (2. * momentum2));
  return IsPositiveFinite(sigma) ? sigma : fallback;
}

inline double ValidOrFallbackError(double error, double fallback,
                                   bool allowFallback, bool& valid)
{
  if (IsPositiveFinite(error)) return error;
  if (allowFallback && IsPositiveFinite(fallback)) return fallback;
  valid = false;
  return 0.;
}

} // namespace CeeAdapterDetail

template <class CeeTrackT>
bool SelectTrackReferenceState(
    const CeeTrackT& track,
    CeeTrackInputMode mode,
    CeeTrackReferenceState& output,
    const CeeAdapterConfig& config = CeeAdapterConfig(),
    CeeAdapterStatus* status = 0)
{
  using namespace CeeAdapterDetail;
  output = CeeTrackReferenceState();

  const bool magneticMode = mode == kCeeMagneticFieldTrackingMode;
  if (!magneticMode && mode != kCeeZeroFieldTofMode) {
    SetStatus(kCeeAdapterInvalidTrackState, status);
    return false;
  }
  if (config.requireConvergedTrackFit && !track.GetIsFitConverged()) {
    SetStatus(kCeeAdapterTrackFitNotConverged, status);
    return false;
  }

  Vec3 seedPosition;
  Vec3 seedMomentum;
  Mat6 seedCovariance;
  const CeeAdapterStatus seedStatus = CopyVertexSeed(
    track, seedPosition, seedMomentum, seedCovariance, config);
  const bool validSeed = seedStatus == kCeeAdapterSuccess;

  if (validSeed) {
    output.positionCm = seedPosition;
    output.direction = seedMomentum;
    output.positionCovarianceCm2 = PositionBlock(seedCovariance);
    output.sigmaDirectionRad = DirectionSigmaFromMomentumCovariance(
      seedMomentum, seedCovariance, config.fallbackDirectionSigmaRad);
    output.usedVertexSeed = true;
    output.hasTrackingCovariance = true;
    output.trackingCovariance = seedCovariance;
  }
  else {
    if (magneticMode || !config.allowLegacyStateInZeroField) {
      SetStatus(seedStatus, status);
      return false;
    }
    // Fall back as one complete reference state.  Never combine the persisted
    // seed position/direction with a covariance assigned to GetX/Y/Z.
    output.positionCm = Vec3(track.GetX(), track.GetY(), track.GetZ());
    output.direction = Vec3(track.GetPx(), track.GetPy(), track.GetPz());
    if (!FillFallbackPositionCovariance(config,
                                        output.positionCovarianceCm2)) {
      SetStatus(kCeeAdapterInvalidTrackingCovariance, status);
      return false;
    }
    output.sigmaDirectionRad = config.fallbackDirectionSigmaRad;
    output.usedVertexSeed = false;
    output.hasTrackingCovariance = false;
  }

  if (!IsFinite(output.positionCm) || !IsFinite(output.direction) ||
      !(Norm(output.direction) > config.minimumMomentumGeV) ||
      !IsPositiveFinite(output.sigmaDirectionRad)) {
    SetStatus(kCeeAdapterInvalidTrackState, status);
    return false;
  }

  SetStatus(output.usedVertexSeed ? kCeeAdapterSuccessVertexSeed
                                  : kCeeAdapterSuccessLegacyFallback,
            status);
  return true;
}

template <class CeeTrackT>
bool MakeTrackInput(const CeeTrackT& track,
                    int reconstructedTrackId,
                    CeeTrackInputMode mode,
                    const CeeTofMeasurement& tof,
                    unsigned int hypothesisMask,
                    TrackInput& output,
                    const CeeAdapterConfig& config = CeeAdapterConfig(),
                    CeeAdapterStatus* status = 0)
{
  using namespace CeeAdapterDetail;
  output = TrackInput();

  if (hypothesisMask == 0u) {
    SetStatus(kCeeAdapterInvalidHypothesisMask, status);
    return false;
  }
  const bool magneticMode = mode == kCeeMagneticFieldTrackingMode;
  CeeTrackReferenceState reference;
  CeeAdapterStatus referenceStatus = kCeeAdapterInvalidTrackState;
  if (!SelectTrackReferenceState(track, mode, reference, config,
                                 &referenceStatus)) {
    SetStatus(referenceStatus, status);
    return false;
  }

  output.id = reconstructedTrackId;
  output.hypothesisMask = hypothesisMask;
  output.lambdaRoleDcaBoundaryCm = config.lambdaRoleDcaBoundaryCm;
  output.qualityChi2 = std::isfinite(track.GetFitChi2()) &&
                       track.GetFitChi2() >= 0. ? track.GetFitChi2() : 0.;
  const double fitNdf = track.GetFitNdf();
  if (std::isfinite(fitNdf) && fitNdf >= 1.) {
    const double maximumInt = static_cast<double>(
      std::numeric_limits<int>::max());
    output.qualityNdf = static_cast<int>(
      std::min(maximumInt, std::floor(fitNdf + 0.5)));
  }
  else {
    output.qualityNdf = 1;
  }

  output.positionCm = reference.positionCm;
  output.direction = reference.direction;
  output.positionCovarianceCm2 = reference.positionCovarianceCm2;
  output.sigmaDirectionRad = reference.sigmaDirectionRad;

  if (magneticMode) {
    const double fittedCharge = track.GetFitCharge();
    output.chargeKnown = std::isfinite(fittedCharge) &&
                         std::fabs(fittedCharge) >=
                           config.fittedChargeThreshold;
    output.charge = output.chargeKnown ? (fittedCharge > 0. ? 1 : -1) : 0;
    if (!output.chargeKnown) {
      SetStatus(kCeeAdapterMissingCharge, status);
      return false;
    }
  }
  else {
    // With B=0 the sign returned by the track fit is not a measured charge.
    // Mark it unknown so it cannot silently veto a p/pi mass hypothesis.
    output.chargeKnown = false;
    output.charge = 0;
  }

  if (magneticMode) {
    output.hasTrackingMomentum = true;
    output.trackingMomentumGeV = reference.direction;
    output.hasTrackingCovariance = true;
    output.trackingCovariance = reference.trackingCovariance;
  }
  else {
    // Deliberately ignore the magnitude returned by GetP* or the seed. In
    // zero field it contains direction information only; TOF supplies |p|.
    output.hasTrackingMomentum = false;
    output.hasTrackingCovariance = false;
  }

  if (tof.available) {
    if (!IsFinite(tof.positionCm) || !std::isfinite(tof.timeNs)) {
      SetStatus(kCeeAdapterInvalidTof, status);
      return false;
    }
    output.hasTof = true;
    output.tofPositionCm = tof.positionCm;
    output.tofTimeNs = tof.timeNs;
    output.sigmaTofTimeNs = IsPositiveFinite(tof.sigmaTimeNs)
      ? tof.sigmaTimeNs : config.fallbackTofTimeSigmaNs;
    output.sigmaPathLengthCm = IsPositiveFinite(tof.sigmaPathLengthCm)
      ? tof.sigmaPathLengthCm : config.fallbackPathLengthSigmaCm;
    if (!IsPositiveFinite(output.sigmaTofTimeNs) ||
        !IsPositiveFinite(output.sigmaPathLengthCm)) {
      SetStatus(kCeeAdapterInvalidTof, status);
      return false;
    }
  }
  else if (!magneticMode) {
    SetStatus(kCeeAdapterMissingTof, status);
    return false;
  }

  SetStatus(referenceStatus, status);
  return true;
}

template <class CeeVertexT>
bool MakePrimaryVertex(const CeeVertexT& vertex,
                       PrimaryVertex& output,
                       std::size_t* selectedIndex = 0,
                       const CeeAdapterConfig& config = CeeAdapterConfig(),
                       CeeAdapterStatus* status = 0)
{
  using namespace CeeAdapterDetail;
  output = PrimaryVertex();

  const auto& x = vertex.GetX();
  const auto& y = vertex.GetY();
  const auto& z = vertex.GetZ();
  const auto& rank = vertex.GetRank();
  const std::size_t coordinateCount = std::min(x.size(),
    std::min(y.size(), z.size()));
  if (coordinateCount == 0u) {
    SetStatus(kCeeAdapterNoValidVertex, status);
    return false;
  }

  bool found = false;
  std::size_t bestIndex = 0u;
  int bestRank = std::numeric_limits<int>::min();
  for (std::size_t i = 0; i < coordinateCount; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) ||
        !std::isfinite(z[i])) continue;
    const bool hasRank = i < rank.size();
    if (!hasRank && !config.allowMissingVertexRank) continue;
    const int candidateRank = hasRank ? rank[i]
                                      : std::numeric_limits<int>::min();
    if (!found || candidateRank > bestRank) {
      found = true;
      bestIndex = i;
      bestRank = candidateRank;
    }
  }
  if (!found) {
    SetStatus(kCeeAdapterNoValidVertex, status);
    return false;
  }

  const auto& xError = vertex.GetXErr();
  const auto& yError = vertex.GetYErr();
  const auto& zError = vertex.GetZErr();
  bool validErrors = true;
  const double sx = ValidOrFallbackError(
    bestIndex < xError.size() ? xError[bestIndex] : 0.,
    config.fallbackVertexSigmaXcm, config.allowVertexErrorFallback,
    validErrors);
  const double sy = ValidOrFallbackError(
    bestIndex < yError.size() ? yError[bestIndex] : 0.,
    config.fallbackVertexSigmaYcm, config.allowVertexErrorFallback,
    validErrors);
  const double sz = ValidOrFallbackError(
    bestIndex < zError.size() ? zError[bestIndex] : 0.,
    config.fallbackVertexSigmaZcm, config.allowVertexErrorFallback,
    validErrors);
  if (!validErrors) {
    SetStatus(kCeeAdapterInvalidVertexUncertainty, status);
    return false;
  }

  output.positionCm = Vec3(x[bestIndex], y[bestIndex], z[bestIndex]);
  output.covarianceCm2 = Mat3();
  output.covarianceCm2.m[0][0] = sx * sx;
  output.covarianceCm2.m[1][1] = sy * sy;
  output.covarianceCm2.m[2][2] = sz * sz;
  if (selectedIndex) *selectedIndex = bestIndex;
  SetStatus(kCeeAdapterSuccess, status);
  return true;
}

// Utility for external state converters that already store a 6x6 covariance
// as 21 lower-triangle values. It does not assert covariance validity; the
// strict MakeTrackInput path above validates CeeTpcTrack vertex-seed states.
inline void SetTrackingCovarianceFromLowerTriangle(
    TrackInput& track, const double lowerTriangle[21])
{
  track.hasTrackingCovariance = true;
  int index = 0;
  for (int row = 0; row < 6; ++row) {
    for (int column = 0; column <= row; ++column, ++index) {
      track.trackingCovariance.m[row][column] = lowerTriangle[index];
      track.trackingCovariance.m[column][row] = lowerTriangle[index];
    }
  }
}

} // namespace CeeV0Lite

#endif
