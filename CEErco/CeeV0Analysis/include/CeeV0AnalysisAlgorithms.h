#ifndef CEE_V0_ANALYSIS_ALGORITHMS_H
#define CEE_V0_ANALYSIS_ALGORITHMS_H

#include <cstddef>
#include <vector>

namespace CeeV0Analysis
{

enum TofDetector
{
  kNoTofDetector = 0,
  kItofDetector = 1,
  kEtofDetector = 2
};

struct ProjectionTrack
{
  ProjectionTrack();

  int id;
  double positionCm[3];
  double direction[3];
};

struct TofHitData
{
  TofHitData();

  int id;
  TofDetector detector;
  // CeeiTOFHit::GetDetectorID()/100. Ignored for eTOF.
  int module;
  double positionCm[3];
  double rawTime;
};

struct TofMatch
{
  TofMatch();

  std::size_t trackIndex;
  std::size_t hitIndex;
  int trackId;
  int hitId;
  TofDetector detector;
  double distanceCm;
};

/**
 * Spatial residual between a straight TPC track and a TOF hit.
 *
 * The iTOF module-plane convention is copied from the existing CEE AnaMaker:
 * side modules 0--8 and 15--24 use an x plane (residual in yz), while forward
 * modules 9--14 use a z plane (residual in xy). eTOF always uses a z plane.
 * A negative line parameter is allowed because a zero-field fit may store the
 * otherwise equivalent opposite direction. The finder later orients the track
 * toward the assigned TOF hit.
 */
double ProjectionDistanceCm(const ProjectionTrack& track,
                            const TofHitData& hit);

/**
 * Build all accepted track-hit edges and find a globally unique assignment.
 * The optimization first maximizes the number of matched tracks and then
 * minimizes the sum of distance/maximumDistance. Normalization prevents the
 * raw cm scale of one TOF detector from dominating the other detector.
 * A track and a physical hit can each appear at most once.
 */
std::vector<TofMatch> AssignUniqueTofHits(
    const std::vector<ProjectionTrack>& tracks,
    const std::vector<TofHitData>& hits,
    double maximumItofDistanceCm,
    double maximumEtofDistanceCm);

enum T0Source
{
  kFixedT0 = 0,
  kT0Hit = 1
};

enum T0Aggregation
{
  kFirstT0 = 0,
  kEarliestT0 = 1,
  kMeanT0 = 2,
  kMedianT0 = 3
};

struct T0SelectionConfig
{
  T0SelectionConfig();

  T0Source source;
  T0Aggregation aggregation;
  double fixedTimeNs;
  double rawToNsScale;
  double offsetNs;
  double sigmaNs;
};

struct T0SelectionResult
{
  T0SelectionResult();

  bool valid;
  int numberOfHits;
  double selectedRawTime;
  double timeNs;
  double sigmaNs;
};

T0SelectionResult SelectEventT0(const std::vector<double>& rawHitTimes,
                                const T0SelectionConfig& config);

const char* TofDetectorName(TofDetector detector);

} // namespace CeeV0Analysis

#endif
