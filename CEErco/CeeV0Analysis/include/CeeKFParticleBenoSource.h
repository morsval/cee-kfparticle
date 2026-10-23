#ifndef CEE_KFP_BENO_SOURCE_H
#define CEE_KFP_BENO_SOURCE_H
#include "CeeKFParticleEvent.h"
#include "CeeV0Analyzer.h"
#include <functional>
class TClonesArray;
class CeeBenoTpcTrack;
class CeeBenoTpcTrackCovMatrix;
namespace CeeV0Analysis {
struct BenoInputOptions {
  BenoInputOptions()
      : rejectSharedHits(true), requireUniqueMatch(false), preferItof(true),
        maximumMatchChi2(1.e30), geometryPositionToleranceCm(0.10),
        geometryDirectionToleranceRad(0.003) {}
  bool rejectSharedHits, requireUniqueMatch, preferItof;
  double maximumMatchChi2;
  // Cross-check the documented beam-axis perigee against origin()/globMom().
  // Allows quantization, but rejects mislabeled upstream GenFit fallback data.
  double geometryPositionToleranceCm, geometryDirectionToleranceRad;
};
struct BenoInputStatistics {
  BenoInputStatistics()
      : entries(0), tracksConsidered(0), tracksAccepted(0), invalidGeometry(0),
        missingMatch(0), sharedMatch(0) {}
  unsigned long long entries, tracksConsidered, tracksAccepted;
  unsigned long long invalidGeometry, missingMatch, sharedMatch;
};
// ROOT/list reader for already-produced BenoDst. Uses the EDM dictionaries,
// not FairRootManager / CeeBenoDstReader, and never reruns detector matching.
class CeeKFParticleBenoSource : public CeeKFParticleEventSource {
public:
  using TrackIndexSelector =
      std::function<std::vector<int>(long long, const TClonesArray &)>;
  using StateConverter = std::function<bool(
      const CeeBenoTpcTrack &, const CeeBenoTpcTrackCovMatrix &,
      CeeV0Lite::TrackInput &, std::string &)>;
  explicit CeeKFParticleBenoSource(
      const AnalysisConfig &, const BenoInputOptions & = BenoInputOptions());
  ~CeeKFParticleBenoSource();
  void SetTrackIndexSelector(TrackIndexSelector);
  // Override only for a documented alternative parameterization / nuclei.
  // Return geometry covariance for TOF, physical momentum + covariance for
  // tracking.
  void SetStateConverter(StateConverter);
  bool Open(std::string &);
  ReadStatus Next(AnalysisEvent &, std::string &);
  const BenoInputStatistics &GetStatistics() const;

private:
  CeeKFParticleBenoSource(const CeeKFParticleBenoSource &) = delete;
  CeeKFParticleBenoSource &operator=(const CeeKFParticleBenoSource &) = delete;
  struct Implementation;
  Implementation *fImplementation;
};
} // namespace CeeV0Analysis
#endif
