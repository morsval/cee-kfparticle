#ifndef CEE_KFP_EVENT_H
#define CEE_KFP_EVENT_H
#include "CeeV0LiteTypes.h"
#include <limits>
#include <string>
#include <vector>

namespace CeeV0Analysis {
// Format-independent value objects: never retain pointers into an input TTree.
struct AnalysisTrack {
  AnalysisTrack() : nHits(-1), fitConverged(false),
    chi2Ndf(std::numeric_limits<double>::quiet_NaN()),
    dEdx(std::numeric_limits<double>::quiet_NaN()),
    nSigmaProton(std::numeric_limits<double>::quiet_NaN()),
    nSigmaPion(std::numeric_limits<double>::quiet_NaN()),
    nSigmaKaon(std::numeric_limits<double>::quiet_NaN()),
    tofMass2GeV2(std::numeric_limits<double>::quiet_NaN()) {}
  CeeV0Lite::TrackInput input;
  int nHits;
  bool fitConverged;
  double chi2Ndf, dEdx;
  // Optional calibrated PID. Missing values are NaN, never assumed to be 0.
  double nSigmaProton, nSigmaPion, nSigmaKaon, tofMass2GeV2;
};
struct AnalysisEvent {
  AnalysisEvent() : eventId(0), inputTrackCount(0) {}
  long long eventId; // source entry number, unique across one input chain
  CeeV0Lite::EventInput input;
  std::size_t inputTrackCount;
  std::vector<AnalysisTrack> tracks;
};
class CeeKFParticleEventSource {
 public:
  enum ReadStatus { kEvent, kEnd, kError };
  virtual ~CeeKFParticleEventSource() {}
  virtual ReadStatus Next(AnalysisEvent& event, std::string& error) = 0;
};
// Both default ROOT output and completely user-owned TTree implementations
// use this interface. Called once for each accepted event, including empty ones.
class CeeKFParticleCandidateSink {
 public:
  virtual ~CeeKFParticleCandidateSink() {}
  virtual bool Begin(std::string& error) = 0;
  virtual bool Write(const AnalysisEvent& event,
                     const std::vector<CeeV0Lite::Candidate>& mothers,
                     std::string& error) = 0;
  virtual bool Finish(std::string& error) = 0;
};
} // namespace CeeV0Analysis
#endif
