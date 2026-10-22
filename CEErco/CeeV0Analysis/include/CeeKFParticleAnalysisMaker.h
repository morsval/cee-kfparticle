#ifndef CEE_KFP_ANALYSIS_MAKER_H
#define CEE_KFP_ANALYSIS_MAKER_H
#include "CeeKFParticleEvent.h"
#include "CeeKFParticleInterface.h"
#include <functional>
namespace CeeV0Analysis {
struct MakerTrackCuts {
  MakerTrackCuts() : minimumHits(0), maximumChi2Ndf(std::numeric_limits<double>::infinity()),
    requireConverged(false) {}
  int minimumHits;
  double maximumChi2Ndf;
  bool requireConverged;
};
struct MakerStatistics {
  MakerStatistics() : eventsRead(0), eventsAccepted(0), tracksSupplied(0),
    tracksSelected(0), mothersReconstructed(0), mothersAccepted(0) {}
  unsigned long long eventsRead, eventsAccepted, tracksSupplied, tracksSelected,
    mothersReconstructed, mothersAccepted;
};
class CeeKFParticleAnalysisMaker {
 public:
  typedef std::function<bool(const AnalysisEvent&)> EventSelector;
  typedef std::function<bool(const AnalysisTrack&,const AnalysisEvent&)> TrackSelector;
  typedef std::function<unsigned(const AnalysisTrack&,const AnalysisEvent&)> PIDSelector;
  typedef std::function<bool(const CeeV0Lite::Candidate&,const AnalysisEvent&)> MotherSelector;
  // Optional caller-provided indices into event.tracks; avoids scanning the
  // entire event in Make(). Empty list deliberately means no selected tracks.
  typedef std::function<std::vector<std::size_t>(const AnalysisEvent&)> TrackListSelector;
  explicit CeeKFParticleAnalysisMaker(const CeeV0Lite::FinderConfig& config = CeeV0Lite::FinderConfig());
  virtual ~CeeKFParticleAnalysisMaker() {}
  void SetTrackCuts(const MakerTrackCuts& cuts) { fTrackCuts=cuts; }
  void SetEventSelector(EventSelector f) { fEventSelector=f; }
  void SetTrackSelector(TrackSelector f) { fTrackSelector=f; }
  void SetPIDSelector(PIDSelector f) { fPIDSelector=f; }
  void SetMotherSelector(MotherSelector f) { fMotherSelector=f; }
  void SetTrackListSelector(TrackListSelector f) { fTrackListSelector=f; }
  CeeKFParticleInterface& GetInterface() { return fInterface; }
  const CeeKFParticleInterface& GetInterface() const { return fInterface; }
  // Baseline Finder cuts run inside reconstruction; AcceptMother runs AFTER
  // reconstruction and can inspect both mother (7D) and motherTopo (8D).
  virtual bool AcceptEvent(const AnalysisEvent&) const;
  virtual bool AcceptTrack(const AnalysisTrack&,const AnalysisEvent&) const;
  virtual unsigned IdentifyTrack(const AnalysisTrack&,const AnalysisEvent&) const;
  virtual bool AcceptMother(const CeeV0Lite::Candidate&,const AnalysisEvent&) const;
  bool Make(const AnalysisEvent& event,std::string& error);
  bool Run(CeeKFParticleEventSource& source,CeeKFParticleCandidateSink& sink,
           std::string& error,long long maximumEvents=-1);
  const std::vector<CeeV0Lite::Candidate>& GetCandidates() const { return fCandidates; }
  const std::vector<CeeV0Lite::TrackInput>& GetSelectedTracks() const { return fSelectedTracks; }
  const MakerStatistics& GetStatistics() const { return fStatistics; }
  bool EventAccepted() const { return fEventAccepted; }
  void ResetStatistics() { fStatistics=MakerStatistics(); }
 private:
  CeeKFParticleInterface fInterface;
  MakerTrackCuts fTrackCuts;
  MakerStatistics fStatistics;
  EventSelector fEventSelector;
  TrackSelector fTrackSelector;
  PIDSelector fPIDSelector;
  MotherSelector fMotherSelector;
  TrackListSelector fTrackListSelector;
  bool fEventAccepted;
  std::vector<CeeV0Lite::TrackInput> fSelectedTracks;
  std::vector<CeeV0Lite::Candidate> fCandidates;
};
}
#endif
