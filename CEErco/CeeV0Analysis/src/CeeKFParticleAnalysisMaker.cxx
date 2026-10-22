#include "CeeKFParticleAnalysisMaker.h"
#include <cmath>
#include <exception>
#include <set>
namespace CeeV0Analysis {
CeeKFParticleAnalysisMaker::CeeKFParticleAnalysisMaker(const CeeV0Lite::FinderConfig& c)
  : fInterface(c), fEventAccepted(false) {}
bool CeeKFParticleAnalysisMaker::AcceptEvent(const AnalysisEvent& e) const {
  return !fEventSelector || fEventSelector(e);
}
bool CeeKFParticleAnalysisMaker::AcceptTrack(const AnalysisTrack& t,const AnalysisEvent& e) const {
  if (fTrackCuts.minimumHits>0 && t.nHits<fTrackCuts.minimumHits) return false;
  if (std::isfinite(fTrackCuts.maximumChi2Ndf) &&
      (!std::isfinite(t.chi2Ndf) || t.chi2Ndf<0. || t.chi2Ndf>fTrackCuts.maximumChi2Ndf)) return false;
  if (fTrackCuts.requireConverged && !t.fitConverged) return false;
  return !fTrackSelector || fTrackSelector(t,e);
}
unsigned CeeKFParticleAnalysisMaker::IdentifyTrack(const AnalysisTrack& t,const AnalysisEvent& e) const {
  return fPIDSelector ? fPIDSelector(t,e) : t.input.hypothesisMask;
}
bool CeeKFParticleAnalysisMaker::AcceptMother(const CeeV0Lite::Candidate& c,const AnalysisEvent& e) const {
  return !fMotherSelector || fMotherSelector(c,e);
}
bool CeeKFParticleAnalysisMaker::Make(const AnalysisEvent& e,std::string& error) {
  error.clear(); fEventAccepted=false; fSelectedTracks.clear(); fCandidates.clear(); fInterface.Clear();
  ++fStatistics.eventsRead;
  if (fTrackCuts.minimumHits<0 || std::isnan(fTrackCuts.maximumChi2Ndf) || fTrackCuts.maximumChi2Ndf<0.) {
    error="invalid Maker track cuts"; return false;
  }
  if (!AcceptEvent(e)) return true;
  std::vector<std::size_t> indices;
  if (fTrackListSelector) indices=fTrackListSelector(e);
  else { indices.reserve(e.tracks.size()); for (std::size_t i=0;i<e.tracks.size();++i) indices.push_back(i); }
  std::set<std::size_t> seen;
  for (std::size_t j=0;j<indices.size();++j) {
    const std::size_t i=indices[j];
    if (i>=e.tracks.size() || !seen.insert(i).second) {
      fSelectedTracks.clear(); error="track-list selector returned an out-of-range or duplicate index"; return false;
    }
    const AnalysisTrack& t=e.tracks[i];
    if (!AcceptTrack(t,e)) continue;
    const unsigned pid=IdentifyTrack(t,e);
    if (!pid) continue;
    CeeV0Lite::TrackInput selected=t.input; selected.hypothesisMask=pid;
    fSelectedTracks.push_back(selected);
  }
  if (!fInterface.Reconstruct(e.input,fSelectedTracks,error)) { fSelectedTracks.clear(); return false; }
  const std::vector<CeeV0Lite::Candidate>& all=fInterface.GetCandidates();
  for (std::size_t i=0;i<all.size();++i) if (AcceptMother(all[i],e)) fCandidates.push_back(all[i]);
  fEventAccepted=true; ++fStatistics.eventsAccepted;
  fStatistics.tracksSupplied+=e.tracks.size(); fStatistics.tracksSelected+=fSelectedTracks.size();
  fStatistics.mothersReconstructed+=all.size(); fStatistics.mothersAccepted+=fCandidates.size();
  return true;
}
bool CeeKFParticleAnalysisMaker::Run(CeeKFParticleEventSource& source,CeeKFParticleCandidateSink& sink,
    std::string& error,long long maximumEvents) {
  error.clear();
  if (maximumEvents < -1) { error="maximumEvents must be -1 or nonnegative"; return false; }
  if (!sink.Begin(error)) return false;
  bool ok=true;
  try {
    for (long long i=0;maximumEvents<0 || i<maximumEvents;++i) {
      AnalysisEvent e;
      const CeeKFParticleEventSource::ReadStatus status=source.Next(e,error);
      if (status==CeeKFParticleEventSource::kEnd) break;
      if (status==CeeKFParticleEventSource::kError || !Make(e,error)) { ok=false; break; }
      if (fEventAccepted && !sink.Write(e,fCandidates,error)) { ok=false; break; }
    }
  } catch (const std::exception& ex) { error=ex.what(); ok=false; }
  // Close even on read/selection/write errors; preserve the original error.
  std::string finishError;
  const bool finished=sink.Finish(finishError);
  if (ok && !finished) error=finishError;
  return ok && finished;
}
}
