#include "CeeKFParticleInterface.h"
#include <set>
namespace CeeV0Analysis {
CeeKFParticleInterface::CeeKFParticleInterface(const CeeV0Lite::FinderConfig& c) : fFinder(c) {}
void CeeKFParticleInterface::Clear() { fCandidates.clear(); fDiagnostics.clear(); }
bool CeeKFParticleInterface::Reconstruct(const CeeV0Lite::EventInput& event,
    const std::vector<CeeV0Lite::TrackInput>& tracks, std::string& error) {
  Clear(); error.clear();
  std::set<int> ids;
  const unsigned supported = CeeV0Lite::kPionHypothesis | CeeV0Lite::kProtonHypothesis;
  for (std::size_t i=0; i<tracks.size(); ++i) {
    if (tracks[i].id < 0 || !ids.insert(tracks[i].id).second) {
      error="selected track IDs must be nonnegative and unique"; return false;
    }
    if (!tracks[i].hypothesisMask || (tracks[i].hypothesisMask & ~supported)) {
      error="selected track has an empty or unsupported PID hypothesis mask"; return false;
    }
  }
  fFinder.Reconstruct(tracks,event,fCandidates,&fDiagnostics);
  return true;
}
}
