#ifndef CEE_KFP_INTERFACE_H
#define CEE_KFP_INTERFACE_H
#include "CeeV0LiteFinder.h"
#include <string>
namespace CeeV0Analysis {
class CeeKFParticleInterface {
 public:
  explicit CeeKFParticleInterface(const CeeV0Lite::FinderConfig& config = CeeV0Lite::FinderConfig());
  void SetConfig(const CeeV0Lite::FinderConfig& config) { fFinder.SetConfig(config); }
  const CeeV0Lite::FinderConfig& GetConfig() const { return fFinder.GetConfig(); }
  // Receives ONLY user-selected tracks with PID hypothesis masks already set.
  // IDs must be unique and nonnegative within an event. No MC truth is used.
  bool Reconstruct(const CeeV0Lite::EventInput& event,
                   const std::vector<CeeV0Lite::TrackInput>& tracks,
                   std::string& error);
  void Clear();
  const std::vector<CeeV0Lite::Candidate>& GetCandidates() const { return fCandidates; }
  const std::vector<CeeV0Lite::TrackReferenceDiagnostic>& GetReferenceDiagnostics() const { return fDiagnostics; }
 private:
  CeeV0Lite::Finder fFinder;
  std::vector<CeeV0Lite::Candidate> fCandidates;
  std::vector<CeeV0Lite::TrackReferenceDiagnostic> fDiagnostics;
};
}
#endif
