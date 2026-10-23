#ifndef CEE_KFP_ROOT_SOURCE_H
#define CEE_KFP_ROOT_SOURCE_H
#include "CeeKFParticleEvent.h"
#include "CeeV0Analyzer.h"
#include <functional>
class TClonesArray;
namespace CeeV0Analysis {
// The only new component coupled to CeeTpcTrack / TOF branch classes.
class CeeKFParticleRootSource : public CeeKFParticleEventSource {
 public:
  typedef std::function<std::vector<int>(long long,const TClonesArray&)> TrackIndexSelector;
  explicit CeeKFParticleRootSource(const AnalysisConfig& config);
  ~CeeKFParticleRootSource();
  // Optional original array indices BEFORE adaptation and TOF assignment.
  // Returned indices must be unique/in range; empty means select no tracks.
  void SetTrackIndexSelector(TrackIndexSelector selector);
  bool Open(std::string& error);
  ReadStatus Next(AnalysisEvent& event,std::string& error);
 private:
  CeeKFParticleRootSource(const CeeKFParticleRootSource&);
  CeeKFParticleRootSource& operator=(const CeeKFParticleRootSource&);
  struct Implementation;
  Implementation* fImplementation;
};
}
#endif
