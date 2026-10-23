#ifndef CEE_KFP_CHANNELS_H
#define CEE_KFP_CHANNELS_H
#include "CeeKFParticleEvent.h"
#include "CeeV0LiteFinder.h"
#include <functional>
namespace CeeV0Analysis {
// Charge is in units of e, momentum is PHYSICAL GeV/c (not rigidity).
struct ParticleHypothesis {
  int pdg, charge;
  double massGeV;
  ParticleHypothesis(int p = 0, int q = 0, double m = 0)
      : pdg(p), charge(q), massGeV(m) {}
};
struct DaughterSpec {
  std::string role, sourceChannel;
  ParticleHypothesis particle;
  // Empty sourceChannel: a selected input track. Otherwise a composite.
  static DaughterSpec Track(const std::string &role,
                            const ParticleHypothesis &p);
  static DaughterSpec Composite(const std::string &role,
                                const std::string &channel);
};
struct ChannelDaughter {
  std::string role;
  int pdg = 0, charge = 0, trackId = -1, candidateId = -1;
  CeeV0Lite::MotherState inputState, fittedState;
  double dcaToPvCm = std::numeric_limits<double>::quiet_NaN();
};
struct ChannelCandidate {
  int id = -1, pdg = 0, charge = 0, ndf = 0;
  std::string channel;
  CeeV0Lite::MotherState mother; // unconstrained decay state, used in cascades
  CeeV0Lite::ProductionState motherTopo; // separate optional collision-PV copy
  CeeV0Lite::FlightObservables flight;
  std::vector<ChannelDaughter> daughters;
  std::vector<int> leafTrackIds;
  double massGeV = 0, sigmaMassGeV = 0, chi2 = 0,
         maximumDaughterSeparationCm = 0;
  const ChannelDaughter &Daughter(const std::string &role) const;
};
struct DecayChannel {
  std::string name;
  int motherPdg = 0;
  std::vector<DaughterSpec> daughters;
  bool fitCollisionPv = false;
  // Applied before this candidate is made available to downstream channels.
  std::function<bool(const ChannelCandidate &, const AnalysisEvent &)> select;
};
class CeeKFParticleChannelSink {
public:
  virtual ~CeeKFParticleChannelSink() {}
  virtual bool Begin(std::string &) = 0;
  virtual bool Write(const AnalysisEvent &,
                     const std::vector<ChannelCandidate> &, std::string &) = 0;
  virtual bool Finish(std::string &) = 0;
};
// Explicit opt-in tracking-momentum engine. The legacy TOF engine is retained.
// User PID callback is invoked BEFORE combinations; no implicit p/pi PID.
class CeeKFParticleChannels {
public:
  typedef std::function<bool(const AnalysisTrack &, const ParticleHypothesis &,
                             const AnalysisEvent &)>
      PIDSelector;
  explicit CeeKFParticleChannels(
      const CeeV0Lite::FinderConfig &c = CeeV0Lite::FinderConfig())
      : config(c) {}
  CeeV0Lite::FinderConfig config;
  std::vector<DecayChannel> channels;
  PIDSelector pid;
  std::size_t maximumCombinations = 1000000;
  bool Reconstruct(const AnalysisEvent &,
                   const std::vector<std::size_t> &selected, std::string &);
  const std::vector<ChannelCandidate> &GetCandidates() const {
    return candidates;
  }
  void Clear() { candidates.clear(); }

private:
  std::vector<ChannelCandidate> candidates;
};
} // namespace CeeV0Analysis
#endif
