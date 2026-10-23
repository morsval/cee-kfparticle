#include "MakerTestFixtures.h"
#include <random>
int main() {
  CeeKFParticleAnalysisMaker m(Configuration(true));
  auto &en = m.Channels();
  DecayChannel c;
  c.name = "lambda";
  c.motherPdg = 3122;
  c.daughters = {DaughterSpec::Track("p", {2212, 1, .9382720813}),
                 DaughterSpec::Track("pi", {-211, -1, .13957039})};
  en.channels = {c};
  en.pid = [](const AnalysisTrack &t, const ParticleHypothesis &p,
              const AnalysisEvent &) {
    return (t.input.id == 100 && p.pdg == 2212) ||
           (t.input.id == 101 && p.pdg == -211);
  };
  std::mt19937 rng(64027);
  std::normal_distribution<double> normal;
  double sum = 0, sum2 = 0, chi2 = 0;
  int good = 0;
  std::string err;
  for (int n = 0; n < 240; ++n) {
    auto e = MakeEvent(true);
    for (auto &t : e.tracks) {
      TrackState ref, at, obs;
      ref.positionCm = t.input.positionCm;
      ref.momentumGeV = t.input.trackingMomentumGeV;
      ref.covariance = t.input.trackingCovariance;
      ref.charge = t.input.charge;
      assert(PropagateState(ref, -20., .5, TransportConfig(), at));
      at.covariance = Mat6();
      for (int i = 0; i < 3; ++i) {
        at.positionCm[i] += .02 * normal(rng);
        at.momentumGeV[i] += .001 * normal(rng);
        at.covariance.m[i][i] = .0004;
        at.covariance.m[i + 3][i + 3] = .000001;
      }
      assert(PropagateState(at, 20., .5, TransportConfig(), obs));
      t.input.positionCm = obs.positionCm;
      t.input.trackingMomentumGeV = obs.momentumGeV;
      t.input.trackingCovariance = obs.covariance;
      Normalize(obs.momentumGeV, t.input.direction);
    }
    assert(m.MakeChannels(e, err));
    if (m.GetChannelCandidates().empty())
      continue;
    auto &a = m.GetChannelCandidates()[0];
    double pull = (a.massGeV - 1.115683) / a.sigmaMassGeV;
    sum += pull;
    sum2 += pull * pull;
    chi2 += a.chi2;
    ++good;
  }
  double mean = sum / good, width = std::sqrt(sum2 / good - mean * mean);
  std::cout << "toy Lambda: accepted=" << good
            << "/240, mass pull mean=" << mean << ", width=" << width
            << ", mean chi2=" << chi2 / good << "\n";
  assert(good > 230 && std::fabs(mean) < .25 && width > .75 && width < 1.25 &&
         chi2 / good > .65 && chi2 / good < 1.35);
}
