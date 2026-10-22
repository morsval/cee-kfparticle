#include "MakerTestFixtures.h"
#include <algorithm>
AnalysisTrack makeTrack(int id, int q, double mass, const Vec3 &v,
                        const Vec3 &p) {
  AnalysisTrack a;
  TrackState t, o;
  t.id = id;
  t.charge = q;
  t.massGeV = mass;
  t.positionCm = v;
  t.momentumGeV = p;
  for (int i = 0; i < 6; ++i)
    t.covariance.m[i][i] = i < 3 ? 1e-4 : 1e-6;
  assert(PropagateState(t, 10., .5, TransportConfig(), o));
  auto &x = a.input;
  x.id = id;
  x.charge = q;
  x.chargeKnown = true;
  x.hasTrackingMomentum = x.hasTrackingCovariance = true;
  x.positionCm = o.positionCm;
  x.trackingMomentumGeV = o.momentumGeV;
  x.trackingCovariance = o.covariance;
  Normalize(o.momentumGeV, x.direction);
  return a;
}
int main() {
  CeeKFParticleAnalysisMaker maker(Configuration(true));
  auto &engine = maker.Channels();
  DecayChannel l;
  l.name = "lambda";
  l.motherPdg = 3122;
  l.daughters = {DaughterSpec::Track("p", {2212, 1, .9382720813}),
                 DaughterSpec::Track("pi", {-211, -1, .13957039})};
  engine.pid = [](const AnalysisTrack &t, const ParticleHypothesis &p,
                  const AnalysisEvent &) {
    return (t.input.id == 100 && p.pdg == 2212) ||
           (t.input.id == 101 && p.pdg == -211) ||
           (t.input.id == 102 && p.pdg == -321);
  };
  engine.channels = {l};
  auto e = MakeEvent(true);
  std::string error;
  assert(maker.MakeChannels(e, error));
  assert(maker.GetChannelCandidates().size() == 1);
  auto first = maker.GetChannelCandidates()[0];
  assert(std::fabs(first.massGeV - 1.115683) < 1e-6);
  assert(Norm(first.mother.positionCm - Vec3(0, 0, 5)) < 1e-5);
  assert(first.ndf == 1 && first.sigmaMassGeV > 0);
  assert(first.Daughter("p").trackId == 100);
  std::reverse(e.tracks.begin(), e.tracks.end());
  assert(maker.MakeChannels(e, error));
  assert(std::fabs(maker.GetChannelCandidates()[0].massGeV - first.massGeV) <
         1e-10);
  // Lambda decays at z=5, Omega at z=2; the composite Lambda must be
  // transported backwards without treating its decay vertex as Omega's SV.
  e.tracks.push_back(
      makeTrack(102, -1, .493677, Vec3(0, 0, 2), Vec3(.3, .1, .7)));
  DecayChannel omega;
  omega.name = "omega";
  omega.motherPdg = 3334;
  omega.daughters = {DaughterSpec::Composite("lambda", "lambda"),
                     DaughterSpec::Track("kaon", {-321, -1, .493677})};
  omega.fitCollisionPv = true;
  engine.channels = {omega, l};
  assert(maker.MakeChannels(e, error));
  assert(maker.GetChannelCandidates().size() == 2);
  auto o = maker.GetChannelCandidates()[1];
  assert(o.channel == "omega" && o.daughters[0].candidateId == 0 &&
         o.leafTrackIds.size() == 3);
  assert(Norm(o.mother.positionCm - Vec3(0, 0, 2)) < 1e-4);
  assert(o.mother.covarianceValid);
  assert(o.Daughter("kaon").trackId == 102);
  // True three-body common vertex, not a cascade.
  DecayChannel three;
  three.name = "three";
  three.motherPdg = 1010010030;
  three.daughters = {DaughterSpec::Track("d", {1000010020, 1, 1.8756129}),
                     DaughterSpec::Track("p", {2212, 1, .9382720813}),
                     DaughterSpec::Track("pi", {-211, -1, .13957039})};
  e.tracks = {makeTrack(1, 1, 1.8756129, Vec3(1, 2, 5), Vec3(.2, .1, .8)),
              makeTrack(2, 1, .9382720813, Vec3(1, 2, 5), Vec3(-.2, .2, .5)),
              makeTrack(3, -1, .13957039, Vec3(1, 2, 5), Vec3(.1, -.2, .3))};
  engine.channels = {three};
  engine.pid = [](const AnalysisTrack &t, const ParticleHypothesis &p,
                  const AnalysisEvent &) {
    return (t.input.id == 1 && p.pdg == 1000010020) ||
           (t.input.id == 2 && p.pdg == 2212) ||
           (t.input.id == 3 && p.pdg == -211);
  };
  assert(maker.MakeChannels(e, error));
  assert(maker.GetChannelCandidates().size() == 1);
  auto a = maker.GetChannelCandidates()[0];
  assert(a.ndf == 3 && a.daughters.size() == 3);
  assert(Norm(a.mother.positionCm - Vec3(1, 2, 5)) < 1e-4);
  assert(IsPositiveSemidefinite(a.mother.covariance));
  std::reverse(engine.channels[0].daughters.begin(),
               engine.channels[0].daughters.end());
  assert(maker.MakeChannels(e, error));
  assert(std::fabs(maker.GetChannelCandidates()[0].massGeV - a.massGeV) < 1e-9);
  // Leaf reuse must reject a bachelor identical to a Lambda daughter.
  engine.channels = {l, omega};
  engine.channels[1].daughters[1] =
      DaughterSpec::Track("pi", {-211, -1, .13957039});
  e = MakeEvent(true);
  engine.pid = [](const AnalysisTrack &t, const ParticleHypothesis &p,
                  const AnalysisEvent &) {
    return (t.input.id == 100 && p.pdg == 2212) ||
           (t.input.id == 101 && p.pdg == -211);
  };
  assert(maker.MakeChannels(e, error));
  assert(maker.GetChannelCandidates().size() == 1);
  engine.channels[0].daughters[0] = DaughterSpec::Composite("cycle", "omega");
  assert(!maker.MakeChannels(e, error));
  assert(error.find("cyclic") != std::string::npos);
  assert(maker.GetChannelCandidates().empty());
  engine.channels = {l};
  engine.maximumCombinations = 0;
  assert(!maker.MakeChannels(e, error));
  assert(maker.GetChannelCandidates().empty());
  engine.maximumCombinations = 100;
  engine.config.momentumMode = kTofMomentum;
  assert(!maker.MakeChannels(e, error));
  assert(error.find("tracking momentum") != std::string::npos);
  std::cout << "generic channels: two/three-body, cascade, covariance, order, "
               "leaf reuse, errors passed\n";
}
