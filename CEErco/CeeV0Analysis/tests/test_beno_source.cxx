#include "CeeBenoEtofHit.h"
#include "CeeBenoEvent.h"
#include "CeeBenoItofHit.h"
#include "CeeBenoTpcEtofMatchInfo.h"
#include "CeeBenoTpcItofMatchInfo.h"
#include "CeeBenoTpcTrack.h"
#include "CeeBenoTpcTrackCovMatrix.h"
#include "CeeKFParticleInput.h"
#include "MakerTestFixtures.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <unistd.h>
AnalysisEvent BenoFixture() {
  auto e = MakeEvent(false);
  auto rotate = [](const Vec3 &v) {
    const double a = .3;
    return Vec3(std::cos(a) * v[0] - std::sin(a) * v[1],
                std::sin(a) * v[0] + std::cos(a) * v[1], v[2]);
  };
  for (auto &t : e.tracks) {
    auto &x = t.input;
    x.positionCm = rotate(x.positionCm) + Vec3(2., 1., 0.);
    x.direction = rotate(x.direction);
    x.trackingMomentumGeV = rotate(x.trackingMomentumGeV);
    x.tofPositionCm = rotate(x.tofPositionCm) + Vec3(2., 1., 0.);
  }
  return e;
}
void WriteBeno(const std::string &path, int variant = 0,
               bool tracking = false) {
  TFile file(path.c_str(), "RECREATE");
  TTree tree("ceesim", "synthetic data, actual Beno EDM");
  TClonesArray *events = new TClonesArray("CeeBenoEvent"),
               *tracks = new TClonesArray("CeeBenoTpcTrack"),
               *covs = new TClonesArray("CeeBenoTpcTrackCovMatrix");
  TClonesArray *matches = new TClonesArray("CeeBenoTpcItofMatchInfo"),
               *hits = new TClonesArray("CeeBenoItofHit");
  tree.Branch("Event", &events);
  tree.Branch("TpcTrack", &tracks);
  if (variant != 9)
    tree.Branch("TpcTrackCovMatrix", &covs);
  tree.Branch("TpcItofMatchInfo", &matches);
  tree.Branch("ItofHit", &hits);
  TClonesArray *em = new TClonesArray("CeeBenoTpcEtofMatchInfo"),
               *eh = new TClonesArray("CeeBenoEtofHit");
  if (variant == 7 || variant == 8) {
    tree.Branch("TpcEtofMatchInfo", &em);
    tree.Branch("EtofHit", &eh);
  }
  // Empty upstream event must be skipped, with eventId retaining chain entry.
  tree.Fill();
  auto *e = new ((*events)[0]) CeeBenoEvent();
  e->setPrimaryVertexIndex(0);
  e->setPrimaryVertex(2., 1., 0.);
  e->setPrimaryVertexError(0.01, 0.01, 0.01);
  const auto fixture = BenoFixture();
  for (int i = 0; i < 3; ++i)
    new ((*hits)[i])
        CeeBenoItofHit(); // raw coordinates invalid: use MatchInfo instead
  for (int i = 0; i < 2; ++i) {
    const auto &x = fixture.tracks[i].input;
    Vec3 u;
    assert(Normalize(x.direction, u));
    double step = -(x.positionCm[0] * u[0] + x.positionCm[1] * u[1]) /
                  (u[0] * u[0] + u[1] * u[1]);
    Vec3 r = x.positionCm + step * u;
    double phi = std::atan2(u[1], u[0]), theta = std::acos(u[2]);
    double d = r[0] * std::sin(phi) - r[1] * std::cos(phi);
    auto *t = new ((*tracks)[i]) CeeBenoTpcTrack();
    t->setTrackId(100 + i * 7);
    t->setNHitsFit(i == 0 ? 30 : -30);
    t->setOrigin(r[0], r[1], r[2]);
    double p = tracking ? Norm(x.trackingMomentumGeV) : 100.;
    t->setGlobalMomentum(-p * u[0], -p * u[1],
                         -p * u[2]); // upstream GetEntryMomentum convention
    t->setCovMatrixIndex(1 - i);
    t->setItofBestMatchIndex(i + 1);
    t->setNSigmaProton(i == 0 ? 0. : 10.);
    t->setNSigmaPion(i == 1 ? 0. : 10.);
    auto *c = new ((*covs)[1 - i]) CeeBenoTpcTrackCovMatrix();
    c->setTpcTrackIndex(i);
    c->setParams(d, r[2], phi, theta, tracking ? (i == 0 ? 1. : -1.) / p : 0.);
    c->setSigmas({{.02f, .02f, .0005f, .0005f, tracking ? .002f : 0.f}});
    std::array<float, 10> corr{};
    corr[1] = 0.2;
    corr[4] = -0.15;
    c->setCorrelations(corr);
    if (variant == 5 && i == 0)
      c->setTpcTrackIndex(1);
    if (variant == 6 && i == 0)
      c->setDcaXY(20.); // mislabeled parameter/reference check
    auto *m = new ((*matches)[1 - i]) CeeBenoTpcItofMatchInfo();
    m->setTpcTrackIndex(i);
    m->setItofHitIndex(i + 1);
    m->setTpcItofMatchFlag(2);
    m->setTpcItofMatchChi2(1.);
    m->setItof(x.tofTimeNs);
    m->setItofHitPositionXYZ(x.tofPositionCm[0], x.tofPositionCm[1],
                             x.tofPositionCm[2]);
    if (variant == 1 && i == 0)
      m->setHitShared(true);
    if (variant == 2 && i == 0)
      m->setItof(-999.);
    if (variant == 3 && i == 0)
      m->setItofHitIndex(50);
    if (variant == 7 || variant == 8) {
      new ((*eh)[i]) CeeBenoEtofHit();
      t->setEtofBestMatchIndex(i);
      auto *me = new ((*em)[i]) CeeBenoTpcEtofMatchInfo();
      me->setTpcTrackIndex(i);
      me->setEtofHitIndex(i);
      me->setTpcEtofMatchFlag(1);
      me->setTpcEtofMatchChi2(2.);
      me->setEtof(x.tofTimeNs + .5);
      me->setEtofHitPositionXYZ(x.tofPositionCm[0], x.tofPositionCm[1],
                                x.tofPositionCm[2]);
      if (variant == 8)
        m->setItof(-999.);
    }
    if (variant == 4 && i == 0)
      t->setItofBestMatchIndex(0); // do not substitute a different match
  }
  tree.Fill();
  tree.Write();
  delete events;
  delete tracks;
  delete covs;
  delete matches;
  delete hits;
  delete em;
  delete eh;
}
int main() {
  const std::string stem = "/tmp/cee_beno_" + std::to_string(getpid());
  for (int i = 0; i <= 9; ++i)
    WriteBeno(stem + std::to_string(i) + ".root", i);
  WriteBeno(stem + "tracking.root", 0, true);
  AnalysisConfig c;
  ConfigureDataMode(c, kRealData);
  c.finder = Configuration(false);
  c.finder.propagateSharedTimeUncertainty = true;
  c.finder.additionalDecayTimeSigmaNs = 0.;
  c.tofTimeSigmaNs = .02;
  c.tofPathLengthSigmaCm = .05;
  std::string error;
  assert(ConfigureInputModes(c, "benodst", "tof", 0., error));
  c.inputPath = stem + "0.root";
  AnalysisEvent event;
  double mass = 0., variance = 0.;
  {
    CeeKFParticleBenoSource s(c);
    assert(s.Open(error));
    assert(s.Next(event, error) == CeeKFParticleEventSource::kEvent);
    assert(event.eventId == 1 && event.tracks.size() == 2 &&
           event.input.startTimeNs == 0. && event.input.sigmaStartTimeNs == .1);
    assert(event.tracks[0].input.id == 100 && event.tracks[1].input.id == 107);
    auto &x = event.tracks[0].input;
    // Built-in conversion is invariant under the arbitrary zero-field |p|
    // placeholder.assert(x.hasGeometryCovariance&&!x.hasTrackingMomentum&&!x.chargeKnown);
    assert(std::isnan(event.tracks[0].chi2Ndf) &&
           !event.tracks[0].fitConverged);
    auto original = BenoFixture();
    assert(std::fabs(x.tofTimeNs - original.tracks[0].input.tofTimeNs) < 1.e-5);
    assert(Norm(x.tofPositionCm - original.tracks[0].input.tofPositionCm) <
           .02);
    assert(std::fabs(x.geometryCovariance.m[1][3]) > 1.e-10);
    // Direction reversal must reverse the position/momentum covariance block.
    Mat6 forward = BuildTofStateCovariance(x, x.direction, 2., .1),
         reverse = BuildTofStateCovariance(x, -1. * x.direction, 2., .1);
    assert(std::fabs(forward.m[1][3] + reverse.m[1][3]) < 1.e-12);
    CeeKFParticleAnalysisMaker maker(c.finder);
    maker.SetPIDSelector(PID);
    assert(maker.Make(event, error));
    if (maker.GetCandidates().empty())
      throw std::runtime_error("Beno TOF produced no Lambda");
    mass = maker.GetCandidates()[0].massGeV;
    variance = maker.GetCandidates()[0].sigmaMassGeV;
    assert(std::fabs(mass - 1.115683) < .001);
    assert(s.Next(event, error) == CeeKFParticleEventSource::kEnd);
  }
  // Shared timing uncertainty changes errors, not central values.
  {
    auto z = c;
    z.t0.sigmaNs = 0.;
    CeeKFParticleBenoSource s(z);
    assert(s.Next(event, error) == CeeKFParticleEventSource::kEvent);
    CeeKFParticleAnalysisMaker maker(z.finder);
    maker.SetPIDSelector(PID);
    assert(maker.Make(event, error));
    assert(std::fabs(maker.GetCandidates()[0].massGeV - mass) < 1.e-12);
    assert(maker.GetCandidates()[0].sigmaMassGeV < variance);
  }
  for (int v : {1, 2, 4, 5, 6}) {
    c.inputPath = stem + std::to_string(v) + ".root";
    CeeKFParticleBenoSource s(c);
    assert(s.Next(event, error) == CeeKFParticleEventSource::kEvent &&
           event.tracks.size() == 1);
  }
  c.inputPath = stem + "3.root";
  {
    CeeKFParticleBenoSource s(c);
    assert(s.Next(event, error) == CeeKFParticleEventSource::kError &&
           event.tracks.empty());
  }
  c.inputPath = stem + "0.root";
  {
    BenoInputOptions o;
    o.requireUniqueMatch = true;
    CeeKFParticleBenoSource s(c, o);
    assert(s.Next(event, error) == CeeKFParticleEventSource::kEvent &&
           event.tracks.empty());
  }
  {
    CeeKFParticleBenoSource s(c);
    s.SetTrackIndexSelector(
        [](long long, const TClonesArray &) { return std::vector<int>{1}; });
    assert(s.Next(event, error) == CeeKFParticleEventSource::kEvent &&
           event.tracks.size() == 1 && event.tracks[0].input.id == 107);
  }
  {
    CeeKFParticleBenoSource s(c);
    s.SetTrackIndexSelector(
        [](long long, const TClonesArray &) { return std::vector<int>{0, 0}; });
    assert(s.Next(event, error) == CeeKFParticleEventSource::kError);
  }
  {
    auto bad = c;
    bad.tofOffsetNs = 1.;
    CeeKFParticleBenoSource s(bad);
    assert(!s.Open(error));
  }
  {
    auto bad = c;
    bad.t0.source = kT0Hit;
    CeeKFParticleBenoSource s(bad);
    assert(!s.Open(error));
  }
  {
    auto bad = c;
    bad.maximumTracksPerEvent = 1;
    CeeKFParticleBenoSource s(bad);
    assert(s.Next(event, error) == CeeKFParticleEventSource::kError);
  }
  for (int v : {7, 8}) {
    c.inputPath = stem + std::to_string(v) + ".root";
    for (bool prefer : {false, true}) {
      BenoInputOptions o;
      o.preferItof = prefer;
      CeeKFParticleBenoSource source(c, o);
      assert(source.Next(event, error) == CeeKFParticleEventSource::kEvent &&
             event.tracks.size() == 2);
      double expected = BenoFixture().tracks[0].input.tofTimeNs +
                        ((v == 8 || !prefer) ? .5 : 0.);
      assert(std::fabs(event.tracks[0].input.tofTimeNs - expected) < 1.e-5);
    }
  }
  assert(!ConfigureInputModes(c, "benodst", "tof", .5, error));
  assert(!ConfigureInputModes(c, "unknown", "tof", 0., error));
  assert(!ConfigureInputModes(c, "raw", "bad", 0., error));
  {
    std::ofstream f(stem + ".list");
    f << stem << "0.root\n" << stem << "0.root\n";
  }
  c.inputPath = stem + ".list";
  {
    CeeKFParticleBenoSource s(c);
    assert(s.Next(event, error) == CeeKFParticleEventSource::kEvent &&
           event.eventId == 1);
    assert(s.Next(event, error) == CeeKFParticleEventSource::kEvent &&
           event.eventId == 3);
  }
  {
    std::ofstream f(stem + ".list");
    f << stem << "0.root\n" << stem << "9.root\n";
  }
  {
    CeeKFParticleBenoSource s(c);
    assert(!s.Open(error));
  }
  c.inputPath = stem + "tracking.root";
  assert(ConfigureInputModes(c, "benodst", "tracking", .5, error));
  {
    CeeKFParticleBenoSource s(c);
    assert(s.Next(event, error) == CeeKFParticleEventSource::kEvent &&
           event.tracks.size() == 2);
    assert(event.tracks[0].input.hasTrackingCovariance &&
           !event.tracks[0].input.hasTof);
    const auto &tr = event.tracks[0].input;
    const double phi = std::atan2(tr.direction[1], tr.direction[0]);
    const Vec3 normal(std::cos(phi), std::sin(phi), 0.);
    double normalVariance = 0.;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        normalVariance += normal[i] * tr.trackingCovariance.m[i][j] * normal[j];
    assert(std::fabs(normalVariance) < 1.e-10); // fixed local reference plane

    CeeKFParticleAnalysisMaker maker(c.finder);
    maker.SetPIDSelector(PID);
    assert(maker.Make(event, error));
    assert(maker.GetCandidates().size() == 1 &&
           std::fabs(maker.GetCandidates()[0].massGeV - 1.115683) < .001);
  }
  std::cout << "Beno actual-EDM ROOT roundtrip, zero-field Lambda mass=" << mass
            << ", sigmaM=" << variance
            << "; shared-T0, covariance/sign, IDs, subset, list, matching and "
               "failure checks passed\n";
  if (std::getenv("CEE_BENO_KEEP_FIXTURES"))
    return 0;
  for (int i = 0; i <= 9; ++i)
    std::remove((stem + std::to_string(i) + ".root").c_str());
  std::remove((stem + "tracking.root").c_str());
  std::remove((stem + ".list").c_str());
}
