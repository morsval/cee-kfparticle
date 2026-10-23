#include "CeeKFParticleBenoSource.h"
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#ifdef CEE_WITH_BENODST
#include "CeeBenoEvent.h"
#include "CeeBenoTpcEtofMatchInfo.h"
#include "CeeBenoTpcItofMatchInfo.h"
#include "CeeBenoTpcTrack.h"
#include "CeeBenoTpcTrackCovMatrix.h"
#include "CeeKFParticleInputFiles.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#endif
namespace CeeV0Analysis {
using namespace CeeV0Lite;
#ifdef CEE_WITH_BENODST
namespace {
bool Finite6(const Mat6 &c) {
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      if (!std::isfinite(c.m[i][j]))
        return false;
  return true;
}
Vec3 Vector(const TVector3 &v) { return Vec3(v.X(), v.Y(), v.Z()); }
double Pid(double v) {
  return v == CeeBenoConst::kInvalidValue
             ? std::numeric_limits<double>::quiet_NaN()
             : v;
}
// Upstream CEETPC_Tracker/CeeTrackingTask::fillPar5AtVertex:
// d=(x*py-y*px)/pt, z0=z-pz*(x*px+y*py)/pt^2, phi,theta,q/p.
// At the transverse foot: r=(d sin(phi),-d cos(phi),z0).
bool PerigeeState(const CeeBenoTpcTrack &t, const CeeBenoTpcTrackCovMatrix &c,
                  bool tracking, const BenoInputOptions &options,
                  TrackInput &out) {
  if (!c.isValid())
    return false;
  const auto a = c.getParams(), sigma = c.getSigmas();
  const auto rho = c.getCorrelations();
  const double d = a[0], phi = a[2], theta = a[3];
  const double sp = std::sin(phi), cp = std::cos(phi), st = std::sin(theta),
               ct = std::cos(theta);
  if (!(st > 1.e-8))
    return false;
  Vec3 r(d * sp, -d * cp, a[1]);
  const Vec3 u(st * cp, st * sp, ct);
  Vec3 storedDirection;
  const Vec3 origin = Vector(t.origin()), p = Vector(t.globMom());
  if (!IsFinite(origin) || !Normalize(p, storedDirection))
    return false;
  const double transverse2 = storedDirection[0] * storedDirection[0] +
                             storedDirection[1] * storedDirection[1];
  if (!(transverse2 > 1.e-16))
    return false;
  const double step =
      -(origin[0] * storedDirection[0] + origin[1] * storedDirection[1]) /
      transverse2;
  const Vec3 foot = origin + step * storedDirection;
  // CeeTrackingTask::GetEntryMomentum can negate globMom; angles/qop and
  // VertexSeed retain the fitted orientation. Compare the same unoriented line.
  if (Norm(foot - r) > options.geometryPositionToleranceCm ||
      std::min(Norm(u - storedDirection), Norm(u + storedDirection)) >
          options.geometryDirectionToleranceRad)
    return false;
  double C[5][5] = {{0.}}, J[6][5] = {{0.}};
  int k = 0;
  for (int i = 0; i < 5; ++i) {
    C[i][i] = sigma[i] * sigma[i];
    for (int j = 0; j < i; ++j) {
      C[i][j] = C[j][i] = rho[k++] * sigma[i] * sigma[j];
    }
  }
  J[0][0] = sp;
  J[0][2] = d * cp;
  J[1][0] = -cp;
  J[1][2] = d * sp;
  J[2][1] = 1.;
  J[3][2] = -st * sp;
  J[3][3] = ct * cp;
  J[4][2] = st * cp;
  J[4][3] = ct * sp;
  J[5][3] = -st;
  double magnitude = 1.;
  if (tracking) {
    // Upstream computes C5 by perturbing a FIXED GenFit plane, whose axes are
    // dir.Cross(beamZ) and beamZ. Restore this plane's position fluctuations;
    // rotating the position onto a new perigee while keeping p unchanged would
    // introduce a spurious phase shift in a field. The saved origin gives the
    // plane's transverse axis (at d=0 the limiting correction vanishes).
    const double transverse = std::hypot(origin[0], origin[1]);
    double alpha = phi;
    if (transverse > 1.e-8) {
      alpha = std::atan2(origin[0], -origin[1]);
      if (std::cos(phi - alpha) < 0.)
        alpha += std::acos(-1.);
    }
    const double cosine = std::cos(phi - alpha),
                 tangent = std::tan(phi - alpha);
    if (!(cosine > 0.99))
      return false;
    const double h = d / cosine, sa = std::sin(alpha), ca = std::cos(alpha);
    r = Vec3(h * sa, -h * ca, a[1] - d * ct / st * tangent);
    if (Norm(r - origin) > options.geometryPositionToleranceCm)
      return false;
    J[0][0] = sa / cosine;
    J[1][0] = -ca / cosine;
    J[0][2] = h * tangent * sa;
    J[1][2] = -h * tangent * ca;
    J[2][0] = -ct / st * tangent;
    J[2][2] = -d * ct / st / (cosine * cosine);
    J[2][3] = d * tangent / (st * st);
    // Current Beno payload stores charge SIGN. Default converter supports
    // |q|=1.
    if (!std::isfinite(a[4]) || std::fabs(a[4]) < 1.e-12 || !(sigma[4] > 0.) ||
        (t.charge() != 1 && t.charge() != -1) || a[4] * t.charge() <= 0.)
      return false;
    magnitude = 1. / std::fabs(a[4]);
    if (std::fabs(magnitude - Norm(p)) > 0.05 * Norm(p) + 0.001)
      return false;
    for (int i = 3; i < 6; ++i) {
      for (int j = 0; j < 4; ++j)
        J[i][j] *= magnitude;
      J[i][4] = -magnitude / a[4] * u[i - 3];
    }
  }
  Mat6 covariance;
  const int dim =
      tracking ? 5 : 4; // q/p and all its cross terms are ignored at B=0/TOF.
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      for (int a = 0; a < dim; ++a)
        for (int b = 0; b < dim; ++b)
          covariance.m[i][j] += J[i][a] * C[a][b] * J[j][b];
  if (!Finite6(covariance) || !IsPositiveSemidefinite(covariance))
    return false;
  out.positionCm = r;
  out.direction = u;
  out.positionCovarianceCm2 = PositionBlock(covariance);
  out.sigmaDirectionRad = std::sqrt(
      std::max(0., (sigma[3] * sigma[3] + st * st * sigma[2] * sigma[2]) / 2.));
  if (tracking) {
    out.hasTrackingMomentum = out.hasTrackingCovariance = out.chargeKnown =
        true;
    out.trackingMomentumGeV = magnitude * u;
    out.trackingCovariance = covariance;
    out.charge = t.charge();
  } else {
    out.hasGeometryCovariance = true;
    out.geometryCovariance = covariance;
  }
  return true;
}
struct Match {
  int track, hit, detector, flag;
  double chi2, time;
  Vec3 position;
  bool shared;
};
} // namespace
#endif
struct CeeKFParticleBenoSource::Implementation {
  AnalysisConfig config;
  BenoInputOptions options;
  BenoInputStatistics stats;
  TrackIndexSelector selector;
  StateConverter converter;
  bool opened = false, failed = false;
  std::string failure;
  Implementation(const AnalysisConfig &c, const BenoInputOptions &o)
      : config(c), options(o) {}
#ifdef CEE_WITH_BENODST
  TChain *chain = nullptr;
  TClonesArray *events = nullptr, *tracks = nullptr, *covs = nullptr,
               *itof = nullptr, *etof = nullptr;
  TClonesArray *itofHits = nullptr, *etofHits = nullptr;
  long long entry = 0;
  ~Implementation() {
    delete chain;
    delete events;
    delete tracks;
    delete covs;
    delete itof;
    delete etof;
    delete itofHits;
    delete etofHits;
  }
  bool Tracking() const {
    return config.finder.momentumMode == kTrackingMomentum;
  }
  void Bind(const char *name, const char *type, TClonesArray *&a) {
    a = new TClonesArray(type);
    chain->SetBranchStatus(name, 1);
    chain->SetBranchStatus((std::string(name) + "*").c_str(), 1);
    if (chain->SetBranchAddress(name, &a) < 0)
      throw std::runtime_error(std::string("cannot bind Beno branch: ") + name);
  }
  void Open() {
    if (opened)
      return;
    if (config.finder.momentumMode != kTrackingMomentum &&
        config.finder.momentumMode != kTofMomentum)
      throw std::runtime_error("select tof or tracking explicitly for BenoDst");
    if (!std::isfinite(config.finder.magneticFieldByTesla) ||
        (!Tracking() && std::fabs(config.finder.magneticFieldByTesla) >
                            config.finder.zeroFieldThresholdTesla))
      throw std::runtime_error("Beno TOF mode requires zero field");
    if (config.trackInputMode !=
        (Tracking() ? kCeeMagneticFieldTrackingMode : kCeeZeroFieldTofMode))
      throw std::runtime_error("Beno trackInputMode and momentumMode disagree");
    if (config.maximumEvents < -1 || !config.maximumTracksPerEvent)
      throw std::runtime_error("invalid input limits");
    // No NDF/convergence flags exist in this schema. Do not silently reject all
    // tracks.
    if (config.maximumTrackChi2Ndf < 1.e29 ||
        config.adapter.requireConvergedTrackFit)
      throw std::runtime_error(
          "BenoDst has no NDF/convergence flag; disable these source cuts");
    if (!Tracking()) {
      std::string error;
      if (!ValidateDataModeTiming(config, &error))
        throw std::runtime_error(error);
      if (config.t0.source != kFixedT0 || config.t0.fixedTimeNs != 0. ||
          config.tofRawToNsScale != 1. || config.tofOffsetNs != 0.)
        throw std::runtime_error(
            "Beno MatchInfo TOF is already calibrated ns minus T0: use fixed "
            "origin 0, scale 1, offset 0");
      if (!(config.tofTimeSigmaNs > 0.) ||
          !std::isfinite(config.tofTimeSigmaNs) ||
          !(config.tofPathLengthSigmaCm > 0.) ||
          !std::isfinite(config.tofPathLengthSigmaCm))
        throw std::runtime_error("invalid TOF time/path resolutions");
      if (config.t0.sigmaNs > 0. &&
          !config.finder.propagateSharedTimeUncertainty)
        throw std::runtime_error(
            "nonzero Beno shared T0 error requires shared-time propagation");
    }
    if (!(options.geometryPositionToleranceCm > 0.) ||
        !std::isfinite(options.geometryPositionToleranceCm) ||
        !(options.geometryDirectionToleranceRad > 0.) ||
        !std::isfinite(options.geometryDirectionToleranceRad) ||
        !(options.maximumMatchChi2 >= 0.))
      throw std::runtime_error("invalid Beno options");
    const auto files = InputFiles::ResolveInputFiles(config.inputPath,
                                                     config.maximumInputFiles);
    if (files.empty())
      throw std::runtime_error("no Beno ROOT inputs");
    bool hasI = false, hasE = false;
    for (std::size_t i = 0; i < files.size(); ++i) {
      std::unique_ptr<TFile> f(TFile::Open(files[i].c_str(), "READ"));
      if (!f || f->IsZombie())
        throw std::runtime_error("cannot open " + files[i]);
      TTree *t = dynamic_cast<TTree *>(f->Get(config.treeName.c_str()));
      if (!t || !t->GetBranch("Event") || !t->GetBranch("TpcTrack") ||
          !t->GetBranch("TpcTrackCovMatrix"))
        throw std::runtime_error(
            "missing Beno Event/TpcTrack/TpcTrackCovMatrix in " + files[i]);
      if (!Tracking()) {
        bool hi = t->GetBranch("TpcItofMatchInfo"),
             he = t->GetBranch("TpcEtofMatchInfo");
        if ((!hi && !he) || (hi && !t->GetBranch("ItofHit")) ||
            (he && !t->GetBranch("EtofHit")))
          throw std::runtime_error("missing Beno TOF MatchInfo/hit branch in " +
                                   files[i]);
        if (i && (hi != hasI || he != hasE))
          throw std::runtime_error(
              "Beno TOF schema differs across input files");
        hasI = hi;
        hasE = he;
      }
    }
    chain = new TChain(config.treeName.c_str());
    for (const auto &f : files)
      if (!chain->Add(f.c_str()))
        throw std::runtime_error("cannot chain " + f);
    chain->SetBranchStatus("*", 0);
    Bind("Event", "CeeBenoEvent", events);
    Bind("TpcTrack", "CeeBenoTpcTrack", tracks);
    Bind("TpcTrackCovMatrix", "CeeBenoTpcTrackCovMatrix", covs);
    if (hasI) {
      Bind("TpcItofMatchInfo", "CeeBenoTpcItofMatchInfo", itof);
      Bind("ItofHit", "CeeBenoItofHit", itofHits);
    }
    if (hasE) {
      Bind("TpcEtofMatchInfo", "CeeBenoTpcEtofMatchInfo", etof);
      Bind("EtofHit", "CeeBenoEtofHit", etofHits);
    }
    opened = true;
  }
  std::vector<Match> Matches() const {
    std::vector<Match> result;
    for (int i = 0; itof && i < itof->GetEntriesFast(); ++i) {
      const auto *m = dynamic_cast<CeeBenoTpcItofMatchInfo *>(itof->At(i));
      if (!m)
        continue;
      result.push_back({m->tpcTrackIndex(), m->itofHitIndex(), 0,
                        m->tpcItofMatchFlag(), m->tpcItofMatchChi2(), m->itof(),
                        Vector(m->itofHitPos()), m->hitShared()});
    }
    for (int i = 0; etof && i < etof->GetEntriesFast(); ++i) {
      const auto *m = dynamic_cast<CeeBenoTpcEtofMatchInfo *>(etof->At(i));
      if (!m)
        continue;
      result.push_back({m->tpcTrackIndex(), m->etofHitIndex(), 1,
                        m->tpcEtofMatchFlag(), m->tpcEtofMatchChi2(), m->etof(),
                        Vector(m->etofHitPos()), m->hitShared()});
    }
    return result;
  }
  ReadStatus Next(AnalysisEvent &out) {
    while (entry < chain->GetEntries() &&
           (config.maximumEvents < 0 || entry < config.maximumEvents)) {
      const long long index = entry++;
      ++stats.entries;
      out = AnalysisEvent();
      if (chain->GetEntry(index) <= 0)
        throw std::runtime_error("Beno event read failed at " +
                                 std::to_string(index));
      if (events->GetEntriesFast() == 0)
        continue; // upstream skipped event
      const auto *e = dynamic_cast<CeeBenoEvent *>(events->At(0));
      if (!e || e->primaryVertexIndex() < 0)
        continue;
      const Vec3 pv = Vector(e->primaryVertex()),
                 err = Vector(e->primaryVertexError());
      if (!IsFinite(pv) || !IsFinite(err) || err[0] <= 0. || err[1] <= 0. ||
          err[2] <= 0.)
        continue;
      out.eventId = index;
      out.inputTrackCount = tracks->GetEntriesFast();
      out.input.primaryVertex.positionCm = pv;
      out.input.primaryVertex.covarianceCm2 = Mat3();
      for (int i = 0; i < 3; ++i)
        out.input.primaryVertex.covarianceCm2.m[i][i] = err[i] * err[i];
      out.input.startTimeNs = 0.;
      out.input.sigmaStartTimeNs = Tracking() ? 0. : config.t0.sigmaNs;
      std::vector<int> indices;
      if (selector)
        indices = selector(index, *tracks);
      else
        for (int i = 0; i < tracks->GetEntriesFast(); ++i)
          indices.push_back(i);
      // Build indexes once per event; no track x hit search or geometric
      // rematching.
      const auto matches = Matches();
      std::map<int, std::vector<const Match *>> byTrack;
      std::map<std::pair<int, int>, std::set<int>> claims;
      for (const auto &m : matches) {
        const auto *hits = m.detector == 0 ? itofHits : etofHits;
        if (m.track < 0 || m.track >= tracks->GetEntriesFast() || !hits ||
            m.hit < 0 || m.hit >= hits->GetEntriesFast())
          throw std::runtime_error(
              "Beno MatchInfo contains an out-of-range track/hit index");
        byTrack[m.track].push_back(&m);
        claims[{m.detector, m.hit}].insert(m.track);
      }
      std::set<int> seen, ids;
      for (int i : indices) {
        if (i < 0 || i >= tracks->GetEntriesFast() || !seen.insert(i).second)
          throw std::runtime_error("invalid/duplicate Beno selected index");
        ++stats.tracksConsidered;
        const auto *t = dynamic_cast<CeeBenoTpcTrack *>(tracks->At(i));
        if (!t)
          continue;
        if (t->trackId() < 0 || !ids.insert(t->trackId()).second)
          throw std::runtime_error("invalid/duplicate Beno trackId");
        if (t->nHitsFit() < config.minimumTpcHits)
          continue;
        const int ci = t->covMatrixIndex();
        const auto *cov =
            ci >= 0 && ci < covs->GetEntriesFast()
                ? dynamic_cast<CeeBenoTpcTrackCovMatrix *>(covs->At(ci))
                : nullptr;
        AnalysisTrack a;
        std::string reason;
        if (!cov || cov->tpcTrackIndex() != i ||
            !(converter
                  ? converter(*t, *cov, a.input, reason)
                  : PerigeeState(*t, *cov, Tracking(), options, a.input))) {
          ++stats.invalidGeometry;
          continue;
        }
        auto &x = a.input;
        if (!IsFinite(x.positionCm) || !IsFinite(x.direction) ||
            !(Norm(x.direction) > 0.) ||
            (Tracking() &&
             (!x.hasTrackingMomentum || !x.hasTrackingCovariance ||
              !x.chargeKnown || !IsFinite(x.trackingMomentumGeV) ||
              !(Norm(x.trackingMomentumGeV) > 0.) ||
              !IsPositiveSemidefinite(x.trackingCovariance))) ||
            (!Tracking() &&
             (!x.hasGeometryCovariance || !Finite6(x.geometryCovariance) ||
              !IsPositiveSemidefinite(x.geometryCovariance)))) {
          ++stats.invalidGeometry;
          continue;
        }
        x.id = t->trackId();
        x.lambdaRoleDcaBoundaryCm = config.adapter.lambdaRoleDcaBoundaryCm;
        x.qualityChi2 = Pid(t->trackChi2());
        x.qualityNdf = 0;
        if (!Tracking()) {
          x.hasTrackingMomentum = x.hasTrackingCovariance = x.chargeKnown =
              false;
          x.charge = 0;
          const Match *best = nullptr;
          bool rejectedShared = false;
          for (const Match *m : byTrack[i]) {
            const int partner = m->detector == 0 ? t->itofBestMatchIndex()
                                                 : t->etofBestMatchIndex();
            if (m->hit != partner || (m->flag != 1 && m->flag != 2) ||
                (options.requireUniqueMatch && m->flag != 1))
              continue;
            if (!std::isfinite(m->chi2) || m->chi2 < 0. ||
                m->chi2 > options.maximumMatchChi2 || !std::isfinite(m->time) ||
                m->time <= 0. || !IsFinite(m->position))
              continue;
            if (options.rejectSharedHits &&
                (m->shared || claims[{m->detector, m->hit}].size() > 1)) {
              rejectedShared = true;
              continue;
            }
            const int preferred = options.preferItof ? 0 : 1;
            if (!best ||
                (m->detector == preferred && best->detector != preferred) ||
                (m->detector == best->detector && m->chi2 < best->chi2))
              best = m;
          }
          if (!best) {
            if (rejectedShared)
              ++stats.sharedMatch;
            else
              ++stats.missingMatch;
            continue;
          }
          x.hasTof = true;
          x.tofPositionCm = best->position;
          x.tofTimeNs = best->time;
          x.sigmaTofTimeNs = config.tofTimeSigmaNs;
          x.sigmaPathLengthCm = config.tofPathLengthSigmaCm;
        }
        a.nHits = t->nHitsFit();
        a.dEdx = Pid(t->dEdx());
        a.nSigmaProton = Pid(t->nSigmaProton());
        a.nSigmaPion = Pid(t->nSigmaPion());
        a.nSigmaKaon = Pid(t->nSigmaKaon());
        // status() is PV association, not convergence; chi2 is not chi2/NDF.
        out.tracks.push_back(a);
        ++stats.tracksAccepted;
        if (out.tracks.size() > config.maximumTracksPerEvent)
          throw std::runtime_error("Beno maximumTracksPerEvent exceeded");
      }
      return kEvent;
    }
    out = AnalysisEvent();
    return kEnd;
  }
#else
  void Open() {
    throw std::runtime_error(
        "BenoDst support is not built: make WITH_BENODST=1 in a CeeRoot "
        "environment containing libCeeBenoDst");
  }
  ReadStatus Next(AnalysisEvent &) { return kError; }
#endif
};
CeeKFParticleBenoSource::CeeKFParticleBenoSource(const AnalysisConfig &c,
                                                 const BenoInputOptions &o)
    : fImplementation(new Implementation(c, o)) {}
CeeKFParticleBenoSource::~CeeKFParticleBenoSource() { delete fImplementation; }
void CeeKFParticleBenoSource::SetTrackIndexSelector(TrackIndexSelector s) {
  fImplementation->selector = s;
}
void CeeKFParticleBenoSource::SetStateConverter(StateConverter s) {
  fImplementation->converter = s;
}
const BenoInputStatistics &CeeKFParticleBenoSource::GetStatistics() const {
  return fImplementation->stats;
}
bool CeeKFParticleBenoSource::Open(std::string &error) {
  error.clear();
  if (fImplementation->failed) {
    error = fImplementation->failure;
    return false;
  }
  try {
    fImplementation->Open();
    return true;
  } catch (const std::exception &e) {
    error = e.what();
    fImplementation->failed = true;
    fImplementation->failure = error;
    return false;
  }
}
CeeKFParticleEventSource::ReadStatus
CeeKFParticleBenoSource::Next(AnalysisEvent &out, std::string &error) {
  out = AnalysisEvent();
  if (!Open(error))
    return kError;
  try {
    return fImplementation->Next(out);
  } catch (const std::exception &e) {
    error = e.what();
    fImplementation->failed = true;
    fImplementation->failure = error;
    out = AnalysisEvent();
    return kError;
  }
}
} // namespace CeeV0Analysis
