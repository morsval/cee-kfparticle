#include "CeeKFParticleChannels.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
namespace CeeV0Analysis {
using namespace CeeV0Lite;
DaughterSpec DaughterSpec::Track(const std::string &r,
                                 const ParticleHypothesis &p) {
  DaughterSpec d;
  d.role = r;
  d.particle = p;
  return d;
}
DaughterSpec DaughterSpec::Composite(const std::string &r,
                                     const std::string &c) {
  DaughterSpec d;
  d.role = r;
  d.sourceChannel = c;
  return d;
}
const ChannelDaughter &ChannelCandidate::Daughter(const std::string &r) const {
  for (const auto &d : daughters)
    if (d.role == r)
      return d;
  throw std::out_of_range("unknown daughter role: " + r);
}
namespace {
struct Input {
  MotherState state;
  int pdg = 0, charge = 0, track = -1, candidate = -1;
  double mass = 0;
  std::vector<int> leaves;
};
double get(const MotherState &s, int i) {
  return i < 3 ? s.positionCm[i] : i < 6 ? s.momentumGeV[i - 3] : s.energyGeV;
}
void set(MotherState &s, int i, double v) {
  if (i < 3)
    s.positionCm[i] = v;
  else if (i < 6)
    s.momentumGeV[i - 3] = v;
  else
    s.energyGeV = v;
}
TrackState track(const Input &x) {
  TrackState t;
  t.id = x.track;
  t.pdg = x.pdg;
  t.charge = x.charge;
  t.massGeV = x.mass;
  t.positionCm = x.state.positionCm;
  t.momentumGeV = x.state.momentumGeV;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      t.covariance.m[i][j] = x.state.covariance.m[i][j];
  return t;
}
bool transport(const Input &x, const Vec3 &v, const FinderConfig &c,
               MotherState &out) {
  TrackState t;
  double path, J[6][6];
  if (!TransportToClosestApproach(track(x), v, c.magneticFieldByTesla,
                                  c.transportConfig, t, path, 0, J))
    return false;
  out = x.state;
  out.positionCm = t.positionCm;
  out.momentumGeV = t.momentumGeV;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j)
      out.covariance.m[i][j] = t.covariance.m[i][j];
    double e = 0;
    for (int k = 0; k < 6; ++k)
      e += J[i][k] * x.state.covariance.m[k][6];
    out.covariance.m[i][6] = out.covariance.m[6][i] = e;
  }
  return true;
}
struct Local {
  MotherState state;
  Vec3 h[2];
  double w[2][2];
  double gain[7][2];
};
bool localize(const std::vector<Input> &inputs, const Vec3 &point,
              const FinderConfig &c, std::vector<Local> &local, Vec3 &vertex,
              Mat3 &vc) {
  local.resize(inputs.size());
  Mat3 info;
  Vec3 rhs;
  for (std::size_t n = 0; n < inputs.size(); ++n) {
    Local &l = local[n];
    if (!transport(inputs[n], point, c, l.state))
      return false;
    Vec3 u;
    if (!Normalize(l.state.momentumGeV, u))
      return false;
    MakeTransverseBasis(u, l.h[0], l.h[1]);
    double r[2][2] = {{0, 0}, {0, 0}};
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b) {
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            r[a][b] += l.h[a][i] * l.state.covariance.m[i][j] * l.h[b][j];
        if (a == b)
          r[a][b] += c.vertexFitterConfig.measurementVarianceFloorCm2;
      }
    if (!InvertSymmetric2(r, l.w) || r[0][0] <= 0 || r[1][1] <= 0)
      return false;
    Mat3 information;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int a = 0; a < 2; ++a)
          for (int b = 0; b < 2; ++b)
            information.m[i][j] += l.h[a][i] * l.w[a][b] * l.h[b][j];
    info = info + information;
    rhs = rhs + Multiply(information, l.state.positionCm);
    for (int i = 0; i < 7; ++i)
      for (int a = 0; a < 2; ++a) {
        l.gain[i][a] = 0;
        for (int j = 0; j < 3; ++j)
          for (int b = 0; b < 2; ++b)
            l.gain[i][a] += l.state.covariance.m[i][j] * l.h[b][j] * l.w[b][a];
      }
  }
  if (!InvertPositiveDefinite3(info, vc))
    return false;
  vertex = Multiply(vc, rhs);
  return IsFinite(vertex);
}
// Iterated batch common-vertex solve against the original independent priors.
// All daughters participate simultaneously; a three-body vertex is NOT a
// sequence of two artificial two-body decay vertices.
bool centralFit(const std::vector<Input> &in, const FinderConfig &c,
                ChannelCandidate &out) {
  Mat3 a;
  Vec3 b;
  for (const auto &x : in) {
    Vec3 u;
    if (!Normalize(x.state.momentumGeV, u))
      return false;
    Mat3 p = Mat3::Identity();
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        p.m[i][j] -= u[i] * u[j];
    a = a + p;
    b = b + Multiply(p, x.state.positionCm);
  }
  Mat3 inv;
  if (!InvertPositiveDefinite3(a, inv))
    return false;
  Vec3 v = Multiply(inv, b);
  Mat3 vc;
  std::vector<Local> l;
  bool converged = false;
  for (int i = 0; i < c.vertexFitterConfig.maximumIterations; ++i) {
    Vec3 next;
    if (!localize(in, v, c, l, next, vc))
      return false;
    double delta = Norm(next - v);
    v = next;
    if (delta < c.vertexFitterConfig.convergenceCm) {
      converged = true;
      break;
    }
  }
  if (!converged)
    return false;
  Vec3 next;
  if (!localize(in, v, c, l, next, vc))
    return false;
  v = next;
  out.mother = MotherState();
  out.mother.positionCm = v;
  out.chi2 = 0;
  out.ndf = 2 * int(in.size()) - 3;
  out.maximumDaughterSeparationCm = 0;
  out.daughters.resize(in.size());
  for (std::size_t n = 0; n < in.size(); ++n) {
    const Local &t = l[n];
    double residual[2] = {Dot(t.h[0], v - t.state.positionCm),
                          Dot(t.h[1], v - t.state.positionCm)};
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        out.chi2 += residual[i] * t.w[i][j] * residual[j];
    MotherState d = t.state;
    for (int i = 3; i < 7; ++i)
      set(d, i,
          get(d, i) + t.gain[i][0] * residual[0] + t.gain[i][1] * residual[1]);
    d.positionCm = v;
    if (in[n].track >= 0)
      d.energyGeV = std::sqrt(in[n].mass * in[n].mass + Norm2(d.momentumGeV));
    if (!std::isfinite(d.energyGeV) || d.energyGeV <= Norm(d.momentumGeV))
      return false;
    // Marginal at fixed linearization, including the uncertain shared vertex.
    double A[7][3] = {{0}};
    for (int i = 0; i < 3; ++i)
      A[i][i] = 1;
    for (int i = 3; i < 7; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 2; ++k)
          A[i][j] += t.gain[i][k] * t.h[k][j];
    Mat7 cond;
    for (int i = 3; i < 7; ++i)
      for (int j = 3; j < 7; ++j) {
        cond.m[i][j] = t.state.covariance.m[i][j];
        for (int k = 0; k < 2; ++k)
          for (int q = 0; q < 3; ++q)
            cond.m[i][j] -=
                t.gain[i][k] * t.h[k][q] * t.state.covariance.m[q][j];
      }
    for (int i = 0; i < 7; ++i)
      for (int j = 0; j < 7; ++j)
        for (int k = 0; k < 3; ++k)
          for (int q = 0; q < 3; ++q)
            cond.m[i][j] += A[i][k] * vc.m[k][q] * A[j][q];
    d.covariance = cond;
    Symmetrize(d.covariance);
    d.covarianceValid = true;
    d.covarianceStatus = kMotherCovarianceValid;
    out.daughters[n].fittedState = d;
    out.mother.momentumGeV = out.mother.momentumGeV + d.momentumGeV;
    out.mother.energyGeV += d.energyGeV;
    for (std::size_t j = 0; j < n; ++j)
      out.maximumDaughterSeparationCm =
          std::max(out.maximumDaughterSeparationCm,
                   Norm(t.state.positionCm - l[j].state.positionCm));
  }
  return true;
}
bool fit(const std::vector<Input> &in, const FinderConfig &c,
         ChannelCandidate &out) {
  if (!centralFit(in, c, out))
    return false;
  // Full-fit numerical propagation, including moving POCA and vertex, not
  // just the fixed-linearization covariance. Composite energy stays a 7th
  // independent coordinate; no Lambda mass constraint is silently imposed.
  Mat7 cov;
  std::vector<Mat7> dcov(in.size());
  for (std::size_t n = 0; n < in.size(); ++n) {
    double J[7][7] = {{0}};
    std::vector<Mat7> dj(in.size());
    for (int col = 0; col < 7; ++col) {
      if (in[n].track >= 0 && col == 6)
        continue;
      double step =
          col < 3
              ? c.vertexFitterConfig.motherJacobianPositionStepCm
              : std::max(
                    c.vertexFitterConfig.motherJacobianMinimumMomentumStepGeV,
                    std::fabs(get(in[n].state, col)) *
                        c.vertexFitterConfig
                            .motherJacobianMomentumRelativeStep);
      auto plus = in, minus = in;
      set(plus[n].state, col, get(in[n].state, col) + step);
      set(minus[n].state, col, get(in[n].state, col) - step);
      if (in[n].track >= 0) {
        plus[n].state.energyGeV = std::sqrt(in[n].mass * in[n].mass +
                                            Norm2(plus[n].state.momentumGeV));
        minus[n].state.energyGeV = std::sqrt(in[n].mass * in[n].mass +
                                             Norm2(minus[n].state.momentumGeV));
      }
      ChannelCandidate p, m;
      if (!centralFit(plus, c, p) || !centralFit(minus, c, m))
        return false;
      for (int row = 0; row < 7; ++row) {
        J[row][col] = (get(p.mother, row) - get(m.mother, row)) / (2 * step);
        for (std::size_t d = 0; d < in.size(); ++d)
          dj[d].m[row][col] = (get(p.daughters[d].fittedState, row) -
                               get(m.daughters[d].fittedState, row)) /
                              (2 * step);
      }
    }
    for (int i = 0; i < 7; ++i)
      for (int j = 0; j < 7; ++j)
        for (int k = 0; k < 7; ++k)
          for (int q = 0; q < 7; ++q) {
            const double ci = in[n].state.covariance.m[k][q];
            cov.m[i][j] += J[i][k] * ci * J[j][q];
            for (std::size_t d = 0; d < in.size(); ++d)
              dcov[d].m[i][j] += dj[d].m[i][k] * ci * dj[d].m[j][q];
          }
  }
  Symmetrize(cov);
  if (!IsFinite(cov) || !IsPositiveSemidefinite(cov))
    return false;
  out.mother.covariance = cov;
  out.mother.covarianceValid = true;
  out.mother.covarianceStatus = kMotherCovarianceValid;
  for (std::size_t d = 0; d < in.size(); ++d) {
    Symmetrize(dcov[d]);
    out.daughters[d].fittedState.covariance = dcov[d];
  }
  const double m2 = out.mother.energyGeV * out.mother.energyGeV -
                    Norm2(out.mother.momentumGeV);
  if (!(m2 > 0))
    return false;
  out.massGeV = std::sqrt(m2);
  double g[7] = {0, 0, 0, 0, 0, 0, out.mother.energyGeV / out.massGeV};
  for (int i = 0; i < 3; ++i)
    g[i + 3] = -out.mother.momentumGeV[i] / out.massGeV;
  const double vm = QuadraticForm(g, cov);
  if (vm < -1e-12)
    return false;
  out.sigmaMassGeV = std::sqrt(std::max(0., vm));
  return true;
}
std::string key(const DaughterSpec &d) {
  if (!d.sourceChannel.empty())
    return "c:" + d.sourceChannel;
  std::ostringstream out;
  out << "t:" << d.particle.pdg << ":" << d.particle.charge << ":"
      << std::setprecision(17) << d.particle.massGeV;
  return out.str();
}
} // namespace
bool CeeKFParticleChannels::Reconstruct(
    const AnalysisEvent &event, const std::vector<std::size_t> &selected,
    std::string &error) {
  candidates.clear();
  error.clear();
  try {
    const auto &fc = config.vertexFitterConfig;
    if (fc.maximumIterations < 1 || !(fc.convergenceCm > 0) ||
        !(fc.measurementVarianceFloorCm2 > 0) ||
        !(fc.motherJacobianPositionStepCm > 0) ||
        !(fc.motherJacobianMomentumRelativeStep > 0) ||
        !(fc.motherJacobianMinimumMomentumStepGeV > 0))
      throw std::runtime_error("invalid generic fitter configuration");
    if (!pid)
      throw std::runtime_error(
          "channel reconstruction requires an explicit PID selector");
    if (config.momentumMode != kTrackingMomentum)
      throw std::runtime_error(
          "generic channels require tracking momentum; zero-field TOF outer "
          "iteration is supported by the legacy two-body engine only");
    if (!std::isfinite(config.magneticFieldByTesla) ||
        config.transportConfig.positionProcessNoiseCm2PerCm != 0 ||
        config.transportConfig.angularProcessNoiseRad2PerCm != 0)
      throw std::runtime_error(
          "invalid field or unsupported process noise in full-fit covariance");
    std::map<std::string, std::size_t> names;
    for (std::size_t i = 0; i < channels.size(); ++i) {
      const auto &ch = channels[i];
      if (ch.name.empty() || !ch.motherPdg || ch.daughters.size() < 2 ||
          !names.emplace(ch.name, i).second)
        throw std::runtime_error(
            "invalid/duplicate channel or fewer than two daughters");
      std::set<std::string> roles;
      for (const auto &d : ch.daughters) {
        if (d.role.empty() || !roles.insert(d.role).second)
          throw std::runtime_error("daughter role must be nonempty and unique");
        if (d.sourceChannel.empty() &&
            (!d.particle.pdg || !(d.particle.massGeV > 0) ||
             !std::isfinite(d.particle.massGeV)))
          throw std::runtime_error("invalid daughter hypothesis");
      }
    }
    std::vector<int> order, mark(channels.size(), 0);
    std::function<void(int)> visit = [&](int i) {
      if (mark[i] == 1)
        throw std::runtime_error("cyclic channel dependencies");
      if (mark[i] == 2)
        return;
      mark[i] = 1;
      for (const auto &d : channels[i].daughters)
        if (!d.sourceChannel.empty()) {
          auto it = names.find(d.sourceChannel);
          if (it == names.end())
            throw std::runtime_error("unknown composite source channel");
          visit(it->second);
        }
      mark[i] = 2;
      order.push_back(i);
    };
    for (std::size_t i = 0; i < channels.size(); ++i)
      visit(i);
    std::set<std::size_t> seen;
    std::set<int> ids;
    std::vector<std::pair<std::size_t, TrackInput>> tracks;
    for (auto i : selected) {
      if (i >= event.tracks.size() || !seen.insert(i).second)
        throw std::runtime_error("invalid/duplicate track index");
      const auto &t = event.tracks[i].input;
      if (t.id < 0 || !ids.insert(t.id).second)
        throw std::runtime_error("track IDs must be nonnegative and unique");
      if (!t.hasTrackingMomentum || !t.hasTrackingCovariance || !t.chargeKnown)
        throw std::runtime_error("generic channels require charge, physical "
                                 "momentum and full covariance");
      TrackInput prepared;
      if (PrepareTrackReference(t, event.input, config, prepared))
        tracks.emplace_back(i, prepared);
    }
    std::map<std::string, std::vector<Input>> cache;
    std::size_t combinations = 0;
    for (int idx : order) {
      const auto &ch = channels[idx];
      std::vector<std::vector<Input>> pools;
      for (const auto &d : ch.daughters) {
        std::vector<Input> pool;
        if (!d.sourceChannel.empty()) {
          for (const auto &m : candidates)
            if (m.channel == d.sourceChannel) {
              Input x;
              x.state = m.mother;
              x.pdg = m.pdg;
              x.charge = m.charge;
              x.candidate = m.id;
              x.mass = m.massGeV;
              x.leaves = m.leafTrackIds;
              pool.push_back(x);
            }
        } else {
          const auto k = key(d);
          auto cached = cache.find(k);
          if (cached != cache.end())
            pool = cached->second;
          else {
            for (const auto &pair : tracks) {
              const auto &t = pair.second;
              if (t.charge != d.particle.charge ||
                  !pid(event.tracks[pair.first], d.particle, event))
                continue;
              Input x;
              x.pdg = d.particle.pdg;
              x.charge = d.particle.charge;
              x.mass = d.particle.massGeV;
              x.track = t.id;
              x.leaves.push_back(t.id);
              x.state.positionCm = t.positionCm;
              x.state.momentumGeV = t.trackingMomentumGeV;
              x.state.energyGeV =
                  std::sqrt(x.mass * x.mass + Norm2(x.state.momentumGeV));
              for (int a = 0; a < 6; ++a)
                for (int b = 0; b < 6; ++b)
                  x.state.covariance.m[a][b] = t.trackingCovariance.m[a][b];
              for (int a = 0; a < 6; ++a) {
                double ce = 0;
                for (int b = 0; b < 3; ++b)
                  ce += t.trackingCovariance.m[a][b + 3] *
                        x.state.momentumGeV[b] / x.state.energyGeV;
                x.state.covariance.m[a][6] = x.state.covariance.m[6][a] = ce;
              }
              for (int b = 0; b < 3; ++b)
                x.state.covariance.m[6][6] += x.state.covariance.m[6][b + 3] *
                                              x.state.momentumGeV[b] /
                                              x.state.energyGeV;
              x.state.covarianceValid = true;
              x.state.covarianceStatus = kMotherCovarianceValid;
              if (!IsPositiveSemidefinite(x.state.covariance))
                throw std::runtime_error("invalid track covariance");
              pool.push_back(x);
            }
            cache[k] = pool;
          }
        }
        pools.push_back(pool);
      }
      std::vector<Input> combo;
      std::set<int> used;
      std::function<void(std::size_t)> enumerate = [&](std::size_t depth) {
        if (depth == pools.size()) {
          if (++combinations > maximumCombinations)
            throw std::runtime_error(
                "maximumCombinations exceeded; no truncated event is returned");
          ChannelCandidate m;
          if (!fit(combo, config, m))
            return;
          m.channel = ch.name;
          m.pdg = ch.motherPdg;
          m.leafTrackIds.assign(used.begin(), used.end());
          for (std::size_t i = 0; i < combo.size(); ++i) {
            auto &d = m.daughters[i];
            const auto &x = combo[i];
            d.role = ch.daughters[i].role;
            d.pdg = x.pdg;
            d.charge = x.charge;
            d.trackId = x.track;
            d.candidateId = x.candidate;
            d.inputState = x.state;
            m.charge += x.charge;
            TrackState at;
            double path;
            if (TransportToClosestApproach(track(x),
                                           event.input.primaryVertex.positionCm,
                                           config.magneticFieldByTesla,
                                           config.transportConfig, at, path))
              d.dcaToPvCm =
                  Norm(at.positionCm - event.input.primaryVertex.positionCm);
          }
          if (ch.fitCollisionPv) {
            ProductionVertexInput p;
            p.mother = m.mother;
            p.vertex = event.input.primaryVertex;
            if (FitProductionVertex(p, config.productionVertexConfig, m.chi2,
                                    m.ndf, m.motherTopo, m.charge,
                                    config.magneticFieldByTesla))
              GetFlightObservables(m.motherTopo, m.flight);
          }
          if (ch.select && !ch.select(m, event))
            return;
          if (config.maximumCandidates &&
              candidates.size() >= config.maximumCandidates)
            throw std::runtime_error(
                "maximumCandidates exceeded; no truncated event is returned");
          m.id = int(candidates.size());
          candidates.push_back(m);
          return;
        }
        for (const auto &x : pools[depth]) {
          bool overlap = false;
          for (int id : x.leaves)
            if (used.count(id))
              overlap = true;
          if (overlap)
            continue;
          // Identical hypotheses are unordered, regardless of role labels.
          bool duplicate = false;
          for (std::size_t j = 0; j < depth; ++j)
            if (key(ch.daughters[j]) == key(ch.daughters[depth]) &&
                combo[j].leaves >= x.leaves)
              duplicate = true;
          if (duplicate)
            continue;
          for (int id : x.leaves)
            used.insert(id);
          combo.push_back(x);
          enumerate(depth + 1);
          combo.pop_back();
          for (int id : x.leaves)
            used.erase(id);
        }
      };
      enumerate(0);
    }
    return true;
  } catch (const std::exception &e) {
    error = e.what();
    candidates.clear();
    return false;
  }
}
} // namespace CeeV0Analysis
