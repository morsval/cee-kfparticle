#ifndef CEE_BENO_ADAPTER_EXAMPLE_H
#define CEE_BENO_ADAPTER_EXAMPLE_H
// Use inside a FairTask after CeeBenoDstReader::init()/event loading.
// Compile against CeeRoot/CeeBenoDst/{reader,edm}; this is an adapter hook,
// NOT a self-contained ROOT file reader or a guessed 5D->6D conversion.
#include "CeeBenoDstReader.h"
#include "CeeBenoEvent.h"
#include "CeeBenoTpcTrack.h"
#include "CeeBenoTpcTrackCovMatrix.h"
#include "CeeKFParticleEvent.h"
#include <functional>
#include <set>
namespace CeeV0Analysis {
using BenoStateConverter = std::function<bool(
    const CeeBenoTpcTrack &, const CeeBenoTpcTrackCovMatrix &,
    CeeV0Lite::TrackInput &, std::string &)>;
inline bool FillBenoTrackingEvent(CeeBenoDstReader &reader,
                                  long long chainEntry,
                                  const std::vector<int> &selectedBenoIndices,
                                  const BenoStateConverter &convert,
                                  AnalysisEvent &result, std::string &error) {
  result = AnalysisEvent();
  error.clear();
  const auto *e = reader.event();
  if (!e || e->primaryVertexIndex() < 0 || !convert) {
    error = "Beno event/PV/state converter missing";
    return false;
  }
  AnalysisEvent out;
  out.eventId = chainEntry;
  out.inputTrackCount = reader.numberOfTpcTracks();
  const auto v = e->primaryVertex(), s = e->primaryVertexError();
  out.input.primaryVertex.positionCm = CeeV0Lite::Vec3(v.X(), v.Y(), v.Z());
  for (int i = 0; i < 3; ++i) {
    if (!std::isfinite(v[i]) || !std::isfinite(s[i]) || s[i] < 0) {
      error = "invalid Beno PV";
      return false;
    }
    out.input.primaryVertex.covarianceCm2.m[i][i] = s[i] * s[i];
  }
  std::set<int> seen;
  for (int i : selectedBenoIndices) {
    if (i < 0 || !seen.insert(i).second) {
      error = "invalid/duplicate Beno list index";
      return false;
    }
    const auto *t = reader.tpcTrack(i);
    if (!t || !t->hasCovMatrix()) {
      error = "Beno track/covariance missing";
      return false;
    }
    const auto *cov = reader.tpcTrackCovMatrix(t->covMatrixIndex());
    if (!cov || !cov->isValid() || cov->tpcTrackIndex() != i) {
      error = "Beno covariance link invalid";
      return false;
    }
    AnalysisTrack a;
    // Caller supplies a physically validated state and C6=J*C5*J^T at the
    // SAME reference point, including any required rigidity/charge conversion.
    if (!convert(*t, *cov, a.input, error))
      return false;
    if (!a.input.hasTrackingMomentum || !a.input.hasTrackingCovariance ||
        !a.input.chargeKnown) {
      error = "Beno converter omitted tracking state";
      return false;
    }
    a.input.id = t->trackId();
    a.nHits = t->nHitsFit();
    a.dEdx = t->dEdx();
    a.input.qualityChi2 = t->trackChi2();
    a.input.qualityNdf = 0;
    // No NDF / convergence flag in this payload: leave chi2Ndf NaN and
    // fitConverged false; do NOT reinterpret t->status() as convergence.
    auto pid = [](double x) {
      return x == CeeBenoConst::kInvalidValue
                 ? std::numeric_limits<double>::quiet_NaN()
                 : x;
    };
    a.nSigmaProton = pid(t->nSigmaProton());
    a.nSigmaPion = pid(t->nSigmaPion());
    a.nSigmaKaon = pid(t->nSigmaKaon());
    out.tracks.push_back(a);
  }
  result = out;
  return true;
}
} // namespace CeeV0Analysis
#endif
