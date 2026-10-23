#ifndef CEE_KFP_INPUT_H
#define CEE_KFP_INPUT_H
#include "CeeKFParticleBenoSource.h"
#include "CeeKFParticleRootSource.h"
#include <cmath>
#include <memory>
namespace CeeV0Analysis {
// Independent choices: source format and physical momentum source.
inline bool ConfigureInputModes(AnalysisConfig &c, const std::string &format,
                                const std::string &mode, double by,
                                std::string &error) {
  error.clear();
  if (format != "raw" && format != "benodst") {
    error = "inputFormat must be raw or benodst";
    return false;
  }
  if (mode != "tof" && mode != "tracking" && mode != "auto") {
    error = "momentumMode must be tof, tracking or auto";
    return false;
  }
  if (!std::isfinite(by)) {
    error = "By must be finite (tesla)";
    return false;
  }
  const bool field = std::fabs(by) > c.finder.zeroFieldThresholdTesla;
  const bool tracking = mode == "tracking" || (mode == "auto" && field);
  if (!tracking && field) {
    error = "TOF momentum mode requires zero field";
    return false;
  }
  c.finder.magneticFieldByTesla = by;
  c.trackInputMode = tracking ? CeeV0Lite::kCeeMagneticFieldTrackingMode
                              : CeeV0Lite::kCeeZeroFieldTofMode;
  c.finder.momentumMode =
      tracking ? CeeV0Lite::kTrackingMomentum : CeeV0Lite::kTofMomentum;
  // Beno MatchInfo times already use the selected measured T0 as origin.
  // Zeroing the coordinate origin does NOT zero its common uncertainty.
  if (format == "benodst") {
    c.t0.source = kFixedT0;
    c.t0.fixedTimeNs = 0.;
    c.t0.rawToNsScale = 1.;
    c.t0.offsetNs = 0.;
  }
  return true;
}
inline std::unique_ptr<CeeKFParticleEventSource>
OpenCeeInput(const std::string &format, const AnalysisConfig &c,
             std::string &error,
             const BenoInputOptions &options = BenoInputOptions()) {
  error.clear();
  if (format == "raw") {
    std::unique_ptr<CeeKFParticleRootSource> p(new CeeKFParticleRootSource(c));
    if (!p->Open(error))
      return {};
    return std::unique_ptr<CeeKFParticleEventSource>(p.release());
  }
  if (format == "benodst") {
    std::unique_ptr<CeeKFParticleBenoSource> p(
        new CeeKFParticleBenoSource(c, options));
    if (!p->Open(error))
      return {};
    return std::unique_ptr<CeeKFParticleEventSource>(p.release());
  }
  error = "inputFormat must be raw or benodst";
  return {};
}
} // namespace CeeV0Analysis
#endif
