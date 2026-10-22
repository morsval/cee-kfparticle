#include "CeeV0Analyzer.h"

#include <cmath>

namespace CeeV0Analysis
{

void ConfigureDataMode(AnalysisConfig& config, DataMode mode)
{
  config.dataMode = mode;
  config.t0.fixedTimeNs = 0.;
  config.t0.rawToNsScale = 1.;
  config.t0.offsetNs = 0.;
  config.t0.aggregation = kMedianT0;
  if (mode == kSimulationData) {
    // Ordinary MC TOF times use the exact generator collision time as their
    // coordinate origin. This is deliberately not a T0-detector measurement.
    config.enableSimulationTruth = true;
    config.t0.source = kFixedT0;
    config.t0.sigmaNs = 0.;
    config.finder.additionalDecayTimeSigmaNs = 0.;
    config.finder.propagateSharedTimeUncertainty = false;
  }
  else {
    config.enableSimulationTruth = false;
    config.t0.source = kT0Hit;
    // Placeholders retained from the initial analysis configuration. They
    // must be replaced or validated against the real T0 calibration.
    config.t0.sigmaNs = 0.10;
    config.finder.additionalDecayTimeSigmaNs = 0.10;
    config.finder.propagateSharedTimeUncertainty = true;
  }
}

AnalysisConfig::AnalysisConfig()
  : outputPath("CeeV0Lambda.root"), treeName("ceesim"),
    dataMode(kSimulationData), enableSimulationTruth(true), maximumEvents(-1),
    maximumInputFiles(-1), debugLevel(0), tofRawToNsScale(1.),
    tofOffsetNs(0.), tofTimeSigmaNs(0.10), tofPathLengthSigmaCm(0.5),
    maximumItofMatchDistanceCm(10.), maximumEtofMatchDistanceCm(10.),
    minimumTpcHits(0), maximumTrackChi2Ndf(1.e30),
    maximumTracksPerEvent(1000u),
    trackInputMode(CeeV0Lite::kCeeZeroFieldTofMode)
{
  finder.momentumMode = CeeV0Lite::kTofMomentum;
  finder.magneticFieldByTesla = 0.;
  finder.reconstructK0Short = false;
  finder.reconstructLambda = true;
  finder.reconstructAntiLambda = false;
  finder.useLambdaDcaRoleAssignment = true;
  finder.requireTofConvergence = true;
  ConfigureDataMode(*this, kSimulationData);
}

bool ApplyCutPreset(const std::string& name, AnalysisConfig& config)
{
  CeeV0Lite::ChannelCuts& cuts = config.finder.lambdaCuts;
  cuts.enabled = true;
  if (name == "legacy") {
    cuts.minimumMassGeV = 0.7;
    cuts.maximumMassGeV = 2.0;
    cuts.maximumDaughterDcaCm = 1.e30;
    cuts.maximumGeometricChi2Ndf = 1.e30;
    cuts.maximumTopologicalChi2Ndf = 1.e30;
    cuts.minimumDecayLengthCm = 0.;
    cuts.minimumDecayLengthSignificance = 0.;
    cuts.minimumPointingCosine = -1.;
    cuts.minimumDaughterChi2FromPv = 0.;
    cuts.minimumDaughterDeviationFromPv = 0.;
    cuts.requireCompleteMotherJacobian = true;
    return true;
  }
  if (name == "loose") {
    cuts = CeeV0Lite::ChannelCuts();
    cuts.enabled = true;
    cuts.minimumMassGeV = 1.05;
    cuts.maximumMassGeV = 1.25;
    return true;
  }
  if (name == "star") {
    cuts.minimumMassGeV = 1.08;
    cuts.maximumMassGeV = 1.18;
    cuts.maximumDaughterDcaCm = 1.;
    cuts.maximumGeometricChi2Ndf = 10.;
    cuts.maximumTopologicalChi2Ndf = 5.;
    cuts.minimumDecayLengthCm = 5.;
    cuts.minimumDecayLengthSignificance = 5.;
    // The active STAR Interface classifies tracks with raw
    // GetDeviationFromVertex(PV) chi2 >= 18.6 as secondary.  CEE uses the same
    // raw-chi2 convention here, though its reduced covariance model means this
    // numerical threshold must still be validated with simulation pulls.
    cuts.minimumPointingCosine = -1.;
    cuts.minimumDaughterChi2FromPv = 18.6;
    cuts.minimumDaughterDeviationFromPv = 0.;
    cuts.requireCompleteMotherJacobian = true;
    return true;
  }
  return false;
}

bool UsesMeasuredT0(const AnalysisConfig& config)
{
  return config.dataMode == kRealData && config.t0.source == kT0Hit;
}

bool ValidateDataModeTiming(const AnalysisConfig& config,
                            std::string* errorMessage)
{
  if (config.dataMode == kRealData) {
    const bool commonValid = std::isfinite(config.t0.sigmaNs) &&
      config.t0.sigmaNs >= 0. &&
      std::isfinite(config.finder.additionalDecayTimeSigmaNs) &&
      config.finder.additionalDecayTimeSigmaNs >= 0.;
    const bool sourceValid =
      (config.t0.source == kFixedT0 &&
       std::isfinite(config.t0.fixedTimeNs)) ||
      (config.t0.source == kT0Hit &&
       std::isfinite(config.t0.rawToNsScale) &&
       config.t0.rawToNsScale > 0. &&
       std::isfinite(config.t0.offsetNs));
    if (commonValid && sourceValid) return true;
    if (errorMessage) {
      *errorMessage =
        "real mode requires a finite fixed time or a positive measured-T0 "
        "scale, plus finite non-negative common timing uncertainties";
    }
    return false;
  }
  if (config.dataMode != kSimulationData) {
    if (errorMessage) *errorMessage = "unknown analysis data mode";
    return false;
  }
  const bool valid = config.t0.source == kFixedT0 &&
                     config.t0.fixedTimeNs == 0. &&
                     config.t0.rawToNsScale == 1. &&
                     config.t0.offsetNs == 0. &&
                     config.t0.aggregation == kMedianT0 &&
                     config.t0.sigmaNs == 0. &&
                     config.finder.additionalDecayTimeSigmaNs == 0. &&
                     !config.finder.propagateSharedTimeUncertainty;
  if (!valid && errorMessage) {
    *errorMessage =
      "simulation mode uses the exact generator time origin t=0 and must "
      "not read or propagate CeeT0Hit; require t0Source=fixed, fixedT0=0, "
      "T0 scale=1, offset=0, policy=median, sigmaT0=0, "
      "additionalTimeSigma=0 and shared-time propagation off";
  }
  return valid;
}

} // namespace CeeV0Analysis
