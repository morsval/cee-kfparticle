#include "CeeV0Analyzer.h"

#include <cassert>
#include <iostream>
#include <string>

int main()
{
  using namespace CeeV0Analysis;

  AnalysisConfig simulation;
  assert(simulation.dataMode == kSimulationData);
  assert(simulation.t0.source == kFixedT0);
  assert(simulation.t0.fixedTimeNs == 0.);
  assert(simulation.t0.rawToNsScale == 1.);
  assert(simulation.t0.offsetNs == 0.);
  assert(simulation.t0.aggregation == kMedianT0);
  assert(simulation.t0.sigmaNs == 0.);
  assert(simulation.finder.additionalDecayTimeSigmaNs == 0.);
  assert(!simulation.finder.propagateSharedTimeUncertainty);
  assert(!UsesMeasuredT0(simulation));
  assert(ValidateDataModeTiming(simulation));
  // The stop-detector resolution remains active; it is not a T0 uncertainty.
  assert(simulation.tofTimeSigmaNs > 0.);

  AnalysisConfig starPreset = simulation;
  assert(ApplyCutPreset("star", starPreset));
  assert(starPreset.finder.lambdaCuts.maximumTopologicalChi2Ndf == 5.);
  assert(starPreset.finder.lambdaCuts.minimumDaughterChi2FromPv == 18.6);
  assert(starPreset.finder.lambdaCuts.minimumDaughterDeviationFromPv == 0.);
  assert(starPreset.finder.lambdaCuts.requireCompleteMotherJacobian);
  AnalysisConfig legacyPreset = simulation;
  legacyPreset.finder.lambdaCuts.minimumDaughterDeviationFromPv = 9.;
  legacyPreset.finder.lambdaCuts.requireCompleteMotherJacobian = false;
  assert(ApplyCutPreset("legacy", legacyPreset));
  assert(legacyPreset.finder.lambdaCuts.minimumDaughterChi2FromPv == 0.);
  assert(legacyPreset.finder.lambdaCuts.minimumDaughterDeviationFromPv == 0.);
  assert(legacyPreset.finder.lambdaCuts.requireCompleteMotherJacobian);
  legacyPreset.finder.lambdaCuts.minimumDaughterChi2FromPv = 7.;
  legacyPreset.finder.lambdaCuts.minimumDaughterDeviationFromPv = 3.;
  assert(ApplyCutPreset("loose", legacyPreset));
  assert(legacyPreset.finder.lambdaCuts.minimumDaughterChi2FromPv == 0.);
  assert(legacyPreset.finder.lambdaCuts.minimumDaughterDeviationFromPv == 0.);
  assert(legacyPreset.finder.lambdaCuts.requireCompleteMotherJacobian);

  std::string error;
  AnalysisConfig invalid = simulation;
  invalid.t0.source = kT0Hit;
  assert(!ValidateDataModeTiming(invalid, &error));
  assert(!error.empty());
  invalid = simulation;
  invalid.t0.fixedTimeNs = 1.;
  assert(!ValidateDataModeTiming(invalid));
  invalid = simulation;
  invalid.t0.rawToNsScale = 1.e-3;
  assert(!ValidateDataModeTiming(invalid));
  invalid = simulation;
  invalid.t0.offsetNs = 1.;
  assert(!ValidateDataModeTiming(invalid));
  invalid = simulation;
  invalid.t0.aggregation = kMeanT0;
  assert(!ValidateDataModeTiming(invalid));
  invalid = simulation;
  invalid.t0.sigmaNs = 0.1;
  assert(!ValidateDataModeTiming(invalid));
  invalid = simulation;
  invalid.finder.additionalDecayTimeSigmaNs = 0.1;
  assert(!ValidateDataModeTiming(invalid));
  invalid = simulation;
  invalid.finder.propagateSharedTimeUncertainty = true;
  assert(!ValidateDataModeTiming(invalid));

  AnalysisConfig real;
  ConfigureDataMode(real, kRealData);
  assert(real.dataMode == kRealData);
  assert(!real.enableSimulationTruth);
  assert(real.t0.source == kT0Hit);
  assert(real.t0.sigmaNs == 0.1);
  assert(real.finder.additionalDecayTimeSigmaNs == 0.1);
  assert(real.finder.propagateSharedTimeUncertainty);
  assert(UsesMeasuredT0(real));
  assert(ValidateDataModeTiming(real));
  AnalysisConfig invalidReal = real;
  invalidReal.t0.rawToNsScale = 0.;
  assert(!ValidateDataModeTiming(invalidReal));
  invalidReal = real;
  invalidReal.t0.sigmaNs = -1.;
  assert(!ValidateDataModeTiming(invalidReal));

  ConfigureDataMode(real, kSimulationData);
  assert(real.dataMode == kSimulationData);
  assert(real.enableSimulationTruth);
  assert(!UsesMeasuredT0(real));
  assert(ValidateDataModeTiming(real));

  std::cout << "CeeV0Analysis simulation/real timing policy tests passed"
            << std::endl;
  return 0;
}
