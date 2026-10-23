#include "CeeV0Analyzer.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace
{

void Usage(const char* program)
{
  std::cout
    << "Usage: " << program << " --input FILE|LIST|DIR [options]\n\n"
    << "Core options:\n"
    << "  --output FILE                 Output ROOT file (CeeV0Lambda.root)\n"
    << "  --mode simulation|real       Data mode (simulation)\n"
    << "  --momentum tof|tracking      Zero-field TOF or magnetic tracking (tof)\n"
    << "  --by-tesla VALUE             Signed uniform CEE field along y / T (0)\n"
    << "                                Beam is z; +By bends +q moving +z toward -x\n"
    << "  --magnetic-track-reference pv-dca|stored  Default pv-dca; stored for replay\n"
    << "  --reference-max-path-cm VALUE  DCA search half-range (1000 cm)\n"
    << "  --cut-preset legacy|loose|star  First closure defaults to legacy\n"
    << "  --events N                   Maximum events, -1 means all\n"
    << "  --files N                    Maximum input files, -1 means all\n"
    << "  --tree NAME                  Input TTree/TChain name (ceesim)\n"
    << "  --debug LEVEL                0, 1 or 2\n\n"
    << "Timing and matching:\n"
    << "  Simulation uses exact generator t=0: no CeeT0Hit or shared T0 error\n"
    << "  --t0-source fixed|hit        Real-data T0 source (default hit)\n"
    << "  --fixed-t0-ns VALUE          Real-data fixed start time (0 ns)\n"
    << "  --t0-scale-to-ns VALUE       Mandatory explicit value for real+hit\n"
    << "  --t0-offset-ns VALUE         Offset after scale\n"
    << "  --t0-sigma-ns VALUE          Real-data T0 resolution (0.10 ns placeholder)\n"
    << "  --t0-policy first|earliest|mean|median  Multiple-hit rule (median)\n"
    << "  --tof-time-scale-to-ns VALUE Raw TOF time scale (1)\n"
    << "  --tof-time-offset-ns VALUE   TOF time offset (0 ns)\n"
    << "  --tof-time-sigma-ns VALUE    TOF time resolution (0.10 ns)\n"
    << "  --tof-path-sigma-cm VALUE    Path-length uncertainty (0.5 cm)\n"
    << "  --additional-decay-time-sigma-ns VALUE  Extra real common error (0.10 ns placeholder)\n"
    << "  --max-beta VALUE             Largest accepted TOF beta (0.9995)\n"
    << "  --itof-match-cm VALUE        Maximum iTOF residual (10 cm)\n"
    << "  --etof-match-cm VALUE        Maximum eTOF residual (10 cm)\n\n"
    << "Track/finder options:\n"
    << "  --min-tpc-hits N             Minimum fitted TPC hits (0)\n"
    << "  --max-track-chi2ndf VALUE    Enable track quality cut\n"
    << "  --max-tracks N               Event protection limit (1000)\n"
    << "  --max-candidates N           Event protection limit (100000)\n"
    << "  --require-track-converged    Require converged TPC fit\n"
    << "  --strict-pv-errors           Reject missing/invalid stored PV errors\n"
    << "  --disable-legacy-track-fallback  Require a valid vertex seed\n"
    << "  --allow-nonconverged-tof     Keep finder results without TOF convergence\n"
    << "  --allow-incomplete-mother-jacobian\n"
    << "                                Keep central candidates when only a +/-\n"
    << "                                covariance refit fails (QA; default reject)\n"
    << "  --disable-dca-role           Disable zero-field p/pi DCA role assignment\n"
    << "  --lambda-role-dca-cm VALUE   TPC p/pi role boundary (1.2 cm)\n"
    << "  --fallback-track-sigma-{x,y,z}-cm VALUE\n"
    << "  --fallback-direction-sigma-rad VALUE\n"
    << "  --fallback-pv-sigma-{x,y,z}-cm VALUE\n"
    << "  --mass-min VALUE, --mass-max VALUE, --max-daughter-dca VALUE,\n"
    << "  --max-geo-chi2ndf VALUE, --max-topo-chi2ndf VALUE,\n"
    << "  --min-decay-length VALUE, --min-ldl VALUE,\n"
    << "  --min-pointing-cos VALUE, --min-daughter-pv-chi2 VALUE\n"
    << "  --min-daughter-pv-deviation VALUE  Deprecated sqrt(chi2) cut\n"
    << "  --no-production-state       Disable auxiliary neutral 8D PV fit (enabled by default)\n"
    << "  --no-truth                   Do not bind MCTrack in simulation mode\n";
}

bool ParseDouble(const char* text, double& value)
{
  errno = 0;
  char* end = 0;
  value = std::strtod(text, &end);
  return errno == 0 && end && *end == '\0' && std::isfinite(value);
}

bool ParseLongLong(const char* text, long long& value)
{
  errno = 0;
  char* end = 0;
  value = std::strtoll(text, &end, 10);
  return errno == 0 && end && *end == '\0';
}

bool Value(int argc, char** argv, int& index, const char*& output)
{
  if (index + 1 >= argc) return false;
  output = argv[++index];
  return true;
}

} // namespace

int main(int argc, char** argv)
{
  CeeV0Analysis::AnalysisConfig config;
  std::string preset = "legacy";
  std::string mode = "simulation";
  bool t0ScaleExplicit = false;
  bool t0OffsetExplicit = false;
  bool t0PolicyExplicit = false;
  bool dcaRoleExplicitlyDisabled = false;

  // Presets are applied before the second-pass individual overrides, so CLI
  // cut overrides have the same result regardless of argument order.
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--cut-preset" && i + 1 < argc) preset = argv[++i];
    else if (std::string(argv[i]) == "--mode" && i + 1 < argc) mode = argv[++i];
  }
  if (mode == "real") {
    CeeV0Analysis::ConfigureDataMode(config, CeeV0Analysis::kRealData);
  }
  else if (mode != "simulation") {
    std::cerr << "Unknown --mode " << mode << std::endl;
    return 1;
  }
  if (!CeeV0Analysis::ApplyCutPreset(preset, config)) {
    std::cerr << "Unknown --cut-preset " << preset << std::endl;
    return 1;
  }

  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    const char* text = 0;
    double number = 0.;
    long long integer = 0;
    if (option == "--help" || option == "-h") {
      Usage(argv[0]);
      return 0;
    }
    else if (option == "--input") {
      if (!Value(argc, argv, i, text)) return 1;
      config.inputPath = text;
    }
    else if (option == "--output") {
      if (!Value(argc, argv, i, text)) return 1;
      config.outputPath = text;
    }
    else if (option == "--tree") {
      if (!Value(argc, argv, i, text)) return 1;
      config.treeName = text;
    }
    else if (option == "--mode" || option == "--cut-preset") {
      if (!Value(argc, argv, i, text)) return 1;
    }
    else if (option == "--events" || option == "--files" ||
             option == "--debug" || option == "--min-tpc-hits" ||
             option == "--max-tracks" || option == "--max-candidates") {
      if (!Value(argc, argv, i, text) || !ParseLongLong(text, integer)) return 1;
      if (option == "--events") {
        config.maximumEvents = integer;
      }
      else if (option == "--files") {
        if (integer < -1 || integer > std::numeric_limits<int>::max()) return 1;
        config.maximumInputFiles = static_cast<int>(integer);
      }
      else if (option == "--debug") {
        if (integer < 0 || integer > std::numeric_limits<int>::max()) return 1;
        config.debugLevel = static_cast<int>(integer);
      }
      else if (option == "--min-tpc-hits") {
        if (integer < 0 || integer > std::numeric_limits<int>::max()) return 1;
        config.minimumTpcHits = static_cast<int>(integer);
      }
      else if (option == "--max-tracks") {
        if (integer <= 0 || static_cast<unsigned long long>(integer) >
            std::numeric_limits<std::size_t>::max()) return 1;
        config.maximumTracksPerEvent = static_cast<std::size_t>(integer);
      }
      else {
        if (integer <= 0 || static_cast<unsigned long long>(integer) >
            std::numeric_limits<std::size_t>::max()) return 1;
        config.finder.maximumCandidates = static_cast<std::size_t>(integer);
      }
    }
    else if (option == "--momentum") {
      if (!Value(argc, argv, i, text)) return 1;
      const std::string value = text;
      if (value == "tof") {
        config.trackInputMode = CeeV0Lite::kCeeZeroFieldTofMode;
        config.finder.momentumMode = CeeV0Lite::kTofMomentum;
        config.finder.allowTofMomentumInMagneticField = false;
        config.finder.useLambdaDcaRoleAssignment = true;
      }
      else if (value == "tracking") {
        config.trackInputMode = CeeV0Lite::kCeeMagneticFieldTrackingMode;
        config.finder.momentumMode = CeeV0Lite::kTrackingMomentum;
        config.finder.allowTofMomentumInMagneticField = false;
        config.finder.useLambdaDcaRoleAssignment = false;
      }
      else return 1;
    }
    else if (option == "--magnetic-track-reference") {
      if (!Value(argc, argv, i, text)) return 1;
      const std::string value = text;
      if (value == "pv-dca") config.finder.usePrimaryVertexDcaReference = true;
      else if (value == "stored") config.finder.usePrimaryVertexDcaReference = false;
      else { std::cerr << "Expected pv-dca or stored." << std::endl; return 1; }
    }
    else if (option == "--t0-source") {
      if (!Value(argc, argv, i, text)) return 1;
      const std::string value = text;
      if (value == "fixed") config.t0.source = CeeV0Analysis::kFixedT0;
      else if (value == "hit") config.t0.source = CeeV0Analysis::kT0Hit;
      else return 1;
    }
    else if (option == "--t0-policy") {
      if (!Value(argc, argv, i, text)) return 1;
      t0PolicyExplicit = true;
      const std::string value = text;
      if (value == "first") config.t0.aggregation = CeeV0Analysis::kFirstT0;
      else if (value == "earliest") config.t0.aggregation = CeeV0Analysis::kEarliestT0;
      else if (value == "mean") config.t0.aggregation = CeeV0Analysis::kMeanT0;
      else if (value == "median") config.t0.aggregation = CeeV0Analysis::kMedianT0;
      else return 1;
    }
    else if (option == "--no-production-state") config.finder.buildProductionState = false;
    else if (option == "--no-truth") config.enableSimulationTruth = false;
    else if (option == "--require-track-converged")
      config.adapter.requireConvergedTrackFit = true;
    else if (option == "--strict-pv-errors")
      config.adapter.allowVertexErrorFallback = false;
    else if (option == "--disable-legacy-track-fallback")
      config.adapter.allowLegacyStateInZeroField = false;
    else if (option == "--allow-nonconverged-tof")
      config.finder.requireTofConvergence = false;
    else if (option == "--allow-incomplete-mother-jacobian")
      config.finder.lambdaCuts.requireCompleteMotherJacobian = false;
    else if (option == "--disable-dca-role")
      dcaRoleExplicitlyDisabled = true;
    else {
      if (!Value(argc, argv, i, text) || !ParseDouble(text, number)) {
        std::cerr << "Unknown option or invalid value: " << option << std::endl;
        return 1;
      }
      if (option == "--by-tesla") config.finder.magneticFieldByTesla = number;
      else if (option == "--reference-max-path-cm")
        config.finder.referenceTransportConfig.maximumAbsPathCm = number;
      else if (option == "--bz-tesla") {
        if (number != 0.) {
          std::cerr << "--bz-tesla is obsolete: CEE uses B=(0,By,0). "
                    << "Use --by-tesla with the signed y field in tesla."
                    << std::endl;
          return 1;
        }
        // Preserve old zero-field scripts without interpreting a z field as y.
        std::cerr << "Deprecated --bz-tesla 0: use --by-tesla 0."
                  << std::endl;
      }
      else if (option == "--fixed-t0-ns") config.t0.fixedTimeNs = number;
      else if (option == "--t0-scale-to-ns") {
        config.t0.rawToNsScale = number;
        t0ScaleExplicit = true;
      }
      else if (option == "--t0-offset-ns") {
        config.t0.offsetNs = number;
        t0OffsetExplicit = true;
      }
      else if (option == "--t0-sigma-ns") config.t0.sigmaNs = number;
      else if (option == "--tof-time-scale-to-ns") config.tofRawToNsScale = number;
      else if (option == "--tof-time-offset-ns") config.tofOffsetNs = number;
      else if (option == "--tof-time-sigma-ns") config.tofTimeSigmaNs = number;
      else if (option == "--tof-path-sigma-cm") config.tofPathLengthSigmaCm = number;
      else if (option == "--additional-decay-time-sigma-ns") config.finder.additionalDecayTimeSigmaNs = number;
      else if (option == "--max-beta") config.finder.maximumBeta = number;
      else if (option == "--itof-match-cm") config.maximumItofMatchDistanceCm = number;
      else if (option == "--etof-match-cm") config.maximumEtofMatchDistanceCm = number;
      else if (option == "--max-track-chi2ndf") config.maximumTrackChi2Ndf = number;
      else if (option == "--lambda-role-dca-cm") config.adapter.lambdaRoleDcaBoundaryCm = number;
      else if (option == "--fallback-track-sigma-x-cm") config.adapter.fallbackTrackSigmaXcm = number;
      else if (option == "--fallback-track-sigma-y-cm") config.adapter.fallbackTrackSigmaYcm = number;
      else if (option == "--fallback-track-sigma-z-cm") config.adapter.fallbackTrackSigmaZcm = number;
      else if (option == "--fallback-direction-sigma-rad") config.adapter.fallbackDirectionSigmaRad = number;
      else if (option == "--fallback-pv-sigma-x-cm") config.adapter.fallbackVertexSigmaXcm = number;
      else if (option == "--fallback-pv-sigma-y-cm") config.adapter.fallbackVertexSigmaYcm = number;
      else if (option == "--fallback-pv-sigma-z-cm") config.adapter.fallbackVertexSigmaZcm = number;
      else if (option == "--mass-min") config.finder.lambdaCuts.minimumMassGeV = number;
      else if (option == "--mass-max") config.finder.lambdaCuts.maximumMassGeV = number;
      else if (option == "--max-daughter-dca") config.finder.lambdaCuts.maximumDaughterDcaCm = number;
      else if (option == "--max-geo-chi2ndf") config.finder.lambdaCuts.maximumGeometricChi2Ndf = number;
      else if (option == "--max-topo-chi2ndf") config.finder.lambdaCuts.maximumTopologicalChi2Ndf = number;
      else if (option == "--min-decay-length") config.finder.lambdaCuts.minimumDecayLengthCm = number;
      else if (option == "--min-ldl") config.finder.lambdaCuts.minimumDecayLengthSignificance = number;
      else if (option == "--min-pointing-cos") config.finder.lambdaCuts.minimumPointingCosine = number;
      else if (option == "--min-daughter-pv-chi2")
        config.finder.lambdaCuts.minimumDaughterChi2FromPv = number;
      else if (option == "--min-daughter-pv-deviation")
        config.finder.lambdaCuts.minimumDaughterDeviationFromPv = number;
      else {
        std::cerr << "Unknown option: " << option << std::endl;
        return 1;
      }
    }
  }

  if (config.inputPath.empty() || config.outputPath.empty() ||
      config.treeName.empty()) {
    Usage(argv[0]);
    return 1;
  }
  // This explicit veto must win independently of CLI option order.
  if (dcaRoleExplicitlyDisabled)
    config.finder.useLambdaDcaRoleAssignment = false;
  std::string timingError;
  if (!CeeV0Analysis::ValidateDataModeTiming(config, &timingError)) {
    std::cerr << timingError << std::endl;
    return 1;
  }
  if (config.dataMode == CeeV0Analysis::kSimulationData && t0ScaleExplicit) {
    std::cerr << "--t0-scale-to-ns is a measured-T0 option and is not used in "
                 "simulation mode."
              << std::endl;
    return 1;
  }
  if (config.t0.source == CeeV0Analysis::kFixedT0 &&
      (t0ScaleExplicit || t0OffsetExplicit || t0PolicyExplicit)) {
    std::cerr << "T0 scale, offset and multi-hit policy are only meaningful "
                 "with --t0-source hit."
              << std::endl;
    return 1;
  }
  if (config.dataMode == CeeV0Analysis::kRealData &&
      config.t0.source == CeeV0Analysis::kT0Hit && !t0ScaleExplicit) {
    std::cerr << "Real-data CeeT0Hit mode requires an explicit "
                 "--t0-scale-to-ns (use 1 if raw values are already ns)."
              << std::endl;
    return 1;
  }
  if (config.maximumEvents < -1 || config.maximumInputFiles < -1 ||
      config.maximumInputFiles == 0 || config.minimumTpcHits < 0 ||
      config.debugLevel < 0 || config.debugLevel > 2 ||
      !(config.t0.sigmaNs >= 0.) || !(config.tofTimeSigmaNs > 0.) ||
      !(config.tofPathLengthSigmaCm > 0.) ||
      !(config.tofRawToNsScale > 0.) ||
      !(config.finder.additionalDecayTimeSigmaNs >= 0.) ||
      !(config.finder.referenceTransportConfig.maximumAbsPathCm > 0.) ||
      !(config.finder.maximumBeta > 0.) ||
      !(config.finder.maximumBeta < 1.) ||
      !(config.maximumItofMatchDistanceCm >= 0.) ||
      !(config.maximumEtofMatchDistanceCm >= 0.) ||
      !(config.adapter.lambdaRoleDcaBoundaryCm > 0.) ||
      !(config.adapter.fallbackTrackSigmaXcm > 0.) ||
      !(config.adapter.fallbackTrackSigmaYcm > 0.) ||
      !(config.adapter.fallbackTrackSigmaZcm > 0.) ||
      !(config.adapter.fallbackDirectionSigmaRad > 0.) ||
      !(config.adapter.fallbackVertexSigmaXcm > 0.) ||
      !(config.adapter.fallbackVertexSigmaYcm > 0.) ||
      !(config.adapter.fallbackVertexSigmaZcm > 0.) ||
      !(config.maximumTrackChi2Ndf > 0.) ||
      !(config.finder.lambdaCuts.minimumMassGeV >= 0.) ||
      !(config.finder.lambdaCuts.minimumMassGeV <
        config.finder.lambdaCuts.maximumMassGeV) ||
      !(config.finder.lambdaCuts.maximumDaughterDcaCm >= 0.) ||
      !(config.finder.lambdaCuts.maximumGeometricChi2Ndf >= 0.) ||
      !(config.finder.lambdaCuts.maximumTopologicalChi2Ndf >= 0.) ||
      !(config.finder.lambdaCuts.minimumDecayLengthCm >= 0.) ||
      !(config.finder.lambdaCuts.minimumDecayLengthSignificance >= 0.) ||
      !(config.finder.lambdaCuts.minimumDaughterChi2FromPv >= 0.) ||
      !(config.finder.lambdaCuts.minimumDaughterDeviationFromPv >= 0.) ||
      config.finder.lambdaCuts.minimumPointingCosine < -1. ||
      config.finder.lambdaCuts.minimumPointingCosine > 1.) {
    std::cerr << "Invalid numeric option range; run --help for expected units."
              << std::endl;
    return 1;
  }
  if (config.t0.source == CeeV0Analysis::kT0Hit &&
      !(config.t0.rawToNsScale > 0.)) {
    std::cerr << "CeeT0Hit scale must be positive." << std::endl;
    return 1;
  }
  if (config.trackInputMode == CeeV0Lite::kCeeZeroFieldTofMode &&
      std::fabs(config.finder.magneticFieldByTesla) >
        config.finder.zeroFieldThresholdTesla) {
    std::cerr << "--momentum tof requires zero B; use --momentum tracking for "
                 "magnetic-field data." << std::endl;
    return 1;
  }
  if (config.trackInputMode == CeeV0Lite::kCeeMagneticFieldTrackingMode &&
      std::fabs(config.finder.magneticFieldByTesla) <=
        config.finder.zeroFieldThresholdTesla) {
    std::cerr << "--momentum tracking requires a non-zero --by-tesla."
              << std::endl;
    return 1;
  }

  CeeV0Analysis::Analyzer analyzer(config);
  return analyzer.Run();
}
