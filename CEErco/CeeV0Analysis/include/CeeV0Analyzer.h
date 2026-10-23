#ifndef CEE_V0_ANALYZER_H
#define CEE_V0_ANALYZER_H

#include "CeeV0AnalysisAlgorithms.h"
#include "CeeV0LiteCeeAdapter.h"
#include "CeeV0LiteFinder.h"

#include <string>

namespace CeeV0Analysis
{

enum DataMode
{
  kSimulationData = 0,
  kRealData = 1
};

struct AnalysisConfig
{
  AnalysisConfig();

  std::string inputPath;
  std::string outputPath;
  std::string treeName;
  DataMode dataMode;
  bool enableSimulationTruth;
  long long maximumEvents;
  int maximumInputFiles;
  int debugLevel;

  T0SelectionConfig t0;
  double tofRawToNsScale;
  double tofOffsetNs;
  double tofTimeSigmaNs;
  double tofPathLengthSigmaCm;
  double maximumItofMatchDistanceCm;
  double maximumEtofMatchDistanceCm;

  int minimumTpcHits;
  double maximumTrackChi2Ndf;
  std::size_t maximumTracksPerEvent;

  CeeV0Lite::CeeTrackInputMode trackInputMode;
  CeeV0Lite::CeeAdapterConfig adapter;
  CeeV0Lite::FinderConfig finder;
};

/**
 * Apply all mode-dependent timing/truth defaults in one place. Call this
 * before user overrides when building an AnalysisConfig outside the CLI.
 */
void ConfigureDataMode(AnalysisConfig& config, DataMode mode);

/** Apply one of: legacy, loose, star. Returns false for an unknown name. */
bool ApplyCutPreset(const std::string& name, AnalysisConfig& config);

/** True only when real-data analysis intentionally consumes CeeT0Hit. */
bool UsesMeasuredT0(const AnalysisConfig& config);

/**
 * Validate the simulation/real timing contract.  Simulation uses the exact
 * generator time origin t=0 and must not carry a measured/shared T0 term.
 */
bool ValidateDataModeTiming(const AnalysisConfig& config,
                            std::string* errorMessage = 0);

class Analyzer
{
 public:
  explicit Analyzer(const AnalysisConfig& config);
  ~Analyzer();

  /** Returns 0 on success and a non-zero value on setup or I/O failure. */
  int Run();

 private:
  Analyzer(const Analyzer&);
  Analyzer& operator=(const Analyzer&);

  struct Implementation;
  Implementation* fImplementation;
};

} // namespace CeeV0Analysis

#endif
