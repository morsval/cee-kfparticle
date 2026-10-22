// From CEErco/CeeV0Analysis, after top-level make:
// root -l -b -q
// 'macros/runCeeChannels.C("input.list","omega.root",0.5,"omega",true)' Last
// true enables CHARGE-ONLY demonstration, NOT calibrated PID.
#if defined(__CLING__)
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(../CeeV0Lite/include)
R__LOAD_LIBRARY(./libCeeKFParticleAnalysis.so)
#endif
#include "../include/CeeKFParticleAnalysisMaker.h"
#include "../include/CeeKFParticleChannelWriter.h"
#include "../include/CeeKFParticleRootSource.h"
#include <cmath>
#include <stdexcept>
using namespace CeeV0Analysis;

// Replace this callback with calibrated PID. It can use an external table
// keyed by eventId and trackId; returning true for several hypotheses is valid.
bool ChannelPID(const AnalysisTrack &t, const ParticleHypothesis &p,
                const AnalysisEvent &) {
  if (std::abs(p.pdg) == 2212)
    return std::isfinite(t.nSigmaProton) && std::fabs(t.nSigmaProton) < 3.;
  if (std::abs(p.pdg) == 211)
    return std::isfinite(t.nSigmaPion) && std::fabs(t.nSigmaPion) < 3.;
  if (std::abs(p.pdg) == 321)
    return std::isfinite(t.nSigmaKaon) && std::fabs(t.nSigmaKaon) < 3.;
  return false; // Add your He3 and deuteron PID here. Never invent a
                // calibration.
}
void runCeeChannels(const char *input = "input.list",
                    const char *output = "mothers.root", double ByTesla = .5,
                    const char *channel = "omega", bool chargeOnlyDemo = false,
                    long long maximumEvents = -1) {
  using namespace CeeV0Lite;
  AnalysisConfig cfg;
  ConfigureDataMode(cfg, kSimulationData);
  cfg.inputPath = input;
  cfg.enableSimulationTruth = false;
  cfg.trackInputMode = kCeeMagneticFieldTrackingMode;
  cfg.finder.momentumMode = kTrackingMomentum;
  cfg.finder.magneticFieldByTesla = ByTesla;
  CeeKFParticleAnalysisMaker maker(cfg.finder);
  MakerTrackCuts quality;
  quality.minimumHits = 15;
  quality.maximumChi2Ndf = 10.;
  quality.requireConverged = true;
  maker.SetTrackCuts(quality);
  maker.Channels().pid = [chargeOnlyDemo](const AnalysisTrack &t,
                                          const ParticleHypothesis &p,
                                          const AnalysisEvent &e) {
    return chargeOnlyDemo || ChannelPID(t, p, e);
  };
  DecayChannel lambda;
  lambda.name = "lambda";
  lambda.motherPdg = 3122;
  lambda.daughters = {DaughterSpec::Track("proton", {2212, 1, .9382720813}),
                      DaughterSpec::Track("pion", {-211, -1, .13957039})};
  // This callback is the baseline selection BEFORE Lambda enters Omega pairs.
  lambda.select = [](const ChannelCandidate &c, const AnalysisEvent &) {
    return c.massGeV > 1.08 && c.massGeV < 1.16 && c.chi2 / c.ndf < 20.;
  };
  DecayChannel omega;
  omega.name = "omega";
  omega.motherPdg = 3334;
  omega.fitCollisionPv = true;
  omega.daughters = {DaughterSpec::Composite("lambda", "lambda"),
                     DaughterSpec::Track("kaon", {-321, -1, .493677})};
  omega.select = [](const ChannelCandidate &c, const AnalysisEvent &) {
    return c.massGeV > 1.5 && c.massGeV < 1.9 && c.chi2 / c.ndf < 20.;
  };
  DecayChannel h2;
  h2.name = "hypertriton2";
  h2.motherPdg = 1010010030;
  h2.fitCollisionPv = true;
  h2.daughters = {DaughterSpec::Track("he3", {1000020030, 2, 2.8083916}),
                  DaughterSpec::Track("pion", {-211, -1, .13957039})};
  DecayChannel h3;
  h3.name = "hypertriton3";
  h3.motherPdg = 1010010030;
  h3.fitCollisionPv = true;
  h3.daughters = {DaughterSpec::Track("deuteron", {1000010020, 1, 1.8756129}),
                  DaughterSpec::Track("proton", {2212, 1, .9382720813}),
                  DaughterSpec::Track("pion", {-211, -1, .13957039})};
  std::string choice(channel);
  if (choice == "omega")
    maker.Channels().channels = {lambda, omega};
  else if (choice == "lambda")
    maker.Channels().channels = {lambda};
  else if (choice == "hypertriton2")
    maker.Channels().channels = {h2};
  else if (choice == "hypertriton3")
    maker.Channels().channels = {h3};
  else
    throw std::runtime_error("unknown channel");
  CeeKFParticleRootSource source(cfg);
  std::string error;
  if (!source.Open(error))
    throw std::runtime_error(error);
  CeeKFParticleChannelWriter writer(output);
  writer.SaveCovariances(true);
  writer.AddColumn("pvChi2Ndf", [](const ChannelCandidate &c,
                                   const AnalysisEvent &) {
    return c.motherTopo.valid && c.motherTopo.pvConstraintNdf > 0
               ? c.motherTopo.pvConstraintChi2 / c.motherTopo.pvConstraintNdf
               : std::numeric_limits<double>::quiet_NaN();
  });
  if (!maker.RunChannels(source, writer, error, maximumEvents))
    throw std::runtime_error(error);
}
