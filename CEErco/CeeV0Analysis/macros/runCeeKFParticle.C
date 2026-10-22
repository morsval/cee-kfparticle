// Run from CEErco/CeeV0Analysis after `make shared` in the CEE environment:
// root -l -b -q 'macros/runCeeKFParticle.C("input.list","mothers.root",0.5)'
// ByTesla=0 selects zero-field TOF; signed nonzero By selects tracking momentum.
#if defined(__CLING__)
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(../CeeV0Lite/include)
R__LOAD_LIBRARY(./libCeeKFParticleAnalysis.so)
#endif
#include "../include/CeeKFParticleAnalysisMaker.h"
#include "../include/CeeKFParticleRootSource.h"
#include "../include/CeeKFParticleRootWriter.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Edit this small class for your analysis; no input-file or ROOT-tree code here.
class UserCeeAnalysis : public CeeV0Analysis::CeeKFParticleAnalysisMaker {
 public:
  explicit UserCeeAnalysis(const CeeV0Lite::FinderConfig& c) : CeeKFParticleAnalysisMaker(c) {}
  unsigned IdentifyTrack(const CeeV0Analysis::AnalysisTrack& t,
                         const CeeV0Analysis::AnalysisEvent&) const override {
    // Default: retain adapter hypotheses; this is NOT a calibrated PID cut.
    // Replace with your dEdx/PID calibration. Example when nSigma is supplied:
    // unsigned mask=0;
    // if (std::isfinite(t.nSigmaProton) && std::fabs(t.nSigmaProton)<3.)
    //   mask |= CeeV0Lite::kProtonHypothesis;
    // if (std::isfinite(t.nSigmaPion) && std::fabs(t.nSigmaPion)<3.)
    //   mask |= CeeV0Lite::kPionHypothesis;
    // return mask; // 0 rejects the track; both bits retain PID ambiguity.
    return t.input.hypothesisMask;
  }
  bool AcceptMother(const CeeV0Lite::Candidate& c,
                    const CeeV0Analysis::AnalysisEvent&) const override {
    // Example loose analysis cut on original mother. No L/sigmaL cut by default.
    // Additional examples (choose thresholds using your analysis validation):
    // if (!(CeeV0Analysis::DcaProton(c)>0.1)) return false;
    // if (!(CeeV0Analysis::DcaPion(c)>0.1)) return false;
    // if (!c.motherTopo.valid || !(c.motherTopo.pvConstraintChi2/2.<5.)) return false;
    return c.massGeV>1.05 && c.massGeV<1.25;
  }
};

void runCeeKFParticle(const char* input="input.root",const char* output="mothers.root",
                     double ByTesla=0.,long long maximumEvents=-1,
                     bool realData=false,bool overwrite=false) {
  using namespace CeeV0Analysis;
  using namespace CeeV0Lite;
  AnalysisConfig config;
  ConfigureDataMode(config,realData ? kRealData:kSimulationData);
  config.inputPath=input; config.outputPath=output; config.maximumEvents=maximumEvents;
  config.enableSimulationTruth=false; // the Maker reader never uses truth for PID
  // Begin with broad baseline cuts. This avoids accidentally inheriting a
  // tight preset when the user wants to apply mother cuts in UserCeeAnalysis.
  if (!ApplyCutPreset("legacy",config)) throw std::runtime_error("cut preset missing");
  config.finder.lambdaCuts.minimumMassGeV=1.0;
  config.finder.lambdaCuts.maximumMassGeV=1.3;
  config.finder.lambdaCuts.maximumDaughterDcaCm=5.; // example loose baseline
  config.finder.lambdaCuts.maximumGeometricChi2Ndf=50.;
  config.finder.lambdaCuts.minimumDecayLengthSignificance=0.;
  config.finder.magneticFieldByTesla=ByTesla;
  const bool magnetic=std::fabs(ByTesla)>config.finder.zeroFieldThresholdTesla;
  config.trackInputMode=magnetic ? kCeeMagneticFieldTrackingMode:kCeeZeroFieldTofMode;
  config.finder.momentumMode=magnetic ? kTrackingMomentum:kTofMomentum;
  config.finder.useLambdaDcaRoleAssignment=!magnetic; // preserve zero-field legacy convention
  config.finder.usePrimaryVertexDcaReference=true;
  config.finder.buildProductionState=true;
  // If running real zero-field data, replace config.t0 / TOF scales, offsets
  // and errors with YOUR detector calibration before running this macro.
  CeeKFParticleRootSource source(config);
  std::string error;
  if (!source.Open(error)) throw std::runtime_error(error);
  UserCeeAnalysis maker(config.finder);
  MakerTrackCuts cuts; cuts.minimumHits=15; cuts.requireConverged=true;
  cuts.maximumChi2Ndf=10.; maker.SetTrackCuts(cuts); // example, configurable
  CeeKFParticleRootWriter writer(output);
  writer.SetOverwrite(overwrite);
  writer.SetVariables({"candidateType","trackId0","trackId1","massGeV","sigmaMassGeV",
    "ptGeV","rapidity","daughterDcaCm","DcaProton","DcaPion",
    "Chi2PrimaryProton","Chi2PrimaryPion","decayLengthCm","pointingCosine",
    "geometricChi2Ndf","motherTopoValid","motherTopoStatus","topoMassGeV",
    "motherTopoPvChi2Ndf","topoLCm","topoSigmaLCm","LOverSigmaL","flightS","sigmaS","SOverSigmaS"});
  // User-defined branch, without changing the default writer:
  writer.AddColumn("massMeV",[](const Candidate& c,const AnalysisEvent&){return c.massGeV*1000.;},"MeV/c^2");
  writer.SaveCovariances(false); // true adds 7x7, 8x8 and daughter 6x6 arrays
  std::ostringstream metadata;
  metadata<<"v0.7.0; input="<<input<<"; ByTesla="<<ByTesla
          <<"; baseline M=[1.0,1.3],daughterDCA<=5,geoChi2/NDF<=50"
          <<"; analysis M=(1.05,1.25); hits>=15,chi2/NDF<=10,converged"
          <<"; PID=adapter masks (replace with calibration); no flight-significance cut";
  writer.SetMetadata(metadata.str());
  if (!maker.Run(source,writer,error)) throw std::runtime_error(error);
  const MakerStatistics& n=maker.GetStatistics();
  std::cout<<"CEE v0.7.0: events="<<n.eventsAccepted<<", selected tracks="<<n.tracksSelected
           <<", mothers="<<n.mothersAccepted<<" -> "<<output<<std::endl;
}
