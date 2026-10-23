#include "CeeV0Analyzer.h"

#include "CeeMCTrack.h"
#include "CeeT0Hit.h"
#include "CeeEtofHit.h"
#include "CeeTpcTrack.h"
#include "CeeTpcVertex.h"
#include "CeeiTOFHit.h"

#include "TBranch.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TList.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace CeeV0Analysis
{

namespace
{

const unsigned int kLambdaDaughterHypotheses =
  CeeV0Lite::kPionHypothesis | CeeV0Lite::kProtonHypothesis;

std::string Trim(const std::string& input)
{
  const std::string whitespace = " \t\r\n";
  const std::size_t first = input.find_first_not_of(whitespace);
  if (first == std::string::npos) return std::string();
  const std::size_t last = input.find_last_not_of(whitespace);
  return input.substr(first, last - first + 1u);
}

bool EndsWith(const std::string& value, const std::string& suffix)
{
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<std::string> ResolveInputFiles(const std::string& input,
                                           int maximumFiles)
{
  std::vector<std::string> files;
  if (EndsWith(input, ".root")) {
    files.push_back(input);
  }
  else if (EndsWith(input, ".list") || EndsWith(input, ".txt")) {
    std::ifstream stream(input.c_str());
    std::string line;
    while (std::getline(stream, line)) {
      line = Trim(line);
      if (line.empty() || line[0] == '#') continue;
      if (EndsWith(line, ".root")) files.push_back(line);
    }
  }
  else {
    TSystemDirectory directory("ceeV0Input", input.c_str());
    TList* entries = directory.GetListOfFiles();
    if (entries) {
      TIter next(entries);
      while (TSystemFile* entry = static_cast<TSystemFile*>(next())) {
        const std::string name = entry->GetName();
        if (!entry->IsDirectory() && EndsWith(name, ".root"))
          files.push_back(input + "/" + name);
      }
    }
    std::sort(files.begin(), files.end());
  }
  if (maximumFiles > 0 && files.size() > static_cast<std::size_t>(maximumFiles))
    files.resize(static_cast<std::size_t>(maximumFiles));
  return files;
}

ProjectionTrack MakeProjectionTrack(const CeeV0Lite::CeeTrackReferenceState& reference,
                                    int id)
{
  ProjectionTrack output;
  output.id = id;
  for (int coordinate = 0; coordinate < 3; ++coordinate) {
    output.positionCm[coordinate] = reference.positionCm[coordinate];
    output.direction[coordinate] = reference.direction[coordinate];
  }
  return output;
}

double TrackChi2Ndf(const CeeTpcTrack& track)
{
  const double direct = track.GetFitChi2Ndf();
  if (std::isfinite(direct) && direct > 0.) return direct;
  const double ndf = track.GetFitNdf();
  const double chi2 = track.GetFitChi2();
  if (std::isfinite(ndf) && ndf > 0. && std::isfinite(chi2) && chi2 >= 0.)
    return chi2 / ndf;
  return std::numeric_limits<double>::infinity();
}

struct SelectedTrack
{
  SelectedTrack() : track(0), originalIndex(-1), usedVertexSeed(false) {}
  CeeTpcTrack* track;
  int originalIndex;
  bool usedVertexSeed;
};

struct MatchMetadata
{
  MatchMetadata()
    : detector(kNoTofDetector), hitId(-1), distanceCm(-1.),
      rawTime(std::numeric_limits<double>::quiet_NaN()),
      timeNs(std::numeric_limits<double>::quiet_NaN()) {}
  TofDetector detector;
  int hitId;
  double distanceCm;
  double rawTime;
  double timeNs;
  CeeV0Lite::Vec3 positionCm;
};

struct TruthResult
{
  TruthResult()
    : available(-1), matched(-1), sameMotherLambda(-1), roleCorrect(-1),
      motherId(-1), motherPdg(0),
      daughterMcId0(-1), daughterMcId1(-1), daughterMcPdg0(0),
      daughterMcPdg1(0) {}
  int available;
  int matched;
  int sameMotherLambda;
  int roleCorrect;
  int motherId;
  int motherPdg;
  int daughterMcId0;
  int daughterMcId1;
  int daughterMcPdg0;
  int daughterMcPdg1;
};

struct CandidateRow
{
  CandidateRow() { Reset(); }

  void Reset()
  {
    event = -1;
    candidateIndex = -1;
    candidateType = 0;
    trackId0 = trackId1 = -1;
    daughterPdg0 = daughterPdg1 = 0;
    nHitsFit0 = nHitsFit1 = -1;
    fitChi2Ndf0 = fitChi2Ndf1 = std::numeric_limits<double>::quiet_NaN();
    dEdx0 = dEdx1 = std::numeric_limits<double>::quiet_NaN();
    fitCharge0 = fitCharge1 = std::numeric_limits<double>::quiet_NaN();
    usedVertexSeed0 = usedVertexSeed1 = -1;
    tofDetector0 = tofDetector1 = 0;
    tofHitId0 = tofHitId1 = -1;
    tofMatchDistance0Cm = tofMatchDistance1Cm = -1.;
    tofRawTime0 = tofRawTime1 = std::numeric_limits<double>::quiet_NaN();
    tofTime0Ns = tofTime1Ns = std::numeric_limits<double>::quiet_NaN();
    startTimeRaw = std::numeric_limits<double>::quiet_NaN();
    startTimeNs = 0.;
    sigmaStartTimeNs = 0.;
    numberT0Hits = numberInputTracks = numberSelectedTracks = 0;
    numberFinderTracks = numberTofMatches = numberCandidates = 0;
    pvX = pvY = pvZ = 0.;
    pvSigmaX = pvSigmaY = pvSigmaZ = 0.;
    pvUsedStoredErrors = -1;
    svX = svY = svZ = 0.;
    svSigmaX = svSigmaY = svSigmaZ = 0.;
    px = py = pz = 0.;
    energyGeV = 0.;
    massGeV = sigmaMassGeV = 0.;
    motherCovarianceValid = 0;
    motherCovarianceStatus = 0;
    motherJacobianComplete = 0;
    motherJacobianPerturbationFailed = 0;
    motherJacobianFitsAttempted = 0;
    motherJacobianFitsSucceeded = 0;
    motherJacobianFailureCount = 0;
    motherJacobianPlusFailureMask = 0;
    motherJacobianMinusFailureMask = 0;
    motherJacobianFailedColumn = -1;
    motherJacobianFailedSide = 0;
    motherJacobianFirstFailureFitStatus = 0;
    covarianceDependentCutsEvaluated = 0;
    for (int i = 0; i < 28; ++i) motherCovariance[i] = 0.;
    motherTopoValid = 0;
    motherTopoStatus = 0;
    motherTopoCorrelationModel = 0;
    motherTopoAtProductionVertex = 0;
    motherTopoIterations = 0;
    motherTopoCharge = 0;
    motherTopoByTesla = std::numeric_limits<double>::quiet_NaN();
    motherTopoPvNdf = 0;
    motherTopoNdf = 0;
    motherTopoFlightValid = 0;
    motherTopoS = std::numeric_limits<double>::quiet_NaN();
    motherTopoSigmaS = std::numeric_limits<double>::quiet_NaN();
    motherTopoPvChi2 = std::numeric_limits<double>::quiet_NaN();
    motherTopoChi2 = std::numeric_limits<double>::quiet_NaN();
    motherTopoChi2Ndf = std::numeric_limits<double>::quiet_NaN();
    motherTopoMassGeV = std::numeric_limits<double>::quiet_NaN();
    motherTopoSigmaMassGeV = std::numeric_limits<double>::quiet_NaN();
    motherTopoDecayLengthCm = std::numeric_limits<double>::quiet_NaN();
    motherTopoSigmaDecayLengthCm = std::numeric_limits<double>::quiet_NaN();
    motherTopoDecayLengthSignificance = std::numeric_limits<double>::quiet_NaN();
    motherTopoCtauCm = std::numeric_limits<double>::quiet_NaN();
    motherTopoSigmaCtauCm = std::numeric_limits<double>::quiet_NaN();
    for (int i=0; i<8; ++i) motherTopoParameters[i] = std::numeric_limits<double>::quiet_NaN();
    for (int i=0; i<36; ++i) motherTopoCovariance[i] = std::numeric_limits<double>::quiet_NaN();
    daughterDcaCm = geometricChi2 = geometricChi2Ndf = 0.;
    geometricNdf = 0;
    topologicalChi2 = topologicalChi2Ndf = 0.;
    topologicalNdf = 0;
    primaryVertexConstraintChi2 = primaryVertexConstraintChi2Ndf = 0.;
    primaryVertexConstraintNdf = 0;
    decayLengthCm = sigmaDecayLengthCm = decayLengthSignificance = 0.;
    pointingCosine = 0.;
    daughterDcaToPv0Cm = daughterDcaToPv1Cm = 0.;
    daughterPvChi20 = daughterPvChi21 = 0.;
    daughterPvSignificance0 = daughterPvSignificance1 = 0.;
    decayTimeNs = 0.;
    tofIterations = 0;
    tofConverged = 0;
    dataMode = 0;
    truth = TruthResult();
  }

  int event;
  int candidateIndex;
  int candidateType;
  int trackId0;
  int trackId1;
  int daughterPdg0;
  int daughterPdg1;
  int nHitsFit0;
  int nHitsFit1;
  double fitChi2Ndf0;
  double fitChi2Ndf1;
  double dEdx0;
  double dEdx1;
  double fitCharge0;
  double fitCharge1;
  int usedVertexSeed0;
  int usedVertexSeed1;
  int tofDetector0;
  int tofDetector1;
  int tofHitId0;
  int tofHitId1;
  double tofMatchDistance0Cm;
  double tofMatchDistance1Cm;
  double tofRawTime0;
  double tofRawTime1;
  double tofTime0Ns;
  double tofTime1Ns;
  double startTimeRaw;
  double startTimeNs;
  double sigmaStartTimeNs;
  int numberT0Hits;
  int numberInputTracks;
  int numberSelectedTracks;
  int numberFinderTracks;
  int numberTofMatches;
  int numberCandidates;
  double pvX;
  double pvY;
  double pvZ;
  double pvSigmaX;
  double pvSigmaY;
  double pvSigmaZ;
  int pvUsedStoredErrors;
  double svX;
  double svY;
  double svZ;
  double svSigmaX;
  double svSigmaY;
  double svSigmaZ;
  double px;
  double py;
  double pz;
  double energyGeV;
  double massGeV;
  double sigmaMassGeV;
  int motherCovarianceValid;
  int motherCovarianceStatus;
  int motherJacobianComplete;
  int motherJacobianPerturbationFailed;
  int motherJacobianFitsAttempted;
  int motherJacobianFitsSucceeded;
  int motherJacobianFailureCount;
  int motherJacobianPlusFailureMask;
  int motherJacobianMinusFailureMask;
  int motherJacobianFailedColumn;
  int motherJacobianFailedSide;
  int motherJacobianFirstFailureFitStatus;
  int covarianceDependentCutsEvaluated;
  double motherCovariance[28];
  int motherTopoValid;
  int motherTopoStatus;
  int motherTopoCorrelationModel;
  int motherTopoAtProductionVertex;
  int motherTopoIterations;
  int motherTopoCharge;
  double motherTopoByTesla;
  int motherTopoPvNdf;
  int motherTopoNdf;
  int motherTopoFlightValid;
  double motherTopoS;
  double motherTopoSigmaS;
  double motherTopoPvChi2;
  double motherTopoChi2;
  double motherTopoChi2Ndf;
  double motherTopoMassGeV;
  double motherTopoSigmaMassGeV;
  double motherTopoDecayLengthCm;
  double motherTopoSigmaDecayLengthCm;
  double motherTopoDecayLengthSignificance;
  double motherTopoCtauCm;
  double motherTopoSigmaCtauCm;
  double motherTopoParameters[8];
  double motherTopoCovariance[36];
  double daughterDcaCm;
  double geometricChi2;
  int geometricNdf;
  double geometricChi2Ndf;
  double topologicalChi2;
  int topologicalNdf;
  double topologicalChi2Ndf;
  double primaryVertexConstraintChi2;
  int primaryVertexConstraintNdf;
  double primaryVertexConstraintChi2Ndf;
  double decayLengthCm;
  double sigmaDecayLengthCm;
  double decayLengthSignificance;
  double pointingCosine;
  double daughterDcaToPv0Cm;
  double daughterDcaToPv1Cm;
  double daughterPvChi20;
  double daughterPvChi21;
  double daughterPvSignificance0;
  double daughterPvSignificance1;
  double decayTimeNs;
  int tofIterations;
  int tofConverged;
  int dataMode;
  TruthResult truth;
};

void CreateCandidateBranches(TTree& tree, CandidateRow& row)
{
#define CEE_BRANCH_INT(name) tree.Branch(#name, &row.name, #name "/I")
#define CEE_BRANCH_DOUBLE(name) tree.Branch(#name, &row.name, #name "/D")
  CEE_BRANCH_INT(event);
  CEE_BRANCH_INT(candidateIndex);
  CEE_BRANCH_INT(candidateType);
  CEE_BRANCH_INT(trackId0);
  CEE_BRANCH_INT(trackId1);
  CEE_BRANCH_INT(daughterPdg0);
  CEE_BRANCH_INT(daughterPdg1);
  CEE_BRANCH_INT(nHitsFit0);
  CEE_BRANCH_INT(nHitsFit1);
  CEE_BRANCH_DOUBLE(fitChi2Ndf0);
  CEE_BRANCH_DOUBLE(fitChi2Ndf1);
  CEE_BRANCH_DOUBLE(dEdx0);
  CEE_BRANCH_DOUBLE(dEdx1);
  CEE_BRANCH_DOUBLE(fitCharge0);
  CEE_BRANCH_DOUBLE(fitCharge1);
  CEE_BRANCH_INT(usedVertexSeed0);
  CEE_BRANCH_INT(usedVertexSeed1);
  CEE_BRANCH_INT(tofDetector0);
  CEE_BRANCH_INT(tofDetector1);
  CEE_BRANCH_INT(tofHitId0);
  CEE_BRANCH_INT(tofHitId1);
  CEE_BRANCH_DOUBLE(tofMatchDistance0Cm);
  CEE_BRANCH_DOUBLE(tofMatchDistance1Cm);
  CEE_BRANCH_DOUBLE(tofRawTime0);
  CEE_BRANCH_DOUBLE(tofRawTime1);
  CEE_BRANCH_DOUBLE(tofTime0Ns);
  CEE_BRANCH_DOUBLE(tofTime1Ns);
  CEE_BRANCH_DOUBLE(startTimeRaw);
  CEE_BRANCH_DOUBLE(startTimeNs);
  CEE_BRANCH_DOUBLE(sigmaStartTimeNs);
  CEE_BRANCH_INT(numberT0Hits);
  CEE_BRANCH_INT(numberInputTracks);
  CEE_BRANCH_INT(numberSelectedTracks);
  CEE_BRANCH_INT(numberFinderTracks);
  CEE_BRANCH_INT(numberTofMatches);
  CEE_BRANCH_INT(numberCandidates);
  CEE_BRANCH_DOUBLE(pvX);
  CEE_BRANCH_DOUBLE(pvY);
  CEE_BRANCH_DOUBLE(pvZ);
  CEE_BRANCH_DOUBLE(pvSigmaX);
  CEE_BRANCH_DOUBLE(pvSigmaY);
  CEE_BRANCH_DOUBLE(pvSigmaZ);
  CEE_BRANCH_INT(pvUsedStoredErrors);
  CEE_BRANCH_DOUBLE(svX);
  CEE_BRANCH_DOUBLE(svY);
  CEE_BRANCH_DOUBLE(svZ);
  CEE_BRANCH_DOUBLE(svSigmaX);
  CEE_BRANCH_DOUBLE(svSigmaY);
  CEE_BRANCH_DOUBLE(svSigmaZ);
  CEE_BRANCH_DOUBLE(px);
  CEE_BRANCH_DOUBLE(py);
  CEE_BRANCH_DOUBLE(pz);
  CEE_BRANCH_DOUBLE(energyGeV);
  CEE_BRANCH_DOUBLE(massGeV);
  CEE_BRANCH_DOUBLE(sigmaMassGeV);
  CEE_BRANCH_INT(motherCovarianceValid);
  CEE_BRANCH_INT(motherCovarianceStatus);
  CEE_BRANCH_INT(motherJacobianComplete);
  CEE_BRANCH_INT(motherJacobianPerturbationFailed);
  CEE_BRANCH_INT(motherJacobianFitsAttempted);
  CEE_BRANCH_INT(motherJacobianFitsSucceeded);
  CEE_BRANCH_INT(motherJacobianFailureCount);
  CEE_BRANCH_INT(motherJacobianPlusFailureMask);
  CEE_BRANCH_INT(motherJacobianMinusFailureMask);
  CEE_BRANCH_INT(motherJacobianFailedColumn);
  CEE_BRANCH_INT(motherJacobianFailedSide);
  CEE_BRANCH_INT(motherJacobianFirstFailureFitStatus);
  CEE_BRANCH_INT(covarianceDependentCutsEvaluated);
  tree.Branch("motherCovariance", row.motherCovariance,
              "motherCovariance[28]/D");
  CEE_BRANCH_INT(motherTopoValid);
  CEE_BRANCH_INT(motherTopoStatus);
  CEE_BRANCH_INT(motherTopoCorrelationModel);
  CEE_BRANCH_INT(motherTopoAtProductionVertex);
  CEE_BRANCH_INT(motherTopoIterations);
  CEE_BRANCH_INT(motherTopoCharge);
  CEE_BRANCH_DOUBLE(motherTopoByTesla);
  CEE_BRANCH_INT(motherTopoPvNdf);
  CEE_BRANCH_INT(motherTopoNdf);
  CEE_BRANCH_INT(motherTopoFlightValid);
  CEE_BRANCH_DOUBLE(motherTopoS);
  CEE_BRANCH_DOUBLE(motherTopoSigmaS);
  CEE_BRANCH_DOUBLE(motherTopoPvChi2);
  CEE_BRANCH_DOUBLE(motherTopoChi2);
  CEE_BRANCH_DOUBLE(motherTopoChi2Ndf);
  CEE_BRANCH_DOUBLE(motherTopoMassGeV);
  CEE_BRANCH_DOUBLE(motherTopoSigmaMassGeV);
  CEE_BRANCH_DOUBLE(motherTopoDecayLengthCm);
  CEE_BRANCH_DOUBLE(motherTopoSigmaDecayLengthCm);
  CEE_BRANCH_DOUBLE(motherTopoDecayLengthSignificance);
  CEE_BRANCH_DOUBLE(motherTopoCtauCm);
  CEE_BRANCH_DOUBLE(motherTopoSigmaCtauCm);
  tree.Branch("motherTopoParameters", row.motherTopoParameters, "motherTopoParameters[8]/D");
  tree.Branch("motherTopoCovariance", row.motherTopoCovariance, "motherTopoCovariance[36]/D");
  CEE_BRANCH_DOUBLE(daughterDcaCm);
  CEE_BRANCH_DOUBLE(geometricChi2);
  CEE_BRANCH_INT(geometricNdf);
  CEE_BRANCH_DOUBLE(geometricChi2Ndf);
  CEE_BRANCH_DOUBLE(topologicalChi2);
  CEE_BRANCH_INT(topologicalNdf);
  CEE_BRANCH_DOUBLE(topologicalChi2Ndf);
  CEE_BRANCH_DOUBLE(primaryVertexConstraintChi2);
  CEE_BRANCH_INT(primaryVertexConstraintNdf);
  CEE_BRANCH_DOUBLE(primaryVertexConstraintChi2Ndf);
  CEE_BRANCH_DOUBLE(decayLengthCm);
  CEE_BRANCH_DOUBLE(sigmaDecayLengthCm);
  CEE_BRANCH_DOUBLE(decayLengthSignificance);
  CEE_BRANCH_DOUBLE(pointingCosine);
  CEE_BRANCH_DOUBLE(daughterDcaToPv0Cm);
  CEE_BRANCH_DOUBLE(daughterDcaToPv1Cm);
  CEE_BRANCH_DOUBLE(daughterPvChi20);
  CEE_BRANCH_DOUBLE(daughterPvChi21);
  CEE_BRANCH_DOUBLE(daughterPvSignificance0);
  CEE_BRANCH_DOUBLE(daughterPvSignificance1);
  // 0.5.x ROOT-schema compatibility.  These aliases intentionally retain
  // the historical sqrt(chi2) meaning; new analyses should use the explicit
  // daughterPvSignificance* names (or daughterPvChi2* for STAR-style cuts).
  tree.Branch("daughterDeviation0", &row.daughterPvSignificance0,
              "daughterDeviation0/D");
  tree.Branch("daughterDeviation1", &row.daughterPvSignificance1,
              "daughterDeviation1/D");
  CEE_BRANCH_DOUBLE(decayTimeNs);
  CEE_BRANCH_INT(tofIterations);
  CEE_BRANCH_INT(tofConverged);
  CEE_BRANCH_INT(dataMode);
  tree.Branch("truthAvailable", &row.truth.available, "truthAvailable/I");
  tree.Branch("truthMatched", &row.truth.matched, "truthMatched/I");
  tree.Branch("truthSameMotherLambda", &row.truth.sameMotherLambda,
              "truthSameMotherLambda/I");
  tree.Branch("truthRoleCorrect", &row.truth.roleCorrect,
              "truthRoleCorrect/I");
  tree.Branch("truthMotherId", &row.truth.motherId, "truthMotherId/I");
  tree.Branch("truthMotherPdg", &row.truth.motherPdg, "truthMotherPdg/I");
  tree.Branch("truthDaughterMcId0", &row.truth.daughterMcId0,
              "truthDaughterMcId0/I");
  tree.Branch("truthDaughterMcId1", &row.truth.daughterMcId1,
              "truthDaughterMcId1/I");
  tree.Branch("truthDaughterMcPdg0", &row.truth.daughterMcPdg0,
              "truthDaughterMcPdg0/I");
  tree.Branch("truthDaughterMcPdg1", &row.truth.daughterMcPdg1,
              "truthDaughterMcPdg1/I");
#undef CEE_BRANCH_INT
#undef CEE_BRANCH_DOUBLE
}

} // namespace

struct Analyzer::Implementation
{
  explicit Implementation(const AnalysisConfig& input)
    : config(input), chain(0), tracks(0), vertices(0), itofHits(0),
      etofHits(0), t0Hits(0), mcTracks(0), truthBranchBound(false),
      output(0), candidateTree(0), hEventCounter(0), hInputTracks(0),
      hSelectedTracks(0), hTofMatches(0), hTofMatchDistance(0), hT0(0),
      hT0Raw(0), hTofMinusT0(0), hTrackReferenceSource(0),
      hReferenceStatus(0), hAdapterStatus(0), hPvStoredErrors(0),
      hPvDcaReferenceStatus(0), hPvDcaReferencePath(0), hPvDcaReferenceDca(0),
      hGeneratedLambdaCount(0), hGeneratedLambdaPt(0),
      hGeneratedLambdaRapidity(0),
      hMass(0), hMassTruth(0), hMassBackground(0), hDaughterDca(0),
      hDecayLength(0), hPointingCosine(0), hGeometricChi2(0),
      hTopologicalChi2(0), hPrimaryVertexConstraintChi2(0),
      hDaughterPvChi2(0), finder(config.finder)
  {
  }

  ~Implementation()
  {
    delete output;
    delete chain;
    delete tracks;
    delete vertices;
    delete itofHits;
    delete etofHits;
    delete t0Hits;
    delete mcTracks;
  }

  bool HasBranch(const char* name) const
  {
    return chain && chain->GetBranch(name) != 0;
  }

  bool SetupInput()
  {
    const std::vector<std::string> files = ResolveInputFiles(
      config.inputPath, config.maximumInputFiles);
    if (files.empty()) {
      std::cerr << "CeeV0Analysis: no ROOT input files found at "
                << config.inputPath << std::endl;
      return false;
    }

    chain = new TChain(config.treeName.c_str());
    for (std::size_t i = 0; i < files.size(); ++i) {
      if (chain->Add(files[i].c_str()) == 0) {
        std::cerr << "CeeV0Analysis: failed to add " << files[i] << std::endl;
        return false;
      }
    }
    if (chain->GetEntries() < 0) return false;

    const char* required[] = {"CeeTpcTrack", "CeeTpcVertex"};
    for (std::size_t i = 0; i < sizeof(required) / sizeof(required[0]); ++i) {
      if (!HasBranch(required[i])) {
        std::cerr << "CeeV0Analysis: required branch " << required[i]
                  << " is missing" << std::endl;
        return false;
      }
    }
    const bool hasItof = HasBranch("CeeiTOFHit");
    const bool hasEtof = HasBranch("CeeEtofHit");
    if (config.trackInputMode == CeeV0Lite::kCeeZeroFieldTofMode &&
        !hasItof && !hasEtof) {
      std::cerr << "CeeV0Analysis: zero-field TOF mode requires at least one "
                   "of CeeiTOFHit/CeeEtofHit" << std::endl;
      return false;
    }
    const bool useMeasuredT0 = UsesMeasuredT0(config);
    if (useMeasuredT0 && !HasBranch("CeeT0Hit")) {
      std::cerr << "CeeV0Analysis: --t0-source hit requires CeeT0Hit"
                << std::endl;
      return false;
    }

    chain->SetBranchStatus("*", 0);
    tracks = new TClonesArray("CeeTpcTrack");
    vertices = new TClonesArray("CeeTpcVertex");
    chain->SetBranchStatus("CeeTpcTrack", 1);
    chain->SetBranchStatus("CeeTpcTrack*", 1);
    chain->SetBranchStatus("CeeTpcVertex", 1);
    chain->SetBranchStatus("CeeTpcVertex*", 1);
    chain->SetBranchAddress("CeeTpcTrack", &tracks);
    chain->SetBranchAddress("CeeTpcVertex", &vertices);
    if (hasItof &&
        config.trackInputMode == CeeV0Lite::kCeeZeroFieldTofMode) {
      itofHits = new TClonesArray("CeeiTOFHit");
      chain->SetBranchStatus("CeeiTOFHit", 1);
      chain->SetBranchStatus("CeeiTOFHit*", 1);
      chain->SetBranchAddress("CeeiTOFHit", &itofHits);
    }
    if (hasEtof &&
        config.trackInputMode == CeeV0Lite::kCeeZeroFieldTofMode) {
      etofHits = new TClonesArray("CeeEtofHit");
      chain->SetBranchStatus("CeeEtofHit", 1);
      chain->SetBranchStatus("CeeEtofHit*", 1);
      chain->SetBranchAddress("CeeEtofHit", &etofHits);
    }

    if (useMeasuredT0) {
      t0Hits = new TClonesArray("CeeT0Hit");
      chain->SetBranchStatus("CeeT0Hit", 1);
      chain->SetBranchStatus("CeeT0Hit*", 1);
      chain->SetBranchAddress("CeeT0Hit", &t0Hits);
    }

    // This is the only place that enables or binds MCTrack. In real mode this
    // block is never entered, so reconstruction cannot accidentally depend on
    // truth branches or FairLinks.
    if (config.dataMode == kSimulationData && config.enableSimulationTruth) {
      if (HasBranch("MCTrack")) {
        mcTracks = new TClonesArray("CeeMCTrack");
        chain->SetBranchStatus("MCTrack", 1);
        chain->SetBranchStatus("MCTrack*", 1);
        chain->SetBranchAddress("MCTrack", &mcTracks);
        truthBranchBound = true;
      }
      else {
        std::cerr << "CeeV0Analysis: MCTrack is absent; continuing without truth"
                  << std::endl;
      }
    }

    std::cout << "CeeV0Analysis: added " << files.size() << " file(s), "
              << chain->GetEntries() << " event(s)" << std::endl;
    return true;
  }

  bool SetupOutput()
  {
    output = TFile::Open(config.outputPath.c_str(), "RECREATE");
    if (!output || output->IsZombie()) {
      std::cerr << "CeeV0Analysis: cannot create " << config.outputPath
                << std::endl;
      return false;
    }
    output->cd();
    candidateTree = new TTree("LambdaCandidates",
                              "CEE lightweight Lambda candidates");
    CreateCandidateBranches(*candidateTree, row);
    hEventCounter = new TH1D("hEventCounter", "Event processing;stage;events",
                             7, 0.5, 7.5);
    const char* stages[] = {"all", "valid PV", "valid time origin", "selected tracks",
                            "momentum input ready", "finder run", "has Lambda"};
    for (int i = 0; i < 7; ++i) hEventCounter->GetXaxis()->SetBinLabel(i + 1, stages[i]);
    hInputTracks = new TH1D("hInputTracks", "TPC tracks;N_{track};events", 200, -0.5, 199.5);
    hSelectedTracks = new TH1D("hSelectedTracks", "Selected TPC tracks;N_{track};events", 200, -0.5, 199.5);
    hTofMatches = new TH1D("hTofMatches", "Unique TOF matches;N_{match};events", 200, -0.5, 199.5);
    hTofMatchDistance = new TH1D("hTofMatchDistance", "Track-TOF match distance;distance (cm);matches", 200, 0., 20.);
    hT0 = new TH1D("hT0", "Selected event start time;t_{0} (ns);events", 400, -20., 20.);
    hT0Raw = new TH1D("hT0Raw", "Selected raw T0Hit time;raw T0Hit time;events", 1000, -20000., 20000.);
    hTofMinusT0 = new TH1D("hTofMinusT0", "Converted TOF minus event T0;t_{TOF}-t_{0} (ns);matches", 600, -20., 100.);
    hTrackReferenceSource = new TH1D("hTrackReferenceSource", "Track reference state;source;tracks", 3, 0.5, 3.5);
    hTrackReferenceSource->GetXaxis()->SetBinLabel(1, "vertex seed");
    hTrackReferenceSource->GetXaxis()->SetBinLabel(2, "legacy fallback");
    hTrackReferenceSource->GetXaxis()->SetBinLabel(3, "rejected");
    hReferenceStatus = new TH1D("hReferenceAdapterStatus", "Reference-state adapter status;status;tracks", 14, -0.5, 13.5);
    hAdapterStatus = new TH1D("hTrackAdapterStatus", "Final track adapter status;status;tracks", 14, -0.5, 13.5);
    for (int status = 0; status < 14; ++status) {
      const char* label = CeeV0Lite::CeeAdapterStatusMessage(
        static_cast<CeeV0Lite::CeeAdapterStatus>(status));
      hReferenceStatus->GetXaxis()->SetBinLabel(status + 1, label);
      hAdapterStatus->GetXaxis()->SetBinLabel(status + 1, label);
    }
    hPvStoredErrors = new TH1D("hPvStoredErrors", "Primary-vertex stored errors;state;events", 2, -0.5, 1.5);
    hPvStoredErrors->GetXaxis()->SetBinLabel(1, "fallback/missing");
    hPvStoredErrors->GetXaxis()->SetBinLabel(2, "stored errors");
    hPvDcaReferenceStatus = new TH1D("hPvDcaReferenceStatus",
      "Magnetic track reference;status;tracks", 4, -0.5, 3.5);
    const char* referenceLabels[] = {"unchanged", "PV DCA", "invalid input", "transport failed"};
    for (int i = 0; i < 4; ++i)
      hPvDcaReferenceStatus->GetXaxis()->SetBinLabel(i+1, referenceLabels[i]);
    const double referenceLimit = config.finder.referenceTransportConfig.maximumAbsPathCm;
    hPvDcaReferencePath = new TH1D("hPvDcaReferencePath",
      "Reference relocation;signed arc (cm);tracks", 400, -referenceLimit, referenceLimit);
    hPvDcaReferenceDca = new TH1D("hPvDcaReferenceDca",
      "Track to PV DCA at new reference;DCA (cm);tracks", 400, 0., 100.);
    hGeneratedLambdaCount = new TH1D(
      "hGeneratedLambdaCount", "Generated MC #Lambda per event;N_{#Lambda};events",
      50, -0.5, 49.5);
    hGeneratedLambdaPt = new TH1D(
      "hGeneratedLambdaPt", "Generated MC #Lambda;p_{T} (GeV/c);particles",
      250, 0., 5.);
    hGeneratedLambdaRapidity = new TH1D(
      "hGeneratedLambdaRapidity", "Generated MC #Lambda;rapidity;particles",
      300, -3., 3.);
    hMass = new TH1D("hLambdaMass", "All #Lambda candidates;inv. mass (GeV/c^{2});candidates", 520, 0.7, 2.0);
    hMassTruth = new TH1D("hLambdaMassTruth", "Truth-matched #Lambda;inv. mass (GeV/c^{2});candidates", 520, 0.7, 2.0);
    hMassBackground = new TH1D("hLambdaMassBackground", "Unmatched #Lambda;inv. mass (GeV/c^{2});candidates", 520, 0.7, 2.0);
    hDaughterDca = new TH1D("hDaughterDca", "Daughter-pair DCA;DCA (cm);candidates", 200, 0., 10.);
    hDecayLength = new TH1D("hDecayLength", "#Lambda decay length;L (cm);candidates", 300, 0., 60.);
    hPointingCosine = new TH1D("hPointingCosine", "#Lambda pointing cosine;cos(#theta);candidates", 220, -1.1, 1.1);
    hGeometricChi2 = new TH1D("hGeometricChi2Ndf", "Geometric #chi^{2}/ndf;#chi^{2}/ndf;candidates", 250, 0., 50.);
    hTopologicalChi2 = new TH1D("hTopologicalChi2Ndf", "Topological #chi^{2}/ndf;#chi^{2}/ndf;candidates", 250, 0., 50.);
    hPrimaryVertexConstraintChi2 = new TH1D(
      "hPrimaryVertexConstraintChi2Ndf",
      "PV-only constraint #chi^{2}/ndf;#chi^{2}/ndf;candidates",
      250, 0., 50.);
    hDaughterPvChi2 = new TH1D(
      "hDaughterPvChi2", "Daughter-to-PV #chi^{2};#chi^{2};daughters",
      400, 0., 200.);

    std::ostringstream summary;
    summary << "mode=" << (config.dataMode == kSimulationData ? "simulation" : "real")
            << "; t0Source=" << (config.t0.source == kFixedT0 ? "fixed" : "hit")
            << "; t0Scale=" << config.t0.rawToNsScale
            << "; t0OffsetNs=" << config.t0.offsetNs
            << "; t0Aggregation=" << static_cast<int>(config.t0.aggregation)
            << "; t0SigmaNs=" << config.t0.sigmaNs
            << "; fixedT0Ns=" << config.t0.fixedTimeNs
            << "; measuredT0Enabled=" << UsesMeasuredT0(config)
            << "; simulationExactTimeOrigin="
            << (config.dataMode == kSimulationData)
            << "; tofScale=" << config.tofRawToNsScale
            << "; tofOffsetNs=" << config.tofOffsetNs
            << "; tofSigmaNs=" << config.tofTimeSigmaNs
            << "; tofPathSigmaCm=" << config.tofPathLengthSigmaCm
            << "; trackInputMode=" << static_cast<int>(config.trackInputMode)
            << "; momentumMode=" << static_cast<int>(config.finder.momentumMode)
            << "; magneticTrackReference="
            << (config.finder.usePrimaryVertexDcaReference ? "pv-dca" : "stored")
            << "; referenceMaxPathCm=" << config.finder.referenceTransportConfig.maximumAbsPathCm
            << "; referenceToleranceCm=" << config.finder.referenceTransportConfig.closestApproachToleranceCm
            << "; referenceMaxIterations=" << config.finder.referenceTransportConfig.maximumClosestApproachIterations
            << "; referenceCovariance=fixed-nominal-arc; referenceProcessNoise=0"
            << "; ByTesla=" << config.finder.magneticFieldByTesla
            << "; magneticFieldAxis=y; beamAxis=z; fieldConvention=CEE-By-20260920"
            << "; iTOFMatchMaxCm=" << config.maximumItofMatchDistanceCm
            << "; eTOFMatchMaxCm=" << config.maximumEtofMatchDistanceCm
            << "; minTpcHits=" << config.minimumTpcHits
            << "; maxTrackChi2Ndf=" << config.maximumTrackChi2Ndf
            << "; requireTrackConverged=" << config.adapter.requireConvergedTrackFit
            << "; allowLegacyTrackFallback=" << config.adapter.allowLegacyStateInZeroField
            << "; allowPvErrorFallback=" << config.adapter.allowVertexErrorFallback
            << "; fallbackTrackSigmaCm=" << config.adapter.fallbackTrackSigmaXcm
            << "," << config.adapter.fallbackTrackSigmaYcm
            << "," << config.adapter.fallbackTrackSigmaZcm
            << "; fallbackDirectionSigmaRad=" << config.adapter.fallbackDirectionSigmaRad
            << "; fallbackPvSigmaCm=" << config.adapter.fallbackVertexSigmaXcm
            << "," << config.adapter.fallbackVertexSigmaYcm
            << "," << config.adapter.fallbackVertexSigmaZcm
            << "; massMinGeV=" << config.finder.lambdaCuts.minimumMassGeV
            << "; massMaxGeV=" << config.finder.lambdaCuts.maximumMassGeV
            << "; maxDaughterDcaCm=" << config.finder.lambdaCuts.maximumDaughterDcaCm
            << "; maxGeoChi2Ndf=" << config.finder.lambdaCuts.maximumGeometricChi2Ndf
            << "; maxTopoChi2Ndf=" << config.finder.lambdaCuts.maximumTopologicalChi2Ndf
            << "; minDecayLengthCm=" << config.finder.lambdaCuts.minimumDecayLengthCm
            << "; minLdl=" << config.finder.lambdaCuts.minimumDecayLengthSignificance
            << "; minPointingCos=" << config.finder.lambdaCuts.minimumPointingCosine
            << "; minDaughterPvChi2=" << config.finder.lambdaCuts.minimumDaughterChi2FromPv
            << "; legacyMinDaughterPvDeviation="
            << config.finder.lambdaCuts.minimumDaughterDeviationFromPv
            << "; lambdaRoleDcaCm=" << config.adapter.lambdaRoleDcaBoundaryCm
            << "; lambdaDcaRoleEnabled=" << config.finder.useLambdaDcaRoleAssignment
            << "; requireCompleteMotherJacobian="
            << config.finder.lambdaCuts.requireCompleteMotherJacobian
            << "; requireTofConvergence=" << config.finder.requireTofConvergence
            << "; maxTofIterations=" << config.finder.maximumTofIterations
            << "; tofVertexConvergenceCm=" << config.finder.tofVertexConvergenceCm
            << "; tofTimeConvergenceNs=" << config.finder.tofTimeConvergenceNs
            << "; maximumBeta=" << config.finder.maximumBeta
            << "; additionalDecayTimeSigmaNs=" << config.finder.additionalDecayTimeSigmaNs
            << "; propagateSharedTimeUncertainty="
            << config.finder.propagateSharedTimeUncertainty
            << "; sharedTimeDerivativeStepNs="
            << config.finder.sharedTimeDerivativeStepNs
            << "; sharedTimeDerivativeStabilityTolerance="
            << config.finder.sharedTimeDerivativeStabilityTolerance
            << "; motherJacobianPositionStepCm="
            << config.finder.vertexFitterConfig.motherJacobianPositionStepCm
            << "; motherJacobianMomentumRelativeStep="
            << config.finder.vertexFitterConfig.motherJacobianMomentumRelativeStep
            << "; motherJacobianMinimumMomentumStepGeV="
            << config.finder.vertexFitterConfig.motherJacobianMinimumMomentumStepGeV
            << "; maximumTracksPerEvent=" << config.maximumTracksPerEvent
            << "; maximumCandidates=" << config.finder.maximumCandidates
            << "; simulationTruthEnabled=" << config.enableSimulationTruth
            << "; buildProductionState=" << config.finder.buildProductionState
            << "; productionMaximumIterations=" << config.finder.productionVertexConfig.maximumIterations
            << "; productionRelativeStateTolerance=" << config.finder.productionVertexConfig.relativeStateTolerance
            << "; productionConstraintToleranceCm=" << config.finder.productionVertexConfig.constraintToleranceCm
            << "; productionPvDerivativeStepCm=" << config.finder.productionVertexConfig.pvDerivativeStepCm
            << "; productionPvDerivativeStabilityTolerance=" << config.finder.productionVertexConfig.pvDerivativeStabilityTolerance
            << "; productionMaximumAbsFlightPathCm=" << config.finder.productionVertexConfig.maximumAbsFlightPathCm
            << "; productionMaximumAbsBendingAngleRad=" << config.finder.productionVertexConfig.maximumAbsBendingAngleRad
            << "; motherTopoStateOrder=x,y,z,px,py,pz,E,S"
            << "; motherTopoCovariancePacking=lowerTriangleRowMajor36"
            << "; motherTopoModel=uniformByJointPvConstraint; motherTopoUsedForCuts=0"
            << "; motherTopoSUnit=cm/(GeV/c); motherTopoCtauUnit=cm"
            << "; motherTopoInputAssumption=independentRawTracksPvT0"
            << "; motherTopoTofPvResponse=fullSolveLinearized; motherTopoErrorModel=localFirstOrder"
            << "; motherStateOrder=x,y,z,px,py,pz,E"
            << "; motherCovariancePacking=lowerTriangleRowMajor"
            << "; topologicalChi2Definition=(geometricChi2+pvConstraintChi2)/(geometricNdf+2)"
            << "; daughterPvCutDefinition=rawChi2"
            << "; analyticFixedPathJacobian=" << config.finder.transportConfig.useAnalyticFixedPathJacobian
            << "; lambdaOnly=1; tpcTpcOnly=1; outputFormat=6; ceeV0LiteVersion=0.7.0";
    TNamed runConfig("CeeV0AnalysisConfig", summary.str().c_str());
    runConfig.Write();
    return true;
  }

  bool SelectPrimaryVertex(CeeV0Lite::PrimaryVertex& outputVertex,
                           bool& usedStoredErrors)
  {
    bool found = false;
    usedStoredErrors = false;
    int bestRank = std::numeric_limits<int>::min();
    for (int i = 0; vertices && i < vertices->GetEntriesFast(); ++i) {
      CeeTpcVertex* vertex = dynamic_cast<CeeTpcVertex*>(vertices->At(i));
      if (!vertex) continue;
      CeeV0Lite::PrimaryVertex candidate;
      std::size_t selectedIndex = 0;
      CeeV0Lite::CeeAdapterStatus status;
      if (!CeeV0Lite::MakePrimaryVertex(*vertex, candidate, &selectedIndex,
                                        config.adapter, &status)) continue;
      const std::vector<int>& ranks = vertex->GetRank();
      const int rank = selectedIndex < ranks.size()
                     ? ranks[selectedIndex] : std::numeric_limits<int>::min();
      if (!found || rank > bestRank) {
        found = true;
        bestRank = rank;
        outputVertex = candidate;
        const std::vector<double>& xError = vertex->GetXErr();
        const std::vector<double>& yError = vertex->GetYErr();
        const std::vector<double>& zError = vertex->GetZErr();
        usedStoredErrors = selectedIndex < xError.size() &&
                           selectedIndex < yError.size() &&
                           selectedIndex < zError.size() &&
                           std::isfinite(xError[selectedIndex]) && xError[selectedIndex] > 0. &&
                           std::isfinite(yError[selectedIndex]) && yError[selectedIndex] > 0. &&
                           std::isfinite(zError[selectedIndex]) && zError[selectedIndex] > 0.;
      }
    }
    return found;
  }

  T0SelectionResult EventT0() const
  {
    std::vector<double> rawTimes;
    if (UsesMeasuredT0(config) && t0Hits) {
      for (int i = 0; i < t0Hits->GetEntriesFast(); ++i) {
        CeeT0Hit* hit = dynamic_cast<CeeT0Hit*>(t0Hits->At(i));
        if (hit) rawTimes.push_back(hit->GetTime());
      }
    }
    return SelectEventT0(rawTimes, config.t0);
  }

  void SelectTracks(std::vector<SelectedTrack>& selected,
                    std::vector<ProjectionTrack>& projections) const
  {
    for (int i = 0; tracks && i < tracks->GetEntriesFast(); ++i) {
      CeeTpcTrack* track = dynamic_cast<CeeTpcTrack*>(tracks->At(i));
      if (!track) continue;
      if (track->GetnHitsFit() < config.minimumTpcHits) continue;
      const double chi2Ndf = TrackChi2Ndf(*track);
      if (config.maximumTrackChi2Ndf < 1.e29 &&
          (!std::isfinite(chi2Ndf) || chi2Ndf > config.maximumTrackChi2Ndf))
        continue;
      if (config.adapter.requireConvergedTrackFit &&
          !track->GetIsFitConverged()) continue;
      CeeV0Lite::CeeTrackReferenceState reference;
      CeeV0Lite::CeeAdapterStatus referenceStatus;
      if (!CeeV0Lite::SelectTrackReferenceState(
            *track, config.trackInputMode, reference, config.adapter,
            &referenceStatus)) {
        hReferenceStatus->Fill(static_cast<int>(referenceStatus));
        hTrackReferenceSource->Fill(3.);
        continue;
      }
      hReferenceStatus->Fill(static_cast<int>(referenceStatus));
      hTrackReferenceSource->Fill(reference.usedVertexSeed ? 1. : 2.);
      const ProjectionTrack projection = MakeProjectionTrack(reference, i);
      const double norm2 = projection.direction[0] * projection.direction[0] +
                           projection.direction[1] * projection.direction[1] +
                           projection.direction[2] * projection.direction[2];
      if (!std::isfinite(norm2) || norm2 < 1.e-24) continue;
      SelectedTrack entry;
      entry.track = track;
      entry.originalIndex = i;
      entry.usedVertexSeed = reference.usedVertexSeed;
      selected.push_back(entry);
      projections.push_back(projection);
      if (selected.size() >= config.maximumTracksPerEvent) break;
    }
  }

  void CollectTofHits(std::vector<TofHitData>& outputHits) const
  {
    for (int i = 0; itofHits && i < itofHits->GetEntriesFast(); ++i) {
      CeeiTOFHit* hit = dynamic_cast<CeeiTOFHit*>(itofHits->At(i));
      if (!hit) continue;
      TofHitData data;
      data.id = i;
      data.detector = kItofDetector;
      data.module = hit->GetDetectorID() / 100;
      data.positionCm[0] = hit->GetX();
      data.positionCm[1] = hit->GetY();
      data.positionCm[2] = hit->GetZ();
      data.rawTime = hit->GetTime();
      if (!std::isfinite(data.positionCm[0]) || !std::isfinite(data.positionCm[1]) ||
          !std::isfinite(data.positionCm[2]) || !std::isfinite(data.rawTime) ||
          !std::isfinite(data.rawTime * config.tofRawToNsScale + config.tofOffsetNs))
        continue;
      outputHits.push_back(data);
    }
    for (int i = 0; etofHits && i < etofHits->GetEntriesFast(); ++i) {
      CeeEtofHit* hit = dynamic_cast<CeeEtofHit*>(etofHits->At(i));
      if (!hit) continue;
      TofHitData data;
      data.id = i;
      data.detector = kEtofDetector;
      data.module = -1;
      data.positionCm[0] = hit->GetX();
      data.positionCm[1] = hit->GetY();
      data.positionCm[2] = hit->GetZ();
      data.rawTime = hit->GetTime();
      if (!std::isfinite(data.positionCm[0]) || !std::isfinite(data.positionCm[1]) ||
          !std::isfinite(data.positionCm[2]) || !std::isfinite(data.rawTime) ||
          !std::isfinite(data.rawTime * config.tofRawToNsScale + config.tofOffsetNs))
        continue;
      outputHits.push_back(data);
    }
  }

  TruthResult MatchTruth(const CeeV0Lite::Candidate& candidate) const
  {
    TruthResult result;
    if (!truthBranchBound || !mcTracks) return result;
    const int trackIds[2] = {candidate.daughterId[0], candidate.daughterId[1]};
    CeeMCTrack* daughters[2] = {0, 0};
    CeeMCTrack* mothers[2] = {0, 0};
    int motherIds[2] = {-1, -1};
    for (int daughter = 0; daughter < 2; ++daughter) {
      if (trackIds[daughter] < 0 || trackIds[daughter] >= tracks->GetEntriesFast())
        return result;
      CeeTpcTrack* reconstructed = dynamic_cast<CeeTpcTrack*>(tracks->At(trackIds[daughter]));
      if (!reconstructed) return result;
      const int mcId = reconstructed->GetMCTrackID();
      if (daughter == 0) result.daughterMcId0 = mcId;
      else result.daughterMcId1 = mcId;
      if (mcId < 0 || mcId >= mcTracks->GetEntriesFast()) return result;
      daughters[daughter] = dynamic_cast<CeeMCTrack*>(mcTracks->At(mcId));
      if (!daughters[daughter]) return result;
      motherIds[daughter] = daughters[daughter]->GetMotherId();
      if (daughter == 0) result.daughterMcPdg0 = daughters[daughter]->GetPdgCode();
      else result.daughterMcPdg1 = daughters[daughter]->GetPdgCode();
    }
    // Both reconstructed tracks now have resolvable MC particles. This is
    // enough to classify a combination as truth-known background even when a
    // daughter is primary (motherId=-1) or the two mothers differ.
    result.available = 1;
    result.matched = 0;
    result.sameMotherLambda = 0;
    result.roleCorrect = daughters[0]->GetPdgCode() == 2212 &&
                         daughters[1]->GetPdgCode() == -211 ? 1 : 0;
    if (motherIds[0] >= 0 && motherIds[0] == motherIds[1]) {
      if (motherIds[0] >= mcTracks->GetEntriesFast()) {
        result.available = -1;
        result.matched = -1;
        result.sameMotherLambda = -1;
        result.roleCorrect = -1;
        return result;
      }
      mothers[0] = dynamic_cast<CeeMCTrack*>(mcTracks->At(motherIds[0]));
      if (!mothers[0]) {
        result.available = -1;
        result.matched = -1;
        result.sameMotherLambda = -1;
        result.roleCorrect = -1;
        return result;
      }
      result.motherId = motherIds[0];
      result.motherPdg = mothers[0]->GetPdgCode();
      result.sameMotherLambda = result.motherPdg == 3122 ? 1 : 0;
    }
    result.matched = result.sameMotherLambda && result.roleCorrect ? 1 : 0;
    return result;
  }

  void FillGeneratedLambdaTruth()
  {
    if (!truthBranchBound || !mcTracks) return;
    int count = 0;
    for (int i = 0; i < mcTracks->GetEntriesFast(); ++i) {
      CeeMCTrack* particle = dynamic_cast<CeeMCTrack*>(mcTracks->At(i));
      if (!particle || particle->GetPdgCode() != 3122) continue;
      ++count;
      hGeneratedLambdaPt->Fill(particle->GetPt());
      hGeneratedLambdaRapidity->Fill(particle->GetRapidity());
    }
    hGeneratedLambdaCount->Fill(count);
  }

  void FillRow(long long eventIndex, int candidateIndex,
               const CeeV0Lite::Candidate& candidate,
               const CeeV0Lite::EventInput& event,
               const T0SelectionResult& t0Result,
               bool pvUsedStoredErrors,
               int inputTrackCount, int selectedTrackCount,
               int finderTrackCount, int matchCount, int candidateCount,
               const std::map<int, MatchMetadata>& metadata)
  {
    row.Reset();
    row.event = static_cast<int>(eventIndex);
    row.candidateIndex = candidateIndex;
    row.candidateType = static_cast<int>(candidate.type);
    row.trackId0 = candidate.daughterId[0];
    row.trackId1 = candidate.daughterId[1];
    row.daughterPdg0 = candidate.daughterPdg[0];
    row.daughterPdg1 = candidate.daughterPdg[1];
    const int daughterTrackIds[2] = {row.trackId0, row.trackId1};
    for (int daughter = 0; daughter < 2; ++daughter) {
      if (daughterTrackIds[daughter] < 0 || !tracks ||
          daughterTrackIds[daughter] >= tracks->GetEntriesFast()) continue;
      CeeTpcTrack* inputTrack = dynamic_cast<CeeTpcTrack*>(
        tracks->At(daughterTrackIds[daughter]));
      if (!inputTrack) continue;
      CeeV0Lite::CeeTrackReferenceState reference;
      CeeV0Lite::CeeAdapterStatus referenceStatus;
      const bool validReference = CeeV0Lite::SelectTrackReferenceState(
        *inputTrack, config.trackInputMode, reference, config.adapter,
        &referenceStatus);
      if (daughter == 0) {
        row.nHitsFit0 = inputTrack->GetnHitsFit();
        row.fitChi2Ndf0 = TrackChi2Ndf(*inputTrack);
        row.dEdx0 = inputTrack->GetdEdx();
        row.fitCharge0 = inputTrack->GetFitCharge();
        row.usedVertexSeed0 = validReference ? (reference.usedVertexSeed ? 1 : 0) : -1;
      }
      else {
        row.nHitsFit1 = inputTrack->GetnHitsFit();
        row.fitChi2Ndf1 = TrackChi2Ndf(*inputTrack);
        row.dEdx1 = inputTrack->GetdEdx();
        row.fitCharge1 = inputTrack->GetFitCharge();
        row.usedVertexSeed1 = validReference ? (reference.usedVertexSeed ? 1 : 0) : -1;
      }
    }
    const std::map<int, MatchMetadata>::const_iterator first = metadata.find(row.trackId0);
    const std::map<int, MatchMetadata>::const_iterator second = metadata.find(row.trackId1);
    if (first != metadata.end()) {
      row.tofDetector0 = first->second.detector;
      row.tofHitId0 = first->second.hitId;
      row.tofMatchDistance0Cm = first->second.distanceCm;
      row.tofRawTime0 = first->second.rawTime;
      row.tofTime0Ns = first->second.timeNs;
    }
    if (second != metadata.end()) {
      row.tofDetector1 = second->second.detector;
      row.tofHitId1 = second->second.hitId;
      row.tofMatchDistance1Cm = second->second.distanceCm;
      row.tofRawTime1 = second->second.rawTime;
      row.tofTime1Ns = second->second.timeNs;
    }
    row.startTimeRaw = t0Result.selectedRawTime;
    row.startTimeNs = t0Result.timeNs;
    row.sigmaStartTimeNs = t0Result.sigmaNs;
    row.numberT0Hits = t0Result.numberOfHits;
    row.numberInputTracks = inputTrackCount;
    row.numberSelectedTracks = selectedTrackCount;
    row.numberFinderTracks = finderTrackCount;
    row.numberTofMatches = matchCount;
    row.numberCandidates = candidateCount;
    row.pvX = event.primaryVertex.positionCm[0];
    row.pvY = event.primaryVertex.positionCm[1];
    row.pvZ = event.primaryVertex.positionCm[2];
    row.pvSigmaX = std::sqrt(std::max(0., event.primaryVertex.covarianceCm2.m[0][0]));
    row.pvSigmaY = std::sqrt(std::max(0., event.primaryVertex.covarianceCm2.m[1][1]));
    row.pvSigmaZ = std::sqrt(std::max(0., event.primaryVertex.covarianceCm2.m[2][2]));
    row.pvUsedStoredErrors = pvUsedStoredErrors ? 1 : 0;
    row.svX = candidate.decayVertexCm[0];
    row.svY = candidate.decayVertexCm[1];
    row.svZ = candidate.decayVertexCm[2];
    if (candidate.mother.covarianceValid) {
      row.svSigmaX = std::sqrt(std::max(
        0., candidate.decayVertexCovarianceCm2.m[0][0]));
      row.svSigmaY = std::sqrt(std::max(
        0., candidate.decayVertexCovarianceCm2.m[1][1]));
      row.svSigmaZ = std::sqrt(std::max(
        0., candidate.decayVertexCovarianceCm2.m[2][2]));
    }
    else {
      const double invalid = std::numeric_limits<double>::quiet_NaN();
      row.svSigmaX = row.svSigmaY = row.svSigmaZ = invalid;
    }
    row.px = candidate.momentumGeV[0];
    row.py = candidate.momentumGeV[1];
    row.pz = candidate.momentumGeV[2];
    row.energyGeV = candidate.mother.energyGeV;
    row.massGeV = candidate.massGeV;
    row.sigmaMassGeV = candidate.sigmaMassGeV;
    row.motherCovarianceValid = candidate.mother.covarianceValid ? 1 : 0;
    row.motherCovarianceStatus =
      static_cast<int>(candidate.mother.covarianceStatus);
    row.motherJacobianComplete = candidate.motherJacobianComplete ? 1 : 0;
    row.motherJacobianPerturbationFailed =
      candidate.mother.covarianceStatus ==
        CeeV0Lite::kMotherCovarianceJacobianPerturbationFailed ? 1 : 0;
    row.motherJacobianFitsAttempted =
      candidate.mother.jacobianFitsAttempted;
    row.motherJacobianFitsSucceeded =
      candidate.mother.jacobianFitsSucceeded;
    row.motherJacobianFailureCount =
      candidate.mother.jacobianFailureCount;
    row.motherJacobianPlusFailureMask = static_cast<int>(
      candidate.mother.jacobianPlusFailureMask);
    row.motherJacobianMinusFailureMask = static_cast<int>(
      candidate.mother.jacobianMinusFailureMask);
    row.motherJacobianFailedColumn =
      candidate.mother.jacobianFailedColumn;
    row.motherJacobianFailedSide = candidate.mother.jacobianFailedSide;
    row.motherJacobianFirstFailureFitStatus =
      candidate.mother.jacobianFirstFailureFitStatus;
    row.covarianceDependentCutsEvaluated =
      candidate.covarianceDependentCutsEvaluated ? 1 : 0;
    int packedIndex = 0;
    for (int i = 0; i < 7; ++i)
      for (int j = 0; j <= i; ++j)
        row.motherCovariance[packedIndex++] =
          candidate.mother.covariance.m[i][j];
    const CeeV0Lite::ProductionState& topo = candidate.motherTopo;
    row.motherTopoValid = topo.valid;
    row.motherTopoStatus = static_cast<int>(topo.status);
    row.motherTopoCorrelationModel = static_cast<int>(topo.correlationModel);
    row.motherTopoAtProductionVertex = topo.atProductionVertex;
    row.motherTopoIterations = topo.iterations;
    row.motherTopoCharge = topo.charge;
    row.motherTopoByTesla = topo.magneticFieldByTesla;
    row.motherTopoPvNdf = topo.pvConstraintNdf;
    row.motherTopoNdf = topo.totalNdf;
    row.motherTopoFlightValid = candidate.topoFlight.valid;
    for (int i=0; i<3; ++i) {
      row.motherTopoParameters[i] = topo.positionCm[i];
      row.motherTopoParameters[i+3] = topo.momentumGeV[i];
    }
    row.motherTopoParameters[6] = topo.energyGeV;
    row.motherTopoParameters[7] = topo.flightParameterS;
    row.motherTopoS = topo.flightParameterS;
    row.motherTopoSigmaS = topo.valid ? std::sqrt(std::max(0.,topo.covariance.m[7][7]))
      : std::numeric_limits<double>::quiet_NaN();
    int topoIndex = 0;
    for (int i=0; i<8; ++i)
      for (int j=0; j<=i; ++j) row.motherTopoCovariance[topoIndex++] = topo.covariance.m[i][j];
    row.motherTopoPvChi2 = topo.pvConstraintChi2;
    row.motherTopoChi2 = topo.totalChi2;
    row.motherTopoChi2Ndf = topo.valid && topo.totalNdf>0 ? topo.totalChi2/topo.totalNdf
      : std::numeric_limits<double>::quiet_NaN();
    row.motherTopoMassGeV = candidate.topoFlight.massGeV;
    row.motherTopoSigmaMassGeV = candidate.topoFlight.sigmaMassGeV;
    row.motherTopoDecayLengthCm = candidate.topoFlight.decayLengthCm;
    row.motherTopoSigmaDecayLengthCm = candidate.topoFlight.sigmaDecayLengthCm;
    row.motherTopoDecayLengthSignificance = candidate.topoFlight.valid && candidate.topoFlight.sigmaDecayLengthCm>0.
      ? candidate.topoFlight.decayLengthCm/candidate.topoFlight.sigmaDecayLengthCm
      : std::numeric_limits<double>::quiet_NaN();
    row.motherTopoCtauCm = candidate.topoFlight.ctauCm;
    row.motherTopoSigmaCtauCm = candidate.topoFlight.sigmaCtauCm;
    row.daughterDcaCm = candidate.daughterDcaCm;
    row.geometricChi2 = candidate.geometricChi2;
    row.geometricNdf = candidate.geometricNdf;
    row.geometricChi2Ndf = candidate.geometricChi2Ndf;
    row.topologicalChi2 = candidate.topologicalChi2;
    row.topologicalNdf = candidate.topologicalNdf;
    row.topologicalChi2Ndf = candidate.topologicalChi2Ndf;
    row.primaryVertexConstraintChi2 =
      candidate.primaryVertexConstraintChi2;
    row.primaryVertexConstraintNdf =
      candidate.primaryVertexConstraintNdf;
    row.primaryVertexConstraintChi2Ndf =
      candidate.primaryVertexConstraintChi2Ndf;
    row.decayLengthCm = candidate.decayLengthCm;
    row.sigmaDecayLengthCm = candidate.sigmaDecayLengthCm;
    row.decayLengthSignificance = candidate.decayLengthSignificance;
    row.pointingCosine = candidate.pointingCosine;
    row.daughterDcaToPv0Cm = candidate.daughterDcaToPvCm[0];
    row.daughterDcaToPv1Cm = candidate.daughterDcaToPvCm[1];
    row.daughterPvChi20 = candidate.daughterChi2FromPv[0];
    row.daughterPvChi21 = candidate.daughterChi2FromPv[1];
    row.daughterPvSignificance0 = candidate.daughterPvSignificance[0];
    row.daughterPvSignificance1 = candidate.daughterPvSignificance[1];
    row.decayTimeNs = candidate.decayTimeNs;
    row.tofIterations = candidate.tofIterations;
    row.tofConverged = candidate.tofConverged ? 1 : 0;
    row.dataMode = config.dataMode == kSimulationData ? 0 : 1;
    row.truth = MatchTruth(candidate);
  }

  int ProcessEvents()
  {
    const long long chainEntries = chain->GetEntries();
    const long long entries = config.maximumEvents >= 0
                            ? std::min(chainEntries, config.maximumEvents)
                            : chainEntries;
    long long totalCandidates = 0;
    for (long long eventIndex = 0; eventIndex < entries; ++eventIndex) {
      if (chain->GetEntry(eventIndex) <= 0) continue;
      hEventCounter->Fill(1.);
      FillGeneratedLambdaTruth();

      CeeV0Lite::EventInput event;
      bool pvUsedStoredErrors = false;
      if (!SelectPrimaryVertex(event.primaryVertex, pvUsedStoredErrors)) continue;
      hPvStoredErrors->Fill(pvUsedStoredErrors ? 1. : 0.);
      hEventCounter->Fill(2.);

      const T0SelectionResult t0Result = EventT0();
      if (!t0Result.valid) continue;
      event.startTimeNs = t0Result.timeNs;
      event.sigmaStartTimeNs = t0Result.sigmaNs;
      hT0->Fill(event.startTimeNs);
      if (std::isfinite(t0Result.selectedRawTime)) hT0Raw->Fill(t0Result.selectedRawTime);
      hEventCounter->Fill(3.);

      const int inputTrackCount = tracks ? tracks->GetEntriesFast() : 0;
      hInputTracks->Fill(inputTrackCount);
      std::vector<SelectedTrack> selectedTracks;
      std::vector<ProjectionTrack> projections;
      SelectTracks(selectedTracks, projections);
      hSelectedTracks->Fill(selectedTracks.size());
      if (selectedTracks.size() < 2u) continue;
      hEventCounter->Fill(4.);

      std::vector<TofHitData> hitData;
      std::vector<TofMatch> matches;
      if (config.trackInputMode == CeeV0Lite::kCeeZeroFieldTofMode) {
        CollectTofHits(hitData);
        matches = AssignUniqueTofHits(
          projections, hitData, config.maximumItofMatchDistanceCm,
          config.maximumEtofMatchDistanceCm);
      }
      // Magnetic tracking mode deliberately does not attach a hit using a
      // straight tangent projection. Momentum already comes from curvature;
      // a future TOF cross-check needs a dedicated helical matcher.
      hTofMatches->Fill(matches.size());
      std::vector<int> matchForTrack(selectedTracks.size(), -1);
      for (std::size_t i = 0; i < matches.size(); ++i) {
        matchForTrack[matches[i].trackIndex] = static_cast<int>(i);
        hTofMatchDistance->Fill(matches[i].distanceCm);
      }
      if (config.trackInputMode == CeeV0Lite::kCeeZeroFieldTofMode &&
          matches.size() < 2u) continue;
      hEventCounter->Fill(5.);

      std::vector<CeeV0Lite::TrackInput> finderTracks;
      std::map<int, MatchMetadata> metadata;
      for (std::size_t i = 0; i < selectedTracks.size(); ++i) {
        CeeV0Lite::CeeTofMeasurement tof;
        const int matchIndex = matchForTrack[i];
        if (matchIndex >= 0) {
          const TofMatch& match = matches[static_cast<std::size_t>(matchIndex)];
          const TofHitData& hit = hitData[match.hitIndex];
          tof.available = true;
          tof.positionCm = CeeV0Lite::Vec3(hit.positionCm[0], hit.positionCm[1],
                                           hit.positionCm[2]);
          tof.timeNs = hit.rawTime * config.tofRawToNsScale + config.tofOffsetNs;
          if (!std::isfinite(tof.timeNs)) continue;
          tof.sigmaTimeNs = config.tofTimeSigmaNs;
          tof.sigmaPathLengthCm = config.tofPathLengthSigmaCm;
          MatchMetadata info;
          info.detector = match.detector;
          info.hitId = match.hitId;
          info.distanceCm = match.distanceCm;
          info.rawTime = hit.rawTime;
          info.timeNs = tof.timeNs;
          info.positionCm = tof.positionCm;
          metadata[selectedTracks[i].originalIndex] = info;
          hTofMinusT0->Fill(tof.timeNs - event.startTimeNs);
        }
        CeeV0Lite::TrackInput converted;
        CeeV0Lite::CeeAdapterStatus status;
        if (!CeeV0Lite::MakeTrackInput(*selectedTracks[i].track,
              selectedTracks[i].originalIndex, config.trackInputMode, tof,
              kLambdaDaughterHypotheses, converted, config.adapter, &status)) {
          if (config.debugLevel > 1)
            std::cerr << "event " << eventIndex << " track "
                      << selectedTracks[i].originalIndex << " rejected: "
                      << CeeV0Lite::CeeAdapterStatusMessage(status) << std::endl;
          hAdapterStatus->Fill(static_cast<int>(status));
          continue;
        }
        hAdapterStatus->Fill(static_cast<int>(status));
        finderTracks.push_back(converted);
      }
      if (finderTracks.size() < 2u) continue;

      std::vector<CeeV0Lite::Candidate> candidates;
      std::vector<CeeV0Lite::TrackReferenceDiagnostic> referenceDiagnostics;
      finder.Reconstruct(finderTracks, event, candidates, &referenceDiagnostics);
      for (std::size_t i = 0; i < referenceDiagnostics.size(); ++i) {
        const CeeV0Lite::TrackReferenceDiagnostic& info = referenceDiagnostics[i];
        hPvDcaReferenceStatus->Fill(static_cast<int>(info.status));
        if (info.status == CeeV0Lite::kReferenceAtPvDca) {
          hPvDcaReferencePath->Fill(info.closestApproach.signedPathCm);
          hPvDcaReferenceDca->Fill(info.dcaCm);
        } else if (info.status != CeeV0Lite::kReferenceUnchanged && config.debugLevel > 1) {
          std::cerr << "event " << eventIndex << " track " << info.trackId
            << " PV-DCA reference rejected: status=" << info.status
            << " converged=" << info.closestApproach.converged
            << " pathLimit=" << info.closestApproach.hitPathLimit << std::endl;
        }
      }
      hEventCounter->Fill(6.);
      if (!candidates.empty()) hEventCounter->Fill(7.);
      for (std::size_t candidateIndex = 0; candidateIndex < candidates.size();
           ++candidateIndex) {
        const CeeV0Lite::Candidate& candidate = candidates[candidateIndex];
        if (candidate.type != CeeV0Lite::kLambda) continue;
        FillRow(eventIndex, static_cast<int>(candidateIndex), candidate, event,
                t0Result, pvUsedStoredErrors, inputTrackCount,
                static_cast<int>(selectedTracks.size()),
                static_cast<int>(finderTracks.size()),
                static_cast<int>(matches.size()),
                static_cast<int>(candidates.size()), metadata);
        candidateTree->Fill();
        hMass->Fill(candidate.massGeV);
        if (row.truth.matched == 1) hMassTruth->Fill(candidate.massGeV);
        else if (row.truth.matched == 0) hMassBackground->Fill(candidate.massGeV);
        hDaughterDca->Fill(candidate.daughterDcaCm);
        hDecayLength->Fill(candidate.decayLengthCm);
        hPointingCosine->Fill(candidate.pointingCosine);
        hGeometricChi2->Fill(candidate.geometricChi2Ndf);
        if (std::isfinite(candidate.topologicalChi2Ndf))
          hTopologicalChi2->Fill(candidate.topologicalChi2Ndf);
        if (std::isfinite(candidate.primaryVertexConstraintChi2Ndf))
          hPrimaryVertexConstraintChi2->Fill(
            candidate.primaryVertexConstraintChi2Ndf);
        hDaughterPvChi2->Fill(candidate.daughterChi2FromPv[0]);
        hDaughterPvChi2->Fill(candidate.daughterChi2FromPv[1]);
        ++totalCandidates;
      }

      if (config.debugLevel > 0 &&
          (eventIndex < 10 || (eventIndex + 1) % 1000 == 0)) {
        std::cout << "event " << eventIndex << ": inputTracks="
                  << inputTrackCount << ", finderTracks=" << finderTracks.size()
                  << ", uniqueTof=" << matches.size() << ", Lambda="
                  << candidates.size() << std::endl;
      }
    }
    std::cout << "CeeV0Analysis: wrote " << totalCandidates
              << " Lambda candidate(s)" << std::endl;
    return 0;
  }

  int Run()
  {
    std::string timingError;
    if (!ValidateDataModeTiming(config, &timingError)) {
      std::cerr << "CeeV0Analysis: " << timingError << std::endl;
      return 1;
    }
    const bool fieldOn = std::fabs(config.finder.magneticFieldByTesla) >
                         config.finder.zeroFieldThresholdTesla;
    if (config.trackInputMode == CeeV0Lite::kCeeZeroFieldTofMode &&
        (config.finder.momentumMode != CeeV0Lite::kTofMomentum || fieldOn)) {
      std::cerr << "CeeV0Analysis: inconsistent zero-field adapter/finder "
                   "configuration" << std::endl;
      return 1;
    }
    if (config.trackInputMode == CeeV0Lite::kCeeMagneticFieldTrackingMode &&
        (config.finder.momentumMode != CeeV0Lite::kTrackingMomentum || !fieldOn)) {
      std::cerr << "CeeV0Analysis: inconsistent magnetic tracking adapter/finder "
                   "configuration" << std::endl;
      return 1;
    }
    if (!SetupInput()) return 2;
    if (!SetupOutput()) return 3;
    const int status = ProcessEvents();
    output->cd();
    // All trees and histograms are owned by the output directory. One file
    // Write avoids duplicate ROOT key cycles.
    output->Write();
    output->Close();
    return status;
  }

  AnalysisConfig config;
  TChain* chain;
  TClonesArray* tracks;
  TClonesArray* vertices;
  TClonesArray* itofHits;
  TClonesArray* etofHits;
  TClonesArray* t0Hits;
  TClonesArray* mcTracks;
  bool truthBranchBound;
  TFile* output;
  TTree* candidateTree;
  CandidateRow row;
  TH1D* hEventCounter;
  TH1D* hInputTracks;
  TH1D* hSelectedTracks;
  TH1D* hTofMatches;
  TH1D* hTofMatchDistance;
  TH1D* hT0;
  TH1D* hT0Raw;
  TH1D* hTofMinusT0;
  TH1D* hTrackReferenceSource;
  TH1D* hReferenceStatus;
  TH1D* hAdapterStatus;
  TH1D* hPvStoredErrors;
  TH1D* hPvDcaReferenceStatus;
  TH1D* hPvDcaReferencePath;
  TH1D* hPvDcaReferenceDca;
  TH1D* hGeneratedLambdaCount;
  TH1D* hGeneratedLambdaPt;
  TH1D* hGeneratedLambdaRapidity;
  TH1D* hMass;
  TH1D* hMassTruth;
  TH1D* hMassBackground;
  TH1D* hDaughterDca;
  TH1D* hDecayLength;
  TH1D* hPointingCosine;
  TH1D* hGeometricChi2;
  TH1D* hTopologicalChi2;
  TH1D* hPrimaryVertexConstraintChi2;
  TH1D* hDaughterPvChi2;
  CeeV0Lite::Finder finder;
};

Analyzer::Analyzer(const AnalysisConfig& config)
  : fImplementation(new Implementation(config))
{
}

Analyzer::~Analyzer()
{
  delete fImplementation;
}

int Analyzer::Run()
{
  return fImplementation->Run();
}

} // namespace CeeV0Analysis
