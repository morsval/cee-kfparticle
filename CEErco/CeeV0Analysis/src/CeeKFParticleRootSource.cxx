#include "CeeKFParticleRootSource.h"
#include "CeeKFParticleInputFiles.h"
#include "CeeTpcTrack.h"
#include "CeeTpcVertex.h"
#include "CeeiTOFHit.h"
#include "CeeEtofHit.h"
#include "CeeT0Hit.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TList.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>
#include <memory>
namespace CeeV0Analysis {
namespace {
using InputFiles::ResolveInputFiles;

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
  // A stored zero can mean the upstream chi2/NDF field was never filled.
  // Prefer a valid raw chi2/NDF fallback in that case (also handles true zero).
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

}
struct CeeKFParticleRootSource::Implementation {
  explicit Implementation(const AnalysisConfig& c) : config(c),chain(0),tracks(0),vertices(0),
    itofHits(0),etofHits(0),t0Hits(0),entry(0),opened(false),failed(false) {}
  ~Implementation() { delete chain; delete tracks; delete vertices; delete itofHits; delete etofHits; delete t0Hits; }
  AnalysisConfig config;
  TChain* chain;
  TClonesArray *tracks,*vertices,*itofHits,*etofHits,*t0Hits;
  long long entry;
  bool opened,failed;
  std::string failure;
  TrackIndexSelector selector;
  void Bind(const char* name,TClonesArray*& data) {
    data=new TClonesArray(name);
    chain->SetBranchStatus(name,1);
    chain->SetBranchStatus((std::string(name)+"*").c_str(),1);
    if (chain->SetBranchAddress(name,&data)<0) throw std::runtime_error(std::string("cannot bind branch: ")+name);
  }
  void Open() {
    if (opened) return;
    std::string error;
    if (!ValidateDataModeTiming(config,&error)) throw std::runtime_error(error);
    const bool magnetic=config.trackInputMode==CeeV0Lite::kCeeMagneticFieldTrackingMode;
    if ((!magnetic && std::fabs(config.finder.magneticFieldByTesla)>config.finder.zeroFieldThresholdTesla) ||
        (magnetic && config.finder.momentumMode!=CeeV0Lite::kTrackingMomentum) ||
        (!magnetic && config.finder.momentumMode!=CeeV0Lite::kTofMomentum))
      throw std::runtime_error("input mode, field and momentum mode disagree");
    if (config.maximumEvents < -1 || config.maximumTracksPerEvent==0 ||
        !std::isfinite(config.finder.magneticFieldByTesla))
      throw std::runtime_error("invalid input limits or magnetic field");
    const std::vector<std::string> files=ResolveInputFiles(config.inputPath,config.maximumInputFiles);
    if (files.empty()) throw std::runtime_error("no ROOT input files: "+config.inputPath);
    bool hasItof=false,hasEtof=false;
    // Validate EVERY file before the event loop, not just the first TChain tree.
    for (std::size_t i=0;i<files.size();++i) {
      std::unique_ptr<TFile> file(TFile::Open(files[i].c_str(),"READ"));
      if (!file.get() || file->IsZombie()) throw std::runtime_error("cannot open input: "+files[i]);
      TTree* tree=dynamic_cast<TTree*>(file->Get(config.treeName.c_str()));
      if (!tree || !tree->GetBranch("CeeTpcTrack") || !tree->GetBranch("CeeTpcVertex"))
        throw std::runtime_error("missing tree or CeeTpcTrack/CeeTpcVertex in "+files[i]);
      const bool itof=tree->GetBranch("CeeiTOFHit")!=0;
      const bool etof=tree->GetBranch("CeeEtofHit")!=0;
      if (!magnetic && (!itof && !etof)) throw std::runtime_error("zero-field TOF branches missing in "+files[i]);
      if (!magnetic && i>0 && (hasItof!=itof || hasEtof!=etof))
        throw std::runtime_error("TOF branch schema differs across input files");
      if (!magnetic && UsesMeasuredT0(config) && !tree->GetBranch("CeeT0Hit"))
        throw std::runtime_error("measured T0 branch missing in "+files[i]);
      hasItof=itof; hasEtof=etof;
    }
    chain=new TChain(config.treeName.c_str());
    for (std::size_t i=0;i<files.size();++i)
      if (!chain->Add(files[i].c_str())) throw std::runtime_error("cannot add input: "+files[i]);
    chain->SetBranchStatus("*",0);
    Bind("CeeTpcTrack",tracks); Bind("CeeTpcVertex",vertices);
    if (!magnetic && hasItof) Bind("CeeiTOFHit",itofHits);
    if (!magnetic && hasEtof) Bind("CeeEtofHit",etofHits);
    if (!magnetic && UsesMeasuredT0(config)) Bind("CeeT0Hit",t0Hits);
    // Analysis/PID never binds MCTrack. The legacy Analyzer retains truth QA.
    opened=true;
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
    if(config.trackInputMode==CeeV0Lite::kCeeMagneticFieldTrackingMode) {
      T0SelectionResult result;result.valid=true;result.timeNs=0.;result.sigmaNs=0.;
      return result;
    }
    std::vector<double> rawTimes;
    if (UsesMeasuredT0(config) && t0Hits) {
      for (int i = 0; i < t0Hits->GetEntriesFast(); ++i) {
        CeeT0Hit* hit = dynamic_cast<CeeT0Hit*>(t0Hits->At(i));
        if (hit) rawTimes.push_back(hit->GetTime());
      }
    }
    return SelectEventT0(rawTimes, config.t0);
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


  ReadStatus Next(AnalysisEvent& out) {
    Open();
    const long long entries=chain->GetEntries();
    while (entry<entries && (config.maximumEvents<0 || entry<config.maximumEvents)) {
      const long long index=entry++;
      out=AnalysisEvent();
      if (chain->GetEntry(index)<=0) throw std::runtime_error("ROOT event read failed at entry "+std::to_string(index));
      bool stored=false;
      if (!SelectPrimaryVertex(out.input.primaryVertex,stored)) continue;
      const T0SelectionResult t0=EventT0();
      if (!t0.valid) continue;
      out.eventId=index; out.inputTrackCount=tracks->GetEntriesFast();
      out.input.startTimeNs=t0.timeNs; out.input.sigmaStartTimeNs=t0.sigmaNs;
      std::vector<int> indices;
      if (selector) indices=selector(index,*tracks);
      else for (int i=0;i<tracks->GetEntriesFast();++i) indices.push_back(i);
      std::set<int> seen;
      std::vector<SelectedTrack> selected;
      std::vector<ProjectionTrack> projections;
      for (std::size_t j=0;j<indices.size();++j) {
        const int i=indices[j];
        if (i<0 || i>=tracks->GetEntriesFast() || !seen.insert(i).second)
          throw std::runtime_error("input track selector returned duplicate/out-of-range index");
        CeeTpcTrack* t=dynamic_cast<CeeTpcTrack*>(tracks->At(i));
        if (!t || t->GetnHitsFit()<config.minimumTpcHits) continue;
        const double chi2=TrackChi2Ndf(*t);
        if (config.maximumTrackChi2Ndf<1.e29 && (!std::isfinite(chi2) || chi2>config.maximumTrackChi2Ndf)) continue;
        if (config.adapter.requireConvergedTrackFit && !t->GetIsFitConverged()) continue;
        CeeV0Lite::CeeTrackReferenceState state;
        if (!CeeV0Lite::SelectTrackReferenceState(*t,config.trackInputMode,state,config.adapter)) continue;
        SelectedTrack item; item.track=t; item.originalIndex=i; item.usedVertexSeed=state.usedVertexSeed;
        selected.push_back(item); projections.push_back(MakeProjectionTrack(state,i));
        // Never silently truncate physics according to input ordering.
        if (selected.size()>config.maximumTracksPerEvent)
          throw std::runtime_error("maximumTracksPerEvent exceeded; raise limit or preselect input tracks");
      }
      std::vector<TofHitData> hits;
      std::vector<TofMatch> matches;
      if (config.trackInputMode==CeeV0Lite::kCeeZeroFieldTofMode) {
        CollectTofHits(hits);
        matches=AssignUniqueTofHits(projections,hits,config.maximumItofMatchDistanceCm,config.maximumEtofMatchDistanceCm);
      }
      std::vector<int> matchForTrack(selected.size(),-1);
      for (std::size_t i=0;i<matches.size();++i) matchForTrack[matches[i].trackIndex]=static_cast<int>(i);
      for (std::size_t i=0;i<selected.size();++i) {
        CeeV0Lite::CeeTofMeasurement tof;
        if (matchForTrack[i]>=0) {
          const TofHitData& h=hits[matches[matchForTrack[i]].hitIndex];
          tof.available=true;
          tof.positionCm=CeeV0Lite::Vec3(h.positionCm[0],h.positionCm[1],h.positionCm[2]);
          tof.timeNs=h.rawTime*config.tofRawToNsScale+config.tofOffsetNs;
          tof.sigmaTimeNs=config.tofTimeSigmaNs; tof.sigmaPathLengthCm=config.tofPathLengthSigmaCm;
        }
        AnalysisTrack t;
        if (!CeeV0Lite::MakeTrackInput(*selected[i].track,selected[i].originalIndex,
            config.trackInputMode,tof,CeeV0Lite::kPionHypothesis|CeeV0Lite::kProtonHypothesis,t.input,config.adapter)) continue;
        t.nHits=selected[i].track->GetnHitsFit(); t.fitConverged=selected[i].track->GetIsFitConverged();
        t.chi2Ndf=TrackChi2Ndf(*selected[i].track); t.dEdx=selected[i].track->GetdEdx();
        // No calibrated nSigma PID or magnetic TOF-mass2 is invented here.
        out.tracks.push_back(t);
      }
      return kEvent;
    }
    out=AnalysisEvent(); return kEnd;
  }
};
CeeKFParticleRootSource::CeeKFParticleRootSource(const AnalysisConfig& c) : fImplementation(new Implementation(c)) {}
CeeKFParticleRootSource::~CeeKFParticleRootSource() { delete fImplementation; }
void CeeKFParticleRootSource::SetTrackIndexSelector(TrackIndexSelector f) { fImplementation->selector=f; }
bool CeeKFParticleRootSource::Open(std::string& error) {
  error.clear();
  if (fImplementation->failed) { error=fImplementation->failure; return false; }
  try { fImplementation->Open(); return true; }
  catch(const std::exception& ex) { error=ex.what(); fImplementation->failed=true; fImplementation->failure=error; return false; }
}
CeeKFParticleEventSource::ReadStatus CeeKFParticleRootSource::Next(AnalysisEvent& e,std::string& error) {
  e=AnalysisEvent(); error.clear();
  if (!Open(error)) return kError;
  try { return fImplementation->Next(e); }
  catch(const std::exception& ex) { error=ex.what(); fImplementation->failed=true; fImplementation->failure=error; e=AnalysisEvent(); return kError; }
}
}
