#include "CeeKFParticleRootWriter.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TNamed.h"
#include "TTree.h"
#include <cctype>
#include <set>
#include <sstream>
#include <stdexcept>
namespace CeeV0Analysis {
struct CeeKFParticleRootWriter::Implementation {
  Implementation(const std::string& p,const std::string& t) : path(p),treeName(t),file(0),tree(0),
    begun(false),closed(false),overwrite(false),covariances(false),eventId(0),candidateIndex(0) {}
  ~Implementation() { if (file) { file->Close(); delete file; } }
  std::string path,treeName,metadata;
  TFile* file;
  TTree* tree;
  bool begun,closed,overwrite,covariances;
  long long eventId;
  int candidateIndex;
  std::vector<CandidateObservable> columns;
  std::vector<double> values;
  double motherCov[49],topoCov[64],daughterCov0[36],daughterCov1[36];
  void Mutable() const { if (begun) throw std::logic_error("writer configuration is locked after Begin"); }
  void CheckName(const std::string& n) const {
    if (n.empty() || !(std::isalpha(static_cast<unsigned char>(n[0])) || n[0]=='_'))
      throw std::invalid_argument("invalid column name: "+n);
    for (std::size_t i=0;i<n.size();++i)
      if (!(std::isalnum(static_cast<unsigned char>(n[i])) || n[i]=='_')) throw std::invalid_argument("invalid column name: "+n);
    if (n=="eventId" || n=="candidateIndex" || n=="motherCovariance" || n=="motherTopoCovariance" ||
        n=="daughter0Covariance" || n=="daughter1Covariance") throw std::invalid_argument("reserved column name: "+n);
    for (std::size_t i=0;i<columns.size();++i) if (columns[i].name==n) throw std::invalid_argument("duplicate column name: "+n);
  }
};
CeeKFParticleRootWriter::CeeKFParticleRootWriter(const std::string& p,const std::string& t)
  : fImplementation(new Implementation(p,t)) {
  const char* names[]={"candidateType","trackId0","trackId1","massGeV","sigmaMassGeV",
    "daughterDcaCm","DcaProton","DcaPion","decayLengthCm","geometricChi2Ndf",
    "motherTopoValid","topoMassGeV","motherTopoPvChi2Ndf","LOverSigmaL","SOverSigmaS"};
  SetVariables(std::vector<std::string>(names,names+sizeof(names)/sizeof(names[0])));
}
CeeKFParticleRootWriter::~CeeKFParticleRootWriter() { delete fImplementation; }
void CeeKFParticleRootWriter::SetVariables(const std::vector<std::string>& names) {
  Implementation& w=*fImplementation; w.Mutable();
  std::vector<CandidateObservable> next;
  std::set<std::string> seen;
  for (std::size_t i=0;i<names.size();++i) {
    const CandidateObservable* o=FindCandidateObservable(names[i]);
    if (!o || !seen.insert(names[i]).second) throw std::invalid_argument("unknown/duplicate variable: "+names[i]);
    next.push_back(*o);
  }
  w.columns.swap(next);
}
void CeeKFParticleRootWriter::AddColumn(const std::string& name,CandidateObservable::Getter getter,
    const std::string& unit,const std::string& description) {
  Implementation& w=*fImplementation; w.Mutable(); w.CheckName(name);
  if (!getter) throw std::invalid_argument("empty column getter");
  CandidateObservable o={name,unit,"user",description,getter}; w.columns.push_back(o);
}
void CeeKFParticleRootWriter::SaveCovariances(bool e) { fImplementation->Mutable(); fImplementation->covariances=e; }
void CeeKFParticleRootWriter::SetOverwrite(bool e) { fImplementation->Mutable(); fImplementation->overwrite=e; }
void CeeKFParticleRootWriter::SetMetadata(const std::string& t) { fImplementation->Mutable(); fImplementation->metadata=t; }
bool CeeKFParticleRootWriter::Begin(std::string& error) {
  Implementation& w=*fImplementation; error.clear();
  if (w.begun) { error="writer Begin called more than once"; return false; }
  w.begun=true;
  TDirectory::TContext context;
  w.file=TFile::Open(w.path.c_str(),w.overwrite ? "RECREATE":"CREATE");
  if (!w.file || w.file->IsZombie()) { error="cannot create ROOT output: "+w.path; return false; }
  w.file->cd(); w.tree=new TTree(w.treeName.c_str(),"CEE KFParticle v0.7.0 selected mothers");
  w.tree->SetDirectory(w.file);
  w.tree->Branch("eventId",&w.eventId,"eventId/L");
  w.tree->Branch("candidateIndex",&w.candidateIndex,"candidateIndex/I");
  w.values.resize(w.columns.size());
  std::ostringstream schema;
  schema << "eventId\tentry index\tint64\ncandidateIndex\tselected index in event\tint32\n";
  for (std::size_t i=0;i<w.columns.size();++i) {
    const CandidateObservable& o=w.columns[i];
    w.tree->Branch(o.name.c_str(),&w.values[i],(o.name+"/D").c_str());
    schema<<o.name<<'\t'<<o.unit<<'\t'<<o.stage<<'\t'<<o.description<<'\n';
  }
  if (w.covariances) {
    w.tree->Branch("motherCovariance",w.motherCov,"motherCovariance[49]/D");
    w.tree->Branch("motherTopoCovariance",w.topoCov,"motherTopoCovariance[64]/D");
    w.tree->Branch("daughter0Covariance",w.daughterCov0,"daughter0Covariance[36]/D");
    w.tree->Branch("daughter1Covariance",w.daughterCov1,"daughter1Covariance[36]/D");
    schema<<"Covariances: row-major; mother x,y,z,px,py,pz,E; topo adds S; daughters x,y,z,px,py,pz. Invalid mother/topo matrices are NaN.\n";
  }
  TNamed("CeeKFParticleVersion","0.7.0").Write();
  TNamed("CandidateSchema",schema.str().c_str()).Write();
  TNamed("AnalysisConfiguration",w.metadata.c_str()).Write();
  return true;
}
bool CeeKFParticleRootWriter::Write(const AnalysisEvent& e,const std::vector<CeeV0Lite::Candidate>& cs,std::string& error) {
  Implementation& w=*fImplementation; error.clear();
  if (!w.tree || w.closed) { error="ROOT writer is not open"; return false; }
  const double nan=std::numeric_limits<double>::quiet_NaN();
  w.eventId=e.eventId;
  for (std::size_t k=0;k<cs.size();++k) {
    const CeeV0Lite::Candidate& c=cs[k]; w.candidateIndex=static_cast<int>(k);
    // Assign EVERY scalar and matrix for EVERY candidate; no stale branch data.
    for (std::size_t i=0;i<w.columns.size();++i) w.values[i]=w.columns[i].get(c,e);
    if (w.covariances) {
      for (int i=0;i<7;++i) for (int j=0;j<7;++j) w.motherCov[i*7+j]=c.mother.covarianceValid ? c.mother.covariance.m[i][j] : nan;
      for (int i=0;i<8;++i) for (int j=0;j<8;++j) w.topoCov[i*8+j]=c.motherTopo.valid ? c.motherTopo.covariance.m[i][j] : nan;
      for (int i=0;i<6;++i) for (int j=0;j<6;++j) {
        w.daughterCov0[i*6+j]=c.daughterAtVertex[0].covariance.m[i][j];
        w.daughterCov1[i*6+j]=c.daughterAtVertex[1].covariance.m[i][j];
      }
    }
    if (w.tree->Fill()<0) { error="ROOT tree Fill failed"; return false; }
  }
  return true;
}
bool CeeKFParticleRootWriter::Finish(std::string& error) {
  Implementation& w=*fImplementation; error.clear();
  if (w.closed) return true;
  if (!w.file || w.file->IsZombie()) { error="ROOT output was not opened"; return false; }
  TDirectory::TContext context; w.file->cd();
  const bool ok=w.file->Write()>0 && !w.file->TestBit(TFile::kWriteError);
  w.file->Close(); w.closed=true; w.tree=0;
  if (!ok) error="ROOT output Write failed";
  return ok;
}
}
