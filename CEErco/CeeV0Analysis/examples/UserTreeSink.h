#ifndef CEE_USER_TREE_SINK_H
#define CEE_USER_TREE_SINK_H
#include "CeeKFParticleEvent.h"
#include "CeeKFParticleObservables.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
// Minimal user-owned ROOT tree. Use this when AddColumn/SetVariables is not
// enough (for example event-level arrays or experiment-specific object trees).
class UserTreeSink : public CeeV0Analysis::CeeKFParticleCandidateSink {
 public:
  explicit UserTreeSink(const std::string& p) : path(p),file(0),tree(0),eventId(0),mass(0),dcaP(0),dcaPi(0) {}
  ~UserTreeSink() { if(file){file->Close();delete file;} }
  bool Begin(std::string& error) override {
    TDirectory::TContext context;
    file=TFile::Open(path.c_str(),"CREATE");
    if (!file || file->IsZombie()) {error="cannot create user ROOT file";return false;}
    file->cd();tree=new TTree("MyLambda","User-selected Lambda branches");
    tree->SetDirectory(file);
    tree->Branch("eventId",&eventId,"eventId/L");tree->Branch("mass",&mass,"mass/D");
    tree->Branch("DcaProton",&dcaP,"DcaProton/D");tree->Branch("DcaPion",&dcaPi,"DcaPion/D");return true;
  }
  bool Write(const CeeV0Analysis::AnalysisEvent& e,const std::vector<CeeV0Lite::Candidate>& cs,std::string& error) override {
    if (!tree){error="user tree not open";return false;}
    eventId=e.eventId;
    for (const CeeV0Lite::Candidate& c:cs) {
      if(c.type!=CeeV0Lite::kLambda)continue;
      mass=c.massGeV;dcaP=CeeV0Analysis::DcaProton(c);dcaPi=CeeV0Analysis::DcaPion(c);
      if(tree->Fill()<0){error="user tree Fill failed";return false;}
    }
    return true;
  }
  bool Finish(std::string& error) override {
    if(!file || file->IsZombie()){error="user file not open";return false;}
    TDirectory::TContext context;file->cd();
    bool ok=file->Write()>0 && !file->TestBit(TFile::kWriteError);file->Close();tree=0;
    if(!ok) error="user file write failed";
    return ok;
  }
 private:
  std::string path;
  TFile* file;TTree* tree;
  long long eventId;double mass,dcaP,dcaPi;
  UserTreeSink(const UserTreeSink&)=delete;
  UserTreeSink& operator=(const UserTreeSink&)=delete;
};
#endif
