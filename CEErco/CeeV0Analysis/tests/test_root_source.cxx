#include "MakerTestFixtures.h"
#include "CeeKFParticleRootSource.h"
#include "root_fixtures/CeeRootTestTypes.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <fstream>
#include <unistd.h>
void InputFile(const std::string& name,bool withPV=true) {
  TFile f(name.c_str(),"RECREATE");TTree tree("ceesim","synthetic CEE-shaped schema");
  TClonesArray* tracks=new TClonesArray("CeeTpcTrack");TClonesArray* vertices=new TClonesArray("CeeTpcVertex");
  tree.Branch("CeeTpcTrack",&tracks);
  if(withPV)tree.Branch("CeeTpcVertex",&vertices);
  new((*vertices)[0])CeeTpcVertex();
  AnalysisEvent e=MakeEvent(true);
  for(int i=0;i<2;++i){
    CeeTpcTrack* t=new((*tracks)[i])CeeTpcTrack();const TrackInput& x=e.tracks[i].input;
    t->pos.SetXYZ(x.positionCm[0],x.positionCm[1],x.positionCm[2]);
    t->mom.SetXYZ(x.trackingMomentumGeV[0],x.trackingMomentumGeV[1],x.trackingMomentumGeV[2]);t->charge=x.charge;
    for(int a=0;a<6;++a)for(int b=0;b<6;++b)t->cov[a][b]=x.trackingCovariance.m[a][b];
  }
  tree.Fill();tree.Write();delete tracks;delete vertices;
}
int main(){
  const std::string stem="/tmp/cee_kfp_source_"+std::to_string(getpid());
  InputFile(stem+"a.root");InputFile(stem+"b.root");InputFile(stem+"bad.root",false);
  {std::ofstream list(stem+".list");list<<"# ROOT input list\n"<<stem<<"a.root\n"<<stem<<"b.root\n";}
  AnalysisConfig config;config.inputPath=stem+".list";config.finder=Configuration(true);
  config.trackInputMode=kCeeMagneticFieldTrackingMode;
  std::string error;AnalysisEvent e;
  {
    CeeKFParticleRootSource source(config);assert(source.Open(error));
    assert(source.Next(e,error)==CeeKFParticleEventSource::kEvent && e.tracks.size()==2 && e.eventId==0);
    assert(e.tracks[0].input.id==0 && e.tracks[0].dEdx==7. && e.tracks[0].chi2Ndf==1. && std::isnan(e.tracks[0].nSigmaProton));
    CeeKFParticleAnalysisMaker maker(config.finder);
    maker.SetPIDSelector([](const AnalysisTrack& t,const AnalysisEvent&){return t.input.charge>0?kProtonHypothesis:kPionHypothesis;});
    assert(maker.Make(e,error));assert(maker.GetCandidates().size()==1);
    assert(source.Next(e,error)==CeeKFParticleEventSource::kEvent && e.eventId==1);
    assert(source.Next(e,error)==CeeKFParticleEventSource::kEnd && e.tracks.empty());
  }
  {
    CeeKFParticleRootSource source(config);
    source.SetTrackIndexSelector([](long long,const TClonesArray&){return std::vector<int>{1};});
    assert(source.Next(e,error)==CeeKFParticleEventSource::kEvent && e.tracks.size()==1 && e.tracks[0].input.id==1);
    source.SetTrackIndexSelector([](long long,const TClonesArray&){return std::vector<int>{0,0};});
    assert(source.Next(e,error)==CeeKFParticleEventSource::kError && e.tracks.empty());
  }
  {auto real=config;ConfigureDataMode(real,kRealData);
    CeeKFParticleRootSource source(real);
    assert(source.Next(e,error)==CeeKFParticleEventSource::kEvent && e.tracks.size()==2);
    assert(e.input.startTimeNs==0. && e.input.sigmaStartTimeNs==0.); // tracking never needs T0 branch
  }
  config.maximumTracksPerEvent=1;
  {CeeKFParticleRootSource source(config);assert(source.Next(e,error)==CeeKFParticleEventSource::kError);}
  config.maximumTracksPerEvent=1000;
  {std::ofstream list(stem+".list");list<<stem<<"a.root\n"<<stem<<"bad.root\n";}
  {CeeKFParticleRootSource source(config);assert(!source.Open(error));assert(error.find("bad.root")!=std::string::npos);}
  std::remove((stem+"a.root").c_str());std::remove((stem+"b.root").c_str());std::remove((stem+"bad.root").c_str());std::remove((stem+".list").c_str());
  std::cout<<"ROOT source fixture: chained file/list, original IDs, subset, Finder integration, overflow and per-file missing-branch checks passed\n";
}
