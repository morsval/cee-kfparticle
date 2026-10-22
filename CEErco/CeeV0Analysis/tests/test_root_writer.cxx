#include "CeeKFParticleRootWriter.h"
#include "../examples/UserTreeSink.h"
#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <unistd.h>
using namespace CeeV0Analysis;
int main() {
  const std::string path="/tmp/cee_kfp_writer_"+std::to_string(getpid())+".root";
  std::remove(path.c_str()); std::string error;
  AnalysisEvent e; e.eventId=9007199254740993LL; // must not be stored as double
  CeeV0Lite::Candidate c;
  c.daughterPdg[0]=2212;c.daughterPdg[1]=-211;
  c.daughterDcaToPvCm[0]=0.4; c.daughterDcaToPvCm[1]=2.; c.massGeV=1.115;
  c.mother.covarianceValid=true; c.mother.covariance.m[0][0]=0.125;
  c.motherTopo.valid=true;c.motherTopo.flightParameterS=10.;c.motherTopo.covariance.m[7][7]=4.;
  c.topoFlight.valid=true;c.topoFlight.decayLengthCm=5.;c.topoFlight.sigmaDecayLengthCm=1.;
  {
    CeeKFParticleRootWriter out(path);
    out.SetVariables({"massGeV","DcaProton","DcaPion","flightS","sigmaS","LOverSigmaL"});
    out.AddColumn("myMassMeV",[](const CeeV0Lite::Candidate& x,const AnalysisEvent&){return x.massGeV*1000.;},"MeV/c^2");
    out.SaveCovariances();out.SetMetadata("integration test");
    bool thrown=false;try{out.AddColumn("massGeV",CandidateObservable::Getter());}catch(const std::invalid_argument&){thrown=true;}assert(thrown);
    assert(out.Begin(error));
    thrown=false;try{out.SetVariables({"massGeV"});}catch(const std::logic_error&){thrown=true;}assert(thrown);
    assert(out.Write(e,{c},error));
    c.motherTopo.valid=false;c.mother.covarianceValid=false;c.daughterDcaToPvCm[0]=0.9;
    assert(out.Write(e,{c},error));assert(out.Finish(error));assert(out.Finish(error));
    assert(!out.Write(e,{c},error));
  }
  {
    TFile file(path.c_str(),"READ");assert(!file.IsZombie());
    TTree* tree=dynamic_cast<TTree*>(file.Get("MotherCandidates"));assert(tree && tree->GetEntries()==2);
    assert(!tree->GetBranch("pointingCosine"));assert(tree->GetBranch("myMassMeV"));
    assert(dynamic_cast<TNamed*>(file.Get("CandidateSchema")));assert(dynamic_cast<TNamed*>(file.Get("AnalysisConfiguration")));
    double dca=0.,s=0.,sig=0.,ratio=0.,custom=0.,matrix[64],original[49];long long event=0;
    tree->SetBranchAddress("eventId",&event);tree->SetBranchAddress("DcaProton",&dca);
    tree->SetBranchAddress("flightS",&s);tree->SetBranchAddress("sigmaS",&sig);
    tree->SetBranchAddress("LOverSigmaL",&ratio);tree->SetBranchAddress("myMassMeV",&custom);
    tree->SetBranchAddress("motherTopoCovariance",matrix);tree->SetBranchAddress("motherCovariance",original);
    tree->GetEntry(0);assert(event==e.eventId && dca==0.4 && s==10. && sig==2. && ratio==5.);
    assert(custom==1115. && matrix[63]==4. && original[0]==0.125);
    tree->GetEntry(1);assert(dca==0.9 && std::isnan(s) && std::isnan(sig) && std::isnan(ratio));
    for (int i=0;i<64;++i) assert(std::isnan(matrix[i]));for(int i=0;i<49;++i)assert(std::isnan(original[i]));
  }
  {
    CeeKFParticleRootWriter collision(path); assert(!collision.Begin(error)); // CREATE protects existing output
  }
  std::remove(path.c_str());
  std::cout<<"Real ROOT writer: schema, custom columns, int64 IDs, NaN reset, covariance round trip, overwrite protection passed\n";
}
