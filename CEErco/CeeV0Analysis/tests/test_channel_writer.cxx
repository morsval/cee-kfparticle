#include "CeeKFParticleChannelWriter.h"
#include <cassert>
#include <cstdio>
#include <iostream>
#include <unistd.h>
using namespace CeeV0Analysis;
int main() {
  const auto path =
      "/tmp/cee_channel_writer_" + std::to_string(getpid()) + ".root";
  AnalysisEvent e;
  e.eventId = 9007199254740993LL;
  ChannelCandidate lambda;
  lambda.id = 0;
  lambda.channel = "lambda";
  lambda.pdg = 3122;
  lambda.leafTrackIds = {1, 2};
  lambda.daughters.resize(2);
  lambda.daughters[0].role = "p";
  lambda.daughters[0].trackId = 1;
  lambda.daughters[1].role = "pi";
  lambda.daughters[1].trackId = 2;
  ChannelCandidate omega;
  omega.id = 1;
  omega.channel = "omega";
  omega.pdg = 3334;
  omega.leafTrackIds = {1, 2, 3};
  omega.daughters.resize(2);
  omega.daughters[0].role = "lambda";
  omega.daughters[0].candidateId = 0;
  omega.daughters[1].role = "kaon";
  omega.daughters[1].trackId = 3;
  omega.daughters[1].dcaToPvCm = .25;
  omega.mother.covariance.m[6][6] = .0123;
  std::string error;
  {
    CeeKFParticleChannelWriter w(path);
    w.SaveCovariances(true);
    w.AddColumn(
        "myKaonDca", [](const ChannelCandidate &c, const AnalysisEvent &) {
          return c.channel == "omega" ? c.Daughter("kaon").dcaToPvCm : -1.;
        });
    assert(w.Begin(error));
    assert(w.Write(e, {lambda, omega}, error));
    assert(w.Finish(error));
  }
  TFile file(path.c_str());
  auto t = dynamic_cast<TTree *>(file.Get("mothers"));
  assert(t && t->GetEntries() == 2);
  long long id;
  std::vector<int> *links = nullptr;
  std::vector<int> *leaves = nullptr;
  std::vector<double> *cov = nullptr;
  std::vector<std::string> *roles = nullptr;
  double dca = 0;
  t->SetBranchAddress("eventId", &id);
  t->SetBranchAddress("daughterCandidateId", &links);
  t->SetBranchAddress("leafTrackIds", &leaves);
  t->SetBranchAddress("motherCovariance", &cov);
  t->SetBranchAddress("daughterRole", &roles);
  t->SetBranchAddress("myKaonDca", &dca);
  t->GetEntry(1);
  assert(id == e.eventId && links->at(0) == 0 && links->at(1) == -1 &&
         leaves->size() == 3 && roles->at(1) == "kaon" &&
         cov->at(48) == .0123 && dca == .25);
  CeeKFParticleChannelWriter collision(path);
  assert(!collision.Begin(error));
  file.Close();
  std::remove(path.c_str());
  std::cout << "ROOT generic daughter graph and covariance round-trip passed\n";
}
