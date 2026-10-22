#ifndef CEE_KFP_CHANNEL_WRITER_H
#define CEE_KFP_CHANNEL_WRITER_H
#include "CeeKFParticleChannels.h"
#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
namespace CeeV0Analysis {
// Header-only ROOT sink: no custom ROOT dictionary is needed for the standard
// vector<int>, vector<double> and vector<string> branches used here.
class CeeKFParticleChannelWriter : public CeeKFParticleChannelSink {
public:
  using Getter =
      std::function<double(const ChannelCandidate &, const AnalysisEvent &)>;
  explicit CeeKFParticleChannelWriter(const std::string &path) : path_(path) {
    AddColumn("massGeV", [](const ChannelCandidate &c, const AnalysisEvent &) {
      return c.massGeV;
    });
    AddColumn("sigmaMassGeV",
              [](const ChannelCandidate &c, const AnalysisEvent &) {
                return c.sigmaMassGeV;
              });
    AddColumn("chi2Ndf", [](const ChannelCandidate &c, const AnalysisEvent &) {
      return c.chi2 / c.ndf;
    });
  }
  void AddColumn(const std::string &name, Getter g) {
    if (started_)
      throw std::logic_error("configure output before Begin");
    const std::set<std::string> reserved = {"eventId",
                                            "candidateId",
                                            "pdg",
                                            "charge",
                                            "channel",
                                            "leafTrackIds",
                                            "daughterRole",
                                            "daughterPdg",
                                            "daughterCharge",
                                            "daughterTrackId",
                                            "daughterCandidateId",
                                            "daughterDcaToPvCm",
                                            "motherState",
                                            "motherCovariance",
                                            "daughterInputState",
                                            "daughterFittedState",
                                            "daughterInputCovariance",
                                            "daughterFittedCovariance"};
    if (name.empty() || reserved.count(name) || columns_.count(name) || !g)
      throw std::invalid_argument("duplicate/invalid branch name");
    columns_[name] = g;
  }
  void ClearColumns() {
    if (started_)
      throw std::logic_error("configure before Begin");
    columns_.clear();
  }
  void SaveCovariances(bool yes) {
    if (started_)
      throw std::logic_error("configure before Begin");
    covariances_ = yes;
  }
  void SetOverwrite(bool yes) {
    if (started_)
      throw std::logic_error("configure before Begin");
    overwrite_ = yes;
  }
  bool Begin(std::string &error) override {
    error.clear();
    if (started_) {
      error = "writer already started";
      return false;
    }
    started_ = true;
    file_.reset(TFile::Open(path_.c_str(), overwrite_ ? "RECREATE" : "CREATE"));
    if (!file_ || file_->IsZombie()) {
      error = "cannot create channel output";
      return false;
    }
    file_->cd();
    tree_ =
        new TTree("mothers", "Generic decay graph; links scoped by eventId");
    tree_->Branch("eventId", &eventId_, "eventId/L");
    tree_->Branch("candidateId", &id_, "candidateId/I");
    tree_->Branch("pdg", &pdg_, "pdg/I");
    tree_->Branch("charge", &charge_, "charge/I");
    tree_->Branch("channel", &channel_);
    tree_->Branch("leafTrackIds", &leaves_);
    tree_->Branch("daughterRole", &roles_);
    tree_->Branch("daughterPdg", &pdgs_);
    tree_->Branch("daughterCharge", &charges_);
    tree_->Branch("daughterTrackId", &tracks_);
    tree_->Branch("daughterCandidateId", &candidates_);
    tree_->Branch("daughterDcaToPvCm", &dca_);
    tree_->Branch("motherState", &state_);
    tree_->Branch("daughterInputState", &input_);
    tree_->Branch("daughterFittedState", &fitted_);
    if (covariances_) {
      tree_->Branch("motherCovariance", &cov_);
      tree_->Branch("daughterInputCovariance", &inputCov_);
      tree_->Branch("daughterFittedCovariance", &fittedCov_);
    }
    for (const auto &p : columns_)
      tree_->Branch(p.first.c_str(), &values_[p.first],
                    (p.first + "/D").c_str());
    TNamed metadata(
        "schema",
        "v0.7.0 channels; cm, GeV/c, GeV; state=(x,y,z,px,py,pz,E); covariance "
        "row-major 7x7; daughters concatenated in role order; -1=not "
        "applicable; all accepted intermediate candidates retained; PV "
        "constraint is optional separate copy");
    metadata.Write();
    return true;
  }
  bool Write(const AnalysisEvent &e,
             const std::vector<ChannelCandidate> &mothers,
             std::string &error) override {
    if (!tree_) {
      error = "writer not open";
      return false;
    }
    for (const auto &c : mothers) {
      eventId_ = e.eventId;
      id_ = c.id;
      pdg_ = c.pdg;
      charge_ = c.charge;
      channel_ = c.channel;
      leaves_ = c.leafTrackIds;
      roles_.clear();
      pdgs_.clear();
      charges_.clear();
      tracks_.clear();
      candidates_.clear();
      dca_.clear();
      state_.clear();
      cov_.clear();
      input_.clear();
      fitted_.clear();
      inputCov_.clear();
      fittedCov_.clear();
      append(c.mother, state_, cov_);
      for (const auto &d : c.daughters) {
        roles_.push_back(d.role);
        pdgs_.push_back(d.pdg);
        charges_.push_back(d.charge);
        tracks_.push_back(d.trackId);
        candidates_.push_back(d.candidateId);
        dca_.push_back(d.dcaToPvCm);
        append(d.inputState, input_, inputCov_);
        append(d.fittedState, fitted_, fittedCov_);
      }
      for (const auto &p : columns_)
        values_[p.first] = p.second(c, e);
      if (tree_->Fill() < 0) {
        error = "channel tree Fill failed";
        return false;
      }
    }
    return true;
  }
  bool Finish(std::string &error) override {
    if (!file_ || !tree_) {
      error = "writer not open";
      return false;
    }
    file_->cd();
    bool ok = tree_->Write() > 0 && !file_->TestBit(TFile::kWriteError);
    file_->Close();
    tree_ = nullptr;
    if (!ok)
      error = "channel tree Write failed";
    return ok;
  }

private:
  void append(const CeeV0Lite::MotherState &s, std::vector<double> &v,
              std::vector<double> &cov) {
    for (int i = 0; i < 3; ++i)
      v.push_back(s.positionCm[i]);
    for (int i = 0; i < 3; ++i)
      v.push_back(s.momentumGeV[i]);
    v.push_back(s.energyGeV);
    if (covariances_)
      for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 7; ++j)
          cov.push_back(s.covariance.m[i][j]);
  }
  std::string path_, channel_;
  bool started_ = false, overwrite_ = false, covariances_ = false;
  std::unique_ptr<TFile> file_;
  TTree *tree_ = nullptr;
  long long eventId_ = 0;
  int id_ = 0, pdg_ = 0, charge_ = 0;
  std::vector<int> leaves_, pdgs_, charges_, tracks_, candidates_;
  std::vector<std::string> roles_;
  std::vector<double> dca_, state_, cov_, input_, fitted_, inputCov_,
      fittedCov_;
  std::map<std::string, Getter> columns_;
  std::map<std::string, double> values_;
};
} // namespace CeeV0Analysis
#endif
