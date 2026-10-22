#include "MakerTestFixtures.h"
#include "../examples/PicoDstAdapterExample.h"
class Source : public CeeKFParticleEventSource {
 public:
  int i=0; bool fail=false;
  ReadStatus Next(AnalysisEvent& e,std::string& error) {
    if (fail) { error="intentional read failure"; return kError; }
    if (i++==2) return kEnd;
    e=MakeEvent(false); e.eventId=i; if (i==2) e.tracks.clear(); return kEvent;
  }
};
class Sink : public CeeKFParticleCandidateSink {
 public:
  int begins=0,writes=0,finishes=0; bool fail=false;
  bool Begin(std::string&) { ++begins;return true; }
  bool Write(const AnalysisEvent&,const std::vector<Candidate>& c,std::string& error) {
    ++writes; if (fail) { error="intentional write failure"; return false; }
    assert(c.size()==(writes==1?1u:0u));return true;
  }
  bool Finish(std::string&) { ++finishes;return true; }
};
int main() {
  std::string error;
  PicoDstAdapterExample adapter([](AnalysisEvent& e,std::string&){e=MakeEvent(false);return CeeKFParticleEventSource::kEvent;});
  AnalysisEvent adapted; assert(adapter.Next(adapted,error)==CeeKFParticleEventSource::kEvent && adapted.tracks.size()==2);
  for (int b=0;b<2;++b) {
    AnalysisEvent e=MakeEvent(b); CeeKFParticleAnalysisMaker maker(Configuration(b));
    maker.SetPIDSelector(PID);
    assert(maker.Make(e,error)); assert(maker.GetCandidates().size()==1);
    const Candidate c=maker.GetCandidates()[0];
    assert(std::fabs(c.massGeV-1.115683)<0.0001);
    assert(c.motherTopo.valid); assert(c.topoFlight.valid);
    assert(DcaProton(c)==c.daughterDcaToPvCm[0]); assert(DcaPion(c)==c.daughterDcaToPvCm[1]);
    // Same selected input sent directly to the pre-existing Finder is identical.
    std::vector<Candidate> direct; Finder(Configuration(b)).Reconstruct(maker.GetSelectedTracks(),e.input,direct);
    assert(direct.size()==1); assert(direct[0].massGeV==c.massGeV);
    for (int i=0;i<7;++i) for (int j=0;j<7;++j) assert(direct[0].mother.covariance.m[i][j]==c.mother.covariance.m[i][j]);
    // Explicit subset avoids evaluating 200 irrelevant tracks, preserves IDs.
    const AnalysisEvent clean=e;
    for (int i=0;i<200;++i) { AnalysisTrack unused; unused.input.id=1000+i; e.tracks.push_back(unused); }
    int calls=0;
    maker.SetTrackListSelector([](const AnalysisEvent&){return std::vector<std::size_t>{1,0};});
    maker.SetPIDSelector([&calls](const AnalysisTrack& t,const AnalysisEvent& x){++calls;return PID(t,x);});
    assert(maker.Make(e,error)); assert(calls==2); assert(maker.GetSelectedTracks().size()==2);
    assert(maker.GetCandidates().size()==1); assert(maker.GetCandidates()[0].daughterId[0]==100);
    assert(maker.GetCandidates()[0].massGeV==c.massGeV);
    // Post-reconstruction cut sees a complete constrained copy and leaves raw results accessible.
    maker.SetMotherSelector([](const Candidate& x,const AnalysisEvent&){assert(x.motherTopo.valid);return false;});
    assert(maker.Make(clean,error)); assert(maker.GetCandidates().empty()); assert(maker.GetInterface().GetCandidates().size()==1);
    maker.SetTrackListSelector([](const AnalysisEvent&){return std::vector<std::size_t>{0,0};});
    assert(!maker.Make(clean,error)); assert(maker.GetCandidates().empty()); assert(maker.GetInterface().GetCandidates().empty());
    maker.SetTrackListSelector([](const AnalysisEvent&){return std::vector<std::size_t>{500};}); assert(!maker.Make(clean,error));
    maker.SetTrackListSelector([](const AnalysisEvent&){return std::vector<std::size_t>{};});
    assert(maker.Make(clean,error)); assert(maker.GetSelectedTracks().empty()); assert(maker.GetCandidates().empty());
    maker.SetTrackListSelector(CeeKFParticleAnalysisMaker::TrackListSelector());
    maker.SetEventSelector([](const AnalysisEvent&){return false;}); assert(maker.Make(clean,error)); assert(!maker.EventAccepted());
    maker.SetEventSelector(CeeKFParticleAnalysisMaker::EventSelector());
    AnalysisEvent duplicate=clean; duplicate.tracks[1].input.id=100;
    assert(!maker.Make(duplicate,error)); assert(maker.GetSelectedTracks().empty());
    maker.SetPIDSelector([](const AnalysisTrack&,const AnalysisEvent&){return 4u;}); assert(!maker.Make(clean,error));
  }
  Candidate swapped; swapped.daughterPdg[0]=-211; swapped.daughterPdg[1]=-2212;
  swapped.daughterDcaToPvCm[0]=2.;swapped.daughterDcaToPvCm[1]=1.;
  assert(DcaProton(swapped)==1. && DcaPion(swapped)==2.);
  swapped.daughterPdg[1]=211; assert(std::isnan(DcaPion(swapped))); assert(std::isnan(DcaProton(swapped)));
  std::set<std::string> names;
  for (const CandidateObservable& o:CandidateObservables()) assert(names.insert(o.name).second);
  assert(std::isnan(CandidateValue("flightS",swapped,AnalysisEvent())));
  assert(std::isnan(CandidateValue("topoMassGeV",swapped,AnalysisEvent())));
  bool thrown=false;try{CandidateValue("typo",swapped,AnalysisEvent());}catch(const std::invalid_argument&){thrown=true;}assert(thrown);
  CeeKFParticleAnalysisMaker maker(Configuration(false)); maker.SetPIDSelector(PID);
  Source src; Sink sink; assert(maker.Run(src,sink,error));assert(sink.writes==2 && sink.finishes==1); assert(maker.GetCandidates().empty());
  Source bad; bad.fail=true; Sink closing; assert(!maker.Run(bad,closing,error));assert(closing.finishes==1 && error=="intentional read failure");
  Source src2; Sink badSink; badSink.fail=true;assert(!maker.Run(src2,badSink,error));assert(badSink.finishes==1);
  std::cout<<"Maker: zero/field Finder parity, PID, subset, cuts, roles, reset and I/O lifecycle passed\n";
}
