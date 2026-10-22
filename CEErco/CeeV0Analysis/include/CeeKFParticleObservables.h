#ifndef CEE_KFP_OBSERVABLES_H
#define CEE_KFP_OBSERVABLES_H
#include "CeeKFParticleEvent.h"
#include <functional>
namespace CeeV0Analysis {
// Unique |PDG| role; returns -1 for missing/ambiguous roles (e.g. two pions).
int DaughterIndex(const CeeV0Lite::Candidate& c,int absolutePdg);
double DcaProton(const CeeV0Lite::Candidate& c);
double DcaPion(const CeeV0Lite::Candidate& c);
double Chi2PrimaryProton(const CeeV0Lite::Candidate& c);
double Chi2PrimaryPion(const CeeV0Lite::Candidate& c);
struct CandidateObservable {
  typedef std::function<double(const CeeV0Lite::Candidate&,const AnalysisEvent&)> Getter;
  std::string name,unit,stage,description;
  Getter get;
};
const std::vector<CandidateObservable>& CandidateObservables();
const CandidateObservable* FindCandidateObservable(const std::string& name);
// Unknown name throws invalid_argument; unavailable/invalid quantity is NaN.
double CandidateValue(const std::string& name,const CeeV0Lite::Candidate& c,const AnalysisEvent& e);
}
#endif
