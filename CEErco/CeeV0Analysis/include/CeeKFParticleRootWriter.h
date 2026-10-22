#ifndef CEE_KFP_ROOT_WRITER_H
#define CEE_KFP_ROOT_WRITER_H
#include "CeeKFParticleObservables.h"
namespace CeeV0Analysis {
class CeeKFParticleRootWriter : public CeeKFParticleCandidateSink {
 public:
  explicit CeeKFParticleRootWriter(const std::string& path,
                                  const std::string& treeName="MotherCandidates");
  ~CeeKFParticleRootWriter();
  // Configuration is locked by Begin(). Unknown/duplicate names are errors.
  void SetVariables(const std::vector<std::string>& names);
  void AddColumn(const std::string& name,CandidateObservable::Getter getter,
                 const std::string& unit="1",const std::string& description="User observable");
  void SaveCovariances(bool enable=true);
  void SetOverwrite(bool enable=true); // default false (ROOT CREATE)
  void SetMetadata(const std::string& text);
  bool Begin(std::string& error);
  bool Write(const AnalysisEvent&,const std::vector<CeeV0Lite::Candidate>&,std::string& error);
  bool Finish(std::string& error);
 private:
  CeeKFParticleRootWriter(const CeeKFParticleRootWriter&);
  CeeKFParticleRootWriter& operator=(const CeeKFParticleRootWriter&);
  struct Implementation;
  Implementation* fImplementation;
};
}
#endif
