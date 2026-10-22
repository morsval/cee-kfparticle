#ifndef CEE_ROOT_TEST_TYPES_H
#define CEE_ROOT_TEST_TYPES_H
// Synthetic ROOT-streamable fixtures. These are NOT CEE production dictionaries.
#include "TObject.h"
#include "TVector3.h"
#include "TMatrixDSym.h"
#include <vector>
class CeeTpcTrack : public TObject {
 public:
  TVector3 pos,mom; double cov[6][6]={{0.}}; int charge=1,hits=30;
  bool HasVertexSeed() const {return true;}
  const TVector3& GetVertexSeedPos() const {return pos;}
  const TVector3& GetVertexSeedMom() const {return mom;}
  void GetVertexSeedCov(TMatrixDSym& c) const {c.ResizeTo(6,6);for(int i=0;i<6;++i)for(int j=0;j<6;++j)c(i,j)=cov[i][j];}
  bool GetIsFitConverged() const {return true;}
  double GetX()const{return pos.X();} double GetY()const{return pos.Y();} double GetZ()const{return pos.Z();}
  double GetPx()const{return mom.X();} double GetPy()const{return mom.Y();} double GetPz()const{return mom.Z();}
  double GetFitCharge()const{return charge;} double GetFitChi2()const{return 5.;}
  double GetFitChi2Ndf()const{return 0.;} double GetFitNdf()const{return 5.;}
  int GetnHitsFit()const{return hits;} double GetdEdx()const{return 7.;}
  ClassDef(CeeTpcTrack,1);
};
class CeeTpcVertex : public TObject {
 public:
  std::vector<double> x{0.},y{0.},z{0.},dx{0.01},dy{0.01},dz{0.01}; std::vector<int> rank{1};
  const std::vector<double>& GetX()const{return x;} const std::vector<double>& GetY()const{return y;} const std::vector<double>& GetZ()const{return z;}
  const std::vector<double>& GetXErr()const{return dx;}const std::vector<double>& GetYErr()const{return dy;}const std::vector<double>& GetZErr()const{return dz;}
  const std::vector<int>& GetRank()const{return rank;}
  ClassDef(CeeTpcVertex,1);
};
class CeeT0Hit : public TObject {public: double GetTime(){return 0.;} ClassDef(CeeT0Hit,1);};
class CeeiTOFHit : public TObject {public: int GetDetectorID()const{return 0;} double GetX()const{return 0.;} double GetY()const{return 0.;} double GetZ()const{return 0.;} double GetTime()const{return 0.;} ClassDef(CeeiTOFHit,1);};
class CeeEtofHit : public TObject {public: double GetX()const{return 0.;} double GetY()const{return 0.;} double GetZ()const{return 0.;} double GetTime()const{return 0.;} ClassDef(CeeEtofHit,1);};
#endif
