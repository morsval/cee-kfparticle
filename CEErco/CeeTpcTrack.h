#ifndef CEETPCTRACK_H
#define CEETPCTRACK_H

#include <iostream>
#include <vector>
#include <TVector3.h>
#include <TMatrixDSym.h>
#include "FairTask.h"
#include "FairMultiLinkedData.h"

#include <TObject.h>
#include <limits>
#include <assert.h>
#include <cmath>
#include <functional>
#include <fstream>      // std::ifstream
#include <string>

//--------------------------

using namespace std; //for sample purposes
class CeeTpcTrack : public FairMultiLinkedData
{

 public:


    
  // FitParameter[0]=DCAxy (cm, signed, vs beam z-axis), [1]=DCAz (cm, z at xy-foot from origin minus ref z),
  // [2]=Phi (rad, atan2(py,px)), [3]=Theta (rad, arccos(pz/|p|) w.r.t. +z), [4]=q/p from GenFit.
  // Same extrapolation as fX,fY,fZ: extrapolateToLine through (0,0,0) along +z; DCA ref point (0,0,0).
  // Fitcov: symmetric 5x5 covariance (linear map from GenFit 5D on that surface).
  double FitParameter[5];
  double Fitcov[5][5];

  //genfit::Track *fitTrack;
  //genfit::AbsTrackRep *rep;

  
  CeeTpcTrack();

  // Add by liucy
  CeeTpcTrack(double x_, double y_, double z_,
	  double px_, double py_, double pz_,
	  double Fitcharge, double Fitchi2,
	  int Sector, int EventID, int MCTrackID,int nHitsFit, double dEdx);  // CeeTpcTrack(double Fitcharge_, double Fitchi2_, const double FitParameter_[5], 


  CeeTpcTrack(double x_, double y_, double z_, double px_, double py_, double pz_, double Fitcharge_, double Fitchi2_, const double FitParameter_[5],const double Fitcov_[5][5],double Fitindex_, double FitTracklength_);
  // CeeTpcTrack(double Fitcharge_, double Fitchi2_, const double FitParameter_[5], 
  //	   const double Fitcov_[5][5], double Fitindex_, double FitTracklength_, genfit::Track *fitTrack_, genfit::AbsTrackRep *rep_);
   
  /*
  void SetdEdx(Double_t dedx) {fdEdx=dedx;}
  double GetdEdx()  const {return fdEdx;}
  float* GetPoints() {return fPoints;}

  bool     IsSortable() const {return kTRUE;}
  int      GetLabel()   const {return fLab;}
  double   GetMass()    const {return fMass;}
  int      GetNumberOfClusters() const {return fN;}
  */
 
  ~CeeTpcTrack();
  ClassDef(CeeTpcTrack,8)
  void AddHitId(int value);
  void AddHitResidual(TVector3 Res);
  
  // NEW: For vertex reconstruction - save all hit positions
  void AddHitPosition(const TVector3& pos) { fHitPositions.push_back(pos); }
  const std::vector<TVector3>& GetHitPositions() const { return fHitPositions; }
  Int_t GetNHitPositions() const { return fHitPositions.size(); }
  
  // NEW: Save hit covariance for vertex-constrained refit
  void SetHitCovariance(double sigmaX, double sigmaY, double sigmaZ) {
      fHitSigmaX = sigmaX; fHitSigmaY = sigmaY; fHitSigmaZ = sigmaZ;
  }
  void GetHitCovariance(double& sigmaX, double& sigmaY, double& sigmaZ) const {
      sigmaX = fHitSigmaX; sigmaY = fHitSigmaY; sigmaZ = fHitSigmaZ;
  }
  
  // NEW: Track fit convergence status
  Bool_t GetIsFitConverged() const { return fIsFitConverged; }
  void SetIsFitConverged(Bool_t converged) { fIsFitConverged = converged; }

  // NEW: Persisted vertex seed state so VertexingTask can skip refits
  void SetVertexSeed(const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov);
  Bool_t HasVertexSeed() const { return fHasVertexSeed; }
  const TVector3& GetVertexSeedPos() const { return fVertexSeedPos; }
  const TVector3& GetVertexSeedMom() const { return fVertexSeedMom; }
  void GetVertexSeedCov(TMatrixDSym& cov) const;
  /** Copy 5x5 covariance of (DCAxy, DCAz, Phi, Theta, q/p) into cov (resized to 5x5). */
  void GetTrackParameterCovariance(TMatrixDSym& cov) const;
  
  Int_t GetMCTrackID() const {return fMCTrackId;}
  Int_t GetEvent() const {return fEvent;}
  Int_t GetSector() const {return fSector;}
  Int_t GetnHitsFit() const {return fnHitsFit;}
  Int_t GetnHitsMax() const { return fnHitsMax; }
  void SetnHitsMax(Int_t value) { fnHitsMax = value; }
  Int_t GetnHitsDedx() const { return fnHitsDedx; }
  void SetnHitsDedx(Int_t value) { fnHitsDedx = value; }
  Double_t GetFitCharge() const { return fFitCharge; }
  Double_t GetFitChi2() const { return fFitChi2; }
  Double_t GetFitChi2Ndf() const { return fFitChi2Ndf; }
  void SetFitChi2Ndf(Double_t value) { fFitChi2Ndf = value; }
  Double_t GetFitNdf() const { return fFitNdf; }
  void SetFitNdf(Double_t value) { fFitNdf = value; }
  Double_t GetMomentum() const { return sqrt(fPx*fPx+fPy*fPy+fPz*fPz); }
  Double_t GetPx() const { return fPx; }
  Double_t GetPy() const { return fPy; }
  Double_t GetPz() const { return fPz; }
  Double_t GetX() const { return fX; }
  Double_t GetY() const { return fY; }
  Double_t GetZ() const { return fZ; }
  Double_t GetdEdx()  const {return fdEdx;}
  Double_t GetFitindex() const { return fFitindex; }
  Double_t GetFitTracklength() const { return fFitTracklength; }
  Double_t GetFakeRatio() const { return fFakeRatio; }
  Double_t GetDcaX() const { return fDcaX; }
  Double_t GetDcaY() const { return fDcaY; }
  Double_t GetDcaZ() const { return fDcaZ; }
  void SetDca(Double_t x, Double_t y, Double_t z);
  void GetDca(Double_t& x, Double_t& y, Double_t& z) const;
  const std::vector<TVector3>& GetResidual() const { return fResidual; }
  const std::vector<int>& GetHitIdList() const { return fHitIdList;}
  
  private:

  // Private Data Members (Digi Objects) ------------
  double fFakeRatio;  // fake ratio
  double fFitindex;
  double fFitTracklength;
  double fX;
  double fY;
  double fZ;
  double fPx;
  double fPy;
  double fPz;
  double fFitCharge;
  double fFitChi2;
  double fFitChi2Ndf;
  double fFitNdf;
  double fdEdx;
  double fDcaX;
  double fDcaY;
  double fDcaZ;
  std::vector<TVector3> fResidual;
  int fEvent;
  int fMCTrackId;
  int fSector;
  int fnHitsFit;
  int fnHitsMax;
  int fnHitsDedx;
  std::vector<int> fHitIdList;
  
  // NEW: For vertex reconstruction
  std::vector<TVector3> fHitPositions;  // All hit positions for vertex fitting
  double fHitSigmaX;                     // Hit error in X
  double fHitSigmaY;                     // Hit error in Y
  double fHitSigmaZ;                     // Hit error in Z
  Bool_t fIsFitConverged;                // Track fit convergence status
  Bool_t fHasVertexSeed;                 // Whether vertex seed info is filled
  TVector3 fVertexSeedPos;               // Seed position (typically DCA)
  TVector3 fVertexSeedMom;               // Seed momentum at that point
  double fVertexSeedCov[6][6];           // 6x6 covariance in (pos,mom)

};

#endif // CIRCLE_H
