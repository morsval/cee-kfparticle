/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             * 
 *              GNU Lesser General Public Licence (LGPL) version 3,             *  
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#include "CeeTpcTrack.h"
#include <algorithm>
ClassImp(CeeTpcTrack)

CeeTpcTrack::CeeTpcTrack()
    : fFakeRatio(0.), fFitindex(0.), fFitTracklength(0.), fX(0.), fY(0.), fZ(0.),
      fPx(0.), fPy(0.), fPz(0.), fFitCharge(0.), fFitChi2(0.), fFitChi2Ndf(0.), fFitNdf(0.), fdEdx(0.),
      fDcaX(0.), fDcaY(0.), fDcaZ(0.), fEvent(0), fMCTrackId(-1), fSector(-1), fnHitsFit(0),
      fnHitsMax(0), fnHitsDedx(0), fHitSigmaX(0.4), fHitSigmaY(0.5), fHitSigmaZ(0.2),
      fIsFitConverged(kTRUE), fHasVertexSeed(kFALSE)
{
    fVertexSeedPos.SetXYZ(0.,0.,0.);
    fVertexSeedMom.SetXYZ(0.,0.,0.);
    for(int i = 0; i < 6; ++i) {
        for(int j = 0; j < 6; ++j) {
            fVertexSeedCov[i][j] = 0.;
        }
    }
}

//CeeTpcTrack::CeeTpcTrack(double Fitcharge_, double Fitchi2_, const double FitParameter_[5], 
//		   const double Fitcov_[5][5], double Fitindex_, double FitTracklength_, genfit::Track *fitTrack_, genfit::AbsTrackRep *rep_)



CeeTpcTrack::CeeTpcTrack(double x_, double y_, double z_, 
                   double px_, double py_, double pz_,
                   double Fitcharge_, double Fitchi2_, 
                   int Sector, int EventID, int MCTrackID,int nHitsFit, double dEdx)
: fFakeRatio(0.), fFitindex(0.), fFitTracklength(0.), fX(x_), fY(y_), fZ(z_),
 fPx(px_),fPy(py_),fPz(pz_),fFitCharge(Fitcharge_),
 fFitChi2(Fitchi2_),fFitChi2Ndf(0.),fFitNdf(0.),fdEdx(dEdx),
 fDcaX(0.), fDcaY(0.), fDcaZ(0.),
 fEvent(EventID),fMCTrackId(MCTrackID),fSector(Sector),fnHitsFit(nHitsFit),
 fnHitsMax(0), fnHitsDedx(0), fIsFitConverged(kTRUE),
 fHasVertexSeed(kFALSE)  // Vertex seed not set until explicitly stored
{
  // Initialize FitParameter and Fitcov arrays to zero
  for(int i = 0; i < 5; i++) {
    FitParameter[i] = 0.0;
    for(int j = 0; j < 5; j++) {
      Fitcov[i][j] = 0.0;
    }
  }
  
  fVertexSeedPos.SetXYZ(0.,0.,0.);
  fVertexSeedMom.SetXYZ(0.,0.,0.);
  for(int i = 0; i < 6; ++i) {
    for(int j = 0; j < 6; ++j) {
      fVertexSeedCov[i][j] = 0.;
    }
  }
  
  //cout<<" SUCCEED Write"<<endl;
}



CeeTpcTrack::CeeTpcTrack(double x_, double y_, double z_, double px_, double py_, double pz_,double Fitcharge_, double Fitchi2_, const double FitParameter_[5],const double Fitcov_[5][5],double Fitindex_, double FitTracklength_)
: fIsFitConverged(kTRUE),
  fHasVertexSeed(kFALSE)
{
  fFitCharge = Fitcharge_;
  fFitChi2 = Fitchi2_;           
  fFitChi2Ndf = 0.0;
  fFitNdf = 0.0;

  for(int i =0; i<5; i++) {FitParameter[i] = FitParameter_[i];};

  
  for(int i =0; i<5; i++)
    for(int j =0; j<5; j++)
      {
	{Fitcov[i][j] = Fitcov_[i][j];};
      }
  
  fFitindex = Fitindex_;
  fFitTracklength = FitTracklength_;
  fnHitsMax = 0;
  fnHitsDedx = 0;

  // Initialize fFakeRatio (not provided in constructor parameters)
  fFakeRatio = 0.0;

  fX = x_;
  fY = y_;
  fZ = z_;
  fPx = px_;
  fPy = py_;
  fPz = pz_;
  fdEdx = 0.;
  fDcaX = 0.;
  fDcaY = 0.;
  fDcaZ = 0.;

  fVertexSeedPos.SetXYZ(0.,0.,0.);
  fVertexSeedMom.SetXYZ(0.,0.,0.);
  for(int i = 0; i < 6; ++i) {
    for(int j = 0; j < 6; ++j) {
      fVertexSeedCov[i][j] = 0.;
    }
  }

  /*
  fitTrack = fitTrack_;
  rep = rep_;
  */
   //cout<<" SUCCEED Write"<<endl;
}
void CeeTpcTrack::AddHitId(int value) {
    fHitIdList.push_back(value); // 向 vector 添加数据
}
void CeeTpcTrack::AddHitResidual(TVector3 Res) {
    fResidual.push_back(Res); // 向 vector 添加数据
}

void CeeTpcTrack::SetDca(Double_t x, Double_t y, Double_t z)
{
  fDcaX = x;
  fDcaY = y;
  fDcaZ = z;
}

void CeeTpcTrack::GetDca(Double_t& x, Double_t& y, Double_t& z) const
{
  x = fDcaX;
  y = fDcaY;
  z = fDcaZ;
}

CeeTpcTrack::~CeeTpcTrack(){fHitIdList.clear();}

void CeeTpcTrack::SetVertexSeed(const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov)
{
  fVertexSeedPos = pos;
  fVertexSeedMom = mom;
  int n = std::min(6, cov.GetNrows());
  int m = std::min(6, cov.GetNcols());
  for(int i = 0; i < 6; ++i) {
    for(int j = 0; j < 6; ++j) {
      if(i < n && j < m) {
        fVertexSeedCov[i][j] = cov(i, j);
      } else {
        fVertexSeedCov[i][j] = 0.;
      }
    }
  }
  fHasVertexSeed = kTRUE;
}

void CeeTpcTrack::GetVertexSeedCov(TMatrixDSym& cov) const
{
  if(cov.GetNrows() != 6) {
    cov.ResizeTo(6, 6);
  }
  for(int i = 0; i < 6; ++i) {
    for(int j = 0; j < 6; ++j) {
      cov(i, j) = fVertexSeedCov[i][j];
    }
  }
}

void CeeTpcTrack::GetTrackParameterCovariance(TMatrixDSym& cov) const
{
  cov.ResizeTo(5, 5);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      cov(i, j) = Fitcov[i][j];
    }
  }
}
