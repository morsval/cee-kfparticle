// $Id: StKFParticleAnalysisMaker.h,v 1.16 2014/08/06 11:43:53 jeromel Exp $
/*!
 * \class  StKFParticleAnalysisMaker
 * \author Maksym Zyzak
 * \date   2017/10/17
 * \brief  class for analysis of PicoDst
 */                                                                      
#ifndef STAR_StKFParticleAnalysisMaker
#define STAR_StKFParticleAnalysisMaker
//#define __DEVT__
#ifndef StMaker_H
#include "StMaker.h"
#endif
#include "TMVA/Reader.h"

#include "TString.h"
#include <TTree.h>
#include "StThreeVectorF.hh"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLorentzVector.h"

#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdEpInfo.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "TClonesArray.h"
#include "TObject.h"
//#include "StRoot/StPileupUtil/StPileupUtil.h"

#define MaxNumberOfTH1F      200
#define MaxNumberOfTH2F      100
#define MaxNumberOfTH3F      100
#define MaxNumberOfTProfile  200
#define MAXMEC		     20000


class StKFParticleInterface;
class StKFParticlePerformanceInterface;
class KFParticle;
class StPicoDst;
class StMuDst;
class TNtuple;
class TFile;
class TChain;
class StRefMultCorr;


class StPicoDstMaker;
class StPicoDst;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TProfile2D;

class TNtuple;
class TFile;
class TChain;
class StPicoEvent;

class StEpdEpFinder;
class StPicoEpdHit;
class TClonesArray;
class MyStarEvent;

class StKFParticleAnalysisMaker : public StMaker {
 private:
  static const int fNNTuples = 8;
  Char_t                mBeg[1];        //!
  StMuDst                          *fMuDst;
  StPicoDstMaker                   *mPicoDstMaker ;                    //  Make PicoDst pointer available to member functions
  StPicoDst                        *mPicoDst;                          //!
  StKFParticleInterface            *mStKFParticleInterface;            //!
  StKFParticlePerformanceInterface *mStKFParticlePerformanceInterface; //!

  TH1F* histogram[MaxNumberOfTH1F]   ;     //  1D Histograms


	TH3F*	h3f_tpc_RunPhiEta_pos; 
	TH3F*	h3f_tpc_RunPhiEta_neg; 

	TH3F*	h3f_tpc_PtEtaPhi_pos; 
	TH3F*	h3f_tpc_PtEtaPhi_neg; 

	TH1F* hCent;
	TH1F* hVz_before_cut;
	TH1F* hVz_after_trigger_cut;
	TH1F* hVz_after_vertex_cut;
	TH1F* hVz_after_centrality_cut;
	TH1F* hVz_after_badrun_cut;
	TH1F* hVz_after_processEvt_cut;



  TString mOutputFileName;

  TFile* fOutputFile;
  TTree*   fEventTree;
  MyStarEvent* fMyEvent;
  //Subhash
  
  Int_t runnumber;
  Int_t runPointer;
  Int_t nBTofmatch;
  Int_t cent;
  Int_t CurrentEvent_centrality;
  float Pi;
  float twoPi;
 

 
  //Subhash TPC EP related
  Int_t           runnumber_flag;
  Int_t                cent_flag;
  Int_t         runnumberPointer; 
  int centIdx;
  
  //Subhash EPD EP related
  StEpdEpFinder* mEpFinder;
  StEpdGeom* mEpdGeom;

  TH1F *dndy_epd_full;
  TH1F *dndy_epd_cut;

  ULong_t mEventsStarted               ;     //  Number of Events read
  ULong_t mEventsProcessed             ;     //  Number of Events processed and analyze
  Char_t                mEnd[1];        //!
 
  StRefMultCorr *mRefmultCorrUtil;
  
 public: 
  StKFParticleAnalysisMaker(const char *name="KFParticleAnalysis");
  virtual       ~StKFParticleAnalysisMaker();
  virtual Int_t  Init();
  virtual Int_t  InitRun(Int_t runumber);
  void           BookVertexPlots();
  virtual Int_t  Make();
  virtual Int_t  Finish();
  static void    PrintMem(const Char_t *opt = "");

  void BookTree();
  bool AcceptTrack(StPicoTrack*)   ;     //  Function to make cuts on track quality (ie. nhits or eta)

  bool ProtonCheck(Float_t nsigma);
  bool ProtonMCheck(Float_t beta,Float_t p);
  bool ProtonPtCut(float ProtonPt, bool isProton);
  bool PionCheck(Float_t nsigma);
  bool PionMCheck(Float_t beta,Float_t p);

  bool V0DirectionCut(Float_t pLxLxPV); 
  bool ProtonDcaCut(Float_t dca, Float_t pt); 
  bool PionDcaCut(Float_t dca, Float_t pt); 
  bool DaughterDcaCut(Float_t dca, Float_t pt); 
  bool LambdaDcaCut(Float_t dca, Float_t pt); 
  bool LambdaDecayLengthCut(Float_t decaylength, Float_t pt); 
  bool LambdaMassCut(Float_t LambdaMass); 
  bool XiMassCut(Float_t XiMass); 
  bool OmegaMassCut(Float_t XiMass);


  bool PassTpcSelection(Float_t eta_min, Float_t eta_max, StPicoTrack* , TVector3 pV);
  bool PassEpdSelection(Int_t row, Short_t side, Float_t nMip);
  bool PassLambdaSelection(double ld_chi2topo, double ld_chi2ndf, double ld_l, double chi2primary_p, double chi2primary_pi,double ld_mass);

  bool AcceptEvent   (StPicoEvent*)   ;     //  Function to make cuts on event quanitites (ie. vertex position)
  bool AcceptTrigger (StPicoEvent*)   ;     //  Function to make cuts on the trigger words (Cusomize each year)

  double CalculateHelicity(TLorentzVector PLam4, TLorentzVector PPro4);
  double CalculatePolarization(bool isLambda, double Psi1, TLorentzVector PLam4, TLorentzVector PPro4);
  double CalculateXiPolarization(bool isXi, double Psi1, TLorentzVector PXi4, TLorentzVector PLambda4);
  TVector3 ConvertToTV3(StThreeVectorF inVector);

  //Subhash
  Int_t MakeRunPointer(int frunnumber);//required for QA hist
  Int_t Centrality(int gRefMult);
  Double_t FitWeight(Double_t refMult);
  //StPileupUtil *mPileupTool;
  Float_t reweight;
  
  ClassDef(StKFParticleAnalysisMaker,0)   //
};
#endif
// $Log: StKFParticleAnalysisMaker.h,v $
