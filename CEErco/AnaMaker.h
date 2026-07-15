////////////////////////////////////////////////////////////
// AnaMaker Class: reconstruct lambda candidicates using the output of CEE fast
//	simulation
//	
// 2024.12.2 version 1.0 Xionghong He
////////////////////////////////////////////////////////////

#ifndef CENTRALITY_H
#define CENTRALITY_H

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"

class FairEventHeader;
class CeeMCTrack;
class CeeTpcPoint;
class CeeMWDCPoint;
class CeeiTOFPoint;
class CeeeTOFPoint;
class CeeZdcPoint;
class CeeT0Hit;
class CeeTpcHit;
class CeeiTOFHit;
class CeeEtofHit;
class CeeZdcHit;
class CeeTpcTrack;
class CeeTpcVertex;
class CeeMWDC_Hit;
class CeeMWDC_Track;

class TH1F;
class TH2F;
class TProfile;

class AnaMaker {

public: 

	//// constructor
	AnaMaker(const TString& inFileDir, int fileNum=1, int debug=0);
	//// destructor
	~AnaMaker();

	//// initialize
	bool Init(TString inFileDir, int fileNum); // define in and out file
	void InitialHistgram();
	void InitializeArrays();
	void DebugLevel(int level) {fDebug = level;}

	void SetOutputFile(const TString& file) {fOutFileName = file;}
	void SetEventNumber(int num) {fEvtNo = num;}
	void RotatePi() {fRotPi = true;}
	void RotateProton() {fRotP = true;}

	//// access

	////
	void ProcessData();
	void Finish();

private:
	int fDebug;
	TString fOutFileName;
	TFile* fOutFile;
	TChain* fChain;
	TTree* fRecTree;
	
	Int_t evtId;
	Int_t nLambda;
	Float_t psi;
	Int_t pi_pdg[500];
	Float_t pi_dedx[500];
	Float_t pi_dca[500];
	Float_t pi_px[500], pi_py[500], pi_pz[500];
	Float_t mc_pi_px[500], mc_pi_py[500], mc_pi_pz[500];
	Int_t proton_pdg[500];
	Float_t proton_dedx[500];
	Float_t proton_dca[500];
	Float_t proton_px[500], proton_py[500], proton_pz[500];
	Float_t mc_proton_px[500], mc_proton_py[500], mc_proton_pz[500];
	Float_t ld_mass[500];
	Int_t ld_mctruth[500];
	//add 6/4
	Int_t ld_mctruth_hit[500];
        Float_t ld_dca[500];
	Float_t dca1to2[500];
	Float_t decaylength[500];
	Float_t ld_px[500];
	Float_t ld_py[500];
	Float_t ld_pz[500];
	Float_t mc_ld_px[500];
	Float_t mc_ld_py[500];
	Float_t mc_ld_pz[500];
	Float_t mc_proton_dca[500];
  
  //add 6/8
  Int_t proton_source[500];
  Int_t pi_source[500];
  Int_t has_mwdc[500];
  Float_t proton_theta[500];
  Float_t pi_theta[500];
  //mwdc
  //add 6/10
  Float_t mwdc_true_p_mc_dca[500];
  Float_t mwdc_true_pi_mc_dca[500];
  // MWDC track reco vs MC (proton)
   Float_t mwdc_p_dp[500];       // Δ|p| (reco - MC)
   Float_t mwdc_p_dpt[500];      // ΔpT
   Float_t mwdc_p_dtheta[500];   // Δθ (degree)
   Float_t mwdc_p_dr[500];       // Δpos (cm)
   Float_t mwdc_p_ddca[500];     // ΔDCA = reco DCA - MC DCA
   Float_t mwdc_p_dphi[500];
   // MWDC track reco vs MC (pion)
   Float_t mwdc_pi_dp[500];
   Float_t mwdc_pi_dpt[500];
   Float_t mwdc_pi_dtheta[500];
   Float_t mwdc_pi_dr[500];
   Float_t mwdc_pi_ddca[500];  
   Float_t mwdc_pi_dphi[500];
   
        Float_t ld_y[500];
        //add 2/27
        Int_t p_mcid[500];           
        Int_t p_tof_trkidx[500];     
        Int_t p_tof_pntidx[500];    
        Int_t p_tofmatch[500]; 
        Int_t pi_mcid[500];
        Int_t pi_tof_trkidx[500];
        Int_t pi_tof_pntidx[500];
        Int_t pi_tofmatch[500];
  //add 5/11
  // ======================================================
// eTOF matching diagnostic branches for each reconstructed Lambda
// ======================================================
  Float_t p_true_etof_dcaXY[500];
  Float_t pi_true_etof_dcaXY[500];

  Float_t p_nearest_etof_dcaXY[500];
  Float_t pi_nearest_etof_dcaXY[500];

  Int_t p_nearest_etof_pdg[500];
  Int_t pi_nearest_etof_pdg[500];

  Int_t p_nearest_etof_same[500];
  Int_t pi_nearest_etof_same[500];
  
  // ======================================================
// iTOF matching diagnostic branches for each reconstructed Lambda
// ======================================================
  Float_t p_true_itof_dcaXY[500];
  Float_t pi_true_itof_dcaXY[500];

  Float_t p_nearest_itof_dcaXY[500];
  Float_t pi_nearest_itof_dcaXY[500];

  Int_t p_nearest_itof_pdg[500];
  Int_t pi_nearest_itof_pdg[500];

  Int_t p_nearest_itof_same[500];
  Int_t pi_nearest_itof_same[500];

  //add 5/13
  // ======================================================
   // TPC projection residual at eTOF z plane
   // dx/dy/dr = projected TPC position - true-linked eTOF hit position
   // ======================================================
   Float_t p_true_etof_proj_dx[500];
   Float_t p_true_etof_proj_dy[500];
   Float_t p_true_etof_proj_dr[500];
  
   Float_t pi_true_etof_proj_dx[500];
   Float_t pi_true_etof_proj_dy[500];
   Float_t pi_true_etof_proj_dr[500]; 
 
	TH2F* mc_pty;
	TH1F* mc_pt;
	TProfile* mc_v1;
	TProfile* rec_v1;
	TH1F* itofTime_diff;
	TH1F* itofLength_diff;
	TH1F* etofTime_diff;
        TH1F* etofLength_diff;
        //add
        TH1F* itofLength_diff2;
  //Add 5/9
  // ======================================================
// Lambda daughter TPC/eTOF link diagnostic histograms
// ======================================================
  TH1F* hLinkSummary;

  TH1F* hEvtMcLambdaPPi;
  TH1F* hEvtBothTpc;
  TH1F* hEvtBothEtof;
  TH1F* hEvtBothTpcEtof;

  TH1F* hEvtRecoTrueCand;
  TH1F* hEvtRecoTrueUnique;

  TH2F* hEvtRecoUniqueVsBothTpcEtof;

  TH1F* hProtonDauTpcMult;
  TH1F* hPionDauTpcMult;
  TH1F* hProtonDauEtofMult;
  TH1F* hPionDauEtofMult;

  TH2F* hLambdaDaughterDetStatus;
  //====
  //add 5/11
  // ======================================================
// DcaXY and nearest-hit matching diagnostic histograms
// ======================================================
  TH1F* hTruthDetect_pTrueDcaXY;
  TH1F* hTruthDetect_piTrueDcaXY;
  TH1F* hTruthDetect_pNearestDcaXY;
  TH1F* hTruthDetect_piNearestDcaXY;

  TH1F* hTruthDetect_pNearestPdg;
  TH1F* hTruthDetect_piNearestPdg;
  TH1F* hTruthDetect_pNearestWrongPdg;
  TH1F* hTruthDetect_piNearestWrongPdg;

  TH1F* hTruthDetectNearestSummary;
  TH1F* hTruthDetectToFMatchCutSummary;

  TH1F* hRecoTrue_pTrueDcaXY;
  TH1F* hRecoTrue_piTrueDcaXY;
  TH1F* hRecoTrue_pNearestDcaXY;
  TH1F* hRecoTrue_piNearestDcaXY;

  TH1F* hRecoTrue_pNearestPdg;
  TH1F* hRecoTrue_piNearestPdg;
  
  
	TH1F* tofdca;
	TH1F* mc_tofdca;
  //add
  TH1F* mc_dndy;
  TH1F* rec_dndy_all;
  TH1F* rec_dndy_mctruth;
  TH1F* rec_dndy_mctruth_hit;  
 
  //add 5/11
  // ======================================================
// Lambda daughter TPC/iTOF link diagnostic histograms
// ======================================================
TH1F* hLinkSummary_iTOF;

TH1F* hEvtMcLambdaPPi_iTOF;
TH1F* hEvtBothTpc_iTOF;
TH1F* hEvtBothiTof_iTOF;
TH1F* hEvtBothTpcItof_iTOF;

TH1F* hEvtRecoTrueCand_iTOF;
TH1F* hEvtRecoTrueUnique_iTOF;

TH2F* hEvtRecoUniqueVsBothTpcItof_iTOF;

TH1F* hProtonDauTpcMult_iTOF;
TH1F* hPionDauTpcMult_iTOF;
TH1F* hProtonDauiTofMult_iTOF;
TH1F* hPionDauiTofMult_iTOF;

TH2F* hLambdaDaughterDetStatus_iTOF;

// ======================================================
// iTOF DcaXY and nearest-hit matching diagnostic histograms
// ======================================================
TH1F* hTruthDetect_pTrueDcaXY_iTOF;
TH1F* hTruthDetect_piTrueDcaXY_iTOF;
TH1F* hTruthDetect_pNearestDcaXY_iTOF;
TH1F* hTruthDetect_piNearestDcaXY_iTOF;

TH1F* hTruthDetect_pNearestPdg_iTOF;
TH1F* hTruthDetect_piNearestPdg_iTOF;
TH1F* hTruthDetect_pNearestWrongPdg_iTOF;
TH1F* hTruthDetect_piNearestWrongPdg_iTOF;

TH1F* hTruthDetectNearestSummary_iTOF;
TH1F* hTruthDetectToFMatchCutSummary_iTOF;

TH1F* hRecoTrue_pTrueDcaXY_iTOF;
TH1F* hRecoTrue_piTrueDcaXY_iTOF;
TH1F* hRecoTrue_pNearestDcaXY_iTOF;
TH1F* hRecoTrue_piNearestDcaXY_iTOF;

TH1F* hRecoTrue_pNearestPdg_iTOF;
TH1F* hRecoTrue_piNearestPdg_iTOF;

  //add 5/13
  // ======================================================
// TPC projection residual histograms at eTOF z plane
// ======================================================
TH1F* hTruthDetect_pProjDx;
TH1F* hTruthDetect_pProjDy;
TH1F* hTruthDetect_pProjDr;

TH1F* hTruthDetect_piProjDx;
TH1F* hTruthDetect_piProjDy;
TH1F* hTruthDetect_piProjDr;

TH2F* hTruthDetect_pProjDxDy;
TH2F* hTruthDetect_piProjDxDy;

TH1F* hRecoTrue_pProjDx;
TH1F* hRecoTrue_pProjDy;
TH1F* hRecoTrue_pProjDr;

TH1F* hRecoTrue_piProjDx;
TH1F* hRecoTrue_piProjDy;
TH1F* hRecoTrue_piProjDr;

TH2F* hRecoTrue_pProjDxDy;
TH2F* hRecoTrue_piProjDxDy;


// ======================================================
// theta acceptance / efficiency diagnostic histograms
// theta is MC daughter polar angle relative to +z beam axis, in degree
// ======================================================
TH1F* hThetaProtonAll;
TH1F* hThetaProtonTpc;
TH1F* hThetaProtonEtof;
TH1F* hThetaProtonTpcEtof;

TH1F* hThetaPionAll;
TH1F* hThetaPionTpc;
TH1F* hThetaPionEtof;
TH1F* hThetaPionTpcEtof;

TH1F* hEffThetaProtonTpc;
TH1F* hEffThetaProtonEtof;
TH1F* hEffThetaProtonTpcEtof;

TH1F* hEffThetaPionTpc;
TH1F* hEffThetaPionEtof;
TH1F* hEffThetaPionTpcEtof;
  
    
	TH1F* mc_p_dca;
	TH1F* mc_pi_dca;
	TH1F* itofToT[24];
	TH1F* itofTime[24];
	TH1F* it0Time;
	TH2F* h2_dEdx_beta;
	TH2F* h2_dEdx_betagamma;
	TH2F* h2_time_length;
	//add 2/3
  TH1F* itof_dz;
  TH1F* itof_dr;
  TH1F* itof_dL_line;
  TH1F* itof_dL_mcdef;
  TH2F* itof_dLold_vs_dz;
  TH2F* itof_dLline_vs_dz;
  TH1F* h_pTofMis_hitPdg; 
  TH1F* h_pTofMis_hitMomPdg;
  TH1F* h_pTofMis_dca;
  TH2F* h_pTofMis_dcaVsMass; 
  //add 3/27
  TH1F* tofdca2;
  TH1F* mc_tofdca2;  
//  TH1F* tofdca2 = nullptr;
//  TH1F* mc_tofdca2 = nullptr;  
  //add 3/3
  TH1F* h_piTofMis_hitPdg    = nullptr;
  TH1F* h_piTofMis_hitMomPdg = nullptr;
  TH1F* h_piTofMis_dca       = nullptr;

  TH1F* h_pTofMis_hitPdg_mWin    = nullptr;
  TH1F* h_pTofMis_hitMomPdg_mWin = nullptr;

  TH1F* h_piTofMis_hitPdg_mWin    = nullptr;
  TH1F* h_piTofMis_hitMomPdg_mWin = nullptr;
  TH1F* h_mass_true_misTof = nullptr;
  
	FairEventHeader* iEventHeader;
	CeeMCTrack*    iMCTrack;
	CeeMCTrack*    fProtonMCTrack;
	CeeMCTrack*    fPionMCTrack;
	CeeTpcPoint*   iTpcPoint;
	CeeMWDCPoint*  iMWDCPoint;
	CeeiTOFPoint*  iiTOFPoint;
	CeeeTOFPoint*  ieTOFPoint;
	CeeZdcPoint*   iZdcPoint;
	CeeT0Hit*      iT0Hit;
	CeeTpcHit*     iTpcHit;
	CeeEtofHit*    iEtofHit;
	CeeZdcHit*     iZdcHit;
	CeeTpcVertex*  iTpcVertex;
	CeeMWDC_Hit*   iMWDC_Hit;
	CeeMWDC_Track* iMWDC_Track;

  CeeTpcTrack*   fProtonTrack;
  CeeTpcTrack*   fPionTrack;
  CeeiTOFHit*    fProtonTofHit;
  CeeiTOFHit*    fPionTofHit;

//  CeeEtofHit*    fProtonTofHit;
//  CeeEtofHit*    fPionTofHit;
	TClonesArray*     fMCTrack;
	TClonesArray*     fTpcPoint;
	TClonesArray*     fMWDCPoint;
	TClonesArray*     fiTOFPoint;
	TClonesArray*     feTOFPoint;
	TClonesArray*     fZdcPoint;
	TClonesArray*     fT0Hit;
	TClonesArray*     fTpcHit;
	TClonesArray*     fiTOFHit;
	TClonesArray*     fEtofHit;
	TClonesArray*     fZdcHit;
	TClonesArray*     fTpcTrack;
	TClonesArray*     fTpcVertex;
	TClonesArray*     fMWDC_Hit;
	TClonesArray*     fMWDC_Track;

	Int_t fEvtNo;
	Bool_t fRotPi;
	Bool_t fRotP;
	// event vertex
	TVector3 fEvtVtx;

	Double_t GetDca(const TVector3& pos, const TVector3& p4, const TVector3& vertex);
	Double_t GetDcaXY(const TVector3& pos, const TVector3& p4, const TVector3& vertex);
};

#endif
