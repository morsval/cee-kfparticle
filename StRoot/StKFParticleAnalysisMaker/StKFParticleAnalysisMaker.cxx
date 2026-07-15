//*-- Author : Yuri Fisyak 02/02/2016
#include "StKFParticleAnalysisMaker.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TTree.h"
//--- KF particle classes ---
#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "KFParticleTopoReconstructor.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"
//--- Pico classes ---
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoArrays.h"
//--- Mu classes ---
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
//--- TMVA classes ---
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
//--- StRefMult class ---
// #include "StRefMultCorr/StRefMultCorr.h"
// #include "StRefMultCorr/CentralityMaker.h"
//--- EP ep finder class ---
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "TObjArray.h"

#include "MyEvent.h"
#include "run_3.h" //==14.6 GeV, 2023.06.25
// #include "StRoot/StPileupUtil/StPileupUtil.h"
const float fxt_y_mid = -1.045;
const int totalRunNumber = 191;

#define OmegaPdgMass 1.67245
#define OmegaPdg 3334
#define AntiOmegaPdg -3334
#define KaonPdgMass 0.493677
#define KaonMinusPdg -321
#define KaonPlusPdg 321
#define XiPdgMass 1.32171
#define XiPdg 3312
#define AntiXiPdg -3312
#define LambdaPdgMass 1.11568
#define ProtonPdgMass 0.938272
#define PionPdgMass 0.139570
#define LambdaPdg 3122
#define ProtonPdg 2212
#define PiMinusPdg -211
#define AntiLambdaPdg -3122
#define AntiProtonPdg -2212
#define PiPlusPdg 211

ClassImp(StKFParticleAnalysisMaker);

//________________________________________________________________________________
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name) : StMaker(name)
{
  memset(mBeg, 0, mEnd - mBeg + 1);

  mOutputFileName = name;
  fEventTree = 0;
  // mPicoDstMaker       =  maker ;
  mEventsStarted = 0;
  mEventsProcessed = 0; // Zero the Number of Events processed by the maker
}
//________________________________________________________________________________
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker()
{
  SafeDelete(mStKFParticleInterface);
  SafeDelete(mStKFParticlePerformanceInterface);
}

//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Init()
{

  runnumber = -13100001;
  Pi = 3.14159;
  twoPi = 6.28318;
  runnumber_flag = totalRunNumber;
  cent_flag = 10;
  runnumberPointer = 0;
  mOutputFileName += ".root";
  fOutputFile = new TFile(mOutputFileName,"recreate");
  fMyEvent = new MyStarEvent;
	fEventTree = new TTree("EventTree","Event");
  fEventTree->Branch("fMyEvent",&fMyEvent,100000,1);
 // fEventTree->Branch("myEvent",&fMyEvent);
	fMyEvent->fMyTpcTrackRefA.reserve(3000);
	fMyEvent->fMyTpcTrackRefB.reserve(3000);
 	fMyEvent->fMyLambda.reserve(3000);
//	fMyEvent->fMyEpdHit.reserve(50000);
	
	mStKFParticleInterface = new StKFParticleInterface;
	bool storeMCHistograms = false;
  mStKFParticlePerformanceInterface = new StKFParticlePerformanceInterface(mStKFParticleInterface->GetTopoReconstructor(), storeMCHistograms);  

	vector<float> TpcHistogramPtBinEdges = {
		0.15, 0.153773, 0.157681, 0.161729, 0.165922, 0.170266, 0.174765, 0.179426, 0.184254, 0.189256, 0.194437, 0.199803, 0.205362, 0.21112, 0.217085, 0.223264, 0.229665, 0.236295, 
		0.243163, 0.250277, 0.257647, 0.26528, 0.273188, 0.281379, 0.289864, 0.298654, 0.307758, 0.31719, 0.326959, 0.337079, 0.347562, 0.358421, 0.369669, 0.381321, 0.393391, 0.405894, 
		0.418845, 0.432261, 0.446158, 0.460553, 0.475465, 0.490912, 0.506912, 0.523487, 0.540656, 0.558441, 0.576864, 0.595948, 0.615716, 0.636193, 0.657405, 0.679377, 0.702138, 0.725715, 
		0.750138, 0.775437, 0.801643, 0.828789, 0.856909, 0.886038, 0.916211, 0.947467, 0.979843, 1.01338, 1.04812, 1.08411, 1.12139, 1.16, 1.2};
	vector<float> TpcHistogramEtaBinEdges; for( float EtaBinEdge=-1.5; EtaBinEdge<=0.; EtaBinEdge+=0.02 )TpcHistogramEtaBinEdges.push_back(EtaBinEdge);
	vector<float> TpcHistogramPhiBinEdges; for( float PhiBinEdge=-TMath::Pi(); PhiBinEdge<=TMath::Pi(); PhiBinEdge+=2.*TMath::Pi()/200. )TpcHistogramPhiBinEdges.push_back(PhiBinEdge);

	h3f_tpc_PtEtaPhi_pos =  new TH3F("h3f_tpc_PtEtaPhi_pos", "", TpcHistogramPtBinEdges.size()-1,&TpcHistogramPtBinEdges[0],TpcHistogramEtaBinEdges.size()-1,&TpcHistogramEtaBinEdges[0],TpcHistogramPhiBinEdges.size()-1,&TpcHistogramPhiBinEdges[0]); 
	h3f_tpc_PtEtaPhi_neg =  new TH3F("h3f_tpc_PtEtaPhi_neg", "", TpcHistogramPtBinEdges.size()-1,&TpcHistogramPtBinEdges[0],TpcHistogramEtaBinEdges.size()-1,&TpcHistogramEtaBinEdges[0],TpcHistogramPhiBinEdges.size()-1,&TpcHistogramPhiBinEdges[0]); 	
//
//	h3f_tpc_PtEtaPhi_pos  =  new TH3F("h3f_tpc_PtEtaPhi_pos", "", 105, 0.15, 1.2, 100, -2., 0, 200,-1.*TMath::Pi(),TMath::Pi() );
//	h3f_tpc_PtEtaPhi_neg  =  new TH3F("h3f_tpc_PtEtaPhi_neg", "", 105, 0.15, 1.2, 100, -2., 0, 200,-1.*TMath::Pi(),TMath::Pi() );

	hVz_before_cut          = new TH1F("hVz_before_cut","hVz_before_cut;Vz;counts",60,197,203);
	hVz_after_trigger_cut   = new TH1F("hVz_after_trigger_cut","hVz_after_trigger_cut;Vz;counts",60,197,203);
	hVz_after_vertex_cut    = new TH1F("hVz_after_vertex_cut","hVz_after_vertex_cut;Vz;counts",60,197,203);
	hVz_after_centrality_cut= new TH1F("hVz_after_centrality_cut","hVz_after_centrality_cut;Vz;counts",60,197,203);
	hVz_after_badrun_cut    = new TH1F("hVz_after_badrun_cut","hVz_after_badrun_cut;Vz;counts",60,197,203);

	hVz_after_processEvt_cut= new TH1F("hVz_after_processEvt_cut","hVz_after_centrality_cut;Vz;counts",60,197,203);
	hCent = new TH1F("hCent",";centrality;Counts",15,-2,13);
  
	cout<<"extremely good"<<endl;    
	return kStOK;
}
//________________________________________________________________________________
Int_t StKFParticleAnalysisMaker::InitRun(Int_t runumber)
{
  return StMaker::InitRun(runumber);
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::PrintMem(const Char_t *opt)
{
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  cout << opt
       << "\tMemory : Total = " << info.fMemTotal
       << "\tUsed = " << info.fMemUsed
       << "\tFree = " << info.fMemFree
       << "\tSwap Total = " << info.fSwapTotal
       << "\tUsed = " << info.fSwapUsed
       << "\tFree = " << info.fSwapFree << endl;
}

//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Make()
{

  mPicoDst = StPicoDst::instance();
  if (!mPicoDst)
    return kStOK;
  StPicoEvent *picoEvent = mPicoDst->event();
  if (!picoEvent)
    return kStOK;

	int RunId =  picoEvent->runId();
	int RunYear =  floor(RunId / pow(10, 6)); 
  fMyEvent->mRunId = picoEvent->runId();

  fMyEvent->mRunYear = floor(RunId / pow(10, 6)); // actual year data was taken + 1
  fMyEvent->mRunDay  = floor((RunId - RunYear * pow(10, 6)) / pow(10, 3));
  fMyEvent->mEventId = picoEvent->eventId();

  TVector3 PrimaryVertex = picoEvent->primaryVertex();
  fMyEvent->mPVx = PrimaryVertex.X();
  fMyEvent->mPVy = PrimaryVertex.Y();
  fMyEvent->mPVz = PrimaryVertex.Z();

//  Int_t nGlobalTracks = picoEvent->numberOfGlobalTracks();
  Double_t Vx = PrimaryVertex.X();
  Double_t Vy = PrimaryVertex.Y();
  Double_t Vz = PrimaryVertex.Z();

  if ((Vx < 1.0e-5 && Vx > -1.0e-5) && (Vy < 1.0e-5 && Vy > -1.0e-5) && (Vz < 1.0e-5 && Vz > -1.0e-5))
    return kStOK; //== added by q.hu on 2023.06.25

  mEventsStarted++;

	hVz_before_cut->Fill(Vz);

  if (!AcceptTrigger(picoEvent))
    return kStOK; // Skip this event if it doens't pass the trigger cuts
  
	hVz_after_trigger_cut->Fill(Vz);

	if (!AcceptEvent(picoEvent))
    return kStOK; // Skip this event if it doens't pass the 'event' cuts
	
	hVz_after_vertex_cut->Fill(Vz);
  //-------------------------------------------------------------------------
  // Subhash
  //-------------------------------------------------------------------------
  runnumber = picoEvent->runId();
  runPointer = -999;
  runPointer = MakeRunPointer(runnumber);
  if (runPointer == -999)
    return kStOK;
  runnumberPointer = runPointer;
	
	hVz_after_badrun_cut->Fill(Vz);

  nBTofmatch = picoEvent->nBTOFMatch(); // not been used in code

  int countrefmult = 0;
  int nTracks = mPicoDst->numberOfTracks();
  for (Int_t itr = 0; itr < nTracks; itr++)
  {
    const StPicoTrack *ptrk = (StPicoTrack *)mPicoDst->track(itr);
    if (!ptrk)
      continue;
    if (!ptrk->isPrimary())
      continue;
    countrefmult++;
  }

	if(countrefmult>195) return kStOK;

  int cent = Centrality(countrefmult);
  fMyEvent->mCentrality = cent;
  fMyEvent->mRefMult = countrefmult;
  CurrentEvent_centrality = cent + 1;
  reweight = FitWeight(countrefmult);
  
	hCent->Fill(cent);	
	if (cent > 8 || cent < 0) return kStOK;

	hVz_after_centrality_cut->Fill(Vz);

  // Define TPC subA and subB
  float eta_subA_min = -0.5;
  float eta_subA_max = -0.4;
  float eta_subB_min = -0.2;
  float eta_subB_max = -0.1;

  int nTrackRefA=0;
  int nTrackRefB=0;

	double Qx_raw_1st_tpcA=0;
	double Qx_raw_1st_tpcB=0;
	double Qy_raw_1st_tpcA=0;
	double Qy_raw_1st_tpcB=0;

	double Qx_phi_1st_tpcA=0;
	double Qx_phi_1st_tpcB=0;
	double Qy_phi_1st_tpcA=0;
	double Qy_phi_1st_tpcB=0;


	double Qx_raw_2nd_tpcA=0;
	double Qx_raw_2nd_tpcB=0;
	double Qy_raw_2nd_tpcA=0;
	double Qy_raw_2nd_tpcB=0;

	double Qx_phi_2nd_tpcA=0;
	double Qx_phi_2nd_tpcB=0;
	double Qy_phi_2nd_tpcA=0;
	double Qy_phi_2nd_tpcB=0;




  for(int iTrack = 0; iTrack<nTracks; iTrack++){
      StPicoTrack *tpcTrack = (StPicoTrack *)mPicoDst->track(iTrack);
			
			if(PassTpcSelection(eta_subA_min, eta_subA_max, tpcTrack, PrimaryVertex)){
        fMyEvent->fMyTpcTrackRefA[nTrackRefA].mCharge   = tpcTrack->charge();
        fMyEvent->fMyTpcTrackRefA[nTrackRefA].mPt       = tpcTrack->pMom().Perp();
        fMyEvent->fMyTpcTrackRefA[nTrackRefA].mEta      = tpcTrack->pMom().PseudoRapidity();
        fMyEvent->fMyTpcTrackRefA[nTrackRefA].mPhi      = tpcTrack->pMom().Phi();
        fMyEvent->fMyTpcTrackRefA[nTrackRefA].nHitsFit  = tpcTrack->nHitsFit();
        fMyEvent->fMyTpcTrackRefA[nTrackRefA].nHitsRatio= (float)tpcTrack->nHitsFit()/(float)tpcTrack->nHitsMax();
			//		double mmPt  = tpcTrack->pMom().Perp();
			//		double mmPhi = tpcTrack->pMom().Phi();
			//		Qx_raw_1st_tpcA += cos(1.*mmPhi) ;	
			//		Qx_phi_1st_tpcA += mmPt*cos(1.*mmPhi);	

			//		Qy_raw_1st_tpcA += sin(1.*mmPhi) ;	
			//		Qy_phi_1st_tpcA += mmPt*sin(1.*mmPhi);	
	
			//		Qx_raw_2nd_tpcA += cos(2.*mmPhi) ;	
			//		Qx_phi_2nd_tpcA += mmPt*cos(2.*mmPhi);	
      //           
			//		Qy_raw_2nd_tpcA += sin(2.*mmPhi) ;	
			//		Qy_phi_2nd_tpcA += mmPt*sin(2.*mmPhi);	
			//	
					nTrackRefA++;
			}//TPC refA

      if(PassTpcSelection(eta_subB_min, eta_subB_max, tpcTrack, PrimaryVertex)){
        fMyEvent->fMyTpcTrackRefB[nTrackRefB].mCharge   = tpcTrack->charge();
        fMyEvent->fMyTpcTrackRefB[nTrackRefB].mPt       = tpcTrack->pMom().Perp();;
        fMyEvent->fMyTpcTrackRefB[nTrackRefB].mEta      = tpcTrack->pMom().PseudoRapidity();
        fMyEvent->fMyTpcTrackRefB[nTrackRefB].mPhi      = tpcTrack->pMom().Phi();
        fMyEvent->fMyTpcTrackRefB[nTrackRefB].nHitsFit  = tpcTrack->nHitsFit();
        fMyEvent->fMyTpcTrackRefB[nTrackRefB].nHitsRatio= (float)tpcTrack->nHitsFit()/(float)tpcTrack->nHitsMax();
 			//		double mmPt  = tpcTrack->pMom().Perp();
			//		double mmPhi = tpcTrack->pMom().Phi();
			//		Qx_raw_1st_tpcB += cos(1.*mmPhi) ;	
			//		Qx_phi_1st_tpcB += mmPt*cos(1.*mmPhi);	

			//		Qy_raw_1st_tpcB += sin(1.*mmPhi) ;	
			//		Qy_phi_1st_tpcB += mmPt*sin(1.*mmPhi);	
	
			//		Qx_raw_2nd_tpcB += cos(2.*mmPhi) ;	
			//		Qx_phi_2nd_tpcB += mmPt*cos(2.*mmPhi);	
      //           
			//		Qy_raw_2nd_tpcB += sin(2.*mmPhi) ;	
			//		Qy_phi_2nd_tpcB += mmPt*sin(2.*mmPhi);	
	
					nTrackRefB++;
      }//TPC refB
		


  }// loop TPC tracks
			//if(mEventsProcessed == 7542 ) cout<<"nTrackRefA="<<nTrackRefA<<endl;
//	cout<<"nTrackRefB="<<nTrackRefB<<endl;
  fMyEvent->nTpcRefA = nTrackRefA;
  fMyEvent->nTpcRefB = nTrackRefB;
  
//	fMyEvent->QxTpcARaw1st = Qx_raw_1st_tpcA;
//	fMyEvent->QyTpcARaw1st = Qy_raw_1st_tpcA;
//	fMyEvent->QxTpcBRaw1st = Qx_raw_1st_tpcB;
//	fMyEvent->QyTpcBRaw1st = Qy_raw_1st_tpcB;
//
//	fMyEvent->QxTpcAPhi1st = Qx_phi_1st_tpcA;
//	fMyEvent->QyTpcAPhi1st = Qy_phi_1st_tpcA;
//	fMyEvent->QxTpcBPhi1st = Qx_phi_1st_tpcB;
//	fMyEvent->QyTpcBPhi1st = Qy_phi_1st_tpcB;
//	
//
//	fMyEvent->QxTpcARaw2nd = Qx_raw_2nd_tpcA;
//	fMyEvent->QyTpcARaw2nd = Qy_raw_2nd_tpcA;
//	fMyEvent->QxTpcBRaw2nd = Qx_raw_2nd_tpcB;
//	fMyEvent->QyTpcBRaw2nd = Qy_raw_2nd_tpcB;
//
//	fMyEvent->QxTpcAPhi2nd = Qx_phi_2nd_tpcA;
//	fMyEvent->QyTpcAPhi2nd = Qy_phi_2nd_tpcA;
//	fMyEvent->QxTpcBPhi2nd = Qx_phi_2nd_tpcB;
//	fMyEvent->QyTpcBPhi2nd = Qy_phi_2nd_tpcB;
//	

	fMyEvent->fMyTpcTrackRefA.resize(nTrackRefA);
  fMyEvent->fMyTpcTrackRefB.resize(nTrackRefB);

  //-------------------------------------------------------------------------------
  // EPD for FXT 3 GeV
  StEpdGeom my_EpdGeom;
  double epd_eta, epd_phi, n_Mip, TileWeight;
  int TT, PP, ring, tileId;

	double Qx_raw_1st_epd=0; 
	double Qy_raw_1st_epd=0;

	double Qx_phi_1st_epd=0;
	double Qy_phi_1st_epd=0;

	double Qx_raw_2nd_epd=0; 
	double Qy_raw_2nd_epd=0;

	double Qx_phi_2nd_epd=0;
	double Qy_phi_2nd_epd=0;

  int nHitsEPD=0;
  for (int iEpdHit = 0; iEpdHit < (int)mPicoDst->numberOfEpdHits(); iEpdHit++)
  {
    StPicoEpdHit EpdHit = *mPicoDst->epdHit(iEpdHit);
    epd_phi = (my_EpdGeom.RandomPointOnTile(EpdHit.id()) - PrimaryVertex).Phi();
    tileId = EpdHit.id();
    epd_eta = (my_EpdGeom.RandomPointOnTile(tileId) - PrimaryVertex).Eta();
    //dndy_epd_full->Fill(epd_eta);

    if(!PassEpdSelection(EpdHit.row(), EpdHit.side(), EpdHit.nMIP())) continue;

    //dndy_epd_cut->Fill(epd_eta);

    n_Mip = EpdHit.nMIP();
    TileWeight = (n_Mip > 2) ? 2 : n_Mip;
    TT = EpdHit.tile();
    PP = EpdHit.position();
    ring = EpdHit.row();

		Qx_raw_1st_epd +=  cos(1.0 * epd_phi);
		Qy_raw_1st_epd +=  sin(1.0 * epd_phi);
		
		Qx_phi_1st_epd +=  cos(1.0 * epd_phi);
		Qy_phi_1st_epd +=  sin(1.0 * epd_phi);
	

		Qx_raw_2nd_epd +=  cos(2.0 * epd_phi);
		Qy_raw_2nd_epd +=  sin(2.0 * epd_phi);
		
		Qx_phi_2nd_epd +=  cos(2.0 * epd_phi);
		Qy_phi_2nd_epd +=  sin(2.0 * epd_phi);

	//if(mEventsProcessed == 7542 ) cout<<setw(10) <<"Events ID="<<setw(10) << picoEvent->eventId()<<setw(10) <<" Epd phi "<< setw(10) << epd_phi<<setw(10) <<endl;

    //fMyEvent->fMyEpdHit[nHitsEPD].mEpdPhi  = epd_phi;
    //fMyEvent->fMyEpdHit[nHitsEPD].mEpdEta  = epd_eta;
    //fMyEvent->fMyEpdHit[nHitsEPD].mEpdTile = TT;
    //fMyEvent->fMyEpdHit[nHitsEPD].mEpdTileWeight = TileWeight;
    //fMyEvent->fMyEpdHit[nHitsEPD].mEpdPosition   = PP;
    //fMyEvent->fMyEpdHit[nHitsEPD].mEpdRing       = ring;
    //fMyEvent->fMyEpdHit[nHitsEPD].mEpdnMIP       = n_Mip;
		
		nHitsEPD++;
  } // iEpdHit
//			if(mEventsProcessed == 7542 ) cout<<"EPD hits number = "<<nHitsEPD <<endl;
  fMyEvent->nEpdHit = nHitsEPD;
	
	fMyEvent->QxEpdRaw1st = Qx_raw_1st_epd; 
	fMyEvent->QyEpdRaw1st = Qy_raw_1st_epd; 
	                        
	fMyEvent->QxEpdPhi1st = Qx_phi_1st_epd; 
	fMyEvent->QyEpdPhi1st = Qy_phi_1st_epd; 

	fMyEvent->QxEpdRaw2nd = Qx_raw_2nd_epd; 
	fMyEvent->QyEpdRaw2nd = Qy_raw_2nd_epd; 
	                        
	fMyEvent->QxEpdPhi2nd = Qx_phi_2nd_epd; 
	fMyEvent->QyEpdPhi2nd = Qy_phi_2nd_epd; 


	//fMyEvent->fMyEpdHit.resize(nHitsEPD);

  // find max global track index this loop was also used to find Psi2, you can revive that
  int maxGBTrackIndex = -1;
  for (unsigned int iTrack = 0; iTrack < mPicoDst->numberOfTracks(); iTrack++)
  {
    StPicoTrack *gTrack = mPicoDst->track(iTrack);
    if (!gTrack)
      continue;
    int index = gTrack->id();
    if (index > maxGBTrackIndex)
      maxGBTrackIndex = index;
  } // end global track loop

  vector<KFMCTrack> mcTracks(0);
  vector<int> mcIndices(maxGBTrackIndex + 1);
  for (unsigned int iIndex = 0; iIndex < mcIndices.size(); iIndex++)
    mcIndices[iIndex] = -1;

  vector<int> triggeredTracks;

  // Process the event
  if (maxGBTrackIndex > 0)
    mStKFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex + 1);

  bool isGoodKF = mStKFParticleInterface->ProcessEvent(mPicoDst, triggeredTracks); //== here make cuts for nHitsFit et al.,2023.06.21
  if (!isGoodKF)
    return kStOK;

	hVz_after_processEvt_cut ->Fill(Vz);
	

  int centralityBin = -1;
  float centralityWeight = 0.;

 // int eventId = -1;
 // int runId = -1;

  centralityWeight = 1;
  Int_t nevent = 100000;

  mStKFParticlePerformanceInterface->SetMCTracks(mcTracks);
  mStKFParticlePerformanceInterface->SetMCIndexes(mcIndices);
  mStKFParticlePerformanceInterface->SetCentralityBin(centralityBin);
  mStKFParticlePerformanceInterface->SetCentralityWeight(centralityWeight);
  mStKFParticlePerformanceInterface->SetPrintEffFrequency(nevent);
  mStKFParticlePerformanceInterface->PerformanceAnalysis();

  // make a map from trackID to track index in global track array
  // Do not do this step unless you want some information which is not in the KFParticle
  vector<int> trackMap;
  trackMap.resize(maxGBTrackIndex + 1, -1);

  for (unsigned int iTrack = 0; iTrack < mPicoDst->numberOfTracks(); iTrack++)
  {
    StPicoTrack *gTrack = mPicoDst->track(iTrack);
    if (!gTrack)
      continue;
    int index = gTrack->id();
    trackMap[index] = iTrack;
  } // end global track loop

  float ld_chi2topo;
  float ld_chi2ndf;
  float ld_ldl;
  float ld_l;
  float ld_dca;
	float ld_px,ld_py,ld_pz; 
	int nLambda = 0;

  float chi2primary_proton;
  float Proton_m2;
  float DcaProton;
  float nSigmaP_proton;
  int   proton_nhits;
  int   proton_nhitsmax;
  float proton_px, proton_py, proton_pz;

  float chi2primary_pion;
  float Pion_m2;
  float DcaPion;
  float nSigmaPi_pi;
  int   pion_nhits;
  int   pion_nhitsmax;
  float pi_px, pi_py, pi_pz;
  int ProtonTrackIndex, PionTrackIndex;
	float ld_mass;

  // cout<<"start loop KFP"<<endl;
		//	if(mEventsProcessed == 7542 )    cout << "the number of rc particles is " << mStKFParticlePerformanceInterface->GetNReconstructedParticles() << endl;

				// Loop over the KFParticles to save/analyze the Lambdas
  for (int iParticle = 0; iParticle < mStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
  {
    KFParticle particle = mStKFParticleInterface->GetParticles()[iParticle];
		//if(mEventsProcessed == 7542 ) cout<<" particle PDG= "<<particle.GetPDG()<<endl;
    // cout<<"particle PDG= "<<particle.GetPDG()<<endl;
    if (LambdaPdg == particle.GetPDG())
    {
      
			ld_px= particle.GetPx();	
			ld_py= particle.GetPy();	
			ld_pz= particle.GetPz();	
			ld_mass = particle.GetMass();
			
			KFParticleSIMD tempSIMDParticle(particle);
      float_v l,dl;
      KFParticleSIMD pv(mStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
      tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
      tempSIMDParticle.SetProductionVertex(pv);
      ld_dca = tempSIMDParticle.GetDistanceFromVertex(pv)[0];
      ld_ldl = l[0]/dl[0];
      ld_l = l[0];
      ld_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
      ld_chi2ndf = particle.Chi2() / particle.NDF();
			// loop over Lambda daughters
      for (int iDaughter = 0; iDaughter < particle.NDaughters(); iDaughter++)
      {

        const int daughterId = particle.DaughterIds()[iDaughter];
        KFParticle daughter = mStKFParticleInterface->GetParticles()[daughterId];
       const int globalTrackId = daughter.DaughterIds()[0];
        int gTrackIndex = trackMap[globalTrackId];

        // get proton
        if (ProtonPdg == daughter.GetPDG())
        {
					float proton_mass_check=-999.; 
					float proton_massError_check=-999.;
					daughter.GetMass(proton_mass_check,proton_massError_check);
					//cout<<"proton mass="<< daughter.GetMass()<<endl;
					//cout<<"proton mass="<< proton_mass_check << " +- "<< proton_massError_check<<endl;
          
					ProtonTrackIndex = gTrackIndex;
          Proton_m2 = mStKFParticleInterface->GetMass2(daughter.DaughterIds()[0]);
          nSigmaP_proton  = mStKFParticleInterface->GetdEdXNSigmaProton(daughter.DaughterIds()[0]);
          chi2primary_proton = daughter.GetDeviationFromVertex(mStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex()); //== added by q.hu on 2023.06.15
          DcaProton = mStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
          proton_nhits = mStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
          proton_nhitsmax = mStKFParticleInterface->GetNHitsMax(daughter.DaughterIds()[0]);
					
					daughter.SetProductionVertex(particle);
					daughter.TransportToProductionVertex();

					proton_px = daughter.GetPx();
					proton_py = daughter.GetPy();
					proton_pz = daughter.GetPz();

        } // end proton
        // get pion
        if (PiMinusPdg == daughter.GetPDG())
        {
         
					PionTrackIndex = gTrackIndex;
          Pion_m2 = mStKFParticleInterface->GetMass2(daughter.DaughterIds()[0]);
          nSigmaPi_pi  = mStKFParticleInterface->GetdEdXNSigmaPion(daughter.DaughterIds()[0]);
          chi2primary_pion = daughter.GetDeviationFromVertex(mStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex()); //== added by q.hu on 2023.06.15
          DcaPion = mStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
          pion_nhits = mStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
          pion_nhitsmax = mStKFParticleInterface->GetNHitsMax(daughter.DaughterIds()[0]);

					daughter.SetProductionVertex(particle);
					daughter.TransportToProductionVertex();

					pi_px = daughter.GetPx();
					pi_py = daughter.GetPy();
					pi_pz = daughter.GetPz();

				} // end pi minus
      }   // end Lambda daughter loop 
      //if(!PassLambdaSelection(ld_chi2topo,ld_chi2ndf,ld_ldl,chi2primary_proton,chi2primary_pion,ld_mass)) continue; // TODO
			if(ld_mass<1.08 || ld_mass>1.18) continue;
//			if(mEventsProcessed == 1221 ) cout<<"decayed proton px = "<<proton_px<<endl;
//			if(mEventsProcessed == 1221 ) cout<<"Event ID= "<<picoEvent->eventId()<<endl;
      fMyEvent->fMyLambda[nLambda].mLdChiTopo = ld_chi2topo; //TODO
      fMyEvent->fMyLambda[nLambda].mLdChi2ndf = ld_chi2ndf;
      fMyEvent->fMyLambda[nLambda].mLdLdl = ld_ldl;
      fMyEvent->fMyLambda[nLambda].mLdL = ld_l;
      fMyEvent->fMyLambda[nLambda].mLdMass = ld_mass;
      fMyEvent->fMyLambda[nLambda].mLdDca = ld_dca;

      fMyEvent->fMyLambda[nLambda].mLdPx = ld_px;
      fMyEvent->fMyLambda[nLambda].mLdPy = ld_py;
      fMyEvent->fMyLambda[nLambda].mLdPz = ld_pz;

      fMyEvent->fMyLambda[nLambda].mProtonPx = proton_px;
      fMyEvent->fMyLambda[nLambda].mProtonPy = proton_py;
      fMyEvent->fMyLambda[nLambda].mProtonPz = proton_pz;
      fMyEvent->fMyLambda[nLambda].mProtonId = ProtonTrackIndex;
      fMyEvent->fMyLambda[nLambda].mProtonMass2= Proton_m2;
      fMyEvent->fMyLambda[nLambda].mProtonNsigma = nSigmaP_proton;
      fMyEvent->fMyLambda[nLambda].mProtonNhits = proton_nhits;
      fMyEvent->fMyLambda[nLambda].mProtonNhitsMax = proton_nhitsmax;
      fMyEvent->fMyLambda[nLambda].mProtonChi2Primary = chi2primary_proton;
      fMyEvent->fMyLambda[nLambda].mProtonDCA = DcaProton;

      fMyEvent->fMyLambda[nLambda].mPionPx = pi_px;
      fMyEvent->fMyLambda[nLambda].mPionPy = pi_py;
      fMyEvent->fMyLambda[nLambda].mPionPz = pi_pz;
      fMyEvent->fMyLambda[nLambda].mPionId = PionTrackIndex;
      fMyEvent->fMyLambda[nLambda].mPionMass2= Pion_m2;
      fMyEvent->fMyLambda[nLambda].mPionNsigma = nSigmaPi_pi;
      fMyEvent->fMyLambda[nLambda].mPionNhits = pion_nhits;
      fMyEvent->fMyLambda[nLambda].mPionNhitsMax = pion_nhitsmax;
      fMyEvent->fMyLambda[nLambda].mPionChi2Primary = chi2primary_pion;
      fMyEvent->fMyLambda[nLambda].mPionDCA = DcaPion;

      nLambda++;
    }
  } // end loop over KFParticles
	//if(nLambda!=0) cout<<	mEventsProcessed<<endl;
	//if(mEventsProcessed == 1221 ) cout<<"lambda number ="<<nLambda<<endl;
  fMyEvent->nLambda = nLambda;
  fMyEvent->fMyLambda.resize(nLambda);

	fEventTree->Fill();
  mEventsProcessed++;
	
//this part is used to check the 
	//if(mEventsProcessed == 1   ) cout<<setw(10) <<"Events="<<setw(10) << mEventsProcessed<<setw(10) <<" countrefmult = "<< setw(10) <<countrefmult<<setw(10) <<" runId="<<setw(10) <<picoEvent->runId()<<setw(10) <<"nLambda="<< nLambda  <<endl;
//	if(mEventsProcessed == 567 ) cout<<setw(10) <<"Events="<<setw(10) << mEventsProcessed<<setw(10) <<" countrefmult = "<< setw(10) <<countrefmult<<setw(10) <<" runId="<<setw(10) <<picoEvent->runId()<<setw(10) <<"nLambda="<< nLambda <<endl;
//	if(mEventsProcessed == 7543 ) cout<<setw(10) <<"Events="<<setw(10) << mEventsProcessed<<setw(10) <<" countrefmult = "<< setw(10) <<countrefmult<<setw(10) <<" runId="<<setw(10) <<picoEvent->runId()<<setw(10) <<"nLambda="<< nLambda <<endl;
//	if(mEventsProcessed == 10082) cout<<setw(10) <<"Events="<<setw(10) << mEventsProcessed<<setw(10) <<" countrefmult = "<< setw(10) <<countrefmult<<setw(10) <<" runId="<<setw(10) <<picoEvent->runId()<<setw(10) <<"nLambda="<< nLambda <<endl;
                                                                                                                                   
//	if(mEventsProcessed == 1   ) cout<<setw(10) <<"Events="<<setw(10) << mEventsProcessed<<setw(10) <<" "<< setw(10) << "nTpcRefA= "<<setw(10)<<nTrackRefA<<setw(10) <<"nEpdHit="<<setw(10)<<nHitsEPD <<endl;
//	if(mEventsProcessed == 567 ) cout<<setw(10) <<"Events="<<setw(10) << mEventsProcessed<<setw(10) <<" "<< setw(10) << "nTpcRefA= "<<setw(10)<<nTrackRefA<<setw(10) <<"nEpdHit="<<setw(10)<<nHitsEPD <<endl;
//	if(mEventsProcessed == 7543 ) cout<<setw(10) <<"Events="<<setw(10) << mEventsProcessed<<setw(10) <<" "<< setw(10) << "nTpcRefA= "<<setw(10)<<nTrackRefA<<setw(10) <<"nEpdHit="<<setw(10)<<nHitsEPD <<endl;
//	if(mEventsProcessed == 10082 ) cout<<setw(10) <<"Events="<<setw(10) << mEventsProcessed<<setw(10) <<" "<< setw(10) << "nTpcRefA= "<<setw(10)<<nTrackRefA<<setw(10) <<"nEpdHit="<<setw(10)<<nHitsEPD <<endl;


  return kStOK;
}

Int_t StKFParticleAnalysisMaker::Finish()
{

  cout << "Finishing" << endl;
  // epd
  //  mEpFinder->Finish();

  cout << endl
       << endl
       << "Events Started = " << mEventsStarted << "     Events Fully Processed = " << mEventsProcessed << endl
       << endl;

	fOutputFile->cd();	//fOutputFile->Write();	
	fEventTree->Write(); 

	//h3f_tpc_RunPhiEta_pos->Write(); 
	//h3f_tpc_RunPhiEta_neg->Write(); 
	//h3f_tpc_PtEtaPhi_pos ->Write(); 
	//h3f_tpc_PtEtaPhi_neg ->Write(); 

	hVz_before_cut					->Write();  
	hVz_after_trigger_cut   ->Write();
	hVz_after_vertex_cut    ->Write();
	hVz_after_badrun_cut    ->Write();
	hVz_after_centrality_cut->Write();
	hVz_after_processEvt_cut->Write();

	hCent->Write();

  delete mEpdGeom;
  delete mEpFinder;

  return kStOK;
}

//------------------------------------------------------------
bool StKFParticleAnalysisMaker::AcceptTrigger(StPicoEvent *picoEvent)
{
  if (!(picoEvent->isTrigger(620052)))
    return false;
  else
    return true;
}
//----------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::AcceptEvent(StPicoEvent *picoEvent)
{
  // Cut parameters for each event
  Float_t VertexZMin = 198.0;       // distance event occurs along beam direction from center (cm)
  const Float_t VertexZMax = 202.0; // cm

  Float_t vertex[3] = {0};
  vertex[0] = picoEvent->primaryVertex().X();
  vertex[1] = picoEvent->primaryVertex().Y();
  vertex[2] = picoEvent->primaryVertex().Z();
  // histogram[0]->Fill(6);
  //  Cut on Vertex location
  if (vertex[2] < VertexZMin || vertex[2] > VertexZMax)
    return false;

  if (sqrt(vertex[0] * vertex[0] + (vertex[1] + 2.0) * (vertex[1] + 2.0)) > 1.5)
    return false;

  return true;
}
//----------------------------------------------------

bool StKFParticleAnalysisMaker::AcceptTrack(StPicoTrack *track)

{

  // Cut Parameters for individual tracks

  const Float_t dcaCut = 3.;  // cm
  const Float_t PtMin = 0.15; // GeV
  const Float_t PtMax = 10.0; // GeV
  const Float_t EtaMin = -1.5;
  const Float_t EtaMax = 1.5;
  const Float_t FitRatio = 0.52; // Number of hits over number of hits possible
  const Int_t nHitMin = 15;      // 15 is typical but sometimes goes as high as 25
  const Int_t nHitMax = 100;     // 45 pad rows in the TPC and so anything bigger than 45+Silicon is infinite
  const Int_t nHitPossMin = 5;   // Don't bother to fit tracks if # possible hits is too low, also protects / 0.

  TVector3 trackMomentum = track->gMom();

  // histogram[6]->Fill(0);

  // StPicoTrack doesn't have "flag" so I just got rid of this cut. This is an ancient cut
  // if ( track->flag() <  0  )     return false ;          // Track quality
  // histogram[6]->Fill(1);
  // if ( track->dcaGlobal().mag() < dcaCut  )  return false ;   // 3D DCA for global tracks

  if (track->nHitsMax() < nHitPossMin)
    return false; // Minimum number of Possible hits, see above.
  // histogram[6]->Fill(2);

  // if ( trackMomentum.Eta()       <  EtaMin  || trackMomentum.Eta()      >  EtaMax  ) return false ;
  // histogram[6]->Fill(3);

  if (trackMomentum.Pt() < PtMin || trackMomentum.Pt() > PtMax)
    return false;
  // histogram[6]->Fill(4);

  if (track->nHitsFit() < nHitMin || track->nHitsFit() > nHitMax)
    return false;
  // histogram[6]->Fill(5);

  if (((float)track->nHitsFit() / (float)track->nHitsMax()) < FitRatio)
    return false;
  // histogram[6]->Fill(6);

  return true;
}


//--------------------------------------------------------------
TVector3 StKFParticleAnalysisMaker::ConvertToTV3(StThreeVectorF inVector)
{
  // Just convert StThreeVectorF to TVector3
  TVector3 outVector;
  outVector.SetXYZ(inVector.x(), inVector.y(), inVector.z());
  return outVector;
}
//--------------------------------------------------------------
double StKFParticleAnalysisMaker::CalculateHelicity(TLorentzVector PLam4, TLorentzVector PPro4)
{
  // function to calculate helicity, p_{proton}^{*} . p_{Lambda}

  PPro4.Boost(-PLam4.BoostVector());
  TVector3 ProtonInLambdaFrame = PPro4.Vect();
  ProtonInLambdaFrame = ProtonInLambdaFrame * (1. / (ProtonInLambdaFrame.Mag()));
  TVector3 LambdaInLabFrame = PLam4.Vect();
  LambdaInLabFrame = LambdaInLabFrame * (1. / (LambdaInLabFrame.Mag()));

  return ProtonInLambdaFrame.Dot(LambdaInLabFrame);
}
//--------------------------------------------------------------
double StKFParticleAnalysisMaker::CalculatePolarization(bool isLambda, double Psi1, TLorentzVector PLam4, TLorentzVector PPro4)
{
  // function to calculate polarization, pass isLambda = true for a Lambda and isLambda = false for an antiLambda
  // I am not adding the requisite 8/(pi*alpha) prescale or resolution correction.
  // Vectors are all assumed to be given in the lab frame

  PPro4.Boost(-PLam4.BoostVector());
  TVector3 SpinVec = PPro4.Vect();
  if (!isLambda)
    SpinVec = -SpinVec;
  float SpinPhi = SpinVec.Phi();
  if (SpinPhi < 0.)
    SpinPhi += 2. * TMath::Pi();
  double Polarization = sin(Psi1 - SpinPhi);
  return Polarization;
}

//--------------------------------------------------------------

bool StKFParticleAnalysisMaker::PassTpcSelection(Float_t eta_min, Float_t eta_max, StPicoTrack *track, TVector3 pV)
{
  if (!track) return false;
  float eta = track->pMom().PseudoRapidity();
  if(eta<eta_min || eta> eta_max) return false;

  float pt = track->pMom().Perp(); // zero for global tracks
  //if (pt < 0.15 || pt > 2.0) return false;
  if (pt < 0.15 || pt > 1.2) return false; //based on Joe's analysis note

  int q = track->charge();
  if (q != 1 && q != -1) return false;
	
  float ddca = track->gDCA(pV).Mag();
  if (fabs(ddca) >= 3.0) return false;

  int nHitsFit = track->nHitsFit();
  float ratio = (float)nHitsFit / (float)track->nHitsMax();
  if (ratio < 0.52 || ratio >= 1.05) return false;
  if (nHitsFit <= 15) return false;
  return true;
}

bool StKFParticleAnalysisMaker::PassEpdSelection(Int_t row, Short_t side, Float_t nMip){
    if(row < 2  || row >16)   return false;
    if(side > 0)    return false;
    if(nMip< 0.3)   return false;
    return true;
}

bool StKFParticleAnalysisMaker::PassLambdaSelection(double ld_chi2topo, double ld_chi2ndf, double ld_ldl, double chi2primary_p, double chi2primary_pi, double ld_mass){
  if(ld_mass<1.08 || ld_mass>1.2) return false;
  if(ld_chi2topo>5) return false;
  if(ld_chi2ndf>5)  return false;
  if(ld_ldl<2)        return false;
  if(chi2primary_p<5)  return false;
  if(chi2primary_pi<5) return false;
  return true;
}






Int_t StKFParticleAnalysisMaker::MakeRunPointer(int runnumber)
{
  int pointer = -999;
  for (int i = 0; i < totalRunNumber; i++)
  {

    if (runnumber == numbers_3[i])
      pointer = i; //==q.hu, 2022.03.25
  }

  if (pointer == -999)
    cout << "Run number are not found! " << runnumber << endl;
  return pointer;
}

Int_t StKFParticleAnalysisMaker::Centrality(int gRefMult)
{
  int centrality;
  int centFull[9] = {5, 9, 16, 26, 41, 60, 86, 119, 142};
  if (gRefMult >= centFull[8])
    centrality = 8;
  else if (gRefMult >= centFull[7])
    centrality = 7;
  else if (gRefMult >= centFull[6])
    centrality = 6;
  else if (gRefMult >= centFull[5])
    centrality = 5;
  else if (gRefMult >= centFull[4])
    centrality = 4;
  else if (gRefMult >= centFull[3])
    centrality = 3;
  else if (gRefMult >= centFull[2])
    centrality = 2;
  else if (gRefMult >= centFull[1])
    centrality = 1;
  else if (gRefMult >= centFull[0])
    centrality = 0;
  else
    centrality = -1;

  return centrality;
}

Double_t StKFParticleAnalysisMaker::FitWeight(Double_t refMult)
{
  Double_t Weight = 1.0;

  Double_t par0 = 1.31492e+00;
  Double_t par1 = -1.66622e+01;
  Double_t par2 = 2.09257e+00;
  Double_t par3 = -4.35674e+00;
  Double_t par4 = -2.38170e-03;
  Double_t par5 = 8.48711e+02;
  Double_t par6 = 4.95939e-06;

  if (
      refMult != -(par3 / par2) // avoid denominator = 0
      && refMult < 70)
  {
    Weight = par0 + par1 / (par2 * refMult + par3) + par4 * (par2 * refMult + par3) + par5 / pow(par2 * refMult + par3, 2) + par6 * pow(par2 * refMult + par3, 2); // Parametrization of MC/data RefMult ratio
  }

  return Weight;
}
