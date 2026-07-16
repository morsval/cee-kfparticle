//////////////////////////////////////////////////////////////////////////////////
// MyTrack class: track reconstruction from ToF hit and start vertex
//
// v1.0 Xionghong 2024.12.2
//////////////////////////////////////////////////////////////////////////////////

#ifndef MYTRACK_H
//#define MUTRACK_H
#define MYTRACK_H
#include "TVector3.h"
#include "TLorentzVector.h"
//add by jcc
#include "Rtypes.h"

class MyTrack
{
public:
	
	MyTrack(const TVector3& startPos, const TVector3& tofHit, const TVector3& evtVtx, Double_t tofTime, Float_t mass, TVector3 direction = TVector3(0., 0., 0.));
	MyTrack(const MyTrack& trk);
	~MyTrack();

	////
	TVector3 evtVtx() const {return fevtVtx; }

	//// one point in the track
	TVector3 pos() const { return fstartPos; }
	//// tof hit of the track
	TVector3 tofHit() const { return ftofHit; }
	//// tof flight time
	Double_t tofTime() const { return ftofTime; }
	//// momentum of the track
	Double_t px() const {return fPx;}
	Double_t py() const {return fPy;}
	Double_t pz() const {return fPz;}
	Double_t e() const {return fE;}
	TVector3 mom() const { return fmom; }
	//// 4momentum of the track
	TLorentzVector p4() const { return fp4; }
	//// velocity of the track
	Double_t velocity() const { return fvelocity; }
	//// direction of the track
	TVector3 direction() const { return fdirection; }
	//// mass
	Double_t mass() const { return fmass; }
	//// dca
	Double_t dca() const { return  fdca; }
	//// pdgcode
	Int_t pdgCode() const { return fpdg; }
	//// momPdgCode
	Int_t momPdgCode() const { return fmomPdg; }
      
  
	//// modifier
	void setPos(const TVector3& val)    { fstartPos = val; }
	void setMom(const TVector3& val);
	void setMass(Double_t val)     { fmass = val; }
	void setDca(Double_t val)      { fdca = val; }
	void setPdgCode(Int_t val)    { fpdg = val; }
	void setMomPdgCode(Int_t val) { fmomPdg = val; }
	void setTofTime(Double_t val)  { ftofTime = val; }
	void setTofLength(Double_t val)  { ftofLength = val; }
	void setVelocity( Double_t val) {fvelocity = val; }
	void setDirection(const TVector3& val) { fdirection = val; }

	void RotateXY(Double_t angle);

private:
	
	TVector3 fevtVtx;
	TVector3 fdirection;
	TVector3 fstartPos;
	TVector3 ftofHit;
	TVector3 fmom;
	TLorentzVector fp4;
	Double_t fdca;
	Double_t fmass;
	Double_t ftofTime;
	Double_t ftofLength;
	Double_t fvelocity;
	Double_t fPx;
	Double_t fPy;
	Double_t fPz;
	Double_t fE;
	Int_t fpdg;
	Int_t fmomPdg;
        
        //add by jcc
        Bool_t fUseExternalMom;
 
	void reconstructTrack();

};
#endif
