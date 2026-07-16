// root class Header
#include <iostream>
#include "TMatrixD.h"
#include "MyTrack.h"
//add 3/12
#include "TRandom3.h"

//static const bool   kApplyTofTimeSmear = true;   // false 就关闭
static const bool   kApplyTofTimeSmear = false;   // false 就关闭
static const double kTofTimeSigma_ns   = 0.0;   // 100 ps sigma; 若要±50ps宽度就改成0.05
static TRandom3 gTofRnd(0);

const Double_t c_speed = 2.99792458*(1e8);

MyTrack::MyTrack(const TVector3& startPos, const TVector3& tofHit, const TVector3& evtVtx, Double_t tofTime, Float_t mass, TVector3 direction)
	: ftofHit(tofHit),
	fdirection(direction),
	fevtVtx(evtVtx),
	ftofTime(tofTime),
	fvelocity(0.),
	fdca(0.),
	fmass(mass),
	fPx(0.),
	fPy(0.),
	fPz(0.),
	fE(0.),
	fp4(TLorentzVector(0., 0., 0., 0.)),
	fpdg(0),
	ftofLength(0.),
	fmomPdg(0)
{
	setPos(startPos);
	
	if(fdirection.Mag()<0.0001)
		fdirection = ftofHit - fstartPos;

	fdirection.SetMag(1.0);
	
	reconstructTrack();
}

MyTrack::~MyTrack() 
{
}

MyTrack::MyTrack(const MyTrack& trk)
	: fevtVtx(trk.fevtVtx),
	fstartPos(trk.fstartPos), 
	ftofHit(trk.ftofHit),
	fdirection(trk.fdirection),
	fmom(trk.fmom),
	fp4(trk.fp4),
	fdca(trk.fdca),             
	fmass(trk.fmass),           
	fpdg(trk.fpdg),
	ftofTime(trk.ftofTime),
	fmomPdg(trk.fmomPdg),
	ftofLength(trk.ftofLength),
	fvelocity(trk.fvelocity)
{
}

void MyTrack::reconstructTrack()
{
	ftofLength = (ftofHit-fevtVtx).Mag();

	if(ftofHit.Mag() <0.0001){
		Info("MyTrack::reconstructTrack", "Tof Hit position is not correct");
		fmom = TVector3(0, 0, 0);
		return;
	}
	//// fast sim
	//Double_t speed = ftofLength*0.01/ftofTime/c_speed; // 0.01: cm->m
	//// full sim
  //itoftime smear
  
  Double_t tofTimeUsed = ftofTime;
      if (kApplyTofTimeSmear) {
        tofTimeUsed = gTofRnd.Gaus(ftofTime, kTofTimeSigma_ns);
        if (tofTimeUsed <= 1e-6) tofTimeUsed = 1e-6;
       }
   ftofTime = tofTimeUsed;     
  
	Double_t speed = ftofLength*0.01/(1.*ftofTime*(1e-9))/c_speed; // 0.01: cm->m, 1e-9: ns->s
	if(speed>=1.0){
		//Warning("reconstructTrack", "track speed > light speed");
		speed = 0.0;
	}
	setVelocity(speed);

	Double_t gamma = 1.0/sqrt(1.0 - (speed*speed));

	TVector3 v = speed*fdirection;
	TVector3 momentum = gamma*fmass*v;

	setMom(momentum);

	fmom = momentum;
	fE = sqrt(momentum.Mag2()+fmass*fmass);
	fPx = momentum.X();
	fPy = momentum.Y();
	fPz = momentum.Z();
	fp4 = TLorentzVector(momentum, fE);
}


void MyTrack::setMom(const TVector3& val)    
{ 
	fmom = val;
	fE = sqrt(val.Mag2()+fmass*fmass);
	fPx = val.X();
	fPy = val.Y();
	fPz = val.Z();
	fp4 = TLorentzVector(val, fE);

}


void MyTrack::RotateXY(Double_t angle)
{
  // float angle - angle of rotation in XY plane in [rad]

  // Before rotation the center of the coordinat system should be moved to the vertex position; move back after rotation
	TVector3 relativePos = fstartPos - fevtVtx;

	TMatrixD rotationMatrix(3, 3);
	rotationMatrix(0, 0) = cos(angle);
	rotationMatrix(0, 1) = -sin(angle);
	rotationMatrix(0, 2) = 0;
	rotationMatrix(1, 0) = sin(angle);
	rotationMatrix(1, 1) = cos(angle);
	rotationMatrix(1, 2) = 0;
	rotationMatrix(2, 0) = 0;
	rotationMatrix(2, 1) = 0;
	rotationMatrix(2, 2) = 1;

	relativePos = rotationMatrix * relativePos;
	fstartPos = fevtVtx + relativePos;

	fmom = rotationMatrix*fmom;
}
