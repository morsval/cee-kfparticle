#include <iostream>
#include <cmath>
#include <functional>
#include "TRandom.h"

// root class Header
#include "MyTrack.h"
#include "MyLambda.h"

const double c_speed = 2.99792458*(1e+8);
const double tol = 1e-6; // 容差
const int max_iter = 20; //最大迭代次数
const double pmass = 0.938272;
const double pimass = 0.13957;
const double ldmass = 1.11568;

MyLambda::MyLambda(MyTrack* proton, MyTrack* pim, Bool_t bkg)
  : fproton(proton),
    fpim(pim),
    fdecayLength(0.),
    fdecayTime(0.),
    fdca(0.),
    fdca1to2(0.),
    fmom(TVector3(0, 0, 0)),
    fdecayVertex(TVector3(0, 0, 0)),
		fRotPi(bkg)
{
	if(fRotPi) {
		gRandom->SetSeed(0);
		Double_t angle = gRandom->Uniform(150,210)/180*TMath::Pi();
		//std::cout << "1: x=" << fpim->pos().X() << ", px=" << fpim->mom().X() << std::endl;
		fpim->RotateXY(angle);
		//std::cout << "2: x=" << fpim->pos().X() << ", px=" << fpim->mom().X() << std::endl;
	}
  
	calculateDca1to2();
  recVertex();
	calculateMom();
}

MyLambda::~MyLambda()
{
}

void MyLambda::calculateDca1to2()
{ 
  // Direction of daughter particles
  TVector3 dir1 = fproton->direction();
  TVector3 dir2 = fpim->direction();
  TVector3 r2_r1 = fproton->pos() - fpim->pos();
  
  // Longitudinal direction (cross product of directions)
  TVector3 crossProd = dir1.Cross(dir2);
  Double_t crossProdMag = crossProd.Mag();

  // Projection of r2_r1 onto the longitudinal direction (crossProd)
	if(crossProdMag==0) // parellal lines
		fdca1to2 = r2_r1.Perp(dir1);
  fdca1to2 = fabs(r2_r1.Dot(crossProd)) / crossProdMag;

	//std::cout << fdca1to2 << std::endl;
}

void MyLambda::recVertex()
{ 
  // Direction of daughter particles
  TVector3 dir1 = fproton->direction();
  TVector3 dir2 = fpim->direction();
  //TVector3 r2_r1 = fproton->pos() - fpim->pos();
  TVector3 r2_r1 = fpim->pos() - fproton->pos();
  
  // Longitudinal direction (cross product of directions)
  TVector3 crossProd = dir1.Cross(dir2);
  Double_t crossProdMag2 = crossProd.Mag2();
  
  if (crossProdMag2 == 0) {
    // If the lines are parallel, set decay vertex to a midpoint of both tracks
    fdecayVertex = 0.5*(fproton->pos()) + 0.5*(fpim->pos());
    return;
  }
  
  // Calculate parameters t and s to find the closest points on both lines
  Double_t t = r2_r1.Cross(dir2).Dot(crossProd) / crossProdMag2;
  Double_t s = r2_r1.Cross(dir1).Dot(crossProd) / crossProdMag2;
  
	//std::cout << "t=" <<t <<", s=" << s << std::endl;
  // Calculate points on both tracks
  TVector3 point1 = fproton->pos() + t * dir1;
  TVector3 point2 = fpim->pos() + s * dir2;

  
  // Set the decay vertex to the midpoint of these two points
  //fdecayVertex = 0.5*point1 + 0.5*point2; // fast simulation
	fdecayVertex = point1; // full simulation
	//std::cout << "rec proton x=" << point1.X() << ", rec pion x=" << point2.X() <<", rec ld x=" << fdecayVertex.X() << std::endl;
	//std::cout << "rec lambda x=" << fdecayVertex.X() << ", y=" << fdecayVertex.Y() << ", z=" << fdecayVertex.Z() << std::endl;

	fdecayLength = (fdecayVertex-fproton->evtVtx()).Mag();
}

void MyLambda::calculateMom(){

   
	////---- method1: proton and pion momentum are calculated assuming primary track
	TVector3 recMom1 = fproton->mom()+fpim->mom();
	
	//std::cout << "rec1=" << (fproton->p4()+fpim->p4()).M() << ", rec2=";
	
	////---- method2: calculate pion momentum
	Double_t protonLength = (fproton->tofHit()-fdecayVertex).Mag()*0.01;
	Double_t protonV = fproton->velocity();
	Double_t protonTime = protonLength/protonV;
	
	fdecayTime = fproton->tofTime()*1e-9*c_speed - protonTime; // 1e-9: ns -> s
	
	// pion
	Double_t pionTime = fpim->tofTime()*1e-9*c_speed - fdecayTime;
	Double_t pionLength = (fpim->tofHit()-fdecayVertex).Mag()*0.01;
	Double_t pionV = pionLength/pionTime;

	Double_t gamma = 1.0/sqrt(1.0 - (pionV*pionV));

	TVector3 v = pionV*fpim->direction();
	TVector3 momentum = gamma*pimass*v;
	fpim->setMom(momentum);
	
	//std::cout << "rec1 mom=" << fpim->mom().Mag() << ", rec2 mom=" << momentum.Mag();// << std::endl;
	TVector3 recMom2 = fproton->mom()+fpim->mom();
	//std::cout << (fproton->p4()+fpim->p4()).M();
	
	setMom(recMom2);	
	//std::cout << "rec1 mom=" << recMom1.Mag() << ", rec2 mom=" << recMom2.Mag();// << std::endl;
  
//  TVector3 recMom2 = fproton->mom()+fpim->mom();
//  setMom(recMom2);
}

double MyLambda::lorentzFactor(double v){
	return 1.0/sqrt(1.0-v*v);
}
