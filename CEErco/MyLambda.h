//////////////////////////////////////////////////////////////////////////////////
// MyLambda class: lambda candidates
//
// v1.0 Xionghong 2024.12.2
//////////////////////////////////////////////////////////////////////////////////

#ifndef MYLAMBDA_H
#define MYLAMBDA_H

#include "TVector3.h"

class MyTrack;

class MyLambda
{
public:
	
	MyLambda(MyTrack* proton, MyTrack* pim, Bool_t bkg=false);
	~MyLambda();

	//// decay vertex
	TVector3 decayVertex() const { return fdecayVertex; }
	//// decay length
	Double_t decayLength() const { return fdecayLength; }
	//// decay time
	Double_t decayTime() const { return fdecayTime; }
	//// momentum
	TVector3 mom() const { return fmom; }
	//// dca
	Double_t dca() const { return  fdca; }
	//// dca of proton with pion
	Double_t dca1to2() const { return  fdca1to2; }

	//// modifier
	void setDecayVertex(const TVector3& val)    { fdecayVertex = val; }
	void setDecayLength(Double_t val)    { fdecayLength = val; }
	void setMom(const TVector3& val)    { fmom = val; }
	void setDca(Double_t val)      { fdca = val; }
	void setDecayTime(Double_t val) { fdecayTime = val; }

private:

	MyTrack* fproton;
	MyTrack* fpim;
	
	TVector3 fdecayVertex;
	TVector3 fmom;
	Double_t fdca;
	Double_t fdca1to2;
	Double_t fdecayLength;
	Double_t fdecayTime;
	Bool_t fRotPi;

	//// dca with another track
	void calculateDca1to2();
	//// cross point of two tracks
	void recVertex();
	//// reset momemtum of two daugters
	void calculateMom();

	double lorentzFactor(double v);
	//void calculate_v1_v2(double v0, double l0, double l1, double l2, double t1, double t2, double &v1, double &v2);
	//double F(double v0, double m0, double m1, double m2, double l0, double l1, double l2, double t1, double t2);
	//double binary_search(double m0, double m1, double m2, double l0, double l1, double l2, double t1, double t2, double a, double b, double tol);

};

#endif
