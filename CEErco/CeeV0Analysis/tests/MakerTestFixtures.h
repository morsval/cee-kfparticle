#ifndef CEE_MAKER_TEST_FIXTURES_H
#define CEE_MAKER_TEST_FIXTURES_H
#include "CeeKFParticleAnalysisMaker.h"
#include "CeeKFParticleObservables.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
using namespace CeeV0Analysis;
using namespace CeeV0Lite;

AnalysisEvent MakeEvent(bool magnetic) {
  AnalysisEvent e; e.eventId=42; e.inputTrackCount=10;
  e.input.primaryVertex.covarianceCm2=Mat3::Identity(0.0001);
  const double M=1.115683,mp=0.9382720813,mpi=0.13957039,beta=0.6,gamma=1.25;
  const double p=std::sqrt((M*M-(mp+mpi)*(mp+mpi))*(M*M-(mp-mpi)*(mp-mpi)))/(2*M);
  for (int i=0;i<2;++i) {
    AnalysisTrack t; t.nHits=30; t.fitConverged=true; t.chi2Ndf=1.;
    t.nSigmaProton=(i==0 ? 0.5:10.); t.nSigmaPion=(i==1 ? 0.5:10.);
    TrackInput& x=t.input; x.id=100+i;
    const double mass=i==0?mp:mpi;
    const Vec3 momentum(i==0?p:-p,0.,gamma*beta*std::sqrt(mass*mass+p*p));
    const Vec3 vertex(0.,0.,5.); Vec3 u; assert(Normalize(momentum,u));
    x.positionCm=vertex+20.*u; x.direction=u;
    x.positionCovarianceCm2=Mat3::Identity(0.0004); x.sigmaDirectionRad=0.0005;
    x.hypothesisMask=kPionHypothesis|kProtonHypothesis;
    x.chargeKnown=magnetic; x.charge=i==0?1:-1;
    x.hasTrackingMomentum=magnetic; x.trackingMomentumGeV=momentum;
    x.hasTrackingCovariance=magnetic;
    for (int j=0;j<6;++j) x.trackingCovariance.m[j][j]=j<3?0.0004:0.000001;
    x.hasTof=!magnetic; x.tofPositionCm=vertex+100.*u;
    x.tofTimeNs=5./(beta*29.9792458)+100.*std::sqrt(Norm2(momentum)+mass*mass)/(Norm(momentum)*29.9792458);
    x.sigmaTofTimeNs=0.02; x.sigmaPathLengthCm=0.05;
    if (magnetic) {
      TrackState in,out; in.positionCm=vertex; in.momentumGeV=momentum;
      in.charge=x.charge; in.massGeV=mass; in.covariance=x.trackingCovariance;
      assert(PropagateState(in,20.,0.5,TransportConfig(),out));
      x.positionCm=out.positionCm; x.trackingMomentumGeV=out.momentumGeV;
      x.trackingCovariance=out.covariance; assert(Normalize(out.momentumGeV,x.direction));
    }
    e.tracks.push_back(t);
  }
  return e;
}
FinderConfig Configuration(bool magnetic) {
  FinderConfig c; c.reconstructK0Short=false; c.reconstructAntiLambda=false;
  c.momentumMode=magnetic?kTrackingMomentum:kTofMomentum;
  c.magneticFieldByTesla=magnetic?0.5:0.; c.useLambdaDcaRoleAssignment=false;
  c.lambdaCuts.minimumMassGeV=1.; c.lambdaCuts.maximumMassGeV=1.3;
  c.lambdaCuts.maximumDaughterDcaCm=10.; c.lambdaCuts.maximumGeometricChi2Ndf=1.e6;
  c.lambdaCuts.maximumTopologicalChi2Ndf=1.e6;
  c.lambdaCuts.minimumDecayLengthCm=0.; c.lambdaCuts.minimumDecayLengthSignificance=0.;
  c.lambdaCuts.minimumDaughterChi2FromPv=0.; c.lambdaCuts.minimumPointingCosine=-1.;
  return c;
}
unsigned PID(const AnalysisTrack& t,const AnalysisEvent&) {
  unsigned mask=0;
  if (std::isfinite(t.nSigmaProton) && std::fabs(t.nSigmaProton)<3.) mask|=kProtonHypothesis;
  if (std::isfinite(t.nSigmaPion) && std::fabs(t.nSigmaPion)<3.) mask|=kPionHypothesis;
  return mask;
}
#endif
