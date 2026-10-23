#include "CeeV0LiteFinder.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
using namespace CeeV0Lite;

TrackInput Input(const TrackState& s) {
  TrackInput t; t.id=s.id;t.positionCm=s.positionCm;
  t.trackingMomentumGeV=s.momentumGeV;Normalize(s.momentumGeV,t.direction);
  t.trackingCovariance=s.covariance;t.positionCovarianceCm2=PositionBlock(s.covariance);
  t.hasTrackingMomentum=t.hasTrackingCovariance=t.chargeKnown=true;t.charge=s.charge;
  t.hypothesisMask=s.charge>0?kProtonHypothesis:kPionHypothesis;return t;
}
int main() {
  FinderConfig cfg; cfg.magneticFieldByTesla=.5;cfg.momentumMode=kTrackingMomentum;
  EventInput ev;ev.primaryVertex.positionCm=Vec3(.04,-.02,.06);
  ev.primaryVertex.covarianceCm2=Mat3::Identity(.0025);
  TrackState s;s.id=17;s.positionCm=Vec3(3.,2.,30.);s.momentumGeV=Vec3(.2,.09,.8);s.charge=1;
  // A full positive covariance with nonzero position-momentum correlations.
  double a[6][6]={{0.}};
  for(int i=0;i<6;++i){a[i][i]=i<3?.08:.004;if(i>0)a[i][i-1]=.001;}
  for(int i=0;i<6;++i)for(int j=0;j<6;++j)for(int k=0;k<6;++k)
    s.covariance.m[i][j]+=a[i][k]*a[j][k];
  double maxRoundtrip=0.,maxCov=0.;
  for(int charge=-2;charge<=2;++charge)if(charge)for(int sign=-1;sign<=1;sign+=2){
    s.charge=charge;cfg.magneticFieldByTesla=sign*.5;
    TrackInput raw=Input(s),out;TrackReferenceDiagnostic info;
    assert(PrepareTrackReference(raw,ev,cfg,out,&info));assert(info.status==kReferenceAtPvDca);
    assert(info.dcaCm>.01); // DCA selection must not move daughter onto the PV.
    assert(std::fabs(Dot(out.positionCm-ev.primaryVertex.positionCm,out.direction))<1.e-8);
    // Local geometric minimum, checked on either side independently.
    TrackState ref=s;ref.positionCm=out.positionCm;ref.momentumGeV=out.trackingMomentumGeV;ref.covariance=out.trackingCovariance;
    for(int side=-1;side<=1;side+=2){TrackState near;assert(PropagateState(ref,side*.01,cfg.magneticFieldByTesla,cfg.transportConfig,near,false));assert(Norm(near.positionCm-ev.primaryVertex.positionCm)>info.dcaCm-1.e-10);}
    TrackState back;assert(PropagateState(ref,-info.closestApproach.signedPathCm,cfg.magneticFieldByTesla,cfg.transportConfig,back));
    maxRoundtrip=std::max(maxRoundtrip,Norm(back.positionCm-s.positionCm)+Norm(back.momentumGeV-s.momentumGeV));
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)maxCov=std::max(maxCov,std::fabs(back.covariance.m[i][j]-s.covariance.m[i][j])/std::sqrt(s.covariance.m[i][i]*s.covariance.m[j][j]));
    TrackInput twice;TrackReferenceDiagnostic twiceInfo;assert(PrepareTrackReference(out,ev,cfg,twice,&twiceInfo));assert(Norm(twice.positionCm-out.positionCm)<1.e-7);
    EventInput other=ev;other.primaryVertex.covarianceCm2=Mat3::Identity(100.);
    TrackInput noConstraint;assert(PrepareTrackReference(raw,other,cfg,noConstraint));
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)assert(noConstraint.trackingCovariance.m[i][j]==out.trackingCovariance.m[i][j]);
    FinderConfig noisy=cfg;noisy.referenceTransportConfig.positionProcessNoiseCm2PerCm=1.;noisy.referenceTransportConfig.angularProcessNoiseRad2PerCm=1.;
    TrackInput pure;assert(PrepareTrackReference(raw,ev,noisy,pure));
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)assert(pure.trackingCovariance.m[i][j]==out.trackingCovariance.m[i][j]);
  }
  assert(maxRoundtrip<1.e-9&&maxCov<1.e-8);
  cfg.magneticFieldByTesla=.5;s.charge=1;TrackInput raw=Input(s),out;out.id=987;TrackReferenceDiagnostic info;
  FinderConfig limited=cfg;limited.referenceTransportConfig.maximumAbsPathCm=.001;
  assert(!PrepareTrackReference(raw,ev,limited,out,&info));assert(out.id==987);assert(info.status==kReferenceTransportFailed);assert(info.closestApproach.hitPathLimit);
  TrackInput missing=raw;missing.hasTrackingCovariance=false;
  assert(!PrepareTrackReference(missing,ev,cfg,out,&info));assert(info.status==kReferenceInvalidInput);
  TrackInput nan=raw;nan.trackingCovariance.m[0][0]=std::numeric_limits<double>::quiet_NaN();assert(!PrepareTrackReference(nan,ev,cfg,out));
  FinderConfig zero=cfg;zero.magneticFieldByTesla=0.;assert(PrepareTrackReference(missing,ev,zero,out,&info));assert(info.status==kReferenceUnchanged);assert(out.positionCm[2]==missing.positionCm[2]);
  // End-to-end: automatic once-per-track preparation equals explicit preparation
  // followed by stored-reference reconstruction, with same IDs and output C.
  std::vector<TrackInput> inputs,anchored;
  for(int d=0;d<2;++d){TrackState at=s;at.id=d;at.charge=d?-1:1;at.positionCm=Vec3(0.,0.,5.);at.momentumGeV=d?Vec3(-.1,0.,.13):Vec3(.1,0.,.71);TrackState far;assert(PropagateState(at,20.,.5,cfg.transportConfig,far));inputs.push_back(Input(far));TrackInput t;assert(PrepareTrackReference(inputs.back(),ev,cfg,t));anchored.push_back(t);}
  cfg.useLambdaDcaRoleAssignment=false;cfg.lambdaCuts.maximumGeometricChi2Ndf=1.e8;cfg.lambdaCuts.maximumTopologicalChi2Ndf=1.e8;cfg.lambdaCuts.minimumPointingCosine=-1.;cfg.lambdaCuts.minimumDecayLengthCm=0.;
  std::vector<Candidate> automatic,manual;std::vector<TrackReferenceDiagnostic> diag;
  Finder(cfg).Reconstruct(inputs,ev,automatic,&diag);assert(diag.size()==2);assert(automatic.size()==1);
  FinderConfig stored=cfg;stored.usePrimaryVertexDcaReference=false;
  Finder(stored).Reconstruct(anchored,ev,manual);assert(manual.size()==automatic.size());
  assert(automatic[0].massGeV==manual[0].massGeV);
  for(int i=0;i<7;++i)for(int j=0;j<7;++j)assert(automatic[0].mother.covariance.m[i][j]==manual[0].mother.covariance.m[i][j]);
  cfg.magneticFieldByTesla=stored.magneticFieldByTesla=0.;
  Finder(cfg).Reconstruct(inputs,ev,automatic);Finder(stored).Reconstruct(inputs,ev,manual);
  assert(automatic.size()==manual.size());assert(automatic.size()==1);assert(automatic[0].massGeV==manual[0].massGeV);
  inputs.push_back(missing);cfg.magneticFieldByTesla=.5;
  Finder(cfg).Reconstruct(inputs,ev,automatic,&diag);assert(diag.size()==3);assert(diag[2].status==kReferenceInvalidInput);assert(automatic.size()==1);
  std::cout<<"PV DCA reference: polarity/charge, stationarity, reversibility, full covariance, PV/noise independence, failures, Finder integration, zero-field bypass passed; max roundtrip="<<maxRoundtrip<<" scaled C="<<maxCov<<"\n";
}
