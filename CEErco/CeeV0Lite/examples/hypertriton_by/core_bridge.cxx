#include "CeeV0LiteVertexFitter.h"
#include "CeeV0LiteFinder.h"
#include "CeeV0LiteProductionVertex.h"
#include <cmath>
#include <limits>
using namespace CeeV0Lite;

// Input per candidate: two Cartesian states (12), two full 6x6 C (72).
// Output: status, chi2, daughter DCA, mother q7, C7, mass, sigma mass.
extern "C" void fit_batch(int n, const double* x, const double* cov,
                           double by, double* out) {
  #pragma omp parallel for schedule(dynamic,32)
  for (int k=0; k<n; ++k) {
    double* o=out+61*k;
    for (int j=0;j<61;++j) o[j]=std::numeric_limits<double>::quiet_NaN();
    TrackState t[2];
    for(int d=0;d<2;++d) {
      t[d].charge=d==0 ? 2 : -1;
      t[d].massGeV=d==0 ? 2.80839160743 : 0.13957039;
      for(int i=0;i<3;++i) {t[d].positionCm[i]=x[k*12+d*6+i];t[d].momentumGeV[i]=x[k*12+d*6+3+i];}
      for(int i=0;i<6;++i) for(int j=0;j<6;++j) t[d].covariance.m[i][j]=cov[k*72+d*36+i*6+j];
    }
    VertexFitter fitter;
    VertexFitResult f;
    if(!fitter.Fit(t[0],t[1],by,f)) {o[0]=-int(f.status);continue;}
    o[1]=f.chi2;o[2]=f.daughterDcaCm;
    for(int i=0;i<3;++i) {o[3+i]=f.mother.positionCm[i];o[6+i]=f.mother.momentumGeV[i];}
    o[9]=f.mother.energyGeV;
    if(!fitter.RecalculateMotherCovariance(t[0],t[1],by,f.mother)) {o[0]=-100-int(f.mother.covarianceStatus);continue;}
    for(int i=0;i<7;++i) for(int j=0;j<7;++j) o[10+7*i+j]=f.mother.covariance.m[i][j];
    if(!MotherMassAndError(f.mother,o[59],o[60])) {o[0]=-200;continue;}
    o[0]=1;
  }
}

extern "C" int propagate(const double* x, int charge, double length, double by, double* out) {
  TrackState a,b;
  a.charge=charge;
  for(int i=0;i<3;++i){a.positionCm[i]=x[i];a.momentumGeV[i]=x[i+3];}
  if(!PropagateState(a,length,by,TransportConfig(),b,false))return 0;
  for(int i=0;i<3;++i){out[i]=b.positionCm[i];out[i+3]=b.momentumGeV[i];}
  return 1;
}

extern "C" int neutral_reference(const double* q, const double* c,
                                  const double* pv, double pvvar,
                                  double* out) {
  ProductionVertexInput in;
  for(int i=0;i<3;++i){in.mother.positionCm[i]=q[i];in.mother.momentumGeV[i]=q[i+3];in.vertex.positionCm[i]=pv[i];in.vertex.covarianceCm2.m[i][i]=pvvar;}
  in.mother.energyGeV=q[6];in.mother.covarianceValid=true;
  for(int i=0;i<7;++i) for(int j=0;j<7;++j)in.mother.covariance.m[i][j]=c[7*i+j];
  ProductionState s;
  if(!FitProductionVertex(in,ProductionVertexConfig(),0.,1,s,0))return 0;
  TransportToDecayVertex(s);
  for(int i=0;i<3;++i){out[i]=s.positionCm[i];out[i+3]=s.momentumGeV[i];}
  out[6]=s.energyGeV;out[7]=s.flightParameterS;
  for(int i=0;i<8;++i)for(int j=0;j<8;++j)out[8+8*i+j]=s.covariance.m[i][j];
  out[72]=s.pvConstraintChi2;return 1;
}

// Formal v0.6.3 charged-mother fit. The example uses independent tracking/PV
// inputs; callers with correlated inputs must fill the public cross block.
// Return at the decay endpoint to compare the same mass with the Python toy.
extern "C" void production_fit_batch(int n, const double* q, const double* cov,
                                     const double* pv, const double* pvc,
                                     int charge, double by, int maximumIterations,
                                     double* out) {
  #pragma omp parallel for schedule(dynamic,32)
  for (int k=0;k<n;++k) {
    double* row=out+75*k;
    for (int i=0;i<75;++i)row[i]=std::numeric_limits<double>::quiet_NaN();
    ProductionVertexInput in;
    for (int i=0;i<3;++i) {
      in.mother.positionCm[i]=q[k*7+i];in.mother.momentumGeV[i]=q[k*7+i+3];
      in.vertex.positionCm[i]=pv[k*3+i];
      for(int j=0;j<3;++j)in.vertex.covarianceCm2.m[i][j]=pvc[k*9+i*3+j];
    }
    in.mother.energyGeV=q[k*7+6];in.mother.covarianceValid=true;
    for(int i=0;i<7;++i)for(int j=0;j<7;++j)in.mother.covariance.m[i][j]=cov[k*49+i*7+j];
    ProductionVertexConfig config;config.maximumIterations=maximumIterations;
    ProductionState state;
    const bool ok=FitProductionVertex(in,config,0.,1,state,charge,by);
    row[73]=state.status;row[74]=state.iterations;
    if(!ok)continue;
    if(!TransportToDecayVertex(state)){row[73]=kProductionUnphysicalState;continue;}
    for(int i=0;i<3;++i){row[i]=state.positionCm[i];row[i+3]=state.momentumGeV[i];}
    row[6]=state.energyGeV;row[7]=state.flightParameterS;
    for(int i=0;i<8;++i)for(int j=0;j<8;++j)row[8+8*i+j]=state.covariance.m[i][j];
    row[72]=state.pvConstraintChi2;
  }
}

// Use exactly the formal Finder reference preparation, including its defaults.
// Per daughter diagnostic: status, signed reference arc, DCA, stationarity.
extern "C" void prepare_reference_batch(int n, const double* x, const double* c,
    const double* pv, double by, double* xo, double* co, double* diagnostics) {
  #pragma omp parallel for schedule(dynamic,32)
  for(int k=0;k<n;++k) {
    FinderConfig cfg;cfg.momentumMode=kTrackingMomentum;cfg.magneticFieldByTesla=by;
    EventInput event;for(int i=0;i<3;++i)event.primaryVertex.positionCm[i]=pv[3*k+i];
    for(int d=0;d<2;++d) {
      TrackInput in,out;TrackReferenceDiagnostic info;
      in.id=2*k+d;in.hasTrackingMomentum=in.hasTrackingCovariance=in.chargeKnown=true;
      in.charge=d==0?2:-1;
      for(int i=0;i<3;++i){in.positionCm[i]=x[12*k+6*d+i];in.trackingMomentumGeV[i]=x[12*k+6*d+3+i];}
      for(int i=0;i<6;++i)for(int j=0;j<6;++j)in.trackingCovariance.m[i][j]=c[72*k+36*d+6*i+j];
      const bool ok=PrepareTrackReference(in,event,cfg,out,&info);
      double* di=diagnostics+8*k+4*d;
      di[0]=info.status;di[1]=info.closestApproach.signedPathCm;di[2]=info.dcaCm;di[3]=info.closestApproach.stationarityResidualCm;
      for(int i=0;i<6;++i)xo[12*k+6*d+i]=ok?(i<3?out.positionCm[i]:out.trackingMomentumGeV[i-3]):std::numeric_limits<double>::quiet_NaN();
      for(int i=0;i<6;++i)for(int j=0;j<6;++j)co[72*k+36*d+6*i+j]=ok?out.trackingCovariance.m[i][j]:std::numeric_limits<double>::quiet_NaN();
    }
  }
}
