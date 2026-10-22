#include "CeeV0LiteProductionVertex.h"
#include "CeeV0LiteTransport.h"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace CeeV0Lite;
namespace {
void Close(double a,double b,double tolerance=1.e-8) {
  assert(std::isfinite(a) && std::isfinite(b));
  if (std::fabs(a-b)>tolerance*std::max(1.,std::max(std::fabs(a),std::fabs(b)))) {
    std::cerr<<"Mismatch "<<std::setprecision(17)<<a<<" vs "<<b<<"\n";
    assert(false);
  }
}
double Get(const ProductionState& s,int i) {
  return i<3 ? s.positionCm[i] : i<6 ? s.momentumGeV[i-3] : i==6 ? s.energyGeV : s.flightParameterS;
}
void Set(ProductionState& s,int i,double v) {
  if (i<3) s.positionCm[i]=v;
  else if (i<6) s.momentumGeV[i-3]=v;
  else if (i==6) s.energyGeV=v;
  else s.flightParameterS=v;
}
ProductionVertexInput Fixture(int charge,double by,double flight,bool perturb=true) {
  TrackState production,decay;
  production.positionCm=Vec3(.015,-.03,.02);
  production.momentumGeV=Vec3(.31,.17,1.4); production.charge=charge;
  assert(PropagateState(production,flight*Norm(production.momentumGeV),by,
                        TransportConfig(),decay,false));
  ProductionVertexInput in;
  in.mother.positionCm=decay.positionCm;
  in.mother.momentumGeV=decay.momentumGeV;
  in.mother.energyGeV=std::sqrt(Norm2(decay.momentumGeV)+2.99131*2.99131);
  in.mother.covarianceValid=true; in.mother.covarianceStatus=kMotherCovarianceValid;
  in.vertex.positionCm=production.positionCm;
  if (perturb) {
    in.mother.positionCm=in.mother.positionCm+Vec3(.04,-.035,.05);
    in.mother.momentumGeV=in.mother.momentumGeV+Vec3(.001,-.0015,.002);
    in.mother.energyGeV+=.001;
    in.vertex.positionCm=in.vertex.positionCm+Vec3(-.007,.003,.008);
  }
  const double sigma[10]={.10,.12,.15,.012,.015,.023,.015,.035,.04,.045};
  double lower[10][10]={{0.}},c[10][10]={{0.}};
  for (int i=0;i<10;++i) {
    lower[i][i]=sigma[i];
    for (int j=0;j<i;++j) lower[i][j]=sigma[i]*.035*std::sin(double(2*i+j));
  }
  for (int i=0;i<10;++i) for (int j=0;j<10;++j)
    for (int k=0;k<10;++k) c[i][j]+=lower[i][k]*lower[j][k];
  for (int i=0;i<7;++i) for (int j=0;j<7;++j) in.mother.covariance.m[i][j]=c[i][j];
  for (int i=0;i<3;++i) for (int j=0;j<3;++j) in.vertex.covarianceCm2.m[i][j]=c[i+7][j+7];
  for (int i=0;i<7;++i) for (int j=0;j<3;++j) in.motherVertexCovariance[i][j]=c[i][j+7];
  return in;
}
void Print(const ProductionVertexInput& in,const ProductionState& state,int charge,double by) {
  double y[10],c[10][10]={{0.}};
  for (int i=0;i<3;++i) {y[i]=in.mother.positionCm[i];y[i+3]=in.mother.momentumGeV[i];y[i+7]=in.vertex.positionCm[i];}
  y[6]=in.mother.energyGeV;
  for (int i=0;i<7;++i) for(int j=0;j<7;++j)c[i][j]=in.mother.covariance.m[i][j];
  for (int i=0;i<3;++i) for(int j=0;j<3;++j)c[i+7][j+7]=in.vertex.covarianceCm2.m[i][j];
  for (int i=0;i<7;++i) for(int j=0;j<3;++j)c[i][j+7]=c[j+7][i]=in.motherVertexCovariance[i][j];
  std::cout<<std::setprecision(17)<<"{\"charge\":"<<charge<<",\"by\":"<<by<<",\"input\":[";
  for (int i=0;i<10;++i) std::cout<<(i?",":"")<<y[i];
  std::cout<<"],\"covariance\":[";
  for (int i=0;i<100;++i) std::cout<<(i?",":"")<<c[i/10][i%10];
  std::cout<<"],\"state\":[";
  for (int i=0;i<8;++i) std::cout<<(i?",":"")<<Get(state,i);
  std::cout<<"],\"state_covariance\":[";
  for (int i=0;i<64;++i) std::cout<<(i?",":"")<<state.covariance.m[i/8][i%8];
  ProductionState decay=state;assert(TransportToDecayVertex(decay));
  std::cout<<"],\"decay_state\":[";
  for (int i=0;i<8;++i) std::cout<<(i?",":"")<<Get(decay,i);
  std::cout<<"],\"chi2\":"<<state.pvConstraintChi2<<"}\n";
}
void EndpointCheck(int charge,double by,double flight) {
  ProductionState start;
  start.valid=true;start.status=kProductionValid;start.atProductionVertex=true;
  start.charge=charge;start.magneticFieldByTesla=by;
  start.positionCm=Vec3(.2,-.1,.5);start.momentumGeV=Vec3(.31,.17,1.4);
  start.energyGeV=3.4;start.flightParameterS=flight;start.covariance=Mat8();
  double lower[8][8]={{0.}};const double sd[8]={.1,.1,.1,.012,.012,.018,.016,.7};
  for(int i=0;i<8;++i){lower[i][i]=sd[i];for(int j=0;j<i;++j)lower[i][j]=.1*sd[i]*std::cos(i+j);}
  for(int i=0;i<8;++i)for(int j=0;j<8;++j)for(int k=0;k<8;++k)start.covariance.m[i][j]+=lower[i][k]*lower[j][k];
  ProductionState moved=start;
  assert(TransportToDecayVertex(moved));
  TrackState track,atDecay;track.positionCm=start.positionCm;track.momentumGeV=start.momentumGeV;track.charge=charge;
  assert(PropagateState(track,flight*Norm(track.momentumGeV),by,TransportConfig(),atDecay,false));
  for(int i=0;i<3;++i){Close(moved.positionCm[i],atDecay.positionCm[i]);Close(moved.momentumGeV[i],atDecay.momentumGeV[i]);}
  double j[8][8]={{0.}};
  for(int k=0;k<8;++k) {
    ProductionState plus=start,minus=start;const double step=1.e-5*std::max(1.,std::fabs(Get(start,k)));
    Set(plus,k,Get(start,k)+step);Set(minus,k,Get(start,k)-step);
    assert(TransportToDecayVertex(plus) && TransportToDecayVertex(minus));
    for(int i=0;i<8;++i)j[i][k]=(Get(plus,i)-Get(minus,i))/(2.*step);
  }
  for(int i=0;i<8;++i)for(int k=0;k<8;++k) {
    double numeric=0.;for(int a=0;a<8;++a)for(int b=0;b<8;++b)numeric+=j[i][a]*start.covariance.m[a][b]*j[k][b];
    const double scale=std::sqrt(moved.covariance.m[i][i]*moved.covariance.m[k][k]);
    assert(std::fabs(numeric-moved.covariance.m[i][k])<2.e-7*scale);
  }
  FlightObservables before,after;assert(GetFlightObservables(start,before)&&GetFlightObservables(moved,after));
  Close(before.massGeV,after.massGeV);Close(before.sigmaMassGeV,after.sigmaMassGeV);
  Close(before.ctauCm,after.ctauCm);Close(before.sigmaCtauCm,after.sigmaCtauCm);
  assert(TransportToDecayVertex(moved)); // idempotent
  assert(TransportToProductionVertex(moved));
  for(int i=0;i<8;++i){Close(Get(start,i),Get(moved,i));for(int k=0;k<8;++k)Close(start.covariance.m[i][k],moved.covariance.m[i][k]);}
}
}

int main(int argc,char** argv) {
  ProductionVertexConfig cfg;ProductionState out;
  const bool print=argc>1 && std::string(argv[1])=="--reference";
  for(int charge=-2;charge<=2;++charge)if(charge)for(int sign=-1;sign<=1;sign+=2)for(int path=-1;path<=1;path+=2) {
    const double by=sign*.5,flight=path*35.;
    const ProductionVertexInput in=Fixture(charge,by,flight);
    assert(FitProductionVertex(in,cfg,.7,1,out,charge,by));
    assert(out.valid&&out.atProductionVertex&&out.charge==charge&&out.magneticFieldByTesla==by);
    assert(out.pvConstraintNdf==2&&out.totalNdf==3&&IsPositiveSemidefinite(out.covariance));
    Close(out.totalChi2,out.pvConstraintChi2+.7);
    if(print)Print(in,out,charge,by);
    else EndpointCheck(charge,by,flight);
  }
  if(print)return 0;
  // Exact curved trajectory closes at the original PV; a straight pointing
  // hypothesis has a large residual even with the very same physical state.
  ProductionVertexInput in=Fixture(1,.5,80.,false);
  assert(FitProductionVertex(in,cfg,0.,1,out,1,.5));
  Close(out.flightParameterS,80.,1.e-7);Close(out.pvConstraintChi2,0.);
  for(int i=0;i<3;++i)Close(out.positionCm[i],in.vertex.positionCm[i]);
  ProductionState straight;
  assert(FitProductionVertex(in,cfg,0.,1,straight,0));
  assert(straight.pvConstraintChi2>1.);
  // Common PV/SV translation must not inflate the flight-parameter error.
  const ProductionState independent=out;
  for(int i=0;i<3;++i){in.mother.covariance.m[i][i]+=.04;in.vertex.covarianceCm2.m[i][i]+=.04;in.motherVertexCovariance[i][i]+=.04;}
  assert(FitProductionVertex(in,cfg,0.,1,out,1,.5));
  Close(out.covariance.m[7][7],independent.covariance.m[7][7]);
  for(int i=0;i<3;++i)Close(out.covariance.m[i][i],independent.covariance.m[i][i]+.04);
  // Exact-known PV is a legitimate singular 10D input.
  in=Fixture(1,.5,10.);in.vertex.covarianceCm2=Mat3();
  for(int i=0;i<7;++i)for(int j=0;j<3;++j)in.motherVertexCovariance[i][j]=0.;
  assert(FitProductionVertex(in,cfg,0.,1,out,1,.5));
  for(int i=0;i<3;++i)Close(out.positionCm[i],in.vertex.positionCm[i],1.e-7);
  const ProductionState knownPv=out;
  assert(TransportToDecayVertex(out)&&TransportToProductionVertex(out));
  for(int i=0;i<8;++i){Close(Get(out,i),Get(knownPv,i));for(int j=0;j<8;++j)Close(out.covariance.m[i][j],knownPv.covariance.m[i][j]);}
  // A momentum parallel to By does not curve, including for a charged mother.
  in=Fixture(1,.5,8.,false);in.mother.momentumGeV=Vec3(0.,1.,0.);
  in.mother.positionCm=in.vertex.positionCm+Vec3(0.,8.,0.);
  assert(FitProductionVertex(in,cfg,0.,1,out,1,.5));
  Close(out.flightParameterS,8.);Close(out.pvConstraintChi2,0.);
  for(int i=0;i<3;++i)Close(out.momentumGeV[i],in.mother.momentumGeV[i]);
  // Explicit zero field reproduces the neutral solver exactly; tiny field is continuous.
  in=Fixture(1,0.,15.);assert(FitProductionVertex(in,cfg,0.,1,straight));
  assert(FitProductionVertex(in,cfg,0.,1,out,1,0.));
  for(int i=0;i<8;++i){Close(Get(straight,i),Get(out,i),1.e-12);for(int j=0;j<8;++j)Close(straight.covariance.m[i][j],out.covariance.m[i][j],1.e-12);}
  assert(FitProductionVertex(in,cfg,0.,1,out,1,1.e-10));
  for(int i=0;i<8;++i){Close(Get(straight,i),Get(out,i),1.e-7);for(int j=0;j<8;++j)Close(straight.covariance.m[i][j],out.covariance.m[i][j],1.e-7);}
  assert(!FitProductionVertex(in,cfg,0.,1,out,1));
  assert(out.status==kProductionMissingMagneticField&&!out.valid&&std::isnan(out.flightParameterS));
  assert(!FitProductionVertex(in,cfg,0.,1,out,1,std::numeric_limits<double>::quiet_NaN()));
  assert(out.status==kProductionInvalidConfiguration);
  in=Fixture(1,.5,20.);cfg.maximumAbsFlightPathCm=1.;
  assert(!FitProductionVertex(in,cfg,0.,1,out,1,.5)&&out.status==kProductionFlightRangeExceeded);
  cfg=ProductionVertexConfig();cfg.maximumIterations=1;
  assert(!FitProductionVertex(in,cfg,0.,1,out,1,.5)&&out.status==kProductionNotConverged);
  cfg=ProductionVertexConfig();in.motherVertexCovariance[0][0]=100.;
  assert(!FitProductionVertex(in,cfg,0.,1,out,1,.5)&&out.status==kProductionInvalidJointCovariance);
  std::cout<<"Charged By PV fit, signs, correlated/known PV, S Jacobian, endpoint covariance and neutral limits passed\n";
}
