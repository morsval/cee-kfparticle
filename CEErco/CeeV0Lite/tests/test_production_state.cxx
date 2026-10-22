#include "CeeV0LiteProductionVertex.h"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using namespace CeeV0Lite;
namespace {
void Close(double a, double b, double tol=1.e-9) {
  assert(std::isfinite(a) && std::isfinite(b));
  assert(std::fabs(a-b) < tol*std::max(1.,std::max(std::fabs(a),std::fabs(b))));
}
ProductionVertexInput Axial(double z=10.) {
  ProductionVertexInput in;
  in.mother.positionCm=Vec3(.3,-.2,z);
  in.mother.momentumGeV=Vec3(0.,0.,2.);
  in.mother.energyGeV=std::sqrt(5.);
  in.mother.covarianceValid=true;
  in.mother.covarianceStatus=kMotherCovarianceValid;
  for(int i=0;i<3;++i) in.mother.covariance.m[i][i]=.04;
  in.mother.covariance.m[6][6]=.0001;
  in.vertex.covarianceCm2=Mat3::Identity(.01);
  return in;
}
ProductionVertexInput General() {
  ProductionVertexInput in;
  in.mother.positionCm=Vec3(.4,-.12,4.2);
  in.mother.momentumGeV=Vec3(.08,-.025,1.1);
  in.mother.energyGeV=std::sqrt(Norm2(in.mother.momentumGeV)+1.115*1.115);
  in.vertex.positionCm=Vec3(.015,-.01,.03);
  in.mother.covarianceValid=true;
  const double sigma[10]={.1,.12,.15,.02,.015,.025,.025,.03,.05,.06};
  double lower[10][10]={{0.}}, c[10][10]={{0.}};
  for(int i=0;i<10;++i) {
    lower[i][i]=sigma[i];
    for(int j=0;j<i;++j) lower[i][j]=sigma[i]*.035*std::cos(double(i+j));
  }
  for(int i=0;i<10;++i) for(int j=0;j<10;++j)
    for(int k=0;k<10;++k)c[i][j]+=lower[i][k]*lower[j][k];
  for(int i=0;i<7;++i)for(int j=0;j<7;++j)in.mother.covariance.m[i][j]=c[i][j];
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)in.vertex.covarianceCm2.m[i][j]=c[i+7][j+7];
  for(int i=0;i<7;++i)for(int j=0;j<3;++j)in.motherVertexCovariance[i][j]=c[i][j+7];
  return in;
}
void PrintReference(const ProductionVertexInput& in, const ProductionState& out) {
  double x[10], q[8];
  for(int i=0;i<3;++i) {
    x[i]=in.mother.positionCm[i];x[i+3]=in.mother.momentumGeV[i];x[i+7]=in.vertex.positionCm[i];
    q[i]=out.positionCm[i];q[i+3]=out.momentumGeV[i];
  }
  x[6]=in.mother.energyGeV;q[6]=out.energyGeV;q[7]=out.flightParameterS;
  double c[10][10]={{0.}};
  for(int i=0;i<7;++i)for(int j=0;j<7;++j)c[i][j]=in.mother.covariance.m[i][j];
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)c[i+7][j+7]=in.vertex.covarianceCm2.m[i][j];
  for(int i=0;i<7;++i)for(int j=0;j<3;++j)c[i][j+7]=c[j+7][i]=in.motherVertexCovariance[i][j];
  std::cout<<std::setprecision(17)<<"{\"input\":[";
  for(int i=0;i<10;++i)std::cout<<(i?",":"")<<x[i];
  std::cout<<"],\"covariance\":[";
  for(int i=0;i<100;++i)std::cout<<(i?",":"")<<c[i/10][i%10];
  std::cout<<"],\"state\":[";
  for(int i=0;i<8;++i)std::cout<<(i?",":"")<<q[i];
  std::cout<<"],\"state_covariance\":[";
  for(int i=0;i<64;++i)std::cout<<(i?",":"")<<out.covariance.m[i/8][i%8];
  std::cout<<"],\"chi2\":"<<out.pvConstraintChi2<<"}\n";
}
}

int main(int argc,char** argv) {
  ProductionVertexConfig cfg;
  ProductionState state;
  ProductionVertexInput in=General();
  assert(FitProductionVertex(in,cfg,1.,1,state));
  if(argc>1 && std::string(argv[1])=="--reference") {PrintReference(in,state);return 0;}
  assert(state.valid && state.atProductionVertex && state.pvConstraintNdf==2 && state.totalNdf==3);
  Close(state.totalChi2,state.pvConstraintChi2+1.);
  const ProductionState saved=state;
  FlightObservables before,after;
  assert(GetFlightObservables(state,before));
  assert(TransportToDecayVertex(state));
  assert(!state.atProductionVertex);
  assert(GetFlightObservables(state,after));
  Close(before.ctauCm,after.ctauCm);Close(before.sigmaCtauCm,after.sigmaCtauCm);
  assert(TransportToDecayVertex(state)); // idempotent
  assert(TransportToProductionVertex(state));
  for(int i=0;i<3;++i)Close(saved.positionCm[i],state.positionCm[i]);
  for(int i=0;i<8;++i)for(int j=0;j<8;++j)Close(saved.covariance.m[i][j],state.covariance.m[i][j]);
  Close(saved.flightParameterS,state.flightParameterS);
  assert(IsPositiveSemidefinite(state.covariance));

  // Independent analytic result with known momentum. x/y are weighted means;
  // z has no pointing information. PV-z uncertainty anti-correlates with S.
  in=Axial();
  const MotherState original=in.mother;
  assert(FitProductionVertex(in,cfg,0.,1,state));
  Close(state.pvConstraintChi2,(.3*.3+.2*.2)/.05);
  Close(state.positionCm[0],.3*.01/.05);Close(state.positionCm[1],-.2*.01/.05);
  Close(state.positionCm[2],0.);Close(state.covariance.m[0][0],.04*.01/.05);
  Close(state.flightParameterS,5.);Close(state.covariance.m[7][7],.05/4.);
  Close(state.covariance.m[2][7],-.01/2.);
  assert(GetFlightObservables(state,after));
  Close(after.decayLengthCm,10.);Close(after.sigmaDecayLengthCm,std::sqrt(.05));
  Close(after.ctauCm,5.);Close(after.sigmaCtauCm,std::sqrt(.05/4.+25.*5.*.0001));
  assert(TransportToDecayVertex(state));
  Close(state.positionCm[2],10.);Close(state.covariance.m[2][7],.04/2.);
  Close(in.mother.positionCm[0],original.positionCm[0]);
  Close(in.mother.covariance.m[0][0],original.covariance.m[0][0]);
  ProductionState more;
  cfg.maximumIterations=100;
  assert(FitProductionVertex(in,cfg,0.,1,more));
  assert(TransportToDecayVertex(more));
  for(int i=0;i<8;++i)for(int j=0;j<8;++j)Close(more.covariance.m[i][j],state.covariance.m[i][j]);

  // Shared translation of PV and SV must cancel from their separation. This
  // checks cross-block signs and prevents counting PV uncertainty twice.
  in=Axial();
  for(int i=0;i<3;++i) {
    in.mother.covariance.m[i][i]+=.01;
    in.motherVertexCovariance[i][i]=.01;
  }
  assert(FitProductionVertex(in,cfg,0.,1,state));
  Close(state.covariance.m[7][7],.04/4.);
  Close(state.pvConstraintChi2,(.3*.3+.2*.2)/.04);
  in=Axial(-10.);
  assert(FitProductionVertex(in,cfg,0.,1,state));
  assert(GetFlightObservables(state,after));
  Close(state.flightParameterS,-5.);Close(after.ctauCm,-5.);

  in=Axial();in.vertex.covarianceCm2=Mat3();
  assert(FitProductionVertex(in,cfg,0.,1,state));
  Close(Norm(state.positionCm-in.vertex.positionCm),0.);
  assert(!FitProductionVertex(in,cfg,0.,1,state,-1));
  assert(state.status==kProductionMissingMagneticField && !state.valid && std::isnan(state.flightParameterS));
  in=Axial();in.motherVertexCovariance[0][0]=10.;
  assert(!FitProductionVertex(in,cfg,0.,1,state));
  assert(state.status==kProductionInvalidJointCovariance);
  in=General();cfg.maximumIterations=1;
  assert(!FitProductionVertex(in,cfg,0.,1,state));
  assert(state.status==kProductionNotConverged && !state.valid);
  assert(!GetFlightObservables(state,after) && std::isnan(after.ctauCm));
  std::cout<<"Neutral 8D constraint, analytic covariance, correlated PV, lifetime and endpoint tests passed\n";
}
