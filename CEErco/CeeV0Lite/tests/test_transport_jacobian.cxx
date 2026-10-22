#include "CeeV0LiteTransport.h"
#include "CeeV0LiteVertexFitter.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>

using namespace CeeV0Lite;
namespace {
// Independent symmetric Jacobi eigenvalue check. The core's unpivoted
// semidefinite Cholesky is not a reliable oracle for rank-five POCA matrices.
bool PositiveSemidefinite(const Mat6& c) {
  double a[6][6];
  for(int i=0;i<6;++i)for(int j=0;j<6;++j) {
    if(c.m[i][i]<0. || c.m[j][j]<0.)return false;
    const double scale=std::sqrt(c.m[i][i]*c.m[j][j]);
    a[i][j]=scale>0.?c.m[i][j]/scale:0.;
  }
  for(int iter=0;iter<300;++iter) {
    int p=0,q=1;double largest=0.;
    for(int i=0;i<6;++i)for(int j=i+1;j<6;++j)
      if(std::fabs(a[i][j])>largest){largest=std::fabs(a[i][j]);p=i;q=j;}
    if(largest<1.e-14)break;
    const double phi=.5*std::atan2(2.*a[p][q],a[q][q]-a[p][p]);
    const double cs=std::cos(phi),sn=std::sin(phi),ap=a[p][p],aq=a[q][q],off=a[p][q];
    for(int i=0;i<6;++i)if(i!=p&&i!=q) {
      const double ip=a[i][p],iq=a[i][q];
      a[i][p]=a[p][i]=cs*ip-sn*iq;a[i][q]=a[q][i]=sn*ip+cs*iq;
    }
    a[p][p]=cs*cs*ap-2.*sn*cs*off+sn*sn*aq;
    a[q][q]=sn*sn*ap+2.*sn*cs*off+cs*cs*aq;
    a[p][q]=a[q][p]=0.;
  }
  for(int i=0;i<6;++i)if(a[i][i]<-1.e-9)return false;
  return true;
}
double CovarianceDifference(const Mat6& a, const Mat6& b) {
  double worst=0.;
  for (int i=0;i<6;++i) for (int j=0;j<6;++j) {
    const double scale=std::sqrt(std::max(a.m[i][i],b.m[i][i])*
                                 std::max(a.m[j][j],b.m[j][j]));
    worst=std::max(worst,std::fabs(a.m[i][j]-b.m[i][j])/std::max(1.e-15,scale));
  }
  return worst;
}
TrackState State(std::mt19937& rng) {
  std::normal_distribution<double> normal;
  std::uniform_real_distribution<double> uniform(0.,1.);
  TrackState t;
  t.positionCm=Vec3(20.*normal(rng),20.*normal(rng),20.*normal(rng));
  t.momentumGeV=Vec3(normal(rng),normal(rng),normal(rng));
  t.momentumGeV=t.momentumGeV*(std::pow(10.,-1.7+2.3*uniform(rng))/Norm(t.momentumGeV));
  t.charge=(rng()%2)?1:-1;
  t.massGeV=0.13957039;
  double l[6][6]={{0.}};
  for(int i=0;i<6;++i) for(int j=0;j<=i;++j)
    l[i][j]=(i<3?.02:.001)*(i==j?1.:.15*normal(rng));
  for(int i=0;i<6;++i) for(int j=0;j<6;++j)
    for(int k=0;k<6;++k) t.covariance.m[i][j]+=l[i][k]*l[j][k];
  return t;
}
}
int main() {
  TransportConfig fast, reference;
  reference.useAnalyticFixedPathJacobian=false;
  // Tight roots let this test compare mathematical derivatives instead of
  // differentiating the stopping-error of the original root solver.
  fast.closestApproachToleranceCm=reference.closestApproachToleranceCm=1.e-10;
  fast.maximumClosestApproachIterations=reference.maximumClosestApproachIterations=30;
  std::mt19937 rng(20261014);
  std::uniform_real_distribution<double> uniform(-1.,1.);
  double worstFixed=0.,worstPoca=0.;int accepted=0;
  for(int trial=0;trial<1200;++trial) {
    TrackState in=State(rng),a,b;
    if(trial%11==0)in.charge=0;
    if(trial%13==0)in.momentumGeV=Vec3(0.,Norm(in.momentumGeV),0.);
    const double field=trial%3==0?0.:(trial%3==1?.5:-.5);
    const double path=150.*uniform(rng);
    assert(PropagateState(in,path,field,fast,a));
    assert(PropagateState(in,path,field,reference,b));
    assert(Norm(a.positionCm-b.positionCm)==0.);
    assert(Norm(a.momentumGeV-b.momentumGeV)==0.);
    worstFixed=std::max(worstFixed,CovarianceDifference(a.covariance,b.covariance));
    const Vec3 target=a.positionCm+Vec3(.1*uniform(rng),.1*uniform(rng),.1*uniform(rng));
    double sa=0.,sb=0.;ClosestApproachResult da,db;
    const bool oka=TransportToClosestApproach(in,target,field,fast,a,sa,&da);
    const bool okb=TransportToClosestApproach(in,target,field,reference,b,sb,&db);
    assert(oka==okb);assert(sa==sb);
    assert(da.converged==db.converged && da.hitPathLimit==db.hitPathLimit);
    if(!oka)continue;
    ++accepted;
    assert(Norm(a.positionCm-b.positionCm)==0.);
    assert(Norm(a.momentumGeV-b.momentumGeV)==0.);
    worstPoca=std::max(worstPoca,CovarianceDifference(a.covariance,b.covariance));
    assert(PositiveSemidefinite(a.covariance));
  }
  // Bounded/nonconverged diagnostic states and multi-turn ties must retain FD.
  TrackState line=State(rng),out;line.positionCm=Vec3();line.momentumGeV=Vec3(.05,0.,0.);
  double s;ClosestApproachResult d;
  TransportConfig bounded=fast;bounded.maximumAbsPathCm=1.;
  bounded.requireClosestApproachConvergence=false;
  assert(TransportToClosestApproach(line,Vec3(5.,1.,0.),0.,bounded,out,s,&d));
  assert(d.hitPathLimit);
  line.charge=1;
  const double radius=.05/(.00299792458*.5);
  assert(TransportToClosestApproach(line,Vec3(0.,0.,2.*radius),.5,fast,out,s,&d));
  assert(d.converged);
  // Low-momentum orthogonal straight lines have a valid exact crossing.
  TrackState other=line;other.positionCm=Vec3(1.,1.,0.);other.momentumGeV=Vec3(0.,.1,0.);
  VertexFitResult fit;assert(VertexFitter().Fit(line,other,0.,fit));
  assert(Norm(fit.vertexCm-Vec3(1.,0.,0.))<1.e-8);
  assert(fit.daughterDcaCm<1.e-8);
  other.momentumGeV=Vec3(.1,0.,0.);
  assert(!VertexFitter().Fit(line,other,0.,fit));
  assert(fit.status==kVertexFitDegenerateGeometry);
  std::cout<<"transport derivative comparison: accepted="<<accepted

    <<" max scaled covariance difference fixed="<<worstFixed
    <<" POCA="<<worstPoca<<std::endl;

  assert(worstFixed<2.e-5);
  assert(worstPoca<2.e-4);
}
