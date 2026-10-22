#include "CeeV0LiteProductionVertex.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
using namespace CeeV0Lite;
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kMagneticRotation = 0.00299792458; // cm, GeV/c, Tesla

// Backward flight from decay to production:
// rP=rD-A(S)pD, pP=R(S)pD, R=dA/dS. The signs follow dp/ds=q v x By.
// This S is arc length / |p|; unlike fixed-length daughter transport, the
// rotation angle k*S does not contain another division by momentum.
bool HelixMatrices(double s, int charge, double by,
                   Mat3& a, Mat3& r, Mat3& derivative)
{
  const double k=kMagneticRotation*double(charge)*by, angle=k*s;
  if (!std::isfinite(s) || !std::isfinite(k) || !std::isfinite(angle)) return false;
  const double sn=std::sin(angle), co=std::cos(angle);
  double sinc, cosc;
  if (std::fabs(angle)<1.e-4) {
    const double u=angle*angle;
    sinc=1.-u/6.+u*u/120.-u*u*u/5040.;
    cosc=angle*(.5-u/24.+u*u/720.-u*u*u/40320.);
  } else {
    sinc=sn/angle;
    const double half=std::sin(.5*angle);
    cosc=2.*half*half/angle;
  }
  a=Mat3(); r=Mat3(); derivative=Mat3();
  a.m[0][0]=a.m[2][2]=s*sinc; a.m[1][1]=s;
  a.m[0][2]=s*cosc; a.m[2][0]=-a.m[0][2];
  r.m[0][0]=r.m[2][2]=co; r.m[1][1]=1.;
  r.m[0][2]=sn; r.m[2][0]=-sn;
  derivative.m[0][0]=derivative.m[2][2]=-k*sn;
  derivative.m[0][2]=k*co; derivative.m[2][0]=-k*co;
  for (int i=0;i<3;++i) for (int j=0;j<3;++j)
    if (!std::isfinite(a.m[i][j]) || !std::isfinite(r.m[i][j]) ||
        !std::isfinite(derivative.m[i][j])) return false;
  return true;
}

Vec3 Product(const Mat3& a, const Vec3& p)
{
  Vec3 result;
  for (int i=0;i<3;++i) for (int j=0;j<3;++j) result[i]+=a.m[i][j]*p[j];
  return result;
}

// Three geometrical equations with one unmeasured flight parameter. Project
// the LINEARIZED system perpendicular to d(Ap)/dS, removing that nuisance.
// This needs only a 2x2 inverse, also allowing a singular 10D prior (e.g. a
// known PV). No covariance inverse or fictitious error for S is introduced.
bool ChargedConstraint(const double x[10], double s, int charge, double by,
                       double g[3], double h[3][10], double ht[2][10],
                       Vec3 basis[2], Vec3& velocity)
{
  Mat3 a,r,dr;
  if (!HelixMatrices(s,charge,by,a,r,dr)) return false;
  const Vec3 p(x[3],x[4],x[5]);
  const Vec3 residual=Vec3(x[0]-x[7],x[1]-x[8],x[2]-x[9])-Product(a,p);
  velocity=Product(r,p);
  Vec3 direction;
  if (!IsFinite(residual) || !Normalize(velocity,direction)) return false;
  MakeTransverseBasis(direction,basis[0],basis[1]);
  for (int i=0;i<3;++i) {
    g[i]=residual[i];
    for (int j=0;j<10;++j) h[i][j]=0.;
    h[i][i]=1.; h[i][i+7]=-1.;
    for (int j=0;j<3;++j) h[i][j+3]=-a.m[i][j];
  }
  for (int i=0;i<2;++i) for (int j=0;j<10;++j) {
    ht[i][j]=0.;
    for (int m=0;m<3;++m) ht[i][j]+=basis[i][m]*h[m][j];
  }
  return true;
}

template<int N> bool ValidCovariance(const double (&c)[N][N], double (*factor)[N]=0)
{
  double scale[N], a[N][N] = {{0.}}, l[N][N] = {{0.}};
  const double tolerance = 1.e-10;
  for (int i=0; i<N; ++i) {
    if (!std::isfinite(c[i][i]) || c[i][i]<0.) return false;
    scale[i]=std::sqrt(c[i][i]);
  }
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
    if (!std::isfinite(c[i][j])) return false;
    if (scale[i]==0. || scale[j]==0.) {
      if (c[i][j]!=0. || c[j][i]!=0.) return false;
    } else {
      a[i][j]=c[i][j]/(scale[i]*scale[j]);
      if (std::fabs(a[i][j]-c[j][i]/(scale[i]*scale[j]))>tolerance)
        return false;
    }
  }
  for (int i=0; i<N; ++i) for (int j=0; j<=i; ++j) {
    double v=a[i][j];
    for (int k=0; k<j; ++k) v-=l[i][k]*l[j][k];
    if (i==j) {
      if (v < -tolerance) return false;
      l[i][j]=v>tolerance ? std::sqrt(v) : 0.;
    } else if (l[j][j]>0.) l[i][j]=v/l[j][j];
    else if (std::fabs(v)>tolerance) return false;
  }
  if (factor) for (int i=0;i<N;++i) for (int j=0;j<N;++j) factor[i][j]=scale[i]*l[i][j];
  return true;
}

// n,e1,e2 are fixed from the prior momentum. The plane-intersection residual
// g=e.d - (n.d)/(n.p) e.p vanishes precisely when d is parallel to p.
// Unlike differentiating a freshly chosen transverse basis, this gives an
// explicit smooth Jacobian even away from the constraint surface.
bool Constraint(const double x[10], const Vec3& n, const Vec3 basis[2],
                double g[2], double h[2][10])
{
  const Vec3 p(x[3],x[4],x[5]), d(x[0]-x[7],x[1]-x[8],x[2]-x[9]);
  const double np=Dot(n,p), momentum=Norm(p);
  if (!(momentum>1.e-12) || !(np>1.e-4*momentum)) return false;
  const double s=Dot(n,d)/np;
  for (int a=0; a<2; ++a) {
    const double ep=Dot(basis[a],p);
    g[a]=Dot(basis[a],d)-s*ep;
    const Vec3 w=basis[a]-(ep/np)*n;
    for (int j=0; j<10; ++j) h[a][j]=0.;
    for (int i=0; i<3; ++i) {
      h[a][i]=w[i]; h[a][i+3]=-s*w[i]; h[a][i+7]=-w[i];
    }
    if (!std::isfinite(g[a])) return false;
  }
  return true;
}

bool Gain(const double c[10][10], const double h[2][10],
          double k[10][2], double inverse[2][2])
{
  double cht[10][2]={{0.}}, w[2][2]={{0.}};
  for (int i=0; i<10; ++i) for (int a=0; a<2; ++a)
    for (int j=0; j<10; ++j) cht[i][a]+=c[i][j]*h[a][j];
  for (int a=0; a<2; ++a) for (int b=0; b<2; ++b)
    for (int i=0; i<10; ++i) w[a][b]+=h[a][i]*cht[i][b];
  const double scale=std::max(w[0][0],w[1][1]);
  if (!(w[0][0]>0.) || !(w[1][1]>0.) || !std::isfinite(scale) ||
      !(w[0][0]*w[1][1]-w[0][1]*w[1][0]>1.e-14*scale*scale) ||
      !InvertSymmetric2(w,inverse)) return false;
  for (int i=0; i<10; ++i) for (int a=0; a<2; ++a)
    k[i][a]=cht[i][0]*inverse[0][a]+cht[i][1]*inverse[1][a];
  return true;
}

bool TransportEndpoint(ProductionState& state, bool production)
{
  double lower[8][8];
  if (!state.valid || !IsFinite(state.positionCm) ||
      !IsFinite(state.momentumGeV) || !std::isfinite(state.energyGeV) ||
      !std::isfinite(state.flightParameterS) ||
      !std::isfinite(state.magneticFieldByTesla) ||
      !ValidCovariance(state.covariance.m,lower)) return false;
  if (state.atProductionVertex==production) return true;
  const double sign=production ? 1. : -1.;
  Mat3 a,r,dr;
  if (!HelixMatrices(sign*state.flightParameterS,state.charge,
                     state.magneticFieldByTesla,a,r,dr)) return false;
  const Vec3 momentum=Product(r,state.momentumGeV);
  const Vec3 momentumDerivative=Product(dr,state.momentumGeV);
  double j[8][8]={{0.}};
  for (int i=0; i<8; ++i) j[i][i]=1.;
  for (int i=0; i<3; ++i) {
    for (int b=0;b<3;++b) {
      j[i][b+3]=-a.m[i][b]; j[i+3][b+3]=r.m[i][b];
    }
    j[i][7]=-sign*momentum[i];
    j[i+3][7]=sign*momentumDerivative[i];
  }
  // C=L L^T -> (J L)(J L)^T. The Gram form avoids negative position
  // variances from cancellation when an endpoint is exactly known.
  double jl[8][8]={{0.}};
  for (int i=0;i<8;++i) for (int a=0;a<8;++a) for (int b=0;b<8;++b)
    jl[i][a]+=j[i][b]*lower[b][a];
  Mat8 c;
  for (int i=0;i<8;++i) for (int k=0;k<8;++k) for (int a=0;a<8;++a)
    c.m[i][k]+=jl[i][a]*jl[k][a];
  Symmetrize(c);
  const Vec3 position=state.positionCm-Product(a,state.momentumGeV);
  if (!IsFinite(position) || !IsFinite(momentum) || !IsPositiveSemidefinite(c)) return false;
  state.positionCm=position;
  state.momentumGeV=momentum;
  state.covariance=c;
  state.atProductionVertex=production;
  return true;
}

bool Error(const double gradient[8], const Mat8& c, double& error)
{
  double value=0., scale=0.;
  for (int i=0; i<8; ++i) for (int j=0; j<8; ++j) {
    const double term=gradient[i]*c.m[i][j]*gradient[j];
    value+=term; scale+=std::fabs(term);
  }
  if (!std::isfinite(value) || value < -64.*std::numeric_limits<double>::epsilon()*scale)
    return false;
  error=std::sqrt(std::max(0.,value));
  return true;
}
} // namespace

namespace CeeV0Lite {

ProductionVertexConfig::ProductionVertexConfig()
  : maximumIterations(30), relativeStateTolerance(1.e-8),
    constraintToleranceCm(1.e-7), pvDerivativeStepCm(1.e-3),
    pvDerivativeStabilityTolerance(0.05), maximumAbsFlightPathCm(1000.),
    maximumAbsBendingAngleRad(3.14159265358979323846) {}

ProductionVertexInput::ProductionVertexInput()
{
  for (int i=0; i<7; ++i) for (int j=0; j<3; ++j)
    motherVertexCovariance[i][j]=0.;
}

ProductionState::ProductionState()
  : valid(false), status(kProductionNotRun),
    correlationModel(kProductionUserSuppliedJointCovariance),
    atProductionVertex(false), iterations(0),
    charge(0), magneticFieldByTesla(kNaN),
    positionCm(kNaN,kNaN,kNaN), momentumGeV(kNaN,kNaN,kNaN),
    energyGeV(kNaN), flightParameterS(kNaN), pvConstraintChi2(kNaN),
    pvConstraintNdf(0), totalChi2(kNaN), totalNdf(0)
{
  for (int i=0; i<8; ++i) for (int j=0; j<8; ++j) covariance.m[i][j]=kNaN;
}

FlightObservables::FlightObservables()
  : valid(false), massGeV(kNaN), sigmaMassGeV(kNaN), decayLengthCm(kNaN),
    sigmaDecayLengthCm(kNaN), ctauCm(kNaN), sigmaCtauCm(kNaN) {}

static bool FitNeutralProductionVertex(const ProductionVertexInput& input,
                         const ProductionVertexConfig& config,
                         double geometricChi2, int geometricNdf,
                         ProductionState& output, int motherCharge)
{
  output=ProductionState();
  if (motherCharge!=0) { output.status=kProductionUnsupportedCharge; return false; }
  if (config.maximumIterations<1 || !(config.relativeStateTolerance>0.) ||
      !std::isfinite(config.relativeStateTolerance) ||
      !(config.constraintToleranceCm>0.) || !std::isfinite(config.constraintToleranceCm) ||
      !(geometricChi2>=0.) || !std::isfinite(geometricChi2) || geometricNdf<0) {
    output.status=kProductionInvalidConfiguration; return false;
  }
  if (!input.mother.covarianceValid || !IsFinite(input.mother.positionCm) ||
      !IsFinite(input.mother.momentumGeV) || !IsFinite(input.vertex.positionCm) ||
      !(input.mother.energyGeV>0.) || !std::isfinite(input.mother.energyGeV) ||
      !(input.mother.energyGeV*input.mother.energyGeV>Norm2(input.mother.momentumGeV))) {
    output.status=kProductionInvalidMother; return false;
  }
  double x0[10], x[10], c[10][10]={{0.}};
  for (int i=0; i<3; ++i) {
    x0[i]=input.mother.positionCm[i]; x0[i+3]=input.mother.momentumGeV[i];
    x0[i+7]=input.vertex.positionCm[i];
  }
  x0[6]=input.mother.energyGeV;
  for (int i=0; i<10; ++i) x[i]=x0[i];
  for (int i=0; i<7; ++i) for (int j=0; j<7; ++j) c[i][j]=input.mother.covariance.m[i][j];
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) c[i+7][j+7]=input.vertex.covarianceCm2.m[i][j];
  for (int i=0; i<7; ++i) for (int j=0; j<3; ++j)
    c[i][j+7]=c[j+7][i]=input.motherVertexCovariance[i][j];
  if (!ValidCovariance(c)) { output.status=kProductionInvalidJointCovariance; return false; }
  Vec3 n, basis[2];
  if (!Normalize(input.mother.momentumGeV,n)) { output.status=kProductionInvalidMother; return false; }
  MakeTransverseBasis(n,basis[0],basis[1]);
  double h[2][10], g[2], inverse[2][2], k[10][2];
  bool converged=false;
  for (int iteration=0; iteration<config.maximumIterations; ++iteration) {
    output.iterations=iteration+1;
    if (!Constraint(x,n,basis,g,h) || !Gain(c,h,k,inverse)) {
      output.status=kProductionSingularConstraint; return false;
    }
    double innovation[2]={g[0],g[1]};
    for (int a=0; a<2; ++a) for (int i=0; i<10; ++i)
      innovation[a]+=h[a][i]*(x0[i]-x[i]);
    double next[10], change=0.;
    for (int i=0; i<10; ++i) {
      next[i]=x0[i]-k[i][0]*innovation[0]-k[i][1]*innovation[1];
      if (!std::isfinite(next[i])) { output.status=kProductionUnphysicalState; return false; }
      change=std::max(change,std::fabs(next[i]-x[i])/(1.+std::fabs(x[i])));
      x[i]=next[i];
    }
    if (!Constraint(x,n,basis,g,h)) { output.status=kProductionSingularConstraint; return false; }
    if (change<config.relativeStateTolerance &&
        std::hypot(g[0],g[1])<config.constraintToleranceCm) { converged=true; break; }
  }
  if (!converged) { output.status=kProductionNotConverged; return false; }
  if (!Gain(c,h,k,inverse)) { output.status=kProductionSingularConstraint; return false; }

  // Joseph form with an exact geometrical constraint (R=0) avoids subtracting
  // nearly equal covariances. Every iteration above used c, the original prior.
  double a[10][10]={{0.}}, ac[10][10]={{0.}}, posterior[10][10]={{0.}};
  for (int i=0; i<10; ++i) for (int j=0; j<10; ++j) {
    a[i][j]=(i==j ? 1. : 0.)-k[i][0]*h[0][j]-k[i][1]*h[1][j];
  }
  // a must be complete before computing A C.
  for (int i=0; i<10; ++i) for (int j=0; j<10; ++j) {
    ac[i][j]=0.;
    for (int b=0; b<10; ++b) ac[i][j]+=a[i][b]*c[b][j];
  }
  for (int i=0; i<10; ++i) for (int j=0; j<10; ++j)
    for (int b=0; b<10; ++b) posterior[i][j]+=ac[i][b]*a[j][b];
  const Vec3 p(x[3],x[4],x[5]), d(x[0]-x[7],x[1]-x[8],x[2]-x[9]);
  const double p2=Norm2(p), s=Dot(p,d)/p2;
  if (!(p2>1.e-24) || !(x[6]>0.) || !(x[6]*x[6]>p2) || !std::isfinite(s)) {
    output.status=kProductionUnphysicalState; return false;
  }
  // First form q at decay. Its first seven variables are the fitted mother;
  // the S derivative uses BOTH endpoints of the correlated joint posterior.
  double t[8][10]={{0.}};
  for (int i=0; i<7; ++i) t[i][i]=1.;
  for (int i=0; i<3; ++i) {
    t[7][i]=p[i]/p2; t[7][i+7]=-p[i]/p2;
    t[7][i+3]=(d[i]-2.*s*p[i])/p2;
  }
  ProductionState state;
  state.charge=0; state.magneticFieldByTesla=0.;
  state.positionCm=Vec3(x[0],x[1],x[2]); state.momentumGeV=p;
  state.energyGeV=x[6]; state.flightParameterS=s; state.covariance=Mat8();
  for (int i=0; i<8; ++i) for (int j=0; j<8; ++j)
    for (int b=0; b<10; ++b) for (int e=0; e<10; ++e)
      state.covariance.m[i][j]+=t[i][b]*posterior[b][e]*t[j][e];
  Symmetrize(state.covariance);
  double innovation[2]={g[0],g[1]};
  for (int b=0; b<2; ++b) for (int i=0; i<10; ++i)
    innovation[b]+=h[b][i]*(x0[i]-x[i]);
  state.pvConstraintChi2=innovation[0]*(inverse[0][0]*innovation[0]+inverse[0][1]*innovation[1])+
    innovation[1]*(inverse[1][0]*innovation[0]+inverse[1][1]*innovation[1]);
  state.pvConstraintNdf=2; state.totalChi2=geometricChi2+state.pvConstraintChi2;
  state.totalNdf=geometricNdf+2; state.iterations=output.iterations;
  state.status=kProductionValid; state.valid=true;
  FlightObservables flight;
  if (!(state.pvConstraintChi2>=0.) || !std::isfinite(state.totalChi2) ||
      !GetFlightObservables(state,flight) || !TransportToProductionVertex(state)) {
    output.status=kProductionUnphysicalState; return false;
  }
  output=state;
  return true;
}

static bool FitChargedProductionVertex(const ProductionVertexInput& input,
                                       const ProductionVertexConfig& config,
                                       double geometricChi2, int geometricNdf,
                                       ProductionState& output, int charge, double by)
{
  output=ProductionState();
  if (config.maximumIterations<1 || !(config.relativeStateTolerance>0.) ||
      !std::isfinite(config.relativeStateTolerance) ||
      !(config.constraintToleranceCm>0.) || !std::isfinite(config.constraintToleranceCm) ||
      !(config.maximumAbsFlightPathCm>0.) || !std::isfinite(config.maximumAbsFlightPathCm) ||
      !(config.maximumAbsBendingAngleRad>0.) ||
      !(config.maximumAbsBendingAngleRad<=3.14159265358979323846) ||
      !(geometricChi2>=0.) || !std::isfinite(geometricChi2) || geometricNdf<0) {
    output.status=kProductionInvalidConfiguration; return false;
  }
  if (!input.mother.covarianceValid || !IsFinite(input.mother.positionCm) ||
      !IsFinite(input.mother.momentumGeV) || !IsFinite(input.vertex.positionCm) ||
      !(input.mother.energyGeV>0.) || !std::isfinite(input.mother.energyGeV) ||
      !(input.mother.energyGeV*input.mother.energyGeV>Norm2(input.mother.momentumGeV))) {
    output.status=kProductionInvalidMother; return false;
  }
  double x0[10],x[10],c[10][10]={{0.}};
  for (int i=0;i<3;++i) {
    x0[i]=input.mother.positionCm[i]; x0[i+3]=input.mother.momentumGeV[i];
    x0[i+7]=input.vertex.positionCm[i];
  }
  x0[6]=input.mother.energyGeV;
  for (int i=0;i<10;++i) x[i]=x0[i];
  for (int i=0;i<7;++i) for (int j=0;j<7;++j) c[i][j]=input.mother.covariance.m[i][j];
  for (int i=0;i<3;++i) for (int j=0;j<3;++j) c[i+7][j+7]=input.vertex.covarianceCm2.m[i][j];
  for (int i=0;i<7;++i) for (int j=0;j<3;++j)
    c[i][j+7]=c[j+7][i]=input.motherVertexCovariance[i][j];
  if (!ValidCovariance(c)) { output.status=kProductionInvalidJointCovariance; return false; }
  const double momentum2=Norm2(input.mother.momentumGeV);
  if (!(momentum2>1.e-24)) { output.status=kProductionInvalidMother; return false; }
  // Straight projection is only an initial guess, never a charged constraint.
  double s=Dot(input.mother.positionCm-input.vertex.positionCm,
               input.mother.momentumGeV)/momentum2;
  const double curvature=kMagneticRotation*double(charge)*by;
  double g[3],h[3][10],ht[2][10],gain[10][2],inverse[2][2];
  Vec3 basis[2],velocity;
  bool converged=false;
  for (int iteration=0;iteration<config.maximumIterations;++iteration) {
    output.iterations=iteration+1;
    const double p=Norm(Vec3(x[3],x[4],x[5]));
    if (!std::isfinite(s) || !std::isfinite(p) ||
        std::fabs(s)*p>config.maximumAbsFlightPathCm ||
        std::fabs(curvature*s)>=config.maximumAbsBendingAngleRad) {
      output.status=kProductionFlightRangeExceeded; return false;
    }
    if (!ChargedConstraint(x,s,charge,by,g,h,ht,basis,velocity) ||
        !Gain(c,ht,gain,inverse)) { output.status=kProductionSingularConstraint; return false; }
    double b[3]={g[0],g[1],g[2]}, bt[2]={0.,0.};
    for (int a=0;a<3;++a) for (int i=0;i<10;++i) b[a]+=h[a][i]*(x0[i]-x[i]);
    for (int a=0;a<2;++a) for (int i=0;i<3;++i) bt[a]+=basis[a][i]*b[i];
    double delta[10],next[10],change=0.;
    for (int i=0;i<10;++i) {
      delta[i]=-gain[i][0]*bt[0]-gain[i][1]*bt[1];
      next[i]=x0[i]+delta[i];
      if (!std::isfinite(next[i])) { output.status=kProductionUnphysicalState; return false; }
      change=std::max(change,std::fabs(next[i]-x[i])/(1.+std::fabs(x[i])));
    }
    Vec3 remaining(b[0],b[1],b[2]);
    for (int a=0;a<3;++a) for (int i=0;i<10;++i) remaining[a]+=h[a][i]*delta[i];
    const double ds=Dot(velocity,remaining)/Norm2(velocity);
    const double scaledDs=std::fabs(ds)/(1.+std::fabs(s));
    s+=ds;
    for (int i=0;i<10;++i) x[i]=next[i];
    if (!std::isfinite(s) || std::fabs(s)*Norm(Vec3(x[3],x[4],x[5]))>config.maximumAbsFlightPathCm ||
        std::fabs(curvature*s)>=config.maximumAbsBendingAngleRad) {
      output.status=kProductionFlightRangeExceeded; return false;
    }
    if (!ChargedConstraint(x,s,charge,by,g,h,ht,basis,velocity)) {
      output.status=kProductionSingularConstraint; return false;
    }
    if (change<config.relativeStateTolerance && scaledDs<config.relativeStateTolerance &&
        Norm(Vec3(g[0],g[1],g[2]))<config.constraintToleranceCm) { converged=true; break; }
  }
  if (!converged) { output.status=kProductionNotConverged; return false; }
  if (!Gain(c,ht,gain,inverse)) { output.status=kProductionSingularConstraint; return false; }
  const Vec3 momentum(x[3],x[4],x[5]);
  if (!(x[6]>0.) || !(x[6]*x[6]>Norm2(momentum))) {
    output.status=kProductionUnphysicalState; return false;
  }
  // Jx maps prior fluctuations to the constrained 10D state. The S row follows
  // H dx-v dS=0: Js=v^T H Jx / |v|^2, retaining all PV/mother correlations.
  double jx[10][10]={{0.}}, t[8][10]={{0.}};
  for (int i=0;i<10;++i) for (int j=0;j<10;++j)
    jx[i][j]=(i==j ? 1. : 0.)-gain[i][0]*ht[0][j]-gain[i][1]*ht[1][j];
  for (int i=0;i<7;++i) for (int j=0;j<10;++j) t[i][j]=jx[i][j];
  const double v2=Norm2(velocity);
  for (int j=0;j<10;++j) for (int a=0;a<3;++a) for (int i=0;i<10;++i)
    t[7][j]+=velocity[a]*h[a][i]*jx[i][j]/v2;
  // Form the production endpoint directly from the same prior Jacobian.
  // Its position is the fitted PV, avoiding a subtraction of correlated
  // coordinates when that PV has exactly zero measurement variance.
  Mat3 arc,rotation,rotationDerivative;
  if (!HelixMatrices(s,charge,by,arc,rotation,rotationDerivative)) {
    output.status=kProductionUnphysicalState; return false;
  }
  const Vec3 productionMomentum=Product(rotation,momentum);
  const Vec3 dpds=Product(rotationDerivative,momentum);
  double tp[8][10]={{0.}};
  for (int j=0;j<10;++j) {
    for (int i=0;i<3;++i) {
      tp[i][j]=jx[i+7][j];
      tp[i+3][j]=dpds[i]*t[7][j];
      for (int a=0;a<3;++a) tp[i+3][j]+=rotation.m[i][a]*jx[a+3][j];
    }
    tp[6][j]=jx[6][j]; tp[7][j]=t[7][j];
  }
  ProductionState state;
  state.charge=charge; state.magneticFieldByTesla=by; state.atProductionVertex=true;
  state.positionCm=Vec3(x[7],x[8],x[9]); state.momentumGeV=productionMomentum;
  state.energyGeV=x[6]; state.flightParameterS=s; state.covariance=Mat8();
  for (int i=0;i<8;++i) for (int j=0;j<8;++j)
    for (int a=0;a<10;++a) for (int b=0;b<10;++b)
      state.covariance.m[i][j]+=tp[i][a]*c[a][b]*tp[j][b];
  Symmetrize(state.covariance);
  double residual[3]={g[0],g[1],g[2]}, innovation[2]={0.,0.};
  for (int a=0;a<3;++a) for (int i=0;i<10;++i) residual[a]+=h[a][i]*(x0[i]-x[i]);
  for (int a=0;a<2;++a) for (int i=0;i<3;++i) innovation[a]+=basis[a][i]*residual[i];
  state.pvConstraintChi2=innovation[0]*(inverse[0][0]*innovation[0]+inverse[0][1]*innovation[1])+
    innovation[1]*(inverse[1][0]*innovation[0]+inverse[1][1]*innovation[1]);
  state.pvConstraintNdf=2; state.totalChi2=geometricChi2+state.pvConstraintChi2;
  state.totalNdf=geometricNdf+2; state.iterations=output.iterations;
  state.status=kProductionValid; state.valid=true;
  FlightObservables flight;
  if (!(state.pvConstraintChi2>=0.) || !std::isfinite(state.totalChi2) ||
      !GetFlightObservables(state,flight) || !TransportToProductionVertex(state)) {
    output.status=kProductionUnphysicalState; return false;
  }
  output=state;
  return true;
}

bool FitProductionVertex(const ProductionVertexInput& input,
                         const ProductionVertexConfig& config,
                         double geometricChi2, int geometricNdf,
                         ProductionState& output, int motherCharge)
{
  if (motherCharge!=0) {
    output=ProductionState(); output.status=kProductionMissingMagneticField; return false;
  }
  return FitProductionVertex(input,config,geometricChi2,geometricNdf,output,0,0.);
}

bool FitProductionVertex(const ProductionVertexInput& input,
                         const ProductionVertexConfig& config,
                         double geometricChi2, int geometricNdf,
                         ProductionState& output, int motherCharge, double magneticFieldByTesla)
{
  if (!std::isfinite(magneticFieldByTesla)) {
    output=ProductionState(); output.status=kProductionInvalidConfiguration; return false;
  }
  if (motherCharge==0 || magneticFieldByTesla==0.) {
    const bool ok=FitNeutralProductionVertex(input,config,geometricChi2,geometricNdf,output,0);
    if (ok) { output.charge=motherCharge; output.magneticFieldByTesla=magneticFieldByTesla; }
    return ok;
  }
  return FitChargedProductionVertex(input,config,geometricChi2,geometricNdf,
                                    output,motherCharge,magneticFieldByTesla);
}

bool TransportToDecayVertex(ProductionState& state) { return TransportEndpoint(state,false); }
bool TransportToProductionVertex(ProductionState& state) { return TransportEndpoint(state,true); }

bool GetFlightObservables(const ProductionState& state, FlightObservables& output)
{
  output=FlightObservables();
  const double p=Norm(state.momentumGeV), e=state.energyGeV, s=state.flightParameterS;
  const double m2=e*e-p*p;
  if (!state.valid || !(p>1.e-12) || !(e>0.) || !(m2>0.) ||
      !std::isfinite(m2) || !std::isfinite(s) || !IsPositiveSemidefinite(state.covariance)) return false;
  FlightObservables result;
  result.massGeV=std::sqrt(m2); result.decayLengthCm=s*p; result.ctauCm=s*result.massGeV;
  double gm[8]={0.}, gl[8]={0.}, gt[8]={0.};
  for (int i=0; i<3; ++i) {
    gm[i+3]=-state.momentumGeV[i]/result.massGeV;
    gl[i+3]=s*state.momentumGeV[i]/p;
    gt[i+3]=s*gm[i+3];
  }
  gm[6]=e/result.massGeV; gl[7]=p;
  gt[6]=s*gm[6]; gt[7]=result.massGeV;
  if (!Error(gm,state.covariance,result.sigmaMassGeV) ||
      !Error(gl,state.covariance,result.sigmaDecayLengthCm) ||
      !Error(gt,state.covariance,result.sigmaCtauCm) ||
      !std::isfinite(result.ctauCm) || !std::isfinite(result.decayLengthCm)) return false;
  result.valid=true; output=result; return true;
}
} // namespace CeeV0Lite
