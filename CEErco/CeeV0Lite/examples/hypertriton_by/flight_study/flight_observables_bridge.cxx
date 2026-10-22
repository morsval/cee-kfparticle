#include "CeeV0LiteProductionVertex.h"
#include <limits>
using namespace CeeV0Lite;
extern "C" void flight_observables_batch(int n,const double* q,const double* c,double* out){
  for(int k=0;k<n;++k){
    for(int i=0;i<6;++i)out[k*6+i]=std::numeric_limits<double>::quiet_NaN();
    ProductionState s;s.valid=true;s.atProductionVertex=false;s.charge=1;s.magneticFieldByTesla=.5;
    for(int i=0;i<3;++i){s.positionCm[i]=q[k*8+i];s.momentumGeV[i]=q[k*8+3+i];}
    s.energyGeV=q[k*8+6];s.flightParameterS=q[k*8+7];
    for(int i=0;i<8;++i)for(int j=0;j<8;++j)s.covariance.m[i][j]=c[k*64+i*8+j];
    FlightObservables f;if(!GetFlightObservables(s,f))continue;
    out[k*6]=f.massGeV;out[k*6+1]=f.sigmaMassGeV;out[k*6+2]=f.decayLengthCm;
    out[k*6+3]=f.sigmaDecayLengthCm;out[k*6+4]=f.ctauCm;out[k*6+5]=f.sigmaCtauCm;
  }
}
