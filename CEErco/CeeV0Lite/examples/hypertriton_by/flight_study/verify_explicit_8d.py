#!/usr/bin/env python3
"""Explicit eight-parameter GLS vs profiled q7 and the formal C++ output.

Parameters: production xyz, production momentum xyz, energy, S.
Ten measured components: decay q7 and PV xyz. Full covariance is whitened.
Uses independent scipy.optimize and numerical Jacobians, not the 2D Kalman gain.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import argparse,json,subprocess
import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
from profiled_7d import fit7

def explicit_fit(y,C,charge,by):
    chol=np.linalg.cholesky(C);curvature=.00299792458*charge*by
    def model(z):
        s=z[7];ang=curvature*s
        a=s*np.sinc(ang/np.pi);b=.5*s*ang*np.sinc(ang/(2*np.pi))**2
        px,py,pz=z[3:6]
        r=z[:3]+[a*px-b*pz,s*py,a*pz+b*px]
        p=[np.cos(ang)*px-np.sin(ang)*pz,py,np.sin(ang)*px+np.cos(ang)*pz]
        return np.r_[r,p,z[6],z[:3]]
    s0=np.dot(y[:3]-y[7:],y[3:6])/np.dot(y[3:6],y[3:6]);a0=curvature*s0
    start=np.r_[y[7:],np.cos(a0)*y[3]+np.sin(a0)*y[5],y[4],-np.sin(a0)*y[3]+np.cos(a0)*y[5],y[6],s0]
    fit=least_squares(lambda z:solve_triangular(chol,model(z)-y,lower=True),start,jac='3-point',x_scale='jac',xtol=1e-13,ftol=1e-13,gtol=1e-12,max_nfev=600)
    if not fit.success:raise RuntimeError(fit.message)
    covariance=np.linalg.inv(fit.jac.T@fit.jac)
    def decay(z):return np.r_[model(z)[:7],z[7]]
    J=np.empty((8,8))
    for j in range(8):
        h=1e-5*(1+abs(fit.x[j]));zp=fit.x.copy();zm=fit.x.copy();zp[j]+=h;zm[j]-=h
        J[:,j]=(decay(zp)-decay(zm))/(2*h)
    return decay(fit.x),J@covariance@J.T,2*fit.cost,fit.x,covariance

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cee-dir',type=Path,required=True);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,default=Path('output'));args=ap.parse_args()
    out=args.out.resolve();out.mkdir(exist_ok=True,parents=True);core=args.cee_dir.resolve()/'CeeV0Lite'
    exe=core/'obj/test_charged_production'
    if not exe.exists():subprocess.run(['make','obj/test_charged_production'],cwd=core,check=True)
    fixtures=[json.loads(l) for l in subprocess.check_output([str(exe),'--reference'],text=True).splitlines()]
    d=dict(np.load(args.input));eligible=d['base']&d['topo_valid']&(d['chi']<10)
    rng=np.random.default_rng(20261014);ids=[]
    for label in [True,False]:
        ii=np.flatnonzero(eligible&(d['label']==label));ids.extend(rng.choice(ii,48,replace=False).tolist())
    # Include flight-selection boundary candidates, where numerical differences matter.
    ratio=d['L']/d['sigma_L'];ii=np.flatnonzero(eligible);ids.extend(ii[np.argsort(abs(ratio[ii]-3))[:16]].tolist());ids=sorted(set(ids))
    checks=[]
    for i in ids:
        y=np.r_[d['raw'][i,3:10],d['pv'][i]];C=np.zeros((10,10));C[:7,:7]=d['raw'][i,10:59].reshape(7,7);C[7:,7:]=np.eye(3)*.05**2
        checks.append((y,C,1,.5,d['q8'][i],d['c8'][i],float(d['chi'][i]),'toy',int(d['id'][i])))
    # Add signs, field strengths, and nonzero mother/PV cross covariances from fixtures.
    fixture_errors=[]
    for j,f in enumerate(fixtures):
        y=np.array(f['input']);C=np.array(f['covariance']).reshape(10,10)
        qd,cd,chi,qp,cp=explicit_fit(y,C,f['charge'],f['by']);scale=np.sqrt(np.diag(cp))
        fixture_errors.append({'state_sigma':float(np.max(abs(qp-np.array(f['state']))/scale)),'cov_scaled':float(np.max(abs(cp-np.array(f['state_covariance']).reshape(8,8))/np.outer(scale,scale))),'chi2':abs(chi-f['chi2'])})
        checks.append((y,C,f['charge'],f['by'],qd,cd,chi,'fixture',j))
    maximum={'state_in_sigma':0.,'covariance_scaled':0.,'mass_GeV':0.,'chi2':0.};rows=[]
    for y,C,q,by,target,tc,tchi,kind,ident in checks:
        e,ec,echi,_,_=explicit_fit(y,C,q,by)
        a,ac,achi,ok,_=fit7(y[None,:7],C[None,:7,:7],y[None,7:],C[None,7:,7:],by,q,C[None,:7,7:])
        assert ok[0],(kind,ident)
        scale=np.sqrt(np.diag(ec)[:7]);mass=lambda z:np.sqrt(z[6]**2-np.dot(z[3:6],z[3:6]))
        errs={'state_in_sigma':float(np.max(abs(a[0]-e[:7])/scale)),'covariance_scaled':float(np.max(abs(ac[0]-ec[:7,:7])/np.outer(scale,scale))),'mass_GeV':abs(mass(a[0])-mass(e)),'chi2':abs(float(achi[0])-echi)}
        if kind=='toy':
            errs['state_in_sigma']=max(errs['state_in_sigma'],float(np.max(abs(target[:7]-e[:7])/scale)))
            errs['covariance_scaled']=max(errs['covariance_scaled'],float(np.max(abs(tc[:7,:7]-ec[:7,:7])/np.outer(scale,scale))))
            errs['chi2']=max(errs['chi2'],abs(tchi-echi))
        for k in maximum:maximum[k]=max(maximum[k],errs[k])
        rows.append({'kind':kind,'id':ident,**errs})
    limits={'state_in_sigma':2e-5,'covariance_scaled':2e-5,'mass_GeV':1e-7,'chi2':1e-6}
    result={'toy_candidates':len(ids),'synthetic_correlated_and_signed_cases':len(fixtures),'total_cases':len(checks),'maximum_differences':maximum,'limits':limits,'fixture_Cpp_production_endpoint_maximum':{k:max(x[k] for x in fixture_errors) for k in fixture_errors[0]},'passed':all(maximum[k]<v for k,v in limits.items()),'comparison':'Independent explicit eight-parameter GLS, independent profiled q7, formal C++ full 8D output','rows':rows}
    fixture_limits={'state_sigma':2e-5,'cov_scaled':2e-5,'chi2':1e-6}
    result['fixture_limits']=fixture_limits
    result['passed']=result['passed'] and all(result['fixture_Cpp_production_endpoint_maximum'][k]<v for k,v in fixture_limits.items())
    (out/'explicit_8d_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2),flush=True)
    assert result['passed'],maximum
if __name__=='__main__':main()
