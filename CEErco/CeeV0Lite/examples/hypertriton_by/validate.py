#!/usr/bin/env python3
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import ctypes,json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
from charged_topology import fit,helix_matrices,endpoint,KAPPA,observables
from simulate import forward,core_fit,PV_SIGMA

HERE=Path(__file__).resolve().parent
def main():
    path=HERE/'output/candidates.npz'
    if not path.exists():path=HERE/'pilot/candidates.npz'
    d=np.load(path);ii=np.flatnonzero(d['paired']&d['label'])[:32]
    q=d['raw'][ii,3:10];c=d['raw'][ii,10:59].reshape(-1,7,7);pv=d['pv'][ii]
    pvc=np.broadcast_to(np.eye(3)*PV_SIGMA**2,(len(ii),3,3))
    r,c8,chi,ok,_=fit(q,c,pv,pvc,.5)
    assert ok.all()
    state_err=[];cov_err=[];chi_err=[]
    for k in range(len(q)):
        C=np.zeros((10,10));C[:7,:7]=c[k];C[7:,7:]=pvc[k]
        prior=np.r_[q[k],pv[k]];L=np.linalg.cholesky(C)
        def model(u):
            A,R,_=helix_matrices([u[7]],.5)
            m=np.r_[u[:7],u[:3]-A[0]@u[3:6]]
            J=np.zeros((10,8));J[:7,:7]=np.eye(7);J[7:,:3]=np.eye(3);J[7:,3:6]=-A[0];J[7:,7]=-R[0]@u[3:6]
            return np.linalg.solve(L,m-prior),np.linalg.solve(L,J)
        start=np.r_[q[k],np.dot(q[k,:3]-pv[k],q[k,3:6])/np.dot(q[k,3:6],q[k,3:6])]
        ls=least_squares(lambda u:model(u)[0],start,jac=lambda u:model(u)[1],x_scale='jac',xtol=1e-13,ftol=1e-13,gtol=1e-11,max_nfev=200)
        cr=np.linalg.inv(ls.jac.T@ls.jac)
        scale=np.sqrt(np.diag(cr));state_err.append(np.max(np.abs(ls.x-r[k])/scale))
        cov_err.append(np.max(np.abs(cr-c8[k])/np.outer(scale,scale)))
        chi_err.append(abs(2*ls.cost-chi[k]))
    prod,pcov=endpoint(r,c8,.5,True);back,bcov=endpoint(prod,pcov,.5,False)
    rtstate=float(np.max(np.abs(back-r)));scale=np.sqrt(np.diagonal(c8,axis1=1,axis2=2))
    rtcov=float(np.max(np.abs(bcov-c8)/(scale[:,:,None]*scale[:,None,:])))
    # Zero field dispatch must agree with the neutral branch. Independent neutral GLS is in tests/.
    zero,zC,zchi,zok,_=fit(q,c,pv,pvc,0)
    lib=ctypes.CDLL(str(HERE/'libtoy_core.so'));P=ctypes.POINTER(ctypes.c_double)
    lib.neutral_reference.argtypes=[P,P,P,ctypes.c_double,P];lib.neutral_reference.restype=ctypes.c_int
    ze=[];zc=[]
    for k in range(len(q)):
        qq=np.ascontiguousarray(q[k]);cc=np.ascontiguousarray(c[k]);pp=np.ascontiguousarray(pv[k]);o=np.zeros(73)
        assert lib.neutral_reference(qq.ctypes.data_as(P),cc.ctypes.data_as(P),pp.ctypes.data_as(P),PV_SIGMA**2,o.ctypes.data_as(P))==1
        ze.append(float(np.max(np.abs(o[:8]-zero[k]))));sc=np.sqrt(np.diag(zC[k]));zc.append(float(np.max(abs(o[8:72].reshape(8,8)-zC[k])/np.outer(sc,sc))))
    # Independent Lorentz-force integration, including |q|=2 and field/charge signs.
    lib.propagate.argtypes=[P,ctypes.c_int,ctypes.c_double,ctypes.c_double,P];lib.propagate.restype=ctypes.c_int
    maxode=0.;maxformula=0.
    for charge in [-2,-1,1,2]:
      for by in [-.5,.5]:
       for length in [-50.,12.,100.]:
        x=np.array([.2,-.5,1.,.35,.21,.87]);out=np.zeros(6)
        assert lib.propagate(x.ctypes.data_as(P),charge,length,by,out.ctypes.data_as(P))
        def rhs(t,u):
            p=u[3:];pn=np.linalg.norm(p)
            return np.r_[p/pn,KAPPA*charge*np.cross(p,np.array([0.,by,0.]))/pn]
        sol=solve_ivp(rhs,[0,length],x,rtol=2e-12,atol=1e-13).y[:,-1]
        maxode=max(maxode,float(np.max(np.abs(out-sol))))
        maxformula=max(maxformula,float(np.max(np.abs(out-forward(x[None,:],length,charge,by)[0]))))
    # S storage alone is an identity operation on the pre-constraint mass.
    qstore=np.column_stack((q,np.ones(len(q))));storage=np.sqrt(qstore[:,6]**2-np.sum(qstore[:,3:6]**2,axis=1))
    identity=float(np.max(np.abs(storage-d['raw'][ii,59])))
    # Confirm simulated inputs were sent to exactly the bundled core fitter.
    rr=core_fit(d['qa_x'][:8],d['qa_cov'][:8]);core_repeat=float(np.nanmax(abs(rr-d['raw'][:8])))
    results=dict(independent_GLS_cases=len(ii),max_GLS_state_difference_in_sigma=max(state_err),max_GLS_covariance_difference_scaled=max(cov_err),max_GLS_chi2_difference=max(chi_err),python_reference_endpoint_roundtrip_state_max=rtstate,python_reference_endpoint_roundtrip_covariance_scaled_max=rtcov,zero_field_vs_neutral_branch_state_max=max(ze),zero_field_vs_neutral_branch_covariance_scaled_max=max(zc),Lorentz_ODE_vs_daughter_core_max=maxode,generator_vs_daughter_core_max=maxformula,mass_change_from_storing_S_GeV=identity,repeated_core_max=core_repeat)
    assert max(state_err)<1e-4 and max(cov_err)<1e-4 and max(chi_err)<1e-7
    assert rtstate<1e-9 and rtcov<1e-8 and max(ze)<1e-5 and max(zc)<1e-4
    assert maxode<1e-8 and maxformula<1e-9 and identity<1e-10 and core_repeat<1e-12
    results['passed']=True
    (HERE/'output/validation.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))

if __name__=='__main__':main()
