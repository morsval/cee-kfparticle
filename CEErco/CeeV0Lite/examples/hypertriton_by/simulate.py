#!/usr/bin/env python3
"""Hypertriton integration toy using the formal C++ v0.6.5 core in CEE By."""
import os
os.environ.setdefault('OMP_NUM_THREADS','4')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse, ctypes, json, time
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from charged_topology import KAPPA,fit,observables

HERE=Path(__file__).resolve().parent
MH=2.99131; MHE=2.80839160743; MPI=.13957039; C_LIGHT=29.9792458
BY=.5; PV_SIGMA=.05; POS_SIGMA=.08; ANGLE_SIGMA=.002
TOF_SIGMA=.100; T0_SIGMA=.030; PATH_SIGMA=.20

def forward(x,length,charge,by=BY):
    x=np.asarray(x);p=x[:,3:6];pn=np.linalg.norm(p,axis=1)
    s=np.asarray(length)/pn;k=KAPPA*charge*by;a=k*s
    si=s*np.sinc(a/np.pi);om=s*(a/2)*np.sinc(a/(2*np.pi))**2
    out=x.copy()
    out[:,0]+=si*p[:,0]-om*p[:,2];out[:,1]+=s*p[:,1];out[:,2]+=si*p[:,2]+om*p[:,0]
    out[:,3]=np.cos(a)*p[:,0]-np.sin(a)*p[:,2];out[:,5]=np.cos(a)*p[:,2]+np.sin(a)*p[:,0]
    return out

def parent_sample(n,rng):
    T=json.loads((HERE/'data/spectrum_fit.json').read_text())['T_GeV']
    grid=np.linspace(.3,4.,10000);density=grid*np.exp(-(np.sqrt(MH**2+grid**2)-MH)/T)
    cdf=cumulative_trapezoid(density,grid,initial=0);cdf/=cdf[-1]
    pt=np.interp(rng.random(n),cdf,grid);yc=rng.uniform(-.25,0,n)
    yboost=np.arccosh(3./(2*.9382720813));yl=yc+yboost
    phi=rng.uniform(0,2*np.pi,n);p=np.column_stack((pt*np.cos(phi),pt*np.sin(phi),np.sqrt(MH**2+pt**2)*np.sinh(yl)))
    ct=rng.exponential(.223*C_LIGHT,n);s=ct/MH;pn=np.linalg.norm(p,axis=1)
    mother0=np.column_stack((np.zeros((n,3)),p));mother=forward(mother0,s*pn,1)
    q=np.sqrt((MH**2-(MHE+MPI)**2)*(MH**2-(MHE-MPI)**2))/(2*MH)
    z=rng.uniform(-1,1,n);az=rng.uniform(0,2*np.pi,n)
    rest=q*np.column_stack((np.sqrt(1-z*z)*np.cos(az),np.sqrt(1-z*z)*np.sin(az),z))
    P=mother[:,3:6];E=np.sqrt(MH*MH+pn*pn);gamma=E/MH;beta=P/E[:,None]
    ds=[]
    for r,m in [(rest,MHE),(-rest,MPI)]:
        dot=np.sum(r*beta,axis=1)
        f=gamma**2/(gamma+1)*dot+gamma*np.sqrt(m*m+q*q)
        ds.append(np.column_stack((mother[:,:3],r+f[:,None]*beta)))
    decaytime=s*E/C_LIGHT
    return np.stack(ds,axis=1),s,decaytime,pt,yc,mother

def smear(truth,rng):
    n=len(truth);reco=np.empty_like(truth);cov=np.zeros((n,2,6,6));relall=[]
    for d,m,z in [(0,MHE,2),(1,MPI,-1)]:
        reference=forward(truth[:,d],np.full(n,20.),z)
        p=reference[:,3:6];pn=np.linalg.norm(p,axis=1);u=p/pn[:,None]
        beta=pn/np.sqrt(pn*pn+m*m);rig=pn/abs(z)
        rel=np.sqrt(.020**2+(.005*rig)**2+(.006/(beta*rig))**2)
        radial=rng.normal(size=n)*pn*rel
        tangent=rng.normal(size=(n,3));tangent-=np.sum(tangent*u,axis=1)[:,None]*u
        reco[:,d,:3]=reference[:,:3]+rng.normal(0,POS_SIGMA,(n,3))
        reco[:,d,3:6]=p+radial[:,None]*u+(pn*ANGLE_SIGMA)[:,None]*tangent
        uu=u[:,:,None]*u[:,None,:]
        cov[:,d,:3,:3]=np.eye(3)*POS_SIGMA**2
        cov[:,d,3:,3:]=(pn*rel)[:,None,None]**2*uu+(pn*ANGLE_SIGMA)[:,None,None]**2*(np.eye(3)-uu)
        relall.append(rel)
    return reco,cov,np.array(relall).T

def core_fit(x,c):
    lib=ctypes.CDLL(str(HERE/'libtoy_core.so'));P=ctypes.POINTER(ctypes.c_double)
    lib.fit_batch.argtypes=[ctypes.c_int,P,P,ctypes.c_double,P]
    x=np.ascontiguousarray(x);c=np.ascontiguousarray(c);out=np.empty((len(x),61))
    t=time.time()
    for a in range(0,len(x),4000):
        b=min(a+4000,len(x));xx=x[a:b];cc=c[a:b];oo=out[a:b]
        lib.fit_batch(b-a,xx.ctypes.data_as(P),cc.ctypes.data_as(P),BY,oo.ctypes.data_as(P))
        print(f'core reconstruction {b}/{len(x)}; {time.time()-t:.1f}s',flush=True)
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--signal',type=int,default=6000);ap.add_argument('--background',type=int,default=120000);ap.add_argument('--seed',type=int,default=20260929);ap.add_argument('--out',default='output');ap.add_argument('--track-reference',choices=['pv-dca','stored'],default='pv-dca');args=ap.parse_args()
    outdir=HERE/args.out;outdir.mkdir(exist_ok=True)
    rng=np.random.default_rng(args.seed);ns=args.signal;nb=args.background
    signal,s,td,pt,yc,mother=parent_sample(ns,rng)
    # Independent one-daughter marginals; there is no injected background mass PDF.
    ba,_,ta,_,_,_=parent_sample(nb,rng);bb,_,tb,_,_,_=parent_sample(nb,rng)
    background=np.stack((ba[:,0],bb[:,1]),axis=1)
    # Fixed, declared mixture: 80% prompt-prompt, 20% independent displaced tracks.
    prompt=rng.random(nb)<.8;background[prompt,:,:3]=0;ta[prompt]=0;tb[prompt]=0
    truth=np.concatenate((signal,background));labels=np.arange(ns+nb)<ns
    reco,cov,rel=smear(truth,rng);pv=rng.normal(0,PV_SIGMA,(ns+nb,3))
    # STAR-style tracking momentum; TOF is a PID consistency test, not p from TOF.
    # A difference of daughter times eliminates common t0 and the unknown decay time.
    ptrue=np.linalg.norm(truth[:,:,3:6],axis=2);pr=np.linalg.norm(reco[:,:,3:6],axis=2)
    masses=np.array([MHE,MPI]);truepath=rng.uniform(220,280,(ns+nb,2));path=truepath+rng.normal(0,PATH_SIGMA,(ns+nb,2))
    times=np.column_stack((np.r_[td,ta],np.r_[td,tb]))+truepath*np.sqrt(1+(masses/ptrue)**2)/C_LIGHT
    times+=rng.normal(0,T0_SIGMA,(ns+nb,1))+rng.normal(0,np.sqrt(TOF_SIGMA**2-T0_SIGMA**2),(ns+nb,2))
    expect=path*np.sqrt(1+(masses/pr)**2)/C_LIGHT
    dt=(times[:,0]-times[:,1])-(expect[:,0]-expect[:,1])
    dtdp=-path*masses**2/(C_LIGHT*pr**2*np.sqrt(pr**2+masses**2))
    sigdt=np.sqrt(2*(TOF_SIGMA**2-T0_SIGMA**2)+np.sum((dtdp*pr*rel)**2+(PATH_SIGMA*np.sqrt(1+(masses/pr)**2)/C_LIGHT)**2,axis=1))
    pid=np.abs(dt)<3*sigdt
    accept=(pr[:,1]>.10)&(pr[:,1]<2.)&(pr[:,0]>.3)&(pr[:,0]<10.)&pid
    idx=np.flatnonzero(accept);n=len(idx)
    fit_x=reco[idx];fit_c=cov[idx]
    reference_diag=np.full((n,2,4),np.nan);reference_diag[:,:,0]=0
    if args.track_reference=='pv-dca':
        from pv_dca_reference import prepare_reference
        fit_x,fit_c,reference_diag=prepare_reference(fit_x,fit_c,pv[idx],BY)
    reference_ok=np.all(reference_diag[:,:,0]<=1,axis=1)
    raw=np.full((n,61),np.nan);raw[:,0]=-300  # rejected input reference
    raw[reference_ok]=core_fit(fit_x[reference_ok],fit_c[reference_ok])
    q=raw[:,3:10];c7=raw[:,10:59].reshape(-1,7,7)
    base=(raw[:,0]==1)&(raw[:,1]<10)&(raw[:,2]<1.)
    fitidx=np.flatnonzero(base)
    q8=np.full((n,8),np.nan);c8=np.full((n,8,8),np.nan);chi=np.full(n,np.nan);valid=np.zeros(n,bool)
    for a in range(0,len(fitidx),5000):
        ii=fitidx[a:a+5000];r=fit(q[ii],c7[ii],pv[idx[ii]],np.broadcast_to(np.eye(3)*PV_SIGMA**2,(len(ii),3,3)),BY)
        q8[ii],c8[ii],chi[ii],valid[ii],_=r
    mo,sm,L,sL,ct,sct=observables(q8,c8)
    sigS=np.sqrt(np.maximum(c8[:,7,7],0));sel=base&valid&(chi<10)&(q8[:,7]>3*sigS)
    # Same-ID sample isolates refitting from candidate rejection.
    paired=base&valid
    strue=np.r_[s,np.full(nb,np.nan)][idx]
    np.savez_compressed(outdir/'candidates.npz',id=idx,label=labels[idx],raw=raw,q8=q8,c8=c8,pv=pv[idx],base=base,paired=paired,selected=sel,topo_valid=valid,chi=chi,mass=mo,sigma_mass=sm,L=L,sigma_L=sL,ct=ct,sigma_ct=sct,S_true=strue,signal_pt=pt,signal_y=yc,signal_mother=mother,
                        # A subset of measured inputs allows exact numerical rechecks.
                        qa_x=fit_x[:200],qa_cov=fit_c[:200],reference_diagnostics=reference_diag)
    counters=dict(generated_signal=ns,generated_background=nb,after_acceptance_signal=int(labels[idx].sum()),after_acceptance_background=int((~labels[idx]).sum()),fit_status={str(int(k)):int(v) for k,v in zip(*np.unique(raw[:,0],return_counts=True))})
    for name,mask in [('base',base),('paired',paired),('selected',sel)]:counters[name]=dict(signal=int((mask&labels[idx]).sum()),background=int((mask&~labels[idx]).sum()))
    config=dict(track_reference=args.track_reference,reference_transport="fixed nominal arc; full 6D C; no PV constraint; no process noise",reference_search_half_range_cm=1000.,reference_tolerance_cm=1e-8,reference_max_iterations=64,reference_status_counts={str(int(k)):int(v) for k,v in zip(*np.unique(reference_diag[:,:,0],return_counts=True))},seed=args.seed,mass_GeV=MH,lifetime_ps=223,field_By_T=BY,mother_charge=1,daughter_charges=[2,-1],daughter_masses_GeV=masses.tolist(),y_cm_range=[-.25,0],pt_GeV_range=[.3,4],y_boost=float(np.arccosh(3/(2*.9382720813))),pv_sigma_cm=PV_SIGMA,reference_position_sigma_cm=POS_SIGMA,direction_sigma_rad=ANGLE_SIGMA,smearing_reference_arc_from_true_SV_cm=20,momentum_resolution='sqrt(0.020^2+(0.005*R)^2+(0.006/(beta*R))^2), R=|p|/|Z| in GeV/c',tof_total_sigma_ns=TOF_SIGMA,common_t0_sigma_ns=T0_SIGMA,tof_path_sigma_cm=PATH_SIGMA,tof_path_range_cm=[220,280],background_prompt_fraction=.8,counts=counters)
    (outdir/'config.json').write_text(json.dumps(config,indent=2));print(json.dumps(counters,indent=2),flush=True)

if __name__=='__main__':main()
