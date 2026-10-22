#!/usr/bin/env python3
"""Controlled PV+chi2 selection vs an additional L/sigmaL cut.
Usage: python run_flight_study.py --cee-dir /path/CEErco --input candidates.npz --out output
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','4')
from pathlib import Path
import argparse,ctypes,json,subprocess,sys,hashlib
import numpy as np
from profiled_7d import fit7

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cee-dir',type=Path,required=True);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,default=Path('output'));a=ap.parse_args()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=True);here=Path(__file__).resolve().parent
    core=a.cee_dir.resolve()/'CeeV0Lite';d=np.load(a.input);label=d['label'];raw=d['raw'];q8=d['q8'];c8=d['c8'];n=len(label)
    build=['g++','-O3','-std=c++11','-fPIC','-shared','-I'+str(core/'include'),str(here/'flight_observables_bridge.cxx')]
    build += [str(core/'src'/('CeeV0Lite'+x+'.cxx')) for x in ['Math','Types','Transport','ProductionVertex']]
    build += ['-o',str(out/'libflight_observables.so')];subprocess.run(build,check=True)
    lib=ctypes.CDLL(str(out/'libflight_observables.so'));P=ctypes.POINTER(ctypes.c_double)
    lib.flight_observables_batch.argtypes=[ctypes.c_int,P,P,P]
    qq=np.ascontiguousarray(q8);cc=np.ascontiguousarray(c8);obs=np.empty((n,6))
    lib.flight_observables_batch(n,qq.ctypes.data_as(P),cc.ctypes.data_as(P),obs.ctypes.data_as(P))
    mass,sm,L,sL,ct,sct=obs.T
    base=d['base'];valid=d['topo_valid']&np.isfinite(obs).all(axis=1)&(sL>0)
    left=base&valid&(d['chi']<10);right=left&(L/sL>3)
    old=left&(q8[:,7]>3*np.sqrt(c8[:,7,7]))
    p=np.linalg.norm(q8[:,3:6],axis=1);u=q8[:,3:6]/p[:,None];S=q8[:,7]
    varp=np.einsum('ni,nij,nj->n',u,c8[:,3:6,3:6],u)
    covSp=np.einsum('ni,ni->n',c8[:,7,3:6],u)
    vL=p*p*c8[:,7,7]+S*S*varp+2*S*p*covSp
    check=base&valid
    rel=float(np.max(abs(vL[check]-sL[check]**2)/sL[check]**2))
    assert rel<1e-9,rel
    # Independent seven-dimensional output fit, same measured prior and PV.
    q7=np.full((n,7),np.nan);c7=np.full((n,7,7),np.nan);chi7=np.full(n,np.nan);ok7=np.zeros(n,bool)
    ii=np.flatnonzero(base)
    for begin in range(0,len(ii),4000):
        jj=ii[begin:begin+4000];pvc=np.broadcast_to(np.eye(3)*.05**2,(len(jj),3,3))
        q7[jj],c7[jj],chi7[jj],ok7[jj],_=fit7(raw[jj,3:10],raw[jj,10:59].reshape(-1,7,7),d['pv'][jj],pvc)
        print('profiled7D',min(begin+4000,len(ii)),len(ii),flush=True)
    both=base&ok7&d['topo_valid'];scale=np.sqrt(np.diagonal(c8[both,:7,:7],axis1=1,axis2=2))
    mass7=np.sqrt(np.maximum(q7[:,6]**2-np.sum(q7[:,3:6]**2,axis=1),0))
    cmp={'inputs':int(base.sum()),'both_valid':int(both.sum()),'only_7d_valid':int((base&ok7&~d['topo_valid']).sum()),'only_8d_valid':int((base&~ok7&d['topo_valid']).sum()),'both_invalid':int((base&~ok7&~d['topo_valid']).sum()),'max_state_difference_in_sigma':float(np.max(abs(q7[both]-q8[both,:7])/scale)),'max_covariance_difference_scaled':float(np.max(abs(c7[both]-c8[both,:7,:7])/(scale[:,:,None]*scale[:,None,:]))),'max_mass_difference_GeV':float(np.max(abs(mass7[both]-mass[both]))),'max_chi2_difference':float(np.max(abs(chi7[both]-d['chi'][both]))),'pv_chi2_selection_changed':int(np.sum((ok7&(chi7<10)&base)!=(d['topo_valid']&(d['chi']<10)&base))),'C7_definition':'Marginal covariance, no conditioning on S=0','method_7d':'Independent NumPy weighted rank-one elimination, q7 output only','method_8d':'Formal C++ transverse constraints and full 8D output, same decay endpoint'}
    cmp['passed']=bool(cmp['only_7d_valid']==0 and cmp['only_8d_valid']==0 and cmp['pv_chi2_selection_changed']==0 and cmp['max_state_difference_in_sigma']<2e-5 and cmp['max_covariance_difference_scaled']<2e-5 and cmp['max_chi2_difference']<1e-5)
    (out/'comparison_7d_8d.json').write_text(json.dumps(cmp,indent=2)+'\n')
    print(json.dumps(cmp,indent=2),flush=True)
    assert cmp['passed'],'Independent comparison failed; do not claim equivalence'
    def metrics(mask):
        sg=mask&label;bg=mask&~label;window=abs(mass-2.99131)<.006
        ns=int((sg&window).sum());nb=int((bg&window).sum())
        return {'signal_total':int(sg.sum()),'background_total':int(bg.sum()),'window_signal':ns,'window_background':nb,'window_S_over_B':ns/nb,'window_significance':ns/np.sqrt(ns+nb),'signal_mass_RMS_MeV':float(np.std(mass[sg])*1000)}
    stats={'core_version':'0.6.4','field_By_T':.5,'input_pairs':n,'base_pairs':int(base.sum()),'common_PV_chi2_over_2_max':5,'right_L_over_sigmaL_min':3,'sigmaL_formula':'p^2 Var(S)+S^2 Var(p)+2 S p Cov(S,p)','sigmaL_formula_max_relative_variance_difference':rel,'left':metrics(left),'right':metrics(right),'old_S_selection':metrics(old),'cut_switch_changed_candidates':int((right!=old).sum()),'cut_switch_changed_signal':int(((right!=old)&label).sum()),'cut_switch_changed_background':int(((right!=old)&~label).sum()),'retained_signal_fraction':float((right&label).sum()/(left&label).sum()),'retained_background_fraction':float((right&~label).sum()/(left&~label).sum()),'retention_denominator':'PV refit valid and chi2PV/2 < 5 (left panel), no flight selection','mass_window_GeV':[2.98531,2.99731],'candidate_source_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest()}
    bins=np.arange(2.950,3.06001,.0005);bins=np.round(bins,7)
    hist={'edges':bins.tolist(),'centers':((bins[1:]+bins[:-1])/2).tolist()}
    for name,mask in [('left',left),('right',right)]:hist[name]={k:np.histogram(mass[mask&m],bins)[0].tolist() for k,m in [('all',np.ones(n,bool)),('signal',label),('background',~label)]}
    hist['max_y']=int(np.ceil(max(max(hist['left']['all']),max(hist['right']['all']))/50)*50+50)
    (out/'histograms.json').write_text(json.dumps(hist)+'\n');(out/'summary.json').write_text(json.dumps(stats,indent=2)+'\n')
    jj=np.flatnonzero(base)
    cols=np.column_stack([d['id'][jj],label[jj].astype(int),d['topo_valid'][jj].astype(int),left[jj].astype(int),right[jj].astype(int),old[jj].astype(int),mass[jj],L[jj],sL[jj],L[jj]/sL[jj],S[jj],np.sqrt(c8[jj,7,7]),d['chi'][jj],chi7[jj],mass7[jj]])
    np.savetxt(out/'candidate_comparison.csv',cols,delimiter=',',fmt='%.15g',header='id,is_signal,pv_valid,left_PV_chi2,right_LsigmaL,old_SsigmaS,mass_GeV,L_cm,sigmaL_cm,L_over_sigmaL,S,sigmaS,chi2PV_8d,chi2PV_7d,mass7d_GeV',comments='')
    print(json.dumps(stats,indent=2),flush=True)
    # Standalone scientific export; native PPT charts use the same bins/data.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),sharex=True,sharey=True)
    for ax,name,title in zip(axes,['left','right'],['PV refit + chi2PV / 2 < 5','Same selection + L / sigmaL > 3']):
        for k,color in [('all','#243746'),('signal','#B33A3A'),('background','#2462A6')]:ax.stairs(hist[name][k],bins,label=k,color=color,lw=1.3)
        ax.set(title=title,xlim=(2.95,3.06),ylim=(0,hist['max_y']),xlabel=r'$M(^{3}{\rm He}+\pi^-)$ [GeV/$c^2$]');ax.legend(frameon=False)
    axes[0].set_ylabel('Candidates / (0.5 MeV/c^2)');fig.tight_layout();fig.savefig(out/'controlled_mass_spectra.pdf');fig.savefig(out/'controlled_mass_spectra.png',dpi=180);plt.close(fig)
if __name__=='__main__':main()
