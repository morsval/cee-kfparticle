"""Same PV and L/sigmaL cuts, independent q7-derived observable vs formal q8.

Input is the existing v0.6.4 toy candidates.npz. No changed reconstruction,
extra mass/lifetime constraint, or copied q8 flight output on the q7 branch.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from equivalent_7d import fit7

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,default=Path('results'));args=ap.parse_args()
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    d=dict(np.load(args.input));raw=d['raw'];base=d['base'];n=len(base)
    q7=np.full((n,7),np.nan);c7=np.full((n,7,7),np.nan);chi7=np.full(n,np.nan);ok7=np.zeros(n,bool)
    l7=np.full(n,np.nan);dl7=np.full(n,np.nan);ii=np.flatnonzero(base)
    for begin in range(0,len(ii),4000):
        jj=ii[begin:begin+4000];pvc=np.broadcast_to(np.eye(3)*.05**2,(len(jj),3,3))
        q7[jj],c7[jj],chi7[jj],ok7[jj],_,l7[jj],dl7[jj]=fit7(raw[jj,3:10],raw[jj,10:59].reshape(-1,7,7),d['pv'][jj],pvc,return_flight=True)
        print('Independent q7 with derived L/sigmaL:',min(begin+4000,len(ii)),len(ii),flush=True)
    # The right branch is the original formal C++ q8 fit. Propagate its complete
    # C8 with the exact same scalar definition used by GetFlightObservables.
    q8=d['q8'];c8=d['c8'];p8=np.linalg.norm(q8[:,3:6],axis=1);g=np.zeros((n,8))
    g[:,3:6]=q8[:,7,None]*q8[:,3:6]/p8[:,None];g[:,7]=p8
    l8=q8[:,7]*p8;dl8=np.sqrt(np.einsum('ni,nij,nj->n',g,c8,g));ok8=d['topo_valid'];chi8=d['chi']
    m7=np.sqrt(q7[:,6]**2-np.sum(q7[:,3:6]**2,axis=1));m8=np.sqrt(q8[:,6]**2-np.sum(q8[:,3:6]**2,axis=1))
    ratio7=l7/dl7;ratio8=l8/dl8
    pv7=base&ok7&(chi7<10);pv8=base&ok8&(chi8<10)
    left=pv7&(dl7>0)&(ratio7>3);right=pv8&(dl8>0)&(ratio8>3)
    both=base&ok7&ok8
    scale=np.sqrt(np.diagonal(c8[both,:7,:7],axis1=1,axis2=2))
    diff={'validity_mismatches':int(np.sum((base&ok7)!=(base&ok8))),
      'PV_cut_mismatches':int(np.sum(pv7!=pv8)),'final_cut_mismatches':int(np.sum(left!=right)),
      'max_L_absolute_cm':float(np.max(abs(l7[both]-l8[both]))),
      'max_sigmaL_relative':float(np.max(abs(dl7[both]-dl8[both])/dl8[both])),
      'max_L_over_sigmaL_absolute':float(np.max(abs(ratio7[both]-ratio8[both]))),
      'max_mass_absolute_GeV':float(np.max(abs(m7[both]-m8[both]))),
      'max_state_in_sigma':float(np.max(abs(q7[both]-q8[both,:7])/scale)),
      'max_covariance_scaled':float(np.max(abs(c7[both]-c8[both,:7,:7])/(scale[:,:,None]*scale[:,None,:]))),
      'max_chi2_absolute':float(np.max(abs(chi7[both]-chi8[both])))}
    label=d['label'];bins=np.round(np.arange(2.95,3.06001,.0005),7)
    hist={'edges':bins.tolist()};stats={}
    for key,m,mask in [('left',m7,left),('right',m8,right)]:
        win=abs(m-2.99131)<.006;ns=int(np.sum(mask&win&label));nb=int(np.sum(mask&win&~label))
        stats[key]={'signal_total':int(np.sum(mask&label)),'background_total':int(np.sum(mask&~label)),
          'signal_window':ns,'background_window':nb,'S_over_B':ns/nb,'S_over_sqrt_SB':ns/np.sqrt(ns+nb)}
        hist[key]={k:np.histogram(m[mask&s],bins)[0].tolist() for k,s in [('all',np.ones(n,bool)),('signal',label),('background',~label)]}
    same_hist=all(hist['left'][k]==hist['right'][k] for k in ['all','signal','background'])
    passed=all(diff[k]==0 for k in ['validity_mismatches','PV_cut_mismatches','final_cut_mismatches']) and diff['max_sigmaL_relative']<1e-8 and diff['max_L_over_sigmaL_absolute']<1e-8 and same_hist
    result={'method_left':'Independent q7 state; S locally eliminated; L and sigmaL derived directly from original input covariance and fitted Jacobians; no 8x8 state/covariance',
      'method_right':'Formal C++ q8 state and complete C8 flight propagation','shared_cuts':'PV fit valid, chi2PV/2<5, L/sigmaL>3',
      'base_candidates':int(base.sum()),'both_valid':int(both.sum()),'mass_window_GeV':[2.98531,2.99731],
      'same_histogram_bins':same_hist,'passed':passed,'differences':diff,**stats,
      'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest()}
    (out/'equivalent_selection_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'histograms.json').write_text(json.dumps(hist)+'\n')
    jj=np.flatnonzero(base)
    cols=np.column_stack([d['id'][jj],label[jj],ok7[jj],ok8[jj],left[jj],right[jj],m7[jj],m8[jj],l7[jj],l8[jj],dl7[jj],dl8[jj],ratio7[jj],ratio8[jj],chi7[jj],chi8[jj]])
    np.savetxt(out/'equivalent_candidates.csv',cols,delimiter=',',fmt='%.15g',header='id,is_signal,valid7,valid8,selected7,selected8,mass7_GeV,mass8_GeV,L7_cm,L8_cm,sigmaL7_cm,sigmaL8_cm,ratio7,ratio8,chi2PV7,chi2PV8',comments='')
    print(json.dumps(result,indent=2),flush=True)
    assert passed,'Do not report agreement: equivalent selection validation failed'

if __name__=='__main__':main()
