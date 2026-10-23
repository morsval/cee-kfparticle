#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;OUT=HERE/'output';MH=2.99131
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':12,'axes.titlesize':14,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})
BLUE='#2168a5';ORANGE='#cf6b27';RED='#ae3446';GRAY='#bfc5cc'

def main():
    d=np.load(OUT/'candidates.npz');label=d['label'];raw=d['raw'][:,59];topo=d['mass']
    base=d['base'];pair=d['paired'];sel=d['selected'];conf=json.loads((OUT/'config.json').read_text())
    sigpair=pair&label
    def width(x):return float(np.diff(np.quantile(x,[.158655,.841345]))[0]/2*1000)
    def metrics(m,mask):
        sig=mask&label;bg=mask&~label;win=np.abs(m-MH)<.006
        ns=int((sig&win).sum());nb=int((bg&win).sum())
        return dict(signal_total=int(sig.sum()),background_total=int(bg.sum()),signal_in_6MeV_window=ns,background_in_6MeV_window=nb,S_over_B=ns/nb if nb else None,S_over_sqrt_SB=ns/np.sqrt(ns+nb),signal_mean_MeV=float(np.mean(m[sig])*1000),signal_RMS_MeV=float(np.std(m[sig])*1000),signal_sigma68_MeV=width(m[sig]))
    stats={'baseline':metrics(raw,base),'paired_before':metrics(raw,pair),'paired_after':metrics(topo,pair),'selected_after':metrics(topo,sel)}
    before=raw[sigpair];after=topo[sigpair]
    rng=np.random.default_rng(9881);est=[]
    for _ in range(500):
        ids=rng.integers(len(before),size=len(before))
        est.append(100*(1-np.std(after[ids])/np.std(before[ids])))
    stats['same_candidates_RMS_improvement_percent']=float(100*(1-np.std(after)/np.std(before)))
    stats['RMS_improvement_bootstrap_sigma_percent']=float(np.std(est))
    stats['signal_retained_fraction']=float((sel&label).sum()/(base&label).sum())
    stats['background_retained_fraction']=float((sel&~label).sum()/(base&~label).sum())
    s_p=(d['q8'][sigpair,7]-d['S_true'][sigpair])/np.sqrt(d['c8'][sigpair,7,7])
    mp=(topo[sigpair]-MH)/d['sigma_mass'][sigpair]
    stats['S_pull_mean_std']=[float(s_p.mean()),float(s_p.std())]
    stats['S_pull_core68_width']=float(np.diff(np.quantile(s_p,[.158655,.841345]))[0]/2)
    stats['S_pull_abs_gt5']=int(np.sum(np.abs(s_p)>5))
    stats['S_pull_abs_gt5_candidate_ids']=d['id'][sigpair][np.abs(s_p)>5].tolist()
    stats['mass_pull_mean_std']=[float(mp.mean()),float(mp.std())]
    stats['mean_PV_chi2_ndf2']=float(np.mean(d['chi'][sigpair]))
    bins=np.arange(2.950,3.06001,.0005);cent=(bins[1:]+bins[:-1])/2
    def draw(ax,m,mask,title,stage):
        count,_=np.histogram(m[mask],bins);sg,_=np.histogram(m[mask&label],bins);bg,_=np.histogram(m[mask&~label],bins)
        ax.stairs(bg,bins,fill=True,color=GRAY,alpha=.65,label='Combinatorial background (MC truth)')
        ax.stairs(sg,bins,color=RED,lw=1.7,label='Hypertriton signal (MC truth)')
        ax.errorbar(cent,count,yerr=np.sqrt(count),fmt='o',ms=2.8,elinewidth=.7,color=BLUE,label='Signal + background',zorder=4)
        ax.axvline(MH,color='#6a7077',lw=1,ls=':',zorder=0)
        ax.set(xlim=(bins[0],bins[-1]),ylim=(0,max(count)*1.45),xlabel=r'$M(^{3}{\rm He}+\pi^-)$  [GeV/$c^2$]',ylabel=r'Candidates / (0.5 MeV/$c^2$)',title=title)
        ax.grid(axis='y',alpha=.16);ax.legend(loc='upper right',frameon=False,fontsize=8.6)
        t=stats[stage]
        lines=[f"Injected mass: {MH:.5f} GeV/$c^2$", f"Signal: {t['signal_total']:,}; background: {t['background_total']:,}", f"Within $\\pm$6 MeV: S/B = {t['S_over_B']:.2f}"]
        for row,txt in enumerate(lines):
            ax.text(.035,.955-row*.043,txt,transform=ax.transAxes,va='top',fontsize=9.2)
        return count,sg,bg
    descriptions=[('mass_before',raw,base,'7D secondary-vertex fit','baseline','Common cuts: daughter vertex $\\chi^2<10$, DCA $<1$ cm; same TPC / TOF PID'),
                  ('mass_after_S',topo,sel,'8D charged-mother PV fit + flight selection','selected_after',r'Additional cuts: $\chi^2_{\rm PV}/2<5$, $S/\sigma_S>3$; no mass constraint')]
    histcols=[cent];headers=['mass_center_GeV']
    for stem,m,mask,title,stage,desc in descriptions:
        fig,ax=plt.subplots(figsize=(10,6.4));fig.subplots_adjust(left=.11,right=.97,top=.81,bottom=.17)
        fig.suptitle(r'$^{3}_{\Lambda}{\rm H}\rightarrow{}^{3}{\rm He}+\pi^-$  |  Toy MC',x=.11,ha='left',y=.98,fontsize=19,fontweight='bold')
        fig.text(.11,.909,r'STAR 3 GeV spectrum shape; CEE field $B_y=+0.5$ T; approximate detector response',fontsize=10,color='#4b5159')
        h=draw(ax,m,mask,title,stage);histcols.extend(h);headers.extend([stem+'_all',stem+'_signal',stem+'_background'])
        fig.text(.11,.075,desc,fontsize=10)
        fig.text(.11,.032,'Simulation only. Chosen signal/background mixture; no STAR acceptance or efficiency normalization.',fontsize=9,color='#6a7077')
        fig.savefig(OUT/(stem+'.png'));fig.savefig(OUT/(stem+'.svg'));plt.close(fig)
    # Same-event comparison proves any width change comes from the fit, not selection.
    fig,axs=plt.subplots(1,2,figsize=(12,4.7));fig.subplots_adjust(left=.08,right=.98,bottom=.18,top=.82,wspace=.30)
    rb=np.arange(-12,12.001,.4)
    axs[0].hist((before-MH)*1000,bins=rb,histtype='step',lw=1.8,color=BLUE,label=f"7D: RMS {np.std(before)*1000:.3f} MeV")
    axs[0].hist((after-MH)*1000,bins=rb,histtype='step',lw=1.8,color=ORANGE,label=f"8D + PV: RMS {np.std(after)*1000:.3f} MeV")
    axs[0].set(xlabel=r'$M-M_{\rm injected}$ [MeV/$c^2$]',ylabel='Truth-matched signal / 0.4 MeV',title=f'Identical {len(before):,} signal candidates')
    axs[0].legend(frameon=False,fontsize=9)
    axs[1].hist(s_p,bins=np.arange(-6,28.01,.4),density=True,histtype='step',lw=1.8,color=ORANGE,label=f'Full mean {s_p.mean():.3f}, std {s_p.std():.3f}')
    gx=np.linspace(-6,28,500);axs[1].plot(gx,np.exp(-gx*gx/2)/np.sqrt(2*np.pi),ls='--',color='#737980',label='Unit Gaussian')
    axs[1].set(yscale='log',ylim=(1e-5,1),xlim=(-6,28))
    axs[1].text(.97,.69,f"Core 68% width: {stats['S_pull_core68_width']:.3f}\n{stats['S_pull_abs_gt5']} / {len(s_p)} outside $\\pm$5",transform=axs[1].transAxes,ha='right',fontsize=9)
    axs[1].set(xlabel=r'$(S_{\rm fit}-S_{\rm true})/\sigma_S$',ylabel='Probability density',title='Flight-parameter uncertainty check')
    axs[1].legend(frameon=False,fontsize=9)
    fig.suptitle('Refit effect and uncertainty calibration | Toy MC',x=.08,ha='left',fontweight='bold')
    fig.text(.08,.03,'No flight-significance or PV-chi2 cut in these checks. Storing S alone changes the mass by exactly zero.',fontsize=10)
    fig.savefig(OUT/'same_candidates_and_S_pull.png');plt.close(fig)
    np.savetxt(OUT/'mass_histograms.csv',np.column_stack(histcols),delimiter=',',header=','.join(headers),comments='')
    np.savetxt(OUT/'candidate_summary.csv',np.column_stack([d['id'],label.astype(int),base.astype(int),pair.astype(int),sel.astype(int),raw,topo,d['q8'][:,7],np.sqrt(d['c8'][:,7,7]),d['chi'],d['ct'],d['sigma_ct']]),delimiter=',',header='id,is_signal,base,paired,selected,mass_SV_GeV,mass_PV_GeV,S_cm_per_GeV,sigma_S,chi2_PV_ndf2,ct_cm,sigma_ct_cm',comments='')
    (OUT/'summary.json').write_text(json.dumps(stats,indent=2));print(json.dumps(stats,indent=2))

if __name__=='__main__':main()
