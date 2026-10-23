"""Old scientific plot style, equal 7D/8D selections, explicit mass window."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,default=Path('results'));a=ap.parse_args();out=a.results
    hist=json.loads((out/'histograms.json').read_text());stats=json.loads((out/'equivalent_selection_summary.json').read_text());assert stats['passed']
    edges=np.array(hist['edges']);centers=(edges[:-1]+edges[1:])/2;lo,hi=stats['mass_window_GeV'];M0=2.99131
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.labelsize':13,'axes.titlesize':15,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'pdf.fonttype':42})
    blue='#2168a5';red='#ae3446';gray='#cdd2d8';orange='#bf7226'
    def draw(ax,key):
        panel=hist[key];bg=np.array(panel['background']);sg=np.array(panel['signal']);all_=np.array(panel['all'])
        ax.axvspan(lo,hi,facecolor='#fff3d9',alpha=.9,zorder=0)
        ax.stairs(bg,edges,fill=True,color=gray,alpha=.95,label='Background (MC truth)',zorder=1)
        ax.stairs(sg,edges,color=red,lw=1.7,label='Signal (MC truth)',zorder=3)
        ax.errorbar(centers,all_,yerr=np.sqrt(all_),fmt='o',ms=2.6,elinewidth=.65,color=blue,label='Signal + background',zorder=4)
        for x in [lo,hi]:ax.axvline(x,color=orange,lw=1.25,ls=(0,(5,3)),zorder=5)
        ax.axvline(M0,color='#7d848d',lw=.8,ls=':',zorder=2)
        ax.annotate('',xy=(lo,560),xytext=(hi,560),arrowprops={'arrowstyle':'|-|','color':orange,'lw':1.1,'mutation_scale':5})
        ax.text(M0,577,r'$M_0\pm6$ MeV/$c^2$',ha='center',va='bottom',fontsize=11,color=orange)
        title=r'7D state + derived $L/\sigma_L$' if key=='left' else r'8D state + stored covariance'
        ax.set(xlim=(2.95,3.06),ylim=(0,650),xlabel=r'$M(^{3}\mathrm{He}+\pi^-)$ [GeV/$c^2$]',ylabel=r'Candidates / (0.5 MeV/$c^2$)',title=title)
        ax.set_xticks([2.96,2.98,3.00,3.02,3.04,3.06]);ax.set_yticks(np.arange(0,651,100))
        ax.grid(axis='y',color='#e7e9ed',lw=.7,zorder=-1)
        ax.legend(loc='upper right',frameon=False,fontsize=9.2,handlelength=1.5,labelspacing=.45,borderaxespad=.6)
    for key in ['left','right']:
        fig,ax=plt.subplots(figsize=(6.6,4.6));fig.subplots_adjust(left=.14,right=.95,bottom=.17,top=.9)
        draw(ax,key);fig.savefig(out/f'mass_equivalent_{key}.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(14,6.4));fig.subplots_adjust(left=.075,right=.965,bottom=.24,top=.77,wspace=.25)
    for ax,key in zip(axes,['left','right']):draw(ax,key)
    fig.text(.075,.95,r'$^{3}_{\Lambda}\mathrm{H}\rightarrow{}^{3}\mathrm{He}+\pi^-$   Equivalent 7D / 8D selections',fontsize=20,fontweight='bold')
    fig.text(.075,.89,r'Both: charged-mother PV refit, $\chi^2_{\mathrm{PV}}/2<5$, $L/\sigma_L>3$; $B_y=+0.5$ T',fontsize=13)
    fig.text(.075,.84,r'Counting window: $2.98531<M<2.99731$ GeV/$c^2$ ($M_0\pm6$ MeV/$c^2$)',fontsize=13,color=orange)
    for xpos,key in [(.075,'left'),(.58,'right')]:
        m=stats[key];fig.text(xpos,.14,f"In window: signal {m['signal_window']:,}, background {m['background_window']:,}",fontsize=14)
        fig.text(xpos,.09,f"S/B = {m['S_over_B']:.2f};  S/sqrt(S+B) = {m['S_over_sqrt_SB']:.2f}",fontsize=13)
    fig.text(.075,.035,'Same selected candidate IDs and histogram bins. 7D flight errors are calculated independently from the original input covariance.',fontsize=11,color='#59636e')
    fig.savefig(out/'CEE_Hypertriton_7D_8D_Equivalent.pdf');plt.close(fig)

if __name__=='__main__':main()
