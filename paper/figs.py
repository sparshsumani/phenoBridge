import sys, types, json
sys.modules['gradio']=types.ModuleType('gradio'); sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.abspath(__file__))))
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator, FixedFormatter
import mechanistic_mvp as M
from pathlib import Path
R=str(Path(__file__).resolve().parent.parent)+'/'; res=json.load(open('results.json'))
pos,xyz,pl=M.load_alphafold_ca(R+'AF-Q6UB99-F1-model_v6.pdb')
fp=pd.read_csv('fp_all.csv'); fpm=fp.set_index(fp.aa_pos.astype(int)); segs=res['segments']
_,df_mut,_=M.load_kbgs_json(R+'ANKRD11_KBGS_individuals.json')
dfv=df_mut.copy(); dfv['aa_pos']=dfv.alleleKey.map(M.allelekey_to_aa_pos); dfv=dfv.dropna(subset=['aa_pos']); dfv['aa_pos']=dfv.aa_pos.astype(int)
dfv['Semi']=dfv.aa_pos.map(fpm.SemiCoreImpact); dfv['pl']=dfv.aa_pos.map(fpm.seed_pLDDT)
plt.rcParams.update({'font.family':'serif','font.serif':['FreeSerif'],'mathtext.fontset':'stix','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'svg.fonttype':'path'})
C1,C2,C3='#2b5c8a','#c0504d','#7f7f7f'
fig,ax=plt.subplots(3,1,figsize=(7.0,2.9),sharex=True,gridspec_kw=dict(height_ratios=[1,1.2,.8],hspace=.14))
ax[0].fill_between(pos,pl,color=C1,alpha=.25,lw=0); ax[0].plot(pos,pl,color=C1,lw=.5); ax[0].axhline(70,color=C3,lw=.5,ls='--')
ax[0].set_ylabel('pLDDT'); ax[0].set_ylim(0,100)
ax[1].plot(fp.aa_pos,fp.SemiCoreImpact,color=C1,lw=.55,label='SemiCoreImpact'); ax[1].plot(fp.aa_pos,fp.CoreImpact,color=C2,lw=.6,label='CoreImpact')
ax[1].set_ylabel('RWR mass'); ax[1].set_ylim(0,1.32); ax[1].set_yticks([0,.5,1])
ax[1].legend(frameon=False,loc='upper center',ncol=2,fontsize=7,bbox_to_anchor=(.5,1.04))
ax[2].hist(dfv.aa_pos,bins=80,color=C3); ax[2].set_ylabel('Variants'); ax[2].set_xlabel('ANKRD11 residue (UniProt Q6UB99)'); ax[2].set_xlim(1,2663)
for a in ax:
    for s,e in segs: a.axvspan(s,e,color='#f2c14e',alpha=.22,lw=0)
for a,t in zip(ax,'abc'): a.text(.006,.97,f'({t})',transform=a.transAxes,fontweight='bold',va='top',ha='left')
fig.savefig('fig/fig1.svg',bbox_inches='tight'); plt.close()

fig,ax=plt.subplots(1,2,figsize=(3.4,1.6),gridspec_kw=dict(width_ratios=[1.1,1],wspace=.6))
keys=['CoreImpact','SemiCoreImpact','DisorderImpact','CommBackboneImpact']; x=np.arange(4)
ax[0].bar(x-.2,[res['fp600'][k] for k in keys],.4,color=C1,label='res. 600')
ax[0].bar(x+.2,[res['fp2400'][k] for k in keys],.4,color=C2,label='res. 2400')
ax[0].set_xticks(x); ax[0].set_xticklabels(['Core','Semi','Dis.','Bb.'],fontsize=7); ax[0].set_ylabel('RWR mass'); ax[0].set_ylim(0,1.5); ax[0].set_yticks([0,.5,1])
ax[0].legend(frameon=False,fontsize=6.5,loc='upper center',ncol=2,bbox_to_anchor=(.5,1.05),handlelength=1,columnspacing=.8)
ax[0].set_title('(a) example fingerprints',fontsize=8)
ax[1].scatter(dfv.pl,np.log10(dfv.Semi),s=5,color=C1,alpha=.55,lw=0); ax[1].axvline(70,color=C3,lw=.5,ls='--')
ax[1].set_xlabel('seed pLDDT'); ax[1].set_ylabel(r'$\log_{10}$ SemiCoreImpact'); ax[1].set_title('(b) cohort variants',fontsize=8)
fig.savefig('fig/fig2.svg',bbox_inches='tight'); plt.close()

B=pd.read_csv('assocB.csv').head(14).iloc[::-1]; yy=np.arange(len(B))
fig,ax=plt.subplots(figsize=(3.4,2.45))
ax.errorbar(B.OR,yy,xerr=[B.OR-B.lo,B.hi-B.OR],fmt='o',ms=3,color=C1,elinewidth=.7,capsize=1.5)
ax.axvline(1,color=C3,lw=.6,ls='--'); ax.set_xscale('log')
ax.xaxis.set_major_locator(FixedLocator([.2,.5,1,2,5])); ax.xaxis.set_major_formatter(FixedFormatter(['0.2','0.5','1','2','5'])); ax.xaxis.set_minor_locator(NullLocator())
short={'Attention deficit hyperactivity disorder':'ADHD','Delayed speech and language development':'Speech/language delay','Delayed skeletal maturation':'Delayed bone age'}
ax.set_yticks(yy); ax.set_yticklabels([f"{short.get(l,l)} (n={n})" for l,n in zip(B.label,B.n)],fontsize=6.5)
ax.set_xlabel('OR per 10-fold increase in SemiCoreImpact')
fig.savefig('fig/fig3.svg',bbox_inches='tight'); plt.close()
