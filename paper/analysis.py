import sys, types, json, re
sys.modules['gradio'] = types.ModuleType('gradio')  # app UI not needed
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.abspath(__file__))))
import numpy as np, pandas as pd, networkx as nx
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import mechanistic_mvp as M

from pathlib import Path
R = str(Path(__file__).resolve().parent.parent)+'/'
out = {}
model = M.build_model(R+'AF-Q6UB99-F1-model_v6.pdb', None, 8.0, 0.3)
pos, xyz, pl = M.load_alphafold_ca(R+'AF-Q6UB99-F1-model_v6.pdb')
out['L'] = len(pos); out['plddt'] = dict(min=pl.min(), med=np.median(pl), mean=pl.mean(), max=pl.max(),
  q50=np.quantile(pl,.5), q90=np.quantile(pl,.9), frac70=(pl>=70).mean(), frac50=(pl<50).mean())
out['G_all'] = (model.G_all.number_of_nodes(), model.G_all.number_of_edges())
out['G_core'] = (model.G_core.number_of_nodes(), model.G_core.number_of_edges())
out['G_core_cc'] = nx.number_connected_components(model.G_core)
out['n_core'], out['n_semi_only'], out['n_dis'] = len(model.core_nodes), len(model.semi_nodes-model.core_nodes), len(model.disorder_nodes)
out['backbone'] = len(model.comm_backbone)
bb = sorted(model.comm_backbone); out['backbone_range'] = (bb[0], bb[-1])

# structured segments (pLDDT>=70 runs >= 20 aa)
segs=[]; s=None
for p,c in zip(pos,pl):
    if c>=70 and s is None: s=p
    if c<70 and s is not None:
        if p-s>=20: segs.append((s,p-1))
        s=None
if s is not None and pos[-1]-s>=19: segs.append((s,pos[-1]))
out['segments']=segs

# fingerprint for every residue
fp_all = pd.DataFrame([M.fingerprint_for_variant(model,int(p)) for p in pos])
fp_all.to_csv('fp_all.csv', index=False)

# cohort
df_pheno, df_mut, labels = M.load_kbgs_json(R+'ANKRD11_KBGS_individuals.json')
out['n_ind']=len(df_pheno); out['n_hpo']=df_pheno.shape[1]-1; out['n_allele']=df_mut.alleleKey.nunique()
dfv = df_mut.copy(); dfv['aa_pos']=dfv.alleleKey.map(M.allelekey_to_aa_pos)
out['n_sv']=int(dfv.alleleKey.str.contains('_SV_').sum())
dfv = dfv.dropna(subset=['aa_pos']); dfv['aa_pos']=dfv.aa_pos.astype(int)
dfv['TRUNC_like']=dfv.alleleKey.map(M.is_truncating).astype(int)
out['n_mapped_ind']=dfv.individualId.nunique(); out['n_mapped_pos']=dfv.aa_pos.nunique()
out['frac_trunc']=dfv.groupby('individualId').TRUNC_like.max().mean()
fpm = fp_all.set_index(fp_all.aa_pos.astype(int))
dfv['Semi']=dfv.aa_pos.map(fpm.SemiCoreImpact); dfv['Core']=dfv.aa_pos.map(fpm.CoreImpact)
dfv['BB']=dfv.aa_pos.map(fpm.CommBackboneImpact); dfv['pl']=dfv.aa_pos.map(fpm.seed_pLDDT)
dfv['Dis']=dfv.aa_pos.map(fpm.DisorderImpact)
pat = dfv.groupby('individualId').agg(Semi=('Semi','max'),Core=('Core','max'),BB=('BB','max'),
        Trunc=('TRUNC_like','mean'),aa=('aa_pos','min'),pl=('pl','max')).reset_index()
pat['logSemi']=np.log10(pat.Semi+1e-15)
out['semi_desc']=pat.Semi.describe().to_dict()
out['frac_semi_gt_0.5']=(pat.Semi>0.5).mean()
out['frac_seed_lowconf']=(pat.pl<50).mean()
out['n_var_in_struct']=int((pat.pl>=70).sum())
# variants per segment
out['var_hist']=np.histogram(dfv.aa_pos, bins=np.linspace(1,2663,21).astype(int))[0].tolist()

D = df_pheno.merge(pat, on='individualId', how='inner')
# association (A) pipeline as in hackathon: NA coded as 0 ; (B) explicit annotation only
def assoc(D, na_as_zero):
    rows=[]
    for h in [c for c in df_pheno.columns if c.startswith('HP:')]:
        y = D[h].fillna(0.0) if na_as_zero else D[h]
        m = y.isin([0,1]); y=y[m].astype(int)
        npos=int(y.sum()); n=len(y)
        if npos<15 or n-npos<15: continue
        X=sm.add_constant(D.loc[m,['logSemi','Trunc']])
        try: r=sm.Logit(y,X).fit(disp=0)
        except Exception: continue
        ci=r.conf_int().loc['logSemi']
        rows.append(dict(hpo=h,label=labels[h],n=n,npos=npos,OR=np.exp(r.params.logSemi),
             lo=np.exp(ci[0]),hi=np.exp(ci[1]),p=r.pvalues.logSemi))
    T=pd.DataFrame(rows); T['q']=multipletests(T.p,method='fdr_bh')[1]
    return T.sort_values('p').reset_index(drop=True)
A=assoc(D,True); B=assoc(D,False)
A.to_csv('assocA.csv',index=False); B.to_csv('assocB.csv',index=False)
out['A_ntests']=len(A); out['B_ntests']=len(B)
out['A_nom']=int((A.p<.05).sum()); out['B_nom']=int((B.p<.05).sum())
out['A_minq']=A.q.min(); out['B_minq']=B.q.min()
# missingness
miss = df_pheno.drop(columns='individualId').isna().mean()
out['miss_median']=miss.median(); out['miss_max']=miss.max()

# robustness of fingerprint to alpha / cutoff at cohort positions
uniq = sorted(dfv.aa_pos.unique())
base = np.array([fpm.SemiCoreImpact[p] for p in uniq])
from scipy.stats import spearmanr
rob={}
for a in [0.15,0.5]:
    m2=M.build_model(R+'AF-Q6UB99-F1-model_v6.pdb',None,8.0,a)
    v=np.array([M.fingerprint_for_variant(m2,int(p))['SemiCoreImpact'] for p in uniq]); rob[f'alpha={a}']=spearmanr(base,v).correlation
for c in [6.0,10.0]:
    m2=M.build_model(R+'AF-Q6UB99-F1-model_v6.pdb',None,c,0.3)
    v=np.array([M.fingerprint_for_variant(m2,int(p))['SemiCoreImpact'] for p in uniq]); rob[f'cutoff={c}']=spearmanr(base,v).correlation
out['robust']=rob
# fingerprints of demo positions
out['fp600']=M.fingerprint_for_variant(model,600); out['fp2400']=M.fingerprint_for_variant(model,2400)

# ---------- figures ----------
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.spines.top':False,'axes.spines.right':False,
  'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6})
C1,C2,C3='#2b5c8a','#c0504d','#7f7f7f'
fig,ax=plt.subplots(3,1,figsize=(7.0,3.6),sharex=True,gridspec_kw=dict(height_ratios=[1,1.2,0.8],hspace=0.12))
ax[0].fill_between(pos,pl,color=C1,alpha=.25,lw=0); ax[0].plot(pos,pl,color=C1,lw=.5)
ax[0].axhline(70,color=C3,lw=.5,ls='--'); ax[0].set_ylabel('pLDDT'); ax[0].set_ylim(0,100)
ax[1].plot(fp_all.aa_pos,fp_all.SemiCoreImpact,color=C1,lw=.6,label='SemiCoreImpact')
ax[1].plot(fp_all.aa_pos,fp_all.CoreImpact,color=C2,lw=.6,label='CoreImpact')
ax[1].set_ylabel('RWR mass'); ax[1].legend(frameon=False,loc='upper left',ncol=2,fontsize=7)
ax[2].hist(dfv.aa_pos,bins=80,color=C3); ax[2].set_ylabel('Variants'); ax[2].set_xlabel('ANKRD11 residue (UniProt Q6UB99)')
ax[2].set_xlim(1,2663)
for a in ax:
    for s,e in segs: a.axvspan(s,e,color='#f2c14e',alpha=.18,lw=0)
fig.savefig('fig/fig1.pdf',bbox_inches='tight'); plt.close()

fig,ax=plt.subplots(1,2,figsize=(3.4,1.9),gridspec_kw=dict(width_ratios=[1.1,1],wspace=0.55))
keys=['CoreImpact','SemiCoreImpact','DisorderImpact','CommBackboneImpact']; lab=['Core','Semi','Dis.','Backb.']
x=np.arange(4)
ax[0].bar(x-.2,[out['fp600'][k] for k in keys],.4,color=C1,label='res. 600')
ax[0].bar(x+.2,[out['fp2400'][k] for k in keys],.4,color=C2,label='res. 2400')
ax[0].set_xticks(x); ax[0].set_xticklabels(lab,fontsize=7); ax[0].set_ylabel('RWR mass'); ax[0].legend(frameon=False,fontsize=6.5)
ax[0].set_title('(a) fingerprints',fontsize=8)
ax[1].scatter(pat.pl,pat.logSemi,s=5,color=C1,alpha=.6,lw=0)
ax[1].set_xlabel('seed pLDDT'); ax[1].set_ylabel(r'$\log_{10}$ SemiCoreImpact'); ax[1].set_title('(b) cohort variants',fontsize=8)
fig.savefig('fig/fig2.pdf',bbox_inches='tight'); plt.close()

# forest plot for B top 12 + A marker
T=B.copy().head(14).iloc[::-1]
fig,ax=plt.subplots(figsize=(3.4,2.6))
yy=np.arange(len(T))
ax.errorbar(T.OR,yy,xerr=[T.OR-T.lo,T.hi-T.OR],fmt='o',ms=3,color=C1,ecolor=C1,elinewidth=.7,capsize=1.5)
ax.axvline(1,color=C3,lw=.6,ls='--'); ax.set_xscale('log')
ax.set_yticks(yy); ax.set_yticklabels([f"{l[:30]} ({n})" for l,n in zip(T.label,T.n)],fontsize=6.3)
ax.set_xlabel('OR per 10-fold increase in SemiCoreImpact')
fig.savefig('fig/fig3.pdf',bbox_inches='tight'); plt.close()

def conv(o):
    if isinstance(o,(np.integer,)): return int(o)
    if isinstance(o,(np.floating,)): return float(o)
    if isinstance(o,tuple): return list(o)
    return str(o)
json.dump(out,open('results.json','w'),indent=1,default=conv)
print(json.dumps(out,indent=1,default=conv))
print(A.head(10).to_string()); print(B.head(14).to_string())
