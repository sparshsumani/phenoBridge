"""Isolates why mechanistic_mvp.py and ANKRD_mechanistic.ipynb disagree.

Refits the app's own fit_hpo_models_from_cohort() twice with no change to its
code: once as shipped (all 333 individuals, the 77 structural-variant carriers
pinned at logSemiCore = -15) and once restricted to the 256 individuals whose
variant actually maps to an amino-acid position. Also prints the SemiCoreImpact
region-definition mismatch between the two codebases. Read-only.

Run:  python audit/diagnose_divergence.py
"""
import os, sys, tempfile, textwrap
import matplotlib; matplotlib.use("Agg")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# mechanistic_mvp.py imports gradio at module scope purely for its UI. Stub it so
# the model-fitting code can be imported without installing the UI stack.
_stub = tempfile.mkdtemp()
with open(os.path.join(_stub, "gradio.py"), "w") as fh:
    fh.write(textwrap.dedent('''
        class _X:
            def __init__(self, *a, **k): pass
            def __enter__(self): return self
            def __exit__(self, *a): return False
            def click(self, *a, **k): pass
            def launch(self, *a, **k): pass
        def __getattr__(name): return _X
        Blocks = Row = Markdown = Number = Button = Code = Textbox = Plot = Dataframe = _X
    '''))
sys.path.insert(0, _stub)

import importlib.util, numpy as np, pandas as pd, statsmodels.api as sm, json, re
from statsmodels.stats.multitest import multipletests
spec = importlib.util.spec_from_file_location("mvp", os.path.join(ROOT, "mechanistic_mvp.py"))
mvp = importlib.util.module_from_spec(spec); spec.loader.exec_module(mvp)
JS = os.path.join(ROOT, "ANKRD11_KBGS_individuals.json")
m = mvp.build_model(os.path.join(ROOT, "AF-Q6UB99-F1-model_v6.pdb"), None, 8.0, 0.3)

# --- feature definition mismatch: mvp SemiCoreImpact (semi MINUS core) vs notebook (semi INCLUDING core)
semi_only=set(m.semi_nodes)-set(m.core_nodes)
print("mvp semi_only nodes:",len(semi_only),"| notebook semicore nodes: 1331")
for aa in [600,2400,1197]:
    p=mvp.rwr(m.G_all,aa,alpha=0.3)
    a=sum(p.get(n,0.) for n in semi_only); b=sum(p.get(n,0.) for n in m.semi_nodes)
    print(f"  aa={aa}: mvp SemiCoreImpact={a:.6g}  notebook-style(incl core)={b:.6g}  ratio={b/a if a else float('nan'):.3f}")

# --- rebuild the mvp design matrix, then refit on 333 vs 256
df_pheno,df_mut,lblmap=mvp.load_kbgs_json(JS)
dfv=df_mut.copy(); dfv["aa_pos"]=dfv["alleleKey"].map(mvp.allelekey_to_aa_pos)
dfv=dfv.dropna(subset=["aa_pos"]); dfv["aa_pos"]=dfv["aa_pos"].astype(int)
dfv["TRUNC_like"]=dfv["alleleKey"].map(mvp.is_truncating).astype(int)
fps=dfv["aa_pos"].map(lambda p: mvp.fingerprint_for_variant(m,int(p)))
dfv["SemiCoreImpact"]=fps.map(lambda d:0.0 if d is None else float(d["SemiCoreImpact"]))
mp=dfv.groupby("individualId")[["SemiCoreImpact"]].max().reset_index()
tp=dfv.groupby("individualId")["TRUNC_like"].mean().reset_index().rename(columns={"TRUNC_like":"TruncFrac"})
dm=df_pheno.merge(mp,on="individualId",how="left").merge(tp,on="individualId",how="left").fillna(0.0)
dm["logSemiCore"]=np.log10(dm["SemiCoreImpact"]+1e-15)
mapped=set(mp.individualId)
dm["is_mapped"]=dm.individualId.isin(mapped)
print("\nlogSemiCore distribution:")
print("  unmapped (SV) patients n=%d -> logSemiCore == -15 exactly: %d"%((~dm.is_mapped).sum(),int((dm.loc[~dm.is_mapped,'logSemiCore']==-15).sum())))
print("  mapped patients n=%d, logSemiCore range %.3f..%.3f"%(dm.is_mapped.sum(),dm.loc[dm.is_mapped,'logSemiCore'].min(),dm.loc[dm.is_mapped,'logSemiCore'].max()))

hpo=[c for c in df_pheno.columns if c!="individualId"]
def fit(frame,tag):
    rows=[]
    for h in hpo:
        yv=frame[h].dropna()
        if not len(yv): continue
        prev=float((yv==1.0).mean())
        if not (0.05<=prev<=0.95): continue
        y=frame[h]; msk=y.isin([0.0,1.0]); y=y[msk].astype(int)
        X=frame.loc[msk,["logSemiCore"]].copy(); X["TruncFrac"]=frame.loc[msk,"TruncFrac"].astype(float)
        X=sm.add_constant(X,has_constant="add")
        try: r=sm.Logit(y,X).fit(disp=0)
        except Exception: continue
        rows.append({"phenotype":h,"label":lblmap.get(h,""),"p":float(r.pvalues.get("logSemiCore",np.nan)),"n":len(y)})
    d=pd.DataFrame(rows).sort_values("p")
    q=multipletests(d.p.values,method="fdr_bh")[1] if len(d) else []
    print(f"\n=== {tag}: {len(d)} models | raw p<0.05: {int((d.p<0.05).sum())} | BH q<0.05: {int((q<0.05).sum())}")
    d2=d.copy(); d2["q"]=q
    print(d2.head(8).to_string(index=False))
    return d2
a=fit(dm,"APP AS SHIPPED (all 333, SV patients pinned at logSemiCore=-15)")
b=fit(dm[dm.is_mapped].reset_index(drop=True),"SAME CODE, restricted to the 256 aa-mapped patients")
