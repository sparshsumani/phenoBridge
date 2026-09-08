"""Verification harness: linear top-to-bottom re-execution of ANKRD_mechanistic.ipynb.

Runs the notebook's logic against the two raw files in the repository root and
prints every quantity the notebook recorded, so the recorded outputs can be
diffed against a fresh run. Read-only: it changes no pipeline code and writes
no files. Cell numbers in the [cNN] prefixes refer to ANKRD_mechanistic.ipynb.
"""
import json, re, random, math
from pathlib import Path
import numpy as np, pandas as pd, networkx as nx
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parent.parent
JSON=str(ROOT/"ANKRD11_KBGS_individuals.json")
PDB=str(ROOT/"AF-Q6UB99-F1-model_v6.pdb")
data=json.loads(Path(JSON).read_text())

# ---- cell 2 ----
df_hpo_headers=pd.DataFrame(data["hpoHeaders"]); hpo_ids=df_hpo_headers["hpoId"].tolist()
ind_rows=[];ind_dis=[]
for r in data["rows"]:
    ind=r.get("individualData",{})
    ind_rows.append({k:ind.get(k) for k in ["individualId","sex","pmid","title","comment","ageOfOnset","ageAtLastEncounter","deceased"]})
    for did in r.get("diseaseIdList",[]) or []: ind_dis.append({"individualId":ind.get("individualId"),"diseaseId":did})
df_individuals=pd.DataFrame(ind_rows).drop_duplicates(subset=["individualId"]).reset_index(drop=True)
df_individual_diseases=pd.DataFrame(ind_dis).drop_duplicates().reset_index(drop=True)
allele_rows=[]
for r in data["rows"]:
    i=(r.get("individualData") or {}).get("individualId")
    for k,c in (r.get("alleleCountMap") or {}).items(): allele_rows.append({"individualId":i,"alleleKey":k,"alleleCount":int(c)})
df_alleles=pd.DataFrame(allele_rows)
df_allele_matrix=df_alleles.pivot_table(index="individualId",columns="alleleKey",values="alleleCount",aggfunc="sum",fill_value=0).reset_index()
TYPE_TO_VALUE={"Observed":1,"Excluded":0,"Na":np.nan,"NA":np.nan,None:np.nan}
pheno_rows=[]
for r in data["rows"]:
    i=(r.get("individualData") or {}).get("individualId"); hd=r.get("hpoData",[]) or []
    hd=(hd+[{"type":"Na"}]*(len(hpo_ids)-len(hd)))[:len(hpo_ids)]
    row={"individualId":i}
    for k,cell in enumerate(hd): row[hpo_ids[k]]=TYPE_TO_VALUE.get((cell or {}).get("type"),np.nan)
    pheno_rows.append(row)
df_pheno=pd.DataFrame(pheno_rows)
df_one_disease=df_individual_diseases.sort_values(["individualId","diseaseId"]).groupby("individualId",as_index=False).first()
df_model=df_individuals.merge(df_one_disease,on="individualId",how="left").merge(df_allele_matrix,on="individualId",how="left").merge(df_pheno,on="individualId",how="left")
print("[c2] rows:",len(data['rows']),"pheno matrix:",df_pheno.shape,"model:",df_model.shape)

# ---- cell 3/4 ----
df_mut=df_alleles.dropna(subset=["alleleKey"]).drop_duplicates().copy(); df_mut["alleleKey"]=df_mut["alleleKey"].astype(str)
def parse_protein_pos(s):
    m=re.search(r"p\.?\(?[A-Za-z]{1,3}(\d+)[A-Za-z]{1,3}\)?",str(s)); return int(m.group(1)) if m else np.nan
def parse_cdna_pos(s):
    m=re.search(r"c\.?(\d+)",str(s)); return int(m.group(1)) if m else np.nan
df_mut2=df_mut.copy()
df_mut2["aa_pos_protein"]=df_mut2["alleleKey"].apply(parse_protein_pos)
df_mut2["cdna_pos"]=df_mut2["alleleKey"].apply(parse_cdna_pos)
df_mut2["aa_pos_cdna"]=df_mut2["cdna_pos"].apply(lambda c: np.nan if pd.isna(c) else int(np.ceil(c/3.0)))
df_mut2["aa_pos"]=df_mut2["aa_pos_protein"].combine_first(df_mut2["aa_pos_cdna"])
def is_trunc(s): return bool(re.search(r"(stop|ter|\*|nonsense|frameshift|fs|del|dup|ins)",str(s),re.I))
def is_miss(s):
    s=str(s); return bool(re.search(r"p\.?\(?[A-Za-z]{1,3}\d+[A-Za-z]{1,3}\)?",s) and not re.search(r"(ter|\*|fs|stop)",s,re.I))
df_mut2["MISSENSE_like"]=df_mut2["alleleKey"].apply(is_miss).astype(int)
df_mut2["TRUNC_like"]=df_mut2["alleleKey"].apply(is_trunc).astype(int)
print("[c4] df_mut2:",df_mut2.shape,"aa_pos present:",int(df_mut2.aa_pos.notna().sum()),"/",len(df_mut2),
      "range:",df_mut2.aa_pos.min(),df_mut2.aa_pos.max(),"MISSENSE_like sum:",int(df_mut2.MISSENSE_like.sum()))

# ---- cell 7/8 ----
chain=next(next(PDBParser(QUIET=True).get_structure("AF",PDB).get_models()).get_chains())
ca_xyz=[];plddt=[];unip_pos=[]
for res in chain.get_residues():
    if res.get_id()[0]!=" " or "CA" not in res: continue
    ca_xyz.append(res["CA"].coord.astype(float)); plddt.append(float(res["CA"].bfactor)); unip_pos.append(int(res.get_id()[1]))
ca_xyz=np.asarray(ca_xyz);plddt=np.asarray(plddt);unip_pos=np.asarray(unip_pos)
print(f"[c8] CA residues: {len(plddt)}  pos {unip_pos.min()}-{unip_pos.max()}  pLDDT min/med/mean/max: {plddt.min()} {np.median(plddt)} {plddt.mean():.6f} {plddt.max()}")

# ---- cell 9 ----
def build_graph(ca,pos,pl,cut=8.0,minp=0.0):
    keep=pl>=minp; idx=np.where(keep)[0]; xyz=ca[keep]
    G=nx.Graph()
    for oi,up,c in zip(idx,pos[keep],pl[keep]): G.add_node(int(oi),unip_pos=int(up),plddt=float(c))
    for a,b in cKDTree(xyz).query_pairs(r=cut):
        i,j=int(idx[a]),int(idx[b]); G.add_edge(i,j,dist=float(np.linalg.norm(ca[i]-ca[j])))
    return G
G_all=build_graph(ca_xyz,unip_pos,plddt,8.0,0.0); G_core=build_graph(ca_xyz,unip_pos,plddt,8.0,70.0)
print("[c9] G_all",G_all.number_of_nodes(),G_all.number_of_edges()," G_core",G_core.number_of_nodes(),G_core.number_of_edges())
print("[c10] structured>=70:",sum(1 for _,a in G_all.nodes(data=True) if a['plddt']>=70),
      " disordered<=50:",sum(1 for _,a in G_all.nodes(data=True) if a['plddt']<=50))

def shortest_path_backbone(G,sources,targets,max_pairs=2000,top_k=200,seed=0):
    rng=random.Random(seed); sources=list(sources);targets=list(targets)
    if not sources or not targets: return set()
    pairs=[(s,t) for s in sources for t in targets]
    if len(pairs)>max_pairs: pairs=rng.sample(pairs,k=max_pairs)
    counts={n:0 for n in G.nodes()}
    for s,t in pairs:
        try: path=nx.shortest_path(G,s,t)
        except nx.NetworkXNoPath: continue
        for n in path[1:-1]: counts[n]+=1
    return set(sorted(counts,key=lambda n:counts[n],reverse=True)[:top_k])
Na_={n for n,a in G_core.nodes(data=True) if a["unip_pos"]<=200}; Ca_={n for n,a in G_core.nodes(data=True) if a["unip_pos"]>=2463}
print("[c11] comm_backbone size:",len(shortest_path_backbone(G_core,Na_,Ca_,2000,200,0)))

def rwr(G,seed_node,alpha=0.3,max_iter=200,tol=1e-10):
    nodes=list(G.nodes()); idx={n:i for i,n in enumerate(nodes)}; n=len(nodes)
    nbrs=[list(G.neighbors(x)) for x in nodes]; deg=np.array([max(1,len(nbrs[i])) for i in range(n)],float)
    p0=np.zeros(n); p0[idx[seed_node]]=1.0; p=p0.copy()
    for _ in range(max_iter):
        pn=np.zeros_like(p)
        for i in range(n):
            s=0.0
            for nb in nbrs[i]: s+=p[idx[nb]]/deg[idx[nb]]
            pn[i]=(1-alpha)*s+alpha*p0[i]
        if np.linalg.norm(pn-p,1)<tol: p=pn;break
        p=pn
    return {nodes[i]:float(p[i]) for i in range(n)}

# ---- cell 23: FINAL region definitions ----
q90=float(np.quantile(plddt,0.90)); q50=float(np.quantile(plddt,0.50))
core_nodes={n for n,a in G_all.nodes(data=True) if a["plddt"]>=q90}
disorder_nodes={n for n,a in G_all.nodes(data=True) if a["plddt"]<=q50}
semicore_nodes=set(G_all.nodes())-disorder_nodes
print(f"[c23] q50={q50} q90={q90} core={len(core_nodes)} semicore={len(semicore_nodes)} disorder={len(disorder_nodes)}")
# ---- cell 24: backbone used for logBackbone ----
G_coreQ=G_all.subgraph(core_nodes).copy()
print("[c24] G_coreQ",G_coreQ.number_of_nodes(),G_coreQ.number_of_edges())
if G_coreQ.number_of_edges()<50: G_coreQ=G_all.subgraph(semicore_nodes).copy()
comm_backbone=shortest_path_backbone(G_coreQ,{n for n,a in G_coreQ.nodes(data=True) if a["unip_pos"]<=200},
                                      {n for n,a in G_coreQ.nodes(data=True) if a["unip_pos"]>=2463},1500,200,0)
print("[c24] backbone nodes:",len(comm_backbone))

pos_to_node={a["unip_pos"]:n for n,a in G_all.nodes(data=True)}
cache={}
def fp(aa):
    aa=int(aa)
    if aa not in cache:
        s=pos_to_node.get(aa)
        if s is None: cache[aa]=None
        else:
            infl=rwr(G_all,s,0.3)
            cache[aa]={"CoreImpact":sum(infl.get(n,0.) for n in core_nodes),
                       "SemiCoreImpact":sum(infl.get(n,0.) for n in semicore_nodes),
                       "DisorderImpact":sum(infl.get(n,0.) for n in disorder_nodes),
                       "CommBackboneImpact":sum(infl.get(n,0.) for n in comm_backbone)}
    return cache[aa]

# ---- cells 15/18/19/20/21/26 ----
mech0=(df_mut2.dropna(subset=["aa_pos"]).assign(aa_pos=lambda d:d.aa_pos.astype(int)))
df_model2=df_model.merge(mech0.groupby("individualId")[["alleleCount"]].max().reset_index().drop(columns=["alleleCount"]),on="individualId",how="left").fillna(0.0)
variant_cols=[c for c in df_model2.columns if c.startswith("c") or c.startswith("ANKRD11_SV_")]
long_var=df_model2[["individualId"]+variant_cols].melt(id_vars=["individualId"],var_name="variant_key",value_name="present")
long_var=long_var[long_var["present"]==1].copy()
print("[c19] patient-variant rows:",len(long_var),"| variant_cols included:",len(variant_cols))
map_df=df_mut2[["alleleKey","aa_pos"]].dropna().copy(); map_df["aa_pos"]=map_df["aa_pos"].astype(int)
long_var2=long_var.merge(map_df,left_on="variant_key",right_on="alleleKey",how="left")
print("[c20] long_var2 rows:",len(long_var2),"mapped frac:",round(long_var2.aa_pos.notna().mean(),6))
sub2=long_var2.dropna(subset=["aa_pos"]).copy(); sub2["aa_pos"]=sub2["aa_pos"].astype(int)
print("[c22] sub rows:",len(sub2),"unique aa_pos:",sub2.aa_pos.nunique(),"individuals mapped:",sub2.individualId.nunique(),"total:",df_model2.individualId.nunique())
f=sub2["aa_pos"].map(fp)
for col in ["CoreImpact","SemiCoreImpact","DisorderImpact","CommBackboneImpact"]: sub2[col]=f.map(lambda d: None if d is None else d[col])
mech_patient=sub2.groupby("individualId")[["CoreImpact","SemiCoreImpact","DisorderImpact","CommBackboneImpact"]].max().reset_index()
df_model3=df_model2.drop(columns=["CoreImpact","SemiCoreImpact","DisorderImpact","CommBackboneImpact"],errors="ignore").merge(mech_patient,on="individualId",how="left").fillna(0.0)
EPS=1e-15
for a,b in [("logCore","CoreImpact"),("logBackbone","CommBackboneImpact"),("logSemiCore","SemiCoreImpact")]:
    df_model3[a]=np.log10(df_model3[b]+EPS)
print("[c31] corr(SemiCore,Disorder):",round(float(np.corrcoef(df_model3.SemiCoreImpact,df_model3.DisorderImpact)[0,1]),10))
print("[c32] backbone zeros:",round(float((df_model3.CommBackboneImpact==0).mean()),10),"core zeros:",round(float((df_model3.CoreImpact==0).mean()),10))
df_mech=df_model3[df_model3.individualId.isin(set(mech_patient.individualId))].copy()
print("[c36] mechanistic cohort size:",len(df_mech))
pheno_cols=[c for c in df_mech.columns if c.startswith("HP:")]
valid=[c for c in pheno_cols if 15<=int(df_mech[c].sum())<=len(df_mech)-15]
print("[c34/37] pheno cols:",len(pheno_cols),"valid:",len(valid))
tf=(sub2.groupby("individualId")["aa_pos"].min()/2663.0).reset_index(name="TruncFrac")
df_mech=df_mech.merge(tf,on="individualId",how="left"); df_mech["TruncFrac"]=df_mech["TruncFrac"].fillna(1.0)
hp=pd.DataFrame({"hpo":pheno_cols,"n_pos":[int(df_mech[c].sum()) for c in pheno_cols]})
hp["prevalence"]=hp.n_pos/len(df_mech); hp=hp.sort_values(["prevalence","n_pos"],ascending=False)
hpv=hp[hp.hpo.isin(valid)]
print("[c41] top10 prevalence:"); [print(f"   {r.hpo}: {r.n_pos}/{len(df_mech)} ({100*r.prevalence:.1f}%)") for r in hpv.head(10).itertuples()]
top5=hpv.head(5)["hpo"].tolist()
print("[c42] unadjusted:")
for ph in top5:
    m=sm.Logit(df_mech[ph],sm.add_constant(df_mech[["logSemiCore"]])).fit(disp=0)
    print(f"   {ph} coef={m.params['logSemiCore']:.6f} OR={np.exp(m.params['logSemiCore']):.6f} p={m.pvalues['logSemiCore']:.6f}")
print("[c43] adjusted +TruncFrac:")
for ph in top5:
    m=sm.Logit(df_mech[ph],sm.add_constant(df_mech[["logSemiCore","TruncFrac"]])).fit(disp=0)
    print(f"   {ph} coef={m.params['logSemiCore']:.6f} OR={np.exp(m.params['logSemiCore']):.6f} p={m.pvalues['logSemiCore']:.6f}")
res=[]
for ph in valid:
    try:
        m=sm.Logit(df_mech[ph],sm.add_constant(df_mech[["logSemiCore","TruncFrac"]])).fit(disp=0)
        res.append({"phenotype":ph,"coef":m.params["logSemiCore"],"OR":np.exp(m.params["logSemiCore"]),"pval":m.pvalues["logSemiCore"]})
    except Exception: continue
dr=pd.DataFrame(res); dr["p_fdr"]=multipletests(dr.pval,method="fdr_bh")[1]
print("[c45] all valid + FDR (top 10 by p_fdr):"); print(dr.sort_values("p_fdr").head(10).to_string(index=False))
m=sm.Logit(df_mech["HP:0001249"],sm.add_constant(df_mech[["logSemiCore","logBackbone","TruncFrac"]])).fit(disp=0)
print("[c49] ID ~ logSemiCore+logBackbone+TruncFrac:")
print("   n=",int(m.nobs),"pseudoR2=",round(float(m.prsquared),5),"LLR p=",round(float(m.llr_pvalue),6))
print(pd.DataFrame({"coef":m.params,"p":m.pvalues,"OR":np.exp(m.params)}).round(6).to_string())
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
X=df_mech[["logSemiCore","logBackbone","TruncFrac"]]; y=df_mech["HP:0001249"]
Xtr,Xte,ytr,yte=train_test_split(X,y,test_size=0.3,random_state=42)
t=DecisionTreeClassifier(max_depth=3,min_samples_leaf=15,random_state=42).fit(Xtr,ytr)
print("[c53] decision tree:"); print(classification_report(yte,t.predict(Xte)))
print("   majority-class baseline acc =",round(float(max(yte.mean(),1-yte.mean())),4))
print("[c54] importances:",dict(zip(X.columns,t.feature_importances_.round(6))))
