"""EXPLORATORY (not in ANALYSIS_PLAN.md): robustness of the one secondary-sensitivity hit
(S3: DisorderImpact vs intellectual disability). Output: out/exploratory_disorder_ID.csv"""
import re
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
from common import OUT
from stats_utils import lrt_ml, lrt_firth

F = pd.read_csv(OUT / "individuals.csv").merge(pd.read_csv(OUT / "clusters.csv")[
    ["individualId", "source", "keep_one_per_family", "keep_one_per_family_or_dup"]], on="individualId")
pt = F[F.is_point == 1].copy()
pt["zD"] = (pt.DisorderImpact - pt.DisorderImpact.mean()) / pt.DisorderImpact.std(ddof=0)
pt["zS"] = (pt.log10_semi - pt.log10_semi.mean()) / pt.log10_semi.std(ddof=0)
h = "HP:0001249"
rows = []

def run(name, d, extra=(), expo="zD"):
    d = d[d[h].notna()]
    y = d[h].values.astype(float)
    base = [np.ones(len(d)), d.is_trunc.values, d.rel_pos.values] + [d[c].values for c in extra]
    X0 = np.column_stack([b for i, b in enumerate(base) if i == 0 or np.ptp(b) > 0])
    X1 = np.column_stack([X0, d[expo].values])
    s, _, p, (b, cov), conv = lrt_ml(X0, X1, y)
    fs, fp, fb, _ = lrt_firth(X1, y, [X1.shape[1] - 1])
    rows.append(dict(analysis=name, n=len(d), n_exc=int((y == 0).sum()), OR_perSD=np.exp(b[-1]), p_ml=p,
                     p_firth=fp, converged=conv))

run("as S3 (all point variants)", pt)
run("one per family", pt[pt.keep_one_per_family == 1])
run("one per family or likely dup", pt[pt.keep_one_per_family_or_dup == 1])
run("truncating only", pt[pt.is_trunc == 1])
# adjust for the largest sources (indicator for each source with >=10 ID-annotated individuals)
d = pt[pt[h].notna()].copy()
big = d.source.value_counts()
big = big[big >= 10].index.tolist()
for s in big:
    d[f"src_{s}"] = (d.source == s).astype(float)
run(f"adjusted for source ({len(big)} indicators)", d, extra=[f"src_{s}" for s in big])
run("adjusted for log10 SemiCore", pt.assign(lsemi=pt.zS), extra=["lsemi"])
run("adjusted for c.1903_1907del hotspot", pt, extra=["hotspot_1903_1907del"])
run("excluding c.1903_1907del carriers", pt[pt.hotspot_1903_1907del == 0])
import patsy
for df_ in (4, 6, 8):
    spl = np.asarray(patsy.dmatrix(f"bs(x, df={df_}, degree=3) - 1", {"x": pt.rel_pos.values}))
    q = pt.copy()
    for j in range(spl.shape[1]):
        q[f"bs{j}"] = spl[:, j]
    run(f"position as cubic B-spline df={df_}", q, extra=[f"bs{j}" for j in range(spl.shape[1])])
# leave-one-source-out
for s in big:
    run(f"leave out {s}", pt[pt.source != s])
# GEE / cluster-robust by family cluster
d = pt[pt[h].notna()].merge(pd.read_csv(OUT / "clusters.csv")[["individualId", "family_cluster"]], on="individualId")
X = sm.add_constant(d[["is_trunc", "rel_pos", "zD"]].astype(float))
r = sm.GLM(d[h].astype(float), X, family=sm.families.Binomial()).fit(cov_type="cluster",
                                                                   cov_kwds={"groups": pd.factorize(d.family_cluster)[0]})
rows.append(dict(analysis="cluster-robust SE by family cluster (Wald)", n=len(d), OR_perSD=np.exp(r.params.zD),
                 p_ml=r.pvalues.zD))
# rank-based, no model
dd = pt[pt[h].notna()]
u = stats.mannwhitneyu(dd.DisorderImpact[dd[h] == 1], dd.DisorderImpact[dd[h] == 0])
rows.append(dict(analysis="Mann-Whitney DisorderImpact by ID status (unadjusted)", n=len(dd), p_ml=u.pvalue,
                 n_exc=int((dd[h] == 0).sum())))
R = pd.DataFrame(rows)
R.to_csv(OUT / "exploratory_disorder_ID.csv", index=False)
print(R.to_string())
print(dd.groupby(h).DisorderImpact.describe())
# print(dd[dd[h] == 0][["individualId", "source", "vclass", "aa_pos", "DisorderImpact"]].sort_values("aa_pos").to_string())

# EXPLORATORY 2: MC2023 short-stature/exon-9 claim with CNVs counted as "outside exon 9"
F2 = pd.read_csv(OUT / "individuals.csv")
F2["grp"] = F2.exon9.map({1.0: "exon9 point", 0.0: "non-exon9 point"}).fillna("CNV")
d = F2[F2["HP:0004322"].notna()]
print(d.groupby("grp")["HP:0004322"].agg(["sum", "count", "mean"]))
a, b = d[d.grp == "exon9 point"]["HP:0004322"], d[d.grp != "exon9 point"]["HP:0004322"]
print("exon 9 vs everything else (incl. CNV):", round(a.mean(), 3), round(b.mean(), 3),
      stats.fisher_exact([[a.sum(), len(a) - a.sum()], [b.sum(), len(b) - b.sum()]]))
