"""EXPLORATORY (post hoc, written after the pre-registered results were seen).

Question: within variants where AlphaMissense is in its own 'ambiguous' class
(0.34 <= am < 0.564; Cheng et al. 2023), does neighbouring-label evidence separate P from B,
and does the propagation component reach ClinGen evidence strengths there?
Models are the pre-registered fits (trained on all TRAIN); only the subset is new.
"""
import json
import numpy as np, pandas as pd
import run_vus3d as R

OUT = R.OUT
rng = np.random.default_rng(R.SEED + 7)
F = pd.read_parquet(OUT / "features.parquet")
best = json.load(open(OUT / "tuned_params.json"))
tr = F[F.set == "train"].reset_index(drop=True)
te = F[(F.set == "test") & (F.n_lab_res >= 5)].reset_index(drop=True)
models = {"M0": None, "M1": None, "M2-1D": best["w1d"], "M2-3D": best["r3d"], "M2-RWR": best["rwr"]}
fits = {m: R.fit(tr, m, p) for m, p in models.items()}

def band(df):
    return (df.am_pathogenicity >= 0.34) & (df.am_pathogenicity < 0.564)

rows, contr = [], []
for name, sel in (("am_ambiguous", band(te)), ("am_confident", ~band(te))):
    t = te[sel].reset_index(drop=True)
    preds = {m: fits[m].predict_proba(R.design(t, m, p))[:, 1] for m, p in models.items()}
    for m in models:
        rows.append(dict(subset=name, model=m, n=len(t), n_P=int(t.y.sum()), n_genes=t.gene.nunique(),
                         strat_auroc=R.strat_auc(t.y.values, preds[m], t.gene.values),
                         auroc=R.roc_auc_score(t.y, preds[m])))
    B = R.boot_gene(t, preds, [("M2-RWR", "M1"), ("M2-3D", "M1"), ("M2-1D", "M1"), ("M2-RWR", "M2-1D")],
                    2000, rng, metric="strat")
    B.insert(0, "subset", name); contr.append(B)
pd.DataFrame(rows).to_csv(OUT / "explore_am_band_metrics.csv", index=False)
pd.concat(contr).to_csv(OUT / "explore_am_band_contrasts.csv", index=False)

# Calibration of the propagation component inside the AM-ambiguous band.
# TRAIN: out-of-fold (gene-grouped 5-fold) component, as in the pre-registered calibration.
par = best["rwr"]
genes = tr.gene.unique(); frng = np.random.default_rng(R.SEED)
fold = dict(zip(genes, frng.permutation(len(genes)) % 5)); f = tr.gene.map(fold).values
comp = np.empty(len(tr))
for i in range(5):
    k = f == i
    m2, m1 = R.fit(tr[~k], "M2-RWR", par), R.fit(tr[~k], "M1")
    comp[k] = m2.decision_function(R.design(tr[k], "M2-RWR", par)) - m1.decision_function(R.design(tr[k], "M1"))
comp_te = fits["M2-RWR"].decision_function(R.design(te, "M2-RWR", par)) - fits["M1"].decision_function(R.design(te, "M1"))
bt, bte = band(tr).values, band(te).values
edges = [-np.inf, -1, -0.5, 0.5, 1, np.inf]
THR = {"Supporting": 2.406, "Moderate": 5.790, "Strong": 33.53}

def lr(y, c, lo, hi):
    s = (c > lo) & (c <= hi)
    nP, nB = int(y[s].sum()), int((1 - y[s]).sum())
    return nP, nB, ((nP + .5) / (y.sum() + .5)) / ((nB + .5) / ((1 - y).sum() + .5))

cal = []
ytr, yte = tr.y.values[bt], te.y.values[bte]; ctr, cte = comp[bt], comp_te[bte]
gtr = tr.gene.values[bt]; ug, gi = np.unique(gtr, return_inverse=True)
for lo, hi in zip(edges[:-1], edges[1:]):
    a = lr(ytr, ctr, lo, hi); b = lr(yte, cte, lo, hi)
    # one-sided 95% lower (or upper, for benign direction) gene-bootstrap bound on TRAIN LR
    bs = []
    for _ in range(1000):
        c = np.bincount(rng.integers(0, len(ug), len(ug)), minlength=len(ug))[gi]
        s = (ctr > lo) & (ctr <= hi)
        P, Bn = (c * ytr).sum(), (c * (1 - ytr)).sum()
        bs.append(((c * ytr * s).sum() + .5) / (P + .5) / (((c * (1 - ytr) * s).sum() + .5) / (Bn + .5)))
    lo5, hi95 = np.percentile(bs, 5), np.percentile(bs, 95)
    strength = "none"
    for nm, t in THR.items():
        if lo5 >= t:
            strength = f"PM1-like {nm}"
        if hi95 <= 1 / t:
            strength = f"benign {nm}"
    cal.append(dict(bin=f"({lo},{hi}]", train_P=a[0], train_B=a[1], train_LR=a[2], train_LR_lo5=lo5,
                    train_LR_hi95=hi95, strength=strength, test_P=b[0], test_B=b[1], test_LR=b[2]))
pd.DataFrame(cal).to_csv(OUT / "explore_am_band_calibration.csv", index=False)
pd.set_option("display.width", 250)
print(pd.DataFrame(rows).round(4).to_string(index=False))
print(pd.concat(contr)[["subset", "contrast", "delta_auc", "ci_lo", "ci_hi", "p_boot"]].round(4).to_string(index=False))
print(pd.DataFrame(cal).round(3).to_string(index=False))
