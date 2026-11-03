"""Allelic-disorder analysis, implementing allelic/PLAN.md (+ Amendment 1).

Unit: labelled residue (variants collapsed to residues). Genes are the units of global
inference. Run after the data pipeline (a0*_*.py) has written derived_allelic/.

    python3 run_allelic.py [--perms 1000] [--calib 200]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse.csgraph import shortest_path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402

OUT = HERE / "out"
SEED = 20261105
MIN_RES = 5

POS_CTRL = ["LMNA", "PTPN11", "FGFR3", "RET", "KCNQ1", "SCN4A", "GJB2", "GCK", "CASR", "TP63",
            "SMCHD1", "STAT1", "FLNA", "SCN9A", "PIEZO1", "AARS1", "CLCN7", "TRPV4", "TGFBI"]
NULL_CTRL = ["FGFR2", "PRNP", "VCP", "ABCA4"]


# ------------------------------------------------------------------ data
def load_variants():
    v = pd.read_parquet(C.OUT / "variants.parquet")
    return v


def residue_labels(v, majority=False):
    """One row per (gene, residue) with a disease label.
    Primary: all single-label variants at the residue agree. majority=True: majority label,
    ties dropped."""
    s = v[v.single_label].copy()
    g = s.groupby(["gene", "uniprot", "pos"]).disease_unit.agg(list).reset_index()
    out = []
    for r in g.itertuples(index=False):
        vals = pd.Series(r.disease_unit).value_counts()
        if len(vals) == 1:
            lab, conflict = vals.index[0], False
        elif majority and vals.iloc[0] > vals.iloc[1]:
            lab, conflict = vals.index[0], True
        else:
            lab, conflict = None, True
        out.append(dict(gene=r.gene, uniprot=r.uniprot, pos=int(r.pos), label=lab, conflict=conflict,
                        n_variants=int(vals.sum())))
    R = pd.DataFrame(out)
    return R


def eligible(R, min_res=MIN_RES):
    R = R[R.label.notna()]
    keep = []
    for g, d in R.groupby("gene"):
        vc = d.label.value_counts()
        ok = vc[vc >= min_res].index
        if len(ok) >= 2:
            keep.append(d[d.label.isin(ok)])
    return pd.concat(keep, ignore_index=True) if keep else R.iloc[:0]


# ------------------------------------------------------------------ distances
def distances(uniprot, positions, alphas=(0.3,)):
    """Return dict name -> (n x n) distance matrix among the given residues."""
    pos_all, xyz, _ = C.ca_coords(uniprot)
    idx = {p: k for k, p in enumerate(pos_all)}
    ii = np.array([idx[p] for p in positions])
    P = np.asarray(positions, float)
    D = {"seq": np.abs(P[:, None] - P[None, :])}
    X = xyz[ii]
    D["3d"] = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    _, _, A = C.contact_graph(uniprot)
    hop = shortest_path(A, method="D", unweighted=True, directed=False, indices=ii)
    D["hop"] = hop[:, ii]
    for a in alphas:
        M = np.vstack([C.rwr_profile(uniprot, p, alpha=a).to_numpy() for p in positions])
        Mn = M / M.max(axis=1, keepdims=True)
        S = Mn[:, ii]
        S = 0.5 * (S + S.T)
        D[f"rwr{a}"] = 1.0 - S
    return D


# ------------------------------------------------------------------ kNN
def knn_predict(Dm, y, k, rng):
    """Leave-one-out kNN. Neighbours: all j != i with d_ij <= k-th smallest distance (ties
    included). Majority vote; label ties -> label whose closest member is nearest; then random."""
    n = len(y)
    pred = np.empty(n, dtype=object)
    for i in range(n):
        d = Dm[i].copy()
        d[i] = np.inf
        kth = np.partition(d, k - 1)[k - 1]
        nb = np.where(d <= kth)[0]
        labs, cnt = np.unique(y[nb], return_counts=True)
        top = labs[cnt == cnt.max()]
        if len(top) > 1:
            best = np.array([d[nb][y[nb] == t].min() for t in top])
            top = top[best == best.min()]
        pred[i] = top[0] if len(top) == 1 else top[rng.integers(len(top))]
    return pred


def balanced_accuracy(y, pred):
    labs = np.unique(y)
    return float(np.mean([(pred[y == l] == l).mean() for l in labs]))


def ba_all(D, y, k, rng):
    return {m: balanced_accuracy(y, knn_predict(Dm, y, k, rng)) for m, Dm in D.items()}


# ------------------------------------------------------------------ primary
def gene_scores(E, k=3, alphas=(0.3,), perms=1000, rng=None, keep_D=False):
    rows, Dstore = [], {}
    for g, d in E.groupby("gene"):
        d = d.sort_values("pos")
        uni = d.uniprot.iloc[0]
        pos = d.pos.to_numpy(int)
        y = d.label.to_numpy(object)
        D = distances(uni, pos, alphas)
        if keep_D:
            Dstore[g] = (D, y, pos)
        obs = ba_all(D, y, k, rng)
        null = {m: np.empty(perms) for m in D}
        for b in range(perms):
            yp = rng.permutation(y)
            for m, Dm in D.items():
                null[m][b] = balanced_accuracy(yp, knn_predict(Dm, yp, k, rng))
        rec = dict(gene=g, uniprot=uni, n_res=len(y), n_dis=len(np.unique(y)))
        for m in D:
            rec[f"ba_{m}"] = obs[m]
            rec[f"chance_{m}"] = float(null[m].mean())
            rec[f"p_{m}"] = float((1 + (null[m] >= obs[m]).sum()) / (1 + perms))
        rows.append(rec)
    return pd.DataFrame(rows), Dstore


def paired_test(G, a, b, rng, nboot=5000):
    d = (G[f"ba_{a}"] - G[f"ba_{b}"]).to_numpy()
    nz = d[d != 0]
    p = float(stats.wilcoxon(nz).pvalue) if len(nz) >= 5 else np.nan
    boot = np.array([np.median(rng.choice(d, len(d))) for _ in range(nboot)])
    return dict(contrast=f"{a} - {b}", n_genes=len(d), median_diff=float(np.median(d)),
                mean_diff=float(d.mean()), ci_lo=float(np.percentile(boot, 2.5)),
                ci_hi=float(np.percentile(boot, 97.5)), n_better=int((d > 0).sum()),
                n_worse=int((d < 0).sum()), n_tie=int((d == 0).sum()), p_wilcoxon=p)


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0
    for r, i in enumerate(o):
        run = max(run, (m - r) * p[i])
        adj[i] = min(1, run)
    return adj


def primary_tests(G, rwr="rwr0.3", rng=None):
    T = [paired_test(G, "3d", "seq", rng), paired_test(G, rwr, "3d", rng), paired_test(G, rwr, "hop", rng)]
    T = pd.DataFrame(T)
    T.insert(0, "hypothesis", ["H1", "H2", "H2b"])
    T["p_holm"] = holm(T.p_wilcoxon.values)
    return T


def calibration(Dstore, k, ncal, rng):
    """Full label permutations within every gene; distribution of median differences and
    Wilcoxon p for H1/H2/H2b."""
    out = []
    for b in range(ncal):
        rows = []
        for g, (D, y, _) in Dstore.items():
            yp = rng.permutation(y)
            rows.append(dict(gene=g, **{f"ba_{m}": v for m, v in ba_all(D, yp, k, rng).items()}))
        Gp = pd.DataFrame(rows)
        for name, a, bb in (("H1", "3d", "seq"), ("H2", "rwr0.3", "3d"), ("H2b", "rwr0.3", "hop")):
            d = (Gp[f"ba_{a}"] - Gp[f"ba_{bb}"]).to_numpy()
            nz = d[d != 0]
            out.append(dict(rep=b, hypothesis=name, median_diff=float(np.median(d)),
                            p=float(stats.wilcoxon(nz).pvalue) if len(nz) >= 5 else np.nan))
    return pd.DataFrame(out)


# ------------------------------------------------------------------ secondary
FEATS = {"plddt": "plddt", "core": "core_a0.3", "semi": "semi_a0.3", "disorder": "disorder_a0.3",
         "degree": "degree", "burial": "burial", "nbr1_plddt": "nbr1_plddt", "am": "am_mean", "rel_pos": "rel_pos"}


def residue_features(E, v):
    """Attach residue-level features (residue tables) and residue-mean AlphaMissense."""
    am = v.groupby(["gene", "pos"]).am_pathogenicity.mean().rename("am_mean").reset_index() \
        if "am_pathogenicity" in v.columns else None
    parts = []
    for uni, d in E.groupby("uniprot"):
        rf = C._residue_table(uni)
        parts.append(d.merge(rf, on="pos", how="left"))
    F = pd.concat(parts, ignore_index=True)
    if am is not None:
        F = F.merge(am, on=["gene", "pos"], how="left")
    return F


def kw_p(x, y):
    m = ~np.isnan(x)
    groups = [x[m & (y == l)] for l in np.unique(y)]
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) < 2 or np.ptp(np.concatenate(groups)) == 0:
        return np.nan
    return float(stats.kruskal(*groups).pvalue)


def spline_resid(x, pos):
    from patsy import dmatrix
    m = ~np.isnan(x)
    r = np.full(len(x), np.nan)
    if m.sum() >= 12 and len(np.unique(pos[m])) >= 8:
        X = np.asarray(dmatrix("cr(p, df=4)", {"p": pos[m] / pos.max()}, return_type="matrix"))
    else:
        X = np.c_[np.ones(m.sum()), pos[m]]
    b, *_ = np.linalg.lstsq(X, x[m], rcond=None)
    r[m] = x[m] - X @ b
    return r


def secondary_features(F, perms, rng):
    rows, null_counts, rot_counts = [], {f: np.zeros(perms, int) for f in FEATS}, {f: np.zeros(perms, int) for f in FEATS}
    for g, d in F.groupby("gene"):
        d = d.sort_values("pos")
        y = d.label.to_numpy(object)
        pos = d.pos.to_numpy(float)
        for f, col in FEATS.items():
            if col not in d.columns:
                continue
            x = d[col].to_numpy(float)
            p = kw_p(x, y)
            xr = spline_resid(x, pos) if f != "rel_pos" else None
            pr = kw_p(xr, y) if xr is not None else np.nan
            for b in range(perms):
                yp = rng.permutation(y)
                null_counts[f][b] += (kw_p(x, yp) < 0.05)
                if xr is not None:
                    yr = np.roll(y, rng.integers(1, len(y)))
                    rot_counts[f][b] += (kw_p(xr, yr) < 0.05)
            rows.append(dict(gene=g, feature=f, p=p, p_resid=pr))
    T = pd.DataFrame(rows)
    S = []
    for f in FEATS:
        t = T[T.feature == f]
        if not len(t):
            continue
        o, orr = int((t.p < .05).sum()), int((t.p_resid < .05).sum())
        S.append(dict(feature=f, n_genes=len(t), obs_p05=o, null_mean=float(null_counts[f].mean()),
                      perm_p=float((1 + (null_counts[f] >= o).sum()) / (1 + perms)),
                      obs_resid_p05=orr, rot_null_mean=float(rot_counts[f].mean()),
                      rot_perm_p=float((1 + (rot_counts[f] >= orr).sum()) / (1 + perms)) if f != "rel_pos" else np.nan))
    return T, pd.DataFrame(S)


def ad_ar(F, diseases):
    inh = diseases.set_index("disease_unit").inheritance_clean.to_dict()
    F = F.assign(inh=F.label.map(inh))
    rows = []
    for g, d in F.groupby("gene"):
        ad, ar = d[d.inh == "AD"], d[d.inh == "AR"]
        if len(ad) >= MIN_RES and len(ar) >= MIN_RES:
            rec = dict(gene=g, n_ad=len(ad), n_ar=len(ar))
            for f in ["burial", "degree", "plddt", "core_a0.3", "am_mean"]:
                if f in d.columns:
                    rec[f] = float(np.nanmedian(ar[f]) - np.nanmedian(ad[f]))
            rows.append(rec)
    T = pd.DataFrame(rows)
    S = []
    for f in ["burial", "degree", "plddt", "core_a0.3", "am_mean"]:
        if f in T.columns and T[f].notna().sum() >= 5:
            x = T[f].dropna().to_numpy()
            nz = x[x != 0]
            S.append(dict(feature=f, n_genes=len(x), median_AR_minus_AD=float(np.median(x)), n_pos=int((x > 0).sum()),
                          p_wilcoxon=float(stats.wilcoxon(nz).pvalue) if len(nz) >= 5 else np.nan,
                          p_sign=float(stats.binomtest(int((x > 0).sum()), int((x != 0).sum())).pvalue)))
    return T, pd.DataFrame(S)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--calib", type=int, default=200)
    ap.add_argument("--sens-perms", type=int, default=200)
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    v = load_variants()
    diseases = pd.read_parquet(C.OUT / "diseases.parquet")
    log = {}

    R = residue_labels(v)
    E = eligible(R)
    log["residues_labelled"] = int(R.label.notna().sum())
    log["residues_conflict"] = int(R.conflict.sum())
    log["eligible_genes"] = int(E.gene.nunique())
    log["eligible_residues"] = int(len(E))
    E.to_csv(OUT / "eligible_residues.csv", index=False)

    # primary
    G, Dstore = gene_scores(E, k=3, alphas=(0.3,), perms=args.perms, rng=rng, keep_D=True)
    G.to_csv(OUT / "gene_scores_primary.csv", index=False)
    T = primary_tests(G, rng=rng)
    T.to_csv(OUT / "primary_tests.csv", index=False)
    print(T.to_string(index=False))
    sa = []
    for m in ["seq", "3d", "hop", "rwr0.3"]:
        sa.append(dict(metric=m, median_ba=float(G[f"ba_{m}"].median()), median_chance=float(G[f"chance_{m}"].median()),
                       n_genes_p05=int((G[f"p_{m}"] < .05).sum()),
                       fisher_p=float(stats.combine_pvalues(G[f"p_{m}"], method="fisher").pvalue)))
    pd.DataFrame(sa).to_csv(OUT / "signal_at_all.csv", index=False)

    cal = calibration(Dstore, 3, args.calib, rng)
    cal.to_csv(OUT / "calibration.csv", index=False)
    log["calibration"] = cal.groupby("hypothesis").agg(median_of_median_diff=("median_diff", "median"),
                                                       frac_p05=("p", lambda s: float((s < .05).mean()))).to_dict("index")

    # controls
    G[G.gene.isin(POS_CTRL + NULL_CTRL)].to_csv(OUT / "controls.csv", index=False)

    # secondary features and AD/AR
    F = residue_features(E, v)
    Tf, Sf = secondary_features(F, args.sens_perms, rng)
    Tf.to_csv(OUT / "features_by_gene.csv", index=False)
    Sf.to_csv(OUT / "features_summary.csv", index=False)
    Ta, Sa = ad_ar(F, diseases)
    Ta.to_csv(OUT / "ad_ar_by_gene.csv", index=False)
    Sa.to_csv(OUT / "ad_ar_summary.csv", index=False)

    # sensitivity
    sens = []

    def run_sens(tag, EE, k=3, alphas=(0.3,), rwr="rwr0.3"):
        if EE.gene.nunique() < 5:
            sens.append(dict(analysis=tag, note="fewer than 5 genes"))
            return
        Gs, _ = gene_scores(EE, k=k, alphas=alphas, perms=args.sens_perms, rng=rng)
        Ts = primary_tests(Gs, rwr=rwr, rng=rng)
        Ts.insert(0, "analysis", tag)
        sens.extend(Ts.to_dict("records"))

    run_sens("k1", E, k=1)
    run_sens("k5", E, k=5)
    run_sens("alpha0.15", E, alphas=(0.15,), rwr="rwr0.15")
    run_sens("alpha0.5", E, alphas=(0.5,), rwr="rwr0.5")
    run_sens("majority_labels", eligible(residue_labels(v, majority=True)))
    if "stars" in v.columns:
        run_sens("clinvar_ge2_stars", eligible(residue_labels(v[(v.source == "ppkt") | (v.stars >= 2)])))
    run_sens("clinvar_only", eligible(residue_labels(v[v.source.isin(["clinvar", "both"])])))
    gm = pd.read_parquet(C.OUT / "genes.parquet")
    if "mean_plddt" in gm.columns:
        hi = set(gm.loc[gm.mean_plddt >= 70, "gene"])
        run_sens("mean_pLDDT_ge70", E[E.gene.isin(hi)])
    # leakage 1: drop single-site diseases (all residues within 5 aa)
    keep = []
    for (g, l), d in E.groupby(["gene", "label"]):
        if d.pos.max() - d.pos.min() > 5:
            keep.append(d)
    run_sens("drop_single_site_diseases", eligible(pd.concat(keep)))
    S = pd.DataFrame(sens)
    S.to_csv(OUT / "sensitivity.csv", index=False)
    json.dump(log, open(OUT / "log.json", "w"), indent=1, default=str)
    print(json.dumps(log, indent=1, default=str))


if __name__ == "__main__":
    main()
