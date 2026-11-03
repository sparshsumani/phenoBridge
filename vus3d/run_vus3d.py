"""VUS-3D evaluation, implementing vus3d/PLAN.md.

    python3 run_vus3d.py [--boot 2000]

Every score for a query residue is computed from TRAIN (Jan 2023) labels only, excluding
labels at the query residue itself (primary) or including them (sensitivity).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss, average_precision_score

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402  (vus3d/common.py: DATA paths, contact_graph, ca_coords)

OUT = HERE / "out"
SEED = 20261106
EPS = 1e-4
ALPHAS = (0.15, 0.3, 0.5)
SIGMAS = (4.0, 6.0, 8.0, 12.0)
LAMBDAS = (2.0, 5.0, 10.0, 20.0)
DEFAULT = {"rwr": 0.3, "r3d": 6.0, "w1d": 5.0}
VCEP_PM1 = ["PTEN", "TP53", "HRAS", "KRAS", "NRAS", "RAF1", "PTPN11", "SOS1", "SOS2", "BRAF", "MAP2K1", "MAP2K2",
            "RIT1", "SHOC2", "CBL", "LZTR1", "RRAS2", "MYH7", "MYBPC3", "CDH1", "RUNX1"]


# ------------------------------------------------------------------ data
def load():
    v = pd.read_parquet(C.OUT / "variants.parquet")
    v = v[v.af_pos.notna() & v.am_pathogenicity.notna() & v.label.isin(["P", "B"])].copy()
    v["pos"] = v.af_pos.astype(int)
    v["y"] = (v.label == "P").astype(int)
    am = v.am_pathogenicity.clip(1e-4, 1 - 1e-4)
    v["am_logit"] = np.log(am / (1 - am))
    return v


def train_residue_labels(tr):
    """Residue -> 'P'/'B' if all TRAIN labels at it agree; mixed residues dropped."""
    g = tr.groupby(["uniprot", "pos"]).y.agg(["min", "max"]).reset_index()
    g = g[g["min"] == g["max"]]
    g["lab"] = np.where(g["max"] == 1, "P", "B")
    return g[["uniprot", "pos", "lab"]]


# ------------------------------------------------------------------ kernels
class Protein:
    def __init__(self, uni):
        self.uni = uni
        self.pos, self.W, self.A = C.contact_graph(uni)
        _, self.xyz, _ = C.ca_coords(uni)
        self.idx = {p: k for k, p in enumerate(self.pos)}
        self._lu = {}

    def lu(self, a):
        if a not in self._lu:
            n = len(self.pos)
            self._lu[a] = splu((sp.identity(n, format="csc") - (1 - a) * self.W).tocsc())
        return self._lu[a]

    def rwr_mass(self, q_idx, src_mask, a):
        """m(q) = sum_{j in src} p_q(j) for all query indices q, where
        p_q = a (I - (1-a) W)^-1 e_q. Uses one transposed solve per source set, plus the
        diagonal p_q(q) for the leave-own-residue-out correction."""
        lu = self.lu(a)
        x = a * lu.solve(src_mask.astype(float), trans="T")   # x[q] = 1_src^T M e_q
        diag = np.empty(len(q_idx))
        for t, q in enumerate(q_idx):
            e = np.zeros(len(self.pos)); e[q] = 1.0
            diag[t] = a * lu.solve(e)[q]
        return x[q_idx], diag

    def r3d_w(self, q_idx, s):
        d2 = ((self.xyz[q_idx][:, None, :] - self.xyz[None, :, :]) ** 2).sum(-1)
        w = np.exp(-d2 / (2 * s * s))
        return w / w.sum(1, keepdims=True)

    def w1d(self, q_idx, lam):
        d = np.abs(self.pos[q_idx][:, None] - self.pos[None, :]).astype(float)
        w = np.exp(-d / lam)
        return w / w.sum(1, keepdims=True)


def protein_features(prot, labs, queries, include_self=False):
    """Label-mass features for query residues (array of residue numbers) in one protein.
    Returns DataFrame indexed by query residue with columns for every kernel/parameter."""
    n = len(prot.pos)
    P = np.zeros(n, bool); B = np.zeros(n, bool)
    for p, l in zip(labs.pos, labs.lab):
        k = prot.idx.get(int(p))
        if k is not None:
            (P if l == "P" else B)[k] = True
    q_idx = np.array([prot.idx[int(q)] for q in queries])
    selfP, selfB = P[q_idx].astype(float), B[q_idx].astype(float)
    out = {}
    for a in ALPHAS:
        mP, dg = prot.rwr_mass(q_idx, P, a)
        mB, _ = prot.rwr_mass(q_idx, B, a)
        if not include_self:
            mP = mP - dg * selfP; mB = mB - dg * selfB
        out[f"rwr{a}_P"], out[f"rwr{a}_B"] = mP, mB
    for name, params, fn in (("r3d", SIGMAS, prot.r3d_w), ("w1d", LAMBDAS, prot.w1d)):
        for s in params:
            w = fn(q_idx, s)
            if not include_self:
                w[np.arange(len(q_idx)), q_idx] = 0.0
            out[f"{name}{s}_P"], out[f"{name}{s}_B"] = w @ P, w @ B
    F = pd.DataFrame(out, index=queries)
    return F.clip(lower=0)


def build_features(v, include_self=False):
    tr = v[v.set == "train"]
    labs = train_residue_labels(tr)
    rows = []
    for uni, d in v.groupby("uniprot"):
        lab_u = labs[labs.uniprot == uni]
        if len(lab_u) == 0:
            continue
        try:
            prot = Protein(uni)
        except Exception as e:  # missing graph
            print("skip", uni, e)
            continue
        qs = np.array(sorted(set(d.pos) & set(prot.pos)))
        F = protein_features(prot, lab_u, qs, include_self)
        F["uniprot"] = uni
        F["pos"] = F.index
        # gene prior, leave-own-residue-out for every query residue
        nP, nB = int((lab_u.lab == "P").sum()), int((lab_u.lab == "B").sum())
        own = lab_u.set_index("pos").lab.reindex(qs)
        oP = (own == "P").astype(int).values if not include_self else 0
        oB = (own == "B").astype(int).values if not include_self else 0
        F["n_lab_res"] = nP + nB - (oP + oB if not include_self else 0)
        pr = (nP - oP + 1) / (nP + nB - oP - oB + 2)
        F["gene_prior"] = np.log(pr / (1 - pr))
        rows.append(F.reset_index(drop=True))
    F = pd.concat(rows, ignore_index=True)
    return v.merge(F, on=["uniprot", "pos"], how="inner")


# ------------------------------------------------------------------ models
def design(df, model, par=None):
    X = [df.am_logit.values]
    if model != "M0":
        X.append(df.gene_prior.values)
    if model.startswith("M2"):
        k = {"M2-RWR": f"rwr{par}", "M2-3D": f"r3d{par}", "M2-1D": f"w1d{par}"}[model]
        X += [np.log(df[f"{k}_P"].values + EPS), np.log(df[f"{k}_B"].values + EPS)]
    return np.column_stack(X)


def fit(df, model, par=None):
    m = LogisticRegression(C=1e6, max_iter=2000)
    m.fit(design(df, model, par), df.y.values)
    return m


def grouped_cv_auc(tr, model, par, k=5, seed=SEED):
    genes = tr.gene.unique()
    rng = np.random.default_rng(seed)
    fold = dict(zip(genes, rng.permutation(len(genes)) % k))
    f = tr.gene.map(fold).values
    pred = np.empty(len(tr))
    for i in range(k):
        te = f == i
        m = fit(tr[~te], model, par)
        pred[te] = m.predict_proba(design(tr[te], model, par))[:, 1]
    return roc_auc_score(tr.y, pred)


def tune(tr):
    best, cv = {}, []
    for model, key, grid in (("M2-RWR", "rwr", ALPHAS), ("M2-3D", "r3d", SIGMAS), ("M2-1D", "w1d", LAMBDAS)):
        aucs = {p: grouped_cv_auc(tr, model, p) for p in grid}
        best[key] = max(aucs, key=aucs.get)
        cv += [dict(model=model, param=p, cv_auc=a) for p, a in aucs.items()]
    for model in ("M0", "M1"):
        cv.append(dict(model=model, param=None, cv_auc=grouped_cv_auc(tr, model, None)))
    return best, pd.DataFrame(cv)


# ------------------------------------------------------------------ evaluation
def strat_auc(y, s, g):
    """Gene-stratified AUROC: concordance over P-B pairs within the same gene."""
    num = den = 0.0
    for gg in np.unique(g):
        m = g == gg
        yy, ss = y[m], s[m]
        if yy.min() == yy.max():
            continue
        p, b = ss[yy == 1], ss[yy == 0]
        r = stats_rank(np.r_[p, b])
        nP, nB = len(p), len(b)
        u = r[:nP].sum() - nP * (nP + 1) / 2
        num += u; den += nP * nB
    return num / den if den else np.nan


def stats_rank(x):
    from scipy.stats import rankdata
    return rankdata(x)


def boot_gene(te, preds, pairs, nboot, rng, metric="strat"):
    genes = te.gene.unique()
    gidx = {g: np.where(te.gene.values == g)[0] for g in genes}
    out = {p: [] for p in pairs}
    y = te.y.values
    gl = te.gene.values

    def met(ii, s, gg):
        if metric == "strat":
            return strat_auc(y[ii], s[ii], gg)
        return roc_auc_score(y[ii], s[ii])
    for _ in range(nboot):
        pick = rng.choice(genes, len(genes))
        ii = np.concatenate([gidx[g] for g in pick])
        gg = np.concatenate([np.full(len(gidx[g]), k) for k, g in enumerate(pick)])  # resampled genes stay distinct
        if y[ii].min() == y[ii].max():
            continue
        for a, b in pairs:
            out[(a, b)].append(met(ii, preds[a], gg) - met(ii, preds[b], gg))
    res = []
    allii = np.arange(len(y))
    for (a, b), d in out.items():
        d = np.array(d)
        obs = met(allii, preds[a], gl) - met(allii, preds[b], gl)
        p = 2 * min((d <= 0).mean(), (d >= 0).mean())
        res.append(dict(metric=metric, contrast=f"{a} - {b}", delta_auc=obs, ci_lo=np.percentile(d, 2.5),
                        ci_hi=np.percentile(d, 97.5), ci90_lo=np.percentile(d, 5), ci90_hi=np.percentile(d, 95),
                        p_boot=max(p, 1 / len(d))))
    return pd.DataFrame(res)


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0
    for r, i in enumerate(o):
        run = max(run, (m - r) * p[i]); adj[i] = min(1, run)
    return adj


def within_gene_auc(te, pred, min_each=3):
    a = []
    for g, d in te.groupby("gene"):
        if d.y.sum() >= min_each and (1 - d.y).sum() >= min_each:
            a.append(roc_auc_score(d.y, pred[d.index]))
    return (float(np.mean(a)) if a else np.nan), len(a)


def evaluate(v, F, best, nboot, rng, tag, min_lab=5, test_filter=None):
    tr = F[F.set == "train"]
    te = F[(F.set == "test") & (F.n_lab_res >= min_lab)]
    if test_filter is not None:
        te = te[test_filter(te)]
    te = te.reset_index(drop=True)
    models = {"M0": None, "M1": None, "M2-1D": best["w1d"], "M2-3D": best["r3d"], "M2-RWR": best["rwr"]}
    preds, rows = {}, []
    for m, par in models.items():
        fm = fit(tr, m, par)
        preds[m] = fm.predict_proba(design(te, m, par))[:, 1]
        wg, ng = within_gene_auc(te, pd.Series(preds[m]))
        rows.append(dict(analysis=tag, model=m, param=par, n_test=len(te), n_P=int(te.y.sum()),
                         n_genes=te.gene.nunique(), strat_auroc=strat_auc(te.y.values, preds[m], te.gene.values),
                         auroc=roc_auc_score(te.y, preds[m]),
                         auprc=average_precision_score(te.y, preds[m]), brier=brier_score_loss(te.y, preds[m]),
                         within_gene_auroc=wg, n_genes_within=ng))
    # H3 compares the 3D model chosen on TRAIN CV
    cvbest = "M2-RWR" if best["_cv_rwr"] >= best["_cv_r3d"] else "M2-3D"
    pairs = [("M2-RWR", "M1"), ("M2-RWR", "M2-3D"), (cvbest, "M2-1D"), ("M2-3D", "M1"), ("M2-1D", "M1"), ("M1", "M0")]
    B = boot_gene(te, preds, pairs, nboot, rng, metric="strat")
    B["hypothesis"] = ["H1", "H2", "H3", "", "", ""]
    m = B.hypothesis != ""
    B.loc[m, "p_holm"] = holm(B.loc[m, "p_boot"].values)
    Bp = boot_gene(te, preds, pairs, max(200, nboot // 5), rng, metric="pooled")
    Bp["hypothesis"] = ["H1-pooled", "H2-pooled", "H3-pooled", "", "", ""]
    B = pd.concat([B, Bp], ignore_index=True)
    B.insert(0, "analysis", tag)
    return pd.DataFrame(rows), B, te, preds


def calibration_bins(F, best, te, preds):
    """PM1-like evidence: LR of P vs B by bins of the propagation component (M2-RWR log-odds
    minus M1 log-odds), estimated on TRAIN (out-of-fold) and checked on TEST."""
    tr = F[F.set == "train"].reset_index(drop=True)
    par = best["rwr"]
    genes = tr.gene.unique(); rng = np.random.default_rng(SEED)
    fold = dict(zip(genes, rng.permutation(len(genes)) % 5)); f = tr.gene.map(fold).values
    comp = np.empty(len(tr))
    for i in range(5):
        t = f == i
        m2, m1 = fit(tr[~t], "M2-RWR", par), fit(tr[~t], "M1")
        comp[t] = m2.decision_function(design(tr[t], "M2-RWR", par)) - m1.decision_function(design(tr[t], "M1"))
    m2, m1 = fit(tr, "M2-RWR", par), fit(tr, "M1")
    comp_te = m2.decision_function(design(te, "M2-RWR", par)) - m1.decision_function(design(te, "M1"))
    edges = [-np.inf, -2, -1, -0.5, 0.5, 1, 2, np.inf]
    rows = []
    priorP_tr, priorP_te = tr.y.mean(), te.y.mean()
    for lo, hi in zip(edges[:-1], edges[1:]):
        def lr(y, c, prior):
            s = (c > lo) & (c <= hi)
            nP, nB = int(y[s].sum()), int((1 - y[s]).sum())
            if nP + nB == 0:
                return nP, nB, np.nan
            return nP, nB, ((nP + .5) / (y.sum() + .5)) / ((nB + .5) / ((1 - y).sum() + .5))
        a = lr(tr.y.values, comp, priorP_tr); b = lr(te.y.values, comp_te, priorP_te)
        rows.append(dict(bin=f"({lo},{hi}]", train_P=a[0], train_B=a[1], train_LR=a[2],
                         test_P=b[0], test_B=b[1], test_LR=b[2]))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=2000)
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    v = load()
    F = build_features(v)
    F.to_parquet(OUT / "features.parquet")
    tr = F[F.set == "train"]
    best, cv = tune(tr[tr.n_lab_res >= 5])
    cv.to_csv(OUT / "train_cv.csv", index=False)
    best["_cv_rwr"] = cv[(cv.model == "M2-RWR") & (cv.param == best["rwr"])].cv_auc.iloc[0]
    best["_cv_r3d"] = cv[(cv.model == "M2-3D") & (cv.param == best["r3d"])].cv_auc.iloc[0]
    json.dump(best, open(OUT / "tuned_params.json", "w"), indent=1)

    allm, allb = [], []
    M, B, te, preds = evaluate(v, F, best, args.boot, rng, "primary")
    allm.append(M); allb.append(B)
    print(M.to_string(index=False)); print(B.to_string(index=False))
    calibration_bins(F, best, te, preds).to_csv(OUT / "calibration_bins.csv", index=False)
    pd.DataFrame({"gene": te.gene, "y": te.y, **{k: p for k, p in preds.items()}}).to_csv(OUT / "test_predictions.csv.gz", index=False)

    # gain by data density
    dens = []
    for lo, hi in ((5, 9), (10, 19), (20, 49), (50, 10**6)):
        s = (te.n_lab_res >= lo) & (te.n_lab_res <= hi)
        if s.sum() > 20 and te.y[s].nunique() == 2:
            dens.append(dict(bin=f"{lo}-{hi}", n=int(s.sum()), auc_M1=roc_auc_score(te.y[s], preds["M1"][s]),
                             auc_RWR=roc_auc_score(te.y[s], preds["M2-RWR"][s])))
    pd.DataFrame(dens).to_csv(OUT / "gain_by_density.csv", index=False)

    sb = max(200, args.boot // 5)
    fixed = dict(best); fixed.update(DEFAULT)
    for tag, kw in (("ge2_stars", dict(test_filter=lambda t: t.stars_2026 >= 2)),
                    ("prior_vus_only", dict(test_filter=lambda t: t.status_2023 == "VUS")),
                    ("genes_ge20_labels", dict(min_lab=20)),
                    ("strict_temporal", dict(test_filter=lambda t: t.last_evaluated_2026_after_cutoff.fillna(False).astype(bool))),
                    ("evaluated_after_2024_09", dict(test_filter=lambda t: pd.to_datetime(t.last_evaluated_2026, errors="coerce") >= "2024-09-01")),
                    ("evaluated_before_2024_09", dict(test_filter=lambda t: pd.to_datetime(t.last_evaluated_2026, errors="coerce") < "2024-09-01")),
                    ("non_vcep_pm1_genes", dict(test_filter=lambda t: ~t.gene.isin(VCEP_PM1))),
                    ("vcep_pm1_genes", dict(test_filter=lambda t: t.gene.isin(VCEP_PM1)))):
        M, B, _, _ = evaluate(v, F, best, sb, rng, tag, **kw); allm.append(M); allb.append(B)
    M, B, _, _ = evaluate(v, F, fixed, sb, rng, "fixed_params"); allm.append(M); allb.append(B)
    Fs = build_features(v, include_self=True)
    M, B, _, _ = evaluate(v, Fs, best, sb, rng, "include_same_residue"); allm.append(M); allb.append(B)
    pd.concat(allm).to_csv(OUT / "model_metrics.csv", index=False)
    pd.concat(allb).to_csv(OUT / "contrasts.csv", index=False)


if __name__ == "__main__":
    main()
