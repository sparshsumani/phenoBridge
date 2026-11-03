"""Step 5 (plan 5d): simulation power analysis for the primary per-term LRT, using each term's real
covariate rows (Trunc, RelPos, z(log10 SemiCore)) and its fitted baseline model.
Output: out/power_grid.csv, out/power_mde.csv, out/power_summary.json"""
import json

import numpy as np
import pandas as pd
from scipy import stats

from common import OUT
from stats_utils import logit_ml

SEED = 20261103
rng = np.random.default_rng(SEED)
NSIM = 1000
GRID = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.75, 2.0, 2.25, 2.5, 3.0, 4.0])
ALPHAS = (0.05, 0.005)

F = pd.read_csv(OUT / "individuals.csv")
P = pd.read_csv(OUT / "primary.csv")
meta = json.load(open(OUT / "meta.json"))
pt = F[F.is_point == 1].copy()
sd = pt.log10_semi.std(ddof=0)
pt["z"] = (pt.log10_semi - pt.log10_semi.mean()) / sd


def batch_ll(X, Y, iters=30):
    """Vectorised Newton logistic fits of many outcome vectors Y (S,n) on one design X (n,p)."""
    S, n = Y.shape
    p = X.shape[1]
    B = np.zeros((S, p))
    for _ in range(iters):
        eta = np.clip(B @ X.T, -30, 30)
        mu = 1 / (1 + np.exp(-eta))
        W = np.clip(mu * (1 - mu), 1e-9, None)
        g = (Y - mu) @ X
        Hm = np.einsum("sn,ni,nj->sij", W, X, X) + 1e-8 * np.eye(p)
        step = np.linalg.solve(Hm, g[..., None])[..., 0]
        step = np.clip(step, -5, 5)
        B = B + step
    eta = np.clip(B @ X.T, -30, 30)
    return np.sum(Y * eta - np.logaddexp(0, eta), axis=1)


rows = []
for h in P.hpo:
    d = pt[pt[h].notna()]
    y = d[h].values.astype(float)
    cols = [np.ones(len(d))]
    for c in ("is_trunc", "rel_pos"):
        if np.ptp(d[c].values) > 0:
            cols.append(d[c].values.astype(float))
    X0 = np.column_stack(cols)
    z = d.z.values - d.z.values.mean()
    X1 = np.column_stack([X0, z])
    b0, _, _, _ = logit_ml(X0, y)
    lp0 = X0 @ b0
    for OR in GRID:
        mu = 1 / (1 + np.exp(-(lp0 + np.log(OR) * z)))
        Y = (rng.random((NSIM, len(d))) < mu).astype(float)
        lr = np.maximum(0, 2 * (batch_ll(X1, Y) - batch_ll(X0, Y)))
        pv = stats.chi2.sf(lr, 1)
        rows.append(dict(hpo=h, label=meta["labels"][h], n=len(d), n_obs=int(y.sum()), OR_per_SD=OR,
                         **{f"power_{a}": float((pv < a).mean()) for a in ALPHAS}))
G = pd.DataFrame(rows)
G.to_csv(OUT / "power_grid.csv", index=False)


def mde(g, col, target=0.8):
    x, yv = g.OR_per_SD.values, g[col].values
    above = np.where(yv >= target)[0]
    if len(above) == 0:
        return np.nan
    i = above[0]
    if i == 0:
        return x[0]
    return float(x[i - 1] + (target - yv[i - 1]) * (x[i] - x[i - 1]) / (yv[i] - yv[i - 1]))


md = []
for h, g in G.groupby("hpo", sort=False):
    r = dict(hpo=h, label=g.label.iloc[0], n=int(g.n.iloc[0]), n_obs=int(g.n_obs.iloc[0]),
             typeI_005=float(g.loc[g.OR_per_SD == 1, "power_0.05"].iloc[0]))
    for a in ALPHAS:
        m = mde(g, f"power_{a}")
        r[f"MDE_OR_perSD_{a}"] = m
        r[f"MDE_OR_per10fold_{a}"] = m ** (1 / sd) if np.isfinite(m) else np.nan
    md.append(r)
MD = pd.DataFrame(md).sort_values("n", ascending=False)
MD.to_csv(OUT / "power_mde.csv", index=False)
def med_inf(x):   # terms whose MDE exceeds the grid maximum count as +inf (not dropped)
    v = np.where(np.isnan(x), np.inf, x)
    m = float(np.median(v))
    return m if np.isfinite(m) else f">{GRID[-1]}"


summ = dict(nsim=NSIM, grid=GRID.tolist(), sd_log10_semi=float(sd), n_terms=len(MD), median_n=float(MD.n.median()),
            mean_typeI_at_005=float(MD.typeI_005.mean()), mean_typeI_note="empirical rejection rate at OR=1, alpha=0.05")
for a in ALPHAS:
    col = MD[f"MDE_OR_perSD_{a}"]
    summ[f"median_MDE_perSD_{a}"] = med_inf(col.values)
    summ[f"min_MDE_perSD_{a}"] = float(col.min())
    summ[f"n_terms_MDE_gt_{GRID[-1]}_{a}"] = int(col.isna().sum())
    summ[f"n_terms_MDE_le_1.5_{a}"] = int((col <= 1.5).sum())
    summ[f"n_terms_MDE_le_2_{a}"] = int((col <= 2).sum())
    m = med_inf(col.values)
    summ[f"median_MDE_per10fold_{a}"] = m ** (1 / sd) if isinstance(m, float) else m
    # power at fixed effect sizes, averaged over terms
    for OR in (1.5, 2.0):
        summ[f"mean_power_OR{OR}_{a}"] = float(G.loc[np.isclose(G.OR_per_SD, OR), f"power_{a}"].mean())
(OUT / "power_summary.json").write_text(json.dumps(summ, indent=1))
print(MD.to_string())
print(json.dumps(summ, indent=1))
