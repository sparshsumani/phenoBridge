"""Step 4: pre-specified analyses (see ANALYSIS_PLAN.md). Sections are labelled with the plan numbering.
Outputs in out/: primary.csv, primary_sens.csv, omnibus.json, baseline_*.csv, mc2023_reproduction.csv,
missense_table.csv, draft_reproduction.csv, analysis_summary.json
"""
import json

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy import stats

from common import OUT
from stats_utils import bh, logit_firth, logit_ml, lrt_firth, lrt_ml, profile_ci_firth

SEED = 20261103
rng = np.random.default_rng(SEED)
F = pd.read_csv(OUT / "individuals.csv")
C = pd.read_csv(OUT / "clusters.csv")
F = F.merge(C[["individualId", "keep_one_per_family", "keep_one_per_family_or_dup", "keep_excl_likely_dups"]],
            on="individualId", validate="1:1")
meta = json.load(open(OUT / "meta.json"))
H, LAB = meta["hpo"], meta["labels"]
MIN = 10
summary = {}


def eligible(df, min_n=MIN):
    return [h for h in H if (df[h] == 1).sum() >= min_n and (df[h] == 0).sum() >= min_n]


def drop_const(X, names):
    keep = [0] + [j for j in range(1, X.shape[1]) if np.ptp(X[:, j]) > 0]
    return X[:, keep], [names[j] for j in keep]


# ============================================================ 5a primary
pt = F[F.is_point == 1].copy()
sd_semi = pt.log10_semi.std(ddof=0)
pt["z"] = (pt.log10_semi - pt.log10_semi.mean()) / sd_semi
summary["n_point"] = len(pt)
summary["sd_log10_semi"] = float(sd_semi)
TERMS = eligible(pt)
summary["n_terms_primary"] = len(TERMS)


def nested(df, terms, base_cols, test_cols, firth_too=True, spline=None):
    rows = []
    for h in terms:
        d = df[df[h].notna()]
        y = d[h].values.astype(float)
        Xb = [np.ones(len(d))] + [d[c].values.astype(float) for c in base_cols]
        names = ["const"] + list(base_cols)
        if spline is not None:
            S = spline.loc[d.index].values
            Xb += [S[:, j] for j in range(S.shape[1])]
            names += [f"bs{j}" for j in range(S.shape[1])]
        X0, n0 = drop_const(np.column_stack(Xb), names)
        X1 = np.column_stack([X0] + [d[c].values.astype(float) for c in test_cols])
        stat, df_, p, (b1, cov1), conv = lrt_ml(X0, X1, y)
        j = X1.shape[1] - 1
        r = dict(hpo=h, label=LAB[h], n=len(d), n_obs=int(y.sum()), n_exc=int(len(y) - y.sum()),
                 n_nontrunc=int((d.is_trunc == 0).sum()), base_terms="+".join(n0[1:]),
                 ml_converged=bool(conv), LR=stat, p=p, beta=b1[j], se=np.sqrt(cov1[j, j]))
        if firth_too or not conv:
            fs, fp, fb, fcov = lrt_firth(X1, y, [j])
            r.update(firth_LR=fs, firth_p=fp, firth_beta=fb[j], firth_se=np.sqrt(fcov[j, j]))
            if not conv:   # plan: fall back to Firth for this term
                r.update(LR=fs, p=fp, beta=fb[j], se=np.sqrt(fcov[j, j]), used="firth")
        r.setdefault("used", "ml")
        rows.append(r)
    T = pd.DataFrame(rows)
    T["q"] = bh(T.p)
    if "firth_p" in T:
        T["firth_q"] = bh(T.firth_p)
    T["OR"] = np.exp(T.beta)
    T["lo"] = np.exp(T.beta - 1.96 * T.se)
    T["hi"] = np.exp(T.beta + 1.96 * T.se)
    return T


P = nested(pt, TERMS, ["is_trunc", "rel_pos"], ["z"])
P["OR_per10fold"] = np.exp(P.beta / sd_semi)
P["lo_per10fold"] = np.exp((P.beta - 1.96 * P.se) / sd_semi)
P["hi_per10fold"] = np.exp((P.beta + 1.96 * P.se) / sd_semi)
P = P.sort_values("p").reset_index(drop=True)
P.to_csv(OUT / "primary.csv", index=False)
summary["primary"] = dict(n_terms=len(P), min_p=float(P.p.min()), min_q=float(P.q.min()),
                          n_p05=int((P.p < .05).sum()), n_q05=int((P.q < .05).sum()),
                          top=P.head(5)[["label", "n", "OR", "lo", "hi", "p", "q"]].to_dict("records"),
                          n_firth_fallback=int((P.used == "firth").sum()),
                          firth_min_p=float(P.firth_p.min()), firth_min_q=float(P.firth_q.min()),
                          n_firth_p05=int((P.firth_p < .05).sum()))
print(P[["label", "n", "n_obs", "n_exc", "OR", "lo", "hi", "p", "q", "firth_p", "used"]].to_string())

# ---- omnibus permutation (same permutation for all terms; within truncating / non-truncating strata)
B = 5000
idx_by_stratum = [np.where(pt.is_trunc.values == s)[0] for s in (0, 1)]
zs = pt.z.values
term_data = []
for h in TERMS:
    m = pt[h].notna().values
    y = pt[h].values[m].astype(float)
    X0, _ = drop_const(np.column_stack([np.ones(m.sum()), pt.is_trunc.values[m], pt.rel_pos.values[m]]),
                       ["c", "t", "r"])
    _, ll0, _, _ = logit_ml(X0, y)
    term_data.append((m, y, X0, ll0))


def total_stat(zv):
    T = 0.0
    for m, y, X0, ll0 in term_data:
        X1 = np.column_stack([X0, zv[m]])
        _, ll1, _, _ = logit_ml(X1, y)
        T += max(0.0, 2 * (ll1 - ll0))
    return T


T_obs = total_stat(zs)
T_perm = np.empty(B)
for b in range(B):
    zp = zs.copy()
    for ix in idx_by_stratum:
        zp[ix] = zs[rng.permutation(ix)]
    T_perm[b] = total_stat(zp)
omni = dict(T_obs=float(T_obs), B=B, p=float((1 + (T_perm >= T_obs).sum()) / (B + 1)),
            perm_mean=float(T_perm.mean()), perm_q95=float(np.quantile(T_perm, .95)), n_terms=len(TERMS),
            chi2_df_naive_p=float(stats.chi2.sf(T_obs, len(TERMS))))
np.save(OUT / "omnibus_perm.npy", T_perm)
(OUT / "omnibus.json").write_text(json.dumps(omni, indent=1))
summary["omnibus"] = omni
print("omnibus:", omni)

# ---- sensitivity S2-S6
sens = []


def add_sens(name, T):
    sens.append(dict(analysis=name, n_terms=len(T), min_p=float(T.p.min()), min_q=float(T.q.min()),
                     n_p05=int((T.p < .05).sum()), n_q05=int((T.q < .05).sum()),
                     top_term=T.sort_values("p").label.iloc[0],
                     top_OR=float(T.sort_values("p").OR.iloc[0])))
    T.assign(analysis=name).to_csv(OUT / f"sens_{name}.csv", index=False)


add_sens("primary", P)
add_sens("S1_firth", P.assign(p=P.firth_p, q=P.firth_q, OR=np.exp(P.firth_beta)))
spl = pd.DataFrame(np.asarray(patsy.dmatrix("bs(x, df=4, degree=3) - 1", {"x": pt.rel_pos.values})),
                   index=pt.index)
add_sens("S2_spline_position", nested(pt, TERMS, ["is_trunc"], ["z"], firth_too=False, spline=spl))
for comp, col in [("log10_core", "log10_core"), ("backbone", "CommBackboneImpact"),
                  ("disorder", "DisorderImpact"), ("seed_plddt", "seed_pLDDT")]:
    pt[f"z_{comp}"] = (pt[col] - pt[col].mean()) / pt[col].std(ddof=0)
    add_sens(f"S3_{comp}", nested(pt, TERMS, ["is_trunc", "rel_pos"], [f"z_{comp}"], firth_too=False))
for keep in ("keep_one_per_family", "keep_one_per_family_or_dup"):
    sub = pt[pt[keep] == 1]
    add_sens(f"S4_{keep}", nested(sub, eligible(sub), ["is_trunc", "rel_pos"], ["z"], firth_too=False))
tr = pt[pt.is_trunc == 1]
add_sens("S5_truncating_only", nested(tr, eligible(tr), ["rel_pos"], ["z"], firth_too=False))
S = pd.DataFrame(sens)
S.to_csv(OUT / "primary_sens.csv", index=False)
summary["sensitivity"] = S.to_dict("records")
print(S.to_string())


# ---- S6: draft model reproduction (y ~ log10 semi + trunc), terms >=15/>=15 as in the draft
def draft_model(df, trunc_col, na_zero):
    rows = []
    for h in H:
        y = df[h].fillna(0.0) if na_zero else df[h]
        m = y.isin([0, 1])
        yy = y[m].astype(int)
        if yy.sum() < 15 or len(yy) - yy.sum() < 15:
            continue
        X = sm.add_constant(df.loc[m, ["log10_semi", trunc_col]].astype(float), has_constant="add")
        r = sm.Logit(yy, X).fit(disp=0)
        rows.append(dict(hpo=h, label=LAB[h], n=len(yy), OR=np.exp(r.params.log10_semi), p=r.pvalues.log10_semi))
    T = pd.DataFrame(rows)
    T["q"] = bh(T.p)
    return T


dr = []
for tc in ("regex_trunc", "is_trunc"):
    for na0 in (True, False):
        T = draft_model(pt, tc, na0)
        hits = T[T.p < .05].sort_values("p")
        dr.append(dict(trunc_flag=tc, encoding="A (NA=0)" if na0 else "B (explicit)", n_terms=len(T),
                       n_p05=len(hits), min_q=float(T.q.min()),
                       hits="; ".join(f"{l} OR={o:.2f} p={p:.3f}" for l, o, p in zip(hits.label, hits.OR, hits.p))))
DR = pd.DataFrame(dr)
DR.to_csv(OUT / "draft_reproduction.csv", index=False)
summary["draft_reproduction"] = dr
print(DR.to_string())


# ============================================================ 5b baseline associations (Firth)
def firth_assoc(df, terms, expo, covars=(), binary=True, scale=1.0):
    rows = []
    for h in terms:
        d = df[df[h].notna() & df[expo].notna()]
        y = d[h].values.astype(float)
        x = d[expo].values.astype(float) * scale
        X = np.column_stack([np.ones(len(d)), x] + [d[c].values.astype(float) for c in covars])
        stat, p, b, cov = lrt_firth(X, y, [1])
        _, lo, hi = profile_ci_firth(X, y, 1)
        r = dict(hpo=h, label=LAB[h], n=len(d), OR=np.exp(b[1]), lo=np.exp(lo), hi=np.exp(hi), p=p)
        if binary:
            a = int(((x == 1) & (y == 1)).sum()); bb = int(((x == 1) & (y == 0)).sum())
            c = int(((x == 0) & (y == 1)).sum()); dd = int(((x == 0) & (y == 0)).sum())
            r.update(exp_pos=a, exp_n=a + bb, unexp_pos=c, unexp_n=c + dd,
                     prev_exp=a / max(1, a + bb), prev_unexp=c / max(1, c + dd),
                     fisher_p=stats.fisher_exact([[a, bb], [c, dd]])[1])
        rows.append(r)
    T = pd.DataFrame(rows)
    T["q"] = bh(T.p)
    return T.sort_values("p").reset_index(drop=True)


base = {}
allF = F.copy()
base["B1_SV_vs_point"] = firth_assoc(allF, eligible(allF), "is_sv")
dele = allF[(allF.is_point == 1) | (allF.sv_type == "deletion")]
base["B1s_SVdeletion_vs_point"] = firth_assoc(dele, eligible(dele), "is_sv")
base["B2_nontrunc_vs_trunc"] = firth_assoc(pt, TERMS, "is_nontrunc_point")
ptc = pt[pt.is_ptc == 1].copy()
base["B3i_PTC_position_per10pct"] = firth_assoc(ptc, eligible(ptc), "ptc_rel", binary=False, scale=10.0)
base["B3ii_PTC_loses_Cterm_region"] = firth_assoc(ptc, eligible(ptc), "loses_deg_region")
bsum = []
for k, T in base.items():
    T.to_csv(OUT / f"baseline_{k}.csv", index=False)
    top = T.iloc[0]
    bsum.append(dict(family=k, n_terms=len(T), n_p05=int((T.p < .05).sum()), n_q05=int((T.q < .05).sum()),
                     min_p=float(T.p.min()), min_q=float(T.q.min()), top=top.label, top_OR=float(top.OR),
                     top_lo=float(top.lo), top_hi=float(top.hi)))
    print(k); print(T.head(8)[["label", "n", "OR", "lo", "hi", "p", "q"] +
                            (["prev_exp", "prev_unexp"] if "prev_exp" in T else [])].to_string())
summary["baseline"] = bsum
summary["n_ptc"] = len(ptc)
summary["n_ptc_keep_deg"] = int((ptc.loses_deg_region == 0).sum())


# ============================================================ B4 MC2023 reproduction
def composite(df, terms):
    v = df[terms]
    pos = (v == 1).any(axis=1)
    neg = ~pos & (v == 0).any(axis=1)
    return np.where(pos, 1.0, np.where(neg, 0.0, np.nan))


F["NDD_composite"] = composite(F, ["HP:0001249", "HP:0007018", "HP:0000729"])
COMORB = ["HP:0000407", "HP:0000405", "HP:0000410", "HP:0000403", "HP:0000028", "HP:0006695", "HP:0001655",
          "HP:0001629", "HP:0001631", "HP:0001642", "HP:0001647", "HP:0011968", "HP:0002373", "HP:0007359",
          "HP:0002197", "HP:0002121"]
F["comorbidity_composite"] = composite(F, COMORB)
LAB["NDD_composite"] = "ID/ADHD/ASD (any)"
LAB["comorbidity_composite"] = "Any comorbidity"
pt2 = F[F.is_point == 1]
rep = []


def rep_row(claim, reported, T):
    r = T.iloc[0].to_dict()
    r.update(claim=claim, reported=reported)
    rep.append(r)


rep_row("a. Triangular face: sequence variant vs CNV", "71.1% vs 42.5%",
        firth_assoc(F.assign(is_seq=F.is_point), ["HP:0000325"], "is_seq"))
rep_row("b. Short stature: exon 9 vs outside exon 9 (point variants)", "62.5% vs 27.8%",
        firth_assoc(pt2, ["HP:0004322"], "exon9"))
rep_row("c. ID/ADHD/ASD: c.1903_1907del vs other point variants", "70.4% vs 89.4%",
        firth_assoc(pt2, ["NDD_composite"], "hotspot_1903_1907del"))
rep_row("c'. ID alone: c.1903_1907del vs other point variants", "(not reported)",
        firth_assoc(pt2, ["HP:0001249"], "hotspot_1903_1907del"))
sdel = F[F.sv_type == "deletion"].copy()
sdel["log2_size"] = np.log2(sdel.sv_size_kb)
for h, claim in [("HP:0001572", "d. Macrodontia increases with deletion size"),
                 ("comorbidity_composite", "d. Comorbidities increase with deletion size"),
                 ("HP:0001155", "d. Hand anomalies increase with SMALLER deletions")]:
    rep_row(claim + " (OR per doubling)", "associated", firth_assoc(sdel, [h], "log2_size", binary=False))
REP = pd.DataFrame(rep)
REP.to_csv(OUT / "mc2023_reproduction.csv", index=False)
summary["mc2023"] = REP[["claim", "reported", "n", "OR", "lo", "hi", "p"] +
                        [c for c in ("prev_exp", "prev_unexp", "exp_n", "unexp_n", "fisher_p") if c in REP]].to_dict("records")
print(REP[["claim", "reported", "n", "OR", "lo", "hi", "p", "prev_exp", "prev_unexp", "exp_n", "unexp_n", "fisher_p"]].to_string())

# ============================================================ 5c missense
MS = F[F.vclass == "missense"][["individualId", "hgvs_p", "aa_pos", "vep_exon", "seed_pLDDT", "SemiCoreImpact",
                               "CoreImpact", "am_score", "am_class", "n_annotated"]]
MS.to_csv(OUT / "missense_table.csv", index=False)
summary["n_missense"] = len(MS)
summary["n_inframe"] = int((F.vclass == "inframe").sum())

# descriptive: NMD and C-terminal
summary["nmd"] = dict(n_ptc=int(F.is_ptc.sum()),
                      escape50=int(F.nmd_escape_50nt.sum()), escape55=int(F.nmd_escape_55nt.sum()),
                      long_exon=int(F.ptc_in_long_exon_407.sum()),
                      ptc_exon_counts=F.ptc_exon.value_counts().to_dict(),
                      loses_deg_region=int(F.loses_deg_region.sum()),
                      loses_all_CT_conf=int(F.loses_all_CT_conf.sum()))
(OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=1, default=str))
print("done")
