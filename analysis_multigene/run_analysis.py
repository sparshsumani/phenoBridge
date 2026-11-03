#!/usr/bin/env python3
"""Pre-specified multi-gene analysis (see PLAN.md).

Reads the derived phenopacket-store tables built by multigene/ and writes
results to analysis_multigene/out/.

Usage: python3 run_analysis.py [--perms 1000] [--data /home/user/data_ext/derived]
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import fisher_exact, hypergeom, norm, rankdata

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
MIN_OBS = 5          # >= 5 observed and >= 5 excluded per (gene, term) test
RNG_SEED = 20261103


# ---------------------------------------------------------------- helpers
def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def storey_pi1(p, lam=0.5):
    p = np.asarray(p, float)
    if len(p) == 0:
        return np.nan
    pi0 = min(1.0, np.mean(p > lam) / (1 - lam))
    return 1.0 - pi0


def mwu_p_matrix(X, y):
    """Two-sided Mann-Whitney U p-values (normal approx., tie-corrected).

    X: (B, n) feature values (each row one permutation); y: (n,) 0/1 groups.
    Returns p (B,) and rank-biserial r (B,) for group y==1 vs y==0.
    """
    X = np.atleast_2d(X)
    n = X.shape[1]
    n1 = int(y.sum())
    n0 = n - n1
    R = rankdata(X, axis=1)
    U = R[:, y == 1].sum(axis=1) - n1 * (n1 + 1) / 2.0
    mu = n1 * n0 / 2.0
    # tie correction
    ties = np.zeros(X.shape[0])
    for b in range(X.shape[0]):
        _, c = np.unique(X[b], return_counts=True)
        ties[b] = (c ** 3 - c).sum()
    var = n1 * n0 / 12.0 * ((n + 1) - ties / (n * (n - 1)))
    var = np.where(var <= 0, np.nan, var)
    z = (U - mu) / np.sqrt(var)
    p = 2 * norm.sf(np.abs(z))
    p = np.where(np.isnan(p), 1.0, p)
    r = 2 * U / (n1 * n0) - 1
    return p, r


def fisher_p_lookup(n, n_case, k_marg):
    """Two-sided Fisher p for every possible 'a' given margins (hypergeometric).

    n = test size, n_case = #observed (y==1), k_marg = #feature==1.
    Returns array p[a] for a in 0..min(n_case,k_marg).
    """
    lo = max(0, k_marg + n_case - n)
    hi = min(n_case, k_marg)
    a = np.arange(lo, hi + 1)
    pmf = hypergeom.pmf(a, n, k_marg, n_case)
    out = np.ones(hi + 1)
    for i, ai in enumerate(a):
        out[ai] = min(1.0, pmf[pmf <= pmf[i] * (1 + 1e-7)].sum())
    return out


def neighbourhood_plddt(pdb_gz, positions):
    """Mean pLDDT over 1-hop and 2-hop Calpha contact neighbourhoods (8 A)."""
    pos, xyz, pl = [], [], []
    with gzip.open(pdb_gz, "rt") as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                pos.append(int(line[22:26]))
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                pl.append(float(line[60:66]))
    pos = np.array(pos)
    xyz = np.array(xyz)
    pl = np.array(pl)
    idx = {p: i for i, p in enumerate(pos)}
    tree = cKDTree(xyz)
    nbrs = tree.query_ball_point(xyz, r=8.0)
    out = {}
    for p in positions:
        if p not in idx:
            continue
        i = idx[p]
        h1 = set(nbrs[i])
        h2 = set(h1)
        for j in h1:
            h2.update(nbrs[j])
        out[p] = (pl[list(h1)].mean(), pl[list(h2)].mean())
    return out


# ---------------------------------------------------------------- data
def load(data, alpha_suffix=""):
    d = Path(data)
    genes = pd.read_parquet(d / "genes.parquet")
    genes = genes[genes.selected].copy()
    ind = pd.read_parquet(d / "individuals.parquet")
    var = pd.read_parquet(d / "variants.parquet")
    var = var[var.gene.isin(genes.gene) & var.mappable.fillna(False).astype(bool)]
    return genes, ind, var


def build_carriers(genes, ind, var, biallelic=False):
    """One row per included individual with features of their variant."""
    keep_ind = ind[ind.gene_selected.fillna(False).astype(bool)]
    if biallelic:
        keep_ind = keep_ind[keep_ind.biallelic == True]  # noqa: E712
    else:
        keep_ind = keep_ind[(keep_ind.biallelic == False) & (keep_ind.n_causative_alleles == 1)]  # noqa: E712
    v = var[var.phenopacket_id.isin(keep_ind.phenopacket_id)].copy()
    if not biallelic:
        v = v[v.n_alleles.fillna(1) == 1]
        v = v.drop_duplicates("phenopacket_id", keep=False)  # exactly one mappable allele
    return v


def add_s6(v, genes):
    s1, s2 = np.full(len(v), np.nan), np.full(len(v), np.nan)
    v = v.reset_index(drop=True)
    for g, sub in v.groupby("gene"):
        f = genes.loc[genes.gene == g, "structure_file"].iloc[0]
        pos = sub.structure_pos.dropna().astype(int).unique().tolist()
        nb = neighbourhood_plddt(f, pos)
        for i, sp in zip(sub.index, sub.structure_pos):
            if pd.notna(sp) and int(sp) in nb:
                s1[i], s2[i] = nb[int(sp)]
    v["nbr1_plddt"] = s1
    v["nbr2_plddt"] = s2
    return v


def feature_table(v, a=""):
    eps = 1e-15
    f = pd.DataFrame({
        "phenopacket_id": v.phenopacket_id.values,
        "gene": v.gene.values,
        "variant_class": v.variant_class.values,
        "rel_pos": v.rel_protein_pos.astype(float).values,
        "S1_logSemi": np.log10(v[f"res_SemiCoreImpact{a}"].astype(float).values + eps),
        "S2_pLDDT": v.res_plddt.astype(float).values,
        "S3_logCore": np.log10(v[f"res_CoreImpact{a}"].astype(float).values + eps),
        "S4_degree": v.res_degree_8A.astype(float).values,
        "S5_AlphaMissense": v.am_pathogenicity.astype(float).values,
        "S6a_nbr1_pLDDT": v.nbr1_plddt.values,
        "S6b_nbr2_pLDDT": v.nbr2_plddt.values,
        "S7_in_domain": v.res_in_domain.astype(float).values,
        "S8_logDisorder": np.log10(v[f"res_DisorderImpact{a}"].astype(float).values + eps),
    })
    return f


CONT = ["S1_logSemi", "S2_pLDDT", "S3_logCore", "S4_degree", "S5_AlphaMissense",
        "S6a_nbr1_pLDDT", "S6b_nbr2_pLDDT", "S8_logDisorder"]


def residualise(f, cols):
    """Within gene, regress each feature on rel_pos (linear OLS); keep residuals."""
    f = f.copy()
    for c in cols:
        r = np.full(len(f), np.nan)
        for g, idx in f.groupby("gene").groups.items():
            sub = f.loc[idx, [c, "rel_pos"]].dropna()
            if len(sub) < 3:
                continue
            X = np.c_[np.ones(len(sub)), sub.rel_pos.values]
            beta, *_ = np.linalg.lstsq(X, sub[c].values, rcond=None)
            r[f.index.get_indexer(sub.index)] = sub[c].values - X @ beta
        f[c + "_res"] = r
    return f


def pheno_matrix(pheno, ids):
    p = pheno[pheno.phenopacket_id.isin(ids)]
    return p.pivot_table(index="phenopacket_id", columns="hpo_id", values="value", aggfunc="max")


# ---------------------------------------------------------------- tests
def run_tests(f, P, feats, perms, rng, label_map):
    """Analysis A: per (gene, term) MWU (continuous) or Fisher (binary) on the
    position-residualised feature, plus a joint within-gene permutation null."""
    rows = []
    perm_counts = {ft: np.zeros(perms, int) for ft in feats}
    for g, fg in f.groupby("gene"):
        fg = fg.reset_index(drop=True)
        ids = fg.phenopacket_id.values
        Pg = P.reindex(ids)
        N = len(fg)
        perm_idx = np.array([rng.permutation(N) for _ in range(perms)]) if perms else None
        for h in Pg.columns:
            y_all = Pg[h].values
            m_ann = ~np.isnan(y_all)
            n_obs = int(np.nansum(y_all[m_ann] == 1))
            n_exc = int(np.nansum(y_all[m_ann] == 0))
            if n_obs < MIN_OBS or n_exc < MIN_OBS:
                continue
            for ft in feats:
                col = ft if ft in ("P0_rel_pos", "S7_in_domain") else ft + "_res"
                x_all = fg[col].values if ft != "P0_rel_pos" else fg.rel_pos.values
                m = m_ann & ~np.isnan(x_all)
                y = y_all[m].astype(int)
                if y.sum() < MIN_OBS or (len(y) - y.sum()) < MIN_OBS:
                    continue
                x = x_all[m]
                if ft == "S7_in_domain":
                    k = int(x.sum())
                    if k == 0 or k == len(x):
                        continue
                    a = int(((x == 1) & (y == 1)).sum())
                    odds, p = fisher_exact([[a, k - a], [y.sum() - a, len(y) - k - y.sum() + a]])
                    eff = np.log(odds) if np.isfinite(odds) and odds > 0 else np.nan
                    if perms:
                        # permuted x values on the annotated subset (joint within gene)
                        xp = x_all[perm_idx][:, m]          # (B, n); nan possible
                        ok = ~np.isnan(xp).any(axis=1)
                        cnt = 0
                        lut_cache = {}
                        pp = np.ones(perms)
                        for b in np.where(ok)[0]:
                            kb = int(xp[b].sum())
                            if kb == 0 or kb == len(y):
                                continue
                            if kb not in lut_cache:
                                lut_cache[kb] = fisher_p_lookup(len(y), int(y.sum()), kb)
                            ab = int(((xp[b] == 1) & (y == 1)).sum())
                            pp[b] = lut_cache[kb][ab]
                        perm_counts[ft] += (pp < 0.05)
                else:
                    p, r = mwu_p_matrix(x[None, :], y)
                    p, eff = float(p[0]), float(r[0])
                    if perms:
                        xp = x_all[perm_idx][:, m]
                        # rows containing nan (feature missing for a permuted-in carrier) are
                        # filled with the observed median, a conservative choice
                        if np.isnan(xp).any():
                            xp = np.where(np.isnan(xp), np.nanmedian(x), xp)
                        pp, _ = mwu_p_matrix(xp, y)
                        perm_counts[ft] += (pp < 0.05)
                rows.append(dict(gene=g, hpo_id=h, hpo_label=label_map.get(h, ""), feature=ft,
                                 n=len(y), n_obs=int(y.sum()), n_exc=int(len(y) - y.sum()),
                                 effect=eff, p=p))
    res = pd.DataFrame(rows)
    summ = []
    for ft in feats:
        r = res[res.feature == ft]
        if len(r) == 0:
            continue
        q = bh(r.p.values)
        res.loc[r.index, "q_feature"] = q
        obs = int((r.p < 0.05).sum())
        pc = perm_counts[ft]
        summ.append(dict(feature=ft, n_tests=len(r), n_genes=r.gene.nunique(),
                         n_p05=obs, expected_p05=0.05 * len(r),
                         perm_mean_p05=float(pc.mean()) if perms else np.nan,
                         perm_p=float((1 + (pc >= obs).sum()) / (1 + perms)) if perms else np.nan,
                         n_q10=int((q < 0.1).sum()), n_genes_q10=int(r.loc[q < 0.1, "gene"].nunique()),
                         pi1=storey_pi1(r.p.values)))
    return res, pd.DataFrame(summ)


def run_class_tests(f, P, perms, rng, label_map):
    """Analysis B: missense vs PTC Fisher per (gene, term), joint within-gene permutation."""
    rows = []
    perm_count = np.zeros(perms, int)
    f = f[f.variant_class.isin(["missense", "nonsense", "frameshift"])].copy()
    f["is_mis"] = (f.variant_class == "missense").astype(int)
    for g, fg in f.groupby("gene"):
        fg = fg.reset_index(drop=True)
        if fg.is_mis.sum() < 5 or (1 - fg.is_mis).sum() < 5:
            continue
        Pg = P.reindex(fg.phenopacket_id.values)
        N = len(fg)
        perm_idx = np.array([rng.permutation(N) for _ in range(perms)]) if perms else None
        xa = fg.is_mis.values
        for h in Pg.columns:
            ya = Pg[h].values
            m = ~np.isnan(ya)
            y = ya[m].astype(int)
            x = xa[m]
            if y.sum() < MIN_OBS or (len(y) - y.sum()) < MIN_OBS:
                continue
            k = int(x.sum())
            if k < 1 or k > len(x) - 1:
                continue
            a = int(((x == 1) & (y == 1)).sum())
            odds, p = fisher_exact([[a, k - a], [y.sum() - a, len(y) - k - y.sum() + a]])
            rows.append(dict(gene=g, hpo_id=h, hpo_label=label_map.get(h, ""), n=len(y),
                             n_obs=int(y.sum()), n_missense=k, log_or=np.log(odds) if 0 < odds < np.inf else np.nan, p=p))
            if perms:
                xp = xa[perm_idx][:, m]
                kb = xp.sum(axis=1)
                ab = ((xp == 1) & (y[None, :] == 1)).sum(axis=1)
                pp = np.ones(perms)
                for kk in np.unique(kb):
                    if kk == 0 or kk == len(y):
                        continue
                    lut = fisher_p_lookup(len(y), int(y.sum()), int(kk))
                    sel = kb == kk
                    pp[sel] = lut[ab[sel]]
                perm_count += (pp < 0.05)
    res = pd.DataFrame(rows)
    res["q"] = bh(res.p.values)
    obs = int((res.p < 0.05).sum())
    summ = dict(feature="B_missense_vs_PTC", n_tests=len(res), n_genes=res.gene.nunique(), n_p05=obs,
                expected_p05=0.05 * len(res), perm_mean_p05=float(perm_count.mean()),
                perm_p=float((1 + (perm_count >= obs).sum()) / (1 + perms)),
                n_q10=int((res.q < 0.1).sum()), n_genes_q10=int(res.loc[res.q < 0.1, "gene"].nunique()),
                pi1=storey_pi1(res.p.values))
    return res, summ


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--data", default="/home/user/data_ext/derived")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    genes, ind, var = load(args.data)
    pheno_prop = pd.read_parquet(Path(args.data) / "phenotypes_propagated.parquet")
    pheno_expl = pd.read_parquet(Path(args.data) / "phenotypes_explicit.parquet")
    pheno_expl = pheno_expl.assign(hpo_id=pheno_expl.hpo_id_current.fillna(pheno_expl.hpo_id))
    labels = dict(zip(pheno_prop.hpo_id, pheno_prop.hpo_label))
    labels.update(dict(zip(pheno_expl.hpo_id, pheno_expl.hpo_label)))

    v = build_carriers(genes, ind, var)
    v = v[v.has_structure_features.fillna(False).astype(bool) | (v.variant_class != "missense")]
    v = add_s6(v.reset_index(drop=True), genes)
    feats_all = ["S1_logSemi", "S2_pLDDT", "S3_logCore", "S4_degree", "S5_AlphaMissense",
                 "S6a_nbr1_pLDDT", "S6b_nbr2_pLDDT", "S7_in_domain", "S8_logDisorder", "P0_rel_pos"]

    def analysis_A(vv, pheno, tag, perms, suffix=""):
        f = feature_table(vv[vv.variant_class == "missense"], suffix)
        f = f[f.S2_pLDDT.notna()]
        f = residualise(f, CONT)
        P = pheno_matrix(pheno, f.phenopacket_id)
        res, summ = run_tests(f, P, feats_all, perms, rng, labels)
        res.to_csv(OUT / f"A_tests_{tag}.csv.gz", index=False)
        summ.insert(0, "analysis", tag)
        return f, res, summ

    summaries = []
    counts = {}
    # primary
    fA, resA, sA = analysis_A(v, pheno_prop, "primary", args.perms)
    summaries.append(sA)
    counts["primary_missense_individuals"] = int(len(fA))
    counts["primary_genes_with_tests"] = int(resA.gene.nunique())
    counts["primary_terms_tested_S1"] = int((resA.feature == "S1_logSemi").sum())
    # reference B
    fB = feature_table(v, "")
    P_B = pheno_matrix(pheno_prop, fB.phenopacket_id)
    resB, sB = run_class_tests(fB, P_B, args.perms, rng, labels)
    resB.to_csv(OUT / "B_tests_primary.csv.gz", index=False)
    sBdf = pd.DataFrame([sB])
    sBdf.insert(0, "analysis", "primary")
    summaries.append(sBdf)
    # sensitivity (fewer permutations)
    sp = max(200, args.perms // 5)
    _, _, s = analysis_A(v, pheno_expl, "sens_explicit", sp); summaries.append(s)
    _, _, s = analysis_A(v, pheno_prop, "sens_alpha015", sp, "_a15"); summaries.append(s)
    _, _, s = analysis_A(v, pheno_prop, "sens_alpha05", sp, "_a5"); summaries.append(s)
    hi = genes.loc[genes.af_mean_plddt >= 50, "gene"]
    _, _, s = analysis_A(v[v.gene.isin(hi)], pheno_prop, "sens_meanpLDDT50", sp); summaries.append(s)
    vb = build_carriers(genes, ind, var, biallelic=True)
    vb = vb[(vb.variant_class == "missense") & vb.has_structure_features.fillna(False).astype(bool)]
    num = ["res_SemiCoreImpact", "res_CoreImpact", "res_DisorderImpact", "res_plddt", "res_degree_8A",
           "am_pathogenicity", "rel_protein_pos", "res_in_domain"]
    vb[num] = vb[num].astype(float)
    agg = vb.groupby("phenopacket_id").agg({**{c: "mean" for c in num}, "gene": "first",
                                            "structure_pos": "first", "variant_class": "first"}).reset_index()
    agg = agg[vb.groupby("phenopacket_id").size().reindex(agg.phenopacket_id).values == 2]
    agg["res_in_domain"] = (agg.res_in_domain >= 0.5).astype(float)
    agg = add_s6(agg.reset_index(drop=True), genes)
    _, _, s = analysis_A(agg, pheno_prop, "sens_biallelic_mean", sp); summaries.append(s)

    S = pd.concat(summaries, ignore_index=True)
    S.to_csv(OUT / "summary.csv", index=False)

    # positive / null controls: best test per gene and feature
    ctrl = ["TBX5", "LMNA", "SPTAN1", "SETD2", "IKZF1", "STXBP1", "POGZ", "EHMT1", "ANKRD11"]
    rows = []
    for g in ctrl:
        for ft in feats_all:
            r = resA[(resA.gene == g) & (resA.feature == ft)]
            if len(r):
                b = r.loc[r.p.idxmin()]
                rows.append(dict(gene=g, analysis="A", feature=ft, n_tests=len(r), n_p05=int((r.p < 0.05).sum()),
                                 best_term=b.hpo_label, best_p=b.p, best_q_feature=b.q_feature, n=b.n))
        r = resB[resB.gene == g]
        if len(r):
            b = r.loc[r.p.idxmin()]
            rows.append(dict(gene=g, analysis="B", feature="missense_vs_PTC", n_tests=len(r),
                             n_p05=int((r.p < 0.05).sum()), best_term=b.hpo_label, best_p=b.p,
                             best_q_feature=b.q, n=b.n))
    pd.DataFrame(rows).to_csv(OUT / "controls.csv", index=False)

    # per-gene summary for the primary feature set
    pg = (resA.groupby(["gene", "feature"]).agg(n_tests=("p", "size"), n_p05=("p", lambda s: int((s < 0.05).sum())),
                                                 min_p=("p", "min")).reset_index())
    pg.to_csv(OUT / "per_gene_A.csv", index=False)
    counts["B_individuals"] = int(len(fB))
    json.dump(counts, open(OUT / "counts.json", "w"), indent=1)
    pd.set_option("display.width", 200)
    print(S.to_string(index=False))
    print(json.dumps(counts, indent=1))


if __name__ == "__main__":
    main()
