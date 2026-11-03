"""Why does the individual-level test mis-calibrate? (exploratory; run after extra_analyses.py)

A structural feature is a property of the variant, so it is constant among carriers of
the same variant. If carriers of the same variant also resemble each other in a
phenotype (true variant effect, same family, same publication), the effective sample
size is smaller than the number of individuals. Kish's design effect for a correlation
between a cluster-constant x and a clustered y is DE = 1 + (m~ - 1) * ICC_y, with
m~ = sum(m_k^2) / N. A rank test that ignores clustering then has a null rejection
rate of about 2 * (1 - Phi(1.96 / sqrt(DE))) at alpha = 0.05.

We compute DE for every primary S1 test and compare the predicted rate with the
variant-level permutation null rate measured in extra_analyses.py.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_analysis as ra  # noqa: E402

OUT = HERE / "out"
DATA = "/home/user/data_ext/derived"


def icc_anova(y, groups):
    """One-way ANOVA estimator of the intraclass correlation (truncated at 0)."""
    df = pd.DataFrame({"y": y, "g": groups})
    k = df.g.nunique()
    N = len(df)
    if k < 2 or k == N:
        return 0.0
    m = df.groupby("g").size().values
    means = df.groupby("g").y.mean()
    grand = df.y.mean()
    msb = (m * (means.values - grand) ** 2).sum() / (k - 1)
    msw = ((df.y - df.g.map(means)) ** 2).sum() / (N - k)
    m0 = (N - (m ** 2).sum() / N) / (k - 1)
    den = msb + (m0 - 1) * msw
    return float(max((msb - msw) / den, 0.0)) if den > 0 else 0.0


def main():
    genes, ind, var = ra.load(DATA)
    pheno = pd.read_parquet(Path(DATA) / "phenotypes_propagated.parquet")
    v = ra.build_carriers(genes, ind, var)
    v = v[v.has_structure_features.fillna(False).astype(bool) & (v.variant_class == "missense")].copy()
    v["nbr1_plddt"] = np.nan
    v["nbr2_plddt"] = np.nan
    f = ra.feature_table(v, "")
    f = f[f.S2_pLDDT.notna() & f.S1_logSemi.notna()]
    P = ra.pheno_matrix(pheno, f.phenopacket_id)
    resA = pd.read_csv(OUT / "A_tests_primary.csv.gz")
    t = resA[resA.feature == "S1_logSemi"]

    rows = []
    for (g, h), _ in t.groupby(["gene", "hpo_id"]):
        fg = f[f.gene == g]
        y = P.reindex(fg.phenopacket_id)[h].values
        m = ~np.isnan(y)
        yy, kk = y[m], fg.vkey.values[m]
        icc = icc_anova(yy, kk)
        sizes = pd.Series(kk).value_counts().values
        mt = (sizes ** 2).sum() / sizes.sum()
        de = 1 + (mt - 1) * icc
        rows.append(dict(gene=g, hpo_id=h, n=int(m.sum()), n_variants=len(sizes), m_tilde=mt, icc=icc, de=de,
                         pred_rate=2 * (1 - norm.cdf(1.96 / np.sqrt(de)))))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "X_design_effect_tests.csv", index=False)
    G = T.groupby("gene").agg(n_tests=("de", "size"), mean_de=("de", "mean"), median_icc=("icc", "median"),
                              pred_rate=("pred_rate", "mean")).reset_index()
    cal = pd.read_csv(OUT / "X_per_gene_calibration.csv")
    G = G.merge(cal[["gene", "carriers_per_variant", "null_rate_variant_S1_logSemi",
                     "null_rate_individual_S1_logSemi"]], on="gene")
    G.to_csv(OUT / "X_design_effect_genes.csv", index=False)
    ok = G[G.n_tests >= 1]
    rho, p = spearmanr(ok.mean_de, ok.null_rate_variant_S1_logSemi)
    rho2, p2 = spearmanr(ok.pred_rate, ok.null_rate_variant_S1_logSemi)

    # publications behind recurrent variants (genes with tests)
    pm = ind.set_index("phenopacket_id").pmid
    ft = f[f.gene.isin(t.gene.unique())].assign(pmid=lambda d: d.phenopacket_id.map(pm))
    by = ft.groupby("vkey").agg(carriers=("pmid", "size"), pmids=("pmid", "nunique"))
    rec = by[by.carriers > 1]
    # same-variant pairs from the same publication
    same_pub = tot = 0
    for _, s in ft.groupby("vkey"):
        n = len(s)
        if n < 2:
            continue
        c = s.pmid.value_counts().values
        same_pub += (c * (c - 1) // 2).sum()
        tot += n * (n - 1) // 2
    top = by.sort_values("carriers", ascending=False).head(8)
    out = dict(
        tests=int(len(T)), median_icc=float(T.icc.median()), frac_icc_gt_0p1=float((T.icc > .1).mean()),
        median_de=float(T.de.median()), mean_de=float(T.de.mean()), max_de=float(T.de.max()),
        pooled_pred_p05=float(T.pred_rate.sum()),
        pooled_perm_variant_mean_p05=float(np.load(OUT / "A_permnull_primary.npz")["S1_logSemi"].mean()),
        spearman_gene_de_vs_variant_null_rate=[float(rho), float(p)],
        spearman_gene_pred_vs_variant_null_rate=[float(rho2), float(p2)],
        recurrent_variants=int(len(rec)), recurrent_variants_one_pmid=int((rec.pmids == 1).sum()),
        individuals_in_recurrent=int(rec.carriers.sum()),
        same_variant_pairs=int(tot), same_variant_pairs_same_pmid=int(same_pub),
        frac_same_variant_pairs_same_pmid=float(same_pub / tot),
        n_pmids_total=int(ft.pmid.nunique()),
        top_variants=[dict(vkey=k, carriers=int(r.carriers), pmids=int(r.pmids)) for k, r in top.iterrows()],
    )
    json.dump(out, open(OUT / "X_design_effect.json", "w"), indent=1)
    print(json.dumps(out, indent=1))
    print(G.sort_values("mean_de").to_string())


if __name__ == "__main__":
    main()
