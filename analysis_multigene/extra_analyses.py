"""Follow-up analyses for the paper (run after run_analysis.py).

1. Variant recurrence in the primary missense set.
2. Per-gene calibration: how many tests reach p < 0.05 under individual-level vs
   variant-level permutation, against how often carriers share a variant.
3. Exploratory check: S1 at alpha = 0.5 residualised on a spline of position.
4. Control genes and the class test, scored with the per-test variant-level permutation p.

Everything labelled exploratory here was not in PLAN.md.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_analysis as ra  # noqa: E402

OUT = HERE / "out"
DATA = "/home/user/data_ext/derived"
PERMS_GENE = 500
PERMS_EXPL = 1000


def main():
    rng = np.random.default_rng(ra.RNG_SEED + 1)
    genes, ind, var = ra.load(DATA)
    pheno = pd.read_parquet(Path(DATA) / "phenotypes_propagated.parquet")
    labels = dict(zip(pheno.hpo_id, pheno.hpo_label))

    v = ra.build_carriers(genes, ind, var)
    v = v[v.has_structure_features.fillna(False).astype(bool) | (v.variant_class != "missense")]
    v = ra.add_s6(v.reset_index(drop=True), genes)
    vm = v[v.variant_class == "missense"]

    f = ra.feature_table(vm, "")
    f = f[f.S2_pLDDT.notna()]
    f = ra.residualise(f, ra.CONT)
    P = ra.pheno_matrix(pheno, f.phenopacket_id)
    resA = pd.read_csv(OUT / "A_tests_primary.csv.gz")
    tested = sorted(resA.gene.unique())

    out = {}
    # ---------------------------------------------------------------- 1. recurrence
    ft = f[f.gene.isin(tested)]
    cpv = ft.groupby("vkey").size()
    out["recurrence"] = dict(
        individuals=int(len(ft)), distinct_variants=int(len(cpv)),
        variants_with_1_carrier=int((cpv == 1).sum()),
        individuals_sharing_variant=int(cpv[cpv > 1].sum()),
        frac_individuals_sharing=float(cpv[cpv > 1].sum() / len(ft)),
        max_carriers=int(cpv.max()), max_variant=str(cpv.idxmax()),
        top10=[(k, int(n)) for k, n in cpv.sort_values(ascending=False).head(10).items()],
    )
    cpv.rename("carriers").to_csv(OUT / "X_carriers_per_variant.csv")

    # ---------------------------------------------------------------- 2. per-gene calibration
    rows = []
    feats = ["S1_logSemi", "P0_rel_pos"]
    for g in tested:
        fg = f[f.gene == g]
        rec = {"gene": g, "n": len(fg), "n_variants": fg.vkey.nunique()}
        for level in ("individual", "variant"):
            _, s, pc = ra.run_tests(fg, P, feats, PERMS_GENE, rng, labels, level=level)
            for ftn in feats:
                srow = s[s.feature == ftn]
                if len(srow) == 0:
                    continue
                nt = int(srow.n_tests.iloc[0])
                rec[f"n_tests_{ftn}"] = nt
                rec[f"null_rate_{level}_{ftn}"] = float(pc[ftn].mean() / nt)
        rows.append(rec)
    pg = pd.DataFrame(rows)
    pg["carriers_per_variant"] = pg.n / pg.n_variants
    pg.to_csv(OUT / "X_per_gene_calibration.csv", index=False)
    ok = pg.dropna(subset=["null_rate_variant_S1_logSemi", "null_rate_individual_S1_logSemi"])
    from scipy.stats import spearmanr
    rho, p_rho = spearmanr(ok.carriers_per_variant, ok.null_rate_variant_S1_logSemi)
    out["per_gene_calibration"] = dict(
        genes=int(len(ok)),
        median_null_rate_individual_S1=float(ok.null_rate_individual_S1_logSemi.median()),
        median_null_rate_variant_S1=float(ok.null_rate_variant_S1_logSemi.median()),
        genes_variant_null_rate_gt_0p10=int((ok.null_rate_variant_S1_logSemi > 0.10).sum()),
        spearman_rho_cpv_vs_variant_null_rate=float(rho), spearman_p=float(p_rho),
    )

    # ---------------------------------------------------------------- 3. exploratory alpha 0.5 + spline
    f5 = ra.feature_table(vm, "_a5")
    f5 = f5[f5.S2_pLDDT.notna()]
    f5 = ra.residualise(f5, ra.CONT, spline=True)
    P5 = ra.pheno_matrix(pheno, f5.phenopacket_id)
    _, s5, _ = ra.run_tests(f5, P5, ["S1_logSemi", "S5_AlphaMissense", "P0_rel_pos"], PERMS_EXPL, rng,
                            labels, level="variant")
    s5.insert(0, "analysis", "explor_alpha05_spline")
    s5.to_csv(OUT / "X_alpha05_spline_summary.csv", index=False)
    out["alpha05_spline"] = s5[["feature", "n_tests", "n_p05", "perm_mean_p05", "perm_p"]].to_dict("records")

    # ---------------------------------------------------------------- 4. controls with permutation p
    ctrl = ["TBX5", "LMNA", "SPTAN1", "SETD2", "IKZF1", "STXBP1", "POGZ", "EHMT1", "ANKRD11"]
    ind_res = pd.read_csv(OUT / "A_tests_primary_individual_perm.csv.gz")
    crow = []
    for g in ctrl:
        r = resA[resA.gene == g]
        ri = ind_res[ind_res.gene == g]
        for ftn in ["S1_logSemi", "S2_pLDDT", "S5_AlphaMissense", "P0_rel_pos"]:
            rr = r[r.feature == ftn]
            rri = ri[ri.feature == ftn]
            if len(rr) == 0:
                continue
            crow.append(dict(gene=g, feature=ftn, n_tests=len(rr),
                             n_carriers=int(rr.n.max()),
                             n_p05_analytic=int((rr.p < 0.05).sum()),
                             n_p05_perm_individual=int((rri.p_perm < 0.05).sum()),
                             n_p05_perm_variant=int((rr.p_perm < 0.05).sum()),
                             min_p_analytic=float(rr.p.min()),
                             min_p_perm_variant=float(rr.p_perm.min())))
    C = pd.DataFrame(crow)
    C.to_csv(OUT / "X_controls_perm.csv", index=False)
    nv = f[f.gene.isin(ctrl)].groupby("gene").agg(n=("vkey", "size"), n_variants=("vkey", "nunique"))
    out["controls_missense_n_variants"] = nv.reset_index().to_dict("records")

    B = pd.read_csv(OUT / "B_tests_primary.csv.gz")
    Bi = pd.read_csv(OUT / "B_tests_individual_perm.csv.gz")
    qv = ra.bh(B.p_perm.values)
    qi = ra.bh(Bi.p_perm.values)
    out["class_test"] = dict(
        tests=int(len(B)), genes=int(B.gene.nunique()),
        bh_analytic_q10=int((B.q < 0.1).sum()), bh_analytic_q10_genes=int(B.loc[B.q < 0.1, "gene"].nunique()),
        bh_perm_variant_q10=int((qv < 0.1).sum()), bh_perm_variant_q10_genes=int(B.loc[qv < 0.1, "gene"].nunique()),
        bh_perm_individual_q10=int((qi < 0.1).sum()),
        bh_perm_individual_q10_genes=int(Bi.loc[qi < 0.1, "gene"].nunique()),
        lost_genes=sorted(set(B.loc[B.q < 0.1, "gene"]) - set(B.loc[qv < 0.1, "gene"])),
        kept_genes=sorted(set(B.loc[qv < 0.1, "gene"])),
    )

    json.dump(out, open(OUT / "X_extra.json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
