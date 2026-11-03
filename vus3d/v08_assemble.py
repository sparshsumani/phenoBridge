#!/usr/bin/env python3
"""Step 8: final tables in derived_vus3d/: variants.parquet, genes.parquet, counts.json
(README.md is maintained by hand next to them).

variants.parquet = one row per usable train/test VariationID (v05 'usable': mapped to MANE, single-fragment
AlphaFold v6 model, AF position mapped, reference aa matches the model).  Residue features of feature
proteins are attached by (uniprot, af_pos) with prefix res_.  AlphaMissense exact substitution from v06.
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from common import INTER, OUT, RFDIR, AA1TO3, CUTOFF_DATE, save_json, load_json

RES_DROP = {"gene", "uniprot", "pos", "ca_x", "ca_y", "ca_z"}


def main():
    u = pd.read_parquet(INTER / "variants_usable.parquet")
    g = pd.read_parquet(INTER / "genes.parquet")
    ex = pd.read_parquet(INTER / "am_exact.parquet")
    meta = pd.read_parquet(INTER / "structure_meta.parquet")
    feat = set(meta.uniprot)
    out = pd.DataFrame(dict(
        variation_id=u.variation_id, allele_id=u.allele_id, gene=u.gene, hgnc_id=u.hgnc_id, uniprot=u.uniprot,
        mane_np=u.mane_np, mane_enst=u.mane_enst, protein_pos=u.mane_pos.astype(int), af_pos=u.af_pos.astype(int),
        uniprot_pos=u.uniprot_pos.astype("Int64"), ref=u.mane_ref, alt=u.mane_alt))
    out["hgvs_p"] = ["p." + AA1TO3[r] + str(int(p)) + AA1TO3[a] for r, p, a in zip(out.ref, out.protein_pos, out.alt)]
    out["clinvar_name"] = u.name; out["clinvar_hgvs_p"] = u.hgvs_p; out["clinvar_hgvs_c"] = u.hgvs_c
    out["clinvar_name_transcript"] = u.name_tx; out["map_method"] = u.map_method
    out["chrom"], out["pos_vcf"], out["ref_vcf"], out["alt_vcf"] = u.chrom, u.pos_vcf.astype(int), u.ref_vcf, u.alt_vcf
    out["rsid"] = u.rsid.replace({"-1": None})
    out["set"], out["label"] = u.set, u.label
    out["prior_status_2023"] = u.prior_status_2023
    for c in ["status_2023", "category_2023", "stars_2023", "review_status_2023", "last_evaluated_2023", "n_submitters_2023", "status_2023_lookup",
              "status_2026", "category_2026", "stars_2026", "review_status_2026", "last_evaluated_2026", "n_submitters_2026", "status_2026_lookup"]:
        out[c] = u[c].to_numpy()
    out["stars_2023"] = pd.array(out.stars_2023, dtype="Int64"); out["stars_2026"] = pd.array(out.stars_2026, dtype="Int64")
    out["n_submitters_2023"] = pd.array(out.n_submitters_2023, dtype="Int64"); out["n_submitters_2026"] = pd.array(out.n_submitters_2026, dtype="Int64")
    out["n_submitters"] = np.where(out.set == "train", out.n_submitters_2023, out.n_submitters_2026)
    out["n_submitters"] = pd.array(out.n_submitters, dtype="Int64")
    out["last_evaluated_2023"] = pd.to_datetime(out.last_evaluated_2023); out["last_evaluated_2026"] = pd.to_datetime(out.last_evaluated_2026)
    out["last_evaluated_2026_after_cutoff"] = pd.array(np.where(out.last_evaluated_2026.isna(), None, out.last_evaluated_2026 > pd.Timestamp(CUTOFF_DATE)), dtype="boolean")
    out["residue_key"], out["protein_change_key"] = u.residue_key, u.protein_change_key
    out["test_residue_has_train_label"] = pd.array(u.test_residue_has_train_label, dtype="boolean")
    out["test_same_protein_change_in_train"] = pd.array(u.test_same_protein_change_in_train, dtype="boolean")
    out["n_train_variants_same_residue"] = u.n_train_variants_same_residue
    out["n_vids_same_protein_change_same_set"] = u.n_vids_same_protein_change_same_set
    gi = g.set_index("gene")
    out["gene_n_train_residues"] = out.gene.map(gi.n_train_residues).astype(int)
    for k in (5, 10, 20):
        out[f"gene_evaluable_{k}"] = out.gene.map(gi[f"evaluable_{k}"]).astype(bool)
    out["protein_length_af"] = u.protein_length_af.astype(int)
    out = out.merge(ex, on="variation_id", how="left")
    out["has_residue_features"] = out.uniprot.isin(feat)
    # residue features
    parts = []
    for acc in sorted(feat & set(out.uniprot)):
        rf = pd.read_parquet(RFDIR / f"{acc}.parquet")
        rf = rf[[c for c in rf.columns if c not in RES_DROP or c in ("uniprot", "pos")]]
        parts.append(rf)
    rf = pd.concat(parts, ignore_index=True)
    rf = rf.rename(columns={c: "res_" + c for c in rf.columns if c not in ("uniprot", "pos")}).rename(columns={"pos": "af_pos"})
    n0 = len(out)
    out = out.merge(rf, on=["uniprot", "af_pos"], how="left")
    assert len(out) == n0
    out["res_aa_matches_ref"] = pd.array(np.where(out.has_residue_features, out.res_aa_af == out.ref, None), dtype="boolean")
    out = out.sort_values(["set", "gene", "af_pos", "alt", "variation_id"]).reset_index(drop=True)
    out.to_parquet(OUT / "variants.parquet", index=False)

    g["has_residue_features"] = g.uniprot.isin(feat)
    mcols = ["L", "mean_plddt", "frac_plddt_lt50", "frac_plddt_ge70", "q_core", "q_semi", "G_all_edges", "n_edges", "am_accession", "am_mapping", "source"]
    mi = meta.set_index("uniprot")
    for c in mcols:
        if c in mi:
            g["struct_" + c] = g.uniprot.map(mi[c])
    g.to_parquet(OUT / "genes.parquet", index=False)

    counts = {}
    for k, f in [("downloads", "downloads.json"), ("clinvar_parse", "clinvar_parse_log.json"), ("sets", "sets_log.json"),
                 ("protein_mapping", "protein_mapping_log.json"), ("genes", "genes_log.json"), ("alphamissense", "am_log.json"),
                 ("features", "features_log.json")]:
        p = INTER / f
        if p.exists():
            counts[k] = load_json(p)
    for k in ("clinvar_2023", "clinvar_2026"):
        counts["downloads"].get(k, {}).pop("sha256_note", None)
    # strip bulky header list
    tr, te = out[out.set == "train"], out[out.set == "test"]
    counts["final"] = dict(
        variants=len(out), by_set_label={f"{a}|{b}": int(n) for (a, b), n in out.groupby(["set", "label"]).size().items()},
        test_by_label_prior={f"{a}|{b}": int(n) for (a, b), n in te.groupby(["label", "prior_status_2023"]).size().items()},
        train_by_label_2026_category={f"{a}|{b}": int(n) for (a, b), n in tr.groupby(["label", "category_2026"]).size().items()},
        test_stars_2026={f"{a}|{b}": int(n) for (a, b), n in te.groupby(["label", "stars_2026"]).size().items()},
        genes=int(out.gene.nunique()), proteins_with_residue_features=int(len(feat)),
        variants_with_residue_features={f"{a}|{b}": int(n) for (a, b), n in out[out.has_residue_features].groupby(["set", "label"]).size().items()},
        am_status={f"{a}|{b}": int(n) for (a, b), n in out.groupby(["set", "am_status"]).size().items()},
        am_coverage={s_: float((d.am_status == "ok").mean()) for s_, d in out.groupby("set")},
        test_last_evaluated_2026_after_cutoff=dict(n=len(te), with_date=int(te.last_evaluated_2026.notna().sum()),
                                                  after=int((te.last_evaluated_2026_after_cutoff == True).sum()),
                                                  by_prior={p: dict(n=len(d), after=int((d.last_evaluated_2026_after_cutoff == True).sum()),
                                                                    no_date=int(d.last_evaluated_2026.isna().sum())) for p, d in te.groupby("prior_status_2023")}),
        res_aa_matches_ref=f"{int((out.res_aa_matches_ref == True).sum())}/{int(out.res_aa_matches_ref.notna().sum())}")
    save_json(counts, OUT / "counts.json")
    print(json.dumps(counts["final"], indent=1, default=str))


if __name__ == "__main__":
    main()
