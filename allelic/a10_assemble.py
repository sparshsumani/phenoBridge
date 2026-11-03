#!/usr/bin/env python3
"""Step 10: final tables in derived_allelic/: variants.parquet, genes.parquet, diseases.parquet,
gene_disease_counts.parquet, counts.json.  (README.md is maintained by hand next to them.)
Residue features are attached to variants of feature proteins by (uniprot, af_pos) with prefix res_.
AlphaMissense exact substitution: feature proteins -> residue am_pos (AF->AM mapping) + alt;
other proteins -> exact '<uniprot> <ref><uniprot_pos><alt>' row (ref must match AM by construction).
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from common import INTER, OUT, RFDIR, GRAPHDIR, AA1TO3, save_json, load_json


def main():
    v = pd.read_parquet(INTER / "variants_positioned.parquet")
    g = pd.read_parquet(INTER / "genes.parquet")
    d = pd.read_parquet(INTER / "diseases.parquet")
    gd = pd.read_parquet(INTER / "gene_disease_counts.parquet")
    am = pd.read_parquet(INTER / "alphamissense_subset.parquet")
    meta = pd.read_parquet(INTER / "structure_meta.parquet")
    feat_prot = set(meta.uniprot)
    dn = dict(zip(d.disease_unit, d.disease_name)); dm = dict(zip(d.disease_unit, d.mondo))
    di = dict(zip(d.disease_unit, d.inheritance_clean)); dc = dict(zip(d.disease_unit, d.inheritance_classes))
    db = dict(zip(d.disease_unit, d.broad_unit))

    out = pd.DataFrame(dict(
        gene=v.gene, hgnc_id=v.hgnc_id, uniprot=v.uniprot, mane_enst=v.hgnc_id.map(g.set_index("hgnc_id").mane_enst),
        mane_np=v.hgnc_id.map(g.set_index("hgnc_id").mane_np),
        protein_pos=v.mane_pos, af_pos=v.af_pos, uniprot_pos=v.uniprot_pos, ref=v.mane_ref, alt=v.mane_alt))
    out["hgvs_p"] = ["p." + AA1TO3[r] + str(int(p)) + AA1TO3[a] for r, p, a in zip(out.ref, out.protein_pos, out.alt)]
    out["variant_key"] = out.gene + ":" + out.ref + out.protein_pos.astype(str) + out.alt
    out["source"] = v.source
    out["disease_units"] = v.clinvar_units
    out["n_disease_units"] = v.clinvar_n_units
    out["single_label"] = v.clinvar_single_label
    out["multi_label"] = v.clinvar_n_units > 1
    lab = lambda us: " | ".join(str(dn.get(u)) for u in us) or None
    out["disease_label"] = v.clinvar_units.map(lab)
    out["disease_mondo"] = v.clinvar_units.map(lambda us: ";".join(str(dm.get(u)) for u in us) or None)
    out["disease_unit"] = [us[0] if len(us) == 1 else None for us in v.clinvar_units]
    out["inheritance"] = [di.get(us[0]) if len(us) == 1 else None for us in v.clinvar_units]
    out["inheritance_per_unit"] = v.clinvar_units.map(lambda us: ";".join(f"{u}={di.get(u)}" for u in us) or None)
    out["inheritance_classes"] = v.clinvar_units.map(lambda us: sorted({c for u in us for c in (dc.get(u) if dc.get(u) is not None else [])}))
    out["unit_is_broad"] = [bool(db.get(us[0])) if len(us) == 1 else None for us in v.clinvar_units]
    out["omim_units"] = v.omim_units
    out["omim_single_label"] = v.omim_single_label
    out["combined_units"] = v.combined_units
    out["combined_n_units"] = v.combined_n_units
    out["combined_single_label"] = v.combined_single_label
    out["combined_disease_unit"] = [us[0] if len(us) == 1 else None for us in v.combined_units]
    out["combined_inheritance"] = [di.get(us[0]) if len(us) == 1 else None for us in v.combined_units]
    out["labels_agree"] = v.labels_agree
    for c in ["clinvar_variation_ids", "n_clinvar_records", "clinvar_hgvs_p", "clinvar_hgvs_c", "clinvar_mapping_method",
              "clinvar_classification", "stars", "review_status", "n_submitters", "n_submitters_sum", "last_evaluated", "rsids",
              "clinvar_conditions_kept", "clinvar_conditions_all", "clinvar_has_generic_condition",
              "clinvar_has_nonspecific_condition", "clinvar_has_grouping_condition", "clinvar_dropped_cancer_condition",
              "clinvar_cancer_flag", "ppkt_hgvs_p", "ppkt_n_carriers", "ppkt_n_pmids", "ppkt_n_carriers_single_allele",
              "ppkt_n_pmids_single_allele", "ppkt_n_carriers_hemizygous", "ppkt_n_carriers_biallelic", "ppkt_n_pmids_biallelic",
              "ppkt_units_all", "ppkt_units_single_allele", "ppkt_units_biallelic", "ppkt_disease_labels", "ppkt_cohorts", "ppkt_pmids"]:
        out[c] = v[c]
    out["stars"] = out.stars.astype("Int64"); out["n_submitters"] = out.n_submitters.astype("Int64")
    out["n_submitters_sum"] = out.n_submitters_sum.astype("Int64"); out["n_clinvar_records"] = out.n_clinvar_records.astype("Int64")
    out["af_n_fragments"] = v.af_n_fragments.astype("Int64")
    out["af_single_fragment"] = v.af_single_fragment
    out["ref_matches_af"] = v.ref_matches_af
    out["ref_matches_uniprot"] = v.ref_matches_uniprot
    out["protein_length"] = v.protein_length_af
    out["rel_protein_pos"] = out.af_pos / out.protein_length
    out["countable"] = v.countable
    out["gene_eligible"] = out.hgnc_id.map(g.set_index("hgnc_id").eligible).fillna(False)
    out["gene_eligible_any_definition"] = out.hgnc_id.map(g.set_index("hgnc_id").eligible_any_definition).fillna(False)
    out["gene_ad_ar_eligible"] = out.hgnc_id.map(g.set_index("hgnc_id").ad_ar_eligible).fillna(False)
    # unit-level >=5 flag for the variant's (single) unit within its gene
    prim = gd[gd.definition == "clinvar"].set_index(["hgnc_id", "disease_unit"]).n_variants_countable
    out["unit_n_variants_in_gene"] = [prim.get((h, u)) if isinstance(u, str) else None for h, u in zip(out.hgnc_id, out.disease_unit)]
    out["unit_n_variants_in_gene"] = out.unit_n_variants_in_gene.astype("Int64")
    out["unit_ge5_in_gene"] = out.unit_n_variants_in_gene >= 5
    out["has_residue_features"] = out.uniprot.isin(feat_prot) & out.countable

    # residue features
    rf = pd.concat([pd.read_parquet(RFDIR / f"{a}.parquet") for a in sorted(feat_prot)], ignore_index=True)
    rf = rf.drop(columns=["gene"]).rename(columns={"pos": "af_pos"})
    rf["af_pos"] = rf.af_pos.astype("Int64")
    rfr = rf.rename(columns={c: "res_" + c for c in rf.columns if c not in ("uniprot", "af_pos")})
    out = out.merge(rfr, on=["uniprot", "af_pos"], how="left")
    # AlphaMissense exact substitution
    amk = am.set_index(["uniprot_id", "pos", "alt"])[["am_pathogenicity", "am_class", "ref"]]
    amk = amk[~amk.index.duplicated()]
    path, cls, st = [], [], []
    amacc = dict(zip(meta.uniprot, meta.am_accession))
    for a, ap, up, ref, alt, haf in zip(out.uniprot, out.res_am_pos, out.uniprot_pos, out.ref, out.alt, out.has_residue_features):
        if haf:
            acc, p = amacc.get(a), ap
        else:
            acc, p = a, up
        if acc is None or pd.isna(p) or (acc, int(p), alt) not in amk.index:
            path.append(np.nan); cls.append(None); st.append("not_found"); continue
        r = amk.loc[(acc, int(p), alt)]
        if r.ref != ref:
            path.append(np.nan); cls.append(None); st.append("ref_mismatch"); continue
        path.append(float(r.am_pathogenicity)); cls.append(r.am_class); st.append("ok")
    out["am_pathogenicity"] = path; out["am_class"] = cls; out["am_status"] = st
    out = out.sort_values(["gene", "protein_pos", "alt"]).reset_index(drop=True)
    out.to_parquet(OUT / "variants.parquet", index=False)

    # genes
    gg = g.drop(columns=[c for c in ("af_pdb",) if c in g.columns]).copy()
    gg["af_pdb"] = g.af_pdb
    gg["has_residue_features"] = gg.uniprot.isin(feat_prot)
    gg["residue_features_file"] = [str(RFDIR / f"{a}.parquet") if a in feat_prot else None for a in gg.uniprot]
    gg["contact_graph_file"] = [str(GRAPHDIR / f"{a}.parquet") if a in feat_prot else None for a in gg.uniprot]
    sm = meta.set_index("uniprot")
    for c in ("L", "mean_plddt", "frac_plddt_lt50", "frac_plddt_ge70", "q_core", "q_semi", "G_all_edges", "backbone_size", "am_accession", "am_mapping"):
        gg["struct_" + c] = gg.uniprot.map(sm[c]) if c in sm else None
    gg.to_parquet(OUT / "genes.parquet", index=False)

    # diseases
    dd = d.copy()
    pr = gd[gd.definition == "clinvar"]
    dd["n_single_label_variants_countable"] = dd.disease_unit.map(pr.groupby("disease_unit").n_variants_countable.sum()).fillna(0).astype(int)
    dd["n_genes"] = dd.disease_unit.map(pr[pr.n_variants_countable > 0].groupby("disease_unit").size()).fillna(0).astype(int)
    dd["genes"] = dd.disease_unit.map(pr[pr.n_variants_countable > 0].groupby("disease_unit").gene.agg(lambda s: ";".join(sorted(set(s)))))
    dd["inheritance_terms"] = dd.inheritance_terms.map(list)
    dd["inheritance_classes"] = dd.inheritance_classes.map(list)
    dd.to_parquet(OUT / "diseases.parquet", index=False)
    gd.to_parquet(OUT / "gene_disease_counts.parquet", index=False)

    logs = {}
    for k, f in [("downloads", "downloads.json"), ("clinvar_filter", "clinvar_filter_log.json"),
                 ("disease_labels", "disease_label_log.json"), ("protein_mapping", "protein_mapping_log.json"),
                 ("ppkt_merge", "ppkt_merge_log.json"), ("inheritance", "inheritance_log.json"), ("genes", "genes_log.json"),
                 ("alphamissense", "am_log.json"), ("features", "features_log.json"), ("sanity", "sanity.json")]:
        if (INTER / f).exists():
            logs[k] = load_json(INTER / f)
    if "clinvar_filter" in logs:
        logs["clinvar_filter"].pop("header_columns", None)
        logs["clinvar_filter"].pop("germline_classification_counts_GRCh38_SNV", None)
    if "downloads" in logs:
        for r in logs["downloads"].values():
            r.pop("sha256", None)
    logs["final"] = dict(
        variants=len(out), variants_by_source={k: int(x) for k, x in out.source.value_counts().items()},
        variants_single_label=int(out.single_label.sum()), variants_multi_label=int(out.multi_label.sum()),
        variants_no_unit=int((out.n_disease_units == 0).sum()),
        variants_countable=int(out.countable.sum()), variants_single_label_countable=int((out.single_label & out.countable).sum()),
        variants_with_residue_features=int(out.has_residue_features.sum()),
        variants_with_residue_features_single_label=int((out.has_residue_features & out.single_label).sum()),
        variants_in_eligible_genes_single_label_unit_ge5=int((out.gene_eligible & out.single_label & out.countable & (out.unit_ge5_in_gene == True)).sum()),
        am_status={k: int(x) for k, x in out.am_status.value_counts().items()},
        am_status_feature_proteins={k: int(x) for k, x in out[out.has_residue_features].am_status.value_counts().items()},
        genes=len(gg), feature_proteins=len(feat_prot), diseases=len(dd),
        residue_rows=len(rf))
    save_json(logs, OUT / "counts.json")
    print(json.dumps(logs["final"], indent=1))


if __name__ == "__main__":
    main()
