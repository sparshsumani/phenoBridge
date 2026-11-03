#!/usr/bin/env python3
"""Step 5: phenopacket-store missense variants and merge with ClinVar at the protein level.

Variant key (both sources): (HGNC id, MANE Select protein position, ref aa, alt aa).

ClinVar: mapped variants (a04 map_status == ok) x disease units (a03).  Several ClinVar VariationIDs
  (different codon changes) can give the same protein change; they are collapsed: union of disease
  units (then ancestor-redundant MONDO units removed again), max stars, max/sum of submitters,
  latest LastEvaluated, all VariationIDs/conditions concatenated.
phenopacket-store (v0.1.27, multigene derived tables): rows with variant_class == missense (VEP on
  MANE Select; protein_pos = MANE position, vep_p_ref/alt).  Each carrier is classed by the
  individual's genotype: single_allele = heterozygous with no second causative allele (zygosity het)
  or hemizygous; biallelic = homozygous or compound heterozygous (individual.biallelic).  Disease =
  individual's disease_id (all OMIM in v0.1.27).  Counts: carriers (phenopackets) and distinct PMIDs.
Merge: outer join on the key.  source = clinvar | ppkt | both.
  ppkt_units_all = set of OMIM ids over all carriers; labels_agree = identical | overlap | disjoint
  (ClinVar units vs ppkt_units_all; NA if either side empty).
  combined_units = ClinVar units U ppkt_units_all; combined_single_label = |combined_units| == 1.
Outputs intermediate/variants_merged.parquet, intermediate/ppkt_merge_log.json
"""
import sys, json
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from common import INTER, DERIVED_PPKT, save_json

KEY = ["hgnc_id", "mane_pos", "mane_ref", "mane_alt"]


def setcmp(a, b):
    a, b = set(a or []), set(b or [])
    if not a or not b:
        return None
    if a == b:
        return "identical"
    return "overlap" if a & b else "disjoint"


def main():
    log = Counter()
    mp = pd.read_parquet(INTER / "clinvar_mapped.parquet")
    cv = pd.read_parquet(INTER / "clinvar_missense.parquet")
    vu = pd.read_parquet(INTER / "clinvar_variant_units.parquet")
    mi = pd.read_parquet(INTER / "mondo_index.parquet")
    du = pd.read_parquet(INTER / "disease_units_clinvar.parquet")
    gm = pd.read_parquet(INTER / "gene_protein_map.parquet")
    ok = mp[mp.map_status == "ok"].merge(
        cv[["variation_id", "clinvar_classification", "review_status", "stars", "n_submitters", "last_evaluated",
            "rsid", "rcv", "origin_simple", "phenotype_list"]], on="variation_id").merge(vu, on="variation_id", how="left")
    log["clinvar_mapped_variation_ids"] = len(ok)
    # MONDO ancestors for re-pruning after the union
    parents = {m: (p.split(";") if isinstance(p, str) else []) for m, p in zip(mi.mondo, mi.parents)}
    memo = {}
    def anc(t):
        if t in memo:
            return memo[t]
        memo[t] = set(); out = set()
        for p in parents.get(t, []):
            out.add(p); out |= anc(p)
        memo[t] = out
        return out
    unit_mondo = dict(zip(du.disease_unit, du.mondo))

    def prune(units):
        units = sorted(set(units))
        drop = {u for u in units if unit_mondo.get(u) and any(
            w != u and unit_mondo.get(w) and unit_mondo[u] in anc(unit_mondo[w]) for w in units)}
        return [u for u in units if u not in drop]

    agg = []
    for k, g in ok.groupby(KEY, sort=False):
        units = prune([u for us in g.disease_units for u in (list(us) if us is not None else [])])
        le = pd.to_datetime(g.last_evaluated, errors="coerce", format="mixed")
        agg.append(dict(zip(KEY, k), gene=g.gene.iat[0], clinvar_variation_ids=";".join(map(str, g.variation_id)),
                        n_clinvar_records=len(g), clinvar_hgvs_p=" | ".join(dict.fromkeys(g.hgvs_p)),
                        clinvar_hgvs_c=" | ".join(dict.fromkeys(g.hgvs_c)),
                        clinvar_mapping_method=";".join(sorted(set(g.map_method))),
                        clinvar_classification=" | ".join(sorted(set(g.clinvar_classification))),
                        stars=int(g.stars.max()), review_status=" | ".join(sorted(set(g.review_status))),
                        n_submitters=int(pd.to_numeric(g.n_submitters).max()), n_submitters_sum=int(pd.to_numeric(g.n_submitters).sum()),
                        last_evaluated=(le.max().date().isoformat() if le.notna().any() else None),
                        rsids=";".join(sorted(set(g.rsid.astype(str)) - {"-1", "nan", "None"})) or None,
                        clinvar_units=units, clinvar_n_units=len(units),
                        clinvar_has_generic_condition=bool(g.has_generic_condition.any()),
                        clinvar_has_nonspecific_condition=bool(g.has_nonspecific_condition.any()),
                        clinvar_has_grouping_condition=bool(g.has_grouping_condition.any()),
                        clinvar_dropped_cancer_condition=bool(g.dropped_cancer_condition.any()),
                        clinvar_cancer_flag=bool(g.cancer_flag.any()),
                        clinvar_conditions_kept=" | ".join(dict.fromkeys(x for s in g.clinvar_condition_names_kept.fillna("") for x in s.split(" | ") if x)),
                        clinvar_conditions_all=" | ".join(dict.fromkeys(x for s in g.clinvar_condition_names_all.fillna("") for x in s.split(" | ") if x))))
    cva = pd.DataFrame(agg)
    log["clinvar_distinct_protein_variants"] = len(cva)
    log["clinvar_protein_variants_from_multiple_records"] = int((cva.n_clinvar_records > 1).sum())

    # ---- phenopacket-store ----
    pv = pd.read_parquet(DERIVED_PPKT / "variants.parquet")
    ind = pd.read_parquet(DERIVED_PPKT / "individuals.parquet")[["phenopacket_id", "disease_id", "disease_label", "pmid", "zygosity", "biallelic", "cohort"]]
    m = pv[pv.variant_class == "missense"].copy()
    log["ppkt_missense_rows"] = len(m)
    m = m.merge(ind.drop(columns=["cohort"]), on="phenopacket_id", how="left")
    m["mane_pos"] = pd.to_numeric(m.protein_pos, errors="coerce").astype("Int64")
    m["mane_ref"], m["mane_alt"] = m.vep_p_ref, m.vep_p_alt
    # MANE version check: ppkt VEP MANE ENSP vs MANE v1.5 ENSP of the same gene
    ensp = dict(zip(gm.hgnc_id, gm.mane_ensp.astype(str).str.split(".").str[0]))
    m["ppkt_ensp_is_mane_v15"] = [ensp.get(h) == (e.split(".")[0] if isinstance(e, str) else None) for h, e in zip(m.hgnc_id, m.ensembl_protein)]
    log["ppkt_missense_rows_ensp_not_mane_v15"] = int((~m.ppkt_ensp_is_mane_v15).sum())
    m["carrier_class"] = np.where(m.biallelic, "biallelic", np.where(m.zygosity == "hemi", "hemizygous", "single_allele_het"))
    m["single_allele"] = ~m.biallelic.astype(bool)
    m = m.dropna(subset=["mane_pos", "mane_ref", "mane_alt", "hgnc_id"])
    log["ppkt_missense_rows_keyed"] = len(m)
    log["ppkt_carrier_class_rows"] = {k: int(v) for k, v in m.carrier_class.value_counts().items()}
    prow = []
    for k, g in m.groupby(KEY, sort=False):
        sa, bi = g[g.single_allele], g[~g.single_allele]
        prow.append(dict(zip(KEY, k), ppkt_gene=g.gene.iat[0], ppkt_hgvs_p=" | ".join(dict.fromkeys(g.vep_hgvsp.dropna())),
                         ppkt_n_carriers=g.phenopacket_id.nunique(), ppkt_n_pmids=g.pmid.nunique(),
                         ppkt_n_carriers_single_allele=sa.phenopacket_id.nunique(), ppkt_n_pmids_single_allele=sa.pmid.nunique(),
                         ppkt_n_carriers_hemizygous=g[g.carrier_class == "hemizygous"].phenopacket_id.nunique(),
                         ppkt_n_carriers_biallelic=bi.phenopacket_id.nunique(), ppkt_n_pmids_biallelic=bi.pmid.nunique(),
                         ppkt_units_all=sorted(set(g.disease_id.dropna())),
                         ppkt_units_single_allele=sorted(set(sa.disease_id.dropna())),
                         ppkt_units_biallelic=sorted(set(bi.disease_id.dropna())),
                         ppkt_disease_labels=" | ".join(dict.fromkeys(g.disease_label.dropna())),
                         ppkt_cohorts=";".join(sorted(set(g.cohort.dropna()))),
                         ppkt_pmids=";".join(sorted(set(g.pmid.dropna())))))
    pa = pd.DataFrame(prow)
    log["ppkt_distinct_missense_variants"] = len(pa)
    log["ppkt_distinct_single_allele_variants"] = int((pa.ppkt_n_carriers_single_allele > 0).sum())
    log["ppkt_distinct_biallelic_variants"] = int((pa.ppkt_n_carriers_biallelic > 0).sum())
    log["ppkt_variants_both_single_and_biallelic_carriers"] = int(((pa.ppkt_n_carriers_single_allele > 0) & (pa.ppkt_n_carriers_biallelic > 0)).sum())

    # ---- merge ----
    for d in (cva, pa):
        d["mane_pos"] = d.mane_pos.astype("Int64")
    v = cva.merge(pa, on=KEY, how="outer", indicator=True)
    v["source"] = v._merge.map({"left_only": "clinvar", "right_only": "ppkt", "both": "both"}).astype(str)
    v = v.drop(columns="_merge")
    v["gene"] = v.gene.fillna(v.ppkt_gene)
    # gene symbol: MANE symbol from gene map
    gsym = dict(zip(gm.hgnc_id, gm.gene))
    v["gene"] = [gsym.get(h, g) for h, g in zip(v.hgnc_id, v.gene)]
    lst = lambda x: list(x) if isinstance(x, (list, np.ndarray)) else []
    v["clinvar_units"] = v.clinvar_units.map(lst)
    v["clinvar_n_units"] = v.clinvar_units.map(len)
    v["clinvar_single_label"] = v.clinvar_n_units == 1
    for c in ("ppkt_units_all", "ppkt_units_single_allele", "ppkt_units_biallelic"):
        v[c] = v[c].map(lst)
    v["labels_agree"] = [setcmp(a, b) for a, b in zip(v.clinvar_units, v.ppkt_units_all)]
    v["combined_units"] = [prune(list(a) + list(b)) for a, b in zip(v.clinvar_units, v.ppkt_units_all)]
    v["combined_n_units"] = v.combined_units.map(len)
    v["combined_single_label"] = v.combined_n_units == 1
    for c in [c for c in v.columns if c.startswith("ppkt_n_")]:
        v[c] = v[c].fillna(0).astype(int)
    v.to_parquet(INTER / "variants_merged.parquet", index=False)
    log["merged_variants"] = len(v)
    log["source"] = {k: int(x) for k, x in v.source.value_counts().items()}
    log["labels_agree_where_both_labelled"] = {str(k): int(x) for k, x in v.labels_agree.value_counts(dropna=False).items()}
    log["clinvar_single_label"] = int(v.clinvar_single_label.sum())
    log["combined_single_label"] = int(v.combined_single_label.sum())
    save_json(dict(log), INTER / "ppkt_merge_log.json")
    print(json.dumps(dict(log), indent=1))


if __name__ == "__main__":
    main()
