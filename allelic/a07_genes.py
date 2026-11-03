#!/usr/bin/env python3
"""Step 7: position mapping of merged variants onto AlphaFold/UniProt numbering, and per-gene
eligibility COUNTS (no feature is looked at here).

Variant usable for counting ("countable") when: its gene has a single-fragment AlphaFold v6 model
(proteins > 2,700 aa are fragmented in AFDB and excluded), the MANE position maps onto the model
(identity or trusted alignment block), and the MANE reference aa equals the model residue.

Per gene and label definition (each variant contributes to a disease unit only if it is single-label
under that definition):
  clinvar   ClinVar disease units only (PRIMARY; a03 rules).           -> n_*  / eligible
  combined  ClinVar units U phenopacket-store carriers' OMIM ids.      -> comb_* / eligible_combined
  omim      PRIMARY restricted to OMIM-phenotype units (MONDO-only units ignored; a variant whose
            only unit is MONDO is not counted)                         -> omim_* / eligible_omim
  eligible            >= 2 disease units each with >= 5 distinct single-label missense variants
  eligible_nonbroad   PRIMARY eligibility not counting 'broad' umbrella units (diseases.broad_unit)
  ad_ar_eligible      >= 1 unit with inheritance_clean AD and >= 5 variants AND >= 1 unit with
                      inheritance_clean AR and >= 5 variants (PRIMARY counts)
Outputs intermediate/variants_positioned.parquet, intermediate/genes.parquet,
        intermediate/gene_disease_counts.parquet, intermediate/genes_log.json
"""
import sys, json
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from common import INTER, MIN_VARIANTS_PER_DISEASE as K, save_json, mg_import

s06 = mg_import("s06_select_genes")


def main():
    v = pd.read_parquet(INTER / "variants_merged.parquet")
    gm = pd.read_parquet(INTER / "gene_protein_map.parquet")
    d = pd.read_parquet(INTER / "diseases.parquet")
    inh = dict(zip(d.disease_unit, d.inheritance_clean))
    dname = dict(zip(d.disease_unit, d.disease_name))
    gi = gm.drop_duplicates("hgnc_id").set_index("hgnc_id")
    maps = {}
    for h, r in gi.iterrows():
        for tag, tgt in (("af", r.af_seq), ("up", r.uniprot_seq)):
            if isinstance(r.mane_protein_seq, str) and isinstance(tgt, str):
                maps[(h, tag)] = None if r.mane_protein_seq == tgt else s06.align_map(r.mane_protein_seq, tgt)[0]

    def conv(h, p, tag):
        if pd.isna(p) or (h, tag) not in maps:
            return None
        m = maps[(h, tag)]
        return int(p) if m is None else m.get(int(p))
    v["uniprot"] = v.hgnc_id.map(gi.uniprot)
    v["af_pos"] = pd.array([conv(h, p, "af") for h, p in zip(v.hgnc_id, v.mane_pos)], dtype="Int64")
    v["uniprot_pos"] = pd.array([conv(h, p, "up") for h, p in zip(v.hgnc_id, v.mane_pos)], dtype="Int64")
    v["af_n_fragments"] = v.hgnc_id.map(gi.af_n_fragments)
    v["af_single_fragment"] = v.af_n_fragments == 1
    afs, ups = gi.af_seq.to_dict(), gi.uniprot_seq.to_dict()
    def chk(s, p, ref):
        if not isinstance(s, str) or pd.isna(p):
            return None
        return bool(1 <= int(p) <= len(s) and s[int(p) - 1] == ref)
    v["ref_matches_af"] = pd.array([chk(afs.get(h), p, r) for h, p, r in zip(v.hgnc_id, v.af_pos, v.mane_ref)], dtype="boolean")
    v["ref_matches_uniprot"] = pd.array([chk(ups.get(h), p, r) for h, p, r in zip(v.hgnc_id, v.uniprot_pos, v.mane_ref)], dtype="boolean")
    v["protein_length_af"] = v.hgnc_id.map(lambda h: len(afs[h]) if isinstance(afs.get(h), str) else np.nan)
    v["countable"] = (v.af_single_fragment & v.af_pos.notna() & (v.ref_matches_af == True)).astype(bool)
    v["omim_units"] = v.clinvar_units.map(lambda us: [u for u in us if u.startswith("OMIM:")])
    v["omim_single_label"] = v.clinvar_single_label & v.omim_units.map(len).eq(1)
    v.to_parquet(INTER / "variants_positioned.parquet", index=False)

    defs = {"clinvar": ("clinvar_units", "clinvar_single_label"),
            "combined": ("combined_units", "combined_single_label"),
            "omim": ("omim_units", "omim_single_label")}
    gd_rows = []
    for name, (ucol, scol) in defs.items():
        s = v[v[scol]]
        for (h, cnt), sub in s.assign(u=s[ucol].map(lambda x: x[0])).groupby(["hgnc_id", "countable"]):
            for u, n in sub.u.value_counts().items():
                gd_rows.append(dict(definition=name, hgnc_id=h, countable=cnt, disease_unit=u, n_variants=int(n)))
    gd = pd.DataFrame(gd_rows)
    gd = gd.pivot_table(index=["definition", "hgnc_id", "disease_unit"], columns="countable", values="n_variants",
                        aggfunc="sum", fill_value=0).reset_index()
    gd.columns.name = None
    gd = gd.rename(columns={True: "n_variants_countable", False: "n_variants_not_countable"})
    for c in ("n_variants_countable", "n_variants_not_countable"):
        if c not in gd:
            gd[c] = 0
    gd["n_variants_any"] = gd.n_variants_countable + gd.n_variants_not_countable
    gd["gene"] = gd.hgnc_id.map(gi.gene)
    gd["disease_name"] = gd.disease_unit.map(dname)
    gd["inheritance_clean"] = gd.disease_unit.map(inh)
    gd = gd.sort_values(["definition", "gene", "n_variants_countable"], ascending=[True, True, False])
    gd.to_parquet(INTER / "gene_disease_counts.parquet", index=False)

    g = gi.reset_index()[["gene", "hgnc_id", "uniprot", "uniprot_entry", "uniprot_map", "uniprot_length", "mane_enst", "mane_np",
                          "mane_ensp", "af_n_fragments", "af_source", "af_single_fragment", "af_pdb", "mane_eq_af", "mane_af_map",
                          "mane_af_identity", "mane_eq_uniprot", "af_eq_uniprot"]].copy()
    vc = v.groupby("hgnc_id")
    g["n_variants_all_sources"] = g.hgnc_id.map(vc.size()).fillna(0).astype(int)
    g["n_variants_clinvar"] = g.hgnc_id.map(v[v.source != "ppkt"].groupby("hgnc_id").size()).fillna(0).astype(int)
    g["n_variants_ppkt"] = g.hgnc_id.map(v[v.source != "clinvar"].groupby("hgnc_id").size()).fillna(0).astype(int)
    g["n_variants_countable"] = g.hgnc_id.map(v[v.countable].groupby("hgnc_id").size()).fillna(0).astype(int)
    g["n_ref_mismatch_af"] = g.hgnc_id.map(v[v.ref_matches_af == False].groupby("hgnc_id").size()).fillna(0).astype(int)
    pre = {"clinvar": "", "combined": "comb_", "omim": "omim_"}
    for name, (ucol, scol) in defs.items():
        p = pre[name]
        sub = gd[gd.definition == name]
        cs = v[v[scol] & v.countable].groupby("hgnc_id").size()
        g[f"{p}n_single_label"] = g.hgnc_id.map(cs).fillna(0).astype(int)
        g[f"{p}n_multi_label"] = g.hgnc_id.map(v[(v[ucol].map(len) > 1) & v.countable].groupby("hgnc_id").size()).fillna(0).astype(int)
        cnt = sub[sub.n_variants_countable > 0]
        g[f"{p}n_disease_units"] = g.hgnc_id.map(cnt.groupby("hgnc_id").size()).fillna(0).astype(int)
        ge = sub[sub.n_variants_countable >= K]
        g[f"{p}n_disease_units_ge{K}"] = g.hgnc_id.map(ge.groupby("hgnc_id").size()).fillna(0).astype(int)
        g[f"{p}disease_units_ge{K}"] = g.hgnc_id.map(ge.groupby("hgnc_id").apply(
            lambda x: ";".join(f"{u}({n})" for u, n in zip(x.disease_unit, x.n_variants_countable)))).fillna("")
        # same ignoring AF fragmentation/mapping (to show what the AF filter costs)
        ge_any = sub[sub.n_variants_any >= K]
        g[f"{p}n_disease_units_ge{K}_ignoring_af"] = g.hgnc_id.map(ge_any.groupby("hgnc_id").size()).fillna(0).astype(int)
        g[f"eligible{'' if name == 'clinvar' else '_' + name}"] = g[f"{p}n_disease_units_ge{K}"] >= 2
    prim = gd[(gd.definition == "clinvar") & (gd.n_variants_countable >= K)]
    broad = set(d[d.broad_unit].disease_unit)
    nb = prim[~prim.disease_unit.isin(broad)]
    g["n_disease_units_ge5_nonbroad"] = g.hgnc_id.map(nb.groupby("hgnc_id").size()).fillna(0).astype(int)
    g["eligible_nonbroad"] = g.n_disease_units_ge5_nonbroad >= 2
    ad = prim[prim.inheritance_clean == "AD"].groupby("hgnc_id").size()
    ar = prim[prim.inheritance_clean == "AR"].groupby("hgnc_id").size()
    g["n_AD_units_ge5"] = g.hgnc_id.map(ad).fillna(0).astype(int)
    g["n_AR_units_ge5"] = g.hgnc_id.map(ar).fillna(0).astype(int)
    g["ad_ar_eligible"] = (g.n_AD_units_ge5 >= 1) & (g.n_AR_units_ge5 >= 1)
    pc = gd[(gd.definition == "combined") & (gd.n_variants_countable >= K)]
    g["ad_ar_eligible_combined"] = g.hgnc_id.isin(set(pc[pc.inheritance_clean == "AD"].hgnc_id) & set(pc[pc.inheritance_clean == "AR"].hgnc_id))
    g["eligible_any_definition"] = g.eligible | g.eligible_combined | g.eligible_omim | g.ad_ar_eligible | g.ad_ar_eligible_combined
    g["fragmented_would_be_eligible"] = (g.af_n_fragments > 1) & (g[f"n_disease_units_ge{K}_ignoring_af"] >= 2)
    g = g.sort_values("n_single_label", ascending=False)
    g.to_parquet(INTER / "genes.parquet", index=False)
    el = g[g.eligible]
    log = dict(
        genes=len(g), genes_with_variants=int((g.n_variants_all_sources > 0).sum()),
        genes_af_single_fragment=int(g.af_single_fragment.sum()), genes_af_fragmented=int((g.af_n_fragments > 1).sum()),
        genes_af_missing=int((g.af_n_fragments == 0).sum()),
        fragmented_genes=sorted(g[g.af_n_fragments > 1].gene), fragmented_genes_that_would_be_eligible=sorted(g[g.fragmented_would_be_eligible].gene),
        variants_merged=len(v), variants_countable=int(v.countable.sum()),
        variants_clinvar_single_label_countable=int((v.clinvar_single_label & v.countable).sum()),
        variants_combined_single_label_countable=int((v.combined_single_label & v.countable).sum()),
        genes_ge1_unit_ge5=int((g.n_disease_units_ge5 >= 1).sum()),
        eligible=int(g.eligible.sum()), eligible_combined=int(g.eligible_combined.sum()), eligible_omim=int(g.eligible_omim.sum()), eligible_nonbroad=int(g.eligible_nonbroad.sum()),
        ad_ar_eligible=int(g.ad_ar_eligible.sum()), ad_ar_eligible_combined=int(g.ad_ar_eligible_combined.sum()),
        eligible_any_definition=int(g.eligible_any_definition.sum()),
        eligible_genes=el.gene.tolist(), eligible_combined_genes=g[g.eligible_combined].gene.tolist(),
        eligible_omim_genes=g[g.eligible_omim].gene.tolist(),
        ad_ar_eligible_genes=g[g.ad_ar_eligible].gene.tolist(), ad_ar_eligible_combined_genes=g[g.ad_ar_eligible_combined].gene.tolist(),
        eligible_variants_single_label=int(el.n_single_label.sum()),
        eligible_variants_in_units_ge5=int(prim[prim.hgnc_id.isin(el.hgnc_id)].n_variants_countable.sum()),
        eligible_units_ge5=int(prim[prim.hgnc_id.isin(el.hgnc_id)].shape[0]),
        distribution_n_units_ge5={int(k): int(x) for k, x in g.n_disease_units_ge5.value_counts().sort_index().items()})
    save_json(log, INTER / "genes_log.json")
    print(json.dumps({k: x for k, x in log.items() if not k.endswith("_genes") or k.startswith("ad_ar") or k.startswith("frag")}, indent=1))
    print(el[["gene", "uniprot", "uniprot_length", "n_single_label", "n_disease_units", "n_disease_units_ge5", "disease_units_ge5"]].head(40).to_string())


if __name__ == "__main__":
    main()
