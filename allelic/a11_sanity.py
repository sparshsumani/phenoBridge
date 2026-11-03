#!/usr/bin/env python3
"""Step 11: sanity checks (no association / no feature-by-disease comparison).

  * reference-aa match rates (MANE seq, AF model, UniProt), VEP validation agreement
  * ClinVar rows / VariationIDs dropped at each filter
  * top 20 eligible genes by single-label variant count with their disease units (>=5) and counts
  * presence/status of known allelic genes
  * spot-check of 5 random single-label variants of eligible genes against the live ClinVar record
    (NCBI E-utilities esummary db=clinvar, cached in cache/clinvar_esummary/): title/protein change,
    germline classification, review status and trait names vs our table
  * RWR: common.rwr_profile vs mechanistic_mvp.rwr (iterative, repo code) on one feature protein, and
    region masses recomputed from rwr_profile vs the stored CoreImpact/DisorderImpact columns
Output intermediate/sanity.json, also written into derived_allelic/counts.json['sanity'].
"""
import sys, json, random
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import common as C
from common import INTER, OUT, save_json, load_json

KNOWN = ["LMNA", "PTPN11", "FGFR3", "RET", "RYR1", "GJB2", "KCNQ1", "SCN5A", "TBX5", "MYH7", "KRAS", "HRAS", "COL2A1",
         "FGFR2", "PRPH2", "FGFR1", "SCN1A", "TP53", "PIK3CA", "CACNA1A", "KCNH2", "TGFBR2", "COL1A1", "NOTCH3", "TTR", "GNAS"]


def esummary(vid):
    mg = C.mg_import("s03_vep_annotate")  # ensures multigene common is loaded
    mgc = sys.modules["_multigene_common"]
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={vid}&retmode=json"
    return mgc.cached_request("clinvar_esummary", url, lambda: mgc.http_json(url))


def main():
    v = pd.read_parquet(OUT / "variants.parquet")
    g = pd.read_parquet(OUT / "genes.parquet")
    gd = pd.read_parquet(OUT / "gene_disease_counts.parquet")
    d = pd.read_parquet(OUT / "diseases.parquet")
    pm = load_json(INTER / "protein_mapping_log.json")
    cf = load_json(INTER / "clinvar_filter_log.json")
    s = {}
    s["ref_match"] = dict(
        clinvar_variants_on_single_fragment_af=pm["variants_on_single_fragment_af"], af_pos_mapped=pm["variants_af_pos_mapped"],
        ref_matches_af=pm["ref_matches_af"], ref_mismatches_af=pm["ref_mismatches_af"],
        rate_among_mapped=pm["ref_match_rate_af"],
        rate_among_all_on_single_fragment=pm["ref_matches_af"] / pm["variants_on_single_fragment_af"],
        merged_variants_ref_matches_af=f"{int((v.ref_matches_af == True).sum())}/{int(v.ref_matches_af.notna().sum())}",
        merged_variants_ref_matches_uniprot=f"{int((v.ref_matches_uniprot == True).sum())}/{int(v.ref_matches_uniprot.notna().sum())}",
        vep_validation=pm["vep_validation_on_mane_transcript_variants"],
        mismatching_variants=v[v.ref_matches_af == False][["gene", "hgvs_p", "af_pos", "source"]].to_dict("records"))
    s["clinvar_dropped_variation_ids"] = cf["dropped_variation_ids"]
    s["clinvar_variation_ids_after"] = cf["variation_ids_after"]
    el = g[g.eligible].sort_values("n_single_label", ascending=False)
    prim = gd[gd.definition == "clinvar"]
    top = []
    for r in el.head(20).itertuples():
        sub = prim[(prim.hgnc_id == r.hgnc_id) & (prim.n_variants_countable >= 5)].sort_values("n_variants_countable", ascending=False)
        top.append(dict(gene=r.gene, uniprot=r.uniprot, n_single_label=int(r.n_single_label), n_disease_units=int(r.n_disease_units),
                        units_ge5=[f"{u} {n} [{inh}] {nm}" for u, n, inh, nm in zip(sub.disease_unit, sub.n_variants_countable, sub.inheritance_clean, sub.disease_name)]))
    s["top20_eligible"] = top
    kn = {}
    for gene in KNOWN:
        r = g[g.gene == gene]
        if not len(r):
            kn[gene] = "no P/LP >=1-star missense in ClinVar/ppkt"; continue
        r = r.iloc[0]
        sub = prim[(prim.hgnc_id == r.hgnc_id)].sort_values("n_variants_countable", ascending=False)
        kn[gene] = dict(eligible=bool(r.eligible), eligible_combined=bool(r.eligible_combined), ad_ar_eligible=bool(r.ad_ar_eligible),
                        af_n_fragments=int(r.af_n_fragments), n_single_label=int(r.n_single_label), n_units_ge5=int(r.n_disease_units_ge5),
                        units=[f"{u}:{n}" for u, n in zip(sub.disease_unit.head(6), sub.n_variants_countable.head(6))])
    s["known_allelic_genes"] = kn

    # spot-check
    random.seed(7)
    cand = v[v.gene_eligible & v.single_label & v.clinvar_variation_ids.notna() & (v.n_clinvar_records == 1)]
    picks = cand.sample(5, random_state=7)
    sc = []
    for r in picks.itertuples():
        vid = r.clinvar_variation_ids
        res = esummary(vid)
        rec = (res or {}).get("result", {}).get(str(vid), {})
        gc = rec.get("germline_classification", {})
        traits = [t.get("trait_name") for t in gc.get("trait_set", [])]
        xrefs = [x for t in gc.get("trait_set", []) for x in t.get("trait_xrefs", [])]
        omims = sorted({"OMIM:" + x["db_id"] for x in xrefs if x.get("db_source") == "OMIM" and x.get("db_id", "").isdigit()})
        sc.append(dict(variation_id=vid, ours=dict(gene=r.gene, hgnc_p=r.hgvs_p, clinvar_hgvs_p=r.clinvar_hgvs_p,
                                                   classification=r.clinvar_classification, stars=int(r.stars),
                                                   unit=r.disease_unit, label=r.disease_label),
                       clinvar_live=dict(title=rec.get("title"), classification=gc.get("description"), review_status=gc.get("review_status"),
                                         traits=traits, trait_omims=omims),
                       title_matches=bool(rec.get("title")) and (r.clinvar_hgvs_p.split(" | ")[0] in rec.get("title", "")),
                       classification_matches=gc.get("description") == r.clinvar_classification,
                       unit_in_live_traits=r.disease_unit in omims or r.disease_unit.startswith("MONDO:")))
    s["spot_check"] = sc

    # RWR verification on one feature protein
    M = C.mg_import("structfeat").M
    acc = g[g.has_residue_features].sort_values("uniprot_length").uniprot.iloc[len(g[g.has_residue_features]) // 2]
    rf = pd.read_parquet(C.RFDIR / f"{acc}.parquet")
    pos, xyz, pl = C.ca_coords(acc)
    G = M.build_confidence_graph(xyz, pos, pl, cutoff_A=8.0, min_plddt_node=0.0)
    seeds = list(pos[np.linspace(0, len(pos) - 1, 5).astype(int)])
    md, mm = 0.0, 0.0
    core = rf.region.eq("core").to_numpy()
    dis_mask = rf.plddt.to_numpy() <= np.quantile(rf.plddt, 0.5)
    for a in (0.3, 0.15, 0.5):
        for sd in seeds:
            p = C.rwr_profile(acc, int(sd), alpha=a)
            ref = M.rwr(G, int(sd), alpha=a, max_iter=2000, tol=1e-13)
            md = max(md, float(np.abs(p.to_numpy() - np.array([ref[n] for n in pos])).max()))
            sfx = "" if a == 0.3 else f"_a{str(a).replace('0.', '')}"
            i = int(np.where(pos == sd)[0][0])
            mm = max(mm, abs(float(p.to_numpy()[core].sum()) - float(rf["CoreImpact" + sfx].iat[i])),
                     abs(float(p.to_numpy()[dis_mask].sum()) - float(rf["DisorderImpact" + sfx].iat[i])))
    s["rwr_check"] = dict(uniprot=acc, seeds=[int(x) for x in seeds], max_abs_diff_vs_mechanistic_mvp_rwr=md,
                          max_abs_diff_region_mass_vs_table=mm, ok=bool(md < 1e-8 and mm < 1e-8))
    save_json(s, INTER / "sanity.json")
    cj = load_json(OUT / "counts.json")
    cj["sanity"] = s
    save_json(cj, OUT / "counts.json")
    print(json.dumps(s, indent=1, default=str))


if __name__ == "__main__":
    main()
