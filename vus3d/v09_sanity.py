#!/usr/bin/env python3
"""Step 9: sanity checks (descriptive only: no classifier, no P-vs-B comparison of any feature or score).

  * VariationID <-> GRCh38 coordinate agreement between the two releases (from v03)
  * reference-aa match rates (MANE seq / AF model / UniProt; residue table aa vs variant ref)
  * fraction of TEST variants whose 2026 aggregate LastEvaluated is after 2023-01-31 (overall and by 2023 status)
  * spot-check of 5 random TEST variants (2 prior ABSENT, 2 prior VUS, 1 prior CONFLICT; seed 11) against the
    live ClinVar VCV record (E-utilities efetch rettype=vcv, cached in cache/clinvar_vcv/): current aggregate
    germline classification vs our 2026 status, VCV DateCreated, and per-submission (SCV) DateCreated /
    DateLastUpdated / classification / DateLastEvaluated -> whether the record is consistent with its
    2023-01 status (ABSENT: VCV created after 2023-01-05; VUS/CONFLICT: a P/LP or B/LB SCV created or
    updated after 2023-01-05)
  * RWR: common.rwr_profile (vus3d) vs mechanistic_mvp.rwr on one newly computed feature protein, and region
    masses from the profile vs stored CoreImpact / DisorderImpact
Output intermediate/sanity.json, also written into derived_vus3d/counts.json['sanity'].
"""
import sys, json, re
import xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import common as C
from common import INTER, OUT, CUTOFF_DATE, save_json, load_json

RELEASE_2023 = pd.Timestamp("2023-01-05")


def vcv(vid):
    C.mg_import("s03_vep_annotate")
    mgc = sys.modules["_multigene_common"]
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id={vid}&rettype=vcv&is_variationid"
    def fetch():
        import requests
        for _ in range(5):
            r = requests.get(url, timeout=120)
            if r.status_code == 200:
                return {"xml": r.text}
        return None
    return mgc.cached_request("clinvar_vcv", url, fetch)


def parse_vcv(xml):
    root = ET.fromstring(xml)
    va = root.find(".//VariationArchive")
    gc = va.find("./ClassifiedRecord/Classifications/GermlineClassification")
    rec = dict(vcv=va.get("Accession"), version=va.get("Version"), name=va.get("VariationName"), date_created=va.get("DateCreated"),
               date_last_updated=va.get("DateLastUpdated"), aggregate=gc.findtext("Description") if gc is not None else None,
               aggregate_review_status=gc.findtext("ReviewStatus") if gc is not None else None,
               aggregate_last_evaluated=gc.get("DateLastEvaluated") if gc is not None else None)
    scvs = []
    for ca in va.findall(".//ClinicalAssertionList/ClinicalAssertion"):
        acc = ca.find("ClinVarAccession")
        cl = ca.find("Classification")
        scvs.append(dict(scv=acc.get("Accession") if acc is not None else None, submitter=acc.get("SubmitterName") if acc is not None else None,
                         date_created=ca.get("DateCreated"), date_last_updated=ca.get("DateLastUpdated"), submission_date=ca.get("SubmissionDate"),
                         classification=cl.findtext("GermlineClassification") if cl is not None else None,
                         date_last_evaluated=cl.get("DateLastEvaluated") if cl is not None else None,
                         contributes=ca.get("ContributesToAggregateClassification")))
    rec["submissions"] = scvs
    return rec


def labelled(c):
    return c in ("Pathogenic", "Likely pathogenic", "Benign", "Likely benign", "Pathogenic/Likely pathogenic", "Benign/Likely benign")


def main():
    v = pd.read_parquet(OUT / "variants.parquet")
    sl = load_json(INTER / "sets_log.json"); pm = load_json(INTER / "protein_mapping_log.json")
    s = {}
    s["key_agreement"] = {k: x for k, x in sl["key_agreement"].items() if not k.endswith("examples")}
    s["key_agreement"]["snv_vid_diff_coord_examples"] = sl["key_agreement"]["snv_vid_diff_coord_examples"]
    s["key_agreement"]["note"] = ("VariationID -> coordinate: identical for all SNV VariationIDs present in both releases except one with "
                                  "missing 2023 coordinates; coordinate -> VariationID: differing VariationIDs at the same SNV coordinate "
                                  "are 2023 records replaced by new VariationIDs (handled by the coordinate fallback lookup)")
    s["ref_match"] = dict(variants_on_single_fragment_af=pm["variants_on_single_fragment_af"], af_pos_mapped=pm["variants_af_pos_mapped"],
                          ref_matches_af=pm["ref_matches_af"], ref_mismatches_af=pm["ref_mismatches_af"],
                          rate_among_mapped=pm["ref_match_rate_af_among_mapped"], rate_uniprot=pm["ref_match_rate_uniprot"],
                          vep_validation=pm["vep_validation_on_mane_transcript_variants"],
                          residue_table_aa_equals_ref=f"{int((v.res_aa_matches_ref == True).sum())}/{int(v.res_aa_matches_ref.notna().sum())}")
    te = v[v.set == "test"]
    after = te.last_evaluated_2026 > pd.Timestamp(CUTOFF_DATE)
    s["test_last_evaluated_after_cutoff"] = dict(
        cutoff=CUTOFF_DATE, n=len(te), with_date=int(te.last_evaluated_2026.notna().sum()), after=int(after.sum()),
        frac_after_among_dated=float(after.sum() / te.last_evaluated_2026.notna().sum()),
        by_status_2023={k: dict(n=len(d), after=int((d.last_evaluated_2026 > pd.Timestamp(CUTOFF_DATE)).sum()),
                                no_date=int(d.last_evaluated_2026.isna().sum()),
                                year_counts={str(int(y)): int(n) for y, n in d.last_evaluated_2026.dt.year.value_counts().sort_index().items()})
                        for k, d in te.groupby("status_2023")})
    # spot check
    rng = np.random.RandomState(11)
    picks = []
    for st, n in (("ABSENT", 2), ("VUS", 2), ("CONFLICT", 1)):
        sub = te[te.status_2023 == st]
        picks += list(sub.iloc[rng.choice(len(sub), n, replace=False)].itertuples())
    sc = []
    for r in picks:
        res = vcv(r.variation_id)
        if not res:
            sc.append(dict(variation_id=r.variation_id, error="fetch_failed")); continue
        rec = parse_vcv(res["xml"])
        cat = lambda c: "PLP" if c in ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic") else (
            "BLB" if c in ("Benign", "Likely benign", "Benign/Likely benign") else ("VUS" if c == "Uncertain significance" else (
                "CONFLICT" if (c or "").startswith("Conflicting") else "OTHER")))
        lab_scv = [x for x in rec["submissions"] if labelled(x["classification"])]
        new_lab = [x for x in lab_scv if max(pd.Timestamp(x["date_created"] or "1900-01-01"), pd.Timestamp(x["date_last_updated"] or "1900-01-01")) > RELEASE_2023]
        if r.status_2023 == "ABSENT":
            consistent = pd.Timestamp(rec["date_created"]) > RELEASE_2023
        else:
            consistent = len(new_lab) > 0
        sc.append(dict(variation_id=r.variation_id, ours=dict(gene=r.gene, hgvs_p=r.hgvs_p, clinvar_hgvs_p=r.clinvar_hgvs_p, status_2023=r.status_2023,
                                                         status_2023_raw=r.status_2023_raw, status_2026_raw=r.status_2026_raw, stars_2026=int(r.stars_2026),
                                                         last_evaluated_2026=str(r.last_evaluated_2026.date()) if pd.notna(r.last_evaluated_2026) else None),
                       live=dict({k: rec[k] for k in rec if k != "submissions"}, n_submissions=len(rec["submissions"]),
                                 submissions=[{k: x[k] for k in ("scv", "classification", "date_created", "date_last_updated", "date_last_evaluated", "contributes")} for x in rec["submissions"]]),
                       title_contains_hgvs_p=r.clinvar_hgvs_p in (rec["name"] or ""),
                       live_category_equals_2026=cat(rec["aggregate"]) == r.status_2026,
                       history_consistent_with_2023_status=bool(consistent)))
    s["spot_check"] = sc
    s["spot_check_summary"] = dict(n=len(sc), title_ok=sum(x.get("title_contains_hgvs_p", False) for x in sc),
                                   live_category_ok=sum(x.get("live_category_equals_2026", False) for x in sc),
                                   history_ok=sum(x.get("history_consistent_with_2023_status", False) for x in sc))
    # RWR verification on a newly computed protein
    M = C.mg_import("structfeat").M
    meta = pd.read_parquet(INTER / "structure_meta.parquet")
    cand = meta[meta.source == "computed_vus3d"].sort_values("L")
    acc = cand.uniprot.iloc[len(cand) // 2]
    rf = C.residue_table(acc)
    pos, xyz, pl = C.ca_coords(acc)
    G = M.build_confidence_graph(xyz, pos, pl, cutoff_A=8.0, min_plddt_node=0.0)
    seeds = list(pos[np.linspace(0, len(pos) - 1, 5).astype(int)])
    md = mm = 0.0
    core = rf.region.eq("core").to_numpy(); dis_mask = rf.plddt.to_numpy() <= np.quantile(rf.plddt, 0.5)
    for a in (0.3, 0.15, 0.5):
        for sd in seeds:
            p = C.rwr_profile(acc, int(sd), alpha=a)
            ref = M.rwr(G, int(sd), alpha=a, max_iter=2000, tol=1e-13)
            md = max(md, float(np.abs(p.to_numpy() - np.array([ref[n] for n in pos])).max()))
            sfx = "" if a == 0.3 else f"_a{str(a).replace('0.', '')}"
            i = int(np.where(pos == sd)[0][0])
            mm = max(mm, abs(float(p.to_numpy()[core].sum()) - float(rf["CoreImpact" + sfx].iat[i])),
                     abs(float(p.to_numpy()[dis_mask].sum()) - float(rf["DisorderImpact" + sfx].iat[i])))
    pos2, W, A = C.contact_graph(acc)
    s["rwr_check"] = dict(uniprot=acc, L=int(len(pos)), seeds=[int(x) for x in seeds], max_abs_diff_vs_mechanistic_mvp_rwr=md,
                          max_abs_diff_region_mass_vs_table=mm, W_column_sums_max_dev=float(np.abs(np.asarray(W.sum(axis=0)).ravel()[np.asarray(A.sum(axis=0)).ravel() > 0] - 1).max()),
                          ok=bool(md < 1e-8 and mm < 1e-8))
    # all feature proteins have both files and consistent node counts
    bad = []
    for a in meta.uniprot:
        try:
            p_, W_, A_ = C.contact_graph(a)
            if W_.shape[0] != len(p_):
                bad.append(a)
        except Exception as e:
            bad.append(f"{a}:{e}")
        C.contact_graph.cache_clear(); C.residue_table.cache_clear()
    s["graph_files_check"] = dict(proteins=len(meta), failures=bad[:20], n_failures=len(bad))
    save_json(s, INTER / "sanity.json")
    cj = load_json(OUT / "counts.json"); cj["sanity"] = s; save_json(cj, OUT / "counts.json")
    print(json.dumps({k: x for k, x in s.items() if k != "spot_check"}, indent=1, default=str))
    for x in sc:
        print(json.dumps({k: x[k] for k in x if k != "live"}, default=str))


if __name__ == "__main__":
    main()
