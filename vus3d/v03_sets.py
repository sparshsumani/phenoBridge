#!/usr/bin/env python3
"""Step 3: temporal train / test sets from the two ClinVar releases (no feature is looked at here).

Labelled record of a release = missense table of v02 (GRCh38, SNV, not somatic-only, Name parsed,
missense) with category PLP or BLB and review status >= 1 star (allelic/a02 STARS).

  TRAIN   labelled in 2023-01 (PLP -> P, BLB -> B), minus flipped variants.
  FLIPPED labelled in 2023-01 and the opposite category (BLB resp. PLP, any stars) in 2026-11; counted,
          saved to intermediate/flipped.parquet, excluded from train and test.
  TEST    labelled in 2026-11 and, in 2023-01, VUS, CONFLICT or ABSENT.
          2023 status is looked up by VariationID in the full 2023 GRCh38 table (all variant types, all
          stars, no missense/Name requirement); if the VariationID is absent it is looked up by GRCh38
          (chrom, pos, ref, alt) among 2023 SNV rows (status_2023_lookup = 'coordinate'); if neither -> ABSENT.
          Excluded (counted): 2023 category PLP/BLB with 0 stars or not in the 2023 missense table
          (prior_label_not_train), 2023 OTHER (e.g. "not provided", "risk factor", "Pathogenic; other").
  Reference: 2026 status of TRAIN variants (looked up the same way in the 2026 GRCh38 table).

Key agreement: for VariationIDs present in both releases the GRCh38 VCF coordinates are compared, and
for SNV coordinates present in both releases the VariationIDs are compared.
Outputs intermediate/sets.parquet (one row per train/test VariationID), intermediate/flipped.parquet,
        intermediate/excluded_test_candidates.parquet, intermediate/sets_log.json
"""
import sys, json
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from common import INTER, CUTOFF_DATE, save_json

STAT = ["variation_id", "type", "classification", "category", "stars", "review_status", "last_evaluated", "n_submitters", "coord_key"]
PROT = ["variation_id", "allele_id", "name", "name_tx", "name_tx_noversion", "name_gene", "clinvar_gene_symbols", "clinvar_hgnc_ids",
        "hgvs_c", "hgvs_p", "ref_aa", "clinvar_protein_pos", "alt_aa", "chrom", "pos_vcf", "ref_vcf", "alt_vcf", "coord_key", "rsid", "rcv",
        "origin_simple"]


def lookup(vids, coords, table, snv_coord_index):
    """status record of each (vid, coord) in a release: by VariationID, else by SNV coordinate."""
    t = table.set_index("variation_id")
    out = pd.DataFrame(index=range(len(vids)))
    hit = pd.Series(vids).isin(t.index).to_numpy()
    out["lookup"] = np.where(hit, "variation_id", None)
    out["matched_vid"] = np.where(hit, vids, None)
    miss = ~hit
    cvid = pd.Series(coords).map(snv_coord_index)
    use = miss & cvid.notna().to_numpy()
    out.loc[use, "lookup"] = "coordinate"
    out.loc[use, "matched_vid"] = cvid[use].to_numpy()
    out["lookup"] = out.lookup.fillna("absent")
    sub = t.reindex(out.matched_vid.where(out.matched_vid.notna(), "__none__"))
    for c in ["type", "classification", "category", "stars", "review_status", "last_evaluated", "n_submitters", "coord_key"]:
        out[c] = sub[c].to_numpy()
    out.loc[out.lookup == "absent", "category"] = "ABSENT"
    out.loc[out.lookup == "absent", "classification"] = "absent"
    return out


def main():
    log = {}
    g23 = pd.read_parquet(INTER / "clinvar_2023_grch38.parquet", columns=STAT)
    g26 = pd.read_parquet(INTER / "clinvar_2026_grch38.parquet", columns=STAT)
    m23 = pd.read_parquet(INTER / "clinvar_2023_missense.parquet")
    m26 = pd.read_parquet(INTER / "clinvar_2026_missense.parquet")

    # ---- key agreement ----
    both = g23[["variation_id", "coord_key", "type"]].merge(g26[["variation_id", "coord_key", "type"]], on="variation_id", suffixes=("_23", "_26"))
    ok = both.coord_key_23.notna() & both.coord_key_26.notna() & ~both.coord_key_23.str.contains("na|-", regex=True)
    snv = (both.type_23 == "single nucleotide variant") & (both.type_26 == "single nucleotide variant")
    ka = dict(vids_in_both=len(both), vids_only_2023=int((~g23.variation_id.isin(g26.variation_id)).sum()),
              vids_only_2026=int((~g26.variation_id.isin(g23.variation_id)).sum()),
              snv_vids_in_both=int(snv.sum()),
              snv_vid_same_coord=int((snv & (both.coord_key_23 == both.coord_key_26)).sum()),
              snv_vid_diff_coord=int((snv & (both.coord_key_23 != both.coord_key_26)).sum()),
              all_types_vid_same_coord=int((both.coord_key_23 == both.coord_key_26).sum()),
              type_changed=int((both.type_23 != both.type_26).sum()))
    ka["snv_vid_diff_coord_examples"] = both[snv & (both.coord_key_23 != both.coord_key_26)].head(10).to_dict("records")
    s23 = g23[g23.type == "single nucleotide variant"]; s26 = g26[g26.type == "single nucleotide variant"]
    c23 = s23.groupby("coord_key").variation_id.agg(list); c26 = s26.groupby("coord_key").variation_id.agg(list)
    ka["snv_coords_multi_vid_2023"] = int((c23.map(len) > 1).sum()); ka["snv_coords_multi_vid_2026"] = int((c26.map(len) > 1).sum())
    cu23 = c23[c23.map(len) == 1].str[0]; cu26 = c26[c26.map(len) == 1].str[0]
    cc = pd.concat([cu23.rename("v23"), cu26.rename("v26")], axis=1, join="inner")
    ka.update(snv_coords_in_both_unique=len(cc), snv_coord_same_vid=int((cc.v23 == cc.v26).sum()),
              snv_coord_diff_vid=int((cc.v23 != cc.v26).sum()))
    ka["snv_coord_diff_vid_examples"] = cc[cc.v23 != cc.v26].head(10).reset_index().to_dict("records")
    log["key_agreement"] = ka
    print(json.dumps(ka, indent=1, default=str), flush=True)

    # ---- labelled sets ----
    L23 = m23[m23.category.isin(["PLP", "BLB"]) & (m23.stars >= 1)].copy()
    L26 = m26[m26.category.isin(["PLP", "BLB"]) & (m26.stars >= 1)].copy()
    log["labelled_2023"] = {k: int(v) for k, v in L23.category.value_counts().items()}
    log["labelled_2026"] = {k: int(v) for k, v in L26.category.value_counts().items()}

    # train + 2026 status
    st26 = lookup(L23.variation_id.to_numpy(), L23.coord_key.to_numpy(), g26, cu26)
    tr = L23[PROT].reset_index(drop=True)
    tr["label"] = np.where(L23.category.to_numpy() == "PLP", "P", "B")
    for c, src in [("classification", "status_2023"), ("category", "category_2023"), ("stars", "stars_2023"), ("review_status", "review_status_2023"),
                   ("last_evaluated", "last_evaluated_2023"), ("n_submitters", "n_submitters_2023")]:
        tr[src] = L23[c].to_numpy()
    tr["status_2023_lookup"] = "variation_id"
    for c in ["classification", "category", "stars", "review_status", "last_evaluated", "n_submitters"]:
        tr[{"classification": "status_2026"}.get(c, c + "_2026")] = st26[c].to_numpy()
    tr["status_2026_lookup"] = st26.lookup.to_numpy()
    tr["vid_2026"] = st26.matched_vid.to_numpy()
    flip = ((tr.label == "P") & (tr.category_2026 == "BLB")) | ((tr.label == "B") & (tr.category_2026 == "PLP"))
    fl = tr[flip].copy()
    fl.to_parquet(INTER / "flipped.parquet", index=False)
    log["flipped"] = dict(total=int(flip.sum()), P_to_B=int(((tr.label == "P") & flip).sum()), B_to_P=int(((tr.label == "B") & flip).sum()),
                          with_2026_stars_ge1=int((flip & (tr.stars_2026 >= 1)).sum()))
    tr = tr[~flip].copy()
    tr["set"] = "train"
    log["train_2026_status"] = {f"{l}->{c}": int(n) for (l, c), n in tr.groupby(["label", "category_2026"]).size().items()}
    log["train_2026_lookup"] = {k: int(v) for k, v in tr.status_2026_lookup.value_counts().items()}

    # test
    st23 = lookup(L26.variation_id.to_numpy(), L26.coord_key.to_numpy(), g23, cu23)
    te = L26[PROT].reset_index(drop=True)
    te["label"] = np.where(L26.category.to_numpy() == "PLP", "P", "B")
    for c in ["classification", "category", "stars", "review_status", "last_evaluated", "n_submitters"]:
        te[{"classification": "status_2026"}.get(c, c + "_2026")] = L26[c].to_numpy()
    te["status_2026_lookup"] = "variation_id"
    te["vid_2026"] = te.variation_id
    for c in ["classification", "category", "stars", "review_status", "last_evaluated", "n_submitters", "type"]:
        te[{"classification": "status_2023"}.get(c, c + "_2023")] = st23[c].to_numpy()
    te["status_2023_lookup"] = st23.lookup.to_numpy()
    te["vid_2023"] = st23.matched_vid.to_numpy()
    in_train_vid = te.variation_id.isin(L23.variation_id) | te.vid_2023.isin(L23.variation_id)
    te["exclusion"] = None
    te.loc[in_train_vid, "exclusion"] = "labelled_2023 (train or flipped)"
    rest = te.exclusion.isna()
    te.loc[rest & te.category_2023.isin(["PLP", "BLB"]), "exclusion"] = "prior_label_not_train"
    te.loc[te.exclusion.isna() & (te.category_2023 == "OTHER"), "exclusion"] = "prior_OTHER"
    log["test_candidates_2026_labelled"] = len(te)
    log["test_candidate_exclusions"] = {str(k): int(v) for k, v in te.exclusion.value_counts(dropna=False).items()}
    pl = te[te.exclusion == "prior_label_not_train"]
    log["prior_label_not_train_detail"] = dict(
        stars0=int((pl.stars_2023 == 0).sum()), stars_ge1=int((pl.stars_2023 >= 1).sum()),
        by_category=f"{dict(pl.category_2023.value_counts())}",
        stars_ge1_type_2023=f"{dict(pl[pl.stars_2023 >= 1].type_2023.value_counts())}",
        stars_ge1_lookup=f"{dict(pl[pl.stars_2023 >= 1].status_2023_lookup.value_counts())}")
    log["prior_OTHER_classifications"] = {k: int(v) for k, v in te[te.exclusion == "prior_OTHER"].status_2023.value_counts().head(20).items()}
    exc = te[te.exclusion.notna() & (te.exclusion != "labelled_2023 (train or flipped)")]
    exc.to_parquet(INTER / "excluded_test_candidates.parquet", index=False)
    te = te[te.exclusion.isna()].drop(columns=["exclusion"]).copy()
    te["set"] = "test"
    te["prior_status_2023"] = te.category_2023
    log["test"] = {f"{l}|prior_{c}": int(n) for (l, c), n in te.groupby(["label", "prior_status_2023"]).size().items()}
    log["test_status_2023_lookup"] = {k: int(v) for k, v in te.status_2023_lookup.value_counts().items()}
    log["test_last_evaluated_2026_after_cutoff"] = dict(
        n=len(te), with_date=int(te.last_evaluated_2026.notna().sum()),
        after_cutoff=int((te.last_evaluated_2026 > pd.Timestamp(CUTOFF_DATE)).sum()),
        on_or_before_cutoff=int((te.last_evaluated_2026 <= pd.Timestamp(CUTOFF_DATE)).sum()))
    tr["prior_status_2023"] = tr.category_2023
    s = pd.concat([tr, te], ignore_index=True)
    # protein change as written in the other release (consistency check; final protein position comes from v04 mapping)
    o23 = m23.set_index("variation_id")[["name_gene", "hgvs_p", "name_tx"]]; o26 = m26.set_index("variation_id")[["name_gene", "hgvs_p", "name_tx"]]
    s["hgvs_p_2023_name"] = s.variation_id.map(o23.hgvs_p); s["hgvs_p_2026_name"] = s.variation_id.map(o26.hgvs_p)
    s["tx_2023_name"] = s.variation_id.map(o23.name_tx); s["tx_2026_name"] = s.variation_id.map(o26.name_tx)
    bothp = s.hgvs_p_2023_name.notna() & s.hgvs_p_2026_name.notna()
    log["name_protein_change_both_releases"] = dict(n=int(bothp.sum()), same_hgvs_p=int((bothp & (s.hgvs_p_2023_name == s.hgvs_p_2026_name)).sum()),
                                                   same_transcript=int((bothp & (s.tx_2023_name.str.split('.').str[0] == s.tx_2026_name.str.split('.').str[0])).sum()))
    s["label_release"] = np.where(s.set == "train", "2023", "2026")
    s.to_parquet(INTER / "sets.parquet", index=False)
    log["sets"] = {f"{a}|{b}": int(n) for (a, b), n in s.groupby(["set", "label"]).size().items()}
    save_json(log, INTER / "sets_log.json")
    print(json.dumps({k: v for k, v in log.items() if k != "key_agreement"}, indent=1, default=str))


if __name__ == "__main__":
    main()
