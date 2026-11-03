#!/usr/bin/env python3
"""Step 5: per-gene eligibility COUNTS (descriptive only: no feature or score is looked at, and test labels
are only counted, never compared with any train-derived quantity).

Counted variants = 'usable' (v04: mapped to MANE, single-fragment AlphaFold model, AF position mapped,
reference aa matches the model).  A residue = (uniprot, af_pos).
Per gene:
  n_train_P, n_train_B, n_test_P, n_test_B        variants (VariationIDs)
  n_train_residues                                distinct residues carrying >=1 TRAIN label ("TRAIN labelled residues")
  n_train_residues_P / _B / _both                 residues with >=1 train P / >=1 train B / both
  n_test_residues
  n_test_same_residue_as_train                    test variants whose residue carries a train label (PS1/PM5 situations)
  n_test_same_protein_change_as_train             ... with the identical amino-acid substitution in train (PS1-like)
  n_test_same_residue_other_alt                   ... residue labelled in train, but no train variant with the same alt (PM5-like)
  evaluable_k (k = 5, 10, 20)                     >=1 TEST variant and >= k TRAIN labelled residues
  feature_protein                                 evaluable_5 (-> contact graph + residue features in v07)
Outputs intermediate/genes.parquet, intermediate/variants_usable.parquet (with residue-overlap flags), intermediate/genes_log.json
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from common import INTER, THRESHOLDS, MIN_TRAIN_RESIDUES, save_json


def main():
    s = pd.read_parquet(INTER / "sets.parquet")
    mp = pd.read_parquet(INTER / "variants_mapped.parquet")
    gm = pd.read_parquet(INTER / "gene_protein_map.parquet")
    v = s.merge(mp.drop(columns=["set", "chrom", "pos_vcf", "ref_vcf", "alt_vcf", "name_tx", "clinvar_protein_pos", "ref_aa", "alt_aa",
                                 "hgvs_p", "hgvs_c"]), on="variation_id", suffixes=("", "_map"))
    v = v.drop(columns=[c for c in v.columns if c.endswith("_map")])
    u = v[v.usable].copy()
    u["af_pos"] = u.af_pos.astype(int)
    u["residue_key"] = u.uniprot + ":" + u.af_pos.astype(str)
    u["protein_change_key"] = u.residue_key + ":" + u.mane_alt
    tr, te = u[u.set == "train"], u[u.set == "test"]
    tr_res = set(tr.residue_key); tr_pc = set(tr.protein_change_key)
    u["test_residue_has_train_label"] = np.where(u.set == "test", u.residue_key.isin(tr_res), None)
    u["test_same_protein_change_in_train"] = np.where(u.set == "test", u.protein_change_key.isin(tr_pc), None)
    u["n_train_variants_same_residue"] = u.residue_key.map(tr.residue_key.value_counts()).fillna(0).astype(int)
    # duplicates: several VariationIDs (nucleotide changes) -> same protein change within the same set
    u["n_vids_same_protein_change_same_set"] = u.groupby(["set", "protein_change_key"]).variation_id.transform("size")
    u.to_parquet(INTER / "variants_usable.parquet", index=False)

    te = u[u.set == "test"]
    g = u.groupby(["gene", "uniprot"]).agg(hgnc_id=("hgnc_id", "first")).reset_index()
    def cnt(df, name):
        return df.groupby("gene").size().rename(name)
    parts = [cnt(tr[tr.label == "P"], "n_train_P"), cnt(tr[tr.label == "B"], "n_train_B"),
             cnt(te[te.label == "P"], "n_test_P"), cnt(te[te.label == "B"], "n_test_B"),
             tr.groupby("gene").residue_key.nunique().rename("n_train_residues"),
             tr[tr.label == "P"].groupby("gene").residue_key.nunique().rename("n_train_residues_P"),
             tr[tr.label == "B"].groupby("gene").residue_key.nunique().rename("n_train_residues_B"),
             te.groupby("gene").residue_key.nunique().rename("n_test_residues"),
             te[te.test_residue_has_train_label == True].groupby("gene").size().rename("n_test_same_residue_as_train"),
             te[te.test_same_protein_change_in_train == True].groupby("gene").size().rename("n_test_same_protein_change_as_train")]
    g = g.set_index("gene").join(parts).fillna(0)
    for c in [c for c in g.columns if c.startswith("n_")]:
        g[c] = g[c].astype(int)
    g["n_train_residues_both"] = g.n_train_residues_P + g.n_train_residues_B - g.n_train_residues
    g["n_test_same_residue_other_alt"] = g.n_test_same_residue_as_train - g.n_test_same_protein_change_as_train
    g["n_train"] = g.n_train_P + g.n_train_B; g["n_test"] = g.n_test_P + g.n_test_B
    for k in THRESHOLDS:
        g[f"evaluable_{k}"] = (g.n_test >= 1) & (g.n_train_residues >= k)
    g["feature_protein"] = (g.n_test >= 1) & (g.n_train_residues >= MIN_TRAIN_RESIDUES)
    g["test_has_both_labels"] = (g.n_test_P >= 1) & (g.n_test_B >= 1)
    gi = gm.set_index("gene")
    for c in ["mane_enst", "mane_np", "uniprot_length", "af_pdb", "af_n_fragments", "mane_af_map", "mane_eq_af", "uniprot_map"]:
        g[c] = g.index.map(gi[c])
    g = g.reset_index()
    g.to_parquet(INTER / "genes.parquet", index=False)

    log = dict(usable_variants={f"{a}|{b}": int(n) for (a, b), n in u.groupby(["set", "label"]).size().items()},
               genes_with_usable=len(g), genes_with_train=int((g.n_train > 0).sum()), genes_with_test=int((g.n_test > 0).sum()),
               genes_test_has_both_labels=int(g.test_has_both_labels.sum()))
    for k in THRESHOLDS:
        e = g[g[f"evaluable_{k}"]]
        log[f"evaluable_{k}"] = dict(genes=len(e), test_variants=int(e.n_test.sum()), test_P=int(e.n_test_P.sum()), test_B=int(e.n_test_B.sum()),
                                     genes_test_both_labels=int(e.test_has_both_labels.sum()),
                                     test_variants_in_genes_with_both_labels=int(e[e.test_has_both_labels].n_test.sum()))
    log["test_same_residue"] = dict(
        test_total=len(te), same_residue_as_train=int((te.test_residue_has_train_label == True).sum()),
        same_protein_change_as_train=int((te.test_same_protein_change_in_train == True).sum()),
        same_residue_other_alt=int(((te.test_residue_has_train_label == True) & (te.test_same_protein_change_in_train != True)).sum()),
        by_test_label={l: dict(n=int(len(x)), same_residue=int((x.test_residue_has_train_label == True).sum()),
                               same_protein_change=int((x.test_same_protein_change_in_train == True).sum()))
                       for l, x in te.groupby("label")})
    log["duplicate_protein_changes"] = dict(
        train_protein_changes=int(tr.protein_change_key.nunique()), train_vids=len(tr),
        test_protein_changes=int(te.protein_change_key.nunique()), test_vids=len(te),
        protein_changes_with_both_P_and_B_in_train=int((tr.groupby("protein_change_key").label.nunique() > 1).sum()),
        protein_changes_with_both_P_and_B_in_test=int((te.groupby("protein_change_key").label.nunique() > 1).sum()))
    log["feature_proteins"] = int(g.feature_protein.sum())
    log["top20_genes_by_test"] = g.sort_values("n_test", ascending=False).head(20)[["gene", "n_train_P", "n_train_B", "n_train_residues", "n_test_P", "n_test_B"]].to_dict("records")
    save_json(log, INTER / "genes_log.json")
    print(json.dumps(log, indent=1, default=str))


if __name__ == "__main__":
    main()
