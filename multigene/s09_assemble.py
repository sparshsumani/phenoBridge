#!/usr/bin/env python3
"""Step 9: assemble final analysis tables in DERIVED/ and the counts used by README.md.

  individuals.parquet            one row per phenopacket (all 10,377), gene_selected flag
  variants.parquet               one row per phenopacket x causative variant (genomicInterpretation),
                                 VEP/NMD/AlphaMissense annotations + residue features at protein_pos
  phenotypes_explicit.parquet    explicit observed(1)/excluded(0) calls only
  phenotypes_propagated.parquet  true-path propagated calls (see s02)
  genes.parquet                  one row per gene
  residue_features/<UNIPROT>.parquet  (written by s07)
  counts.json                    every count quoted in README.md
"""
import json, gzip, shutil
from collections import Counter
import numpy as np
import pandas as pd
from common import INTER, DERIVED, MIN_INDIVIDUALS

RFDIR = DERIVED / "residue_features"
RF_COLS_DROP = {"gene", "uniprot", "pos"}


def vcd(vc):
    return {str(k): int(x) for k, x in vc.items()}


def mean_plddt(path):
    vals = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                vals.append(float(line[60:66]))
    v = np.array(vals)
    return v.mean(), (v < 50).mean(), (v >= 70).mean()


def main():
    C = {}
    ind = pd.read_parquet(INTER / "individuals_raw.parquet")
    v = pd.read_parquet(INTER / "variants_annot.parquet")
    g = pd.read_parquet(INTER / "gene_summary.parquet")
    am = pd.read_parquet(INTER / "alphamissense_subset.parquet")
    smeta = pd.read_parquet(INTER / "structure_meta.parquet")
    ph = pd.read_parquet(INTER / "phenotypes_explicit.parquet")
    pp = pd.read_parquet(INTER / "phenotypes_propagated.parquet")

    # ---- AlphaMissense match for missense variants ----
    # AlphaMissense (2023) indexes the UniProt canonical sequence it was built on, which for a few
    # proteins differs from today's UniProt/MANE.  Its reference sequence is reconstructed from the
    # table (ref aa per position) and the MANE protein is mapped onto it (identity, or trusted blocks
    # of a global alignment).  If the gene's accession is absent from AlphaMissense, another reviewed
    # entry of the same gene is used (GNAS: P63092 absent -> Q5JWF2, flagged in am_accession).
    from s06_select_genes import align_map
    g_acc = dict(zip(g.hgnc_id, g.uniprot))
    v["uniprot"] = v.hgnc_id.map(g_acc)
    am_seq = {}
    for acc, sub in am.drop_duplicates(["uniprot_id", "pos"]).groupby("uniprot_id"):
        sub = sub.sort_values("pos")
        s_ = ["X"] * int(sub.pos.max())
        for p_, r_ in zip(sub.pos, sub.ref):
            s_[p_ - 1] = r_
        am_seq[acc] = "".join(s_)
    am_score = dict(zip(zip(am.uniprot_id, am.protein_variant), zip(am.am_pathogenicity, am.am_class)))
    mane_seq = dict(zip(g.hgnc_id, g.vep_protein_seq))
    other = dict(zip(g.hgnc_id, g.other_hits))
    gene_am = {}
    for h in v[v.variant_class == "missense"].hgnc_id.dropna().unique():
        acc = g_acc.get(h); ms_ = mane_seq.get(h)
        cand = [acc] if isinstance(acc, str) and acc in am_seq else \
            [x for x in (other.get(h) or "").split(",") if x in am_seq]
        if not cand or not isinstance(ms_, str):
            gene_am[h] = (None, None, None); continue
        best = None
        for c in cand:
            if am_seq[c] == ms_:
                best = (c, None, 1.0); break
            m_, idn = align_map(ms_, am_seq[c])
            if best is None or idn > best[2]:
                best = (c, m_, idn)
        gene_am[h] = best
    score, klass, status, am_acc_col, am_pos_col, am_pv = [], [], [], [], [], []
    for h, c, p, r, a in zip(v.hgnc_id, v.variant_class, v.protein_pos, v.vep_p_ref, v.vep_p_alt):
        if c != "missense":
            for L_ in (score, klass, status, am_acc_col, am_pos_col, am_pv):
                L_.append(None)
            continue
        acc, m_, idn = gene_am.get(h, (None, None, None))
        st, sc, kl, ap, pv = None, None, None, None, None
        if acc is None:
            st = "accession_not_in_alphamissense"
        elif pd.isna(p) or not isinstance(r, str) or not isinstance(a, str):
            st = "no_protein_change_parsed"
        else:
            ap = int(p) if m_ is None else m_.get(int(p))
            if ap is None:
                st = "position_not_alignable"
            elif ap > len(am_seq[acc]) or am_seq[acc][ap - 1] != r:
                st = "ref_aa_mismatch"
            else:
                pv = f"{r}{ap}{a}"
                hit = am_score.get((acc, pv))
                if hit is None:
                    st = "alt_not_found"
                else:
                    sc, kl = hit; st = "matched" if m_ is None else "matched_via_alignment"
        score.append(sc); klass.append(kl); status.append(st); am_acc_col.append(acc); am_pos_col.append(ap); am_pv.append(pv)
    v["am_accession"] = am_acc_col; v["am_pos"] = pd.array(am_pos_col, dtype="Int64"); v["am_protein_variant"] = am_pv
    v["am_pathogenicity"] = pd.to_numeric(pd.Series(score, index=v.index), errors="coerce")
    v["am_class"] = klass; v["am_match_status"] = status
    v["am_accession_is_gene_primary"] = [(x == y) if isinstance(x, str) else None for x, y in zip(v.am_accession, v.uniprot)]

    # ---- gene table ----
    sm = smeta.set_index("uniprot")
    mp = []
    for p, ok in zip(g.af_pdb, g.af_pdb_exists.fillna(False)):
        mp.append(mean_plddt(p) if ok else (np.nan, np.nan, np.nan))
    g["af_mean_plddt"] = [x[0] for x in mp]; g["af_frac_plddt_lt50"] = [x[1] for x in mp]; g["af_frac_plddt_ge70"] = [x[2] for x in mp]
    for c in ["q_core", "q_semi", "n_core", "n_semi_only", "n_disorder", "G_all_edges", "G_core_nodes", "G_core_edges",
              "anchor_rule", "n_N_anchor_nodes", "n_C_anchor_nodes", "backbone_pairs", "backbone_pairs_connected",
              "backbone_size", "backbone_strict_size"]:
        g["struct_" + c] = g.uniprot.map(sm[c]) if c in sm else None
    g["residue_features_file"] = [str(RFDIR / f"{a}.parquet") if s else None for a, s in zip(g.uniprot, g.selected)]
    g["n_individuals_missense_ge20"] = g.n_individuals_missense >= MIN_INDIVIDUALS
    mc = v[v.variant_class == "missense"].groupby("hgnc_id").am_match_status.agg(lambda s: s.str.startswith("matched").mean())
    g["am_missense_match_rate"] = g.hgnc_id.map(mc)
    genes = g.rename(columns={"af_pdb": "structure_file", "sequence": "uniprot_sequence"}).drop(columns=["vep_protein_seq", "af_sequence"])
    genes["af_fragmented"] = genes.af_n_fragments > 1

    # ---- residue features merged at protein_pos ----
    rfs = []
    for acc in genes[genes.selected].uniprot:
        rf = pd.read_parquet(RFDIR / f"{acc}.parquet")
        rfs.append(rf)
    rf = pd.concat(rfs, ignore_index=True)
    rfm = rf.rename(columns={"pos": "protein_pos_int"}).drop(columns=["gene"])
    rfm = rfm.rename(columns={c: ("res_" + c) for c in rfm.columns if c not in ("uniprot", "protein_pos_int")})
    sel_h = set(genes[genes.selected].hgnc_id)
    v["gene_selected"] = v.hgnc_id.isin(sel_h)
    v["protein_pos_int"] = v.structure_pos.astype("Int64")
    v["protein_length"] = v.mane_protein_seq_len          # MANE protein length (VEP coordinates)
    v["uniprot_length"] = v.hgnc_id.map(dict(zip(genes.hgnc_id, genes.uniprot_length)))
    v["rel_protein_pos"] = v.protein_pos / v.protein_length  # relative position on the MANE protein
    v = v.merge(rfm, on=["uniprot", "protein_pos_int"], how="left")
    feat_cols = [c for c in v.columns if c.startswith("res_")]
    # residue features only for mappable variants in selected genes with a matching reference aa
    v["af_length"] = v.uniprot.map(sm["L"])
    usable = v.gene_selected & v.mappable & v.structure_pos.notna() & (v.structure_pos <= v.af_length)
    v.loc[~usable, feat_cols] = np.nan
    v["has_structure_features"] = (usable & v.res_plddt.notna()).astype(bool)
    v = v.drop(columns=["protein_pos_int"])
    # ---- individuals ----
    obs = ph.groupby("phenopacket_id").value.agg(n_hpo_explicit_observed="sum", n_hpo_explicit="count")
    ind = ind.merge(obs, left_on="phenopacket_id", right_index=True, how="left")
    ind["n_hpo_explicit"] = ind.n_hpo_explicit.fillna(0).astype(int)
    ind["n_hpo_explicit_observed"] = ind.n_hpo_explicit_observed.fillna(0).astype(int)
    ind["n_hpo_explicit_excluded"] = ind.n_hpo_explicit - ind.n_hpo_explicit_observed
    vc = v.groupby("phenopacket_id").agg(variant_classes=("variant_class", lambda s: "+".join(sorted(s))),
                                         n_mappable_variants=("mappable", "sum"),
                                         any_missense=("variant_class", lambda s: (s == "missense").any()),
                                         any_ptc=("variant_class", lambda s: s.isin(["nonsense", "frameshift"]).any()))
    ind = ind.merge(vc, left_on="phenopacket_id", right_index=True, how="left")
    ind["hgnc_id"] = ind.phenopacket_id.map(v.groupby("phenopacket_id").hgnc_id.first())
    ind["uniprot"] = ind.hgnc_id.map(g_acc)
    ind["gene_selected"] = ind.hgnc_id.isin(sel_h)
    ind["has_mappable_variant"] = ind.n_mappable_variants.fillna(0) > 0
    # ---- write ----
    ind.to_parquet(DERIVED / "individuals.parquet", index=False)
    v.to_parquet(DERIVED / "variants.parquet", index=False)
    ph.to_parquet(DERIVED / "phenotypes_explicit.parquet", index=False)
    pp.to_parquet(DERIVED / "phenotypes_propagated.parquet", index=False)
    genes.to_parquet(DERIVED / "genes.parquet", index=False)
    # also csv.gz copies of the small tables for convenience
    for name, df in (("individuals", ind), ("genes", genes.drop(columns=["uniprot_sequence"])), ("variants", v)):
        df.to_csv(DERIVED / f"{name}.csv.gz", index=False)

    # ---- counts ----
    vs = v[v.gene_selected]
    C["phenopackets"] = len(ind)
    C["genomic_interpretations"] = len(v)
    C["variant_class_all"] = v.variant_class.value_counts().pipe(vcd)
    C["variant_class_selected_genes"] = vs.variant_class.value_counts().pipe(vcd)
    C["genes_total"] = len(genes)
    C["genes_selected"] = int(genes.selected.sum())
    C["selected_gene_list"] = genes[genes.selected].sort_values("n_individuals_mappable", ascending=False).gene.tolist()
    C["individuals_in_selected_genes"] = int(ind.gene_selected.sum())
    C["individuals_in_selected_genes_with_mappable_variant"] = int((ind.gene_selected & ind.has_mappable_variant).sum())
    C["individuals_with_structure_features"] = int(v[v.has_structure_features].phenopacket_id.nunique())
    C["genes_missense_ge20_all"] = int((genes.n_individuals_missense >= MIN_INDIVIDUALS).sum())
    C["genes_missense_ge20_selected"] = int((genes.selected & (genes.n_individuals_missense >= MIN_INDIVIDUALS)).sum())
    C["genes_missense_ge20_selected_list"] = genes[genes.selected & (genes.n_individuals_missense >= MIN_INDIVIDUALS)].gene.tolist()
    ms = v[v.variant_class == "missense"]
    C["am_missense_rows_all"] = len(ms); C["am_matched_all"] = int(ms.am_match_status.str.startswith("matched").sum())
    C["am_status_all"] = ms.am_match_status.value_counts().pipe(vcd)
    mss = ms[ms.gene_selected]
    C["am_missense_rows_selected"] = len(mss); C["am_matched_selected"] = int(mss.am_match_status.str.startswith("matched").sum())
    C["am_status_selected"] = mss.am_match_status.value_counts().pipe(vcd)
    C["p_pos_compared"] = int(v.p_pos_agree.notna().sum()); C["p_pos_agree"] = int(v.p_pos_agree.sum())
    C["p_pos_agree_by_class"] = v[v.p_pos_agree.notna()].groupby("variant_class").p_pos_agree.mean().round(4).to_dict()
    C["hgvs_p_compared"] = int(v.hgvs_p_string_agree.notna().sum()); C["hgvs_p_agree"] = int(v.hgvs_p_string_agree.sum())
    C["hgvs_c_compared"] = int(v.hgvs_c_string_agree.notna().sum()); C["hgvs_c_agree"] = int(v.hgvs_c_string_agree.sum())
    C["ppkt_transcript_is_mane"] = int(v.ppkt_tx_is_mane.sum()); C["ppkt_transcript_present"] = int(v.transcript_c.notna().sum())
    C["ref_aa_matches_af_selected"] = vs.ref_aa_matches_af.value_counts(dropna=False).pipe(vcd)
    ptc = v[v.variant_class.isin(["nonsense", "frameshift"])]
    C["ptc_all"] = len(ptc); C["nmd_escape_all"] = ptc.nmd_escape.value_counts(dropna=False).pipe(vcd)
    pa = ptc[ptc.gene == "ANKRD11"]
    C["ankrd11_ptc"] = len(pa); C["ankrd11_nmd_escape"] = pa.nmd_escape.value_counts(dropna=False).pipe(vcd)
    C["ankrd11_individuals"] = int((ind.gene == "ANKRD11").sum())
    C["ankrd11_variant_class"] = v[v.gene == "ANKRD11"].variant_class.value_counts().pipe(vcd)
    C["phenotypes_explicit_cells"] = len(ph); C["phenotypes_propagated_cells"] = len(pp)
    C["selection_fail_reasons"] = {
        "n20_mappable": int(genes.passes_n20_mappable.sum()),
        "n20_mappable_but_af_fragmented": genes[genes.passes_n20_mappable & (genes.af_n_fragments > 1)].gene.tolist(),
        "n20_mappable_but_no_af_model": genes[genes.passes_n20_mappable & (genes.af_n_fragments == 0)].gene.tolist(),
        "n20_mappable_af_single_but_lt20_on_structure": genes[genes.passes_n20_mappable & (genes.af_n_fragments == 1) & ~genes.passes_n20_on_structure].gene.tolist(),
    }
    C["fragmented_genes_all"] = genes[genes.af_n_fragments > 1].gene.tolist()
    C["genes_without_uniprot"] = genes[genes.uniprot.isna()].gene.tolist()
    (DERIVED / "counts.json").write_text(json.dumps(C, indent=1, default=str))
    print(json.dumps(C, indent=1, default=str))


if __name__ == "__main__":
    main()
