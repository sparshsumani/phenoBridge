#!/usr/bin/env python3
"""Step 6: merge VEP annotation onto every (individual x variant), fetch MANE protein sequences
(Ensembl REST, cached) and decide the gene set.

Mappable variant classes (protein position defined): missense, nonsense, frameshift,
inframe_indel, start_lost, stop_lost.
Selection: gene has >= 20 individuals carrying >= 1 mappable variant in that gene AND a
single-fragment AlphaFold model (proteome v6) for its UniProt accession AND the AlphaFold /
UniProt sequence equals the MANE Select protein used by VEP (so VEP positions index the model).
Outputs INTER/variants_annot.parquet, INTER/gene_summary.parquet, INTER/selection_log.json
"""
import json
from collections import Counter
import numpy as np
import pandas as pd
from common import INTER, CACHE, http_json, MIN_INDIVIDUALS
from s03_vep_annotate import parse_p

MAPPABLE = {"missense", "nonsense", "frameshift", "inframe_indel", "start_lost", "stop_lost"}
SEQDIR = CACHE / "ens_sequence"; SEQDIR.mkdir(parents=True, exist_ok=True)


def ensp_sequences(ids):
    ids = sorted(set(i for i in ids if isinstance(i, str)))
    todo = [i for i in ids if not (SEQDIR / f"{i}.json").exists()]
    for k in range(0, len(todo), 50):
        b = todo[k:k + 50]
        res = http_json("https://rest.ensembl.org/sequence/id?type=protein", "POST", {"ids": b})
        got = {r["query"]: r for r in res} if isinstance(res, list) else {}
        for i in b:
            (SEQDIR / f"{i}.json").write_text(json.dumps(got.get(i)))
    out = {}
    for i in ids:
        r = json.loads((SEQDIR / f"{i}.json").read_text())
        out[i] = r["seq"].rstrip("*") if r else None
    return out


def af_sequence(path):
    import gzip
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    seq, last = [], None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    seq.append(index_to_one(three_to_index(line[17:20])))
                except Exception:
                    seq.append("X")
    return "".join(seq)


def main():
    var = pd.read_parquet(INTER / "variants_raw.parquet")
    ann = pd.read_parquet(INTER / "vep_annotations.parquet")
    gu = pd.read_parquet(INTER / "gene_uniprot.parquet")
    af = pd.read_parquet(INTER / "af_models.parquet")
    log = Counter()
    v = var.merge(ann.drop(columns=["vep_input"]), on=["chrom", "pos", "ref", "alt", "hgnc_id"], how="left")
    assert len(v) == len(var)
    v.loc[v.is_structural, "variant_class"] = "structural"
    v.loc[v.is_structural, "vep_status"] = "not_run_structural"
    v["variant_class"] = v.variant_class.fillna("unannotated")
    # phenopacket hgvs.p parsing & agreement with VEP
    pp = [parse_p(p) for p in v.hgvs_p]
    v["ppkt_p_pos"] = pd.array([p["pos"] if p else None for p in pp], dtype="Int64")
    v["ppkt_p_ref"] = [p["ref"] if p else None for p in pp]
    v["ppkt_p_alt"] = [p["alt"] if p else None for p in pp]
    v["ppkt_p_norm"] = [h.split(":", 1)[1].replace("(", "").replace(")", "") if isinstance(h, str) else None for h in v.hgvs_p]
    v["vep_p_norm"] = [h.split(":", 1)[1] if isinstance(h, str) else None for h in v.vep_hgvsp]
    v["ppkt_c_norm"] = [h.split(":", 1)[1] if isinstance(h, str) else None for h in v.hgvs_c]
    v["vep_c_norm"] = [h.split(":", 1)[1] if isinstance(h, str) else None for h in v.vep_hgvsc]
    v["ppkt_tx_is_mane"] = [(isinstance(a, str) and isinstance(b, str) and a.split(".")[0] == b.split(".")[0])
                            for a, b in zip(v.transcript_c, v.mane_select_refseq)]
    v["ppkt_tx_is_mane_exact_version"] = [(isinstance(a, str) and a == b) for a, b in zip(v.transcript_c, v.mane_select_refseq)]
    both = v.ppkt_p_pos.notna() & v.vep_p_pos.notna()
    v["p_pos_agree"] = pd.array([bool(a == b) if m else None for a, b, m in zip(v.ppkt_p_pos, v.vep_p_pos, both)], dtype="boolean")
    v["hgvs_p_string_agree"] = pd.array([(a == b) if (isinstance(a, str) and isinstance(b, str)) else None
                                         for a, b in zip(v.ppkt_p_norm, v.vep_p_norm)], dtype="boolean")
    v["hgvs_c_string_agree"] = pd.array([(a == b) if (isinstance(a, str) and isinstance(b, str)) else None
                                         for a, b in zip(v.ppkt_c_norm, v.vep_c_norm)], dtype="boolean")
    # protein position used downstream: VEP (MANE) hgvsp first residue, fallback VEP protein_start
    v["protein_pos"] = v.vep_p_pos.fillna(v.vep_protein_start)
    v["mappable"] = v.variant_class.isin(MAPPABLE) & v.protein_pos.notna()
    # sequences
    seqs = ensp_sequences(v.ensembl_protein.dropna().unique())
    v["mane_protein_seq_len"] = [len(seqs[p]) if isinstance(p, str) and seqs.get(p) else np.nan for p in v.ensembl_protein]
    # gene-level table
    g = gu.merge(af, on=["gene", "uniprot"], how="left")
    mane_ensp = v.dropna(subset=["ensembl_protein"]).query("transcript_source=='MANE_Select'").groupby("hgnc_id").ensembl_protein.agg(lambda s: s.mode().iat[0])
    any_ensp = v.dropna(subset=["ensembl_protein"]).groupby("hgnc_id").ensembl_protein.agg(lambda s: s.mode().iat[0])
    g["vep_protein"] = g.hgnc_id.map(mane_ensp).fillna(g.hgnc_id.map(any_ensp))
    g["vep_transcript"] = g.hgnc_id.map(v.dropna(subset=["ensembl_transcript"]).groupby("hgnc_id").ensembl_transcript.agg(lambda s: s.mode().iat[0]))
    g["vep_protein_seq"] = g.vep_protein.map(lambda p: seqs.get(p) if isinstance(p, str) else None)
    g["af_sequence"] = [af_sequence(p) if isinstance(p, str) and ok else None for p, ok in zip(g.af_pdb, g.af_pdb_exists.fillna(False))]
    g["af_eq_uniprot_seq"] = [(a == b) if isinstance(a, str) and isinstance(b, str) else None for a, b in zip(g.af_sequence, g.sequence)]
    g["mane_eq_uniprot_seq"] = [(a == b) if isinstance(a, str) and isinstance(b, str) else None for a, b in zip(g.vep_protein_seq, g.sequence)]
    g["mane_eq_af_seq"] = [(a == b) if isinstance(a, str) and isinstance(b, str) else None for a, b in zip(g.vep_protein_seq, g.af_sequence)]
    # per-variant ref aa check against AF model sequence
    afseq = dict(zip(g.hgnc_id, g.af_sequence))
    def refchk(h, pos, ref, cls):
        s = afseq.get(h)
        if not isinstance(s, str) or pd.isna(pos) or not isinstance(ref, str) or ref in ("?", "*"):
            return None
        pos = int(pos)
        return bool(1 <= pos <= len(s) and s[pos - 1] == ref)
    v["ref_aa_matches_af"] = pd.array([refchk(h, p, r, c) for h, p, r, c in zip(v.hgnc_id, v.protein_pos, v.vep_p_ref, v.variant_class)], dtype="boolean")
    # counts per gene
    ind = pd.read_parquet(INTER / "individuals_raw.parquet")
    cnt = v.groupby("hgnc_id").agg(n_individuals=("phenopacket_id", "nunique"))
    cnt["n_individuals_mappable"] = v[v.mappable].groupby("hgnc_id").phenopacket_id.nunique()
    cnt["n_individuals_missense"] = v[v.variant_class == "missense"].groupby("hgnc_id").phenopacket_id.nunique()
    cnt["n_individuals_ptc"] = v[v.variant_class.isin(["nonsense", "frameshift"])].groupby("hgnc_id").phenopacket_id.nunique()
    cnt["n_individuals_structural"] = v[v.variant_class == "structural"].groupby("hgnc_id").phenopacket_id.nunique()
    cnt["n_unique_variants"] = v.groupby("hgnc_id").variant_id.nunique()
    cnt["n_unique_mappable_positions"] = v[v.mappable].groupby("hgnc_id").protein_pos.nunique()
    cnt = cnt.fillna(0).astype(int)
    g = g.merge(cnt, left_on="hgnc_id", right_index=True, how="left")
    g["cohort"] = g.hgnc_id.map(v.groupby("hgnc_id").cohort.first())
    g["passes_n20_mappable"] = g.n_individuals_mappable >= MIN_INDIVIDUALS
    g["passes_af_single"] = g.af_n_fragments == 1
    g["passes_seq_concordance"] = g.mane_eq_af_seq.fillna(False).astype(bool)
    g["selected"] = g.passes_n20_mappable & g.passes_af_single & g.passes_seq_concordance
    g["af_mean_plddt"] = np.nan
    v.to_parquet(INTER / "variants_annot.parquet", index=False)
    g.to_parquet(INTER / "gene_summary.parquet", index=False)
    log.update(dict(
        genes_total=len(g), genes_uniprot=int(g.uniprot.notna().sum()),
        genes_n20_mappable=int(g.passes_n20_mappable.sum()),
        genes_n20_mappable_and_af_single=int((g.passes_n20_mappable & g.passes_af_single).sum()),
        genes_n20_mappable_fragmented=g[g.passes_n20_mappable & (g.af_n_fragments > 1)].gene.tolist(),
        genes_n20_mappable_no_af=g[g.passes_n20_mappable & (g.af_n_fragments == 0)].gene.tolist(),
        genes_n20_mappable_af_single_seq_mismatch=g[g.passes_n20_mappable & g.passes_af_single & ~g.passes_seq_concordance].gene.tolist(),
        genes_selected=int(g.selected.sum()),
        genes_missense_n20=int((g.n_individuals_missense >= MIN_INDIVIDUALS).sum()),
        genes_missense_n20_selected=int((g.selected & (g.n_individuals_missense >= MIN_INDIVIDUALS)).sum()),
        p_pos_compared=int(v.p_pos_agree.notna().sum()), p_pos_agree=int(v.p_pos_agree.sum()),
        hgvs_p_compared=int(v.hgvs_p_string_agree.notna().sum()), hgvs_p_agree=int(v.hgvs_p_string_agree.sum()),
        hgvs_c_compared=int(v.hgvs_c_string_agree.notna().sum()), hgvs_c_agree=int(v.hgvs_c_string_agree.sum()),
        ppkt_tx_is_mane=int(v.ppkt_tx_is_mane.sum()), rows=len(v),
    ))
    (INTER / "selection_log.json").write_text(json.dumps(dict(log), indent=1))
    print(json.dumps(dict(log), indent=1))
    print(g[g.selected][["gene", "uniprot", "uniprot_length", "n_individuals", "n_individuals_mappable", "n_individuals_missense"]].sort_values("n_individuals_mappable", ascending=False).to_string())


if __name__ == "__main__":
    main()
