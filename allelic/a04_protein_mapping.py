#!/usr/bin/env python3
"""Step 4: map every ClinVar missense variant to the MANE Select protein, the reviewed UniProt
canonical isoform and the AlphaFold v6 model numbering; check the reference amino acid.

MANE protein change
  * gene = symbol in the ClinVar Name '(GENE)'; MANE Select row from MANE v1.5 summary by symbol
    (fallback: ClinVar HGNC_ID).  MANE protein sequence from MANE.GRCh38.v1.5.refseq_protein.faa.gz.
  * method 'clinvar_tx_is_mane': ClinVar Name transcript (version ignored) == MANE RefSeq_nuc and the
    ClinVar ref aa equals the MANE protein residue at that position -> MANE position = ClinVar position.
  * otherwise method 'vep_mane': Ensembl VEP REST GRCh38 (multigene/s03 run_vep, POST batches of 200,
    cached in cache/vep_region/ keyed by input 'chrom pos . ref alt . . .' from PositionVCF/Ref/AltVCF);
    the MANE Select transcript consequence of the same gene must contain missense_variant; MANE
    position/ref/alt from protein_start + amino_acids.  Variants that are not missense on MANE are dropped.
  * gene table also covers genes of phenopacket-store missense variants (by HGNC id -> MANE row).
  * validation: a random 500 'clinvar_tx_is_mane' variants are also run through VEP (agreement logged).
UniProt / AlphaFold
  * UniProt reviewed human entry (REST stream, raw/uniprot/) whose MANE-Select xref equals the MANE
    Ensembl transcript (version ignored); fallback HGNC xref; fallback primary gene name.
  * AlphaFold: fragments per accession from the v6 proteome tar listing (multigene/s05); single-fragment
    models are extracted (tar -T list) into raw/alphafold/models/.  Fragmented models (proteins >2,700 aa)
    are excluded and listed.
  * MANE pos -> AF pos and -> UniProt pos by identity when sequences are equal, else multigene/s06
    align_map (global alignment; only aligned blocks >=10 aa with >=90% identity are trusted).
  * ref_matches_af: MANE ref aa == AF model residue at the mapped position.
Outputs intermediate/clinvar_mapped.parquet, intermediate/gene_protein_map.parquet,
        intermediate/protein_mapping_log.json
"""
import re, sys, json, subprocess, gzip
from collections import Counter, defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from common import (FILES, INTER, DERIVED_PPKT, AF_TAR, AF_DIR, AF_LISTING, AA3, MAX_AF_SINGLE_FRAGMENT_LEN, save_json,
                    read_fasta_gz, mg_import, af_model_path)

s03 = mg_import("s03_vep_annotate")
s06 = mg_import("s06_select_genes")


def load_mane():
    m = pd.read_csv(FILES["mane_summary"], sep="\t", dtype=str)
    m = m[m.MANE_status == "MANE Select"].copy()
    prot = read_fasta_gz(FILES["mane_protein"])
    m["mane_protein_seq"] = m.RefSeq_prot.map(prot)
    m["nm_nov"] = m.RefSeq_nuc.str.split(".").str[0]
    m["enst_nov"] = m.Ensembl_nuc.str.split(".").str[0]
    return m


def load_uniprot():
    u = pd.read_csv(FILES["uniprot_human"], sep="\t", dtype=str)
    u.columns = ["accession", "entry_name", "gene_primary", "gene_synonyms", "hgnc", "mane_select", "length", "sequence", "seq_version"]
    return u


def vep_pick(res, gene, hgnc):
    if not isinstance(res, dict) or "__error__" in res:
        return None, "vep_error"
    tcs = [t for t in res.get("transcript_consequences", []) if t.get("gene_symbol") == gene or (hgnc and t.get("hgnc_id") == hgnc)]
    pick = next((t for t in tcs if t.get("mane_select")), None)
    if pick is None:
        return None, "vep_no_mane_transcript_for_gene"
    if "missense_variant" not in pick.get("consequence_terms", []):
        return pick, "vep_not_missense_on_mane:" + ",".join(pick.get("consequence_terms", []))
    aa = (pick.get("amino_acids") or "").split("/")
    if len(aa) != 2 or len(aa[0]) != 1 or len(aa[1]) != 1:
        return pick, "vep_amino_acids_unparsed"
    return pick, "ok"


def main():
    log = Counter()
    cv = pd.read_parquet(INTER / "clinvar_missense.parquet")
    mane = load_mane()
    by_sym = mane.drop_duplicates("symbol").set_index("symbol", drop=False)
    by_hgnc = mane.drop_duplicates("HGNC_ID").set_index("HGNC_ID", drop=False)
    rows = []
    for r in cv.itertuples():
        mr = by_sym.loc[r.name_gene] if r.name_gene in by_sym.index else None
        how = "symbol"
        if mr is None:
            for h in str(r.clinvar_hgnc_ids).split(";"):
                if h in by_hgnc.index:
                    mr = by_hgnc.loc[h]; how = "hgnc"; break
        d = dict(variation_id=r.variation_id, gene=r.name_gene, mane_lookup=how if mr is not None else None)
        if mr is None:
            d["map_status"] = "no_mane_select_for_gene"; rows.append(d); continue
        d.update(gene_mane=mr.symbol, hgnc_id=mr.HGNC_ID, mane_nm=mr.RefSeq_nuc,
                 mane_np=mr.RefSeq_prot, mane_enst=mr.Ensembl_nuc, mane_ensp=mr.Ensembl_prot)
        seq = mr.mane_protein_seq
        d["clinvar_tx_is_mane"] = r.name_tx_noversion == mr.nm_nov
        d["clinvar_tx_is_mane_exact_version"] = r.name_tx == mr.RefSeq_nuc
        p = r.clinvar_protein_pos
        if d["clinvar_tx_is_mane"] and isinstance(seq, str) and 1 <= p <= len(seq) and seq[p - 1] == r.ref_aa:
            d.update(map_status="ok", map_method="clinvar_tx_is_mane", mane_pos=p, mane_ref=r.ref_aa, mane_alt=r.alt_aa)
        else:
            d["map_status"] = "needs_vep"
            d["needs_vep_reason"] = "non_mane_transcript" if not d["clinvar_tx_is_mane"] else "ref_mismatch_on_mane_seq"
        rows.append(d)
    mp = pd.DataFrame(rows)
    mp = cv[["variation_id", "chrom", "pos_vcf", "ref_vcf", "alt_vcf", "name_tx", "clinvar_protein_pos", "ref_aa", "alt_aa",
             "hgvs_p", "hgvs_c"]].merge(mp, on="variation_id")
    log.update({f"pre_vep_{k}": int(v) for k, v in mp.map_status.value_counts().items()})
    log.update({f"needs_vep_{k}": int(v) for k, v in mp.needs_vep_reason.value_counts().items()})

    # ---- VEP for non-MANE + random validation sample ----
    need = mp.map_status == "needs_vep"
    rng = np.random.RandomState(0)
    ok_idx = mp.index[mp.map_status == "ok"]
    val_idx = rng.choice(ok_idx, size=min(500, len(ok_idx)), replace=False)
    vcf_ok = mp.pos_vcf.str.fullmatch(r"\d+") & mp.ref_vcf.str.fullmatch("[ACGT]") & mp.alt_vcf.str.fullmatch("[ACGT]")
    mp["vep_input"] = [f"{c} {p} . {a} {b} . . ." if ok else None for c, p, a, b, ok in
                       zip(mp.chrom, mp.pos_vcf, mp.ref_vcf, mp.alt_vcf, vcf_ok)]
    run_idx = sorted(set(mp.index[need & vcf_ok]) | set(val_idx))
    log["vep_inputs"] = len(run_idx)
    log["needs_vep_but_no_vcf_coords"] = int((need & ~vcf_ok).sum())
    vep = s03.run_vep(sorted(set(mp.loc[run_idx, "vep_input"])))
    val = Counter()
    val_set = set(val_idx)
    for i in run_idx:
        r = mp.loc[i]
        res = vep.get(r.vep_input)
        pick, st = vep_pick(res, r.gene, r.hgnc_id)
        if pick is not None:
            mp.loc[i, "vep_hgvsp"] = pick.get("hgvsp"); mp.loc[i, "vep_hgvsc"] = pick.get("hgvsc")
            mp.loc[i, "vep_mane_transcript"] = pick.get("transcript_id")
        vpos = vref = valt = None
        if st == "ok":
            vref, valt = pick["amino_acids"].split("/"); vpos = int(pick["protein_start"])
        if i in val_set and r.map_status == "ok":
            val["validated"] += 1
            val["agree" if (st == "ok" and vpos == r.mane_pos and vref == r.mane_ref and valt == r.mane_alt) else "disagree:" + st] += 1
            continue
        if st == "ok":
            mp.loc[i, ["map_status", "map_method", "mane_pos", "mane_ref", "mane_alt"]] = ["ok", "vep_mane", vpos, vref, valt]
        else:
            mp.loc[i, "map_status"] = st.split(":")[0]
            mp.loc[i, "vep_status_detail"] = st
    mp.loc[need & ~vcf_ok, "map_status"] = "no_vcf_coordinates"
    log["vep_validation_on_mane_transcript_variants"] = dict(val)
    mp["mane_pos"] = pd.to_numeric(mp.mane_pos, errors="coerce").astype("Int64")
    # check VEP-mapped ref against MANE sequence too
    seqs = dict(zip(mane.RefSeq_prot, mane.mane_protein_seq))
    mp["ref_matches_mane_seq"] = pd.array([
        (bool(1 <= int(p) <= len(seqs.get(np_, "") or "") and seqs[np_][int(p) - 1] == ref) if pd.notna(p) and isinstance(np_, str) else None)
        for p, np_, ref in zip(mp.mane_pos, mp.mane_np, mp.mane_ref)], dtype="boolean")
    log.update({f"post_vep_{k}": int(v) for k, v in mp.map_status.value_counts().items()})
    log["map_method"] = {k: int(v) for k, v in mp.map_method.value_counts().items()}
    log["vep_mapped_ref_matches_mane_seq"] = int(mp[mp.map_method == "vep_mane"].ref_matches_mane_seq.sum())
    log["vep_mapped_total"] = int((mp.map_method == "vep_mane").sum())

    # ---- genes -> UniProt -> AlphaFold ----
    up = load_uniprot()
    up["mane_enst_nov"] = up.mane_select.fillna("").map(lambda s: [x.split(".")[0] for x in s.split(";") if x.strip()])
    up["hgnc_list"] = up.hgnc.fillna("").map(lambda s: [x for x in s.split(";") if x.strip()])
    by_enst, by_h, by_g = defaultdict(list), defaultdict(list), defaultdict(list)
    for k, r in up.iterrows():
        for e in r.mane_enst_nov:
            by_enst[e].append(k)
        for h in r.hgnc_list:
            by_h[h].append(k)
        if isinstance(r.gene_primary, str):
            by_g[r.gene_primary].append(k)
    names = AF_LISTING.read_text().split() if AF_LISTING.exists() else subprocess.run(
        ["tar", "-tf", str(AF_TAR)], capture_output=True, text=True, check=True).stdout.split()
    frags = defaultdict(set)
    for n in names:
        m = re.match(r"AF-([A-Z0-9]+)-F(\d+)-model_v(\d+)\.pdb\.gz$", n)
        if m:
            frags[m.group(1)].add(int(m.group(2)))
    genes = mp[mp.map_status == "ok"].drop_duplicates("gene")[["gene", "hgnc_id", "mane_enst", "mane_np", "mane_ensp"]]
    # + genes of phenopacket-store missense variants (multigene derived table), via their HGNC id
    pv = pd.read_parquet(DERIVED_PPKT / "variants.parquet", columns=["hgnc_id", "variant_class"])
    ph = set(pv[pv.variant_class == "missense"].hgnc_id.dropna()) - set(genes.hgnc_id)
    extra = by_hgnc.loc[sorted(h for h in ph if h in by_hgnc.index)].reset_index(drop=True)
    extra = pd.DataFrame(dict(gene=extra.symbol, hgnc_id=extra.HGNC_ID, mane_enst=extra.Ensembl_nuc,
                              mane_np=extra.RefSeq_prot, mane_ensp=extra.Ensembl_prot))
    log["ppkt_only_genes_added"] = len(extra)
    genes = pd.concat([genes, extra], ignore_index=True)
    grows = []
    for g in genes.itertuples():
        e = g.mane_enst.split(".")[0]
        cand, how = by_enst.get(e, []), "uniprot_mane_xref"
        if not cand:
            cand, how = by_h.get(g.hgnc_id, []), "hgnc_xref"
        if not cand:
            cand, how = by_g.get(g.gene, []), "gene_name"
        if not cand:
            grows.append(dict(gene=g.gene, uniprot=None, uniprot_map=None)); continue
        k = cand[0]
        if len(cand) > 1:  # prefer the one whose sequence equals MANE protein
            ms = seqs.get(g.mane_np)
            k = next((c for c in cand if up.loc[c, "sequence"] == ms), cand[0])
        u = up.loc[k]
        nf = len(frags.get(u.accession, ()))
        grows.append(dict(gene=g.gene, hgnc_id=g.hgnc_id, mane_enst=g.mane_enst, mane_np=g.mane_np, mane_ensp=g.mane_ensp,
                          uniprot=u.accession, uniprot_map=how, n_uniprot_candidates=len(cand), uniprot_entry=u.entry_name,
                          uniprot_length=int(u.length), uniprot_seq=u.sequence, uniprot_seq_version=u.seq_version,
                          mane_protein_seq=seqs.get(g.mane_np), af_n_fragments=nf))
    gm = pd.DataFrame(grows)
    gm["af_single_fragment"] = gm.af_n_fragments == 1
    todo = sorted({f"AF-{a}-F1-model_v6.pdb.gz" for a in gm[gm.af_single_fragment].uniprot if not af_model_path(a).exists()})
    log["af_models_to_extract"] = len(todo)
    if todo:
        lst = INTER / "af_extract_list.txt"
        lst.write_text("\n".join(todo) + "\n")
        print("extracting", len(todo), "AF models", flush=True)
        subprocess.run(["tar", "-xf", str(AF_TAR), "-C", str(AF_DIR), "-T", str(lst)], check=True)
    gm["af_pdb"] = [str(af_model_path(a)) if s and af_model_path(a).exists() else None for a, s in zip(gm.uniprot, gm.af_single_fragment)]
    gm["af_seq"] = [s06.af_sequence(p) if isinstance(p, str) else None for p in gm.af_pdb]
    eq = lambda a, b: (a == b) if isinstance(a, str) and isinstance(b, str) else None
    gm["af_eq_uniprot"] = [eq(a, b) for a, b in zip(gm.af_seq, gm.uniprot_seq)]
    gm["mane_eq_uniprot"] = [eq(a, b) for a, b in zip(gm.mane_protein_seq, gm.uniprot_seq)]
    gm["mane_eq_af"] = [eq(a, b) for a, b in zip(gm.mane_protein_seq, gm.af_seq)]
    maps_af, maps_up = {}, {}
    for r in gm.itertuples():
        for tgt, store, col in ((r.af_seq, maps_af, "af"), (r.uniprot_seq, maps_up, "up")):
            if isinstance(r.mane_protein_seq, str) and isinstance(tgt, str):
                store[r.gene] = (None, 1.0) if r.mane_protein_seq == tgt else s06.align_map(r.mane_protein_seq, tgt)
    gm["mane_af_map"] = gm.gene.map(lambda g: None if g not in maps_af else ("identity" if maps_af[g][0] is None else "alignment"))
    gm["mane_af_identity"] = gm.gene.map(lambda g: maps_af[g][1] if g in maps_af else None)
    gm["mane_uniprot_map"] = gm.gene.map(lambda g: None if g not in maps_up else ("identity" if maps_up[g][0] is None else "alignment"))
    gm["mane_uniprot_identity"] = gm.gene.map(lambda g: maps_up[g][1] if g in maps_up else None)
    gm["af_fragmented_excluded"] = gm.af_n_fragments > 1

    def conv(store, g, p):
        if pd.isna(p) or g not in store:
            return None
        m = store[g][0]
        return int(p) if m is None else m.get(int(p))
    mp["af_pos"] = pd.array([conv(maps_af, g, p) for g, p in zip(mp.gene, mp.mane_pos)], dtype="Int64")
    mp["uniprot_pos"] = pd.array([conv(maps_up, g, p) for g, p in zip(mp.gene, mp.mane_pos)], dtype="Int64")
    gi = gm.set_index("gene")
    mp["uniprot"] = mp.gene.map(gi.uniprot)
    afs = gi.af_seq.to_dict(); ups = gi.uniprot_seq.to_dict()
    def chk(s, p, ref):
        if not isinstance(s, str) or pd.isna(p) or not isinstance(ref, str):
            return None
        return bool(1 <= int(p) <= len(s) and s[int(p) - 1] == ref)
    mp["ref_matches_af"] = pd.array([chk(afs.get(g), p, r) for g, p, r in zip(mp.gene, mp.af_pos, mp.mane_ref)], dtype="boolean")
    mp["ref_matches_uniprot"] = pd.array([chk(ups.get(g), p, r) for g, p, r in zip(mp.gene, mp.uniprot_pos, mp.mane_ref)], dtype="boolean")
    mp["af_single_fragment"] = mp.gene.map(gi.af_single_fragment)
    mp.to_parquet(INTER / "clinvar_mapped.parquet", index=False)
    gm.to_parquet(INTER / "gene_protein_map.parquet", index=False)

    ok = mp[mp.map_status == "ok"]
    onaf = ok[ok.af_single_fragment == True]
    log = dict(log)
    log.update(
        variants=len(mp), mapped_to_mane=len(ok), genes_mapped=int(ok.gene.nunique()),
        genes_with_uniprot=int(gm.uniprot.notna().sum()), uniprot_map_method={k: int(v) for k, v in gm.uniprot_map.value_counts().items()},
        genes_af_single=int(gm.af_single_fragment.sum()), genes_af_fragmented=int((gm.af_n_fragments > 1).sum()),
        genes_af_missing=int(((gm.af_n_fragments == 0) & gm.uniprot.notna()).sum()),
        fragmented_genes=sorted(gm[gm.af_n_fragments > 1].gene), af_missing_genes=sorted(gm[(gm.af_n_fragments == 0) & gm.uniprot.notna()].gene),
        genes_mane_ne_af=sorted(gm[gm.mane_eq_af == False].gene), genes_mane_ne_uniprot=sorted(gm[gm.mane_eq_uniprot == False].gene),
        variants_on_single_fragment_af=len(onaf), variants_af_pos_mapped=int(onaf.af_pos.notna().sum()),
        ref_matches_af=int(onaf.ref_matches_af.sum()), ref_mismatches_af=int((onaf.ref_matches_af == False).sum()),
        ref_match_rate_af=float(onaf.ref_matches_af.mean()),
        ref_match_rate_uniprot=float(ok.ref_matches_uniprot.mean()),
        variants_in_fragmented_genes=int((ok.af_single_fragment == False).sum()))
    save_json(dict(log), INTER / "protein_mapping_log.json")
    print(json.dumps({k: v for k, v in log.items() if k not in ("fragmented_genes", "af_missing_genes", "genes_mane_ne_af", "genes_mane_ne_uniprot")}, indent=1))


if __name__ == "__main__":
    main()
