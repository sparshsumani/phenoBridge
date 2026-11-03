#!/usr/bin/env python3
"""Step 4: map every train/test variant to the MANE Select protein, the reviewed UniProt canonical
isoform and the AlphaFold v6 model numbering (allelic/a04 logic; a04's load_mane / load_uniprot /
vep_pick, multigene s03 run_vep (cache/vep_region) and s06 align_map / af_sequence are reused).

  * protein change from the ClinVar Name of the label-defining release (train: 2023-01, test: 2026-11)
  * 'clinvar_tx_is_mane': Name transcript (version ignored) == MANE RefSeq_nuc and the ref aa equals the
    MANE protein residue -> MANE position = ClinVar position
  * else 'vep_mane': Ensembl VEP REST GRCh38 on the VCF coordinates; MANE Select consequence of the same
    gene must be missense_variant (others dropped: map_status vep_not_missense_on_mane / ...)
  * validation: 300 random 'clinvar_tx_is_mane' variants also run through VEP
  * gene -> UniProt (MANE-Select xref, else HGNC xref, else gene name) -> AlphaFold v6 fragments from the
    proteome tar listing; single-fragment models are extracted into raw/alphafold/models/ (tar -T list,
    only the needed AF-*-F1-model_v6.pdb.gz members). Fragmented (>2,700 aa) proteins are excluded.
  * MANE pos -> AF / UniProt pos: identity when sequences are equal, else s06 align_map (blocks >=10 aa,
    >=90 % identity). ref_matches_af: MANE ref aa == AF model residue at the mapped position.
  * protein-level consistency: when the variant has a parsable missense Name in both releases, the
    2023 and 2026 hgvs_p are compared (logged).
Outputs intermediate/variants_mapped.parquet, intermediate/gene_protein_map.parquet, intermediate/protein_mapping_log.json
"""
import re, sys, json, subprocess
from collections import Counter, defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from common import INTER, AF_TAR, AF_DIR, AF_LISTING, save_json, al_import, mg_import, af_model_path

a04 = al_import("a04_protein_mapping")
s03 = mg_import("s03_vep_annotate")
s06 = mg_import("s06_select_genes")


def main():
    log = Counter()
    s = pd.read_parquet(INTER / "sets.parquet")
    mane = a04.load_mane()
    by_sym = mane.drop_duplicates("symbol").set_index("symbol", drop=False)
    by_hgnc = mane.drop_duplicates("HGNC_ID").set_index("HGNC_ID", drop=False)
    seqs = dict(zip(mane.RefSeq_prot, mane.mane_protein_seq))
    rows = []
    for r in s.itertuples():
        mr = by_sym.loc[r.name_gene] if r.name_gene in by_sym.index else None
        how = "symbol"
        if mr is None:
            for h in str(r.clinvar_hgnc_ids).split(";"):
                if h in by_hgnc.index:
                    mr = by_hgnc.loc[h]; how = "hgnc"; break
        d = dict(variation_id=r.variation_id, gene=r.name_gene, mane_lookup=how if mr is not None else None)
        if mr is None:
            d["map_status"] = "no_mane_select_for_gene"; rows.append(d); continue
        d.update(gene_mane=mr.symbol, hgnc_id=mr.HGNC_ID, mane_nm=mr.RefSeq_nuc, mane_np=mr.RefSeq_prot,
                 mane_enst=mr.Ensembl_nuc, mane_ensp=mr.Ensembl_prot)
        seq = mr.mane_protein_seq
        d["clinvar_tx_is_mane"] = r.name_tx_noversion == mr.nm_nov
        p = r.clinvar_protein_pos
        if d["clinvar_tx_is_mane"] and isinstance(seq, str) and 1 <= p <= len(seq) and seq[p - 1] == r.ref_aa:
            d.update(map_status="ok", map_method="clinvar_tx_is_mane", mane_pos=p, mane_ref=r.ref_aa, mane_alt=r.alt_aa)
        else:
            d["map_status"] = "needs_vep"
            d["needs_vep_reason"] = "non_mane_transcript" if not d["clinvar_tx_is_mane"] else "ref_mismatch_on_mane_seq"
        rows.append(d)
    mp = s[["variation_id", "set", "chrom", "pos_vcf", "ref_vcf", "alt_vcf", "name_tx", "clinvar_protein_pos", "ref_aa", "alt_aa",
            "hgvs_p", "hgvs_c"]].merge(pd.DataFrame(rows), on="variation_id")
    assert len(mp) == len(s)
    log.update({f"pre_vep_{k}": int(v) for k, v in mp.map_status.value_counts().items()})
    log.update({f"needs_vep_{k}": int(v) for k, v in mp.needs_vep_reason.value_counts().items()})

    need = mp.map_status == "needs_vep"
    rng = np.random.RandomState(0)
    ok_idx = mp.index[mp.map_status == "ok"]
    val_idx = rng.choice(ok_idx, size=min(300, len(ok_idx)), replace=False)
    vcf_ok = mp.pos_vcf.str.fullmatch(r"\d+") & mp.ref_vcf.str.fullmatch("[ACGT]") & mp.alt_vcf.str.fullmatch("[ACGT]")
    mp["vep_input"] = [f"{c} {p} . {a} {b} . . ." if ok else None for c, p, a, b, ok in zip(mp.chrom, mp.pos_vcf, mp.ref_vcf, mp.alt_vcf, vcf_ok)]
    run_idx = sorted(set(mp.index[need & vcf_ok]) | set(val_idx))
    log["vep_inputs"] = len(run_idx)
    log["needs_vep_but_no_vcf_coords"] = int((need & ~vcf_ok).sum())
    vep = s03.run_vep(sorted(set(mp.loc[run_idx, "vep_input"])))
    val = Counter(); val_set = set(val_idx)
    upd = {}
    for i in run_idx:
        r = mp.loc[i]
        pick, st = a04.vep_pick(vep.get(r.vep_input), r.gene, r.hgnc_id)
        vpos = vref = valt = None
        if st == "ok":
            vref, valt = pick["amino_acids"].split("/"); vpos = int(pick["protein_start"])
        if i in val_set and r.map_status == "ok":
            val["validated"] += 1
            val["agree" if (st == "ok" and vpos == r.mane_pos and vref == r.mane_ref and valt == r.mane_alt) else "disagree:" + st] += 1
            continue
        if st == "ok":
            upd[i] = dict(map_status="ok", map_method="vep_mane", mane_pos=vpos, mane_ref=vref, mane_alt=valt,
                          vep_hgvsp=pick.get("hgvsp"))
        else:
            upd[i] = dict(map_status=st.split(":")[0], vep_status_detail=st, vep_hgvsp=(pick or {}).get("hgvsp"))
    if upd:
        u = pd.DataFrame.from_dict(upd, orient="index")
        for c in u.columns:
            if c not in mp:
                mp[c] = None
            mp.loc[u.index, c] = u[c]
    mp.loc[need & ~vcf_ok, "map_status"] = "no_vcf_coordinates"
    log["vep_validation_on_mane_transcript_variants"] = dict(val)
    mp["mane_pos"] = pd.to_numeric(mp.mane_pos, errors="coerce").astype("Int64")
    mp["ref_matches_mane_seq"] = pd.array([
        (bool(1 <= int(p) <= len(seqs.get(np_, "") or "") and seqs[np_][int(p) - 1] == ref) if pd.notna(p) and isinstance(np_, str) else None)
        for p, np_, ref in zip(mp.mane_pos, mp.mane_np, mp.mane_ref)], dtype="boolean")
    log.update({f"post_vep_{k}": int(v) for k, v in mp.map_status.value_counts().items()})
    log["map_method"] = {k: int(v) for k, v in mp.map_method.value_counts().items()}
    log["map_status_by_set"] = {f"{a}|{b}": int(n) for (a, b), n in mp.groupby(["set", "map_status"]).size().items()}
    log["vep_detail"] = {k: int(v) for k, v in mp.vep_status_detail.value_counts().head(15).items()} if "vep_status_detail" in mp else {}

    # ---- genes -> UniProt -> AlphaFold (a04 logic, all genes of mapped variants) ----
    up = a04.load_uniprot()
    up["mane_enst_nov"] = up.mane_select.fillna("").map(lambda x: [y.split(".")[0] for y in x.split(";") if y.strip()])
    up["hgnc_list"] = up.hgnc.fillna("").map(lambda x: [y for y in x.split(";") if y.strip()])
    by_enst, by_h, by_g = defaultdict(list), defaultdict(list), defaultdict(list)
    for k, r in up.iterrows():
        for e in r.mane_enst_nov:
            by_enst[e].append(k)
        for h in r.hgnc_list:
            by_h[h].append(k)
        if isinstance(r.gene_primary, str):
            by_g[r.gene_primary].append(k)
    names = AF_LISTING.read_text().split()
    frags = defaultdict(set)
    for n in names:
        m = re.match(r"AF-([A-Z0-9]+)-F(\d+)-model_v(\d+)\.pdb\.gz$", n)
        if m:
            frags[m.group(1)].add(int(m.group(2)))
    genes = mp[mp.map_status == "ok"].drop_duplicates("gene")[["gene", "hgnc_id", "mane_enst", "mane_np", "mane_ensp"]]
    grows = []
    for g in genes.itertuples():
        e = g.mane_enst.split(".")[0]
        cand, how = by_enst.get(e, []), "uniprot_mane_xref"
        if not cand:
            cand, how = by_h.get(g.hgnc_id, []), "hgnc_xref"
        if not cand:
            cand, how = by_g.get(g.gene, []), "gene_name"
        if not cand:
            grows.append(dict(gene=g.gene, hgnc_id=g.hgnc_id, uniprot=None, uniprot_map=None, af_n_fragments=0)); continue
        k = cand[0]
        if len(cand) > 1:
            ms = seqs.get(g.mane_np)
            k = next((c for c in cand if up.loc[c, "sequence"] == ms), cand[0])
        u = up.loc[k]
        grows.append(dict(gene=g.gene, hgnc_id=g.hgnc_id, mane_enst=g.mane_enst, mane_np=g.mane_np, mane_ensp=g.mane_ensp,
                          uniprot=u.accession, uniprot_map=how, n_uniprot_candidates=len(cand), uniprot_entry=u.entry_name,
                          uniprot_length=int(u.length), uniprot_seq=u.sequence, uniprot_seq_version=u.seq_version,
                          mane_protein_seq=seqs.get(g.mane_np), af_n_fragments=len(frags.get(u.accession, ()))))
    gm = pd.DataFrame(grows)
    gm["af_source"] = np.where(gm.af_n_fragments > 0, "proteome_tar", None)
    sw = (gm.af_n_fragments == 0) & gm.uniprot.map(lambda a: isinstance(a, str) and af_model_path(a).exists())
    gm.loc[sw, "af_n_fragments"] = 1
    gm.loc[sw, "af_source"] = "swissprot_tar_range_request(multigene)"
    log["af_models_from_swissprot_tar"] = sorted(gm[sw].gene)
    gm["af_single_fragment"] = gm.af_n_fragments == 1
    todo = sorted({f"AF-{a}-F1-model_v6.pdb.gz" for a in gm[gm.af_single_fragment].uniprot if not af_model_path(a).exists()})
    log["af_models_extracted_now"] = len(todo)
    if todo:
        lst = INTER / "af_extract_list.txt"
        lst.write_text("\n".join(todo) + "\n")
        print("extracting", len(todo), "AF models", flush=True)
        subprocess.run(["tar", "-xf", str(AF_TAR), "-C", str(AF_DIR), "-T", str(lst)], check=True)
    gm["af_pdb"] = [str(af_model_path(a)) if s_ and af_model_path(a).exists() else None for a, s_ in zip(gm.uniprot, gm.af_single_fragment)]
    gm["af_seq"] = [s06.af_sequence(p) if isinstance(p, str) else None for p in gm.af_pdb]
    eq = lambda a, b: (a == b) if isinstance(a, str) and isinstance(b, str) else None
    gm["af_eq_uniprot"] = [eq(a, b) for a, b in zip(gm.af_seq, gm.uniprot_seq)]
    gm["mane_eq_uniprot"] = [eq(a, b) for a, b in zip(gm.mane_protein_seq, gm.uniprot_seq)]
    gm["mane_eq_af"] = [eq(a, b) for a, b in zip(gm.mane_protein_seq, gm.af_seq)]
    maps_af, maps_up = {}, {}
    for r in gm.itertuples():
        for tgt, store in ((r.af_seq, maps_af), (r.uniprot_seq, maps_up)):
            if isinstance(r.mane_protein_seq, str) and isinstance(tgt, str):
                store[r.gene] = (None, 1.0) if r.mane_protein_seq == tgt else s06.align_map(r.mane_protein_seq, tgt)
    gm["mane_af_map"] = gm.gene.map(lambda g: None if g not in maps_af else ("identity" if maps_af[g][0] is None else "alignment"))
    gm["mane_af_identity"] = gm.gene.map(lambda g: maps_af[g][1] if g in maps_af else None)
    gm["mane_uniprot_map"] = gm.gene.map(lambda g: None if g not in maps_up else ("identity" if maps_up[g][0] is None else "alignment"))
    gm["mane_uniprot_identity"] = gm.gene.map(lambda g: maps_up[g][1] if g in maps_up else None)

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
    def chk(sq, p, ref):
        if not isinstance(sq, str) or pd.isna(p) or not isinstance(ref, str):
            return None
        return bool(1 <= int(p) <= len(sq) and sq[int(p) - 1] == ref)
    mp["ref_matches_af"] = pd.array([chk(afs.get(g), p, r) for g, p, r in zip(mp.gene, mp.af_pos, mp.mane_ref)], dtype="boolean")
    mp["ref_matches_uniprot"] = pd.array([chk(ups.get(g), p, r) for g, p, r in zip(mp.gene, mp.uniprot_pos, mp.mane_ref)], dtype="boolean")
    mp["af_n_fragments"] = mp.gene.map(gi.af_n_fragments)
    mp["af_single_fragment"] = mp.gene.map(gi.af_single_fragment).fillna(False).astype(bool)
    mp["protein_length_af"] = mp.gene.map(lambda g: len(afs[g]) if isinstance(afs.get(g), str) else np.nan)
    mp["usable"] = ((mp.map_status == "ok") & mp.af_single_fragment & mp.af_pos.notna() & (mp.ref_matches_af == True)).astype(bool)
    # exclusion reason (first failing step)
    reason = np.select([mp.map_status != "ok", mp.uniprot.isna(), mp.af_n_fragments.fillna(0) == 0, mp.af_n_fragments.fillna(0) > 1,
                        mp.af_pos.isna(), mp.ref_matches_af != True],
                       ["map:" + mp.map_status.astype(str), "no_uniprot", "no_af_model", "af_fragmented", "af_pos_unaligned", "ref_mismatch_af"], "usable")
    mp["map_exclusion"] = reason
    mp.drop(columns=["vep_input"]).to_parquet(INTER / "variants_mapped.parquet", index=False)
    gm.to_parquet(INTER / "gene_protein_map.parquet", index=False)

    ok = mp[mp.map_status == "ok"]; onaf = ok[ok.af_single_fragment]
    log = dict(log)
    log.update(variants=len(mp), mapped_to_mane=len(ok), genes_mapped=int(ok.gene.nunique()),
               genes_with_uniprot=int(gm.uniprot.notna().sum()), uniprot_map_method={k: int(v) for k, v in gm.uniprot_map.value_counts().items()},
               genes_af_single=int(gm.af_single_fragment.sum()), genes_af_fragmented=int((gm.af_n_fragments > 1).sum()),
               genes_af_missing=int(((gm.af_n_fragments == 0) & gm.uniprot.notna()).sum()),
               fragmented_genes=sorted(gm[gm.af_n_fragments > 1].gene), af_missing_genes=sorted(gm[(gm.af_n_fragments == 0) & gm.uniprot.notna()].gene),
               n_genes_mane_ne_af=int((gm.mane_eq_af == False).sum()),
               variants_on_single_fragment_af=len(onaf), variants_af_pos_mapped=int(onaf.af_pos.notna().sum()),
               ref_matches_af=int(onaf.ref_matches_af.sum()), ref_mismatches_af=int((onaf.ref_matches_af == False).sum()),
               ref_match_rate_af_among_mapped=float(onaf.ref_matches_af.mean()),
               ref_match_rate_af_among_single_fragment=float(onaf.ref_matches_af.sum() / len(onaf)),
               ref_match_rate_uniprot=float(ok.ref_matches_uniprot.mean()),
               exclusion_by_set_label={f"{a}|{b}|{c}": int(n) for (a, b, c), n in mp.merge(pd.read_parquet(INTER / "sets.parquet", columns=["variation_id", "label"]), on="variation_id").groupby(["set", "label", "map_exclusion"]).size().items()},
               ref_mismatch_examples=mp[(mp.ref_matches_af == False) & mp.af_single_fragment].groupby("gene").size().sort_values(ascending=False).head(15).to_dict())
    save_json(log, INTER / "protein_mapping_log.json")
    print(json.dumps({k: v for k, v in log.items() if k not in ("fragmented_genes", "af_missing_genes")}, indent=1, default=str))


if __name__ == "__main__":
    main()
