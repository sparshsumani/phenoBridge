#!/usr/bin/env python3
"""Step 3: Ensembl VEP (GRCh38 REST) consequence annotation of every vcfRecord variant,
restricted to the gene's MANE Select transcript (fallback: Ensembl canonical), plus
exon structure lookup of those transcripts for the NMD-escape rule.

All REST responses are cached under CACHE/vep_region/ (one file per variant) and
CACHE/ens_lookup/ (one per transcript), so reruns are offline.

Outputs INTER/vep_annotations.parquet (one row per unique chrom/pos/ref/alt/hgnc_id),
        INTER/transcripts.parquet, INTER/vep_log.json
"""
import json, re, hashlib
from collections import Counter
import numpy as np
import pandas as pd
from common import CACHE, INTER, http_json

VEP_URL = ("https://rest.ensembl.org/vep/human/region?mane=1&canonical=1&hgvs=1&numbers=1"
           "&uniprot=1&protein=1&variant_class=1")
LOOKUP_URL = "https://rest.ensembl.org/lookup/id?expand=1"
VEPDIR = CACHE / "vep_region"; VEPDIR.mkdir(parents=True, exist_ok=True)
LKDIR = CACHE / "ens_lookup"; LKDIR.mkdir(parents=True, exist_ok=True)

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U"}


def vkey(s):
    return hashlib.sha1(s.encode()).hexdigest()


def vep_inputs(df):
    return [f"{c} {p} . {r} {a} . . ." for c, p, r, a in zip(df.chrom, df.pos, df.ref, df.alt)]


def run_vep(inputs):
    todo = [s for s in inputs if not (VEPDIR / (vkey(s) + ".json")).exists()]
    print(f"VEP: {len(inputs)} inputs, {len(todo)} not cached")

    def post(batch):
        res = http_json(VEP_URL, "POST", {"variants": batch})
        if isinstance(res, dict) and "__error__" in res:
            if len(batch) == 1:
                (VEPDIR / (vkey(batch[0]) + ".json")).write_text(json.dumps(res))
                return
            for s in batch:  # split to isolate bad record
                post([s])
            return
        got = set()
        for r in res:
            s = r["input"]
            got.add(s)
            (VEPDIR / (vkey(s) + ".json")).write_text(json.dumps(r))
        for s in batch:
            if s not in got:
                (VEPDIR / (vkey(s) + ".json")).write_text(json.dumps({"__error__": "no_result", "input": s}))
    from concurrent.futures import ThreadPoolExecutor
    batches = [todo[i:i + 200] for i in range(0, len(todo), 200)]
    with ThreadPoolExecutor(4) as ex:  # Ensembl allows ~15 req/s; 4 concurrent POSTs is polite
        for k, _ in enumerate(ex.map(post, batches)):
            print("  vep batch", k + 1, "/", len(batches), flush=True)
    return {s: json.loads((VEPDIR / (vkey(s) + ".json")).read_text()) for s in inputs}


def lookup_transcripts(ids):
    ids = sorted(set(i for i in ids if i))
    todo = [i for i in ids if not (LKDIR / f"{i}.json").exists()]
    for k in range(0, len(todo), 50):
        b = todo[k:k + 50]
        res = http_json(LOOKUP_URL, "POST", {"ids": b})
        for i in b:
            (LKDIR / f"{i}.json").write_text(json.dumps(res.get(i)))
    return {i: json.loads((LKDIR / f"{i}.json").read_text()) for i in ids}


def transcript_structure(lk):
    """exon lengths in transcript order, cDNA pos of CDS start, last exon-exon junction (cDNA)."""
    if not lk:
        return None
    strand = lk["strand"]
    exons = sorted(lk["Exon"], key=lambda e: e["start"], reverse=(strand == -1))
    lens = [e["end"] - e["start"] + 1 for e in exons]
    tr = lk.get("Translation")
    cds_start_cdna = None
    if tr:
        g = tr["start"] if strand == 1 else tr["end"]  # genomic coord of first CDS base
        acc = 0
        for e, L in zip(exons, lens):
            if e["start"] <= g <= e["end"]:
                off = (g - e["start"]) if strand == 1 else (e["end"] - g)
                cds_start_cdna = acc + off + 1
                break
            acc += L
    return dict(transcript_id=lk["id"], strand=strand, n_exons=len(exons), exon_lengths=lens,
                tx_length=sum(lens), cds_start_cdna=cds_start_cdna,
                last_junction_cdna=sum(lens[:-1]) if len(lens) > 1 else None,
                protein_length=tr["length"] if tr else None, protein_id=tr["id"] if tr else None,
                gene_id=lk.get("Parent"), display_name=lk.get("display_name"))


# ---- classification ----
CLASS_ORDER = [  # VEP severity order mapped to our classes
    ("transcript_ablation", "structural"),
    ("splice_acceptor_variant", "splice"), ("splice_donor_variant", "splice"),
    ("stop_gained", "nonsense"), ("frameshift_variant", "frameshift"),
    ("stop_lost", "stop_lost"), ("start_lost", "start_lost"),
    ("transcript_amplification", "structural"), ("feature_elongation", "structural"),
    ("feature_truncation", "structural"),
    ("inframe_insertion", "inframe_indel"), ("inframe_deletion", "inframe_indel"),
    ("missense_variant", "missense"), ("protein_altering_variant", "inframe_indel"),
    ("splice_donor_5th_base_variant", "splice"), ("splice_region_variant", "splice_region"),
    ("splice_donor_region_variant", "splice_region"), ("splice_polypyrimidine_tract_variant", "splice_region"),
    ("incomplete_terminal_codon_variant", "other"), ("start_retained_variant", "synonymous"),
    ("stop_retained_variant", "synonymous"), ("synonymous_variant", "synonymous"),
]


def classify(terms):
    if not terms:
        return "other"
    for t, c in CLASS_ORDER:
        if t in terms:
            if c == "splice_region":
                return "splice"  # splice-region hit without any coding change above it
            return c
    return "other"  # UTR / intron / non-coding transcript / upstream / downstream


def parse_p(p):
    """'NP_x:p.(Arg428Gln)' or 'ENSP..:p.Arg428Gln' -> (first_pos, ref1, alt1-or-None, kind, fs_ter)."""
    if not p or not isinstance(p, str):
        return None
    s = p.split(":", 1)[-1].replace("(", "").replace(")", "")
    if s in ("p.?", "p.0", "p.0?", "p.="):
        return None
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)(.*)$", s)
    if not m:
        return None
    ref, pos, rest = m.group(1), int(m.group(2)), m.group(3)
    fs = re.search(r"fsTer(\d+|\?)", rest)
    alt = None
    m2 = re.fullmatch(r"([A-Z][a-z]{2})", rest)
    if m2:
        alt = m2.group(1)
    return dict(pos=pos, ref=AA3.get(ref, "?"), alt=AA3.get(alt) if alt else None, rest=rest,
                fs_ter=(None if not fs or fs.group(1) == "?" else int(fs.group(1))))


def ptc_codon(cls, tc, pp):
    """Codon index (1-based, CDS) of the premature stop for nonsense / frameshift."""
    if pp is None:
        return None
    if cls == "nonsense":
        if pp["rest"] == "Ter":
            return pp["pos"]
        m = re.search(r"delins((?:[A-Z][a-z]{2})*?)Ter", pp["rest"])  # e.g. Arg595delinsSerTer
        if m:
            return pp["pos"] + len(m.group(1)) // 3
        m = re.match(r"_[A-Z][a-z]{2}\d+delins((?:[A-Z][a-z]{2})*?)Ter", pp["rest"])
        if m:
            return pp["pos"] + len(m.group(1)) // 3
        if pp["fs_ter"] is not None:
            return pp["pos"] + pp["fs_ter"] - 1
        return None
    if cls == "frameshift":
        if pp["fs_ter"] is not None:
            return pp["pos"] + pp["fs_ter"] - 1
        if pp["rest"] == "Ter":
            return pp["pos"]
    return None


def main():
    var = pd.read_parquet(INTER / "variants_raw.parquet")
    snv = var[~var.is_structural].drop_duplicates(["chrom", "pos", "ref", "alt", "hgnc_id"]).copy()
    snv["vep_input"] = vep_inputs(snv)
    vep = run_vep(sorted(set(snv.vep_input)))
    log = Counter()
    rows = []
    for _, r in snv.iterrows():
        res = vep[r.vep_input]
        out = dict(chrom=r.chrom, pos=r.pos, ref=r.ref, alt=r.alt, hgnc_id=r.hgnc_id, vep_input=r.vep_input)
        if "__error__" in res:
            out["vep_status"] = "error"
            out["vep_error"] = str(res.get("text", res.get("__error__")))[:300]
            log["vep_error"] += 1
            rows.append(out); continue
        tcs = [t for t in res.get("transcript_consequences", []) if t.get("hgnc_id") == r.hgnc_id]
        pick, src = None, None
        for t in tcs:
            if t.get("mane_select"):
                pick, src = t, "MANE_Select"; break
        if pick is None:
            for t in tcs:
                if t.get("canonical") == 1:
                    pick, src = t, "canonical"; break
        out["vep_most_severe_any"] = res.get("most_severe_consequence")
        out["vep_variant_class"] = res.get("variant_class")
        if pick is None:
            out["vep_status"] = "no_transcript_for_gene"
            log["no_transcript_for_gene"] += 1
            rows.append(out); continue
        out["vep_status"] = "ok"
        out["transcript_source"] = src
        out["ensembl_transcript"] = pick["transcript_id"]
        out["mane_select_refseq"] = pick.get("mane_select")
        out["ensembl_protein"] = pick.get("protein_id")
        out["vep_gene_symbol"] = pick.get("gene_symbol")
        out["biotype"] = pick.get("biotype")
        out["consequence_terms"] = ",".join(pick.get("consequence_terms", []))
        out["impact"] = pick.get("impact")
        out["vep_hgvsc"] = pick.get("hgvsc")
        out["vep_hgvsp"] = pick.get("hgvsp")
        out["vep_protein_start"] = pick.get("protein_start")
        out["vep_protein_end"] = pick.get("protein_end")
        out["vep_amino_acids"] = pick.get("amino_acids")
        out["vep_cds_start"] = pick.get("cds_start")
        out["vep_cdna_start"] = pick.get("cdna_start")
        out["exon"] = pick.get("exon")
        out["intron"] = pick.get("intron")
        sp = pick.get("swissprot") or []
        out["vep_swissprot"] = ",".join(s.split(".")[0] for s in sp)
        out["sift_score"] = pick.get("sift_score")
        out["polyphen_score"] = pick.get("polyphen_score")
        rows.append(out)
    ann = pd.DataFrame(rows)
    ann["variant_class"] = [classify(set(t.split(","))) if isinstance(t, str) else None for t in ann.consequence_terms]
    # size-based rules: >=50 bp variants that cross an exon/intron boundary are exon-level SVs
    ann["indel_length"] = (ann.ref.str.len() - ann.alt.str.len()).abs()
    ann["is_large_indel"] = ann.indel_length >= 50
    span = ann.consequence_terms.fillna("").str.contains("splice_acceptor_variant|splice_donor_variant|transcript_ablation")
    ann.loc[ann.is_large_indel & span, "variant_class"] = "structural"
    ann["structural_reason"] = np.where(ann.is_large_indel & span, "vcf_>=50bp_crossing_exon_boundary", None)
    # malformed / symbolic vcfRecords (e.g. ref 'N', alt 'DEL') are not annotatable -> structural
    sym = ~(ann.ref.str.fullmatch("[ACGTacgt]+") & ann.alt.str.fullmatch("[ACGTacgt]+"))
    ann.loc[sym, "variant_class"] = "structural"
    ann.loc[sym, "vep_status"] = "symbolic_or_malformed_vcf"
    ann.loc[sym, "structural_reason"] = "symbolic_or_malformed_vcf"
    log["symbolic_or_malformed_vcf"] = int(sym.sum())

    # transcripts / exon structure
    lk = lookup_transcripts(ann.ensembl_transcript.dropna().unique())
    tx = pd.DataFrame([s for s in (transcript_structure(lk[i]) for i in lk) if s])
    txd = tx.set_index("transcript_id").to_dict("index")

    # NMD escape
    ptc, dist, esc, last_exon, startprox = [], [], [], [], []
    vpp = [parse_p(p) for p in ann.vep_hgvsp]
    for cls, tid, pp in zip(ann.variant_class, ann.ensembl_transcript, vpp):
        c = ptc_codon(cls, None, pp) if cls in ("nonsense", "frameshift") else None
        t = txd.get(tid)
        if c is None or t is None or t["cds_start_cdna"] is None:
            ptc.append(np.nan); dist.append(np.nan); esc.append(None); last_exon.append(None); startprox.append(None)
            continue
        ptc_cds = 3 * (c - 1) + 1                   # first nt of the premature stop codon, CDS coords
        ptc_cdna = t["cds_start_cdna"] + ptc_cds - 1
        ptc.append(c)
        startprox.append(ptc_cds <= 150)
        if t["last_junction_cdna"] is None:          # single-exon transcript
            dist.append(np.nan); esc.append(True); last_exon.append(True); continue
        d = t["last_junction_cdna"] - ptc_cdna     # nt from PTC start to last exon-exon junction
        dist.append(d)
        last_exon.append(d < 0)
        esc.append(bool(d <= 55))                  # in last exon (d<0) or <=55 nt upstream of last EJC
    ann["ptc_codon"] = ptc
    ann["ptc_dist_to_last_junction_nt"] = dist
    ann["ptc_in_last_exon"] = last_exon
    ann["nmd_escape"] = esc
    ann["ptc_start_proximal_150nt"] = startprox
    ann["vep_p_pos"] = [p["pos"] if p else np.nan for p in vpp]
    ann["vep_p_ref"] = [p["ref"] if p else None for p in vpp]
    ann["vep_p_alt"] = [p["alt"] if p else None for p in vpp]
    for c in ("vep_protein_start", "vep_protein_end", "vep_cds_start", "vep_cdna_start", "ptc_codon", "vep_p_pos"):
        ann[c] = pd.to_numeric(ann[c], errors="coerce").astype("Int64")
    ann.to_parquet(INTER / "vep_annotations.parquet", index=False)
    tx.to_parquet(INTER / "transcripts.parquet", index=False)
    log.update({f"class_{k}": int(v) for k, v in ann.variant_class.value_counts(dropna=False).items()})
    log.update({f"source_{k}": int(v) for k, v in ann.transcript_source.value_counts(dropna=False).items()})
    log["unique_variants"] = len(ann)
    (INTER / "vep_log.json").write_text(json.dumps(dict(log), indent=1))
    print(json.dumps(dict(log), indent=1))


if __name__ == "__main__":
    main()
