"""Step 1: variant annotation.

* Links every individual in ANKRD11_KBGS_individuals.json to its GA4GH phenopacket and verifies the
  match (allele key, sex, all 48 HPO calls).
* Runs Ensembl VEP (REST, cached) on the hg38 VCF records and keeps the MANE Select transcript
  ENST00000301030 (= NM_013275.6).
* Gets the exon structure of ENST00000301030 and computes, for each premature termination codon
  (PTC), its position relative to the last exon-exon junction (NMD rules).
* Parses the structural variants (copy number, size, overlap with ANKRD11 in hg19/hg38).
* Compares with the repo's regex truncation flag and c.-based residue mapping.
Output: out/variants_annotated.csv, out/exons.csv, out/annotation_summary.json
"""
import glob
import json
import re

import numpy as np
import pandas as pd

from common import (CACHE, CDS_LEN, COHORT_JSON, MANE_ENST, OUT, PPKT_DIR, PROT_LEN,
                    cached_get, cached_post, import_mvp)

M = import_mvp()
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*"}

# ---------------------------------------------------------------- cohort + phenopackets
cohort = json.loads(COHORT_JSON.read_text())
H = [h["hpoId"] for h in cohort["hpoHeaders"]]
rows = {r["individualData"]["individualId"]: r for r in cohort["rows"]}
SEX = {"MALE": "M", "FEMALE": "F", "UNKNOWN_SEX": "U"}

recs, mismatch = [], []
for f in sorted(glob.glob(str(PPKT_DIR / "*.json"))):
    pp = json.load(open(f))
    sid = pp["subject"]["id"]
    r = rows.get(sid)
    if r is None:
        mismatch.append((sid, "no JSON row")); continue
    gi = pp["interpretations"][0]["diagnosis"]["genomicInterpretations"]
    assert len(gi) == 1
    vd = gi[0]["variantInterpretation"]["variationDescriptor"]
    akey = list(r["alleleCountMap"])[0]
    feats = {p["type"]["id"]: 0 if p.get("excluded") else 1 for p in pp.get("phenotypicFeatures", [])}
    jv = {h: {"Observed": 1, "Excluded": 0}.get(c["type"]) for h, c in zip(H, r["hpoData"])}
    jv = {h: v for h, v in jv.items() if v is not None}
    ok_allele = vd["id"] == akey
    ok_sex = SEX.get(pp["subject"].get("sex", "UNKNOWN_SEX")) == r["individualData"]["sex"]
    ok_ph = jv == feats
    if not (ok_allele and ok_sex and ok_ph):
        mismatch.append((sid, ok_allele, ok_sex, ok_ph))
    ex = {e["syntax"]: e["value"] for e in vd.get("expressions", [])}
    vcf = vd.get("vcfRecord", {})
    recs.append(dict(individualId=sid, ppkt_file=f.split("/")[-1], sex=r["individualData"]["sex"],
                     pmid=r["individualData"]["pmid"], alleleKey=akey,
                     hgvs_c=ex.get("hgvs.c"), hgvs_g=ex.get("hgvs.g"), hgvs_p=ex.get("hgvs.p"),
                     chrom=vcf.get("chrom"), pos=vcf.get("pos"), ref=vcf.get("ref"), alt=vcf.get("alt"),
                     sv_label=vd.get("label") if not ex else None))
V = pd.DataFrame(recs)
assert len(V) == len(rows) == 333, len(V)
match_summary = dict(n_json=len(rows), n_ppkt=len(V), n_mismatch=len(mismatch), mismatches=mismatch)
print("phenopacket match:", match_summary)

# ---------------------------------------------------------------- exon structure (MANE)
lk = cached_get(f"https://rest.ensembl.org/lookup/id/{MANE_ENST}", {"expand": 1}, name="lookup_mane")
tl = lk["Translation"]
assert lk["strand"] == -1 and tl["length"] == PROT_LEN
cds_g_start, cds_g_end = tl["end"], tl["start"]          # minus strand: CDS starts at the high coord
exons = sorted(lk["Exon"], key=lambda e: -e["end"])
ex_rows, c = [], 0
for i, e in enumerate(exons, 1):
    lo, hi = max(e["start"], tl["start"]), min(e["end"], tl["end"])
    ncod = max(0, hi - lo + 1)
    c_first = c + 1 if ncod else None
    c += ncod
    ex_rows.append(dict(exon=i, id=e["id"], g_start=e["start"], g_end=e["end"], length=e["end"] - e["start"] + 1,
                        coding_nt=ncod, c_first=c_first, c_last=c if ncod else None))
EX = pd.DataFrame(ex_rows)
assert c == CDS_LEN, c
EX["aa_first"] = np.ceil(EX.c_first / 3)
EX["aa_last"] = np.ceil(EX.c_last / 3)
EX.to_csv(OUT / "exons.csv", index=False)
last_coding_exon = int(EX.loc[EX.coding_nt > 0, "exon"].max())
assert last_coding_exon == EX.exon.max(), "stop codon not in last exon"
LAST_JUNCTION_C = int(EX.loc[EX.exon == last_coding_exon - 1, "c_last"].iloc[0])  # c. of last nt before final junction
print(EX.to_string())
print("last exon-exon junction after c.", LAST_JUNCTION_C)


def exon_of_c(cpos):
    if cpos is None or np.isnan(cpos):
        return np.nan
    hit = EX[(EX.c_first <= cpos) & (EX.c_last >= cpos)]
    return int(hit.exon.iloc[0]) if len(hit) else np.nan


# ---------------------------------------------------------------- VEP
pt = V.dropna(subset=["pos"]).copy()
pt["vep_in"] = pt.apply(lambda r: f"{r.chrom.replace('chr', '')} {int(r.pos)} {r.alleleKey} {r.ref} {r.alt} . . .", axis=1)
uniq_in = sorted(pt.vep_in.unique())
vep = []
for i in range(0, len(uniq_in), 150):
    vep += cached_post("https://rest.ensembl.org/vep/homo_sapiens/region", {"variants": uniq_in[i:i + 150]},
                       params={"mane": 1, "hgvs": 1, "numbers": 1, "canonical": 1, "shift_3prime": 1,
                               "transcript_version": 1},
                       name=f"vep_{i}")
vep_by_id = {}
for v in vep:
    tcs = [t for t in v.get("transcript_consequences", []) if t["transcript_id"].startswith(MANE_ENST)]
    vep_by_id[v["id"]] = (v, tcs[0] if tcs else None)
print("VEP records:", len(vep), "unique inputs:", len(uniq_in))


def parse_p(p):
    """Return (first_aa_pos, ptc_aa_pos, kind) from an HGVS p. string."""
    if not isinstance(p, str):
        return np.nan, np.nan, None
    s = p.split(":")[-1].replace("p.", "").strip("()")
    m = re.match(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})?fs(?:Ter|\*)(\d+|\?)", s)
    if m:
        first = int(m.group(2))
        ptc = first + int(m.group(4)) - 1 if m.group(4) != "?" else np.nan
        return first, ptc, "fs"
    m = re.match(r"([A-Z][a-z]{2})(\d+)(Ter|\*)$", s)
    if m:
        return int(m.group(2)), int(m.group(2)), "nonsense"
    m = re.match(r"([A-Z][a-z]{2})(\d+)(?:_[A-Z][a-z]{2}\d+)?delins((?:[A-Z][a-z]{2})*?)Ter$", s)
    if m:  # e.g. Arg595delinsSerTer: stop after the inserted residues
        return int(m.group(2)), int(m.group(2)) + len(m.group(3)) // 3, "nonsense"
    m = re.match(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", s)
    if m:
        return int(m.group(2)), np.nan, "missense"
    m = re.match(r"([A-Z][a-z]{2})(\d+)(?:_([A-Z][a-z]{2})(\d+))?(del|dup|ins|delins)", s)
    if m:
        return int(m.group(2)), np.nan, "inframe"
    m = re.match(r"([A-Z][a-z]{2})(\d+)=", s)
    if m:
        return int(m.group(2)), np.nan, "synonymous"
    return np.nan, np.nan, s


out = []
for _, r in V.iterrows():
    d = r.to_dict()
    d["regex_trunc"] = int(M.is_truncating(r.alleleKey))
    ap = M.allelekey_to_aa_pos(r.alleleKey)
    d["regex_aa_pos"] = ap if not r.alleleKey.startswith("ANKRD11_SV") else np.nan
    if pd.isna(r.pos):
        d["vclass"] = "SV"
        out.append(d); continue
    key = f"{r.chrom.replace('chr', '')} {int(r.pos)} {r.alleleKey} {r.ref} {r.alt} . . ."
    v, tc = vep_by_id[r.alleleKey]
    assert v["input"] == key
    d["vep_most_severe"] = v.get("most_severe_consequence")
    if tc is None:
        d["vclass"] = "no_MANE_consequence"; out.append(d); continue
    cons = tc["consequence_terms"]
    d["vep_consequence"] = "&".join(cons)
    d["vep_transcript"] = tc["transcript_id"]
    d["vep_mane"] = tc.get("mane_select")
    d["vep_hgvsc"] = tc.get("hgvsc")
    d["vep_hgvsp"] = tc.get("hgvsp")
    d["vep_exon"] = tc.get("exon")
    d["vep_intron"] = tc.get("intron")
    d["vep_protein_start"] = tc.get("protein_start")
    d["vep_protein_end"] = tc.get("protein_end")
    d["vep_cds_start"] = tc.get("cds_start")
    d["vep_aa"] = tc.get("amino_acids")
    if "frameshift_variant" in cons:
        vc = "frameshift"
    elif "stop_gained" in cons:
        vc = "nonsense"
    elif "splice_acceptor_variant" in cons or "splice_donor_variant" in cons:
        vc = "splice"
    elif "inframe_deletion" in cons or "inframe_insertion" in cons or "protein_altering_variant" in cons:
        vc = "inframe"
    elif "missense_variant" in cons:
        vc = "missense"
    elif "start_lost" in cons:
        vc = "start_lost"
    elif "stop_lost" in cons:
        vc = "stop_lost"
    else:
        vc = "other:" + "&".join(cons)
    d["vclass"] = vc
    # protein position + PTC from the phenopacket HGVS p. (cross-checked against VEP below)
    first, ptc, kind = parse_p(r.hgvs_p)
    d["pp_p_first"], d["pp_ptc_aa"], d["pp_p_kind"] = first, ptc, kind
    fv, ptcv, kindv = parse_p(tc.get("hgvsp"))
    d["vep_p_first"], d["vep_ptc_aa"], d["vep_p_kind"] = fv, ptcv, kindv
    out.append(d)
A = pd.DataFrame(out)

# protein position: VEP protein_start; for splice variants use the nearest codon of the cDNA position
def cpos_from_hgvsc(h):
    if not isinstance(h, str):
        return np.nan
    m = re.search(r"c\.(-?\d+)([+-]\d+)?", h)
    return float(m.group(1)) if m else np.nan


A["vep_c_anchor"] = A.vep_hgvsc.map(cpos_from_hgvsc)
A["aa_pos"] = A.vep_protein_start.astype(float)
spl = A.vclass.eq("splice") & A.aa_pos.isna()
A.loc[spl, "aa_pos"] = np.ceil(A.loc[spl, "vep_c_anchor"].clip(lower=1) / 3)
A["rel_pos"] = A.aa_pos / PROT_LEN

# PTC position (aa) and NMD prediction
A["ptc_aa"] = A.vep_ptc_aa.where(A.vep_ptc_aa.notna(), A.pp_ptc_aa)
A.loc[A.vclass.eq("nonsense") & A.ptc_aa.isna(), "ptc_aa"] = A.aa_pos
A.loc[~A.vclass.isin(["nonsense", "frameshift"]), "ptc_aa"] = np.nan
A["ptc_c_first"] = 3 * A.ptc_aa - 2                       # first nt of the new stop codon
A["ptc_exon"] = A.ptc_c_first.map(exon_of_c)
A["ptc_dist_to_last_junction"] = LAST_JUNCTION_C - A.ptc_c_first + 1   # nt upstream of last EJ (<=0: last exon)
for rule in (50, 55):
    A[f"nmd_escape_{rule}nt"] = np.where(A.ptc_aa.notna(), (A.ptc_dist_to_last_junction <= rule).astype(float), np.nan)
A["nmd_escape_start_prox150"] = np.where(A.ptc_aa.notna(), (A.ptc_c_first <= 150).astype(float), np.nan)
exlen = dict(zip(EX.exon, EX.length))
A["ptc_exon_len"] = A.ptc_exon.map(exlen)
A["ptc_in_long_exon_407"] = np.where(A.ptc_aa.notna(), (A.ptc_exon_len > 407).astype(float), np.nan)
A["retained_aa"] = A.ptc_aa - 1                            # residues of WT sequence + frameshift tail before PTC
A["trunc_first_aa"] = np.where(A.vclass.eq("frameshift"), A.aa_pos, A.ptc_aa)   # first WT residue lost

# ---------------------------------------------------------------- structural variants
def parse_sv(label):
    s = label.replace(" ", "")
    cn = re.search(r"x(\d)$", s)
    m = re.search(r"\((\d[\d.]*)[-_](\d[\d.]*)\)", s)
    if not m:
        return dict(sv_cn=np.nan, sv_start=np.nan, sv_end=np.nan)
    a, b = m.group(1), m.group(2)
    if "." in a:  # given in Mb
        a, b = float(a) * 1e6, float(b) * 1e6
    return dict(sv_cn=int(cn.group(1)) if cn else np.nan, sv_start=float(a), sv_end=float(b))


GENE = {"hg38": (89267619, 89490706), "hg19": (89334029, 89556969)}   # ENSG00000167522 (Ensembl 38 / 37)
CDS_G = {"hg38": (89268478, 89317019)}
svp = A.loc[A.vclass.eq("SV"), "sv_label"].map(parse_sv).apply(pd.Series)
A = A.join(svp)
A["sv_size_kb"] = (A.sv_end - A.sv_start) / 1e3
for asm, (gs, ge) in GENE.items():
    ov = (np.minimum(A.sv_end, ge) - np.maximum(A.sv_start, gs)).clip(lower=0)
    A[f"sv_gene_overlap_{asm}"] = np.where(A.vclass.eq("SV"), ov / (ge - gs), np.nan)
A["sv_type"] = np.where(A.vclass.eq("SV"), np.select([A.sv_cn == 1, A.sv_cn == 3], ["deletion", "duplication"],
                                                    "unknown"), None)

# ---------------------------------------------------------------- regex flag vs VEP
TRUNC = {"frameshift", "nonsense", "splice"}
A["is_trunc_vep"] = np.where(A.vclass.eq("SV"), np.nan, A.vclass.isin(TRUNC).astype(float))
A["vep_vs_pp_ptc_agree"] = (A.vep_ptc_aa == A.pp_ptc_aa) | (A.vep_ptc_aa.isna() & A.pp_ptc_aa.isna())
A["vep_vs_pp_first_agree"] = (A.vep_p_first == A.pp_p_first) | (A.vep_p_first.isna() & A.pp_p_first.isna())
A["regex_pos_diff"] = A.regex_aa_pos - A.aa_pos

A.to_csv(OUT / "variants_annotated.csv", index=False)

pt = A[A.vclass.ne("SV")]
cmp = pd.crosstab(pt.vclass, pt.regex_trunc)
print(cmp)
summ = dict(
    match=match_summary,
    mane_transcript=lk["id"] + "." + str(lk["version"]),
    vep_mane_select=sorted(pt.vep_mane.dropna().unique().tolist()),
    n_exons=int(len(EX)), coding_exons=EX.loc[EX.coding_nt > 0, "exon"].tolist(),
    largest_exon=int(EX.loc[EX.length.idxmax(), "exon"]), largest_exon_len=int(EX.length.max()),
    largest_exon_coding_frac=float(EX.coding_nt.max() / CDS_LEN),
    last_exon_coding_nt=int(EX.coding_nt.iloc[-1]), last_junction_after_c=LAST_JUNCTION_C,
    n_points=int(len(pt)), n_sv=int(A.vclass.eq("SV").sum()),
    vclass_counts=pt.vclass.value_counts().to_dict(),
    sv_type_counts=A.sv_type.value_counts().to_dict(),
    regex_vs_vep=cmp.to_dict(),
    vep_ptc_disagree_pp=int((~pt.vep_vs_pp_ptc_agree).sum()),
    vep_first_disagree_pp=int((~pt.vep_vs_pp_first_agree).sum()),
)
(OUT / "annotation_summary.json").write_text(json.dumps(summ, indent=1, default=str))
print(json.dumps(summ, indent=1, default=str))
