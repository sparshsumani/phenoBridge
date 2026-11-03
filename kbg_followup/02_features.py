"""Step 2: per-residue fingerprints (repo code, unchanged defaults alpha=0.3, r=8 A) and per-individual
feature table. Output: out/fp_all_residues.csv, out/individuals.csv (one row per individual, 333 rows,
with features + 48 HPO columns coded 1/0/NaN)."""
import json

import numpy as np
import pandas as pd

from common import CACHE, COHORT_JSON, OUT, PDB, PROT_LEN, import_mvp

M = import_mvp()

fp_file = OUT / "fp_all_residues.csv"
if fp_file.exists():
    FP = pd.read_csv(fp_file)
else:
    model = M.build_model(str(PDB), None, 8.0, 0.3)
    pos, _, _ = M.load_alphafold_ca(str(PDB))
    FP = pd.DataFrame([M.fingerprint_for_variant(model, int(p)) for p in pos])
    FP["aa_pos"] = FP.aa_pos.astype(int)
    FP.to_csv(fp_file, index=False)
FP = FP.set_index("aa_pos")
assert len(FP) == PROT_LEN

# confident segments (pLDDT >= 70 for >= 20 residues), as in the draft
pl = FP.seed_pLDDT.values
segs, s = [], None
for p, c in zip(FP.index, pl):
    if c >= 70 and s is None:
        s = p
    if c < 70 and s is not None:
        if p - s >= 20:
            segs.append((s, p - 1))
        s = None
if s is not None and FP.index[-1] - s >= 19:
    segs.append((s, FP.index[-1]))
CT_SEGS = [x for x in segs if x[0] >= 2400]
print("confident segments:", segs, "C-terminal:", CT_SEGS)

# UniProt regions
up = json.loads((CACHE / "uniprot_Q6UB99.json").read_text())
DEG = [(f["location"]["start"]["value"], f["location"]["end"]["value"]) for f in up["features"]
       if f["type"] == "Region" and "degradation" in (f.get("description") or "")][0]
ANK = [(f["location"]["start"]["value"], f["location"]["end"]["value"]) for f in up["features"] if f["type"] == "Repeat"]
print("UniProt degradation region:", DEG, "ANK repeats:", ANK)

# AlphaMissense
am = pd.read_csv(CACHE / "alphamissense_Q6UB99.tsv", sep="\t", comment="#")
am_map = dict(zip(am.protein_variant, am.am_pathogenicity))
am_cls = dict(zip(am.protein_variant, am.am_class))
am_pos_mean = am.assign(pos=am.protein_variant.str[1:-1].astype(int)).groupby("pos").am_pathogenicity.mean()

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}

V = pd.read_csv(OUT / "variants_annotated.csv")
F = V[["individualId", "sex", "alleleKey", "vclass", "hgvs_c", "hgvs_p", "aa_pos", "rel_pos", "ptc_aa",
       "ptc_exon", "nmd_escape_50nt", "nmd_escape_55nt", "ptc_in_long_exon_407", "retained_aa",
       "trunc_first_aa", "vep_exon", "vep_intron", "regex_trunc", "regex_aa_pos", "sv_label", "sv_type",
       "sv_cn", "sv_size_kb", "sv_gene_overlap_hg19", "sv_gene_overlap_hg38"]].copy()
F["is_sv"] = F.vclass.eq("SV").astype(int)
F["is_point"] = 1 - F.is_sv
F["is_trunc"] = F.vclass.isin(["frameshift", "nonsense", "splice"]).astype(int)
F["is_ptc"] = F.vclass.isin(["frameshift", "nonsense"]).astype(int)
F["is_nontrunc_point"] = F.vclass.isin(["missense", "inframe"]).astype(int)
F["class3"] = np.select([F.is_sv == 1, F.is_trunc == 1], ["SV", "truncating"], "missense/inframe")
F["exon9"] = np.where(F.is_point == 1, F.vep_exon.fillna("").str.startswith("9/").astype(int), np.nan)
F["hotspot_1903_1907del"] = np.where(F.is_point == 1, F.hgvs_c.fillna("").str.endswith("c.1903_1907del").astype(int), np.nan)

ap = F.aa_pos
for col in ["seed_pLDDT", "CoreImpact", "SemiCoreImpact", "DisorderImpact", "CommBackboneImpact"]:
    F[col] = ap.map(lambda p: FP.loc[int(p), col] if pd.notna(p) else np.nan)
F["log10_semi"] = np.log10(F.SemiCoreImpact + 1e-15)
F["log10_core"] = np.log10(F.CoreImpact + 1e-15)

# truncation-specific: protein retained, C-terminal regions lost
F["ptc_rel"] = F.ptc_aa / PROT_LEN
F["loses_deg_region"] = np.where(F.is_ptc == 1, (F.ptc_aa <= DEG[0]).astype(float), np.nan)        # PTC before 2369
F["loses_all_CT_conf"] = np.where(F.is_ptc == 1, (F.ptc_aa <= CT_SEGS[0][0]).astype(float), np.nan)  # before 2436
F["retains_all_CT_conf"] = np.where(F.is_ptc == 1, (F.ptc_aa > CT_SEGS[-1][1]).astype(float), np.nan)  # after 2631


def am_lookup(r):
    if r.vclass != "missense" or not isinstance(r.hgvs_p, str):
        return pd.Series(dict(am_score=np.nan, am_class=None))
    s = r.hgvs_p.split("p.(")[-1].rstrip(")")
    key = AA3[s[:3]] + s[3:-3] + AA3[s[-3:]]
    return pd.Series(dict(am_score=am_map.get(key, np.nan), am_class=am_cls.get(key)))


F = F.join(F.apply(am_lookup, axis=1))
F["am_pos_mean"] = ap.map(lambda p: am_pos_mean.get(int(p), np.nan) if pd.notna(p) else np.nan)

# phenotypes
cohort = json.loads(COHORT_JSON.read_text())
H = [h["hpoId"] for h in cohort["hpoHeaders"]]
labels = {h["hpoId"]: h["hpoLabel"] for h in cohort["hpoHeaders"]}
ph = pd.DataFrame([{**{"individualId": r["individualData"]["individualId"]},
                    **{h: {"Observed": 1.0, "Excluded": 0.0}.get(c["type"], np.nan) for h, c in zip(H, r["hpoData"])}}
                   for r in cohort["rows"]])
F = F.merge(ph, on="individualId", how="left", validate="1:1")
F["n_annotated"] = F[H].notna().sum(axis=1)
F.to_csv(OUT / "individuals.csv", index=False)
json.dump(dict(labels=labels, hpo=H, ct_segments=CT_SEGS, segments=segs, deg_region=DEG, ank=ANK),
          open(OUT / "meta.json", "w"), indent=1)

pt = F[F.is_point == 1]
print(F.class3.value_counts())
print("PTC:", int(F.is_ptc.sum()), "loses deg region:", F.loses_deg_region.sum(), "loses all CT conf:",
      F.loses_all_CT_conf.sum(), "retains all CT conf:", F.retains_all_CT_conf.sum())
print("seed pLDDT<50:", (pt.seed_pLDDT < 50).mean(), ">=70:", (pt.seed_pLDDT >= 70).sum())
print(F.loc[F.vclass == "missense", ["individualId", "hgvs_p", "seed_pLDDT", "SemiCoreImpact", "am_score", "am_class"]])
print("corr(log10 semi, rel_pos) point variants:", pt[["log10_semi", "rel_pos"]].corr(method="spearman").iloc[0, 1])
