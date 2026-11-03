#!/usr/bin/env python3
"""Step 2: filter ClinVar variant_summary.txt.gz to germline P/LP, >=1-star, GRCh38 missense SNVs.

Filters (applied in this order; rows and distinct VariationIDs dropped at each are logged):
  F1 Assembly == GRCh38
  F2 Type == "single nucleotide variant"
  F3 ClinicalSignificance (= germline classification) exactly one of
       Pathogenic | Likely pathogenic | Pathogenic/Likely pathogenic
     (excludes conflicting, VUS, benign, risk factor, drug response, and P/LP combined with any
      extra term such as "Pathogenic; risk factor"; somatic clinical impact / oncogenicity are
      separate columns and never enter)
  F4 OriginSimple != "somatic" (pure somatic-origin records excluded; germline/somatic kept)
  F5 review status >= 1 star (drops "no assertion criteria provided", "no classification provided",
     "no classification for the individual variant", "no classifications from unflagged records")
  F6 Name parses as <RefSeq NM_/XM_>(<GENE>):c.<...> (p.<...>)  (gene symbol and p. change present)
  F7 protein change is missense: p.Xaa<n>Yaa with Xaa != Yaa, neither = Ter, and not Met1 (start-loss)
     (synonymous p.=/p.Xaa<n>=, nonsense, stop-loss, start-loss, p.? are dropped and counted)
Output: derived_allelic/intermediate/clinvar_missense.parquet (one row per ClinVar VariationID)
        derived_allelic/intermediate/clinvar_filter_log.json
"""
import re, sys, json, gzip
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd
from common import FILES, INTER, AA3, save_json

STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no assertion provided": 0,
    "no classification for the individual variant": 0,
    "no classifications from unflagged records": 0,
    "no interpretation for the single variant": 0,
}
PLP = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
NAME_RE = re.compile(r"^(?P<tx>[NX]M_\d+(?:\.\d+)?)\((?P<gene>[^)]+)\):(?P<c>c\.\S+)(?: \((?P<p>p\.[^)]+)\))?\s*$")
P_RE = re.compile(r"^p\.(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2}|=|\*|\?)$")
KEEP = ["#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "HGNC_ID", "ClinicalSignificance", "LastEvaluated",
        "RS# (dbSNP)", "RCVaccession", "PhenotypeIDS", "PhenotypeList", "Origin", "OriginSimple", "Assembly",
        "Chromosome", "Start", "ReviewStatus", "NumberSubmitters", "SubmitterCategories", "VariationID",
        "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]


def classify_p(p):
    if p is None or (isinstance(p, float)):
        return "no_protein_change", None
    m = P_RE.match(p)
    if not m:
        return "protein_change_unparsed", None
    ref, pos, alt = m.group("ref"), int(m.group("pos")), m.group("alt")
    if alt == "=" or alt == ref:
        return "synonymous", None
    if alt == "?":
        return "p_unknown", None
    if ref == "Ter":
        return "stop_loss", None
    if alt in ("Ter", "*"):
        return "nonsense", None
    if pos == 1 and ref == "Met":
        return "start_loss", None
    if ref not in AA3 or alt not in AA3:
        return "nonstandard_aa", None
    return "missense", (AA3[ref], pos, AA3[alt])


def main():
    src = FILES["clinvar"]
    with gzip.open(src, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    usecols = [c for c in KEEP if c in header]
    missing = [c for c in KEEP if c not in header]
    log = {"source_file": str(src), "header_columns": header, "missing_expected_columns": missing}
    steps = ["F0_all", "F1_GRCh38", "F2_SNV", "F3_PLP_exact", "F4_not_somatic_only", "F5_stars_ge1",
             "F6_name_parsed", "F7_missense"]
    rows = Counter(); vids = {s: set() for s in steps}
    sig_counts, rs_counts, pclass = Counter(), Counter(), Counter()
    sig_excluded_plp_like = Counter()
    keep = []
    for ch in pd.read_csv(src, sep="\t", dtype=str, usecols=usecols, chunksize=500_000, compression="gzip",
                          quoting=3, keep_default_na=False, na_values=[]):
        def mark(name, d):
            rows[name] += len(d); vids[name].update(d.VariationID)
        mark("F0_all", ch)
        d = ch[ch.Assembly == "GRCh38"]; mark("F1_GRCh38", d)
        d = d[d.Type == "single nucleotide variant"]; mark("F2_SNV", d)
        sig_counts.update(d.ClinicalSignificance)
        x = d[~d.ClinicalSignificance.isin(PLP) & d.ClinicalSignificance.str.contains("athogenic", na=False)]
        sig_excluded_plp_like.update(x.ClinicalSignificance)
        d = d[d.ClinicalSignificance.isin(PLP)]; mark("F3_PLP_exact", d)
        d = d[d.OriginSimple.str.strip().str.lower() != "somatic"]; mark("F4_not_somatic_only", d)
        rs_counts.update(d.ReviewStatus)
        d = d.assign(stars=d.ReviewStatus.map(STARS))
        unk = d.stars.isna()
        if unk.any():
            log.setdefault("unknown_review_status", []).extend(sorted(set(d.ReviewStatus[unk])))
        d = d[d.stars.fillna(0) >= 1]; mark("F5_stars_ge1", d)
        m = d.Name.str.extract(NAME_RE)
        d = pd.concat([d, m.rename(columns={"tx": "name_tx", "gene": "name_gene", "c": "hgvs_c", "p": "hgvs_p"})], axis=1)
        ok = d.name_tx.notna() & d.name_gene.notna() & d.hgvs_p.notna()
        pclass["name_unparsed_or_no_p"] += int((~ok).sum())
        d = d[ok]; mark("F6_name_parsed", d)
        cl = [classify_p(p) for p in d.hgvs_p]
        pclass.update(c for c, _ in cl)
        d = d.assign(p_class=[c for c, _ in cl], ref_aa=[t[0] if t else None for _, t in cl],
                     clinvar_protein_pos=[t[1] if t else None for _, t in cl], alt_aa=[t[2] if t else None for _, t in cl])
        d = d[d.p_class == "missense"]; mark("F7_missense", d)
        keep.append(d)
        print(dict(rows), flush=True)
    cv = pd.concat(keep, ignore_index=True)
    dup = cv.VariationID.duplicated(keep=False)
    log["duplicate_variationid_rows"] = int(dup.sum())
    cv = cv.drop_duplicates("VariationID")
    cv["stars"] = cv.stars.astype(int)
    cv["n_submitters"] = pd.to_numeric(cv.NumberSubmitters, errors="coerce").astype("Int64")
    cv["last_evaluated"] = cv.LastEvaluated.replace({"-": None, "": None})
    cv["name_tx_noversion"] = cv.name_tx.str.split(".").str[0]
    cv["clinvar_protein_pos"] = cv.clinvar_protein_pos.astype(int)
    cv = cv.rename(columns={"#AlleleID": "allele_id", "VariationID": "variation_id", "GeneSymbol": "clinvar_gene_symbols",
                            "HGNC_ID": "clinvar_hgnc_ids", "ClinicalSignificance": "clinvar_classification",
                            "ReviewStatus": "review_status", "RS# (dbSNP)": "rsid", "RCVaccession": "rcv",
                            "PhenotypeIDS": "phenotype_ids", "PhenotypeList": "phenotype_list", "Chromosome": "chrom",
                            "PositionVCF": "pos_vcf", "ReferenceAlleleVCF": "ref_vcf", "AlternateAlleleVCF": "alt_vcf",
                            "OriginSimple": "origin_simple", "SubmitterCategories": "submitter_categories"})
    cv = cv.drop(columns=["Type", "Assembly", "NumberSubmitters", "LastEvaluated", "p_class"])
    cv.to_parquet(INTER / "clinvar_missense.parquet", index=False)
    order = steps
    log["rows_after"] = {s: rows[s] for s in order}
    log["variation_ids_after"] = {s: len(vids[s]) for s in order}
    log["dropped_variation_ids"] = {order[i]: len(vids[order[i - 1]]) - len(vids[order[i]]) for i in range(1, len(order))}
    log["germline_classification_counts_GRCh38_SNV"] = dict(sig_counts.most_common())
    log["excluded_PLP_with_extra_terms_or_conflicting"] = dict(sig_excluded_plp_like.most_common(30))
    log["review_status_after_F4"] = dict(rs_counts.most_common())
    log["protein_change_classes_after_F5"] = dict(pclass.most_common())
    log["final_variation_ids"] = int(len(cv))
    log["final_genes"] = int(cv.name_gene.nunique())
    log["stars_distribution"] = {int(k): int(v) for k, v in cv.stars.value_counts().sort_index().items()}
    log["classification_distribution"] = dict(cv.clinvar_classification.value_counts())
    save_json(log, INTER / "clinvar_filter_log.json")
    print(json.dumps({k: v for k, v in log.items() if k != "header_columns"}, indent=1, default=str))


if __name__ == "__main__":
    main()
