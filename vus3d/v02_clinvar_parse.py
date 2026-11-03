#!/usr/bin/env python3
"""Step 2: parse both ClinVar variant_summary releases (2023-01 and 2026-11-01) with the allelic/a02
filters except the classification filter, and classify every record.

Per release:
  F0 all rows
  F1 Assembly == GRCh38                 -> intermediate/clinvar_<rel>_grch38.parquet (ALL types, slim columns;
                                           used to look up the status of any VariationID / coordinate in that release)
  F2 Type == "single nucleotide variant"
  F4 OriginSimple != "somatic"
  F6 Name parses as <NM_/XM_>(<GENE>):c.<...> (p.<...>)     (allelic/a02 NAME_RE)
  F7 missense (allelic/a02 classify_p: not synonymous / nonsense / stop-loss / start-loss Met1 / p.?)
                                        -> intermediate/clinvar_<rel>_missense.parquet (all classifications, all stars)
  (F5 review status >= 1 star and the classification category are applied in v03; star mapping = a02 STARS.)
Classification category (germline ClinicalSignificance, exact strings):
  PLP      Pathogenic | Likely pathogenic | Pathogenic/Likely pathogenic
  BLB      Benign | Likely benign | Benign/Likely benign
  VUS      Uncertain significance
  CONFLICT starts with "Conflicting" (2023: "Conflicting interpretations of pathogenicity",
           2026: "Conflicting classifications of pathogenicity"; incl. variants with an added "; <term>")
  OTHER    everything else (e.g. "Pathogenic; risk factor", "Pathogenic, low penetrance", "risk factor",
           "drug response", "not provided", "Uncertain significance; risk factor", "other", ...)
Also a02-style sequential counts restricted to labelled classes (PLP or BLB): F3 class, F4, F5, F6, F7.
Output intermediate/clinvar_parse_log.json
"""
import re, sys, json
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd
from common import RELEASES, INTER, save_json, al_import

a02 = al_import("a02_clinvar_filter")
STARS, NAME_RE, classify_p = a02.STARS, a02.NAME_RE, a02.classify_p
PLP = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BLB = {"Benign", "Likely benign", "Benign/Likely benign"}
COLS = ["#AlleleID", "Type", "Name", "GeneSymbol", "HGNC_ID", "ClinicalSignificance", "LastEvaluated", "RS# (dbSNP)",
        "OriginSimple", "Assembly", "Chromosome", "ReviewStatus", "NumberSubmitters", "VariationID", "PositionVCF",
        "ReferenceAlleleVCF", "AlternateAlleleVCF", "RCVaccession"]


def category(sig):
    s = sig.fillna("")
    out = pd.Series("OTHER", index=sig.index)
    out[s.isin(PLP)] = "PLP"
    out[s.isin(BLB)] = "BLB"
    out[s == "Uncertain significance"] = "VUS"
    out[s.str.startswith("Conflicting")] = "CONFLICT"
    return out


def parse(rel):
    src = RELEASES[rel]["file"]
    log = dict(source_file=str(src))
    rows, vids = Counter(), {}
    lab_vids = {}
    def mark(d_, name, store):
        store.setdefault(name, set()).update(d_.VariationID)
    g38, mis = [], []
    unknown_rs = Counter()
    for ch in pd.read_csv(src, sep="\t", dtype=str, usecols=COLS, chunksize=500_000, compression="gzip",
                          quoting=3, keep_default_na=False, na_values=[]):
        rows["F0_all"] += len(ch); mark(ch, "F0_all", vids)
        d = ch[ch.Assembly == "GRCh38"].copy(); rows["F1_GRCh38"] += len(d); mark(d, "F1_GRCh38", vids)
        d["category"] = category(d.ClinicalSignificance)
        d["stars"] = d.ReviewStatus.map(STARS)
        unknown_rs.update(d.ReviewStatus[d.stars.isna()])
        d["stars"] = d.stars.fillna(0).astype(int)
        g38.append(d.drop(columns=["Assembly"]))
        # a02-style sequential counts for labelled classes
        x = d[d.Type == "single nucleotide variant"]; rows["F2_SNV"] += len(x); mark(x, "F2_SNV", vids)
        for lab in ("PLP", "BLB"):
            y = x[x.category == lab]; mark(y, "F3_class", lab_vids.setdefault(lab, {}))
            y = y[y.OriginSimple.str.strip().str.lower() != "somatic"]; mark(y, "F4_not_somatic_only", lab_vids[lab])
            y = y[y.stars >= 1]; mark(y, "F5_stars_ge1", lab_vids[lab])
            m = y.Name.str.extract(NAME_RE); y = y[m.tx.notna() & m.gene.notna() & m.p.notna()]
            mark(y, "F6_name_parsed", lab_vids[lab])
            y = y[[classify_p(p)[0] == "missense" for p in m.loc[y.index, "p"]]]; mark(y, "F7_missense", lab_vids[lab])
        # missense subset, all classes / stars
        x = x[x.OriginSimple.str.strip().str.lower() != "somatic"]; rows["F4_not_somatic_only"] += len(x); mark(x, "F4_not_somatic_only", vids)
        m = x.Name.str.extract(NAME_RE)
        x = pd.concat([x, m.rename(columns={"tx": "name_tx", "gene": "name_gene", "c": "hgvs_c", "p": "hgvs_p"})], axis=1)
        x = x[x.name_tx.notna() & x.name_gene.notna() & x.hgvs_p.notna()]; rows["F6_name_parsed"] += len(x); mark(x, "F6_name_parsed", vids)
        cl = [classify_p(p) for p in x.hgvs_p]
        x = x.assign(p_class=[c for c, _ in cl], ref_aa=[t[0] if t else None for _, t in cl],
                     clinvar_protein_pos=[t[1] if t else None for _, t in cl], alt_aa=[t[2] if t else None for _, t in cl])
        x = x[x.p_class == "missense"]; rows["F7_missense"] += len(x); mark(x, "F7_missense", vids)
        mis.append(x.drop(columns=["Assembly", "p_class"]))
        print(rel, dict(rows), flush=True)
    ren = {"#AlleleID": "allele_id", "VariationID": "variation_id", "GeneSymbol": "clinvar_gene_symbols", "HGNC_ID": "clinvar_hgnc_ids",
           "ClinicalSignificance": "classification", "ReviewStatus": "review_status", "RS# (dbSNP)": "rsid", "RCVaccession": "rcv",
           "Chromosome": "chrom", "PositionVCF": "pos_vcf", "ReferenceAlleleVCF": "ref_vcf", "AlternateAlleleVCF": "alt_vcf",
           "OriginSimple": "origin_simple", "Type": "type", "Name": "name", "NumberSubmitters": "n_submitters", "LastEvaluated": "last_evaluated_raw"}
    def tidy(df):
        df = df.rename(columns=ren)
        df["n_submitters"] = pd.to_numeric(df.n_submitters, errors="coerce").astype("Int64")
        le = df.last_evaluated_raw.replace({"-": None, "": None})
        df["last_evaluated"] = pd.to_datetime(le, format="%b %d, %Y", errors="coerce")
        df["coord_key"] = df.chrom + ":" + df.pos_vcf + ":" + df.ref_vcf + ":" + df.alt_vcf
        return df
    g = tidy(pd.concat(g38, ignore_index=True))
    log["grch38_rows"] = len(g)
    log["grch38_duplicate_variationid_rows"] = int(g.variation_id.duplicated().sum())
    log["last_evaluated_unparsed_nonempty"] = int((g.last_evaluated.isna() & g.last_evaluated_raw.replace({"-": "", None: ""}).ne("")).sum())
    g = g.drop_duplicates("variation_id")
    g.drop(columns=["rcv"]).to_parquet(INTER / f"clinvar_{rel}_grch38.parquet", index=False)
    mm = tidy(pd.concat(mis, ignore_index=True))
    log["missense_duplicate_variationid_rows"] = int(mm.variation_id.duplicated().sum())
    mm = mm.drop_duplicates("variation_id")
    mm["clinvar_protein_pos"] = mm.clinvar_protein_pos.astype(int)
    mm["name_tx_noversion"] = mm.name_tx.str.split(".").str[0]
    mm.to_parquet(INTER / f"clinvar_{rel}_missense.parquet", index=False)
    steps = ["F0_all", "F1_GRCh38", "F2_SNV", "F4_not_somatic_only", "F6_name_parsed", "F7_missense"]
    log["rows_after"] = {s: rows[s] for s in steps}
    log["variation_ids_after"] = {s: len(vids[s]) for s in steps}
    lsteps = ["F3_class", "F4_not_somatic_only", "F5_stars_ge1", "F6_name_parsed", "F7_missense"]
    log["labelled_sequential_variation_ids"] = {lab: {s: len(lab_vids[lab].get(s, ())) for s in lsteps} for lab in ("PLP", "BLB")}
    log["unknown_review_status_strings"] = dict(unknown_rs)
    log["grch38_category_counts"] = {k: int(v) for k, v in g.category.value_counts().items()}
    log["missense_category_counts"] = {k: int(v) for k, v in mm.category.value_counts().items()}
    log["missense_category_by_stars"] = {f"{c}|{s}": int(n) for (c, s), n in mm.groupby(["category", "stars"]).size().items()}
    log["missense_other_classifications_top"] = {k: int(v) for k, v in mm[mm.category == "OTHER"].classification.value_counts().head(25).items()}
    log["missense_conflict_strings"] = {k: int(v) for k, v in mm[mm.category == "CONFLICT"].classification.value_counts().head(10).items()}
    log["review_status_counts_missense"] = {k: int(v) for k, v in mm.review_status.value_counts().items()}
    return log


def main():
    log = {rel: parse(rel) for rel in RELEASES}
    save_json(log, INTER / "clinvar_parse_log.json")
    print(json.dumps(log, indent=1, default=str))


if __name__ == "__main__":
    main()
