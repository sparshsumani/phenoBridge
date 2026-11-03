#!/usr/bin/env python3
"""Step 1: parse all phenopacket-store phenopackets into tidy tables.

Outputs (INTER/):
  individuals_raw.parquet   one row per phenopacket
  variants_raw.parquet      one row per phenopacket x genomicInterpretation (distinct causative variant)
  phenotypes_explicit.parquet  phenopacket x HPO term, value 1 (observed) / 0 (excluded); never imputed
  parse_log.json            counts
"""
import glob, json, re
from collections import Counter
import pandas as pd
from common import PPKT_ROOT, INTER

ISO = re.compile(r"^P(?:(\d+(?:\.\d+)?)Y)?(?:(\d+(?:\.\d+)?)M)?(?:(\d+(?:\.\d+)?)W)?(?:(\d+(?:\.\d+)?)D)?")


def iso_years(s):
    m = ISO.match(s or "")
    if not m or not any(m.groups()):
        return None
    y, mo, w, d = (float(x) if x else 0.0 for x in m.groups())
    return y + mo / 12.0 + w * 7 / 365.25 + d / 365.25


def norm_chrom(c):
    c = str(c)
    if c.startswith("chr"):
        c = c[3:]
    if c in ("M", "MT"):
        return "MT"
    return c


def main():
    files = sorted(glob.glob(str(PPKT_ROOT / "*/phenopackets/*.json")))
    log = Counter()
    ind_rows, var_rows, ph_rows = [], [], []
    conflicts = []
    for f in files:
        d = json.load(open(f))
        log["phenopackets_read"] += 1
        cohort = f.split("/")[-3]
        pid = d["id"]
        subj = d.get("subject", {})
        t = subj.get("timeAtLastEncounter") or {}
        age_raw, age_years = None, None
        if "age" in t:
            age_raw = t["age"].get("iso8601duration")
            age_years = iso_years(age_raw)
        elif "gestationalAge" in t:
            age_raw = "gestational:{}w{}d".format(t["gestationalAge"].get("weeks"), t["gestationalAge"].get("days", 0))
        elif "ontologyClass" in t:
            age_raw = t["ontologyClass"]["id"] + "|" + t["ontologyClass"].get("label", "")
        refs = d.get("metaData", {}).get("externalReferences", [])
        pmids = [r["id"] for r in refs if r.get("id", "").startswith("PMID")]
        dis = d.get("diseases", [])
        interps = d.get("interpretations", [])
        gis = []
        for it in interps:
            gis.extend(it.get("diagnosis", {}).get("genomicInterpretations", []))
        genes, states = [], []
        n_alleles_total = 0
        n_sv = 0
        for k, gi in enumerate(gis):
            vi = gi.get("variantInterpretation", {})
            vd = vi.get("variationDescriptor", {})
            gc = vd.get("geneContext", {})
            ex = {e["syntax"]: e["value"] for e in vd.get("expressions", [])}
            vcf = vd.get("vcfRecord")
            st = vd.get("allelicState", {}).get("label")
            n_all = 2 if st == "homozygous" else 1
            n_alleles_total += n_all
            genes.append(gc.get("symbol"))
            states.append(st)
            is_sv = vcf is None
            n_sv += int(is_sv)
            hc = ex.get("hgvs.c")
            hp = ex.get("hgvs.p")
            var_rows.append(dict(
                phenopacket_id=pid, cohort=cohort, gi_index=k,
                variant_id=vd.get("id"),
                gene=gc.get("symbol"), hgnc_id=gc.get("valueId"),
                interpretation_status=gi.get("interpretationStatus"),
                acmg=vi.get("acmgPathogenicityClassification"),
                allelic_state=st, n_alleles=n_all,
                hgvs_c=hc, transcript_c=hc.split(":")[0] if hc else None,
                hgvs_g=ex.get("hgvs.g"), hgvs_p=hp, protein_acc_p=hp.split(":")[0] if hp else None,
                chrom=norm_chrom(vcf["chrom"]) if vcf else None,
                pos=int(vcf["pos"]) if vcf else None,
                ref=vcf.get("ref") if vcf else None, alt=vcf.get("alt") if vcf else None,
                assembly=vcf.get("genomeAssembly") if vcf else None,
                is_structural=is_sv,
                structural_type=vd.get("structuralType", {}).get("label"),
                sv_label=vd.get("label"),
            ))
        gene_set = sorted({g for g in genes if g})
        if len(gis) == 0:
            zyg = "none"
        elif len(gis) == 1:
            zyg = {"heterozygous": "het", "homozygous": "hom", "hemizygous": "hemi"}.get(states[0], states[0])
        else:
            zyg = "comphet" if all(s == "heterozygous" for s in states) else "multi_" + "+".join(sorted(states))
        ind_rows.append(dict(
            phenopacket_id=pid, cohort=cohort,
            gene=gene_set[0] if len(gene_set) == 1 else ("|".join(gene_set) if gene_set else None),
            n_genes=len(gene_set),
            subject_id=subj.get("id"), sex=subj.get("sex"),
            age_raw=age_raw, age_years=age_years,
            vital_status=(subj.get("vitalStatus") or {}).get("status"),
            pmid=pmids[0] if pmids else None,
            pmid_title=refs[0].get("description") if refs else None,
            disease_id=dis[0]["term"]["id"] if dis else None,
            disease_label=dis[0]["term"].get("label") if dis else None,
            n_diseases=len(dis),
            n_genomic_interpretations=len(gis),
            n_causative_alleles=n_alleles_total,
            allelic_states="+".join(s or "NA" for s in states),
            zygosity=zyg,
            biallelic=(zyg in ("hom", "comphet")) or (len(gis) >= 2),
            n_structural_variants=n_sv,
            file=f.split("/notebooks/")[-1],
        ))
        # phenotypes: explicit only
        pfs = d.get("phenotypicFeatures", [])
        if not pfs:
            log["phenopackets_without_phenotypicFeatures"] += 1
        per = {}
        labels = {}
        for pf in pfs:
            hid = pf["type"]["id"]
            labels[hid] = pf["type"].get("label")
            val = 0 if pf.get("excluded", False) else 1
            per.setdefault(hid, set()).add(val)
            log["pf_entries"] += 1
        for hid, vals in per.items():
            if len(vals) > 1:
                conflicts.append((pid, hid))
                log["explicit_conflict_cells_dropped"] += 1
                continue
            ph_rows.append(dict(phenopacket_id=pid, cohort=cohort, hpo_id=hid, hpo_label=labels[hid], value=vals.pop()))

    ind = pd.DataFrame(ind_rows)
    var = pd.DataFrame(var_rows)
    ph = pd.DataFrame(ph_rows)
    log["individuals"] = len(ind)
    log["genomic_interpretations"] = len(var)
    log["gi_structural_no_vcf"] = int(var.is_structural.sum())
    log["gi_without_gene"] = int(var.gene.isna().sum())
    log["explicit_cells"] = len(ph)
    log["explicit_observed"] = int((ph.value == 1).sum())
    log["explicit_excluded"] = int((ph.value == 0).sum())
    log["duplicate_pf_entries_collapsed"] = log["pf_entries"] - len(ph) - log["explicit_conflict_cells_dropped"]
    for c in ("pos",):
        var[c] = var[c].astype("Int64")
    ind.to_parquet(INTER / "individuals_raw.parquet", index=False)
    var.to_parquet(INTER / "variants_raw.parquet", index=False)
    ph.to_parquet(INTER / "phenotypes_explicit.parquet", index=False)
    (INTER / "parse_log.json").write_text(json.dumps(dict(log, conflict_examples=conflicts[:20]), indent=1))
    print(json.dumps(dict(log), indent=1))
    print(ind.zygosity.value_counts())


if __name__ == "__main__":
    main()
