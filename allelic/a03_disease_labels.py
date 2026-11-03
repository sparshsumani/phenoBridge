#!/usr/bin/env python3
"""Step 3: ClinVar condition parsing and harmonisation to disease units.

ClinVar variant_summary layout: PhenotypeIDS / PhenotypeList are aligned; '|' separates the
condition sets of different RCVs, ';' separates conditions inside one multi-condition set, ','
separates the identifiers of one condition (MONDO:MONDO:x, OMIM:nnnnnn, OMIM:PSnnnnnn, MedGen:...,
Orphanet:..., Human Phenotype Ontology:HP:..., MeSH:...).

Per condition (one name + its IDs):
  * OMIM IDs are typed with OMIM mim2gene.txt: 'phenotype' / 'predominantly phenotypes' MIMs are
    phenotype MIMs; 'gene' MIMs and 'moved/removed' are not disease units; OMIM:PS* = phenotypic series.
  * MONDO term: from the ClinVar ID; else via MONDO xref (any source) to the condition's phenotype MIM.
  * MONDO -> OMIM: MONDO xrefs OMIM:nnnnnn with source "MONDO:equivalentTo" (exact match only).
Condition category (first that applies):
  nonspecific_label  name in {not provided, not specified, see cases, inborn genetic diseases,
                     hereditary disease, ...} or MONDO in a blacklist of top-level terms
  generic_gene_related  name matches '<GENE>-related <disorder|disease|condition>(s)' (or MONDO label does)
                     and no phenotype MIM is given on that condition
  cancer_no_omim     cancer/somatic-type label (regex on name, or MONDO descendant of
                     MONDO:0045024 'cancer or benign tumor' / MONDO:0015356 'hereditary neoplastic syndrome')
                     without a phenotype MIM -> dropped
  no_omim_no_mondo   none of OMIM phenotype MIM / MONDO -> dropped (e.g. MedGen-only, HPO-only)
  grouping_no_omim   no phenotype MIM (directly or as the single MONDO exact OMIM xref) and the MONDO term
                     is a grouping class: in subset 'disease_grouping'/'rare_grouping', or carries an OMIM
                     phenotypic-series xref, or has >1 exact OMIM xref, or has >=1 is_a descendant with an
                     exact OMIM xref (e.g. 'Long QT syndrome', 'Intellectual disability', 'Gaucher disease')
                     -> not a unit (recorded)
  cancer_somatic_type cancer/somatic-type LABEL (CANCER_RE on the ClinVar name) that does carry a phenotype MIM but whose name contains
                     'somatic' or does not look like a hereditary predisposition (regex PREDISP_RE:
                     syndrome|hereditary|familial|susceptib|predispos|named syndromes|'type N'|trailing number)
                     e.g. 'Gastric cancer', 'Acute myeloid leukemia', 'Testicular tumor, somatic' -> dropped
  specific           kept: unit = 'OMIM:<mim>' (direct phenotype MIM; if the condition carries several,
                     each is a unit and the condition is flagged multi_omim), else MONDO exact OMIM
                     equivalent when exactly one, else 'MONDO:<id>'.  cancer-type conditions that have a
                     phenotype MIM are kept and flagged cancer_flag.
Per variant: unit set = union over its conditions of specific units, after ancestor pruning (a
MONDO-based unit whose MONDO term is an is_a ancestor of another unit of the same variant is removed and
recorded as pruned_ancestor_units).  single_label = exactly one unit.

Outputs: intermediate/clinvar_conditions.parquet (one row per variant x condition),
         intermediate/clinvar_variant_units.parquet (one row per VariationID),
         intermediate/disease_units_clinvar.parquet, intermediate/mondo_index.parquet,
         intermediate/disease_label_log.json
"""
import re, sys, json
from collections import Counter, defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd
from common import FILES, INTER, save_json

NONSPEC = {"not provided", "not specified", "see cases", "inborn genetic diseases", "hereditary disease",
           "genetic disease", "disease", "inherited disease", "rare genetic disease", "multiple conditions",
           "variant of unknown significance", "none provided", "other", "unknown"}
NONSPEC_MONDO = {"MONDO:0003847", "MONDO:0000001", "MONDO:0700096", "MONDO:0021199"}
GENERIC_RE = re.compile(r"^\s*\S+-related\s+(disorder|disease|condition|disorders|diseases|conditions)\s*$", re.I)
CANCER_RE = re.compile(r"(carcinoma|cancer|tumou?r|neoplasm|neoplas|lymphoma|leuka?emia|sarcoma|melanoma|(?<!pseudo)glioma|"
                       r"blastoma|adenoma|myeloma|malignant (?:neoplasm|tumou?r)|mesothelioma|myelodysplas|paraganglioma|pheochromocytoma|"
                       r"somatic)", re.I)
PREDISP_RE = re.compile(r"(syndrome|hereditary|familial|susceptib|predispos|neurofibromatosis|lynch|li-fraumeni|"
                        r"tuberous sclerosis|multiple endocrine|polyposis|retinoblastoma|von hippel|cowden|carney|"
                        r"wilms|nevus|\btype\s*\d|\d+[a-z]?$)", re.I)
CANCER_ROOTS = {"MONDO:0045024", "MONDO:0015356", "MONDO:0005070", "MONDO:0004992"}


def parse_mondo(path):
    terms, cur = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = dict(parents=[], omim_exact=[], omim_any=[], omimps=[], orpha=[], subsets=[], obsolete=False, name=None)
                continue
            if line.startswith("[") and line.endswith("]"):
                cur = None; continue
            if cur is None or not line:
                continue
            k, _, v = line.partition(": ")
            if k == "id":
                if v.startswith("MONDO:"):
                    terms[v] = cur
                cur["id"] = v
            elif k == "name":
                cur["name"] = v
            elif k == "is_a":
                cur["parents"].append(v.split()[0])
            elif k == "subset":
                cur["subsets"].append(v.split()[0])
            elif k == "is_obsolete" and v.startswith("true"):
                cur["obsolete"] = True
            elif k == "xref":
                x = v.split()[0]
                exact = "MONDO:equivalentTo" in v
                if re.fullmatch(r"OMIM:\d{6}", x):
                    cur["omim_any"].append(x)
                    if exact:
                        cur["omim_exact"].append(x)
                elif x.startswith("OMIMPS:"):
                    cur["omimps"].append(x)
                elif x.startswith("Orphanet:") and exact:
                    cur["orpha"].append(x.replace("Orphanet:", "ORPHA:"))
    return terms


def ancestors_fn(terms):
    memo = {}
    def anc(t):
        if t in memo:
            return memo[t]
        memo[t] = set()
        out = set()
        for p in terms.get(t, {}).get("parents", []):
            out.add(p); out |= anc(p)
        memo[t] = out
        return out
    return anc


def split_conditions(ids, names):
    """-> list of (rcv_set_index, name, [ids]).  ids/names aligned by '|' then ';'."""
    out, bad = [], 0
    ids_sets = ids.split("|") if isinstance(ids, str) else []
    name_sets = names.split("|") if isinstance(names, str) else []
    if len(ids_sets) != len(name_sets):
        bad += 1
    for k, nm_set in enumerate(name_sets):
        id_set = ids_sets[k] if k < len(ids_sets) else ""
        nms, idl = nm_set.split(";"), id_set.split(";")
        if len(nms) != len(idl):
            # names can contain ';' rarely; fall back: all IDs to the joined name
            nms, idl = [nm_set], [",".join(idl)]
            bad += 1
        for n, i in zip(nms, idl):
            out.append((k, n.strip(), [x for x in i.split(",") if x]))
    return out, bad


def main():
    cv = pd.read_parquet(INTER / "clinvar_missense.parquet")
    terms = parse_mondo(FILES["mondo"])
    anc = ancestors_fn(terms)
    m2g = pd.read_csv(FILES["mim2gene"], sep="\t", comment="#", header=None,
                      names=["mim", "mim_type", "entrez", "symbol", "ensembl"], dtype=str)
    mimtype = dict(zip("OMIM:" + m2g.mim, m2g.mim_type))
    hpoa = pd.read_csv(FILES["hpoa"], sep="\t", comment="#", dtype=str)
    hpoa_names = hpoa.drop_duplicates("database_id").set_index("database_id").disease_name.to_dict()
    omim2mondo = defaultdict(set)
    for t, d in terms.items():
        if d["obsolete"]:
            continue
        for o in d["omim_exact"]:
            omim2mondo[o].add(t)
    omim2mondo_any = defaultdict(set)
    for t, d in terms.items():
        if not d["obsolete"]:
            for o in d["omim_any"]:
                omim2mondo_any[o].add(t)
    children = defaultdict(set)
    for t, d in terms.items():
        for p in d["parents"]:
            children[p].add(t)
    desc_memo = {}
    def has_omim_desc(t):
        if t in desc_memo:
            return desc_memo[t]
        desc_memo[t] = False
        r = any((terms.get(c, {}).get("omim_exact") and not terms[c]["obsolete"]) or has_omim_desc(c) for c in children.get(t, ()))
        desc_memo[t] = r
        return r
    is_cancer_mondo = lambda t: bool(t) and (t in CANCER_ROOTS or bool(anc(t) & CANCER_ROOTS))
    log = Counter()
    crow = []
    for vid, ids, names in zip(cv.variation_id, cv.phenotype_ids, cv.phenotype_list):
        conds, bad = split_conditions(ids, names)
        log["variants_with_misaligned_condition_fields"] += bool(bad)
        seen = set()
        for k, nm, idl in conds:
            key = (nm.lower(), tuple(sorted(idl)))
            if key in seen:
                continue
            seen.add(key)
            omims = [x for x in idl if re.fullmatch(r"OMIM:\d{6}", x)]
            omimps = [x for x in idl if x.startswith("OMIM:PS")]
            mondo = [x.replace("MONDO:MONDO:", "MONDO:") for x in idl if x.startswith("MONDO:")]
            pheno_mims = [o for o in omims if mimtype.get(o) in ("phenotype", "predominantly phenotypes")]
            other_mims = [o for o in omims if o not in pheno_mims]
            mondo_id = mondo[0] if mondo else None
            if mondo_id and terms.get(mondo_id, {}).get("obsolete"):
                log["obsolete_mondo_in_clinvar"] += 1
            if mondo_id is None and len(pheno_mims) == 1:
                cand = omim2mondo.get(pheno_mims[0]) or omim2mondo_any.get(pheno_mims[0]) or set()
                mondo_id = sorted(cand)[0] if len(cand) == 1 else None
            mt = terms.get(mondo_id, {}) if mondo_id else {}
            mondo_name = mt.get("name")
            mondo_omim = [o for o in mt.get("omim_exact", []) if mimtype.get(o) in ("phenotype", "predominantly phenotypes")]
            grouping = bool(mt) and ("disease_grouping" in mt.get("subsets", []) or "rare_grouping" in mt.get("subsets", []))
            cancer = bool(CANCER_RE.search(nm)) or bool(mondo_name and CANCER_RE.search(mondo_name)) or is_cancer_mondo(mondo_id)
            generic = bool(GENERIC_RE.match(nm)) or bool(mondo_name and GENERIC_RE.match(mondo_name))
            units, cat, multi = [], None, False
            if nm.lower() in NONSPEC or (mondo_id in NONSPEC_MONDO and not pheno_mims):
                cat = "nonspecific_label"
            elif generic and not pheno_mims:
                cat = "generic_gene_related"
            elif cancer and not pheno_mims:
                cat = "cancer_no_omim"
            elif CANCER_RE.search(nm) and (re.search("somatic", nm, re.I) or not PREDISP_RE.search(nm)):
                cat = "cancer_somatic_type"   # cancer label with a MIM but not a hereditary predisposition label
            elif pheno_mims:
                cat = "specific"; units = sorted(set(pheno_mims)); multi = len(units) > 1
            elif mondo_id and len(mondo_omim) == 1:
                cat = "specific"; units = mondo_omim
            elif mondo_id and (grouping or omimps or mt.get("omimps") or len(mondo_omim) > 1 or has_omim_desc(mondo_id)):
                cat = "grouping_no_omim"
            elif mondo_id:
                cat = "specific"; units = [mondo_id]
            elif other_mims:
                cat = "omim_gene_or_removed_only"
            else:
                cat = "no_omim_no_mondo"
            log[f"conditions_{cat}"] += 1
            crow.append(dict(variation_id=vid, rcv_set=k, condition_name=nm, condition_ids=",".join(idl),
                             omim_pheno=";".join(pheno_mims) or None, omim_other=";".join(other_mims) or None,
                             omim_ps=";".join(omimps) or None, mondo=mondo_id, mondo_name=mondo_name,
                             mondo_grouping=grouping, cancer_flag=cancer, generic_flag=generic, category=cat,
                             units=";".join(units) or None, multi_omim=multi))
    cd = pd.DataFrame(crow)
    cd.to_parquet(INTER / "clinvar_conditions.parquet", index=False)

    # unit -> MONDO (for ancestor pruning and annotation)
    unit_mondo = {}
    for u, m in zip(cd.units, cd.mondo):
        if isinstance(u, str):
            for x in u.split(";"):
                if x.startswith("MONDO:"):
                    unit_mondo[x] = x
                elif m and x not in unit_mondo and len(u.split(";")) == 1:
                    unit_mondo[x] = m
    for u in list(unit_mondo):
        pass
    for o, ms in omim2mondo.items():
        if len(ms) == 1:
            unit_mondo.setdefault(o, next(iter(ms)))

    vrows = []
    for vid, g in cd.groupby("variation_id", sort=False):
        units = sorted({x for u in g.units.dropna() for x in u.split(";")})
        pruned = []
        for u in units:
            mu = unit_mondo.get(u)
            if not mu:
                continue
            for w in units:
                if w != u and unit_mondo.get(w) and mu in anc(unit_mondo[w]):
                    pruned.append(u); break
        units = [u for u in units if u not in pruned]
        cats = Counter(g.category)
        spec = g[g.category == "specific"]
        vrows.append(dict(variation_id=vid, disease_units=units, n_disease_units=len(units),
                          single_label=len(units) == 1, pruned_ancestor_units=pruned or None,
                          n_conditions=len(g), n_conditions_specific=cats["specific"],
                          has_generic_condition=bool(cats["generic_gene_related"]),
                          has_nonspecific_condition=bool(cats["nonspecific_label"]),
                          has_grouping_condition=bool(cats["grouping_no_omim"]),
                          dropped_cancer_condition=bool(cats["cancer_no_omim"]),
                          cancer_flag=bool(spec.cancer_flag.any()), multi_omim_condition=bool(spec.multi_omim.any()),
                          clinvar_condition_names_kept=" | ".join(dict.fromkeys(spec.condition_name)),
                          clinvar_condition_names_all=" | ".join(dict.fromkeys(g.condition_name))))
    vu = pd.DataFrame(vrows)
    vu.to_parquet(INTER / "clinvar_variant_units.parquet", index=False)

    # unit table
    urows = []
    allu = sorted({u for us in vu.disease_units for u in us})
    names_by_unit = defaultdict(Counter)
    for u, nm in zip(cd.units, cd.condition_name):
        if isinstance(u, str) and len(u.split(";")) == 1:
            names_by_unit[u][nm] += 1
    for u in allu:
        m = unit_mondo.get(u)
        mt = terms.get(m, {})
        urows.append(dict(disease_unit=u, unit_type=u.split(":")[0],
                          disease_name=hpoa_names.get(u) or (names_by_unit[u].most_common(1)[0][0] if names_by_unit[u] else mt.get("name")),
                          clinvar_name=(names_by_unit[u].most_common(1)[0][0] if names_by_unit[u] else None),
                          mondo=m, mondo_name=mt.get("name"),
                          mondo_grouping=("disease_grouping" in mt.get("subsets", [])) if mt else None,
                          mondo_orpha=";".join(mt.get("orpha", [])) or None,
                          omim_mim_type=mimtype.get(u) if u.startswith("OMIM:") else None,
                          cancer_flag=bool(mt and is_cancer_mondo(m)) or bool(CANCER_RE.search(str(names_by_unit[u].most_common(1)[0][0] if names_by_unit[u] else ""))),
                          in_hpoa=u in hpoa_names))
    du = pd.DataFrame(urows)
    du.to_parquet(INTER / "disease_units_clinvar.parquet", index=False)
    mi = pd.DataFrame([dict(mondo=t, name=d["name"], obsolete=d["obsolete"], omim_exact=";".join(d["omim_exact"]) or None,
                            orpha=";".join(d["orpha"]) or None, grouping="disease_grouping" in d["subsets"],
                            parents=";".join(d["parents"]) or None) for t, d in terms.items()])
    mi.to_parquet(INTER / "mondo_index.parquet", index=False)
    log = dict(log)
    log.update(variants=len(vu), variants_with_ge1_unit=int((vu.n_disease_units >= 1).sum()),
               variants_no_unit=int((vu.n_disease_units == 0).sum()),
               variants_single_label=int(vu.single_label.sum()), variants_multi_label=int((vu.n_disease_units > 1).sum()),
               variants_with_pruned_ancestor=int(vu.pruned_ancestor_units.notna().sum()),
               variants_generic_condition=int(vu.has_generic_condition.sum()),
               variants_cancer_flag=int(vu.cancer_flag.sum()),
               disease_units=len(du), disease_units_omim=int((du.unit_type == "OMIM").sum()),
               disease_units_mondo_only=int((du.unit_type == "MONDO").sum()),
               n_units_distribution={int(k): int(v) for k, v in vu.n_disease_units.value_counts().sort_index().items()},
               mondo_terms=len(terms))
    save_json(log, INTER / "disease_label_log.json")
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
