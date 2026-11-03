#!/usr/bin/env python3
"""Step 6: disease-unit table with mode of inheritance from HPO annotations (phenotype.hpoa).

Units = every ClinVar disease unit (a03) U every phenopacket-store disease id among merged variants (a05).
Inheritance terms = all aspect 'I' annotations of the unit in phenotype.hpoa (OMIM units: the OMIM
record; MONDO-only units: the Orphanet record(s) that MONDO declares exact-equivalent, if any).
Each term is assigned a major class from the HPO is_a hierarchy (hp.json, cached raw/hpo/hp.json):
  AD  = HP:0000006 or descendants (incl. AD germline de novo HP:0025352, sex-limited AD, imprinting)
  AR  = HP:0000007 or descendants
  XLR = HP:0001419, XLD = HP:0001423, XL = HP:0001417 (unspecified X-linked)
  MT  = HP:0001427 (mitochondrial), YL = HP:0001450
  SD  = HP:0032113 semidominant
  COMPLEX = HP:0001426 multifactorial and descendants (digenic, oligogenic, polygenic)
  modifiers (sporadic HP:0003745, somatic mosaicism HP:0001442, incomplete penetrance HP:0003829,
  anticipation, contiguous gene syndrome, ...) are recorded in inheritance_terms but get no class.
inheritance_clean (single label when unambiguous): X-linked classes are first collapsed
  (XL+XLR -> XLR, XL+XLD -> XLD, XLR+XLD(+XL) -> XL); then exactly one class remaining -> that class;
  otherwise None (e.g. AD+AR).  inheritance_ad_only / inheritance_ar_only are convenience flags.
broad_unit: the unit's MONDO term has >= 1 is_a descendant with an exact OMIM xref (an umbrella entry
  such as OMIM:268000 'Retinitis pigmentosa'); flagged only, still a unit.
Output intermediate/diseases.parquet, intermediate/inheritance_log.json
"""
import sys, json
from collections import Counter, defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd
from common import FILES, INTER, RAW, save_json

ROOTS = {"HP:0000006": "AD", "HP:0000007": "AR", "HP:0001419": "XLR", "HP:0001423": "XLD", "HP:0001417": "XL",
         "HP:0001427": "MT", "HP:0001450": "YL", "HP:0032113": "SD", "HP:0001426": "COMPLEX"}
PRIORITY = ["XLR", "XLD", "XL", "AD", "AR", "MT", "YL", "SD", "COMPLEX"]  # most specific root first


def hpo_graph():
    j = json.loads((RAW / "hpo/hp.json").read_text())
    g = j["graphs"][0]
    names = {n["id"].split("/")[-1].replace("_", ":"): n.get("lbl") for n in g["nodes"]}
    par = defaultdict(set)
    for e in g["edges"]:
        if e["pred"] == "is_a":
            par[e["sub"].split("/")[-1].replace("_", ":")].add(e["obj"].split("/")[-1].replace("_", ":"))
    return names, par


def main():
    names, par = hpo_graph()
    memo = {}
    def anc(t):
        if t not in memo:
            memo[t] = set()
            s = set()
            for p in par.get(t, ()):
                s |= {p} | anc(p)
            memo[t] = s
        return memo[t]
    def klass(t):
        lin = {t} | anc(t)
        for c in PRIORITY:
            root = next(k for k, v in ROOTS.items() if v == c)
            if root in lin:
                return c
        return None
    hp = pd.read_csv(FILES["hpoa"], sep="\t", comment="#", dtype=str)
    inh = hp[hp.aspect == "I"].groupby("database_id").hpo_id.agg(lambda s: sorted(set(s))).to_dict()
    hname = hp.drop_duplicates("database_id").set_index("database_id").disease_name.to_dict()

    du = pd.read_parquet(INTER / "disease_units_clinvar.parquet")
    v = pd.read_parquet(INTER / "variants_merged.parquet")
    mi = pd.read_parquet(INTER / "mondo_index.parquet")
    omim2mondo = defaultdict(list)
    for m, o, obs in zip(mi.mondo, mi.omim_exact, mi.obsolete):
        if isinstance(o, str) and not obs:
            for x in o.split(";"):
                omim2mondo[x].append(m)
    mname = dict(zip(mi.mondo, mi.name)); morpha = dict(zip(mi.mondo, mi.orpha)); mgroup = dict(zip(mi.mondo, mi.grouping))
    pp_units = sorted({u for us in v.ppkt_units_all for u in us})
    pp_label = {}
    for us, lab in zip(v.ppkt_units_all, v.ppkt_disease_labels.fillna("")):
        if len(us) == 1 and lab and " | " not in lab:
            pp_label.setdefault(us[0], lab)
    cv_units = set(du.disease_unit)
    m2g = pd.read_csv(FILES["mim2gene"], sep="\t", comment="#", header=None, names=["mim", "t", "e", "s", "en"], dtype=str)
    mimtype = dict(zip("OMIM:" + m2g.mim, m2g.t))
    rows = du.to_dict("records")
    for u in pp_units:
        if u in cv_units:
            continue
        ms = omim2mondo.get(u, [])
        m = ms[0] if len(ms) == 1 else None
        rows.append(dict(disease_unit=u, unit_type=u.split(":")[0], disease_name=hname.get(u) or pp_label.get(u),
                         clinvar_name=None, mondo=m, mondo_name=mname.get(m), mondo_grouping=mgroup.get(m),
                         mondo_orpha=morpha.get(m), omim_mim_type=mimtype.get(u), cancer_flag=None, in_hpoa=u in hname))
    d = pd.DataFrame(rows)
    # broad-unit flag: the unit's MONDO term has >=1 is_a descendant with an exact OMIM xref
    # (e.g. OMIM:268000 'Retinitis pigmentosa' = MONDO:0019200, parent of RP1..RPn)
    children = defaultdict(set)
    for m_, p_ in zip(mi.mondo, mi.parents):
        for x in (p_.split(";") if isinstance(p_, str) else []):
            children[x].add(m_)
    has_omim = {m_ for m_, o, obs in zip(mi.mondo, mi.omim_exact, mi.obsolete) if isinstance(o, str) and not obs}
    memo_d = {}
    def n_omim_desc(t):
        if t in memo_d:
            return memo_d[t]
        memo_d[t] = set()
        out = set()
        for c in children.get(t, ()):
            if c in has_omim:
                out.add(c)
            out |= n_omim_desc(c)
        memo_d[t] = out
        return out
    d["mondo_n_omim_descendants"] = [len(n_omim_desc(m_)) if isinstance(m_, str) else 0 for m_ in d.mondo]
    d["broad_unit"] = d.mondo_n_omim_descendants > 0
    d["in_clinvar"] = d.disease_unit.isin(cv_units)
    d["in_ppkt"] = d.disease_unit.isin(pp_units)
    d["ppkt_label"] = d.disease_unit.map(pp_label)
    terms, src = [], []
    for u, orpha in zip(d.disease_unit, d.mondo_orpha):
        t = inh.get(u)
        s = "hpoa_omim" if t else None
        if not t and u.startswith("MONDO:") and isinstance(orpha, str):
            t = sorted({x for o in orpha.split(";") for x in inh.get(o, [])}) or None
            s = "hpoa_orpha_via_mondo" if t else None
        terms.append(t or []); src.append(s)
    d["inheritance_terms"] = terms
    d["inheritance_source"] = src
    d["inheritance_names"] = d.inheritance_terms.map(lambda ts: ";".join(names.get(t, t) for t in ts) or None)
    d["inheritance_classes"] = d.inheritance_terms.map(lambda ts: sorted({c for c in map(klass, ts) if c}))

    def clean(cs):
        cs = set(cs)
        xl = cs & {"XL", "XLR", "XLD"}
        if xl:
            cs -= xl
            cs.add("XLR" if xl <= {"XL", "XLR"} and "XLR" in xl else "XLD" if xl <= {"XL", "XLD"} and "XLD" in xl else "XL")
        return next(iter(cs)) if len(cs) == 1 else None
    d["inheritance_clean"] = d.inheritance_classes.map(clean)
    d["inheritance_ad_only"] = d.inheritance_clean == "AD"
    d["inheritance_ar_only"] = d.inheritance_clean == "AR"
    d.to_parquet(INTER / "diseases.parquet", index=False)
    log = dict(units=len(d), units_clinvar=int(d.in_clinvar.sum()), units_ppkt=int(d.in_ppkt.sum()),
               units_ppkt_only=int((d.in_ppkt & ~d.in_clinvar).sum()),
               with_inheritance=int(d.inheritance_source.notna().sum()),
               inheritance_source={str(k): int(x) for k, x in d.inheritance_source.value_counts(dropna=False).items()},
               inheritance_clean={str(k): int(x) for k, x in d.inheritance_clean.value_counts(dropna=False).items()},
               class_combinations={";".join(k): int(x) for k, x in Counter(tuple(c) for c in d.inheritance_classes).most_common(15)},
               broad_units=int(d.broad_unit.sum()),
               hpoa_version=open(FILES["hpoa"]).read(400).split("#version: ")[1].split("\n")[0])
    save_json(log, INTER / "inheritance_log.json")
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
