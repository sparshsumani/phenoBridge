#!/usr/bin/env python3
"""Step 2: HPO ontology propagation of explicit phenotype calls.

Rules (true-path rule):
  * observed term t  -> t and all its is_a ancestors are observed (1)
  * excluded term t  -> t and all its is_a descendants are excluded (0)
  * a cell (individual, term) that receives both 1 and 0 is a conflict -> dropped, counted.
Terms: obsolete/alt ids are remapped to the current primary id (alt_id / replaced_by) before
propagation. Root HP:0000001 ('All') is not emitted.
To keep the table finite and useful, propagated cells are kept only for terms that are
observed (after propagation) in >=1 individual of the same cohort ("cohort-relevant terms");
excluded-descendant cells for terms nobody in the cohort has are not useful for within-gene
contrasts.  This restriction cannot remove a conflict (a conflict implies an observed call).
Outputs INTER/phenotypes_propagated.parquet, INTER/hpo_terms.parquet, INTER/hpo_log.json
"""
import json
from collections import defaultdict, Counter
import pandas as pd
from common import RAW, INTER

HP_JSON = RAW / "hpo/hp.json"


def cid(iri):
    return iri.rsplit("/", 1)[-1].replace("_", ":")


def load_hpo():
    g = json.load(open(HP_JSON))["graphs"][0]
    version = g.get("meta", {}).get("version")
    labels, obsolete, alt, replaced = {}, set(), {}, {}
    for n in g["nodes"]:
        if not n["id"].startswith("http://purl.obolibrary.org/obo/HP_"):
            continue
        t = cid(n["id"])
        labels[t] = n.get("lbl")
        m = n.get("meta", {})
        if m.get("deprecated"):
            obsolete.add(t)
        for bp in m.get("basicPropertyValues", []):
            if bp["pred"].endswith("hasAlternativeId"):
                alt[bp["val"]] = t
            if bp["pred"].endswith("IAO_0100001"):  # term replaced by
                replaced[t] = bp["val"] if bp["val"].startswith("HP:") else cid(bp["val"])
    parents = defaultdict(set)
    children = defaultdict(set)
    for e in g["edges"]:
        if e["pred"] != "is_a":
            continue
        s, o = cid(e["sub"]), cid(e["obj"])
        if s.startswith("HP:") and o.startswith("HP:"):
            parents[s].add(o)
            children[o].add(s)
    return version, labels, obsolete, alt, replaced, parents, children


def closure(rel):
    memo = {}

    def go(t):
        if t in memo:
            return memo[t]
        out = set()
        stack = list(rel.get(t, ()))
        while stack:
            x = stack.pop()
            if x in out:
                continue
            out.add(x)
            stack.extend(rel.get(x, ()))
        memo[t] = out
        return out
    return go


def main():
    version, labels, obsolete, alt, replaced, parents, children = load_hpo()
    anc = closure(parents)
    desc = closure(children)
    ph = pd.read_parquet(INTER / "phenotypes_explicit.parquet")
    log = Counter()
    log["hpo_version"] = version

    def remap(t):
        if t in alt:
            log["remapped_alt_id"] += 1
            t = alt[t]
        if t in replaced:
            log["remapped_replaced_by"] += 1
            t = replaced[t]
        return t
    ph["hpo_id_current"] = ph.hpo_id.map(remap)
    unknown = ph[~ph.hpo_id_current.isin(labels.keys())]
    log["explicit_cells_unknown_term"] = len(unknown)
    log["explicit_cells_obsolete_unreplaced"] = int(ph.hpo_id_current.isin(obsolete).sum())

    rows = []
    for (pid, coh), sub in ph.groupby(["phenopacket_id", "cohort"], sort=False):
        obs, exc = set(), set()
        for t, v in zip(sub.hpo_id_current, sub.value):
            if t not in labels:
                continue
            if v == 1:
                obs.add(t)
                obs |= anc(t)
            else:
                exc.add(t)
                exc |= desc(t)
        obs.discard("HP:0000001"); exc.discard("HP:0000001")
        conf = obs & exc
        log["conflict_cells_dropped"] += len(conf)
        if conf:
            log["individuals_with_conflict"] += 1
        explicit_ids = set(sub.hpo_id_current)
        for t in obs - conf:
            rows.append((pid, coh, t, 1, t in explicit_ids))
        for t in exc - conf:
            rows.append((pid, coh, t, 0, t in explicit_ids))
    df = pd.DataFrame(rows, columns=["phenopacket_id", "cohort", "hpo_id", "value", "is_explicit"])
    log["propagated_cells_unrestricted"] = len(df)
    log["propagated_observed_unrestricted"] = int((df.value == 1).sum())
    relevant = df[df.value == 1][["cohort", "hpo_id"]].drop_duplicates()
    df = df.merge(relevant, on=["cohort", "hpo_id"], how="inner")
    log["propagated_cells"] = len(df)
    log["propagated_observed"] = int((df.value == 1).sum())
    log["propagated_excluded"] = int((df.value == 0).sum())
    df["hpo_label"] = df.hpo_id.map(labels)
    df["value"] = df.value.astype("int8")
    df.to_parquet(INTER / "phenotypes_propagated.parquet", index=False)
    ph.to_parquet(INTER / "phenotypes_explicit.parquet", index=False)
    pd.DataFrame({"hpo_id": list(labels), "label": [labels[t] for t in labels],
                  "n_ancestors": [len(anc(t)) for t in labels],
                  "n_descendants": [len(desc(t)) for t in labels],
                  "obsolete": [t in obsolete for t in labels]}).to_parquet(INTER / "hpo_terms.parquet", index=False)
    (INTER / "hpo_log.json").write_text(json.dumps(dict(log), indent=1))
    print(json.dumps(dict(log), indent=1))


if __name__ == "__main__":
    main()
