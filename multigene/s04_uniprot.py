#!/usr/bin/env python3
"""Step 4: map each gene (HGNC id from phenopacket geneContext) to its reviewed human UniProt
entry (preferring the entry carrying a MANE-Select cross-reference), fetch sequence and
features, and Pfam domain locations from the InterPro API. Cached under CACHE/uniprot*, CACHE/interpro.

Outputs INTER/gene_uniprot.parquet, INTER/uniprot_features.parquet, INTER/uniprot_log.json
"""
import json, urllib.parse
from collections import Counter
import pandas as pd
from common import INTER, cached_request, http_json

FEATURE_TYPES = ["Domain", "Region", "Repeat", "Zinc finger", "Topological domain", "Transmembrane",
                 "Intramembrane", "Motif", "Coiled coil", "DNA binding", "Compositional bias"]


def uniprot_search(q):
    url = "https://rest.uniprot.org/uniprotkb/search?format=json&size=25&query=" + urllib.parse.quote(q)
    return cached_request("uniprot_search", url, lambda: http_json(url))


def interpro_pfam(acc):
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/?page_size=200"
    def f():
        import requests
        for _ in range(5):
            r = requests.get(url, headers={"Accept": "application/json"}, timeout=120)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 204:
                return {"count": 0, "results": []}
        return None
    return cached_request("interpro_pfam", url, f)


def mane_of(e):
    for x in e.get("uniProtKBCrossReferences", []):
        if x["database"] == "MANE-Select":
            p = {k["key"]: k["value"] for k in x.get("properties", [])}
            return x["id"], p.get("RefSeqNucleotideId"), p.get("RefSeqProteinId"), p.get("ProteinId")
    return None, None, None, None


def main():
    var = pd.read_parquet(INTER / "variants_raw.parquet")
    genes = var.dropna(subset=["hgnc_id"]).drop_duplicates("hgnc_id")[["gene", "hgnc_id"]]
    log = Counter()
    rows, feats = [], []
    for g, h in zip(genes.gene, genes.hgnc_id):
        num = h.split(":")[1]
        res = uniprot_search(f"xref:hgnc-{num} AND reviewed:true AND organism_id:9606")
        hits = res.get("results", []) if res else []
        how = "hgnc_xref"
        if not hits:
            res = uniprot_search(f"gene_exact:{g} AND reviewed:true AND organism_id:9606")
            hits = res.get("results", []) if res else []
            how = "gene_exact"
        if not hits:
            log["no_uniprot"] += 1
            rows.append(dict(gene=g, hgnc_id=h, uniprot=None, map_method="none", n_reviewed_hits=0))
            continue
        with_mane = [e for e in hits if mane_of(e)[0]]
        e = with_mane[0] if with_mane else hits[0]
        if len(hits) > 1:
            log["multiple_reviewed_hits"] += 1
        enst, nm, np_, ensp = mane_of(e)
        acc = e["primaryAccession"]
        rows.append(dict(gene=g, hgnc_id=h, uniprot=acc, map_method=how, n_reviewed_hits=len(hits),
                         other_hits=",".join(x["primaryAccession"] for x in hits if x is not e),
                         uniprot_length=e["sequence"]["length"], sequence=e["sequence"]["value"],
                         uniprot_seq_version=e.get("entryAudit", {}).get("sequenceVersion"),
                         uniprot_mane_enst=enst, uniprot_mane_nm=nm, uniprot_mane_np=np_, uniprot_mane_ensp=ensp,
                         protein_name=e.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")))
        for f in e.get("features", []):
            if f["type"] in FEATURE_TYPES:
                s, t = f["location"]["start"].get("value"), f["location"]["end"].get("value")
                if s is None or t is None:
                    continue
                feats.append(dict(uniprot=acc, source="UniProt", type=f["type"], description=f.get("description"),
                                  start=int(s), end=int(t), feature_id=None))
        pf = interpro_pfam(acc)
        if pf is None:
            log["interpro_fail"] += 1
        else:
            for r in pf.get("results", []):
                md = r["metadata"]
                for p in r.get("proteins", []):
                    for loc in p.get("entry_protein_locations", []):
                        for fr in loc.get("fragments", []):
                            feats.append(dict(uniprot=acc, source="Pfam", type="Pfam", description=md.get("name"),
                                              start=int(fr["start"]), end=int(fr["end"]), feature_id=md["accession"]))
    gu = pd.DataFrame(rows)
    ft = pd.DataFrame(feats)
    gu.to_parquet(INTER / "gene_uniprot.parquet", index=False)
    ft.to_parquet(INTER / "uniprot_features.parquet", index=False)
    log["genes"] = len(gu); log["mapped"] = int(gu.uniprot.notna().sum())
    log["with_mane_xref"] = int(gu.uniprot_mane_enst.notna().sum())
    (INTER / "uniprot_log.json").write_text(json.dumps(dict(log), indent=1))
    print(json.dumps(dict(log), indent=1))


if __name__ == "__main__":
    main()
