#!/usr/bin/env python3
"""Step 5: list the AlphaFold human proteome tar (v6) once, determine per UniProt accession
whether a model exists and how many fragments it has, and extract the PDB of every
single-fragment model for all mapped genes into RAW/alphafold/models/.

Requires RAW/alphafold/UP000005640_9606_HUMAN_v6.tar (downloaded from
https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v6.tar).
Outputs INTER/af_models.parquet, INTER/af_tar_listing.txt
"""
import re, subprocess, json
from collections import defaultdict, Counter
import pandas as pd
from common import INTER, AF_TAR, AF_DIR

LISTING = INTER / "af_tar_listing.txt"


def main():
    if not LISTING.exists():
        subprocess.run(f"tar -tf {AF_TAR} > {LISTING}.tmp && mv {LISTING}.tmp {LISTING}", shell=True, check=True)
    names = LISTING.read_text().split()
    kinds = Counter(re.sub(r"^AF-[A-Z0-9]+-F\d+-", "", n) for n in names)
    frags = defaultdict(set)
    for n in names:
        m = re.match(r"AF-([A-Z0-9]+)-F(\d+)-model_v(\d+)\.pdb\.gz$", n)
        if m:
            frags[m.group(1)].add(int(m.group(2)))
    gu = pd.read_parquet(INTER / "gene_uniprot.parquet")
    rows, want = [], []
    for g, acc in zip(gu.gene, gu.uniprot):
        if acc is None or (isinstance(acc, float)):
            rows.append(dict(gene=g, uniprot=None, af_n_fragments=0, af_fragmented=None, af_pdb=None)); continue
        n = len(frags.get(acc, ()))
        src = "proteome_tar" if n else None
        member = f"AF-{acc}-F1-model_v6.pdb.gz" if n == 1 else None
        if member:
            want.append(member)
        elif n == 0:
            # not in the proteome set (e.g. GNAS P63092): try the Swiss-Prot tar via HTTP range requests
            try:
                from fetch_af_from_swissprot_tar import fetch
                if not (AF_DIR / f"AF-{acc}-F1-model_v6.pdb.gz").exists():
                    fetch(acc)
                member, n, src = f"AF-{acc}-F1-model_v6.pdb.gz", 1, "swissprot_tar_range_request"
            except BaseException as e:  # SystemExit when absent
                print("no AF model for", g, acc, e)
        rows.append(dict(gene=g, uniprot=acc, af_n_fragments=n, af_fragmented=(n > 1) if n else None,
                         af_source=src, af_pdb=str(AF_DIR / member) if member else None))
    todo = [w for w in want if not (AF_DIR / w).exists()]
    if todo:
        lst = INTER / "af_extract_list.txt"
        lst.write_text("\n".join(todo) + "\n")
        subprocess.run(["tar", "-xf", str(AF_TAR), "-C", str(AF_DIR), "-T", str(lst)], check=True)
    df = pd.DataFrame(rows)
    df["af_pdb_exists"] = [bool(p) and (AF_DIR / p.split("/")[-1]).exists() for p in df.af_pdb.fillna("")]
    df.to_parquet(INTER / "af_models.parquet", index=False)
    log = dict(tar_members=len(names), member_kinds=dict(kinds), accessions_in_tar=len(frags),
               genes=len(df), genes_no_model=int((df.af_n_fragments == 0).sum()),
               genes_fragmented=int((df.af_n_fragments > 1).sum()), genes_single=int((df.af_n_fragments == 1).sum()),
               pae_files_in_tar=sum(1 for n in names if "predicted_aligned_error" in n or "pae" in n.lower()))
    (INTER / "af_log.json").write_text(json.dumps(log, indent=1))
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
