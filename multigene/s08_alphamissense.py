#!/usr/bin/env python3
"""Step 8: AlphaMissense (Cheng et al. 2023, aa_substitutions table, UniProt canonical isoforms).
Streams the 1.2 GB gz from Google Storage and keeps only rows for the UniProt accessions of all
mapped genes -> RAW/alphamissense/am_subset.tsv.gz -> INTER/alphamissense_subset.parquet.
(The full file is never stored.)  Columns of source: uniprot_id, protein_variant, am_pathogenicity, am_class.
"""
import subprocess
import pandas as pd
from common import INTER, RAW

URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz"
OUT = RAW / "alphamissense/am_subset.tsv.gz"


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    gu = pd.read_parquet(INTER / "gene_uniprot.parquet")
    accs = sorted(set(gu.uniprot.dropna()))
    lst = OUT.parent / "accessions.txt"
    lst.write_text("\n".join(accs) + "\n")
    if not OUT.exists():
        cmd = (f"set -o pipefail; curl -sS --retry 5 {URL} | zcat | "
               f"awk -F'\\t' 'NR==FNR{{a[$1];next}} /^#/{{next}} $1==\"uniprot_id\"{{print;next}} ($1 in a)' {lst} - "
               f"| gzip > {OUT}.tmp && mv {OUT}.tmp {OUT}")
        subprocess.run(["bash", "-c", cmd], check=True)
    am = pd.read_csv(OUT, sep="\t")
    am["ref"] = am.protein_variant.str[0]
    am["alt"] = am.protein_variant.str[-1]
    am["pos"] = am.protein_variant.str[1:-1].astype(int)
    am.to_parquet(INTER / "alphamissense_subset.parquet", index=False)
    print(len(am), am.uniprot_id.nunique(), "of", len(accs), "accessions present")


if __name__ == "__main__":
    main()
