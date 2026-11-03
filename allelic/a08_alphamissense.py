#!/usr/bin/env python3
"""Step 8: AlphaMissense (Cheng et al. 2023; AlphaMissense_aa_substitutions.tsv.gz, UniProt canonical
isoforms, same source as multigene/s08) downloaded (resumable curl -C -, to a temporary file that is deleted after filtering) and filtered:
  * every row of the UniProt accessions of feature proteins (genes eligible under any definition),
  * plus, for all other variants, only the exact substitution row '<acc>\t<ref><uniprot_pos><alt>'.
The full 1.2 GB file is not kept; the subset is raw/alphamissense_allelic/am_subset_allelic.tsv.gz.
Output intermediate/alphamissense_subset.parquet (uniprot_id, protein_variant, am_pathogenicity, am_class,
ref, pos, alt) and intermediate/am_log.json.  The subset is re-streamed automatically when the filter keys change.
"""
import sys, subprocess, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd
from common import INTER, URLS, FILES, save_json


def main():
    out = FILES["alphamissense"]
    out.parent.mkdir(parents=True, exist_ok=True)
    g = pd.read_parquet(INTER / "genes.parquet")
    v = pd.read_parquet(INTER / "variants_positioned.parquet")
    full = sorted(set(g[g.eligible_any_definition].uniprot.dropna()))
    keys = sorted({f"{a}\t{r}{int(p)}{t}" for a, r, p, t in zip(v.uniprot, v.mane_ref, v.uniprot_pos, v.mane_alt)
                   if isinstance(a, str) and pd.notna(p) and a not in full})
    lst = out.parent / "filter_keys.tsv"
    old_keys = lst.read_text() if lst.exists() else None
    keys_changed = False
    lst.write_text("".join(f"{a}\t*\n" for a in full) + "".join(k + "\n" for k in keys))
    keys_changed = out.exists() and old_keys is not None and old_keys != lst.read_text()
    if keys_changed:
        print("filter keys changed since the subset was streamed -> re-streaming", flush=True)
        out.unlink()
    if not out.exists():
        # resumable download of the full file to a temporary location (the proxy can reset long
        # streams), filter, then delete the full file
        full_gz = out.parent / "_full_tmp" / "AlphaMissense_aa_substitutions.tsv.gz"
        full_gz.parent.mkdir(exist_ok=True)
        for attempt in range(30):
            r = subprocess.run(["curl", "-sS", "-C", "-", "-o", str(full_gz), URLS["alphamissense"]])
            if r.returncode == 0:
                break
            print("download interrupted (curl exit", r.returncode, "), resuming", flush=True)
        else:
            raise SystemExit("AlphaMissense download failed")
        subprocess.run(["gzip", "-t", str(full_gz)], check=True)
        cmd = (f"set -o pipefail; zcat {full_gz} | "
               f"awk -F'\\t' 'NR==FNR{{if($2==\"*\")a[$1]; else b[$1\"\\t\"$2]; next}} /^#/{{next}} "
               f"$1==\"uniprot_id\"{{print;next}} ($1 in a)||(($1\"\\t\"$2) in b)' {lst} - "
               f"| gzip > {out}.tmp && mv {out}.tmp {out}")
        print("filtering AlphaMissense ...", flush=True)
        subprocess.run(["bash", "-c", cmd], check=True)
        full_gz.unlink()
    am = pd.read_csv(out, sep="\t")
    am["ref"] = am.protein_variant.str[0]
    am["alt"] = am.protein_variant.str[-1]
    am["pos"] = am.protein_variant.str[1:-1].astype(int)
    am.to_parquet(INTER / "alphamissense_subset.parquet", index=False)
    log = dict(source=URLS["alphamissense"], full_accessions_requested=len(full),
               full_accessions_found=int(am[am.uniprot_id.isin(full)].uniprot_id.nunique()),
               missing_full_accessions=sorted(set(full) - set(am.uniprot_id)),
               exact_keys_requested=len(keys), exact_rows_found=int((~am.uniprot_id.isin(full)).sum()), rows=len(am), filter_keys_changed_since_stream=bool(keys_changed))
    save_json(log, INTER / "am_log.json")
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
