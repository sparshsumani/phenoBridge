#!/usr/bin/env python3
"""Step 6: AlphaMissense (Cheng et al. 2023, AlphaMissense_aa_substitutions.tsv.gz, UniProt canonical
numbering), streamed once through awk as in allelic/a08 (full file downloaded by v01 with resumable curl):
  * exact substitution row '<acc>\t<ref><uniprot_pos><alt>' for EVERY mapped variant
    -> raw/alphamissense_vus3d/am_subset_vus3d.tsv.gz
  * for the UniProt accessions of feature proteins (v05: >=1 test variant and >=5 train-labelled residues)
    a per-residue summary over the 19 substitutions (ref aa, mean, max, n) computed in the stream
    -> raw/alphamissense_vus3d/am_residue_summary_vus3d.tsv.gz  (keeping all 52 M substitution rows of these
    proteins was avoided for disk economy; they can be re-derived from the full file)
The full 1.2 GB file is deleted after a successful filter; the subset is re-streamed when the filter
keys change (full file needed again: rerun v01).
Exact score per variant: key (uniprot, uniprot_pos, alt); AM ref must equal mane_ref.
Outputs intermediate/am_exact.parquet (variation_id, am_pathogenicity, am_class, am_status),
        intermediate/am_residue_summary.parquet (uniprot_id, pos, ref, am_mean, am_max, am_n_subst),
        intermediate/am_log.json
"""
import sys, subprocess, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from common import INTER, URLS, AM_SUBSET, AM_FULL_TMP, save_json


def main():
    g = pd.read_parquet(INTER / "genes.parquet")
    v = pd.read_parquet(INTER / "variants_mapped.parquet")
    full = sorted(set(g[g.feature_protein].uniprot.dropna()))
    fs = set(full)
    keys = sorted({f"{a}\t{r}{int(p)}{t}" for a, r, p, t in zip(v.uniprot, v.mane_ref, v.uniprot_pos, v.mane_alt)
                   if isinstance(a, str) and pd.notna(p) and isinstance(r, str)})
    lst = AM_SUBSET.parent / "filter_keys.tsv"
    summ = AM_SUBSET.parent / "am_residue_summary_vus3d.tsv.gz"
    new = "".join(f"{a}\t*\n" for a in full) + "".join(k + "\n" for k in keys)
    old = lst.read_text() if lst.exists() else None
    if AM_SUBSET.exists() and old != new:
        print("filter keys changed -> re-streaming", flush=True)
        AM_SUBSET.unlink()
    lst.write_text(new)
    if not (AM_SUBSET.exists() and summ.exists()):
        if not AM_FULL_TMP.exists():
            raise SystemExit("full AlphaMissense file missing: run v01_download.py")
        subprocess.run(["gzip", "-t", str(AM_FULL_TMP)], check=True)
        awk = ("NR==FNR{if($2==\"*\")a[$1]; else b[$1\"\\t\"$2]; next} /^#/{next} "
               "$1==\"uniprot_id\"{print > ex; next} "
               "(($1\"\\t\"$2) in b){print > ex} "
               "($1 in a){p=substr($2,2,length($2)-2); k=$1\"\\t\"p; if(!(k in n)){n[k]=0; s[k]=0; m[k]=-1; r[k]=substr($2,1,1)} "
               "n[k]++; s[k]+=$3; if($3>m[k])m[k]=$3} "
               "END{print \"uniprot_id\\tpos\\tref\\tam_mean\\tam_max\\tam_n_subst\" > sm; "
               "for(k in n) printf \"%s\\t%s\\t%.6f\\t%.4f\\t%d\\n\", k, r[k], s[k]/n[k], m[k], n[k] > sm}")
        cmd = (f"set -o pipefail; zcat {AM_FULL_TMP} | awk -F'\\t' -v ex=/dev/stdout -v sm={summ}.raw '{awk}' {lst} - "
               f"| gzip > {AM_SUBSET}.tmp && gzip -c {summ}.raw > {summ}.tmp && rm {summ}.raw && mv {AM_SUBSET}.tmp {AM_SUBSET} && mv {summ}.tmp {summ}")
        print("filtering AlphaMissense ...", flush=True)
        subprocess.run(["bash", "-c", cmd], check=True)
        AM_FULL_TMP.unlink()
    am = pd.read_csv(AM_SUBSET, sep="\t")
    am["ref"] = am.protein_variant.str[0]
    am["alt"] = am.protein_variant.str[-1]
    am["pos"] = am.protein_variant.str[1:-1].astype(np.int32)
    sm = pd.read_csv(summ, sep="\t")
    sm.sort_values(["uniprot_id", "pos"]).to_parquet(INTER / "am_residue_summary.parquet", index=False)
    # exact substitution per variant
    idx = am.drop_duplicates(["uniprot_id", "pos", "alt"]).set_index(["uniprot_id", "pos", "alt"])[["am_pathogenicity", "am_class", "ref"]]
    acc_in_am = set(am.uniprot_id.astype(str))
    rows = []
    for vid, a, p, r, t in zip(v.variation_id, v.uniprot, v.uniprot_pos, v.mane_ref, v.mane_alt):
        if not isinstance(a, str) or pd.isna(p) or not isinstance(r, str):
            rows.append((vid, np.nan, None, "no_uniprot_position")); continue
        k = (a, int(p), t)
        if k in idx.index:
            x = idx.loc[k]
            if x.ref != r:
                rows.append((vid, np.nan, None, "am_ref_mismatch")); continue
            rows.append((vid, float(x.am_pathogenicity), str(x.am_class), "ok"))
        else:
            rows.append((vid, np.nan, None, "accession_not_in_am" if a not in acc_in_am else "substitution_not_in_am"))
    ex = pd.DataFrame(rows, columns=["variation_id", "am_pathogenicity", "am_class", "am_status"])
    ex.to_parquet(INTER / "am_exact.parquet", index=False)
    vv = v[["variation_id", "set", "usable"]].merge(ex, on="variation_id")
    log = dict(source=URLS["alphamissense"], subset_file=str(AM_SUBSET), residue_summary_file=str(summ),
               full_accessions_requested=len(full), full_accessions_found=int(sm.uniprot_id.nunique()),
               missing_full_accessions=sorted(fs - set(sm.uniprot_id)), residue_summary_rows=len(sm),
               exact_keys_requested=len(keys), exact_rows_found=len(am),
               am_status_all={k: int(n) for k, n in ex.am_status.value_counts().items()},
               am_status_usable_by_set={f"{a}|{b}": int(n) for (a, b), n in vv[vv.usable].groupby(["set", "am_status"]).size().items()})
    save_json(log, INTER / "am_log.json")
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
