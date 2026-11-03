#!/usr/bin/env bash
# Rerun the full multigene pipeline (uses the REST cache in $PHENOBRIDGE_DATA/cache when present).
set -euo pipefail
cd "$(dirname "$0")"
for s in s01_parse_phenopackets s02_hpo_propagate s03_vep_annotate s04_uniprot s05_alphafold_extract \
         s06_select_genes s08_alphamissense s07_residue_features s09_assemble verify_rwr; do
  echo "== $s"; t0=$SECONDS; python3 "$s.py"; echo "   $((SECONDS - t0)) s"
done
