#!/usr/bin/env bash
# Streams the 1.2 GB AlphaMissense table and keeps only ANKRD11 (Q6UB99); the full file is never stored.
set -euo pipefail
cd "$(dirname "$0")/cache"
[ -s uniprot_Q6UB99.json ] || curl -sS "https://rest.uniprot.org/uniprotkb/Q6UB99.json" -o uniprot_Q6UB99.json
[ -s alphamissense_Q6UB99.tsv ] && exit 0
curl -sS https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz | zcat \
  | awk -F'\t' 'NR<=4 || $1=="Q6UB99"' > alphamissense_Q6UB99.tsv.part
mv alphamissense_Q6UB99.tsv.part alphamissense_Q6UB99.tsv
