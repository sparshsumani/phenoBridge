#!/usr/bin/env bash
# Rebuild the allelic-heterogeneity dataset ($PHENOBRIDGE_DATA/derived_allelic, default /home/user/data_ext).
# Requires the multigene pipeline outputs ($PHENOBRIDGE_DATA/derived, intermediate/af_tar_listing.txt) and the
# AlphaFold v6 human proteome tar.  Downloads and REST responses are cached, so a rerun is offline
# except for files deliberately deleted.  Scripts run in isolated mode (python3 -I).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p "${PHENOBRIDGE_DATA:-/home/user/data_ext}/derived_allelic/intermediate"
for s in a01_download a02_clinvar_filter a03_disease_labels a04_protein_mapping a05_ppkt_merge a06_inheritance \
         a07_genes a08_alphamissense a09_residue_features a10_assemble a11_sanity; do
  echo "== $s"; t0=$SECONDS; python3 -I "$s.py"; echo "   $((SECONDS - t0)) s"
done
