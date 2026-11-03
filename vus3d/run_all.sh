#!/usr/bin/env bash
# Rebuild the temporal ClinVar "3D label propagation" VUS dataset ($PHENOBRIDGE_DATA/derived_vus3d,
# default /home/user/data_ext).  Requires the AlphaFold v6 human proteome tar and
# intermediate/af_tar_listing.txt (multigene/s05), MANE + UniProt files (allelic/a01).  Downloads and REST
# responses are cached, so a rerun is offline except for deliberately deleted files.  Scripts run in
# isolated mode (python3 -I) because they read downloaded files.
# NOTE: v06 deletes the full AlphaMissense file after filtering; if the filter keys change (new feature
# proteins / variants) v01 re-downloads it.  No classifier is trained or evaluated here.
set -euo pipefail
cd "$(dirname "$0")"
for s in v01_download v02_clinvar_parse v03_sets v04_protein_mapping v05_genes v06_alphamissense \
         v07_residue_features v08_assemble v09_sanity; do
  echo "== $s"; t0=$SECONDS; python3 -I "$s.py"; echo "   $((SECONDS - t0)) s"
done
