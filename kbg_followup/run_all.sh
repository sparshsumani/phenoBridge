#!/usr/bin/env bash
# Re-runs the whole KBG follow-up. REST/AlphaMissense results are cached in cache/.
set -euo pipefail
cd "$(dirname "$0")"
bash 00_fetch_alphamissense.sh   # AlphaMissense, Q6UB99 only
python3 stats_utils.py            # self-test of ML / Firth implementations
python3 01_annotate.py            # phenopacket match, VEP (MANE), exons, NMD, SV parsing, regex comparison
python3 02_features.py            # fingerprints (all residues, ~2 min first time), per-individual features
python3 03_duplicates.py          # relatedness / duplicate clusters
python3 04_analysis.py            # pre-specified analyses (ANALYSIS_PLAN.md), ~5 min (5,000 permutations)
python3 04b_exploratory_disorder_ID.py   # EXPLORATORY follow-up of the one secondary hit
python3 05_power.py               # simulation power analysis
python3 06_figures.py             # figures -> fig/
python3 07_tables.py             # markdown tables for RESULTS.md
