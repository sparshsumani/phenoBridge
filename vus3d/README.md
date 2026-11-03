# vus3d — temporal ClinVar dataset for "3D label propagation" of missense VUS

Builds a clean, leakage-aware dataset: **TRAIN** = ClinVar missense P/LP and B/LB (≥1★) as of the
**2023-01** release; **TEST** = variants that were VUS, conflicting or absent in 2023-01 and are P/LP or
B/LB (≥1★) in the **2026-11-01** release. All variants are mapped to MANE / UniProt / AlphaFold v6
numbering; AlphaMissense exact-substitution scores are attached; every protein with ≥1 test variant and
≥5 train-labelled residues gets an 8 Å Cα contact graph and per-residue structural features.
**No classifier is trained or evaluated here and no feature or score is compared between P and B test
variants** — only descriptive counts. Dataset documentation (every column, every count):
`$PHENOBRIDGE_DATA/derived_vus3d/README.md` and `counts.json` (default `/home/user/data_ext`).

Built 2026-11-03. Reuses `../allelic/` (a02 filters/star map/Name parser, a04 MANE/UniProt/VEP helpers,
a09 per-protein features + domain annotation, a08 AlphaMissense streaming approach) and through it
`../multigene/` (VEP REST client + cache, align_map, structfeat exact RWR, s07 annotate/am_mapping).

## Run
```bash
./run_all.sh     # v01 … v09, each `python3 -I vNN_*.py`
```
Isolated mode (`-I`): the scripts read downloaded files; code lives here, downloads under `data_ext/raw/`.

| step | script | what | time (this build) |
|---|---|---|---|
| 1 | `v01_download.py` | ClinVar 2023-01 archive + current release, AlphaMissense (resumable curl) + provenance | 1 min cold |
| 2 | `v02_clinvar_parse.py` | both releases: GRCh38 status tables (all types) + missense tables (a02 filters except classification) | 6 min |
| 3 | `v03_sets.py` | key agreement, TRAIN / TEST / FLIPPED / excluded sets | 4 min |
| 4 | `v04_protein_mapping.py` | MANE / VEP / UniProt / AlphaFold mapping, ref-aa check, AF model extraction | 6 min (+11 min cold: VEP + tar) |
| 5 | `v05_genes.py` | per-gene counts, evaluable thresholds, same-residue flags | 10 s |
| 6 | `v06_alphamissense.py` | awk stream of AlphaMissense: exact rows + per-residue summaries | 10 min |
| 7 | `v07_residue_features.py` | contact graphs + residue features for feature proteins (4 processes, resumable) | ~45 min |
| 8 | `v08_assemble.py` | `variants.parquet`, `genes.parquet`, `counts.json` | 1 min |
| 9 | `v09_sanity.py` | sanity checks, ClinVar VCV spot-check, RWR verification → `counts.json["sanity"]` | 1 min |

## Code layout
* `common.py` — paths (`OUT` = `derived_vus3d`, `INTER`, `RFDIR`, `GRAPHDIR`), release files, thresholds;
  loads `allelic/common.py` as `_allelic_common` and re-exports its constants; `al_import()` imports allelic
  `aNN` modules with the right `common`. Analysis helpers (same semantics as `allelic/common.py`, reading
  `derived_vus3d/`):
  ```python
  import sys; sys.path.insert(0, "/home/user/phenoBridge/vus3d"); import common as C
  pos, xyz, plddt = C.ca_coords("P02545")
  pos, W, A = C.contact_graph("P02545")          # W = A D^-1 (D = max(1, degree)), csc
  p = C.rwr_profile("P02545", 482, alpha=0.3)   # Series residue -> stationary prob. (exact solve, sums to 1)
  x = C.rwr_matrix_solve("P02545", b, 0.3)      # alpha (I-(1-alpha)W)^-1 b for any seed vector b (multi-seed)
  rf = C.residue_table("P02545")                # per-residue features
  ```
* `hse_fast.py` — vectorised, value-identical re-implementation of Bio.PDB `HSExposureCA/CB` (the
  pure-Python Biopython loop was ~75–95 % of a09's per-protein time). Verified: 0 mismatches on 64,184
  stored a09 values (12 allelic proteins), and full v07 worker output identical (all 39 columns + edge list)
  to `allelic/a09.one` on 2 proteins.
* `run_vus3d.py`, `PLAN.md`, `notes/` — written by the lead (analysis; not part of the data build).

## Key decisions (details and counts in the derived README)
* Classification = germline `ClinicalSignificance` string → PLP / BLB / VUS / CONFLICT / OTHER (exact
  strings as allelic/a02; "Pathogenic; risk factor", "low penetrance" etc. are OTHER).
* 2023 status of a test candidate is looked up in the **full** 2023 GRCh38 table (any type / stars / Name),
  by VariationID, else by GRCh38 SNV coordinate (2023 records later replaced by a new VariationID), else ABSENT.
* Excluded from TEST (counted): 2023 P/LP/B/LB with 0★ or outside the 2023 missense table, and 2023 OTHER.
* FLIPPED (P↔B between releases) are excluded from both sets.
* Protein change from the label-defining release's Name (train 2023, test 2026); position via MANE
  (VEP for non-MANE Names); only single-fragment AlphaFold v6 models; ref aa must match the model.

## Known limitations
* ~17 % of TEST variants have a 2026 aggregate LastEvaluated on/before 2023-01-31 (classification
  existed before the cutoff but was submitted to ClinVar later, or a VUS-era date carried forward);
  the column `last_evaluated_2026_after_cutoff` allows a stricter test set.
* TEST is dominated by prior-ABSENT benign variants (bulk lab submissions after 2023).
* Fragmented AlphaFold proteins (>2,700 aa; 210 genes incl. TTN, RYR1, FBN1) are excluded.
