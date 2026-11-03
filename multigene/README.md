# multigene — data pipeline (phenopacket-store × AlphaFold structural context)

Builds the multi-gene dataset used to ask whether AlphaFold-derived structural-context
features of a patient's variant explain within-gene phenotypic variability.
No statistics are done here. Large/raw data live outside the repo, under
`$PHENOBRIDGE_DATA` (default `/home/user/data_ext`); the documented outputs are in
`$PHENOBRIDGE_DATA/derived/` (see `derived/README.md` there for every column and count).

## Inputs (not in repo)
* `data_ext/ppkt/.../notebooks/<COHORT>/phenopackets/*.json` — phenopacket-store v0.1.27
* `data_ext/raw/hpo/hp.json` — `curl -sSL -o hp.json https://purl.obolibrary.org/obo/hp.json`
* `data_ext/raw/alphafold/UP000005640_9606_HUMAN_v6.tar` — `curl -sS -C - -o ... https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v6.tar` (4.8 GB)
* AlphaMissense is streamed and filtered by `s08` (full file never stored).
* REST calls (Ensembl VEP/lookup/sequence, UniProt, InterPro) are cached in `data_ext/cache/`; with a warm cache everything reruns offline.

## Run (in this order, from this directory)
```bash
python3 s01_parse_phenopackets.py   # individuals / variants / explicit phenotypes                  (~3 s)
python3 s02_hpo_propagate.py        # true-path propagation                                          (~10 s)
python3 s03_vep_annotate.py         # Ensembl VEP REST on MANE Select + exon structure + NMD rule    (~20 min cold, 4 parallel POSTs; ~5 s cached)
python3 s04_uniprot.py              # UniProt accession / sequence / features + Pfam (InterPro)      (~25 min cold)
python3 s05_alphafold_extract.py    # list tar, fragment status, extract single-fragment PDBs;
                                    #   accessions missing from the proteome tar (GNAS P63092) are range-fetched
                                    #   from swissprot_pdb_v6.tar by fetch_af_from_swissprot_tar.py      (~1-3 min)
python3 s06_select_genes.py         # merge, hgvs.p/hgvs.c cross-check, MANE->UniProt/AF position maps, selection
python3 s08_alphamissense.py        # stream + filter AlphaMissense to our accessions                 (~2 min cold)
python3 s07_residue_features.py     # per-residue features for the selected genes (4 processes, ~1 min; --all = every gene)
python3 s09_assemble.py             # final tables + counts.json                                      (~40 s)
python3 verify_rwr.py               # fast RWR == mechanistic_mvp.rwr on ANKRD11 (exit 1 on mismatch)
```
`s08` must run before `s07` (per-residue AlphaMissense summaries). `run_all.sh` runs the whole chain.

## Files
* `common.py` — paths, REST disk cache, import of `mechanistic_mvp` with `gradio` stubbed.
* `s0*_*.py` — the steps above; `fetch_af_from_swissprot_tar.py` — HTTP-range binary search in the 27 GB Swiss-Prot AF tar.
* `structfeat.py` — per-residue features. Graph/regions/backbone use the repo's own
  functions (`build_confidence_graph`, `region_nodes`, `nodes_by_unip_range`,
  `shortest_path_backbone`); RWR is solved exactly with one sparse LU per protein and α
  (`p_s = α (I − (1−α) A D⁻¹)⁻¹ e_s`), and region masses are obtained for all seeds at once
  from `Mᵀ w = α 1_region`. Verified against `mechanistic_mvp.rwr` (max |Δp| 2e-11, max
  |Δfingerprint| 2e-10 on ANKRD11, α = 0.15/0.3/0.5).
* Backbone anchors: repo rule (N = residues 1–200, C = last 200) when L ≥ 600; for L < 600 the
  first/last ⌊L/3⌋ residues (guard added here).
