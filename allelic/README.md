# allelic — data pipeline for within-gene allelic heterogeneity (ClinVar × phenopacket-store × AlphaFold)

Builds a clean dataset of pathogenic missense variants labelled with the specific disease they cause,
mapped onto AlphaFold v6 models, with PhenoBridge per-residue structural features for every residue
of every eligible protein. **No association test, no between-disease feature comparison and no look at
feature distributions by disease label is done here** — only descriptive counts (variants per gene /
disease / source). Outputs: `$PHENOBRIDGE_DATA/derived_allelic/` (default `/home/user/data_ext`), see
`derived_allelic/README.md` there for every column and `counts.json` for every count.

Built 2026-11-03. Reuses `../multigene/` (VEP REST client + cache `s03`, UniProt/InterPro queries `s04`,
AF sequence/alignment `s06`, domain annotation/AlphaMissense mapping `s07`, `structfeat.py` which uses
`mechanistic_mvp.py` functions) through `common.mg_import` (multigene modules see `multigene/common.py`).

## Run
```bash
./run_all.sh            # a01 … a11, each `python3 -I aNN_*.py`; ~15 min with warm caches
```
Isolated mode (`-I`) is used because the scripts read downloaded files; the scripts live in this repo,
downloads live under `data_ext/raw/` and `data_ext/cache/`.

| step | script | what | time (warm) |
|---|---|---|---|
| 1 | `a01_download.py` | download + provenance (Last-Modified, size, sha256, in-file version) → `intermediate/downloads.json` | 10 s |
| 2 | `a02_clinvar_filter.py` | ClinVar filters F1–F7 | 90 s |
| 3 | `a03_disease_labels.py` | condition parsing, harmonisation to disease units | 80 s |
| 4 | `a04_protein_mapping.py` | MANE / VEP / UniProt / AlphaFold mapping, ref-aa check, AF model extraction | 60 s |
| 5 | `a05_ppkt_merge.py` | phenopacket-store missense + protein-level merge with ClinVar | 150 s |
| 6 | `a06_inheritance.py` | disease table + HPO inheritance | 3 s |
| 7 | `a07_genes.py` | AF/UniProt positions of merged variants, per-gene eligibility counts | 40 s |
| 8 | `a08_alphamissense.py` | AlphaMissense download → filtered subset | 5 s (3 min cold) |
| 9 | `a09_residue_features.py` | residue features + contact graphs for the 99 feature proteins (4 processes) | 7 min |
| 10 | `a10_assemble.py` | final tables + `counts.json` | 10 s |
| 11 | `a11_sanity.py` | sanity checks, ClinVar spot-check, RWR verification → `counts.json["sanity"]` | 5 s |

## Inputs (URL → local file; versions as downloaded 2026-11-03)
| source | URL | local | version |
|---|---|---|---|
| ClinVar | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz | `raw/clinvar/variant_summary.txt.gz` | Last-Modified Sun 01 Nov 2026 14:17:51 GMT, 457,862,838 B |
| MONDO | https://purl.obolibrary.org/obo/mondo.obo | `raw/mondo/mondo.obo` | data-version releases/2026-10-06 |
| HPO annotations | https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa | `raw/hpoa/phenotype.hpoa` | #version 2026-09-02 (hp 2026-09-01) |
| HPO ontology | (multigene) https://purl.obolibrary.org/obo/hp.json | `raw/hpo/hp.json` | 2026-09-01 |
| OMIM mim2gene | https://omim.org/static/omim/data/mim2gene.txt | `raw/omim/mim2gene.txt` | Generated 2026-11-02 (MIM entry types only) |
| MANE | https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz and `…v1.5.refseq_protein.faa.gz` | `raw/mane/` | MANE v1.5 |
| UniProt | `https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=(organism_id:9606) AND (reviewed:true)&fields=accession,id,gene_primary,gene_synonym,xref_hgnc,xref_mane-select,length,sequence,sequence_version` | `raw/uniprot/uniprot_human_reviewed.tsv` | release 2026_03 (20,431 entries) |
| UniProt features / Pfam | rest.uniprot.org search JSON (s04 query), InterPro API `/entry/pfam/protein/uniprot/<acc>` | `cache/uniprot_search`, `cache/interpro_pfam` | live, cached |
| Ensembl VEP | https://rest.ensembl.org/vep/human/region (s03 parameters `mane=1&canonical=1&hgvs=1&numbers=1&uniprot=1&protein=1&variant_class=1`), POST ≤200 variants | `cache/vep_region/` | GRCh38 REST, cached |
| AlphaFold DB | v6 human proteome tar (multigene) | `raw/alphafold/UP000005640_9606_HUMAN_v6.tar` → only needed `AF-*-F1-model_v6.pdb.gz` extracted into `raw/alphafold/models/` | v6 |
| AlphaMissense | https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz | filtered subset `raw/alphamissense_allelic/am_subset_allelic.tsv.gz` (full file deleted after filtering) | Cheng et al. 2023 |
| ClinVar live record (spot check) | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=… | `cache/clinvar_esummary/` | live |
| phenopacket-store | v0.1.27 derived tables from multigene | `derived/variants.parquet`, `derived/individuals.parquet` | v0.1.27 |

## Step details and counts

### a02 — ClinVar filters (rows / distinct VariationIDs; `counts.json → clinvar_filter`)
| filter | VariationIDs after | dropped |
|---|---:|---:|
| F0 all rows of variant_summary | 4,707,240 (9,345,867 rows) | |
| F1 Assembly == GRCh38 | 4,636,853 | 70,387 |
| F2 Type == single nucleotide variant | 4,285,982 | 350,871 |
| F3 ClinicalSignificance (germline) exactly Pathogenic / Likely pathogenic / Pathogenic/Likely pathogenic | 198,070 | 4,087,912 (incl. 176,092 conflicting; P/LP with extra terms such as "; risk factor", "low penetrance" excluded) |
| F4 OriginSimple != somatic | 197,471 | 599 |
| F5 review status ≥ 1 star (drops no assertion criteria / no classification provided) | 180,061 | 17,410 |
| F6 Name parses as `NM_x(GENE):c.… (p.…)` (protein change present) | 135,045 | 45,016 (non-coding/intronic/splice, no p.) |
| F7 missense (excl. synonymous 692, nonsense 72,491, start-loss Met1 1,635, stop-loss 188 rows) | **60,056** | 74,989 |

Stars: 1★ 34,864 · 2★ 21,730 · 3★ 3,454 · 4★ 8. Genes: 3,478. Recorded per variant: stars, review
status, NumberSubmitters, LastEvaluated, classification, RCV, rsID, VCF coordinates.

### a03 — disease labels → disease units
Field layout: `|` separates RCV condition sets, `;` conditions inside one set, `,` IDs of one condition.
Rules (in order) per condition, see the docstring for regexes:
1. `nonspecific_label` – not provided / not specified / See cases / Inborn genetic diseases / Hereditary disease … (28,864 condition links)
2. `generic_gene_related` – `<GENE>-related disorder(s)/disease/condition` without an OMIM phenotype MIM (11,658; the variant is flagged `clinvar_has_generic_condition`)
3. `cancer_no_omim` – cancer/somatic-type label (name regex, or MONDO descendant of cancer-or-benign-tumor / hereditary neoplastic syndrome) without a phenotype MIM, e.g. "Hereditary cancer-predisposing syndrome" (5,806)
4. `cancer_somatic_type` – cancer label with a MIM but not a hereditary-predisposition label (e.g. "Gastric cancer", "Acute myeloid leukemia", "…, somatic") (660)
5. `specific` – has an OMIM **phenotype** MIM (typed with mim2gene: phenotype / predominantly phenotypes) → unit `OMIM:<mim>`; else MONDO with exactly one exact (`MONDO:equivalentTo`) OMIM phenotype xref → that OMIM; else a MONDO term that is not a grouping → unit `MONDO:<id>` (53,704). Cancer-type conditions kept here (hereditary predisposition with a MIM, e.g. Li-Fraumeni 1, Lynch 1, NF1) are flagged `cancer_flag`.
6. `grouping_no_omim` – MONDO grouping (subset disease_grouping/rare_grouping, OMIM phenotypic series, >1 exact OMIM, or any is_a descendant with an exact OMIM xref, e.g. "Long QT syndrome", "RASopathy", "Intellectual disability") → not a unit (17,643)
7. `no_omim_no_mondo` – MedGen/HPO/Orphanet-only conditions (e.g. "Cardiovascular phenotype") → dropped (7,951)

Per variant: union of specific units; a MONDO-based unit that is an is_a ancestor of another unit of
the same variant is pruned (697 variants). Of 60,056 ClinVar variants: 32,494 single-label, 8,168
multi-label, 19,394 with no specific unit. 4,570 ClinVar units (4,388 OMIM, 182 MONDO-only).

### a04 — protein mapping
* MANE Select v1.5 by ClinVar gene symbol (fallback HGNC). 59,945 ClinVar Names are on the MANE RefSeq
  transcript with ref aa = MANE protein residue → position used directly; 103 non-MANE → Ensembl VEP
  (64 missense on MANE, 37 not missense on MANE → dropped, 2 no MANE consequence); 8 genes without MANE
  Select. **60,009 mapped**. VEP validation on 500 random MANE-transcript variants: 499 agree (1 = VEP
  returned no MANE consequence for that symbol).
* UniProt: entry whose MANE-Select xref = MANE ENST (3,485 genes) or HGNC xref (113).
* AlphaFold v6: single-fragment for 3,479 genes (GNAS P63092 from the Swiss-Prot tar fetched by
  multigene), **116 fragmented (>2,700 aa) genes excluded** (list in `counts.json →
  protein_mapping.fragmented_genes`; 8 would otherwise be eligible: CHD7, COL7A1, DYNC1H1, FBN1, KMT2D,
  RYR1, USH2A, VWF), 3 without model (EPM2A, GPX4, SELENON — selenoproteins/absent).
* MANE→AF and MANE→UniProt positions: identity when sequences are equal, else s06 `align_map`
  (blocks ≥10 aa, ≥90 % identity). 350 genes have MANE ≠ AF sequence.
* **Reference aa vs AF model: 54,667 / 54,673 position-mapped ClinVar variants match (99.99 %);
  54,667 / 54,706 (99.93 %) of all ClinVar variants in single-fragment genes** (33 fall in unaligned
  segments). Mismatches: KRAS D153 ×3 (MANE = KRAS4B, UniProt/AF canonical = KRAS4A), DUOX2 S25R, DYSF
  V526L, RPGR E585K. Merged table: 55,221 / 55,227 match.

### a05 — phenopacket-store and merge
phenopacket-store v0.1.27: 5,615 missense rows (VEP on MANE Select; all ENSP equal MANE v1.5) →
2,625 distinct protein variants (1,522 with single-allele carriers, i.e. heterozygous without a second
causative allele or hemizygous; 1,108 with biallelic carriers; 5 with both). Carriers and PMIDs are
counted separately for single-allele and biallelic carriers; diseases are the carriers' OMIM ids.
Key = (HGNC, MANE position, ref, alt). ClinVar 60,009 records → 58,881 protein variants (1,090 from
>1 nucleotide change). **Merged: 60,571 variants — clinvar 57,946, ppkt 1,690, both 935.** Where both
sides have labels (778): identical 524, overlap 190, disjoint 64.

### a06 — inheritance (HPO phenotype.hpoa, aspect I)
4,758 units (4,570 ClinVar + 188 ppkt-only). Terms mapped to classes via the HPO is_a hierarchy (AD,
AR, XLR, XLD, XL, MT, YL, SD, COMPLEX; modifiers recorded but unclassed). Clean label when one class
remains (XL+XLR→XLR, XL+XLD→XLD): AR 2,190 · AD 1,810 · XLR 164 · XLD 52 · XL 49 · COMPLEX 12 · MT 3 ·
YL 2 · none 476 (131 AD+AR; 288 units without HPOA inheritance — MONDO-only units have none because
HPOA Orphanet records carry no inheritance). `broad_unit` flags 109 umbrella units (MONDO term with
OMIM-mapped descendants, e.g. OMIM:268000 Retinitis pigmentosa).

### a07 — genes (counting only)
"Countable" variant: single-fragment AF model, MANE position maps onto the model, ref aa matches
(55,221 of 60,571). PRIMARY definition = ClinVar single-label (28,803 countable).
* genes with ≥1 unit with ≥5 variants: 1,009; **eligible (≥2 units each ≥5 distinct single-label
  missense variants): 92 genes**, 214 units, 4,935 variants in those units.
* variants: combined labels (ClinVar ∪ ppkt) 99 genes; OMIM-only units 91; excluding broad units 81.
* **AD/AR-eligible (≥1 AD-only and ≥1 AR-only unit, each ≥5): 13** — COL3A1, COL4A3, COL4A4, BEST1,
  CAPN3, FH, GJA1, LPL, TWNK, EDAR, TSHR, CLCN7, POLR3B (same 13 with combined labels).

### a08 / a09 — AlphaMissense and residue features (99 feature proteins = eligible under any definition)
AlphaMissense: all 19 substitutions × residues for the 99 accessions (1.85 M rows; 93 identical
numbering to the AF model, 6 mapped by alignment, s07 `am_mapping`) plus the exact substitution row for
every other variant (46,291 / 46,899 keys found). Residue features: `structfeat.compute_residue_features`
(+ added neighbourhood pLDDT, Cα-count burial, half-sphere exposure, coordinates, domains, AlphaMissense
per residue) for 97,133 residues; contact-graph edge lists saved. See derived README.

### a11 — sanity checks (`counts.json → sanity`)
* top-20 eligible genes and units, known allelic genes (eligible: LMNA, RET, GJB2, SCN5A, MYH7, COL2A1,
  FGFR2, FGFR1, PRPH2, SCN1A, CACNA1A, COL1A1; not eligible: PTPN11 (57 multi-label vs 32 single-label,
  NS1 + NSML co-listed), FGFR3, KCNQ1, TBX5, KRAS, HRAS, TP53, KCNH2, TTR, NOTCH3, PIK3CA, TGFBR2, GNAS;
  RYR1 fragmented).
* spot-check of 5 random single-label variants against live ClinVar (E-utilities): 5/5 title, protein
  change, classification and OMIM trait agree.
* `common.rwr_profile` vs `mechanistic_mvp.rwr` (repo iterative code): max |Δp| 1.6e-14; region masses
  from the profile equal the stored CoreImpact/DisorderImpact (3e-16).

## Analysis helpers (`common.py`)
```python
import sys; sys.path.insert(0, "/home/user/phenoBridge/allelic"); import common as C
p = C.rwr_profile("P02545", 482, alpha=0.3)   # pandas Series residue -> stationary prob. (exact solve, sums to 1)
pos, xyz, plddt = C.ca_coords("P02545")
pos, W, A = C.contact_graph("P02545")          # W = A D^-1 (D = max(1, degree)), as mechanistic_mvp.rwr
```

## Known limitations / data-quality notes
* ClinVar variant_summary gives the aggregate (VCV) germline classification; condition links are the
  union over RCVs. Because the aggregate is non-conflicting P/LP, contributing submissions are P/LP, but
  a condition can come from a lower-review-status submission.
* MedGen-only conditions (7,951 links) are dropped even when their names are specific; MONDO-only units
  (182) lack inheritance.
* Multi-label variants are common in classic allelic genes (PTPN11, RET, FGFR3, KCNQ1); they are kept and
  flagged, not used for single-label counts.
* KRAS: MANE (4B) ≠ UniProt/AF canonical (4A) at the C-terminus; positions ≥151 do not map reliably.
* Fragmented AF proteins (116 genes incl. RYR1, FBN1, VWF) are excluded.
