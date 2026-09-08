# ANKRD11 / KBG syndrome codebase audit

**Scope.** Inventory of every dataset, preprocessing step, patient/variant mapping, phenotype
representation, structural-graph construction step, RWR/mechanistic feature, model and evaluation
procedure currently implemented; a reproducible pipeline map; leakage and statistical-validity
findings; and a ledger of exactly which recorded results can be regenerated from the raw files.

**No scientific methodology was changed.** Everything below is observation and verification. The two
scripts in `audit/` are read-only harnesses that re-run the existing code; they do not modify it.

Audited commit: `a6d20d4` ("Hackathon Code"). Reproduction run: 2026-09-08, Python 3.11.15, versions
pinned in `audit/requirements.txt`.

---

## 1. Repository contents

| File | Size | Role |
|---|---|---|
| `ANKRD11_KBGS_individuals.json` | 1.0 MB | Raw cohort (phenotypes + variants) |
| `AF-Q6UB99-F1-model_v6.pdb` | 1.7 MB | Raw structure (AlphaFold ANKRD11) |
| `ANKRD_mechanistic.ipynb` | 772 KB | Exploratory analysis, 57 code cells, 0 markdown cells |
| `mechanistic_mvp.py` | 23 KB | Gradio demo app; refits the phenotype models at startup |

No `requirements.txt`, no README, no tests, no CI, no seeds file, single commit.

---

## 2. Datasets

### 2.1 `ANKRD11_KBGS_individuals.json`

PheTools / phenopacket-store export. `phetoolsSchemaVersion` 0.2, `hpoVersion` 2025-10-22,
`cohortAcronym` KBGS, `cohortType` mendelian. Curated by ORCID 0000-0002-0736-9199 (2025-01-09,
revised 2025-12-31).

| Key | Contents | Used? |
|---|---|---|
| `diseaseList` | 1 entry — OMIM:148050 KBG syndrome, AD (HP:0000006), ANKRD11 / HGNC:21316 / NM_013275.6 | metadata only |
| `hpoHeaders` | 48 HPO terms; list order defines the phenotype column order | **yes** |
| `rows` | 333 individuals, all from **one** publication (PMID:36446582); exactly 1 allele each | **yes** |
| `hgvsVariants` | 164 variants with hg38 coords, ref/alt, `hgvs`, `gHgvs`, **`pHgvs`** | **never read** |
| `structuralVariants` | 73 entries (label, svType, chromosome) | **never read** |
| `intergenicVariants` | 0 | — |

Per-individual fields: `individualId`, `sex`, `pmid`, `title`, `comment`, `ageOfOnset`,
`ageAtLastEncounter`, `deceased`. Only `sex`, `pmid`, `title` are populated — the other four are the
literal string `"na"` for all 333 individuals, so no age or severity covariate exists.

**Phenotype cell census (15,984 cells = 333 × 48):**

| Value | Count | Share |
|---|---|---|
| `Na` (not assessed) | 10,735 | **67.2 %** |
| `Excluded` (absent) | 2,633 | 16.5 % |
| `Observed` (present) | 2,616 | 16.4 % |

### 2.2 `AF-Q6UB99-F1-model_v6.pdb`

AlphaFold DB model v6 for UniProt Q6UB99 (ANKRD11_HUMAN). Single model, single chain A, residues
1–2663, B-factor = pLDDT. Residue numbering equals UniProt numbering.

pLDDT: min 17.77, **median 31.19**, mean 39.208, max 96.44. Only 327 / 2663 residues (12.3 %) reach
pLDDT ≥ 70. ANKRD11 is predominantly intrinsically disordered, and the model reflects that.

No network calls occur in the executed path. Notebook cell 0 imports `requests`, `Bio.SeqIO`,
`PairwiseAligner`, `cdist` and `mannwhitneyu`; none are used.

---

## 3. Pipeline map

### 3.1 Preprocessing (notebook c2–c3; `mechanistic_mvp.py:load_kbgs_json`)

1. **Phenotype matrix** — `hpoData` is a positional list aligned to `hpoHeaders`; padded with
   `{"type":"Na"}` or truncated to length 48. `Observed→1`, `Excluded→0`, `Na/NA/None→NaN`.
   Result 333 × 49.
2. **Allele long table** — from `alleleCountMap`: 333 rows, 232 unique `alleleKey`s.
3. **Allele one-hot matrix** — `pivot_table(aggfunc="sum", fill_value=0)`, 333 × 232.
4. **Individual table** — 8 columns, deduplicated on `individualId` (no duplicates exist).
5. **Merge** → `df_model`, 333 × 289.
6. **`fillna(0.0)` applied to the whole frame** (notebook c15/c21/c26; `mechanistic_mvp.py:404-408`).
   This is the step that collapses `NaN` phenotypes into `0` — see finding **L-1**.

### 3.2 Patient / variant mapping (notebook c4; `mechanistic_mvp.py:350-369`)

`alleleKey`s are non-standard strings such as `c3590_3594del_ANKRD11_NM_013275v6` and
`ANKRD11_SV_16q24_3_89336307_89354085_x1`. The pipeline regex-parses **the key string**, not the
curated annotations:

| Step | Rule | Result on this data |
|---|---|---|
| `parse_protein_pos` | `p\.?\(?[A-Za-z]{1,3}(\d+)[A-Za-z]{1,3}\)?` | matches **0 / 333** (keys are cDNA-style) |
| `parse_cdna_pos` | `c\.?(\d+)` → first cDNA coordinate | 256 / 333 |
| `cdna_to_aa_pos` | `ceil(c / 3)` | — |
| `aa_pos` | protein `combine_first` cDNA | **256 / 333**, range 71–2612, **153 distinct** |
| `TRUNC_like` | `(stop\|ter\|\*\|nonsense\|frameshift\|fs\|del\|dup\|ins)` (`mechanistic_mvp.py` also matches `splice`; identical result here) | 178 / 333 |
| `MISSENSE_like` | requires a `p.` pattern | **0 / 333 — dead feature** |

The 77 unmapped individuals are exactly the structural-variant carriers.

Notebook-only detour (c19–c22): `variant_cols` is selected by `startswith("c")`, which also captures
the `comment` column (233 columns for 232 alleles); the 333-row long table is then merged onto a
**non-unique** `alleleKey → aa_pos` map, inflating it to 1,751 rows (1,674 after `dropna`) — a 5.3×
duplication. Aggregation is `.max()`, so fitted values are unaffected, but c22's printed counts
describe duplicated rows rather than patients.

Patient-level aggregation everywhere: `groupby("individualId").max()` over the impact features.

### 3.3 Phenotype representation

48 binary HPO columns. **No HPO ontology is loaded** — no propagation to ancestor terms, no
information-content weighting, no semantic similarity, no term clustering. The ternary
Observed/Excluded/Na source is collapsed to binary by `fillna(0)`. One derived outcome:
`NeuroScore` = `HP:0001249 + HP:0001263 + HP:0000750` (range 0–3, notebook c46).

### 3.4 Structural graph construction (notebook c7–c11, c23–c24; `mechanistic_mvp.py:53-170`)

1. Parse chain A Cα atoms, skipping heteroatoms and residues without CA → 2663 residues.
2. `build_confidence_graph(cutoff_A=8.0, min_plddt_node)` — nodes filtered by pLDDT; edges from
   `cKDTree.query_pairs(r=8.0)` on Cα coordinates; edge attribute `dist`.
   - **Node identity differs between the two codebases**: the notebook uses the original residue
     index 0…2662, `mechanistic_mvp.py` uses the UniProt position. Both are internally consistent
     (each maps via `unip_pos`), but the two are off by one and are not interchangeable.
   - `G_all` = **2663 nodes / 6911 edges** (mean degree 5.19).
   - `G_core` (pLDDT ≥ 70) = **327 nodes / 1354 edges**.
3. **Region sets — three mutually incompatible definitions appear:**

| Where | core | semi-core | disorder |
|---|---|---|---|
| notebook c10 | ≥ 70 (327) | — | ≤ 50 (2232) |
| notebook c16 | ≥ q95 = 90.994 (134) | — | ≤ q50 = 31.19 (1332) |
| notebook c23 **(final)** | ≥ q90 = 81.19 (268) | all − disorder (**1331, includes core**) | ≤ q50 (1332) |
| `mechanistic_mvp.py` | ≥ q90 (268) | ≥ q50 (1334); scored region is semi − core (**1066**) | ≤ q50 (1332) |

4. `shortest_path_backbone(G, sources, targets, max_pairs, top_k, seed=0)` — all N-anchor
   (pos ≤ 200) × C-anchor (pos ≥ 2463) pairs, subsampled with `random.Random(0)`, unweighted BFS via
   `nx.shortest_path`, count interior-node visits, take `top_k`:

| Where | graph | max_pairs | top_k | feeds anything? |
|---|---|---|---|---|
| notebook c11 | `G_core` | 2000 | 200 | c13/c14 only |
| notebook c24 | subgraph(core q90) | 1500 | 200 | **yes — this is the `logBackbone` used in c49/c53** |
| notebook c27 | subgraph(semicore) | 2000 | 300 | **no — inert** |
| notebook c30 | subgraph(semicore) | 3000 | 600 | **no — inert** |
| `mechanistic_mvp.py` | `G_core` (pLDDT ≥ 70) | 2000 | 200 | yes |

### 3.5 RWR and mechanistic features

- `rwr(G, seed_node, alpha=0.3, max_iter=200, tol=1e-10)`. The notebook uses the pull form, the app
  the push form plus a final L1 renormalisation; they are mathematically equivalent here (no
  isolated nodes, so no mass leaks and the renormalisation is a no-op).
- `rwr_conf_weighted` (pLDDT-weighted transitions, notebook c28) — **defined, never called**.
- **Fingerprint** = region masses of the stationary RWR vector seeded at the variant residue:
  `CoreImpact`, `SemiCoreImpact`, `DisorderImpact`, `CommBackboneImpact`. The app additionally
  returns `aa_pos`, `seed_pLDDT`, `alpha`, `cutoff_A`.
- `logCore` / `logBackbone` / `logSemiCore` = `log10(x + 1e-15)`. **Zeros map to exactly −15.**
- **`TruncFrac` names two different quantities:**
  - notebook c38: `min(aa_pos per patient) / 2663` — a normalised sequence position; unmatched → 1.0
  - `mechanistic_mvp.py:400`: `mean(TRUNC_like per patient)` — a truncating-variant fraction

### 3.6 Models

| ID | Specification | Where |
|---|---|---|
| M1 | `Logit(HPO ~ 1 + logSemiCore)`, top-5 phenotypes by prevalence | c42 |
| M2 | `Logit(HPO ~ 1 + logSemiCore + TruncFrac)`, top-5, then all 31 valid + BH-FDR | c43, c45 |
| M3 | `OLS(NeuroScore ~ logSemiCore + TruncFrac)` | c48 |
| M4 | `Logit(HP:0001249 ~ logSemiCore + logBackbone + TruncFrac)` | c49 |
| M5 | `DecisionTreeClassifier(max_depth=3, min_samples_leaf=15, random_state=42)` on the same 3 features → HP:0001249, single 70/30 split, `random_state=42` | c53 |
| M6 | App: `Logit(HPO ~ 1 + logSemiCore + TruncFrac)` for every term with `0.05 ≤ prevalence ≤ 0.95`, fit on **all 333**; coefficients cached and applied to user input via the sigmoid; ✅ badge at raw *p* < 0.05 | `mechanistic_mvp.py:372-483` |

Phenotype inclusion filters: notebook `15 ≤ n_pos ≤ n − 15` → 31 of 48; app `min_prev` 0.05 → 33 of 48.

### 3.7 Evaluation

Wald *p* and odds ratio on `logSemiCore` (in-sample); BH-FDR **once**, in c45 only; pseudo-R² and
LLR *p* (c49); `classification_report` on one 30 % holdout (c53); feature importances (c54);
`plot_tree` (c55); prevalence tables, a `corrcoef` check (c31) and zero-fraction checks (c32).

**Absent throughout:** cross-validation, ROC/PR/AUC, calibration, bootstrap CIs, permutation tests,
group-aware splitting, external or held-out validation, and any multiplicity control in the
delivered app.

---

## 4. Reproducibility ledger

Verified by `python audit/reproduce_notebook.py`, which re-executes the notebook's logic
top-to-bottom against the two raw files.

### 4.1 Reproduces exactly from raw data

Every one of these matched the recorded output to the precision printed in the notebook:

| Cell | Result |
|---|---|
| c2 | 333 rows, 48 HPO headers, 1 disease; phenotype matrix (333, 49); model table (333, 289) |
| c3 | `df_mut` (333, 3); 333 patients; 232 unique allele keys |
| c4 | (333, 9); aa_pos present 256 / 333; range 71.0–2612.0 |
| c6 | aa_pos min/max 71 / 2612; hotspot bins 533–666: 45, 799–932: 28, 2396–2529: 25, 1065–1198: 21, 400–533: 18 |
| c8 | 2663 Cα; positions 1–2663; pLDDT 17.77 / 31.19 / 39.20781449493053 / 96.44 |
| c9 | `G_all` 2663 / 6911; `G_core` 327 / 1354 |
| c10 | structured ≥ 70: 327; disordered ≤ 50: 2232 |
| c11 | backbone size 200 |
| c13 | `fp(100)` = Core 0.0005393824991717977, Disorder 0.9986534862326728, Backbone 0.0005393824991686175 |
| c14 | per-variant impacts for KBG2/3/4/5/64 |
| c16 | q50 31.19, q95 90.994; core 134; disorder 1332 |
| c19 | 333 patient-variant rows |
| c20 | mapped fraction 0.9560251284980011 |
| c21 | `mech_patient` values (KBG1–KBG5) |
| c22 | 1674 / 1674 / 153 unique aa_pos / 256 individuals / 333 total |
| c23 | q90 81.19; core 268; semi-core 1331; disorder 1332 |
| c24 | `G_coreQ` 268 / 1135; backbone 200 |
| c26 | `df_model3` impact and log columns |
| c27, c30 | backbone 300; backbone 600 |
| c29 | `.describe()` over the four impact columns |
| c31 | corr(SemiCore, Disorder) = −0.2980995913438425 |
| c32 | backbone zeros = core zeros = 0.23723723723723725 |
| c33, c36 | mechanistic cohort size 256 |
| c34, c37 | 48 phenotype columns; 31 valid |
| c39–c41 | prevalence tables and the top-10 listing |
| c42 | all five coef / OR / *p* (HP:0001249 coef 0.669489, OR 1.953239, *p* 0.029525) |
| c43 | all five adjusted coef / OR / *p* (HP:0001249 *p* 0.076536) |
| c45 | full 31-row table including `p_fdr` (minimum 0.639399) |
| c48 | OLS R² 0.003, F 0.3699, all coefficients |
| c49 | n 256, pseudo-R² 0.03086, LLR *p* 0.014661; logSemiCore 0.662383 (*p* 0.041328), logBackbone −0.064033 (*p* 0.043797), TruncFrac 0.658128 (*p* 0.209647) |
| c50 | `exp(params)` |
| c53 | `classification_report` identical (0.45/0.37/0.41 @ 27; 0.69/0.76/0.72 @ 50; accuracy 0.62) |
| c54 | importances 0.566427 / 0.310297 / 0.123276 |

Only three edits are needed to run it: the two hard-coded Drive paths (c2, c7) and removing the
`drive.mount` cell.

### 4.2 Cannot be reproduced

| Cell | Recorded output | Why not |
|---|---|---|
| c5 | "9J0A motif window: 345 to 368; Variants in window: 0" | **`df_pos` is never defined in the notebook.** It came from session state that no longer exists. (The result is a dead end regardless — zero variants fall in the window.) |
| c51–c52 | "Low structural impact probability: 0.6279988986486627 / High: 0.7024763258149385" | **`low` and `high` are never defined.** The scenario probabilities cannot be regenerated. |
| c0, c1 | pip install; Drive mount | Environment-specific (Colab) |
| c44, c55 | boxplot; `plot_tree` figure | Regenerable in principle; not byte-comparable |

Further obstacles to a clean re-run: all 57 cells carry `execution_count: null`, so the **actual**
execution order cannot be verified from the file; there are no markdown cells documenting intent;
and no dependency versions are pinned anywhere in the repository.

### 4.3 Results that exist only in the app

`mechanistic_mvp.py` refits its own models at startup and they appear in no notebook cell. Confirmed
by running its `fit_hpo_models_from_cohort` directly: **33 models, 7 at raw *p* < 0.05**, led by
HP:0001155 *p* = 0.000444, HP:0000343 *p* = 0.001846, HP:0001249 *p* = 0.001911. These are the
numbers the demo badges with ✅. They do not correspond to any notebook result — see **R-2** and
**L-2**.

---

## 5. Findings

Severity: **Critical** = invalidates a reported result; **High** = materially biases or blocks
reproduction; **Medium/Low** = correctness and hygiene.

### L-1 · Critical · Unassessed phenotypes are silently recoded as "absent"

`fillna(0.0)` is applied to the entire merged frame (notebook c15/c21/c26;
`mechanistic_mvp.py:404-408`), so every `Na` HPO cell becomes `0`. **67.2 % of all phenotype cells
are `Na`.** The outcome variable in every model is therefore roughly two-thirds fabricated
negatives:

| HPO term | Prevalence among **assessed** | Prevalence after `Na → 0` |
|---|---|---|
| HP:0001263 Global developmental delay | 95.7 % | 52.9 % |
| HP:0001249 Intellectual disability | 87.0 % | 58.3 % |
| HP:0000750 Delayed speech and language | 86.0 % | 48.0 % |
| HP:0000729 Autistic behavior | 58.3 % | 16.8 % |
| HP:0007018 ADHD | 74.4 % | 18.3 % |

This also silently drives the phenotype-inclusion filters: at its true 95.7 % prevalence,
Global developmental delay would fail the app's `min_prev` bound, but at the coerced 52.9 % it
passes. `Excluded` (a genuine negative) and `Na` (unknown) are indistinguishable downstream even
though the raw file separates them cleanly.

### L-2 · Critical · Structural-variant carriers are pinned at the predictor's extreme *and* are the most under-phenotyped

The 77 SV carriers get no `aa_pos`, so `fillna(0.0)` sets `SemiCoreImpact = 0` and
`logSemiCore = −15` **exactly** — a point mass 12.6 units below the entire observed range for mapped
patients (−2.404 … −0.000). Those same 77 individuals have significantly more `Na` cells than the
256 mapped ones (mean 34.9 vs 31.4 of 48; Mann–Whitney *p* = 0.0038). SV status therefore drives the
predictor to one extreme *and* the coerced outcome toward 0 — a backdoor path that manufactures a
positive `logSemiCore` ↔ phenotype association out of ascertainment.

Verified with `python audit/diagnose_divergence.py`, changing nothing but the row filter:

| `fit_hpo_models_from_cohort` run on | models | raw *p* < 0.05 | BH-FDR *q* < 0.05 | min *q* |
|---|---|---|---|---|
| **All 333 (as shipped)** | 33 | **7** | **4** | 0.0147 |
| **The 256 aa-mapped patients** | 34 | 2 | **0** | 0.659 |

The app's entire significance signal disappears when the 77 pinned patients are removed. Note also
that this assigns whole-gene deletions the *lowest* structural-impact score in the model, inverting
the expected severity ordering.

### S-1 · Critical · The shipped app reports uncorrected *p*-values as evidence

`mechanistic_mvp.py:314` is commented "p only, no q"; `SIG_P_THRESHOLD = 0.05` badges ✅ across 33
simultaneously fitted models with no multiplicity control. The notebook's own c45 *does* apply
BH-FDR over the same family and finds **nothing survives** (minimum *q* = 0.639). The correction was
implemented and then dropped from the delivered artifact.

### S-2 · High · Every reported *p*-value is in-sample; fit statistics are presented as validation

`predict_top10_sig_table` scores user input with coefficients fit on the whole cohort and displays,
in the adjacent column, the *p*-value from that same fit. Outside c53's single tree split there is
no held-out estimate of anything.

### S-3 · High · Model and definition search over one dataset, reported without correction

Region thresholds were tried at q95/q50, then q90/q50, then semi-core; the backbone at `top_k`
200 → 300 → 600 and `max_pairs` 2000 → 1500 → 3000; covariate sets at `{logSemiCore}`,
`{+TruncFrac}`, `{+logBackbone}`; outcomes at top-5, all-31, `NeuroScore`, and a single term. The
headline (c49: logSemiCore *p* = 0.041, logBackbone *p* = 0.044) is the last configuration of a
search whose earlier steps sit in the same notebook, reported with no accounting for that search.
Pseudo-R² is 0.031.

### S-4 · High · Non-independent observations treated as i.i.d.

124 of the 256 mapped patients share a variant with at least one other patient; `c.1903_1907del`
alone covers **34 patients (13 % of the mechanistic cohort)** with an identical feature vector.
Sibling pairs are present (KBG8A/KBG8B and KBG10A/KBG10B share alleles). No clustering, random
effect, or robust standard error is used anywhere.

### L-3 · High · The decision tree's split cuts through identical feature vectors

c53 splits 256 patients 70/30 with no grouping. With only 153 distinct `aa_pos` and one variant
covering 34 patients, the same feature vector appears on both sides of the split, so the tree can
memorise a variant in training and be scored on it in test. Separately, the reported accuracy of
**0.62 is below the majority-class baseline of 0.649** on that holdout — which the notebook does not
state.

### D-1 · High · The curated protein consequences in the file are never read

`hgvsVariants` carries `pHgvs` for all 164 SNV/indel variants — e.g.
`NP_037407.4:p.(Arg1986IlefsTer45)`, `p.(Arg174Ter)`, `p.(Arg2512Gln)` — plus hg38 coordinates,
ref/alt and `gHgvs`. The pipeline ignores all of it and regex-parses the variant *key string*
instead. Consequences:

- `MISSENSE_like` is 0 for all 333 patients because the keys never contain `p.`, while the file
  identifies **14** non-truncating variants.
- `TRUNC_like` misses all **56 stop-gain** variants — their keys read `c520CtoT`, containing no
  `del`/`dup`/`ins`/`fs`/`ter` token — so the covariate labels the clearest truncating class as
  non-truncating. (Ground truth from `pHgvs`: 94 frameshift, 56 stop-gain, 14 other.)
- Off-by-one on some deletions: `c7083del` → `ceil(7083/3)` = 2361 vs the curated
  `p.(Thr2362ProfsTer39)`.
- For the 94 frameshifts, `ceil(c/3)` is the frameshift *start*, not the premature stop, so the RWR
  seed sits at the wrong residue.

The 73 `structuralVariants` entries are likewise unread; their genomic intervals would at minimum
allow the deleted residue range to be modelled instead of dropped.

### R-1 · High · The notebook cannot be run as given

- c1 requires `google.colab.drive.mount`; c2 and c7 hard-code `/content/drive/MyDrive/HackRare/…`.
- c5 references `df_pos`, never defined. c51 references `low` and `high`, never defined.
- 57 code cells, **0 markdown cells**, all `execution_count: null` — execution order is unverifiable.
- No pinned dependencies; no seed except `train_test_split` and the tree.

### R-2 · Medium · The app does not reproduce the notebook

The two compute differently-defined quantities under the same names, on different cohorts:

| | notebook | `mechanistic_mvp.py` |
|---|---|---|
| `SemiCoreImpact` region | semi-core, 1331 nodes, **includes** core | semi − core, **1066** nodes |
| Backbone graph | subgraph(core q90), top 200 | `G_core` (pLDDT ≥ 70), top 200 |
| `TruncFrac` | `min(aa_pos)/2663` | `mean(TRUNC_like)` |
| Regression cohort | 256 aa-mapped | **all 333** |
| Multiplicity control | BH-FDR applied | **none** |
| Node identity | residue index 0…2662 | UniProt position 1…2663 |

Net effect: 0 phenotypes survive FDR in the notebook; the app badges 7. The `SemiCoreImpact`
definitions differ by up to 2.3 % on tested seeds (aa 2400: 0.972058 vs 0.993948).

### R-3 · Medium · Dead and inert code

`rwr_conf_weighted` (c28) is defined and never called. c27 and c30 recompute `comm_backbone` at
`top_k` 300 and 600 *after* `logBackbone` was frozen in c26, so they change no downstream result.
c0's `requests` / `SeqIO` / `PairwiseAligner` / `cdist` / `mannwhitneyu` imports are unused. c12/c17/
c25 redefine `score_influence` and `fingerprint_for_variant`; c16/c23 redefine the region sets;
c13/c15/c21 are superseded by c26.

### D-2 · Medium · Silent 5.3× row inflation in the notebook variant merge

See §3.2. `variant_cols` picks up the `comment` column via `startswith("c")`, and the merge onto a
non-unique `alleleKey → aa_pos` map turns 333 rows into 1,751. Fitted values survive because
aggregation is `.max()`, but c22's diagnostics count duplicated rows, not patients.

### D-3 · Medium · The structural substrate is mostly low-confidence

Median pLDDT is 31.19 and only 12.3 % of residues reach 70. The "semi-core" cut sits at the median —
pLDDT 31.19, inside AlphaFold's *very low* confidence band, where backbone coordinates carry no
reliable three-dimensional information. The 8 Å Cα contact graph over those residues, and every RWR
mass computed across them, largely reflects the geometry of a prediction AlphaFold itself flags as
unreliable. Consistent with this, `G_all` has mean degree 5.19 over a 2663-residue chain, so RWR is
close to one-dimensional diffusion along the sequence — which makes `logSemiCore` substantially a
proxy for sequence position. *Flagged only; no methodology change made.*

### S-5 · Medium · Predictor and covariate are both proxies for sequence position

In the notebook `TruncFrac` is literally `min(aa_pos)/2663`, and per D-3 `logSemiCore` is largely
positional too, so M2 and M4 regress phenotype on two transformations of the same underlying
quantity. Adding `TruncFrac` moves HP:0001249 from *p* = 0.0295 (c42) to *p* = 0.0765 (c43),
consistent with collinearity rather than confounding control.

### S-6 · Low · The app's region sets are not a partition

`region_nodes(plddt_lo=q_semi)` keeps ≥ q50 and `plddt_hi=q_semi` keeps ≤ q50, so the 3 residues at
exactly pLDDT 31.19 belong to both `semi_nodes` and `disorder_nodes`. Immaterial to the results, but
the partition is not one.

### S-7 · Low · Single-source cohort with no adjustable covariates

All 333 individuals come from PMID:36446582. `ageOfOnset`, `ageAtLastEncounter`, `deceased` and
`comment` are `"na"` for everyone, so no age or severity adjustment is possible; `sex` is available
but unused. Ascertainment is publication-driven and uniform, which limits generalisation and makes
the missingness pattern (L-1, L-2) a property of the reporting, not of the patients.

---

## 6. What currently stands

- **The structural half of the pipeline is fully reproducible.** Graph construction, region
  definitions, backbone extraction and the RWR fingerprints regenerate bit-for-bit from the two raw
  files.
- **The notebook's own statistical conclusion is the defensible one:** after BH-FDR across the 31
  tested phenotypes, no association between `logSemiCore` and any phenotype survives (min *q* =
  0.639); the tree performs below its majority-class baseline; the OLS on `NeuroScore` has R² = 0.003.
- **The demo app's headline results do not stand.** Its 7 badged phenotypes are an artifact of L-1
  and L-2 combined and vanish under the cohort restriction shown in §5 L-2.
- **Two recorded results cannot be regenerated at all** (c5, c51–c52) because they depend on
  variables that exist nowhere in the repository.

---

## 7. Running the verification

```bash
python -m venv venv && ./venv/bin/pip install -r audit/requirements.txt
./venv/bin/python audit/reproduce_notebook.py     # re-runs the notebook logic linearly
./venv/bin/python audit/diagnose_divergence.py    # isolates the app/notebook divergence
```

Both scripts are read-only, take no arguments, resolve the raw files relative to the repository
root, and modify no pipeline code.
