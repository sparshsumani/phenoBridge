# Revision notes: variant-recurrence manuscript (2026-11-15)

**Revised from:** `Sumani_Bhatia_variant_recurrence_revised.pdf` (18 pages, uploaded; no source in the repository).
**New editable source:** `body.html` + `style.css`, built by `build.py`.
**Output:** `Sumani_Bhatia_variant_recurrence.pdf` (11 pages including the supplement).

## 1. Summary of changes

### Positioning and title
- **New title:** *Variant recurrence and publication clustering in literature-derived Mendelian disease cohorts: implications for genotype–phenotype inference in rare-disease research.*
  - "Rare-disease research" names the field.
  - The title does not claim that every one of the 44 disorders is formally rare. The Methods and Limitations state that prevalence was not assessed.
- **Central narrative:** literature-derived cohorts are small, recurrent and publication-clustered, so the question is what kind of evidence they supply. The distinction between *patients* and *distinct variants* as units of replication now runs through the abstract, introduction, Table 3 and the conclusions.

### Abstract and front matter
- Structured abstract (Background / Methods / Results / Conclusions), 290 words, with keywords.
- A Key Points box and a separate title page were drafted, then removed when the layout switched to the LaTeX-style two-column format requested on 2026-11-15. The abstract carries the same messages.

### Introduction (rewritten)
- **Why genotype–phenotype inference matters in rare disease:** prognosis and surveillance, counselling, variant interpretation (ACMG/AMP phenotype-specificity evidence) and mechanistic hypotheses.
- **Three properties of literature-derived cohorts:** small size, variant recurrence, and entry through publications. Ascertainment is linked to known differences in penetrance between clinical and population cohorts (Wright 2019).
- **Explicit contrast between two kinds of question:**
  - variant-specific questions, e.g. carriers of *LMNA* p.Arg482Trp, where repeated carriers are informative;
  - relationships across variant properties, where distinct variants are the replication unit.
- **Framing:** pseudoreplication (Hurlbert; Lazic) and the design effect.
- **Stated gap and numbered contributions.** The statement that recurrence alone does not establish dependence is kept.

### Methods
- Reorganised under medical-journal headings.
- **New items:**
  - an ethics statement;
  - an explicit statement that rare-disease prevalence was not assessed;
  - a **Cohort structure** subsection (carriers, distinct variants, publications and largest-publication share per gene), described as added during revision, descriptive and untested;
  - the variant-level scheme framed in standard permutation terms (variant = exchangeability block; Winkler 2014);
  - a short KBG paragraph pointing to the supplement.
- The provenance and deviation list was kept intact, with the new descriptive summaries added to the exploratory list.

### Results
- **Every number, statistical procedure, table and figure from the uploaded version was kept.** An automated comparison of all numbers in the two PDFs found none lost; the only differences were section numbers and citation brackets.
- **New descriptive results** (from `analysis_multigene/cohort_structure.py` → `out/X_cohort_structure*.{json,csv}`):
  - 212 source publications;
  - per gene: median 28.5 carriers (IQR 21–42), 15 distinct variants (8–20) and 3 publications (1–6);
  - the largest publication contributes a median 64% of a gene's carriers, and at least half in 33/44 genes.
- **Headline results are emphasised in the abstract and results subheadings:**
  - 76% recurrence;
  - 4.7% vs 14.9% nominal rejections;
  - 32 → 21 variant-class discoveries.
- **LMNA, SPTAN1 and KBG as practical examples:**
  - *SPTAN1*: spastic paraplegia (driven by p.Arg19Trp) vs motor axonal neuropathy (shared by several variants).
  - *LMNA*: Arg482 lipodystrophy.
  - **KBG** is a new main-text subsection, using existing `kbg_followup/out` results:
    - 333 records citing one review but drawn from 41 sources;
    - 58 individuals in 22 family clusters;
    - 34 carriers of c.1903_1907del;
    - two inconsistent family labels;
    - 20 candidate cross-source duplicate pairs, against ≈14 expected by chance;
    - the one-per-family sensitivity analyses.
- The text says explicitly that the attenuated hotspot associations are not shown to be false.
- **Table 2** now gives individual-level and variant-level counts for all six control genes with tests. *SETD2* and *EHMT1* were filled from `X_controls_perm.csv`.
- **Figure 1C** y-axis relabelled from "Type I error of naive test" to "Null rejection rate, individual-level test", because the old label implied a known ground truth. The data are unchanged.

### Discussion
- New subsections:
  - **Implications for rare-disease registries and phenopacket-based analyses**, including the PMID-as-source problem shown by KBG.
  - **Incorporating these checks into existing workflows**, with Box 1, a reporting checklist explicitly labelled as unvalidated recommendations.
  - **Table 3**, matching the null model to the inferential target, also labelled as reasoning rather than a validated rule.
- The design effect is consistently described as a promising exploratory warning signal, not a validated correction.
- Variant-level permutation is described as one option among several (cluster-robust/mixed models, variant-aggregated tests), not as universally correct.
- **Limitations** keep all six required themes, each under its own heading:
  - exchangeability and ascertainment;
  - publication and family confounding, now including the PMID-proxy weakness;
  - position adjustment and approximate residual permutation;
  - internal assessment of the diagnostic;
  - post hoc decisions and multiplicity;
  - lack of independent validation and limited sensitivity to structural effects.
- The **Conclusions** emphasise reliability and interpretation and make no clinical claims.

### Format
- **LaTeX-like two-column article layout, matching `manuscript_vus3d/`:**
  - centred bold title, then authors and affiliations;
  - indented abstract;
  - numbered sections and subsections, generated by CSS counters;
  - bracketed numeric citations, placed before punctuation;
  - booktabs-style ruled tables;
  - full-width figures and tables placed between column blocks;
  - a compact Vancouver reference list.
- Back matter (Data and code availability, Acknowledgements, Author contributions, Funding, Conflicts of interest) is unnumbered.
- The supplement starts on a new page.
- British spelling is used throughout.

### References
- **Six new references, each verified against Crossref (DOI, title, journal, volume, pages):**
  - Nguengang Wakap 2020 (rare-disease prevalence);
  - Richards 2015 (ACMG/AMP);
  - Wright 2019 (ascertainment and penetrance);
  - Lazic 2010 (pseudoreplication in biomedicine);
  - Winkler 2014 (permutation inference, exchangeability blocks, residual permutation);
  - Breslow & Clayton 1993 (GLMMs).
- **All 35 DOI-bearing existing references re-verified;** Kish 1965 and Donner & Klar 2000 are books.
- Reference 13 of the uploaded version (Donner & Klar) and the other books are kept as before.
- Total: 43 references, numbered automatically in order of first citation.

## 2. Unsupported claims and citations to check before submission

| # | Location | Claim | What is needed |
|---|---|---|---|
| 1 | Introduction ¶1 | Rare diseases affect 3.5–5.9% of the population, and "most have a genetic cause" | The DOI is verified. Confirm in the full text that 3.5–5.9% and ≈72% genetic are Nguengang Wakap's figures. |
| 2 | Introduction ¶3 | Whether a case is published "can depend on its phenotype" | **No citation.** Add evidence on publication/ascertainment bias in case reports, or soften. |
| 3 | Introduction ¶3 | Recurrence arises through founder alleles, hotspots, recurrent *de novo* events and relatives | General knowledge; add a review citation if a reviewer asks. |
| 4 | Introduction ¶3 | Clinically ascertained penetrance and expressivity differ from population estimates (Wright 2019) | The DOI is verified; check that the paraphrase matches the paper's conclusions. |
| 5 | Introduction ¶5 | "To our knowledge … not quantified systematically" | **A novelty claim.** Run a targeted search, including the GPSEA full text (Rekerle 2026) and phenopacket-store papers, for any analysis of recurrence or publication clustering. |
| 6 | Results (class test) | *SPTAN1* p.Arg19Trp and *LMNA* Arg482 "have been described as variant-specific findings" | Check that Van de Vondel 2022 reports p.(Arg19Trp) with spastic paraplegia and that Shackleton 2000 reports R482 in partial lipodystrophy. |
| 7 | Methods (controls) | Control-gene rationales (TBX5 G80R/R237, SETD2 R1740W, IKZF1 N159) | These carry over from the earlier version; confirm each against its cited paper. |
| 8 | Methods (ethics) | "Institutional ethics review was therefore not sought" | **Authors must confirm** this is acceptable to UNSW, Purdue and the target journal; some journals require an explicit exemption statement. |
| 9 | Discussion, Table 3, Box 1 | Reporting recommendations | These are the authors' reasoning, labelled as unvalidated. Consider citing general reporting guidance (e.g. STREGA/STROBE) if a reviewer asks. |
| 10 | Data availability | Release tag and archive | **Missing.** Add a git tag and a Zenodo DOI. |
| 11 | Acknowledgements | HackRare at Harvard University | Authors to confirm the event name, host and year. |
| 12 | Title page | Corresponding author, address, email, CRediT contributions, funding, conflicts of interest | Placeholders. |
| 13 | Methods | Software versions (Python, NumPy, SciPy, pandas, statsmodels, Playwright) | **Not reported.** Add them, or pin them in a lock file. |
| 14 | KBG section | The cohort's single "review" PMID (36446582) | Name the review explicitly (it appears to be the Martinez-Cayuelas 2023 series); confirm. |

## 3. The three largest remaining weaknesses

1. **There is no ground truth for which null is correctly calibrated.**
   - "Three times the nominal rate" is measured *relative to the variant-level null*, which is itself an assumption, and it was chosen after a pilot.
   - We have not shown, under a known generative model, that variant-level permutation controls type I error while individual-level testing does not.
   - We have also not shown how either behaves when the dependence comes from publications or families that span several variants.
   - A reviewer can fairly ask whether the inflation is a property of the data or of the chosen null.
2. **The design-effect diagnostic is assessed in a partly circular way.**
   - The approximation and the variant-level permutation null both depend on the same within-variant phenotype correlation in the same data, so close agreement (73.3 vs 70.6; ρ = 0.92) is expected by construction.
   - It does not show that the diagnostic predicts the true false-positive rate.
   - It was developed and evaluated on one collection, with no external cohort.
3. **Publication and family structure is described but not modelled, and no data are independent.**
   - Publication is proxied by PMID, which can point to an aggregating review (KBG).
   - The planned one-per-family analysis was not done for the multi-gene cohort.
   - No independent replication shows that discoveries surviving the variant-level null replicate better than those that do not.
   - Together with modest novelty against the general pseudoreplication literature, this is the main risk for a methods journal. The paper needs a demonstration that the proposed checks change conclusions *in the right direction*.

## 4. Plan: simulation-based calibration (not performed)

Designed with the ADEMP framework (Morris, White & Crowther, *Stat Med* 2019;38:2074–2102, doi:10.1002/sim.8086, verified). Write and commit the protocol before running.

- **Aims.**
  - Measure type I error and power of each inference scheme under variant, family and publication dependence.
  - Measure the calibration of the design-effect diagnostic.
  - Check that the residual-permutation procedure (position adjustment plus median imputation) behaves as an exact test.
- **Data-generating mechanisms.** Resample the 44 genes' observed structure: carrier counts per variant, the variant–publication incidence and real feature/position vectors. Generate binary phenotypes as
  `logit P(y_i = 1) = β₀ + β·x_v(i) + γ·pos_v(i) + u_v(i) + w_pub(i) + f_fam(i)`,
  where `u ~ N(0, σ_v²)`, `w ~ N(0, σ_p²)`, `f ~ N(0, σ_f²)` and publications may span variants. The grid:

  | Parameter | Values |
  |---|---|
  | σ_v | 0, 0.5, 1, 2 |
  | σ_p | 0, 0.5, 1 |
  | σ_f | 0, 1 (families placed within publications) |
  | Prevalence | 0.1, 0.3, 0.5 |
  | Distinct variants K | 5, 10, 15, 30 |
  | Cluster-size shape | empirical, uniform, one dominant hotspot |
  | β per SD | 0, log 1.5, log 2 |
  | γ | 0 or matched to the observed position effect |

  Add publication-dependent missingness of HPO calls (MCAR vs MNAR by publication).
- **Estimands.**
  - Null: rejection probability at α = 0.05 and 0.005.
  - Non-null: power; the predicted-vs-empirical rejection rate for the design effect; FDR of BH discoveries across a 44-gene panel.
- **Methods compared.**
  - analytic Mann–Whitney;
  - individual-level permutation;
  - variant-level permutation;
  - publication-blocked permutation;
  - a variant-aggregated test (weighted variant-mean phenotype vs feature);
  - a GLMM with variant (and publication) random effects;
  - GEE with cluster-robust SEs by variant and by publication;
  - a DE-adjusted z-test (z/√DE);
  - a cluster bootstrap.

  Each is run with and without position residualisation, and with and without median imputation.
- **Performance measures.**
  - Bias and Monte Carlo SE of rejection rates: ≥2,000 replicates per null scenario gives an MCSE of about 0.5% at 5%.
  - Power curves.
  - Calibration slope and intercept of DE-predicted vs empirical rejection.
  - The lower-tail behaviour seen for S8.
- **Decisive questions.**
  1. With σ_v > 0 and β = 0, is individual-level testing inflated by the amount the DE predicts?
  2. With σ_v = 0 and σ_p > 0, i.e. dependence through publications only, is variant-level permutation still calibrated? If not, publication-blocked or model-based methods are needed.
  3. Is the residual permutation exact, and does median imputation bias it?

## 5. Plan: external validation (not performed)

1. **Temporal holdout within phenopacket-store.** Freeze v0.1.27. Use phenopackets added in later releases, or a split by publication date within genes, as the validation set. Pre-specify:
   - whether recurrence and publication structure replicate;
   - DE-predicted vs permuted rejection agreement in new genes (target calibration slope 0.8–1.2);
   - **the key usefulness test:** do discoveries made in early publications replicate in later publications more often when they survive the variant-level null than when they are significant only analytically? This directly tests whether the stricter null tracks reproducibility.
2. **GPSEA's published associations.** Re-run GPSEA's 253 reported correlations with their own predicates, adding:
   - variant counts per group;
   - largest-variant and largest-publication shares;
   - variant-level permutation;
   - DE.

   Report how many rest on a single recurrent variant. This needs the GPSEA predicate definitions, which should be checked for availability.
3. **Independent case-level resources with family and source information.** Candidates include DECIPHER and gene-specific LOVD databases, subject to access terms (not yet checked). Use them to run the one-per-family analysis and to test whether publication-proxy clustering under-states dependence, as KBG suggests.

## 6. The three highest-priority analyses before journal submission

1. **The simulation study in §4.** It is the single most important addition: it turns "the rejection rate depends on the null" into "this scheme is calibrated under these dependence structures and that one is not". It also gives the design-effect diagnostic a validated range of use.
2. **Temporal/independent replication of discoveries (§5.1).** It is the most direct evidence that the proposed checks improve reproducibility rather than just reducing the number of discoveries.
3. **Model-based and source-aware re-analysis of the real data.** Run the following on the 44-gene cohort, with publication-blocked permutation as a sensitivity analysis:
   - variant-aggregated tests;
   - GLMM/GEE with variant and publication terms;
   - a one-individual-per-family analysis where family labels exist.

   Then report agreement between methods for Analysis A and B discoveries. This answers the expected reviewer request for "an appropriate model" and covers the omitted pre-specified family analysis.

## 7. Submission notes
- The abstract is 290 words. Many genetics journals cap structured abstracts at 250, so trim for the chosen journal.
- The main text is about 5,100 words, plus 3 tables, 2 figures and 1 box; check the journal's limits. Most journals require a double-spaced single-column submission file, so `style.css` would need a submission variant; the two-column PDF suits bioRxiv.
- **Candidate venues** (check current scope and limits):
  - *Genetics in Medicine Open*;
  - *HGG Advances*;
  - *European Journal of Human Genetics*;
  - *Orphanet Journal of Rare Diseases*;
  - *Genetic Epidemiology* (methods).
- For bioRxiv, the current single PDF with the supplement appended is suitable. Most journals will want the supplement as a separate file.
