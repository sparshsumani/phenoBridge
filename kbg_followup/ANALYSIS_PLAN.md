# KBG follow-up: pre-specified analysis plan

Written 2026-11-03, **before** any genotype–phenotype association was computed in this follow-up
(`04_analysis.py` had not been written yet).

**What I had already seen when writing it** (disclosed so the reader can judge how pre-specified this is):
- the draft paper (`paper/paper.html`) and its results: nominal ADHD/triangular-face hits under NA=0 coding,
  none under explicit coding;
- the abstract of Martinez-Cayuelas et al. 2023 (J Med Genet 60:644, PMID 36446582), via Europe PMC;
- marginal counts of each HPO term (observed/excluded) in all 333 individuals;
- covariate-only summaries: variant classes from VEP (174 frameshift, 66 nonsense, 7 splice, 7 missense,
  2 in-frame, 77 SV); PTC positions and NMD rule results; seed pLDDT; Spearman ρ(log10 SemiCore, position) = 0.02;
- duplicate/relatedness clusters (`03_duplicates.py`). Its rules use only genotype, sex and phenotype concordance,
  not any genotype–phenotype association.

## Data and coding
- Individuals: 333 (`ANKRD11_KBGS_individuals.json`). Each is matched 1:1 to its phenopacket; all 333 agree on
  allele, sex and all 48 HPO calls.
- **Phenotype coding**: explicit annotations only. Observed = 1, Excluded = 0, not recorded = missing, and missing
  individuals are dropped term by term. The draft's NA→0 coding ("Encoding A") is used only to reproduce the draft.
- **Variant class** (VEP on ENST00000301030.10 = NM_013275.6, MANE Select). Truncating = frameshift, nonsense or
  canonical splice. Non-truncating point = missense or in-frame indel. SV = copy-number variant
  (70 deletions, 6 duplications, 1 with no copy number given).
- **Position**: protein position from VEP. Splice variants use ⌈c/3⌉ of the exonic anchor. Relative position =
  aa/2663.
- **Fingerprint**: the repo's code is used unchanged (α = 0.3, r = 8 Å). The exposure is log10 SemiCoreImpact,
  standardised (z) within the analysis population. ORs are reported per SD, and per 10-fold for comparison
  with the draft.

## Term eligibility
A term is analysed in a given population if it has ≥10 observed **and** ≥10 excluded among individuals with an
explicit annotation in that population.

## 5a. Primary analysis: does the fingerprint add to the baselines?
- Population: the 256 individuals with a point variant (the fingerprint is undefined for SVs).
- M0: logit P(y=1) = β0 + β1·Trunc + β2·RelPos.  M1: M0 + β3·z(log10 SemiCore).
- Per-term test: likelihood-ratio test of M1 vs M0 (1 df). Benjamini–Hochberg across the eligible terms;
  significance at q < 0.05.
- **Omnibus test**: T = Σ over eligible terms of the LR χ². Null: permute z across individuals within
  variant-class strata (truncating / non-truncating). The same permutation is applied to every term at once,
  which preserves the correlation between terms. B = 5,000 permutations; p = (1 + #{T* ≥ T}) / (B + 1).
- If the ML fit fails to converge or shows separation for a term, the Firth fit is used for that term, and
  this is reported.
- Pre-specified sensitivity analyses (none can rescue the primary result; they are reported regardless):
  - S1: Firth-penalised logistic regression with a penalised LRT (Heinze & Schemper) for every term.
  - S2: relative position modelled flexibly (cubic B-spline, 4 df) instead of linearly.
  - S3: other fingerprint components in place of SemiCore, each separately with BH within component:
    log10 CoreImpact, CommBackboneImpact, DisorderImpact, seed pLDDT. These are secondary.
  - S4: one individual per family cluster, and one per family-or-likely-duplicate cluster (`clusters.csv`).
  - S5: truncating variants only (drops the 9 non-truncating), with M0 = RelPos.
  - S6: the draft's model reproduced with the corrected class flag, under both phenotype codings.
    This is for comparison only.

## 5b. Baseline associations (secondary; BH within each family of tests)
Firth logistic ORs with profile penalised-likelihood 95% CIs and penalised-LRT p-values. Fisher's exact p is also
given for binary exposures.
- B1. SV vs point variant (all 333; eligible terms counted in all 333). Sensitivity: SV deletions only
  (excluding duplications and the SV without a copy number).
- B2. Missense/in-frame vs truncating (point variants). n = 9 non-truncating, so this is expected to be
  underpowered; I report it anyway.
- B3. Truncation position among the 240 PTC variants (frameshift + nonsense):
  - (i) PTC relative position, continuous, OR per 10% of protein length;
  - (ii) PTC that removes the UniProt C-terminal "important for protein degradation" region (PTC ≤ aa 2369)
    vs later PTC.
  - Exon 9 vs outside exon 9 among point variants is analysed in B4b.
- B4. Reproduction of the specific claims in the Martinez-Cayuelas 2023 abstract. These are the same individuals
  (the cohort is curated from that paper), so this checks the data rather than replicating the findings:
  - a. triangular face, sequence variant vs CNV (reported 71.1% vs 42.5%);
  - b. short stature, exon 9 vs outside exon 9 (reported 62.5% vs 27.8%), among point variants;
  - c. ID/ADHD/ASD lower with c.1903_1907del (reported 70.4% vs 89.4%), among point variants.
    Composite = 1 if any of HP:0001249 / HP:0007018 / HP:0000729 is observed; 0 if none is observed and
    ≥1 is excluded. ID alone is also reported.
  - d. macrodontia and comorbidities ↑ with deletion size; hand anomalies ↑ with smaller deletions (SV deletions
    only; OR per doubling of size). "Comorbidities" is operationalised as hearing impairment (any of the 3
    hearing terms), recurrent otitis media, cryptorchidism, any listed heart defect, feeding difficulties or any
    seizure term. These coordinates mix genome assemblies, but size does not depend on assembly.

## 5c. Missense-only subset
If n < 15 (expected: 7 missense), I report counts and a descriptive table only (position, pLDDT, SemiCore,
AlphaMissense) and run no tests.

## 5d. Power analysis (simulation)
- For each eligible primary term, keep that term's real covariate rows (Trunc, RelPos, z). The baseline linear
  predictor comes from M0 fitted to the observed data.
- Simulate y ~ Bernoulli(expit(lp0 + β·(z − z̄))) for OR = e^β on a grid from 1.0 to 4.0 per SD.
  Refit M0/M1 and run the LRT.
- 1,000 simulations per grid point per term. Power is the fraction with p < 0.05 and p < 0.005.
- The minimum detectable OR (80% power) is found by linear interpolation. I report it per term and the median
  across terms, and also per 10-fold change in SemiCore.

## 5e. Sparse-event sensitivity
This is S1 (Firth) for the primary analysis. Firth is also the default for 5b.

## Multiplicity and interpretation
- Confirmatory claims come only from 5a: the per-term BH-adjusted LRT and the omnibus permutation p.
- Everything else is secondary and descriptive. Anything not listed here is labelled **exploratory** in RESULTS.md.
- Random seed: 20261103.
