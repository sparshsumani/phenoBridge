# VUS-3D: plan (written before any data was assembled or any score computed)

**Question.** For a missense variant of uncertain significance, do the *known* pathogenic
and benign variants near it on the AlphaFold structure predict how it will eventually be
classified, beyond AlphaMissense and beyond a gene-level baseline? Does PhenoBridge's
network propagation (RWR on the contact graph) do better than simpler 3D and 1D
neighbourhoods?

This plan was committed before the dataset existed. Anything not listed here is reported
as exploratory.

## Lessons carried over (`manuscript/`, `allelic/RESULTS.md`)
1. Use a question where structure has a strong prior: pathogenic variants cluster in 3D.
2. Use honest baselines. Any gain must exceed AlphaMissense **plus a gene-level prior**,
   so that a model cannot win simply by learning which genes are dangerous.
3. Use leakage-free evaluation. The temporal holdout means no test label is ever used to
   build a score.
4. Genes are the resampling unit for uncertainty.

## Data (built by `vus3d/v0*_*.py`)
- **TRAIN:** ClinVar missense SNVs (GRCh38, ≥1 star) classified P/LP or B/LB in the
  **January 2023** release.
- **TEST:** missense variants that in January 2023 were VUS, conflicting or absent, and
  that are P/LP or B/LB (≥1 star) in the **November 2026** release.
- Variants whose classification flipped between P and B are excluded.
- Only proteins with a single-fragment AlphaFold v6 model are used.

## Scores
All scores for a query residue q use **TRAIN labels only**. Labels at q itself are
excluded (same-residue evidence is ACMG PS1/PM5, not a hotspot).
- **Residue labels.** A residue is P if every TRAIN label at it is P, and B if every
  label is B. Residues with mixed labels are not used as label sources.
- **Label mass.** For a kernel w_q(j), m_P(q) = Σ_{j∈P} w_q(j) and
  m_B(q) = Σ_{j∈B} w_q(j). Model features are log(m_P + ε) and log(m_B + ε), with
  ε = 10⁻⁴ (after normalising w so that Σ_j w_q(j) = 1 over all residues).
- **Kernels:**
  - **RWR** (PhenoBridge): w_q = the exact RWR stationary vector seeded at q on the 8 Å
    Cα contact graph, α ∈ {0.15, 0.3, 0.5};
  - **3D radius:** w_q(j) ∝ exp(−d_qj² / 2σ²), with d the Cα–Cα distance and
    σ ∈ {4, 6, 8, 12} Å;
  - **1D window:** w_q(j) ∝ exp(−|q − j| / λ), λ ∈ {2, 5, 10, 20} residues.
- **Tuning.** Each kernel's parameter is chosen **on TRAIN only**, by 5-fold cross-
  validation with folds grouped by gene. Features for TRAIN variants are computed leaving
  out the variant's own residue.
- **AlphaMissense:** logit of the exact-substitution score.
- **Gene prior:** the smoothed fraction of P among the gene's TRAIN labelled residues,
  (n_P + 1)/(n_P + n_B + 2), on the logit scale.

## Models
All models are logistic regressions fit on TRAIN and evaluated once on TEST.
- **M0:** AlphaMissense.
- **M1:** AlphaMissense + gene prior. This is the baseline.
- **M2-1D:** M1 + 1D-window label masses.
- **M2-3D:** M1 + 3D-radius label masses.
- **M2-RWR:** M1 + RWR label masses.

## Primary hypotheses
Each is tested on TEST variants in genes with ≥5 TRAIN labelled residues. The metric is
pooled AUROC. Uncertainty comes from a paired bootstrap over genes (2,000 resamples);
p is two-sided.
- **H1 (3D network evidence adds):** AUROC(M2-RWR) − AUROC(M1) > 0.
- **H2 (network vs 3D distance):** AUROC(M2-RWR) − AUROC(M2-3D) ≠ 0.
- **H3 (3D vs sequence):** AUROC(best of M2-RWR and M2-3D, chosen on TRAIN CV) −
  AUROC(M2-1D) > 0.

Holm correction is applied across H1–H3.

**What counts as a positive finding:** H1 significant after Holm correction *and*
ΔAUROC ≥ 0.01. A significant but smaller effect is reported as statistically detectable
but practically small.

## Secondary analyses
- **Within-gene AUROC** of each model: the mean over genes with ≥3 P and ≥3 B TEST
  variants. This answers whether the score ranks variants *within* a gene.
- **Gain by data density:** ΔAUROC(M2-RWR − M1) by TRAIN labelled-residue bins
  (5–9, 10–19, 20–49, ≥50).
- **Clinical calibration as hotspot evidence (PM1-like):**
  - From TRAIN only, estimate the likelihood ratio for P versus B in bins of the
    *propagation component* of M2-RWR, i.e. its log-odds minus M1's log-odds.
  - Map these LRs to ClinGen SVI strengths: Supporting LR ≥ 2.08, Moderate ≥ 4.33,
    Strong ≥ 18.7; the agent will verify the thresholds against the literature.
  - On TEST, report the observed LR in each bin and the fraction of P and B variants
    assigned each strength.
- **Precision–recall** and Brier score for M0–M2.

## Sensitivity analyses
- Include same-residue TRAIN labels.
- TEST restricted to variants with ≥2 stars in 2026.
- TEST restricted to variants that were VUS in 2023, excluding absent and conflicting.
- Genes with ≥20 TRAIN labelled residues.
- Fixed pre-chosen kernel parameters (α = 0.3, σ = 6 Å, λ = 5), without tuning.
- Pooled AUROC weighted so that each gene contributes equally.

## Known biases, stated in advance
- **Some TEST classifications may use AlphaMissense** (as PP3 evidence after its 2023
  release). This favours the AlphaMissense baselines and makes the test conservative for
  the 3D score.
- **Some may use hotspot or PM1 arguments based on earlier variants.** This favours
  label-propagation scores. That is the intended clinical use, but it means part of the
  gain could reflect how curators reason rather than biology. We report it as a
  limitation, and the within-gene and same-residue-excluded analyses reduce it.
