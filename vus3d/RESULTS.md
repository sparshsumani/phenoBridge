# VUS-3D results

All numbers are from `vus3d/out/`, produced by `run_vus3d.py --boot 2000` (pre-registered analyses)
and `explore_vus3d.py` (exploratory, post hoc). The plan in `PLAN.md`, with Amendments 1 and 2, was
committed before any model was fit.

## Verdict against the pre-registered criteria

**Negative.** The pre-registered positive finding required H1 to be significant after Holm
correction *and* ΔAUROC ≥ 0.01. Neither condition was met:
- Network label propagation (M2-RWR) improves gene-stratified AUROC over AlphaMissense + gene prior
  (M1) by **+0.0022** (95% CI 0.0002–0.0044; p_boot = 0.026; **Holm p = 0.078**).
- That is about a fifth of the pre-set minimum effect.
- 3D neighbourhoods (RWR or Gaussian distance) are **not better than a 1D sequence window**. The 1D
  window is nominally better, and the difference lies inside the ±0.005 equivalence margin.
- In no bin does the propagation component reach the ClinGen **Supporting** likelihood ratio (2.406).

## Data

| | n |
|---|---|
| Test variants (genes with ≥5 TRAIN-labelled residues) | 55,975 (22,425 P / 33,550 B) |
| Genes | 3,153 |
| Genes contributing within-gene P–B pairs | 1,602 |

## Primary analysis (gene-stratified AUROC; 2,000-resample gene bootstrap)

| Model | Stratified AUROC | Pooled AUROC | AUPRC | Brier |
|---|---|---|---|---|
| M0 AlphaMissense | 0.9629 | 0.9518 | 0.925 | 0.084 |
| M1 + gene prior | 0.9627 | 0.9712 | 0.959 | 0.061 |
| M2-1D (λ = 5) | **0.9660** | 0.9749 | 0.965 | 0.056 |
| M2-3D (σ = 6 Å) | 0.9648 | 0.9756 | 0.966 | 0.054 |
| M2-RWR (α = 0.3) | 0.9649 | 0.9759 | 0.966 | 0.054 |

Kernel parameters were tuned on TRAIN by gene-grouped CV. They equal the pre-chosen defaults, so
the "fixed parameters" sensitivity analysis is identical to the primary analysis.

| Hypothesis | Contrast | Δ | 95% CI | p_boot | Holm p |
|---|---|---|---|---|---|
| H1 | RWR − M1 | +0.0022 | 0.0002, 0.0044 | 0.026 | 0.078 |
| H2 | RWR − 3D | +0.0001 | −0.0005, 0.0007 | 0.62 | 0.62 |
| H3 | RWR − 1D | −0.0011 | −0.0023, 0.0002 | 0.092 | 0.18 |
| (secondary) | 3D − M1 | +0.0021 | 0.0001, 0.0042 | 0.043 | |
| (secondary) | 1D − M1 | +0.0033 | 0.0017, 0.0052 | 0.0005 | |
| (secondary) | M1 − M0 | −0.0002 | −0.0003, −0.0001 | 0.0005 | |

For H3, RWR was the 3D model chosen on TRAIN CV. Its 90% CI is −0.0021 to −0.00003, inside the
±0.005 margin, so **3D is equivalent to 1D** and, if anything, slightly worse.

**Pooled AUROC (secondary) shows why it was demoted.** The gene prior raises pooled AUROC by +0.019
but changes within-gene ranking by −0.0002. Most of the apparent pooled gain is learning which genes
are dangerous. On pooled AUROC, H1 is +0.0046 (95% CI 0.0038–0.0058); that is still below 0.01.

**Mean within-gene AUROC** (721 genes with ≥3 P and ≥3 B): 0.954–0.956 for every model. No model
improves on AlphaMissense at the per-gene level.

## Clinical calibration (PM1-like evidence conditional on AlphaMissense)

The table gives the likelihood ratio for P vs B by bins of the propagation component, i.e. the M2-RWR
log-odds minus the M1 log-odds. TRAIN values are out of fold.

| Component bin | TRAIN LR | TEST LR | TEST n (P/B) |
|---|---|---|---|
| ≤ −2 | 1.25 | 2.03 | 49 / 36 |
| (−2, −1] | 0.75 | 0.59 | 1,045 / 2,656 |
| (−1, −0.5] | 0.60 | 0.55 | 2,171 / 5,922 |
| (−0.5, 0.5] | 0.95 | 1.02 | 14,251 / 20,827 |
| (0.5, 1] | 1.77 | 1.94 | 4,105 / 3,172 |
| (1, 2] | 1.56 | 1.30 | 791 / 912 |
| > 2 | 0.81 | 0.79 | 13 / 25 |

The highest LR is 1.94 on TEST and 1.77 on TRAIN, below the Supporting threshold (2.406). The lowest
is 0.55, which is not below 1/2.406 = 0.416. **No bin qualifies for any evidence strength**, in either
direction. TRAIN and TEST agree, so the component is calibrated; it is just weak.

## Gain by data density (pooled AUROC within bins of TRAIN-labelled residues)

| Labelled residues | n | M1 | M2-RWR | Δ |
|---|---|---|---|---|
| 5–9 | 9,770 | 0.960 | 0.965 | +0.005 |
| 10–19 | 10,758 | 0.963 | 0.968 | +0.005 |
| 20–49 | 15,204 | 0.965 | 0.972 | +0.006 |
| ≥ 50 | 20,243 | 0.977 | 0.981 | +0.004 |

The gain does not grow with label density.

## Sensitivity analyses and stratifications (H1, RWR − M1, gene-stratified)

| Analysis | n test | Δ | 95% CI | Holm p (H1–H3) |
|---|---|---|---|---|
| Primary | 55,975 | +0.0022 | 0.0002, 0.0044 | 0.078 |
| ≥2 stars in 2026 | 10,724 | +0.0045 | 0.0017, 0.0082 | 0.045 |
| VUS in 2023 only | 6,404 | +0.0046 | −0.0005, 0.0096 | 0.32 |
| Genes ≥20 labelled residues | 35,447 | +0.0025 | 0.0002, 0.0047 | 0.060 |
| Strict temporal (evaluated after cutoff) | 47,239 | +0.0023 | −0.0001, 0.0049 | 0.18 |
| Evaluated ≥ 2024-09 | 31,247 | +0.0034 | −0.0000, 0.0072 | 0.17 |
| Evaluated < 2024-09 | 24,055 | +0.0021 | −0.0011, 0.0053 | 0.51 |
| Genes without a VCEP PM1 rule | 54,971 | +0.0019 | −0.0001, 0.0040 | 0.14 |
| **Genes with a VCEP PM1 rule (21)** | 1,004 | **+0.0151** | 0.0037, 0.0542 | 0.0075 |
| Fixed parameters | 55,975 | +0.0022 | 0.0002, 0.0044 | 0.17 |
| Include same-residue labels | 55,975 | +0.0037 | 0.0011, 0.0062 | 0.0075 |

Notes:
- Every analysis except VCEP-PM1 genes gives Δ < 0.005, and the 1D window is never worse than 3D
  except in VCEP-PM1 genes (RWR − 1D +0.0043, 95% CI 0.0007–0.0155).
- **Circularity in VCEP-PM1 genes.** These are the genes whose expert panels award PM1 for
  hotspots, so their 2023–2026 classifications were partly *made* from neighbouring pathogenic
  variants. The +0.015 is the one place the pre-set effect size is reached, and it is also where the
  outcome is most circular. It was a pre-specified stratification, reported and not tested.
- **Same-residue analysis was corrected before interpretation.** The first implementation gave each
  TRAIN variant its own label inside its features. The model over-learned that leaked signal and
  stratified AUROC fell to 0.864. The corrected version uses a residue's label for a TRAIN variant
  only when *other* TRAIN variants supply it; TEST variants use all earlier labels at their residue.
  Only the corrected numbers are reported above. Primary features were never affected.
- **AlphaMissense circularity** (PP3 applied after 2023) would favour AlphaMissense. AlphaMissense is
  still at 0.956 stratified AUROC among variants evaluated before 2024-09, and the propagation gain is
  no larger there.

## Exploratory (post hoc, written after the results above were seen)

**AlphaMissense "ambiguous" band** (0.34 ≤ AM < 0.564; Cheng et al. 2023): 3,620 test variants in
1,165 genes. Here AlphaMissense alone ranks within genes at 0.598.

| Contrast | Δ stratified AUROC | 95% CI |
|---|---|---|
| RWR − M1 | +0.067 | −0.016, 0.137 |
| 3D − M1 | +0.062 | −0.024, 0.135 |
| 1D − M1 | +0.087 | 0.006, 0.161 |
| RWR − 1D | −0.020 | −0.057, 0.020 |

Calibration inside the band:
- **Pathogenic direction:** the (0.5, 1] component bin has TRAIN LR 2.78 and TEST LR 2.75. Its
  gene-bootstrap 5th percentile is 2.37, just below 2.406, so the bin does **not** qualify as
  Supporting under the pre-set rule.
- **Benign direction:** the (−1, −0.5] bin has TRAIN LR 0.38 and TEST LR 0.26. Its 95th percentile
  is 0.44, just above 0.416, so it also does not qualify.

This is the one setting where neighbour evidence could matter clinically, and even there a sequence
window does as well as structure. It is hypothesis-generating only.

## What this means
1. For variants that ClinVar reclassified between 2023 and 2026, AlphaMissense already ranks P above
   B within genes about 96% of the time.
2. Known neighbouring labels add a detectable but clinically negligible ~0.002–0.004 within genes.
3. **Structure is not what carries that small gain:** sequence-window neighbours do as well as or
   better than 3D distance or contact-graph network propagation.
4. Computational "hotspot" evidence built from neighbouring labels does not reach ClinGen Supporting
   strength once AlphaMissense is used. Counting both PP3 and an automated PM1 would risk
   double-counting.
