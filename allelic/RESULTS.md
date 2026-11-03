# Allelic-disorder study: results (full run)

Plan: `PLAN.md` (with Amendment 1), committed before data assembly.

Run settings: `python3 run_allelic.py --perms 1000 --calib 200 --sens-perms 200`.
Data: ClinVar variant_summary (2026-11-01) and phenopacket-store v0.1.27; AlphaFold v6.

## Units
- 3,403 labelled residues; 168 residues with conflicting labels were excluded.
- 71 eligible genes, covering 2,625 residues.
- A gene is eligible with ≥2 specific disease units, each with ≥5 labelled residues.
  Broad MONDO grouping units are excluded.

## Primary hypotheses
Leave-one-residue-out 3-nearest-neighbour balanced accuracy (BA), compared within genes.
Genes are the units of the Wilcoxon signed-rank test.

| Hypothesis | Contrast | Median Δ | Mean Δ | Better / worse / tie | p | Holm p |
|---|---|---|---|---|---|---|
| H1 | 3D − sequence | 0.000 | +0.031 | 34 / 26 / 11 | 0.169 | 0.47 |
| H2 | RWR − 3D | 0.000 | −0.010 | 21 / 29 / 21 | 0.170 | 0.47 |
| H2b | RWR − hop | 0.000 | +0.012 | 25 / 23 / 23 | 0.155 | 0.47 |

- The bootstrap 95% CI of the median Δ is within [−0.013, 0.031] for all three.
- The pre-specified criterion for a positive finding (Holm p < 0.05, median Δ ≥ 0.03,
  ≥ 60% of genes improving) is **not met**.
- **Calibration** (200 complete label permutations): the median Δ was 0 for all three
  hypotheses, and the fraction of p < 0.05 was 0.055 (H1), 0.07 (H2) and 0.035 (H2b),
  i.e. the tests are well calibrated.
- **Sensitivity:**
  - settings varied: k = 1 or 5; RWR α = 0.15 or 0.5; majority labels; ClinVar ≥2 stars
    (23 genes); ClinVar labels only; mean pLDDT ≥ 70 (47 genes); dropping single-site
    diseases;
  - result: no H1/H2/H2b test reaches Holm p < 0.05;
  - smallest raw p: H2 at α = 0.5, p = 0.064, with RWR *worse* than 3D.

## Secondary results
- **S-a. Signal at all.** Nearest-neighbour labels predict disease above chance under
  every distance (Fisher combined p ≤ 7×10⁻⁶), but the effect is small:

  | Distance | Median BA | Chance | Genes with p < 0.05 |
  |---|---|---|---|
  | sequence | 0.495 | 0.461 | 11 of 71 |
  | 3D | 0.502 | 0.463 | 14 |
  | hop | 0.500 | 0.459 | 12 |
  | RWR | 0.500 | 0.462 | 14 |

- **S-b. Features by disease (Kruskal–Wallis; genes with p < 0.05 out of 71).**
  - Unadjusted for position, several features differ by disease more often than the
    null predicts (null mean ≈ 3.3 genes):

    | Feature | Genes (unadjusted) | Perm p |
    |---|---|---|
    | relative position | 12 | 0.005 |
    | AlphaMissense | 11 | 0.005 |
    | Core RWR mass | 11 | 0.005 |
    | contact degree | 9 | 0.005 |
    | burial | 9 | 0.020 |
    | pLDDT | 7 | 0.055 |
    | Semi RWR mass | 6 | 0.16 |

  - After spline adjustment for position, tested against the rotation null, only
    **AlphaMissense** remains (9 genes vs 4.3 expected; p = 0.020). Every structural
    feature has p ≥ 0.26.
- **S-c. AD vs AR.** Only 10 genes are evaluable. No feature differs consistently: for
  burial the median AR − AD is +0.75 Cα, 5 of 10 genes are positive, and p = 0.55.

## Controls
- **Positive:**
  - LMNA (4 diseases) and RET are predicted above chance, but **equally by sequence and
    3D** (LMNA BA 0.43 sequence vs 0.39 3D, chance 0.23; RET 0.85 under all distances).
  - GCK, GJB2, CASR and ABCA4 are at chance.
  - TP63 and CLCN7 are too small.
- **Null control:** FGFR2 (Crouzon vs Pfeiffer) was expected to be at chance but shows
  3D BA 0.75 vs chance 0.47 (p = 0.019).

## Conclusion
With labels and inference that are calibrated, which disease a missense variant causes is
only weakly predictable from where it sits, and AlphaFold 3D or contact-network proximity
(including PhenoBridge's RWR) adds nothing detectable to sequence proximity. The only
feature that adds information beyond position is AlphaMissense. This matches the
phenopacket-store HPO study (`manuscript/`).
