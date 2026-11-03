# Allelic-disorder analysis plan (written before any association was computed)

**Question.** Many rare-disease genes cause several distinct diseases, depending on which
missense variant a patient carries. Does a variant's mechanistic structural context, read
from the AlphaFold model, explain *which* disease it causes? In particular:
- is 3D proximity on the AlphaFold model a better guide than sequence proximity?
- does PhenoBridge's network-propagation (RWR) proximity add to plain 3D distance?

This plan was committed before the dataset was assembled and before any feature, distance
or label comparison was made. Anything not listed here is reported as exploratory.

## Lessons carried over from the phenopacket study
- **The unit is the structural site, not the patient.** Variants are collapsed to *residues*,
  because different substitutions at one residue share every structural feature.
- **Genes are the independent units for global inference.** Within-gene statistics are
  combined across genes with rank tests that treat each gene as one observation.
- **Every structural question has a sequence-only baseline.**

## Data (built by `allelic/a0*_*.py`; see `allelic/README.md`)
- ClinVar P/LP germline missense SNVs, ≥1 review star, linked to exactly one specific
  disease unit (OMIM phenotype MIM; MONDO-harmonised).
- phenopacket-store v0.1.27 missense variants with their carriers' diagnoses.
- Proteins with a single-fragment AlphaFold model (v6).

## Units and eligibility
- **Residue label.** A residue takes a disease label if every single-label variant at that
  residue carries the same disease unit. Residues with conflicting labels are excluded from
  the primary analysis and counted; a sensitivity analysis gives them the majority label,
  dropping ties.
- **Within a gene**, only disease units with ≥5 labelled residues are kept.
- **A gene is eligible** if it then has ≥2 disease units.

## Distances between two variant residues i and j in a gene
- **D_seq:** |i − j| (sequence distance).
- **D_3D:** Euclidean Cα–Cα distance in the AlphaFold model.
- **D_RWR:** 1 − s_ij, where s_ij is the symmetrised RWR affinity
  s_ij = ½ (p_i(j)/max_k p_i(k) + p_j(i)/max_k p_j(k)).
  Here p_i is the exact RWR stationary vector seeded at i on the 8 Å Cα contact graph
  (α = 0.3, as in PhenoBridge), so s_ij is near 1 when j is a strong network neighbour of i.

## Primary statistic: leave-one-residue-out nearest-neighbour prediction
- For each labelled residue i in a gene, predict its disease as the majority label of its
  k = 3 nearest *other* labelled residues under a distance D. Ties go to the label of the
  single nearest residue.
- The per-gene score is **balanced accuracy** (the mean recall over the gene's disease units),
  computed for each of D_seq, D_3D and D_RWR.
- Per-gene chance level is estimated by permuting residue labels 1,000 times.

## Primary hypotheses (two-sided; genes are the units)
- **H1 (3D beats sequence).** Across eligible genes, BA(D_3D) − BA(D_seq) > 0.
  Wilcoxon signed-rank test over genes.
- **H2 (network propagation beats 3D distance).** BA(D_RWR) − BA(D_3D) > 0.
  Wilcoxon signed-rank test over genes.
- Holm correction across H1 and H2.
- Reported with each test: median per-gene difference, the number of genes improving or
  worsening, and a bootstrap 95% CI over genes.

## Secondary analyses
- **S-a. Signal at all.** For each distance, test whether per-gene BA exceeds its
  permutation chance level. Per gene: permutation p. Global: Fisher's combined p across
  genes, plus the count of genes with p < 0.05.
- **S-b. Structural context features.** Within gene, test whether residue-level features
  differ between disease units (Kruskal–Wallis). Features:
  - pLDDT;
  - Core / Semi / Disorder RWR mass (α = 0.3);
  - contact degree;
  - burial proxy;
  - 1-hop neighbourhood pLDDT;
  - AlphaMissense (residue-mean over observed substitutions);
  - relative position.

  Global statistic: the number of genes with p < 0.05, against a label-permutation null
  (labels permuted among residues within each gene, 1,000 permutations).
  "Beyond position": features are residualised on a natural cubic spline of position
  (4 basis functions including the constant), and tested against a **rotation null**:
  disease labels are cyclically shifted along the ordered labelled residues. This keeps the
  sequence clustering of labels but breaks their alignment with the feature.
- **S-c. AD vs AR (Q3).** In genes with ≥5 labelled residues for both an autosomal-dominant
  and an autosomal-recessive disease unit:
  - per-gene difference (AR median − AD median) in the burial proxy, contact degree, pLDDT
    and Core mass;
  - pre-specified direction: AR (loss-of-function) residues are *more* buried and confident
    than AD residues;
  - Wilcoxon signed-rank test over genes, with a sign test as a check.

## Controls
- **Positive controls:** the list in `allelic/notes/literature.md`, fixed before analysis;
  results are reported per gene.
- **Calibration:** with labels permuted within genes, BA differences between distances should
  centre on 0. We report this simulated null distribution of the H1/H2 statistics, using
  200 complete label permutations across all genes.

## Sensitivity
- k ∈ {1, 5};
- α ∈ {0.15, 0.5} for D_RWR;
- majority labels for conflicting residues;
- ClinVar ≥2 stars only;
- ClinVar only (no phenopacket-store);
- excluding genes with mean pLDDT < 70;
- the D_3D threshold version, i.e. shared-neighbourhood vs raw distance (exploratory only).

## What counts as a positive finding
H1 or H2 significant after Holm correction, with the median per-gene improvement in BA
≥ 0.03, and consistent direction in ≥ 60% of genes. Otherwise it is reported as a null,
with the CI on the median difference.
