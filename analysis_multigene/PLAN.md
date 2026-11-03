# Multi-gene analysis plan (written before looking at the multi-gene data)

**Question.** Across Mendelian disease genes in phenopacket-store, do AlphaFold-derived
structural-context features of a patient's variant explain *within-gene* differences in
which HPO phenotypes are present, beyond simple baselines?

The plan was fixed before any association was computed. Anything not listed here is
reported as exploratory.

## Units and inclusion
- **Individual.** Has exactly one causative allele (heterozygous or hemizygous) that can
  be mapped to a protein position. Biallelic individuals are analysed only in a
  sensitivity analysis, where the feature is the mean over the two alleles.
- **Gene.** Has a single-fragment AlphaFold model and at least 20 included individuals.
- **Phenotype cells.** ontology-propagated observed/excluded calls (primary) or explicit
  calls only (sensitivity). Missing calls are never imputed.
- **Test unit (gene g, term h).** Among carriers in the analysis subset, at least 5
  observed **and** at least 5 excluded for h. Terms whose observed carriers are an exact
  superset or subset of another tested term's carriers in the same gene are kept, but
  the global inference uses permutation, which accounts for this dependence.

## Analyses
### A. Missense-only (primary; this is where 3D context should matter)
Subset: missense carriers.

Features, each first residualised on relative sequence position (pos/L) within the gene
using ordinary least-squares with a linear term, so that a structural feature has to add
something beyond "where along the chain":
1. **S1 (primary):** log10 SemiCoreImpact (RWR fingerprint, alpha = 0.3, cutoff 8 Å)
2. S2: seed pLDDT
3. S3: log10 CoreImpact
4. S4: contact degree
5. S5: AlphaMissense pathogenicity (reference feature)
6. P0: relative position itself, not residualised (baseline reference)

Per (g, h) test: two-sided Mann–Whitney U comparing the feature between individuals with
h observed and those with h excluded. Effect size = rank-biserial correlation.

### B. Variant class (reference analysis)
Subset: all mappable individuals. For each (g, h), Fisher exact test of missense vs
protein-truncating (stop_gained or frameshift). This is the genotype–phenotype signal
already known to be common, and it calibrates how much signal the cohorts can show.

### C. Global inference (answers "is there structural signal at all?")
For each feature, pool all (g, h) tests and report:
- the number of BH discoveries at q < 0.1, with BH applied across all tests for that feature;
- Storey's pi1 = 1 − pi0, using lambda = 0.5;
- a **within-gene permutation null**: shuffle feature values among carriers within each
  gene (this keeps the phenotype correlation structure), 1,000 permutations; the global
  statistic is the number of tests with p < 0.05; empirical p = (1 + #perm ≥ obs)/(1 + B).

The primary claim rests on the global permutation p for S1, compared with S2, S5 and P0.

### D. Positive controls
Genes with published structure- or location-driven correlations (list from the literature
review, fixed before running) are reported individually. We check whether the published
direction is recovered. This is a sanity check of sensitivity, not new inference.

### E. KBG / ANKRD11
Reported separately (kbg_followup/), with its own pre-specified plan.

## Sensitivity
- Explicit-only phenotype calls.
- One individual per family/pedigree where identifiable.
- RWR alpha in {0.15, 0.5}.
- Excluding genes with mean pLDDT < 50.

## What would count as a positive finding
S1 global permutation p < 0.01, with pi1 clearly above that of the permuted data, and
the signal not explained by P0 (structural features are residualised on position).
A null would be reported as a null, with the power implied by the sample sizes.
