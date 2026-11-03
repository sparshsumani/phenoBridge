# KBG / ANKRD11 follow-up: results

Code: `kbg_followup/` (run `run_all.sh`; REST and AlphaMissense results are cached in `cache/`).
Plan: `ANALYSIS_PLAN.md`, written before any genotype–phenotype test was run.
Every number below is generated from `out/`, and the tables are produced by `07_tables.py` (`out/tables_md.txt`).

## 1. Data linkage and variant annotation

- All 333 JSON individuals match a phenopacket by subject ID.
- All 333 also agree with their phenopacket on the allele key, the sex, and all 48 HPO calls (0 mismatches).
- VEP (Ensembl REST, GRCh38) was run on the 162 distinct hg38 VCF records.
  - Consequences are read from ENST00000301030.10, which VEP reports as MANE Select (NM_013275.6).
  - VEP's HGVS p. strings and protein positions agree with the phenopacket HGVS p. for all 256 point variants.

**Variant classes (individuals):**

| Class | Individuals |
|---|---|
| frameshift | 174 |
| nonsense | 66 |
| canonical splice | 7 |
| missense | 7 |
| in-frame indel | 2 |
| 16q24.3 CNV | 77 |

The 77 CNVs are **70 deletions (x1), 6 duplications (x3), and 1 without a copy number**. That last one is Youngs2011, "16q24.3(87.92-88.10)". It does not overlap ANKRD11 in either hg19 or hg38.

So 247 of 256 point variants (96.5%) are truncating (frameshift, nonsense or splice).

**Regex "truncating-like" flag vs VEP.** 73 of 256 individuals (28.5%) are misclassified by the regex:
- 65 of 66 nonsense SNVs are called non-truncating (the 1 caught is a delins);
- 6 of 7 splice variants are called non-truncating;
- both in-frame deletions are called truncating.

The regex gives 69.5% truncating; the true figure is 96.5%.

| VEP class (MANE) | regex: not truncating | regex: truncating |
|---|---|---|
| frameshift | 0 | 174 |
| inframe | 0 | 2 |
| missense | 7 | 0 |
| nonsense | 65 | 1 |
| splice | 6 | 1 |

The regex residue mapping (⌈c/3⌉), by contrast, is fine: it agrees with VEP's protein position exactly in 252/256 cases and within 1 residue in all 256. The allele keys never contain p. notation, so every draft position came from ⌈c/3⌉.

## 2. Exon structure and NMD

- ENST00000301030 has **13 exons**. The CDS starts in exon 3.
- **Exon 9 is the giant exon**: 6,578 nt, c.893–7470, aa 298–2490, which is 82.3% of the CDS.
- **The last exon (13) carries only 186 coding nt** (c.7807–7992, aa 2603–2663).
- The last exon–exon junction lies after c.7806.

So ANKRD11 does *not* have a large last exon; the premise that it does is wrong.

**NMD prediction for the 240 PTC variants (frameshift + nonsense).** The PTC position for frameshifts is the new stop codon (p.…fsTerN).

| PTC location | Individuals |
|---|---|
| exon 9 | 224 |
| exon 11 | 7 |
| exon 8 | 3 |
| exon 6 | 2 |
| exon 13 | 2 |
| exon 7 | 1 |
| exon 10 | 1 |

- **50-nt rule** (PTC in the last exon or ≤50 nt upstream of the last junction): **2/240 (0.8%) escape NMD**. Both are p.(Glu2612Ter), c.7834G>T (KBG49, KBG7), in exon 13. The 55-nt rule gives the same 2/240.
- So 238/240 (99.2%) are predicted to trigger NMD by the canonical rule. The draft's statement is consistent with that rule.
- Caveat: 224/240 PTCs lie in the 6.6-kb exon 9. NMD efficiency is reduced for PTCs in long exons (>~400 nt; Lindeboom et al. 2016, Nat Genet). Walz et al. 2015 (Hum Genet; abstract verified) reported that mutant ANKRD11 proteins accumulate aberrantly.
- "Degraded by NMD" is therefore a rule-based prediction, not an established fact, and should be worded that way.
- No PTC is within the first 150 coding nt (the start-proximal escape rule).

**C-terminus.**
- 221/240 PTCs (92%) lie at or before aa 2369, so the protein loses the whole UniProt region 2369–2663 ("important for protein degradation", PMID 25413698).
- 227/240 (95%) lose all three C-terminal confident segments (2436–2455, 2489–2547, 2577–2631).
- None retains all three.
- UniProt Q6UB99 annotates 4 ANK repeats (167–292), many MobiDB-lite disordered regions, and the degradation region. It has **no "RD2/repression domain" feature**. RD1/RD2 come from the older literature (Zhang et al. 2004 JBC); I could not verify their coordinates here.

## 3. Duplicates and relatedness (`03_duplicates.py`)

- The individuals come from 41 sources. The 67 "KBG##" IDs are the new Martinez-Cayuelas 2023 (MC2023) patients.
- After the hotspot rule (alleles recurrent in ≥3 sources cluster only via explicit family labels), **58 individuals fall in 22 multi-member family clusters**. Examples:
  - Murray 2017 families 1 and 3 (3 members each);
  - the Sirmaci 2011 F1 trio (re-reports of Tekin 2004);
  - Ockeloen 2015 P1–P5 (c.7481dup);
  - Crippa 2015 P1–P3 (a duplication);
  - Goldenberg 2016 P1–P4 (the same deletion);
  - Kim 2015 P1–P3;
  - Low 2016 P31–P33;
  - Khalifa 2013 P1A/B;
  - KBG31A/B.
- **Labelling inconsistencies:**
  - KBG10A and KBG10B carry *different* alleles (c.6968_6975del vs c.2329_2332del), as do KBG8A and KBG8B (c.6628G>T vs c.7753C>T). An autosomal-dominant family should share its variant, so either the A/B suffix does not mean "same family" or an allele is mis-assigned.
  - KBG10B shares c.2329_2332del with KBG64.
  - "Low, 2016_P13 (27)" and "Low, 2016_P29 (27)" carry the same original-patient number but different alleles.
- **Cross-source duplicates cannot be identified from these data.**
  - 20 pairs meet the "likely duplicate" rule: same allele, different source, compatible sex, ≥5 shared terms, 0 discordant.
  - But pairs carrying *different* alleles (who cannot be the same person) meet the same phenotype criterion at a rate of 8.0%. That predicts **≈14 false hits** among the eligible same-allele pairs.
  - Most hits involve the c.1903_1907del hotspot.
  - The strongest candidates (≥12 shared terms, 0 discordant) are Gnazzo 2020 P10 vs Parenti 2021 P4 (c.2398_2401del) and Gnazzo 2020 P4 vs Parenti 2021 P18 (c.1903_1907del). Both are Italian series, but both pairs have discordant sex.
  - The rare c.1977C>G is carried by Gnazzo P5, Gnazzo P6 and Parenti P3.
  - These candidates should be checked against the source papers.
- **Sensitivity:** one per family cluster leaves n = 297; one per family-or-likely-duplicate cluster leaves n = 281. Both were run for the primary analysis (S4 below).

## 4. Features

`out/individuals.csv` holds one row per individual (333), with these features:
- variant class and the truncating flag;
- VEP position and relative position;
- PTC position, exon, and NMD-escape (50 and 55 nt) and long-exon flags;
- retained length;
- loss of the 2369–2663 region and of the C-terminal confident segments;
- the four RWR components (repo code unchanged, α = 0.3, r = 8 Å), seed pLDDT, and log10 SemiCore;
- AlphaMissense (missense only) and the per-position mean AlphaMissense;
- SV type, copy number, size, and overlap with ANKRD11 in hg19 and hg38;
- whether the variant is in exon 9, and whether it is the c.1903_1907del hotspot;
- the 48 HPO calls (1/0/NA).

The per-residue fingerprints are in `out/fp_all_residues.csv`.

## 5a. Primary: does log10 SemiCoreImpact add to variant class + relative position?

The analysis covers 256 point-variant individuals and **34 eligible terms** (≥10 observed and ≥10 excluded, explicit annotations only), using nested logistic LRTs with 1 df.

- **No term is significant.**
  - Minimum p = 0.080 (ADHD, OR 0.59 per SD, 95% CI 0.31–1.14).
  - Minimum BH q = 0.789.
  - 0/34 terms reach p < 0.05.
- **Omnibus permutation test** (sum of LR χ² = 30.8; 5,000 within-stratum permutations; permutation mean 36.5; 95th percentile 54.0): **p = 0.70**.
- There is no evidence of any structural signal from the pre-specified fingerprint.
- 6 terms showed ML separation of the class coefficient (all 9 non-truncating individuals share one value). For those, Firth was used, as the plan specifies.

Top 10 terms by p:

| label | n | obs | exc | OR (95% CI) per SD | p | q | Firth p |
|---|---|---|---|---|---|---|---|
| ADHD | 68 | 53 | 15 | 0.59 (0.31–1.14) | 0.080 | 0.789 | 0.123 |
| Microretrognathia | 52 | 16 | 36 | 1.84 (0.85–4.00) | 0.091 | 0.789 | 0.134 |
| Delayed skeletal maturation | 113 | 54 | 59 | 0.70 (0.45–1.09) | 0.110 | 0.789 | 0.118 |
| Intrauterine growth retardation | 62 | 13 | 49 | 1.84 (0.78–4.33) | 0.124 | 0.789 | 0.186 |
| Short stature | 191 | 115 | 76 | 0.79 (0.57–1.09) | 0.138 | 0.789 | 0.153 |
| Bulbous nose | 110 | 74 | 36 | 0.69 (0.41–1.17) | 0.152 | 0.789 | 0.180 |
| Hypotonia | 53 | 21 | 32 | 1.49 (0.83–2.67) | 0.165 | 0.789 | 0.208 |
| Feeding difficulties | 134 | 59 | 75 | 0.80 (0.55–1.16) | 0.230 | 0.789 | 0.244 |
| Abnormality of the hand | 214 | 155 | 59 | 0.83 (0.60–1.14) | 0.236 | 0.789 | 0.236 (Firth used) |
| Intellectual disability | 183 | 158 | 25 | 1.26 (0.86–1.85) | 0.242 | 0.789 | 0.214 |

One SD of log10 SemiCore is 0.519, so one 10-fold change corresponds to OR^(1/0.519).

**Pre-specified sensitivity analyses:**

| analysis | min p | min q | n p<.05 | n q<.05 | top term (OR) |
|---|---|---|---|---|---|
| primary | 0.080 | 0.789 | 0 | 0 | ADHD (0.59) |
| S1 Firth penalised LRT | 0.118 | 0.831 | 0 | 0 | Delayed skeletal maturation (0.72) |
| S2 position as B-spline (4 df) | 0.051 | 0.734 | 0 | 0 | Short stature (0.72) |
| S3 log10 CoreImpact | 0.006 | 0.200 | 3 | 0 | Bulbous nose (1.95) |
| S3 CommBackboneImpact* | 0.008 | 0.279 | 3 | 0 | Delayed skeletal maturation |
| **S3 DisorderImpact** | **2.5e-4** | **0.009** | 1 | **1** | **Intellectual disability (0.36)** |
| S3 seed pLDDT | 0.046 | 0.575 | 1 | 0 | Intellectual disability (2.00) |
| S4 one per family (n = 297) | 0.004 | 0.124 | 1 | 0 | Short stature (0.58) |
| S4 one per family/likely dup (n = 281) | 0.003 | 0.098 | 1 | 0 | Short stature (0.57) |
| S5 truncating only | 0.101 | 0.765 | 0 | 0 | Bulbous nose (0.64) |

\* The backbone impact is mostly ≈0 and very skewed, so its per-SD ORs are ill-conditioned. Do not interpret them.

**The one secondary hit: DisorderImpact vs intellectual disability.**
- OR 0.36 per SD; LRT p = 2.5×10⁻⁴.
- BH q = 0.009 within the DisorderImpact family; q = 0.035 across all 136 S3 tests; q = Bonferroni p = 0.043 across all 170 fingerprint tests (primary + S3).
- ID is recorded as observed in 96% of individuals with DisorderImpact ≤ 0.2, against about 79–82% when it is > 0.5.
- It rests on only 25 ID-excluded individuals. Most of them carry PTCs in the long disordered stretch, aa 570–1520.
- **Exploratory robustness checks (not in the plan; `04b_exploratory_disorder_ID.py`).** The association holds up in all of them:
  - Firth p = 4.1×10⁻⁴;
  - one per family, p = 7×10⁻⁴;
  - adjusted for source, p = 8×10⁻⁴;
  - excluding c.1903_1907del carriers, p = 4×10⁻⁴;
  - position as a cubic B-spline with 4, 6 or 8 df, p ≤ 0.003;
  - leave-one-source-out, p ≤ 0.010 (dropping MC2023 gives p = 0.010);
  - family-clustered SEs, p = 4×10⁻⁴;
  - an unadjusted Mann–Whitney test, p = 8.7×10⁻⁴.
- Interpretation: this is a hypothesis to replicate, not a finding.
  - It is a secondary component, and SemiCore, the pre-specified exposure, shows nothing (OR for ID 1.26, p = 0.24).
  - "ID excluded" in literature cohorts may reflect age at report (age is not available) or ascertainment.
  - No independent cohort was available here.

**Reproducing the draft (comparison only).** The draft's model is y ~ log10 Semi + Trunc, with terms ≥15/≥15.

| Trunc flag | Encoding | Terms | Nominal hits | min q |
|---|---|---|---|---|
| regex | A (NA=0) | 31 | 1 (ADHD, p = 0.023) | 0.57 |
| regex | B (explicit) | 30 | 0 | 0.77 |
| VEP | A (NA=0) | 31 | 1 (ADHD, p = 0.015) | 0.48 |
| VEP | B (explicit) | 30 | 0 | 0.77 |

Triangular face (p = 0.048 in the draft) loses nominal significance when the VEP positions are used (4 individuals shift by 1 residue). The draft's "two nominal hits" is therefore fragile.

## 5b. Baseline associations (Firth ORs, profile 95% CIs; BH within family)

None of the five families of tests has a term with q < 0.05.

| Family (terms) | min p | min q | Strongest term, OR (95% CI); prevalence exposed vs unexposed |
|---|---|---|---|
| B1 CNV vs point (34) | 0.036 | 0.46 | IUGR 3.28 (1.08–10.05); 8/17 vs 13/62 |
| B1s CNV deletion vs point (34) | 0.022 | 0.33 | Abnormality of hand 0.48 (0.26–0.90); 29/52 vs 155/214 |
| B2 missense/in-frame vs truncating (34) | 0.055 | 0.98 | Low anterior hairline 0.18 (0.02–1.03); 1/5 vs 54/83 |
| B3i PTC position, per 10% of length (34) | 0.016 | 0.56 | Cryptorchidism 1.25 (1.04–1.52) |
| B3ii PTC removes 2369–2663 vs later PTC (34) | 0.056 | 0.71 | Cryptorchidism 0.24 (0.04–1.03); 25/65 vs 6/8 |

Other nominal terms:
- B1: microretrognathia (3.69, p = 0.046).
- B1s: triangular face 0.36 (0.15–0.88), p = 0.024, 11/25 vs 68/99; microretrognathia, IUGR and VSD (all p ≈ 0.05).

Only 19 PTCs fall after aa 2369, and only 9 individuals have a non-truncating point variant. B2 and B3ii are therefore very underpowered.

**Reproduction of the MC2023 abstract claims** (same individuals, so this checks the data, not replication):

| Claim | Reported | This cohort | OR (95% CI) | p (Firth) | Fisher p |
|---|---|---|---|---|---|
| Triangular face: sequence vs CNV | 71.1% vs 42.5% | 68.7% vs 51.7% (99/29) | 2.03 (0.88–4.70) | 0.095 | 0.12 |
| (CNV deletions only, B1s) | | 68.7% vs 44.0% (99/25) | 2.7 (1.14–6.7)ᵃ | 0.024 | |
| Short stature: exon 9 vs outside | 62.5% vs 27.8% | 61.7% vs 50.0% (167/24) | 1.60 (0.69–3.76) | 0.27 | 0.37 |
| ID/ADHD/ASD: c.1903_1907del vs others | 70.4% vs 89.4% | 85.7% vs 90.9% (21/176) | 0.54 (0.17–2.22) | 0.36 | 0.43 |
| ID alone: c.1903_1907del vs others | — | 73.7% vs 87.8% (19/164) | 0.37 (0.13–1.19) | 0.093 | 0.15 |
| Macrodontia ↑ with deletion size (per doubling) | associated | n = 45 | 1.22 (0.92–1.67) | 0.15 | |
| Comorbidities ↑ with deletion size (per doubling) | associated | n = 43 | 1.31 (1.02–1.79) | 0.037 | |
| Hand anomalies ↑ with smaller deletions (per doubling) | associated | n = 52 | 0.76 (0.51–1.03) | 0.080 | |

ᵃ This is the inverse of the B1s OR 0.36 (0.15–0.88).

- Directions agree with MC2023 for every claim.
- The magnitudes are smaller for short stature/exon 9 and for c.1903_1907del.
- **Exploratory:** counting CNVs as "outside exon 9" gives 61.7% vs 51.6%, Fisher p = 0.15. That still does not reproduce 27.8%.
- MC2023 used its own full phenotype records and a "phenotypical score". Here only explicit HPO calls are available.

## 5c. Missense subset

- **n = 7 missense** (+2 in-frame). That is below 15, so the plan allows counts only and no tests.
- 6 of 7 missense variants lie in the confident C-terminal domain, aa 2491–2585 (exons 10–12; pLDDT 83–92):
  - p.Arg2512Trp ×2, p.Arg2512Gln, p.Arg2536Gln, p.Arg2585Cys, p.Ile2491Leu (the last at an exon-10 splice-region base);
  - all six have AlphaMissense 0.79–0.98 (likely pathogenic), and all have SemiCore ≤ 0.18 with CoreImpact ≥ 0.82.
- The seventh, p.Thr1721Met (VanDongen 2019 P13, pLDDT 32), has AlphaMissense 0.21 (likely benign) and only 3 annotated terms.
- **So missense variants are not "spread across the gene".** They cluster in the C-terminal structured region.

## 5d. Power (`05_power.py`)

The simulation uses each term's real covariate rows, its fitted baseline model, and 1,000 simulations per grid point.

- **Minimum detectable OR (80% power) per SD of log10 SemiCore:**

  | Term set | α = 0.05 | α = 0.005 |
  |---|---|---|
  | best terms (hand, short stature, macrodontia; n ≈ 190–214) | 1.53–1.59 | 1.74–1.83 |
  | median across the 34 terms (median n = 88) | 1.96 | 2.42 |
  | terms with no 80% power even at OR 4 | 0 | 10 |

- In 10-fold units, the median MDE is 3.6 (α = 0.05) and 5.5 (α = 0.005) per 10-fold change in SemiCore.
- Mean power across terms:

  | True OR per SD | α = 0.05 | α = 0.005 |
  |---|---|---|
  | 1.5 | 41% | 16% |
  | 2.0 | 74% | 49% |

- The empirical type-I error of the ML LRT at α = 0.05 was 6.2% on average (slightly liberal; Firth agrees with it).
- **Conclusion:** the null excludes only large effects, roughly OR ≥ 2 per SD for a typical term. Modest structural effects (OR 1.2–1.5 per SD) would be missed. The study cannot support "no effect", only "no large effect". The individual CIs make the same point (e.g. ID 0.86–1.85, short stature 0.57–1.09).

## 5e. Sparse events
Firth (S1) agrees with ML throughout: minimum Firth p = 0.118, minimum q = 0.83.

## Figures (`fig/`, SVG + PDF, FreeSerif)

1. `fig1_lollipop` (7 in): needle plot of the 256 point variants by class (with marker shapes as secondary encoding), the pLDDT track with confident segments, exons, the ANK and degradation regions, and the 50-nt NMD-escape boundary (aa ≈ 2586). The NMD rule applies to the *PTC* position; for frameshifts that lies downstream of the plotted variant position.
2. `fig2_power` (3.4 in): median power across the 34 terms (band: 10th–90th percentile) vs the true OR per SD, at α = 0.05 and 0.005.
3. `fig3_forest` (7 in):
   - (a) primary per-term ORs per SD of log10 SemiCore;
   - (b) CNV vs point-variant ORs (Firth).

## 7. Errors and overstatements in `paper/paper.html`

| # | Draft says | Correct statement |
|---|---|---|
| 1 | Abstract: "Missense and truncating variants are spread across the gene". | 96.5% of point variants are truncating, and there are only 7 missense individuals. 6 of the 7 missense variants cluster in aa 2491–2585 (exons 10–12). Truncating variants span aa 169–2612, but 93% of PTCs are in exon 9. |
| 2 | §2/Table 1: "77 individuals carry structural variants (16q24.3 deletions)". | These are 70 deletions, 6 **duplications** (x3: Crippa 2015 P1–P3, Bucerzan 2020, Kutkowska-Kazmierczak 2021 P17, Palumbo 2016), and 1 CNV with no copy number (Youngs 2011, "16q24.3(87.92-88.10)"). That last CNV does not overlap ANKRD11 in hg19 or hg38. SV coordinates come from mixed assemblies. |
| 3 | Table 1/§3.4: "Truncating-like among mapped 69.5%". The abstract and Table 2 call the model "adjusted for variant class". | The regex flag misclassifies 73/256 individuals (28.5%): 65/66 nonsense SNVs and 6/7 splice variants are called non-truncating, and both in-frame deletions are called truncating. The true truncating fraction is 96.5% (VEP, MANE). The models were adjusted for a flag that is mostly wrong. Correcting it does not change the null (min q 0.77, Encoding B). |
| 4 | §3.4: "Protein positions are read from p. notation when it is present". | The allele keys never contain p. notation, so every position came from ⌈c/3⌉. This is harmless: positions agree with VEP within 1 residue for all 256, and exactly for 252. The flag lists "fs" and "ter/*", but cDNA keys never contain them, so the flag effectively matches only del/dup/ins. |
| 5 | §3.3: "The three pLDDT regions partition V, so their components sum to one". | Not exactly. The code's Disorder (pLDDT ≤ q0.5) and Semi-core (pLDDT ≥ q0.5) overlap at the 3 residues whose pLDDT equals the median, 31.19 (aa 666, 967, 2105). Component sums reach 1.37 at residue 666 and exceed 1.01 for 33 seeds. No cohort seed is affected materially (max 1.005). |
| 6 | Discussion: "Fingerprints are likely to be more informative for missense variants in the confident C-terminal segments. Testing that subset needs a larger cohort than the 31 individuals available here." | The 31 are individuals with a seed pLDDT ≥ 70, of whom only 6 carry a missense variant. The others are 11 frameshift, 8 nonsense, 5 splice and 1 in-frame. The whole cohort has 7 missense individuals. |
| 7 | Discussion: truncating variants are "probably degraded through NMD". | By the canonical 50–55-nt rule, 238/240 PTCs (99.2%) are predicted to trigger NMD; only the 2 exon-13 PTCs (p.Glu2612Ter) escape. 224/240 lie in the 6,578-nt exon 9, however, where NMD efficiency is reduced (long-exon rule), and Walz et al. 2015 reported accumulation of mutant protein. This should read: "predicted, by the 50-nt rule, to be NMD substrates; not experimentally established". (ANKRD11 does **not** have a large last exon; exon 13 has 186 coding nt.) |
| 8 | Results/Abstract: Encoding A gives "two nominal signals" (ADHD p = 0.026, triangular face p = 0.048). | These are reproducible only with the draft's own ⌈c/3⌉ positions. With VEP positions and either flag, only ADHD is nominal (p = 0.023 with the regex flag, 0.015 with the VEP flag); triangular face is > 0.05. The conclusion is unchanged. |
| 9 | Intro: "Large cohort studies have described … only modest genotype–phenotype correlations [3]". | Acceptable, but MC2023 explicitly reports several correlations: triangular face higher in sequence variants than CNVs; short stature higher with exon 9 variants; fewer ID/ADHD/ASD with c.1903_1907del; deletion size vs macrodontia, comorbidities and hand anomalies. These should be cited specifically. In this cohort they reproduce in direction but not in magnitude, and none reaches q < 0.05. |
| 10 | Conclusion/Discussion treat the null as showing that structure is uninformative. | The cohort can detect only OR ≥ ~2 per SD for a typical term (median MDE 1.96 at α = 0.05, 2.42 at α = 0.005). Power at OR 1.5 per SD is 41% (α = 0.05). The null rules out large effects only. |
| 11 | §2: the cohort is "curated from the literature, primarily the 340-patient series of Martinez-Cayuelas". | All 333 rows carry PMID 36446582, but they derive from 41 original sources. Only 67 are new MC2023 patients. 58 individuals sit in 22 family clusters (relatives are treated as independent). Some families are mislabelled: KBG10A/B and KBG8A/B carry different alleles. Cross-source duplicates cannot be excluded. |

Notes (not errors):
- UniProt Q6UB99 annotates 4 ANK repeats (167–292) and a region 2369–2663 "important for protein degradation" (Walz 2015). It has no RD1/RD2 feature, so any future mention of RD2 needs a literature citation.
- 92% of PTCs remove the 2369–2663 region, so "C-terminal region lost" is almost collinear with "truncating".
- The draft reports ORs per 10-fold change in SemiCore, while the cohort SD is 0.52 log10 units. ORs per SD would be easier to read.

Items 3, 6, 7, 10 and 11 are substantive. Items 4, 5 and 8 are minor.
