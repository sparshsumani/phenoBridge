# Literature notes for the PhenoBridge paper (KBG case study → multi-gene study)

Compiled 2026-11-03. Bibliographic data were checked against PubMed E-utilities (esummary/efetch), Crossref, UniProt, Ensembl REST and publisher pages. Statements from abstracts or full texts that I read are marked **[V]** (verified). Statements I could not check against a primary source are marked **[UNVERIFIED]**. GPSEA per-cohort results come from the authors' public analysis repository (https://github.com/P2GX/gpsea-cs, files `supplement/tex/<GENE>_summary.tex` and `supplement/Table S2.xlsx`), which I cloned and read **[V]**.

---

## 1. Reference verification (the draft's 15 references)

**Summary.** All 15 references are bibliographically correct: authors, title, journal, year, volume and pages all match PubMed or Crossref, and I found **no wrong citation**. Three things should still change:
- Add DOIs, and PMIDs where they exist. At present only ref 3 has a PMID.
- Add an issue number if the house style uses one.
- Fix the content and in-text problems listed after the table.

The corrected lines below use the draft's style with DOI and PMID added.

| # | Status | Corrected citation |
|---|---|---|
| 1 | OK | Sirmaci A, Spiliopoulos M, Brancati F, et al. Mutations in ANKRD11 cause KBG syndrome, characterized by intellectual disability, skeletal malformations, and macrodontia. *Am J Hum Genet* 89(2):289–294, 2011. doi:10.1016/j.ajhg.2011.06.007. PMID 21782149. |
| 2 | OK | Low K, Ashraf T, Canham N, et al. Clinical and genetic aspects of KBG syndrome. *Am J Med Genet A* 170(11):2835–2846, 2016. doi:10.1002/ajmg.a.37842. PMID 27667800. |
| 3 | OK | Martinez-Cayuelas E, Blanco-Kelly F, Lopez-Grondona F, et al. Clinical description, molecular delineation and genotype–phenotype correlation in 340 patients with KBG syndrome: addition of 67 new patients. *J Med Genet* 60(7):644–654, 2023. doi:10.1136/jmg-2022-108632. PMID 36446582. |
| 4 | OK | Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596(7873):583–589, 2021. doi:10.1038/s41586-021-03819-2. PMID 34265844. |
| 5 | OK | Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Res* 50(D1):D439–D444, 2022. doi:10.1093/nar/gkab1061. PMID 34791371. |
| 6 | OK | Tunyasuvunakool K, Adler J, Wu Z, et al. Highly accurate protein structure prediction for the human proteome. *Nature* 596(7873):590–596, 2021. doi:10.1038/s41586-021-03828-1. PMID 34293799. |
| 7 | OK | Akdel M, Pires DEV, Pardo EP, et al. A structural biology community assessment of AlphaFold2 applications. *Nat Struct Mol Biol* 29(11):1056–1067, 2022. doi:10.1038/s41594-022-00849-w. PMID 36344848. |
| 8 | OK | Atilgan AR, Durell SR, Jernigan RL, et al. Anisotropy of fluctuation dynamics of proteins with an elastic network model. *Biophys J* 80(1):505–515, 2001. doi:10.1016/S0006-3495(01)76033-X. PMID 11159421. |
| 9 | OK | Tong H, Faloutsos C, Pan J-Y. Fast random walk with restart and its applications. In *Proc. 6th IEEE Int. Conf. on Data Mining (ICDM'06)*, 613–622, 2006. doi:10.1109/ICDM.2006.70. |
| 10 | OK | Köhler S, Bauer S, Horn D, Robinson PN. Walking the interactome for prioritization of candidate disease genes. *Am J Hum Genet* 82(4):949–958, 2008. doi:10.1016/j.ajhg.2008.02.013. PMID 18371930. |
| 11 | OK | Köhler S, Gargano M, Matentzoglu N, et al. The Human Phenotype Ontology in 2021. *Nucleic Acids Res* 49(D1):D1207–D1217, 2021. doi:10.1093/nar/gkaa1043. PMID 33264411. (Optional: a newer HPO update paper may exist. I did not check.) |
| 12 | OK | Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 57(1):289–300, 1995. doi:10.1111/j.2517-6161.1995.tb02031.x. |
| 13 | OK | Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422–1423, 2009. doi:10.1093/bioinformatics/btp163. PMID 19304878. |
| 14 | OK | Hagberg AA, Schult DA, Swart PJ. Exploring network structure, dynamics, and function using NetworkX. In *Proc. 7th Python in Science Conf. (SciPy 2008)*, 11–15, 2008. doi:10.25080/TCWV9851. |
| 15 | OK | Seabold S, Perktold J. Statsmodels: econometric and statistical modeling with Python. In *Proc. 9th Python in Science Conf.*, 92–96, 2010. doi:10.25080/Majora-92bf1922-011. |

### Content and in-text issues found while verifying (fix these)

1. **Martinez-Cayuelas [3] reports positive genotype–phenotype findings, so "only modest genotype–phenotype correlations" undersells it.** The abstract **[V]** reports:
   - a higher phenotypic score with sequence variants than with CNVs;
   - triangular face in 71.1% (sequence variants) vs 42.5% (CNVs);
   - short stature in 62.5% inside exon 9 vs 27.8% outside;
   - ID/ADHD/ASD less frequent with c.1903_1907del (70.4% vs 89.4%);
   - deletion-size associations.

   The study applied **no multiple-testing correction**. The GPSEA authors state this in their ANKRD11 notebook **[V]**, and only the phenotype-score difference replicated in GPSEA (see §2a). Suggested wording: "nominal correlations with variant class, exon 9 and one recurrent frameshift, of which only the CNV-vs-sequence-variant severity difference has been independently replicated".
2. **Cite GPSEA and phenopacket-store.** They are the data source and the closest prior work, and the draft cites neither (see §2a). This omission would be noticed immediately.
3. **Exon 9 is not the last exon of ANKRD11.** Ensembl canonical ENST00000301030 / NM_013275.6 has **13 exons**. Exon 9 is the largest at 6,578 bp **[V, Ensembl REST]**, but the CDS ends in exon 13. By the 50–55-nt rule, PTVs in exon 9 are therefore expected to trigger NMD. This matches Awamleh et al. 2023 **[V]**: "Variants in exon 9 are predicted to cause ANKRD11 haploinsufficiency via NMD", citing Walz 2015. The draft's Discussion sentence "haploinsufficiency through nonsense-mediated decay probably matters more" is consistent with this. Do **not** write that exon-9 PTVs escape NMD because they sit in the last exon.
4. **pLDDT < 50 and disorder [6,7].** The claim is supported. Consider adding Ruff & Pappu 2021, which cautions about reading pLDDT as a disorder measure (see §2b).
5. **"Annotated ankyrin-repeat domain".** UniProt Q6UB99 gives ANK1–4 at residues 167–196, 200–229, 233–262 and 266–292 **[V]**, consistent with the draft's confident segment 159–285. Cite UniProt.
6. **Problem in the draft's own data (not a literature issue).** The phenopackets were parsed with `python -I` from `/home/user/data_ext/ppkt/.../notebooks/ANKRD11/phenopackets` **[V]**. Of the 256 mapped individuals:

   | Variant class (from `hgvs.p`) | Individuals |
   |---|---|
   | Frameshift | 171 |
   | Nonsense SNV | 69 |
   | In-frame | 2 |
   | Missense | **7** |
   | `p.?` (likely splice) | 7 |

   The draft's string heuristic labels **65 of the 69 nonsense SNVs as non-truncating**, which produces the reported "69.5% truncating-like". The true PTV fraction is about 94%, and **only 7 individuals carry missense variants**. The T covariate in Eq. 3 is therefore badly misclassified. Use `hgvs.p` or VEP consequences. The 31 individuals with seeds at pLDDT ≥ 70 are probably mostly PTVs.

---

## 2. Novelty and closest prior work

### 2a. GPSEA and phenopacket-store (the closest prior work; must be cited and positioned against)

- **Danis D, Bamshad MJ, Bridges Y, et al. A corpus of GA4GH phenopackets: case-level phenotyping for genomic diagnostics and discovery. *HGG Adv* 6(1):100371, 2025. doi:10.1016/j.xhgg.2024.100371. PMID 39394689** (preprint medRxiv 2024, doi:10.1101/2024.05.29.24308104) **[V]**
  - The resource paper. v0.1.19 held 6,668 phenopackets, 475 diseases, 423 genes and 3,834 alleles from 959 publications.
  - It lists genotype–phenotype correlation as an intended use. The resource paper itself does no structural analysis.

- **Rekerle L, Danis D, Rehburg F, et al. GA4GH phenopacket-driven characterization of genotype-phenotype correlations in Mendelian disorders. *Am J Hum Genet* 113(1):57–70, 2026 (online 23 Dec 2025). doi:10.1016/j.ajhg.2025.12.001. PMID 41443197. PMCID PMC12824607** (preprint medRxiv 2025, doi:10.1101/2025.03.05.25323315) **[V, full text read]**

  *What they did:*
  - GPSEA v0.9.11 applied to **85 cohorts** (6,179 individuals; 81 genes; 122 diseases; median 49 per cohort, range 16–462) from phenopacket-store **v0.1.25**.
  - They found **253 significant GPCs, in 48/85 cohorts**. Previously published GPCs existed for 29 of those cohorts. They tested 9,013 hypotheses in total; 217 of the 6,736 categorical tests were significant.

  *Statistics:*
  - Fisher exact tests on observed/excluded HPO terms. Unannotated terms are omitted, which is the draft's "Encoding B"; GPSEA does this by design.
  - HPO annotations are propagated by the true-path rule.
  - **IF-HPO independent filtering** removes:
    - the top HPO levels;
    - terms whose counts equal a child's;
    - terms annotated in < 40% of the cohort;
    - tables too small to reach nominal significance;
    - terms with an empty genotype class.

    This cut tested terms from a mean of 304 to 40.
  - Benjamini–Hochberg correction over the remaining terms.
  - Also Mann–Whitney tests on phenotype scores (counting score, de Vries score), t-tests on measurements, and log-rank tests on onset or survival.

  *Genotype predicates:*
  - Variant effect, specific variant, exon, **UniProt protein region/domain**, SV type, sex, diagnosis, allele count.
  - The "protein features" are UniProt annotations used for visualisation and region predicates.
  - **No AlphaFold, pLDDT, 3D distance, contact network or pathogenicity-score features are used anywhere.** The full text has no mention of AlphaFold or 3D structure.

  *Hypothesis choice:* Hypotheses were chosen per cohort by hand, from the literature or from the lollipop plot. They were not systematic.

  *ANKRD11 cohort* (337 individuals, store v0.1.24, 48 HPO terms) **[V]**:
  - FET for SV vs other, exon 9 vs other, and c.1903_1907del vs other: **0 significant results** in each of 23 tests.
  - The only significant results:
    - a counting phenotype score lower for SV than for other variants (4.62 vs 5.94; Mann–Whitney p = 2.6 × 10⁻⁴), replicating Martinez-Cayuelas;
    - a sex difference (p = 0.0075).

- **Gap our paper fills.** GPSEA tested hand-picked, categorical, 1D genotype classes in 85 cohorts. Our study would test:
  - a *systematic, pre-registered, continuous* structural covariate, derived from AlphaFold, for every cohort of ≥ 30 individuals among ~700 (80 cohorts at ≥ 30 in v0.1.27);
  - with explicit baselines: variant class, linear position or UniProt domain, AlphaMissense;
  - and with calibration on positive and negative controls.

  The question is whether 3D context adds information **beyond** what GPSEA's region and variant-class predicates already capture.

- **"Already done" risk.** A reviewer will say GPSEA already showed widespread GPCs from domain and region predicates. We must show incremental value over a "UniProt domain + variant class" baseline, reuse GPSEA's observed-only encoding and IF-HPO-like filtering (or justify alternatives), and ideally reproduce GPSEA's significant GPCs (Table S2) as positive controls.

### 2b. Protein structure, AlphaFold and spatial-clustering features related to phenotype

**Phenotype differences within a gene or gene family (closest to our question)**

- **Brünger T, Ivaniuk A, Pérez-Palma E, et al.** Conserved missense variant pathogenicity and correlated phenotypes across paralogous genes. *Genome Biol* 26:197, 2025. doi:10.1186/s13059-025-03663-x. PMID 40624578 **[V]**. Across 10 voltage-gated sodium-channel genes and 22 expert-curated disorders, they find cross-paralog phenotype correlation based on **3D spatial position**: LoF-related disorders in SCN1A/2A/5A/8A form overlapping spatial clusters. This is the closest "3D position → phenotype" analysis, but it covers one gene family with curated disorder labels, not HPO terms in phenopackets.
- **Brünger T, Pérez-Palma E, Montanucci L, et al.** Conserved patterns across ion channels correlate with variant pathogenicity and clinical phenotypes. *Brain* 146(3):923–934, 2023. doi:10.1093/brain/awac305. PMID 36036558 **[V]**. They computed 163 structural features for 30 ion channels and a 3D spatial-distance score. Pore-lining residues near the pore axis carry the strongest pathogenic enrichment, and the identified properties are associated with clinical phenotypes and functional data.
- **Brunklaus A, Feng T, Brünger T, et al.** Gene variant effects across sodium channelopathies predict function and guide precision therapy. *Brain* 145(12):4275–4286, 2022. doi:10.1093/brain/awac006. PMID 35037686 **[V]**. Variants at paralogous positions across sodium channels have similar functional effects (35 of 38 pairs). Relevant to "position predicts GoF vs LoF predicts phenotype".
- **Heyne HO, Baez-Nieto D, Iqbal S, et al.** Predicting functional effects of missense variants in voltage-gated sodium and calcium channels. *Sci Transl Med* 12(556):eaay6848, 2020. doi:10.1126/scitranslmed.aay6848. PMID 32801145 **[V title only]**. A GoF/LoF predictor using sequence and structure features. Phenotype follows mechanism in these channels.
- **Badonyi M, Marsh JA.** Prevalence of loss-of-function, gain-of-function and dominant-negative mechanisms across genetic disease phenotypes. *Nat Commun* 16:8392, 2025. doi:10.1038/s41467-025-63234-3. PMID 40998763 **[V]**. A structure-based missense-LoF score applied to 2,837 phenotypes in 1,979 genes. 43% of dominant and 49% of mixed-inheritance genes contain both LoF and non-LoF mechanisms, and combining the score with phenotype semantic similarity prioritises dominant-negative mechanisms. **This is the closest systematic, multi-gene, structure → phenotype work.** It works at the gene/disease level, not on individual patients' HPO terms.
- **Gerasimavicius L, Livesey BJ, Marsh JA.** Loss-of-function, gain-of-function and dominant-negative mutations have profoundly different effects on protein structure. *Nat Commun* 13:3895, 2022. doi:10.1038/s41467-022-31686-6. PMID 35794153 **[V title]**. Non-LoF missense variants are milder structurally (ΔΔG) and cluster in 3D. This motivates using 3D clustering as a mechanism proxy.
- **Ng JK, Chen Y, Akinwe TM, et al.** Proteome-wide assessment of differential missense variant clustering in neurodevelopmental disorders and cancer. *Cell Genom* 5(4):100807, 2025. doi:10.1016/j.xgen.2025.100807. PMID 40073865 (preprint PMID 38352539) **[V]**. 3D-CLUMP on AlphaFold structures finds 18 proteins clustering differentially in NDD vs cancer and 19 the other way round: phenotype-specific 3D clustering at the cohort level, but only two coarse phenotype classes.
- **Gogate N, Jolly A, Rosenfeld JA, et al.** Domain specific phenotypic expansion associated with variants in MACF1. *medRxiv* 2025. doi:10.1101/2025.06.26.25330137. PMID 40666329 **[V]**. HPO-based clustering of 29 individuals separates GAR-domain from non-GAR-domain phenotypes. A single-gene example of domain → HPO cluster.
- **Sun YJ, Velez G, Li A, et al.** Advancing genetic pathogenicity prediction with three-dimensional proteoform-phenotype analysis. *medRxiv* 2026. doi:10.64898/2026.10.07.26364649 **[V, abstract page]**. PRPH2: 46 pathogenic missense variants on AI structure models show distinct spatial patterns for RP vs macular vs pattern dystrophy. A single gene, not peer reviewed.
- **Andriessen E, de Boer E, Lyon GJ, et al.** Functional data strengthen clinical validation of PhenoScore phenotype-guided AI for ANKRD11 missense variants. *Clin Genet* 109(6):1038–1048, 2026. doi:10.1111/cge.70159. PMID 41786677 **[V]**. PhenoScore (facial + HPO) separates KBG from controls (AUC 0.95). It was trained on 79 PTV carriers and tested on 6 functionally confirmed missense carriers; it is not significantly different from AlphaMissense and is complementary to it. Relevant to "phenotype similarity" scoring in KBG.
- **Dingemans AJM, Hinne M, Truijen KMG, et al.** PhenoScore quantifies phenotypic variation for rare genetic diseases by combining facial analysis with other clinical features using a machine-learning framework. *Nat Genet* 55(9):1598–1607, 2023. doi:10.1038/s41588-023-01469-w. PMID 37550531 **[V title]**.

**Spatial or positional clustering of pathogenic variants (genotype only, no patient phenotypes)**

- **Sivley RM, Dou X, Meiler J, et al.** Comprehensive analysis of constraint on the spatial distribution of missense variants in human protein structures. *Am J Hum Genet* 102(3):415–426, 2018. doi:10.1016/j.ajhg.2018.01.017. PMID 29455857 **[V]**. Disease variants cluster in 3D (215 proteins); only 12% of those were detectable in 1D. A key argument that 3D adds information beyond position.
- **Quinodoz M, Peter VG, Cisarova K, et al.** Analysis of missense variants in the human genome reveals widespread gene-specific clustering and improves prediction of pathogenicity. *Am J Hum Genet* 109(3):457–470, 2022. doi:10.1016/j.ajhg.2022.01.006. PMID 35120630 **[V]**. Pathogenic variants cluster within genes in 387 of 840 genes (MutScore).
- **Pérez-Palma E, May P, Iqbal S, et al.** Identification of pathogenic variant enriched regions across genes and gene families. *Genome Res* 30(1):62–71, 2020. doi:10.1101/gr.252601.119. PMID 31871067 **[V title]**.
- **Iqbal S, Pérez-Palma E, Jespersen JB, et al.** Comprehensive characterization of amino acid positions in protein structures reveals molecular effect of missense variants. *PNAS* 117(45):28201–28211, 2020. doi:10.1073/pnas.2002660117. PMID 33106425 **[V title]**.
- **Iqbal S, Brünger T, Pérez-Palma E, et al.** Delineation of functionally essential protein regions for 242 neurodevelopmental genes. *Brain* 146(2):519–533, 2023. doi:10.1093/brain/awac381. PMID 36256779 **[V]**. Essential3D sites on AlphaFold and experimental structures.
- **Lelieveld SH, Wiel L, Venselaar H, et al.** Spatial clustering of de novo missense mutations identifies candidate neurodevelopmental disorder-associated genes. *Am J Hum Genet* 101(3):478–484, 2017. doi:10.1016/j.ajhg.2017.08.004. PMID 28867141 **[V title]**.
- **Turner TN, Douville C, Kim D, et al.** Proteins linked to autosomal dominant and autosomal recessive disorders harbor characteristic rare missense mutation distribution patterns. *Hum Mol Genet* 24(21):5995–6002, 2015. doi:10.1093/hmg/ddv309. PMID 26246501 **[V title]**.
- **Sevim Bayrak C, Stein D, Jain A, et al.** Identification of discriminative gene-level and protein-level features associated with pathogenic gain-of-function and loss-of-function variants. *Am J Hum Genet* 108(12):2301–2318, 2021. doi:10.1016/j.ajhg.2021.10.007. PMID 34762822 **[V title]**.

**AlphaFold, AlphaMissense and disorder (baselines and caveats)**

- **Cheng J, Novati G, Pan J, et al.** Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 381(6664):eadg7492, 2023. doi:10.1126/science.adg7492. PMID 37733863 **[V]**. The baseline score.
- **Ruff KM, Pappu RV.** AlphaFold and implications for intrinsically disordered proteins. *J Mol Biol* 433(20):167208, 2021. doi:10.1016/j.jmb.2021.167208. PMID 34418423 **[V title]**. A caveat on reading pLDDT for IDRs.
- **Hubrich D, Alvarado Valverde J, Lee CY, et al.** Variant characterization in the intrinsically disordered human proteome. *Nat Struct Mol Biol* 33(8):1183–1193, 2026. doi:10.1038/s41594-026-01846-z. PMID 42538434 **[V]**. 37% of annotated missense variants fall in the disordered 25% of the proteome. They modelled 1,300 motif–domain interfaces with AlphaFold and validated deleterious IDR variants that AlphaMissense called benign. Directly relevant to the ANKRD11 negative result: monomer pLDDT/contact graphs miss motif-mediated interactions in IDRs.
- **AlphaFold DB length limit [V, AFDB download page + Tunyasuvunakool 2021].** Human proteins > 2,700 aa are distributed only as overlapping 1,400-aa fragments (F1, F2, …), and only in the human-proteome bulk archive. ANKRD11 (2,663 aa) has a single model. **NF1 (2,839 aa), FBN1 (2,871 aa), ITPR1 (2,758 aa canonical) and KMT2D do not**: a multi-gene pipeline must handle fragments or exclude them.

**Gap statement for (b).** I found no published study that, across many genes, builds patient-level AlphaFold-derived structural descriptors (pLDDT region, 3D neighbourhood, network propagation) and tests them against individual HPO profiles with baselines for variant class, position and AlphaMissense. The closest are:
- Badonyi & Marsh 2025: multi-gene, structure → mechanism → phenotype, at disease level;
- Brünger et al. 2025: 3D position → phenotype, one gene family;
- GPSEA: multi-cohort individual-level HPO tests, no structure.

This rests on PubMed and web searches only (queries: AlphaFold/3D clustering/phenotype/HPO/within-gene). Absence of evidence is not proof; a final check of bioRxiv/medRxiv before submission is advised.

### 2c. Random walk with restart and network propagation on residue contact networks

- **Chennubhotla C, Bahar I.** Markov propagation of allosteric effects in biomolecular systems: application to GroEL–GroES. *Mol Syst Biol* 2:36, 2006. doi:10.1038/msb4100075. PMID 16820777 **[V title]**. A Markov random walk on a residue contact network (hitting and commute times). **The direct methodological precedent for our RWR; cite it.**
- **Chennubhotla C, Bahar I.** Signal propagation in proteins and relation to equilibrium fluctuations. *PLoS Comput Biol* 3(9):e172, 2007 **[UNVERIFIED — from memory; check before citing]**.
- **del Sol A, Fujihashi H, Amoros D, Nussinov R.** Residues crucial for maintaining short paths in network communication mediate signaling in proteins. *Mol Syst Biol* 2:2006.0019, 2006. doi:10.1038/msb4100063. PMID 16738564 **[V title]**. Shortest-path centrality on residue networks, the precedent for our "communication backbone".
- **Atilgan AR, Akan P, Baysal C.** Small-world communication of residues and significance for protein dynamics. *Biophys J* 86(1):85–91, 2004. doi:10.1016/S0006-3495(04)74086-2. PMID 14695252 **[V title]**.
- **Amor BR, Schaub MT, Yaliraki SN, Barahona M.** Prediction of allosteric sites and mediating interactions through bond-to-bond propensities. *Nat Commun* 7:12477, 2016. doi:10.1038/ncomms12477. PMID 27561351 **[V title]**.
- **Wang J, Jain A, McDonald LR, et al. (Dokholyan lab).** Mapping allosteric communications within individual proteins. *Nat Commun* 11:3862, 2020. doi:10.1038/s41467-020-17618-2. PMID 32737291 **[V title]**. Ohm, a perturbation-propagation method on contact networks.
- **Ponzoni L, Bahar I.** Structural dynamics is a determinant of the functional significance of missense variants. *PNAS* 115(16):4164–4169, 2018. doi:10.1073/pnas.1715896115. PMID 29610305 **[V title]**; and **Ponzoni L, Peñaherrera DA, Oltvai ZN, Bahar I.** Rhapsody: predicting the pathogenicity of human missense variants. *Bioinformatics* 36(10):3084–3092, 2020. doi:10.1093/bioinformatics/btaa127. PMID 32101277 **[V title]**. Elastic-network (ENM) features for variant pathogenicity, the closest "network dynamics → variant effect" work.
- **Cowen L, Ideker T, Raphael BJ, Sharan R.** Network propagation: a universal amplifier of genetic associations. *Nat Rev Genet* 18(9):551–562, 2017. doi:10.1038/nrg.2017.38. PMID 28607512 **[V]**. A review of RWR/diffusion. Our [9,10] cover gene-level RWR; cite this review too.
- **Smedley D, Köhler S, Czeschik JC, et al.** Walking the interactome for candidate prioritization in exome sequencing studies of Mendelian diseases. *Bioinformatics* 30(22):3215–3222, 2014. doi:10.1093/bioinformatics/btu508. PMID 25078397 **[V]** (ExomeWalker; gene-level RWR).

**Gap and risk for (c).** I found no published use of residue-level RWR/personalised PageRank on AlphaFold contact graphs as a *patient-phenotype* covariate; searches turned up only gene-level RWR and allostery methods. A reviewer will likely argue two things:
- Personalised PageRank on contact graphs is a known construct (Chennubhotla & Bahar 2006).
- With α = 0.3 the stationary mass is mostly within a few hops, so the fingerprint is roughly a smoothed local pLDDT class. The draft's own Fig. 2b and its parameter-robustness result (ρ ≥ 0.986) support this.

Preempt both by benchmarking against "seed pLDDT", "k-hop mean pLDDT" and "UniProt domain" baselines.

### 2d. KBG / ANKRD11 genotype–phenotype literature

- **Domains [V].**
  - Zhang A, Li CW, Chen JD. Characterization of transcriptional regulatory domains of ankyrin repeat cofactor-1. *Biochem Biophys Res Commun* 358(4):1034–1040, 2007. doi:10.1016/j.bbrc.2007.05.017. PMID 17521611. Defines **RD1 = residues 318–611**, **RD2 = 2369–2663**, and a minimal **activation domain (AD) = 2076–2145** (70 aa).
  - Awamleh et al. 2023 instead quote the AD as 1851–2145 and state that exon 9 encodes RD1, the AD and the start of RD2 **[V]**. Use Zhang's minimal AD (2076–2145) and note the broader definition.
  - UniProt Q6UB99: ANK repeats 167–292; region 2369–2663 "important for protein degradation"; long disordered regions throughout **[V]**.
  - Original identification as a p160 coactivator-interacting protein: Zhang A, et al. *J Biol Chem* 279(32):33799–33805, 2004. doi:10.1074/jbc.M403997200. PMID 15184363 **[V title]**.
- **Sirmaci et al. 2011 [ref 1] [V].** Exome sequencing identified a splice variant c.7570-1G>C (p.Glu2524_Lys2525del), segregating in a family, plus de novo truncating variants. No genotype–phenotype analysis.
- **Walz K, Cohen D, Neilsen PM, et al.** Characterization of ANKRD11 mutations in humans and mice related to KBG syndrome. *Hum Genet* 134(2):181–190, 2015. doi:10.1007/s00439-014-1509-2. PMID 25413698 **[V]**. The C-terminus is required for ANKRD11 degradation (D-box motifs). All 11 variants studied, including the Yoda mouse, affect C-terminal regions, and the mutant proteins accumulate aberrantly. They propose disturbed protein abundance.
- **Barbaric I, Perry MJ, Dear TN, et al.** An ENU-induced mutation in the Ankrd11 gene results in an osteopenia-like phenotype in the mouse mutant Yoda. *Physiol Genomics* 32(3):311–321, 2008. doi:10.1152/physiolgenomics.00116.2007. PMID 17986521 **[V title]**. I did not verify the exact Yoda residue (often quoted as E2502K).
- **Low et al. 2016 [ref 2] [V].** 32 patients; all 20 variants truncating; "**no conclusive phenotype-genotype correlations**".
- **Goldenberg A, Riccardi F, Tessier A, et al.** Clinical and molecular findings in 39 patients with KBG syndrome caused by deletion or mutation of ANKRD11. *Am J Med Genet A* 170(11):2847–2859, 2016. doi:10.1002/ajmg.a.37878. PMID 27605097 **[V]**. 19 16q24.3 deletions and 20 point variants. Deletions never extend toward the 16q telomere. The abstract reports no formal genotype–phenotype correlation.
- **Ockeloen CW, Willemsen MH, de Munnik S, et al.** Further delineation of the KBG syndrome phenotype caused by ANKRD11 aberrations. *Eur J Hum Genet* 23(9):1176–1185, 2015. doi:10.1038/ejhg.2014.253. PMID 25424714 **[V title]**; erratum *Eur J Hum Genet* 23(9):1270 (PMID 26269249).
- **Martinez-Cayuelas et al. 2023 [ref 3] [V].** See §1, point 1. It is the source of most phenopacket-store ANKRD11 phenopackets: the GPSEA notebook says the store parsed its supplementary tables with pyphetools **[V]**.
- **de Boer E, Ockeloen CW, Kampen RA, et al.** Missense variants in ANKRD11 cause KBG syndrome by impairment of stability or transcriptional activity of the encoded protein. *Genet Med* 24(10):2051–2064, 2022. doi:10.1016/j.gim.2022.06.007. PMID 35833929 **[V]**. A second PubMed record with the same title, Genet Med 25(11):100962, 2023 (PMID 37658852), is probably a correction or republication **[UNVERIFIED type]**.
  - 20 missense variants in 29 individuals; missense variants **cluster significantly in RD2 (C-terminus)**.
  - 6 of 10 tested reduced protein stability; one abolished transcriptional activity.
  - Conclusion: LoF/haploinsufficiency, with a KBG phenotype similar to PTVs.
  - Awamleh et al., citing this paper, add that "while the transcript can escape NMD, the protein is unstable" for these C-terminal variants **[V via Awamleh; I did not read de Boer's full text]**.
- **Awamleh Z, Choufani S, Cytrynbaum C, et al.** ANKRD11 pathogenic variants and 16q24.3 microdeletions share an altered DNA methylation signature in patients with KBG syndrome. *Hum Mol Genet* 32(9):1429–1438, 2023. doi:10.1093/hmg/ddac289. PMID 36440975 **[V]**. Exon-9 PTVs and 16q24.3 deletions share one episignature, which supports a single LoF mechanism. It notes that a dominant-negative effect has not been ruled out for C-terminal RD2 variants.
- **Parenti I, Mallozzi MB, Hüning I, et al.** ANKRD11 variants: KBG syndrome and beyond. *Clin Genet* 100(2):187–200, 2021. doi:10.1111/cge.13977. PMID 33955014 **[V title]**. A review. A search summary said it describes NMD escape in patient cell lines **[UNVERIFIED — check full text]**.
- **Su S, Ma J, Zhang Q, et al.** Insights into ANKRD11-related epilepsy from 163 people. *Epilepsia* 67(9):4922–4938, 2026. doi:10.1002/epi.70351. PMID 42322264 **[V]**. A pooled cohort in which missense variants cluster in RD2 and truncating variants outside it. Missense carriers tended toward later seizure onset, normal MRI and better seizure control. Results are descriptive, presented as "may be associated".
- **NMD conclusion.** The published consensus is that exon-9 PTVs act through haploinsufficiency (NMD); exon 9 is not the last exon. C-terminal RD2 missense variants and late PTVs may yield unstable protein. Episignature data unify PTVs and deletions. **No published study shows a robust within-KBG structural or positional phenotype correlation that survives multiple-testing correction.** Our negative result is consistent with GPSEA (0 significant FETs) and with Low 2016.
- **Reviewer "already done" risk for KBG.** Low. GPSEA already reported the KBG negative FET result in a supplement, so state explicitly that we extend it with continuous structural covariates. The very low missense count (7 in the store) means the RD2 missense cluster (de Boer 2022) cannot be tested here; say so.

---

## 3. Positive- and negative-control genes (phenopacket-store v0.1.27, cohort sizes from `/home/user/data_ext/cohort_sizes.txt`)

The variant-class counts come from my own parse of the local store's `hgvs.p` fields (`/home/user/data_ext/ppkt/...`), so they are approximate. "GPSEA" means the result in the gpsea-cs supplement (store v0.1.24/25; BH-adjusted p after IF-HPO). Each row gives n (v0.1.27) and the variant mix.

**General caveat.** In multi-disease cohorts (LMNA, FBN1, ITPR1, SPTAN1, PTPN11, SETD2, TGFBR1, SAMD9L), the diagnosis label is near-deterministic from genotype. HPO differences then partly re-encode disease labels assigned by curators. That is fine for positive controls, but say so.

### Tier 1: clean, published, structure- or position-driven missense signal, replicated by GPSEA

| Rank | Gene (n; mix) | Published correlation (verified) | GPSEA replication | Notes for our pipeline |
|---|---|---|---|---|
| 1 | **TBX5** (150; 80 missense, 59 PTV) | Basson CT, et al. *PNAS* 96(6):2919–2924, 1999. doi:10.1073/pnas.96.6.2919. PMID 10077612 **[V]**: within the T-box, **G80R → cardiac-predominant**; **R237Q/W → severe upper-limb, milder cardiac**. These residues contact different DNA grooves on the Xbra structure. Null alleles affect both organs. | **Yes**, R237Q vs other: VSD 0/17 vs 61/73 (adj p 1.8×10⁻⁹); upper-limb phocomelia 7/22 vs 3/131 (adj p 7×10⁻⁴). Missense vs other: VSD 31/60 vs 30/30. G80R vs other: 0 of 28 tests significant. | Two residues in the same folded domain with opposite phenotypes, a stringent test of 3D resolution. 518 aa, one AF model. |
| 2 | **NF1** (405; 295 missense, 36 in-frame, 58 SV) | Rojnueangnit K, et al. *Hum Mutat* 36(11):1052–1063, 2015. doi:10.1002/humu.22832. PMID 26178382 **[V]**: **p.Arg1809** → no cutaneous/plexiform neurofibromas, Noonan-like features, pulmonic stenosis. Koczkowska M, et al. *Am J Hum Genet* 102(1):69–87, 2018. doi:10.1016/j.ajhg.2017.12.001. PMID 29290338 **[V]**: **codons 844–848** (CSRD) → more plexiform/spinal neurofibromas, optic pathway glioma, malignancy. Upadhyaya M, et al. *Am J Hum Genet* 80(1):140–151, 2007. doi:10.1086/510781. PMID 17160901 **[V]**: **p.Met992del** (in-frame) → no cutaneous neurofibromas. | **Yes**, p.Arg1830 (NP_001035957.1 numbering) vs other: neurofibroma 5/133 vs 153/275 (adj p 6.7×10⁻²⁷), plexiform 0/135 vs 68/255, optic glioma 0/135 vs 39/217. Leu847Pro: 9 significant terms. Met992del: neurofibroma 0/44 vs 158/364. | **Isoform numbering:** GPSEA uses NP_001035957.1 (2,839 aa, the UniProt P21359 canonical length), so Arg1830 is very likely the legacy Arg1809 (exon 23a adds 21 aa) **[UNVERIFIED mapping — check]**. **AlphaFold:** 2,839 aa > 2,700, so only fragment models exist. Met992del is in-frame, not missense. |
| 3 | **LMNA** (259; 225 missense; FPLD2 127, DCM1A 61, EDMD2 41, L-CMD 15, HGPS 15) | Shackleton S, et al. *Nat Genet* 24(2):153–156, 2000. doi:10.1038/72807. PMID 10655060 **[V]**: site-specific missense variants cause FPLD vs EDMD/DCM. Krimm I, et al. *Structure* 10(6):811–823, 2002. doi:10.1016/S0969-2126(02)00777-3. PMID 12057196 **[V]**: Ig-fold (430–545) structure; FPLD hotspot R482W/Q vs EDMD R453W have different molecular consequences. Reports that FPLD variants form a surface patch while muscle variants affect the core **[UNVERIFIED wording — read full text]**. | **Yes**, "Upstream of NLS" (rod) vs other: DCM 35/70 vs 5/102 (adj p 8×10⁻¹¹); lipodystrophy 10/93 vs 132/158 (adj p 6.6×10⁻³⁰); AV block 17/36 vs 8/116. Missense vs other: truncal fat loss 104/104 vs 4/11. | The best 3D test: an Ig-fold surface patch (FPLD) vs core or rod (striated muscle). Exclude HGPS c.1824C>T (p.Gly608=, a cryptic splice) from missense structural features. 664 aa, one model. |
| 4 | **ITPR1** (170; 132 missense; SCA29 104, Gillespie 39, SCA15 27) | McEntagart M, et al. *Am J Hum Genet* 98(5):981–992, 2016. doi:10.1016/j.ajhg.2016.03.018. PMID 27108798 **[V]**: Gillespie syndrome (aniridia) is caused by a **restricted repertoire** (Glu2094, Gly2539 missense; Lys2596 in-frame del); dominant-negative. Gerber S, et al. *Am J Hum Genet* 98(5):971–980, 2016. doi:10.1016/j.ajhg.2016.03.004. PMID 27108797 **[V title]**. Tolonen JP, et al. *Mov Disord* 39(1):141–151, 2024. doi:10.1002/mds.29651. PMID 37964426 **[V]**: SCA29/GLSP missense variants cluster in the IP3-binding domain, the CA8-binding region and the TM channel; within-SCA29 genotype–phenotype associations are "highly variable"; standardised transcript numbering. | **Yes**, GS hotspot vs other: aniridia 14/20 vs 6/48 (adj p 1.6×10⁻⁴). IP3-binding vs other: aniridia 0/14 vs 20/54, more developmental delay. SV deletion (SCA15) vs other: no neurodevelopmental delay (0/19 vs 99/100). | **Numbering:** isoforms differ, so Gly2539 vs Gly2554 etc. Map to UniProt Q14643 carefully. Canonical 2,758 aa > 2,700, so **fragment models only**. |
| 5 | **FBN1** (143; 119 missense; MFS 50, ectopia lentis 44, geleophysic 19, acromicric 13, MFLS 9) | Faivre L, et al. *Am J Hum Genet* 81(3):454–466, 2007. doi:10.1086/520125. PMID 17701892 **[V]**: cysteine-substituting/creating missense → **higher ectopia lentis**; PTC → more severe skeletal/skin; **exons 24–32** → severe, early, complete phenotype. Le Goff C, et al. *Am J Hum Genet* 89(1):7–14, 2011. doi:10.1016/j.ajhg.2011.05.012. PMID 21683322 **[V]**: **TB5 domain (exons 41–42)** → acromicric/geleophysic dysplasia (short stature). Arnaud P, et al. *J Med Genet* 61(5):469–476, 2024 (PMID 38458756) **[V title]**: TB5 variants can also give MFS. | **Yes**, TB domain vs cbEGF: short stature 23/39 vs 0/24 (adj p 8.8×10⁻⁶), ectopia lentis 9/18 vs 48/59. Exon 37 vs other: stiff skin 8/9 vs 0/50 (adj p 8.5×10⁻⁸). Missense vs other: arachnodactyly and thoracic aortic aneurysm lower in missense. | 2,871 aa > 2,700, so fragment models only, and a long cbEGF rod in which AlphaFold inter-domain contacts are unreliable. The signal is domain-level, so it tests whether 3D beats a UniProt-domain baseline. |
| 6 | **SPTAN1** (85; 34 missense, 22 in-frame, 24 PTV; SPG91 28, DEE5 22, DD±epilepsy 21, dHMN11 14) | Van de Vondel L, et al. *Mov Disord* 37(6):1175–1186, 2022. doi:10.1002/mds.28959. PMID 35150594 **[V]**: recurrent **p.Arg19Trp** in 15 spastic-paraplegia patients; ataxia variants. Beijer D, et al. *Brain* 142(9):2605–2616, 2019. doi:10.1093/brain/awz216. PMID 31332438 **[V title]**: nonsense → juvenile HMN. Syrbe S, et al. *Brain* 140(9):2322–2336, 2017. doi:10.1093/brain/awx195. PMID 29050398 **[V]**: variants in the **α/β heterodimerisation domain** (C-terminal repeats) → severe DEE with atrophy; variants outside it are milder. Saitsu H, et al. *Am J Hum Genet* 86(6):881–891, 2010 (PMID 20493457) **[V title]**. | **Yes**, very strong. R19W vs other: spastic paraplegia 21/21 vs 0/20 (adj p 9×10⁻¹¹). Truncating vs other: motor axonal neuropathy 12/12 vs 0/23 (adj p 3×10⁻⁸). Missense vs other: spasticity up, seizures/ID/infantile spasms down. | Combines a class signal (PTV → neuropathy) with a position signal (R19 vs dimerisation domain → DEE). 2,472 aa, one model. In-frame indels need a position rule. |
| 7 | **SETD2** (46; 27 missense; Luscan-Lumish 29, Rabin-Pappas 14) | Rabin R, et al. Genotype-phenotype correlation at codon 1740 of SETD2. *Am J Med Genet A* 182(9):2037–2048, 2020. doi:10.1002/ajmg.a.61724. PMID 32710489 **[V]**: **Arg1740Trp** → microcephaly, profound ID, multiple anomalies. Arg1740Gln → intermediate. LoF → overgrowth/macrocephaly (LLS). | **Yes**, Arg1740Trp vs other: macrocephaly 0/11 vs 19/28, hypertelorism 11/11 vs 5/23, severe GDD 9/9 vs 0/12. Missense vs other: macrocephaly 4/24 vs 15/15. | A single-residue signal. Structure should help only if R1740's neighbourhood differs from other missense sites. A good test that position or AlphaMissense baselines can match the structural features. |
| 8 | **IKZF1** (81; 56 missense; all labelled CVID13) | Boutboul D, et al. *J Clin Invest* 128(7):3071–3087, 2018. doi:10.1172/JCI98164. PMID 29889099 **[V]**: de novo **N159 (DNA-binding ZF)** → dominant-negative combined immunodeficiency, Pneumocystis in all patients, one T-ALL. Kuehn HS, et al. *N Engl J Med* 374(11):1032–1043, 2016. doi:10.1056/NEJMoa1512234. PMID 26981933 **[V title]**: haploinsufficiency → B-cell loss / CVID-like. Hoshino A, et al. *Sci Immunol* 7(69):eabi7160, 2022 (PMID 35333544) **[V title]**: GoF variants → immune dysregulation. | **Yes, weaker.** N159S vs other: recurrent pneumonia 7/9 vs 7/39 (adj p 0.035); T-ALL 3/13 vs 0/69 (adj p 0.047). Missense vs other: 0 of 34 tests significant. | A mechanism (DN vs HI) tied to a 3D DNA-contact residue. A good small positive control. |

### Tier 2: published or GPSEA-supported, but driven by class, linear position or disease label rather than missense 3D context

| Gene (n) | Correlation | Evidence | Comment |
|---|---|---|---|
| **GLI3** (81; 66 PTV) | Johnston JJ, et al. *Am J Hum Genet* 76(4):609–622, 2005. doi:10.1086/429346. PMID 15739154 **[V]**: PTVs in the first third (nt 1–1997) → GCPS; second third (1998–3481) → **PHS**; 3′ third → GCPS. | GPSEA **yes**: mid-region fs/splice → Y-shaped metacarpals 11/23 vs 4/56, anal atresia; C-terminal third → macrocephaly, postaxial polydactyly. | A **position-only control**: structural features of PTV sites should *not* beat a linear-position baseline. If they do, suspect leakage. |
| **PTPN11** (68; NSML 31, NS1 25, metachondromatosis 12) | Tartaglia M, et al. *Am J Hum Genet* 70(6):1555–1563, 2002. doi:10.1086/340847. PMID 11992261 **[V]**: NS missense variants cluster at the **N-SH2/PTP interface** (GoF). Kontaridis MI, et al. *J Biol Chem* 281(10):6785–6792, 2006. doi:10.1074/jbc.M513068200. PMID 16377799 **[V]**: LEOPARD/NSML variants are **catalytically defective** (PTP domain), mechanistically distinct from NS. Sobreira NL, et al. *PLoS Genet* 6(6):e1000991, 2010 (PMID 20577567) **[V title]**: PTPN11 LoF in metachondromatosis. | GPSEA: only missense vs other (driven by metachondromatosis): hypertelorism 37/41 vs 0/12, pulmonic stenosis, webbed neck. **Tyr279Cys vs other: 0 of 8 significant.** NS-vs-NSML structure was not tested directly. | An excellent *a priori* 3D hypothesis (interface vs catalytic pocket → lentigines, HCM, deafness) that GPSEA did not test. If our features recover it, this is a strong selling point; it is not verified in this cohort. |
| **TGFBR1** (42; LDS1 24, MSSE 18) | MSSE = LoF (including truncating/extracellular-domain variants); LDS1 = kinase-domain missense **[UNVERIFIED citation: Goudie DR et al., Nat Genet 2011]**. | GPSEA **yes**: MSSE variants → self-healing epithelioma 18/19 vs 0/21; hypertelorism 0/18 vs 15/19. | Disease-label-driven. |
| **SMARCB1** (36; CSS3 22, rhabdoid predisposition 14) | Coffin–Siris variants are missense/in-frame near the C-terminus; rhabdoid predisposition variants are truncating/deletions (Kosho T, et al. *Am J Med Genet C* 166C(3):262–275, 2014. PMID 25168959 **[V title only]**). | GPSEA: SV vs other → rhabdoid/embryonal tumours (adj p ≤ 6×10⁻⁴). | Class/SV-driven. |
| **UMOD** (199; 155 missense, 43 in-frame) | GPSEA: **cysteine-affecting vs other → hyperuricemia 38/41 vs 26/48** (adj p 2.7×10⁻⁴); EGF-domain vs other → hyperuricemia 14/32 vs 50/57; Cys248Trp vs Gln316Pro → earlier stage-5 CKD (log-rank p 4.1×10⁻⁴). Olinger E, et al. *Kidney Int* 98(3):717–731, 2020. doi:10.1016/j.kint.2020.04.038. PMID 32450155 **[V]** characterises ADTKD-UMOD (gout-free survival), but I did not verify a position-specific claim in it. | GPSEA-derived; no verified published UMOD position → HPO paper. | Large n, nearly all missense: an ideal power test for residue-level features (cysteine/disulfide context). |
| **SATB2** (157; 49 missense, 70 PTV, 30 SV) | Zarate YA, et al. Mutation update for the SATB2 gene. *Hum Mutat* 40(8):1013–1029, 2019. doi:10.1002/humu.23771. PMID 31021519 **[V; abstract mentions "novel genotype-phenotype correlations", specifics not read]**. den Hoed J, et al. *J Med Genet* 61(11):1062–1067, 2024 (PMID 39327041) **[V title]**: p.Gly392 variants have variable functional and clinical effects. | GPSEA: missense vs other → cleft palate 11/49 vs 59/105 (adj p 0.002). | Class-driven. |
| **CHD8** (78; 18 missense, 50 PTV) | Dingemans AJM, et al. *Transl Psychiatry* 12:421, 2022. doi:10.1038/s41398-022-02189-1. PMID 36182950 **[V title/abstract start]**: HPO-based genotype–phenotype analysis. GPSEA reports missense → lower de Vries score. | GPSEA **yes** (Mann–Whitney p 9×10⁻⁴); 0 of 36 FETs significant. | Class/severity-score control. |
| **SAMD9L** (30; all missense; ataxia-pancytopenia 21, SCA49 9) | GPSEA: Arg986Cys vs Ser626Leu → pancytopenia, thrombocytopenia, neutropenia 7/9 vs 0/9 (adj p 0.008). I did not verify the published source. | GPSEA only. | Small; two-residue contrast. |
| **SMAD3** (56) | GPSEA: p.Arg287Trp → osteoarthritis 19/19 vs 7/19 (adj p 8.6×10⁻⁴). GPSEA states this is **not previously published**. | GPSEA only. | A single-variant (probably single-family) effect. Use as an "unexpected" control, not a positive one. |
| **PPP2R1A** (60; store has c. HGVS only, so map with VEP) | Lenaerts L, et al. *Genet Med* 23(2):352–362, 2021. doi:10.1038/s41436-020-00981-2. PMID 33106617 **[V]**: variants that disrupt **B55α binding** (HEAT-repeat interface) → profound ID, epilepsy, corpus callosum hypoplasia, ± microcephaly. Macrocephaly only without the B55α-binding deficit. | GPSEA: Arg182Trp vs other 0 of 19 significant; SV40-ST-binding region 0 of 19. | A strong *a priori* 3D-interface hypothesis that GPSEA did **not** replicate. A useful "hard" case. |

### Negative or null controls (literature and/or GPSEA report no within-gene genotype–phenotype correlation)

| Gene (n) | Evidence |
|---|---|
| **STXBP1** (462; 248 missense, 137 PTV) | Xian J, et al. Assessing the landscape of STXBP1-related disorders in 534 individuals. *Brain* 145(5):1668–1683, 2022. doi:10.1093/brain/awab327. PMID 35190816 **[V]**: hotspots (R406, R292, R551, P139, R190) are **not** associated with distinct syndromes, single features or overall similarity, because "baseline variability … is too high". **Caveat:** PTV/deletion carriers were about twice as likely to have West syndrome, so this is a **null control for position or structure but not for variant class**. GPSEA: missense, Arg406 and exon 14 vs other, 0 of 18 significant each. The largest cohort, so the best-powered null. |
| **ANKRD11** (333) | GPSEA: 0 FET results significant (SV, exon 9, c.1903_1907del); Low 2016: "no conclusive" correlations. Our own case study. |
| **POGZ** (117) | GPSEA: missense/frameshift vs other, 0 of 24 significant; severity score p = 0.43. |
| **EHMT1** (125) | GPSEA: missense, ANK, SET vs other, 0 of 62 significant each (one N-terminal-frameshift ADHD hit). |
| **KMT2D** (64; Kabuki 1) | GPSEA: missense, N-term, p.Glu5425Lys, p.Leu3542Pro, 0 significant. The distinct exon 38–39 missense disorder (Cuvertino S, et al. *Genet Med* 22(5):867–877, 2020. PMID 31949313 **[V title]**) is presumably not in this Kabuki-labelled cohort. Also 5,537 aa, so fragment models only. |
| **TRAF7** (45; all missense) | GPSEA: WD7 and Arg655Gln vs other, 0 of 40 significant. |
| **U2AF2** (48; missense) | GPSEA: R149W, R149/R150, RRM1 vs other, 0 of 177 significant. |
| **FBXO11** (56), **SON** (64), **ZNF462** (39), **KDM6B** (84) | GPSEA: 0 significant for missense vs other, the recurrent SON c.5753_5756del, ablation vs other, and KDM6B JmjC domain / N- vs C-terminal. These are mostly PTV cohorts. |
| **NIPBL** (60) | Gillis LA, et al. *Am J Hum Genet* 75(4):610–623, 2004. doi:10.1086/424698. PMID 15318302 **[V]**: only a *trend* toward milder disease with missense. GPSEA: missense and "disordered" vs other, 0 of 86 significant. A weak or null control. |
| **COL3A1** (41) | Literature: null (haploinsufficient) alleles → later complications and longer survival (Leistritz DF, et al. *Genet Med* 13(8):717–722, 2011. PMID 21637106 **[V title]**; Pepin MG, et al. *Genet Med* 16(12):881–888, 2014. doi:10.1038/gim.2014.72. PMID 24922459 **[V]**). GPSEA: missense, Gly-missense and triple-helix vs other, 0 of 20 significant. The published signal is survival, not HPO, so it is a null for HPO-level tests. |
| **GNAS** (103; PHP1a 80, PPHP 23) | The phenotype is determined by parent-of-origin imprinting, not by variant position **[UNVERIFIED — standard knowledge; cite a review]**. A conceptual negative control: identical alleles give different phenotypes. |

**Recommended control set for the paper:**
- **Positive:** TBX5, NF1, LMNA, ITPR1, FBN1, SPTAN1, SETD2, IKZF1, plus PTPN11 and PPP2R1A as *a priori* 3D hypotheses.
- **Position-only control:** GLI3.
- **Null:** STXBP1, POGZ, EHMT1, TRAF7, U2AF2, ANKRD11.

Note that NF1, FBN1, ITPR1 and KMT2D need AlphaFold fragment handling.

---

## 4. Venue suggestions

### (a) Systematic multi-gene study with an open tool

| Journal | Fit | Length limits |
|---|---|---|
| **HGG Advances** (Cell Press/ASHG, OA) | Phenopacket-store resource paper venue (Danis 2025); natural home | Limits **[UNVERIFIED]**: the author page returned 403. |
| **Bioinformatics** (OUP) | Methods plus open tool | Original Paper "up to 7 pages; approx. 5,000 words excluding figures"; Application Note "up to 4 pages; approx. 2,600 words or 2,000 words plus one figure" **[V, OUP author guidelines]**. |
| **PLOS Computational Biology** | Research, Methods or Software articles | Research articles "can be any length … no restrictions on word count, number of figures"; abstract ≤ 300 words; **Software articles "less than 3500 words"** **[V, PLOS submission guidelines]**. |
| **NAR Genomics and Bioinformatics** (OUP) | Standard or Methods papers; Application Notes for software | Application Notes "max 3000 words, 3 figures/tables or equivalent and 30 references"; no word limit stated for Standard/Methods papers; software must be OSI-licensed open source **[V, Scope & Criteria page]**. |
| **Genome Medicine** (BMC) or **European Journal of Human Genetics** (Nature) | Higher-impact clinical-genomics options | Limits **[UNVERIFIED]**: pages redirected or blocked. |
| **Am J Hum Genet** | Published GPSEA; ambitious | Only if the multi-gene result is clearly positive and novel. Limits **[UNVERIFIED]**. |

### (b) Short negative-result paper (KBG / ANKRD11)

| Journal | Fit | Length limits |
|---|---|---|
| **BMC Research Notes** | Explicitly accepts null results | Research Note: abstract 200 words; introduction + main text + limitations ≤ 2,000 words; "null results" listed among accepted types **[V via search snippet of the official research-note guideline page; direct fetch redirected]**. |
| **Bioinformatics Application Note** | Frame as a tool paper with a candid case study | ≤ 4 pages, ~2,600 words **[V]**. |
| **NAR Genomics and Bioinformatics Application Note** | Same framing | ≤ 3,000 words, 3 figures, 30 references **[V]**. |
| **Clinical Genetics** Short Report (Wiley) | Published the ANKRD11 PhenoScore and Parenti review, so topical | 2,500 words, 3 figures/tables, 20 references, 200-word unstructured abstract **[V only via search-result snippet; Wiley page returned 403, so re-check]**. |
| **PLOS ONE** / **F1000Research** | Soundness-only review | Limits not checked **[UNVERIFIED]**. |

---

## Additional references to add to the draft (minimum set)

1. GPSEA: Rekerle et al. 2026 *AJHG*.
2. Phenopacket-store: Danis et al. 2025 *HGG Adv*.
3. de Boer et al. 2022 *Genet Med* (missense, RD2).
4. Zhang et al. 2007 *BBRC* (RD1/RD2/AD).
5. Walz et al. 2015 *Hum Genet*.
6. Awamleh et al. 2023 *HMG*.
7. Cheng et al. 2023 *Science* (AlphaMissense).
8. Chennubhotla & Bahar 2006 *Mol Syst Biol*.
9. Cowen et al. 2017 *Nat Rev Genet*.
10. Sivley et al. 2018 *AJHG*.
11. Badonyi & Marsh 2025 *Nat Commun*.
12. Brünger et al. 2025 *Genome Biol*.
13. Ruff & Pappu 2021 *JMB* or Hubrich et al. 2026 *NSMB* (IDR caveats).
