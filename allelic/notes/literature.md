# Literature notes: allelic heterogeneity study (PhenoBridge follow-up)

Compiled 2026-11-03 for the planned study of genes in which different missense variants cause **different named diseases**. The unit of analysis is the distinct variant. The three questions:

- **Q1.** Do different-disease variants differ in structural context beyond position, domain and AlphaMissense?
- **Q2.** Are same-disease variants closer in 3D and in RWR "network proximity" than different-disease variants?
- **Q3.** Within mixed-inheritance genes, do AD and AR missense variants differ structurally?

**How sources were checked.**
- Bibliographic data were checked against PubMed E-utilities (esearch/esummary/efetch), Europe PMC REST (including full-text XML), Crossref and publisher pages.
- OMIM disease IDs and inheritance modes come from the HPO annotation files `genes_to_disease.txt` and `phenotype.hpoa` (latest GitHub release, downloaded 2026-11-03).
- Protein lengths come from UniProt REST (reviewed human entries).

**Marking.**
- **[V]**: bibliographic record verified, and I read the abstract or full text for the claim made.
- **[V title]**: record verified, but the specific claim was not checked in the text.
- **[UNVERIFIED]**: not checked against a primary source.

Citations already verified in `paper/notes/literature.md` are reused and marked as such, not re-derived.

---

## 0. Bottom line (novelty verdict)

| Q | Verdict | Closest prior work | "Already done" risk |
|---|---|---|---|
| **Q1** | **Concept established; the formal, calibrated version has not been done.** | Wang 2012 (3D interactome: "disease specificity for different mutations of the same gene can be explained by their location within an interface"); Zhong 2009 and Sahni 2015 (edgetic alleles → distinct disorders); **Badonyi & Marsh 2025** (ClinVar P/LP missense + OMIM phenotype labels + AlphaFold, phenotype-level mLOF = ΔΔG + 3D clustering, 2,837 phenotypes, intragenic mechanistic heterogeneity). | **High.** A reviewer will cite Wang 2012 and Badonyi 2025. Novelty must rest on four things: (i) an explicit within-gene test with a variant-level null; (ii) incremental value over sequence position, UniProt domain and AlphaMissense; (iii) pLDDT-region/RWR/degree/burial descriptors; (iv) calibration on positive and null controls. None of the prior works does (i)–(ii). Badonyi ranks phenotype pairs by mLOF difference (top 5%) without a within-gene significance test or position baseline. Wang 2012 uses interface membership, not continuous structure, and only PDB-resolved interactions. |
| **Q2** | **Not done in the proposed form.** I found no multi-gene study testing whether same-disease variant pairs are closer than different-disease pairs, beyond sequence distance, with a Mantel-type or pair-permutation test. **I found no use of RWR / diffusion "network proximity" on residue contact graphs to compare variants by disease label.** | Brünger 2025: Kendall correlation between phenotype-specific 3D neighbour-count profiles, one gene family (10 VGSC genes), patient counts, Bonferroni, no position control. Ng 2025 (3D-CLUMP): clustering in NDD vs cancer tested separately per class. Guo 2013 and Wang 2012: same interface → same disease. Homburger 2016 and Lemmers 2019: single genes. Nason 2026 (medRxiv): 3DNT case/control neighbourhoods; ClinVar P/LP vs B/LB only, no condition labels. | **Moderate for the 3D-distance part, low for RWR.** Expect three objections: (a) "Wang/Guo/Brünger already showed same location → same disease"; (b) "RWR proximity at α = 0.3 is a smoothed graph distance, so show it adds beyond Euclidean Cα distance and shortest-path hops"; (c) partial Mantel tests are anti-conservative under spatial autocorrelation (Guillot & Rousset 2013). |
| **Q3** | **Largely done at gene or phenotype level; a within-gene, variant-level, position-controlled contrast has not been done.** | Turner 2015 (AD genes more 1D-clustered than AR genes); Gerasimavicius 2022 (AD/non-LoF milder ΔΔG and more 3D-clustered; DN enriched at interfaces; **mixed AD/AR genes were excluded**); Guo 2013 (inheritance mode × interface); **Badonyi 2025** (mixed-inheritance genes: AD phenotypes mostly non-LoF, 9.9% LoF; AR phenotypes in mixed genes have lower ΔΔG, i.e. hypomorphic); Nason 2026 (AD genes have more significant 3D neighbourhoods); MOI predictors DOMINO (gene level), MAVERICK and MOI-Pred (variant level). | **High.** Q3 is publishable only as a variant-level, within-gene test with position/domain/AlphaMissense baselines and pLDDT/RWR features. It also needs explicit "null" mixed genes where AD vs AR reflects zygosity or dosage rather than mechanism (GCK MODY2 vs PNDM1; CASR FHH1 vs NSHPT). |

**Overall framing.** The defensible novelty is **not** "structure explains allelic disease", which is known. It is a *calibrated, pre-registered, multi-gene test of whether AlphaFold-derived context and network proximity add information beyond sequence position, domain and AlphaMissense*. The test treats distinct variants as units and uses explicit positive and null controls. This is the same angle that distinguished the first PhenoBridge paper.

---

## 1. Closest prior work (Task 1)

### 1a. Interactome / "edgetic" models: different interfaces → different diseases

- **Zhong Q, Simonis N, Li QR, et al.** Edgetic perturbation models of human inherited disorders. *Mol Syst Biol* 5:321, 2009. doi:10.1038/msb.2009.80. PMID 19888216 **[V]**
  - *What they did:* Global analysis of ~50,000 Mendelian mutations, plus experimental profiling of mutant alleles. "A large fraction of cases in which a single gene is linked to multiple disorders can be modeled by distinguishing edgetic network perturbations."
  - *Gap for us:* Gene-interaction level, not residue-level structural context. No test against a position baseline.

- **Wang X, Wei X, Thijssen B, Das J, Lipkin SM, Yu H.** Three-dimensional reconstruction of protein networks provides insight into human genetic disease. *Nat Biotechnol* 30(2):159–164, 2012. doi:10.1038/nbt.2106. PMID 22252508 **[V]**
  - *What they did:* 3,949 genes, 62,663 mutations, 3,453 disorders, 4,222 PPIs with atomic-resolution interfaces. In-frame mutations are enriched on interfaces, and "the disease specificity for different mutations of the same gene can be explained by their location within an interface."
  - *This is the most direct "already done" citation for Q1/Q2.*
  - *Gap for us:*
    - limited to PDB/homology-resolved interfaces (pre-AlphaFold);
    - categorical interface assignment;
    - individual-level/mutation-level counting, with no variant-level permutation and no position or AlphaMissense baseline;
    - no monomer-internal context such as pLDDT regions or contact-graph diffusion.

- **Guo Y, Wei X, Das J, et al.** Dissecting disease inheritance modes in a three-dimensional protein network challenges the "guilt-by-association" principle. *Am J Hum Genet* 93(1):78–89, 2013. doi:10.1016/j.ajhg.2013.05.022. PMID 23791107 **[V]**
  - *What they did:* 82,833 mutations classified by inheritance (recessive vs dominant) and type (in-frame vs truncating) on a 3D interactome. Recessive mutations at the same interface tend to cause the same disease; this does **not** hold for dominant mutations.
  - *Gap for us:* Directly relevant to Q2 and Q3, but at the level of interface co-membership. No residue-level distance or RWR, and no within-gene AD-vs-AR structural contrast with position control.

- **Sahni N, Yi S, Taipale M, et al.** Widespread macromolecular interaction perturbations in human genetic disorders. *Cell* 161(3):647–660, 2015. doi:10.1016/j.cell.2015.04.013. PMID 25910212 **[V]**
  - *What they did:* Experimental interaction profiling of thousands of missense alleles. "Different mutations in the same gene leading to different interaction profiles often result in distinct disease phenotypes."
  - *Gap for us:* Experimental, and not a structural-feature test.

- **Mosca R, Tenorio-Laranga J, Olivella R, et al.** dSysMap: exploring the edgetic role of disease mutations. *Nat Methods* 12(3):167–168, 2015. doi:10.1038/nmeth.3289. PMID 25719824 **[V title]**. A resource that maps disease mutations onto structurally resolved interactomes.

- **Meyer MJ, Beltrán JF, Liang S, et al.** Interactome INSIDER: a structural interactome browser for genomic studies. *Nat Methods* 15(2):107–114, 2018. doi:10.1038/nmeth.4540. PMID 29355848 **[V title]**. A source of predicted interfaces for the interface covariate.

- **Cheng F, Zhao J, Wang Y, et al.** Comprehensive characterization of protein–protein interactions perturbed by disease mutations. *Nat Genet* 53(3):342–353, 2021. doi:10.1038/s41588-020-00774-y. PMID 33558758 **[V]**
  - *What they did:* Germline disease variants and somatic mutations are enriched at PPI interfaces; mostly an oncology application.
  - *Gap for us:* No within-gene disease-label analysis.

### 1b. Molecular mechanism (LoF vs GoF/DN) and structure: Marsh lab and related

- **Gerasimavicius L, Livesey BJ, Marsh JA.** Loss-of-function, gain-of-function and dominant-negative mutations have profoundly different effects on protein structure. *Nat Commun* 13:3895, 2022. doi:10.1038/s41467-022-31686-6. PMID 35794153 **[V, full text]**
  - *What they did:*
    - AR mutations are more destabilising than AD; non-LoF mutations are milder.
    - DN mutations are enriched at interfaces.
    - Non-LoF mutations cluster in 3D. This introduced **EDC (Extent of Disease Clustering)**: the ratio of mean distance from non-disease to disease residues over disease to disease residues.
  - **Mixed-inheritance genes were explicitly excluded** from the AD/AR analysis ("excluding those with mixed inheritance").
  - In HGMD genes with both GoF and LoF variants, GoF variants were milder than LoF variants *within the same gene* (Supplementary Fig. 6b).
  - *Gap for us:* Gene-level mechanism labels. The within-gene GoF/LoF contrast uses only ΔΔG, with no position control and no disease-label analysis.

- **Badonyi M, Marsh JA.** Proteome-scale prediction of molecular mechanisms underlying dominant genetic diseases. *PLoS One* 19(8):e0307312, 2024. doi:10.1371/journal.pone.0307312. PMID 39172982 **[V]**. Gene-level DN/GOF/LOF classifiers.

- **Badonyi M, Marsh JA.** Prevalence of loss-of-function, gain-of-function and dominant-negative mechanisms across genetic disease phenotypes. *Nat Commun* 16:8392, 2025. doi:10.1038/s41467-025-63234-3. PMID 40998763. PMCID PMC12462468 **[V, full text read]**. **The closest work to Q1 and Q3.**
  - *Data:*
    - ClinVar P/LP missense (accessed 10 Sep 2024), mapped by VEP to UniProt;
    - OMIM phenotype cross-references for 71.8% of 63,920 variants;
    - AlphaFold models;
    - inheritance from OMIM/HPO.
  - *Method:* A **phenotype-level** mLOF score, computed by kernel density from ΔΔGrank plus EDC, over 2,837 phenotypes in 1,979 genes.
  - *Within-gene results:* "43% of dominant and 49% of mixed-inheritance genes harbour both LoF and non-LoF mechanisms."
    - Table 4 ranks phenotype pairs with the largest mLOF difference (top 5%). Examples: SMCHD1 FSHD2 vs BAMS, KRAS, TP63 RHS vs SHFM4, BRAF, MTOR, AARS1 DEE29 vs CMT2N, CAV3, SCN4A, ACTA1, ADCY5, LHCGR, NPR2, TWNK, POLR3B.
    - Fig. 4 shows the two phenotypes' variants on AlphaFold models.
  - *Mixed-inheritance genes:*
    - AD phenotypes are 9.9% LoF, 53.6% GoF and 36.5% DN.
    - AR phenotypes in mixed genes have lower ΔΔGrank than exclusively AR genes (P = 1.6 × 10⁻³), consistent with hypomorphic alleles.
    - Combining mLOF with HPO semantic similarity of AD–AR phenotype pairs (217 pairs, 103 genes) prioritises DN mechanisms, e.g. CLCN7 OPTA2 vs OPTB4.
  - *Exact gap for us:*
    1. There is no within-gene hypothesis test of a difference between phenotype-specific variant sets, and no permutation null.
    2. There are no baselines for sequence position, UniProt domain or AlphaMissense. ΔΔG and EDC may partly re-encode position or domain.
    3. Inference is at phenotype level, not variant level. There is no pairwise "same-disease closer" analysis (Q2).
    4. pLDDT regions, contact degree, burial and network propagation are not used.
    5. Labels are OMIM only, without phenopacket-store and without harmonisation across ClinVar submitters.

  **We must cite it prominently, reuse its ClinVar→OMIM mapping logic, and ideally reproduce its Table 4 pairs as positive controls.** Its data (OSF doi:10.17605/OSF.IO/AH2UC) include per-phenotype mLOF scores and could serve as a direct benchmark.

- **Backwell L, Marsh JA.** Diverse molecular mechanisms underlying pathogenic protein mutations: beyond the loss-of-function paradigm. *Annu Rev Genomics Hum Genet* 23:475–498, 2022. doi:10.1146/annurev-genom-111221-103208. PMID 35395171 **[V title]**. Review.

- **Bergendahl LT, Gerasimavicius L, Miles J, et al.** The role of protein complexes in human genetic disease. *Protein Sci* 28(8):1400–1411, 2019. doi:10.1002/pro.3667. PMID 31219644 **[V title]**.

- **Hernández Trapero R, Badonyi M, Gerasimavicius L, et al.** Complementary roles of structure and variant effect predictors in RyR1 clinical interpretation. *Hum Mutat* 2025:1834898, 2025. doi:10.1155/humu/1834898. PMID 41080992 **[V]**. Introduces **SPDV (Spatial Proximity to Disease Variants)** for RYR1, a 3D-clustering pathogenicity metric. Single gene; not disease-label specific.

- **Sevim Bayrak C, et al. 2021 AJHG** (GoF vs LoF features) and **Stein D, Kars ME, Wu Y, et al.** Genome-wide prediction of pathogenic gain- and loss-of-function variants from ensemble learning of a diverse feature set (LoGoFunc). *Genome Med* 15:103, 2023. doi:10.1186/s13073-023-01261-9. PMID 38037155 **[V title]**. Variant-level GoF/LoF predictors. Use them as an *additional baseline* for Q3: does structure add beyond LoGoFunc?

### 1c. Spatial / 3D clustering methods (mostly genotype-only)

- **Turner TN, Douville C, Kim D, et al.** Proteins linked to autosomal dominant and autosomal recessive disorders harbor characteristic rare missense mutation distribution patterns. *Hum Mol Genet* 24(21):5995–6002, 2015. doi:10.1093/hmg/ddv309. PMID 26246501 **[V]**. CLUMP (1D clustering by mutation position): AD-disease proteins show more clustering than AR-disease proteins (Wilcoxon P = 5.7 × 10⁻⁴). Gene-level and 1D only.
- **Sivley RM, Dou X, Meiler J, et al.** Comprehensive analysis of constraint on the spatial distribution of missense variants in human protein structures. *Am J Hum Genet* 102(3):415–426, 2018. doi:10.1016/j.ajhg.2018.01.017. PMID 29455857 **[V, methods read on PMC]**.
  - **Univariate Ripley's K** per protein; the null comes from permuting variant positions within the structure (100,000 permutations).
  - The area between observed and null K curves is converted to a Z-score; FDR is 10% (Storey q-values).
  - Two variant classes are compared only by correlating their per-protein Z-scores. **There is no bivariate (cross-K / D) test.** For Q2 we need a two-label method (see §4d).
- **Sivley RM, Sheehan JH, Kropski JA, et al.** Three-dimensional spatial analysis of missense variants in RTEL1 identifies pathogenic variants in patients with Familial Interstitial Pneumonia. *BMC Bioinformatics* 19:18, 2018. doi:10.1186/s12859-018-2010-z. PMID 29361909 **[V]**. A proximity-to-pathogenic vs proximity-to-neutral score (the PathProx idea). A template for a "proximity to same-disease vs other-disease variants" classifier.
- **Kamburov A, Lawrence MS, Polak P, et al.** Comprehensive assessment of cancer missense mutation clustering in protein structures (CLUMPS). *PNAS* 112(40):E5486–E5495, 2015. doi:10.1073/pnas.1516373112. PMID 26392535 **[V]**.
- **Tokheim C, Bhattacharya R, Niknafs N, et al.** Exome-scale discovery of hotspot mutation regions in human cancer using 3D protein structure (HotMAPS). *Cancer Res* 76(13):3719–3731, 2016. doi:10.1158/0008-5472.CAN-15-3190. PMID 27197156 **[V title]**.
- **Meyer MJ, Lapcevic R, Romero AE, et al.** mutation3D: cancer gene prediction through atomic clustering of coding variants in the structural proteome. *Hum Mutat* 37(5):447–456, 2016. doi:10.1002/humu.22963. PMID 26841357 **[V title]**. Cancer driver clustering, not disease labels.
- **Lelieveld SH, Wiel L, Venselaar H, et al.** Spatial clustering of de novo missense mutations identifies candidate neurodevelopmental disorder-associated genes. *Am J Hum Genet* 101(3):478–484, 2017. doi:10.1016/j.ajhg.2017.08.004. PMID 28867141 **[V title]** (already in project notes).
- **Quinodoz M, et al. 2022 AJHG** (MutScore; gene-specific clustering in 387/840 genes) and **Pérez-Palma E, et al. 2020 Genome Res** (pathogenic variant-enriched regions) are already **[V]** in `paper/notes/literature.md`. **Lal D, May P, Perez-Palma E, et al.** Gene family information facilitates variant interpretation and identification of disease-associated genes in neurodevelopmental disorders. *Genome Med* 12:28, 2020. doi:10.1186/s13073-020-00725-6. PMID 32183904 **[V title]**.
- **Homburger JR, Green EM, Caleshu C, et al.** Multidimensional structure-function relationships in human β-cardiac myosin from population-scale genetic variation. *PNAS* 113(24):6701–6706, 2016. doi:10.1073/pnas.1606950113. PMID 27247418 **[V]**. A **3D spatial scan statistic** on MYH7: HCM variants enriched in the converter and "myosin mesa"; patients with variants in the enriched regions have earlier onset. A single-gene precedent for "3D region → clinical phenotype".
- **Tang ZZ, Sliwoski GR, Chen G, et al.** PSCAN: spatial scan tests guided by protein structures improve complex disease gene discovery and signal variant detection. *Genome Biol* 21:217, 2020. doi:10.1186/s13059-020-02121-0. PMID 32847609 **[V title]**.
- **Marceau West R, Lu W, Rotroff DM, et al.** Identifying individual risk rare variants using protein structure guided local tests (POINT). *PLoS Comput Biol* 15(2):e1006722, 2019. doi:10.1371/journal.pcbi.1006722. PMID 30779729 **[V title]**.
- **Jin B, Capra JA, Benchek P, et al.** An association test of the spatial distribution of rare missense variants within protein structures identifies Alzheimer's disease-related patterns. *Genome Res* 32(4):778–790, 2022. doi:10.1101/gr.276069.121. PMID 35210353 **[V title]**. Case–control 3D tests; adaptable to "disease A vs disease B" labels.
- **Ng JK, Chen Y, Akinwe TM, et al.** Proteome-wide assessment of differential missense variant clustering in neurodevelopmental disorders and cancer. *Cell Genom* 5(4):100807, 2025. doi:10.1016/j.xgen.2025.100807. PMID 40073865. PMCID PMC12008811 **[V, full text grepped]**.
  - *What they did:* 3D-CLUMP on AlphaFold models (948 high-quality models, 342 proteins) found 18 proteins clustered in NDD but not cancer and 19 the reverse.
  - *Note:* "Differential" here means **significant in one class and not the other**. It is not a direct two-sample test of spatial separation.
  - *Gap:* Two coarse phenotype classes (germline NDD vs somatic cancer), not named Mendelian diseases. No position baseline.
- **Nason E, Gerges S, Satterstrom FK, et al. (Daly, Finucane).** Systematic identification of disease-associated 3D neighborhoods in protein structures. *medRxiv* 2026. doi:10.64898/2026.05.29.26354366 **[V, full text read via medRxiv; preprint, not peer reviewed]**.
  - *What they did:* The 3D Neighborhood Test runs Fisher tests in 15 Å spheres, with case/control label permutation, on AlphaFold DB and AF3 models.
  - *Results:* In ClinVar P/LP vs B/LB, 872 of 1,121 proteins had an FDR-significant neighbourhood. **AD genes have more FWER-significant neighbourhoods than AR genes.**
  - **No ClinVar condition labels, and no network proximity.** One anecdote (SETD1A) is a schizophrenia neighbourhood different from the ClinVar one. This is the newest competitor and must be cited.
- **Iqbal S, Pérez-Palma E, Jespersen JB, et al.** 2020 *PNAS* (amino-acid position features) and **Iqbal S, Brünger T, et al.** 2023 *Brain* (Essential3D) are already in project notes **[V title / V]**.
- **Hicks M, Bartha I, di Iulio J, Venter JC, Telenti A.** Functional characterization of 3D protein structures informed by human genetic diversity. *PNAS* 116(18):8960–8965, 2019. doi:10.1073/pnas.1820813116. PMID 30988206 **[V]**. 3D intolerance from population data.

### 1d. Phenotype- or disease-specific structural analyses within genes or gene families

- **Brünger T, Ivaniuk A, Pérez-Palma E, et al.** Conserved missense variant pathogenicity and correlated phenotypes across paralogous genes. *Genome Biol* 26:197, 2025. doi:10.1186/s13059-025-03663-x. PMID 40624578. PMCID PMC12235991 **[V, full text read]**.
  - *Data:* 945 affected individuals, 886 unique missense variants, 22 expert-curated phenotypes in 10 VGSC genes, all mapped onto the SCN2A (Nav1.2) structure.
  - *Method:* For each residue, count patients with a variant at that residue or a Cα-proximal neighbour, per phenotype. Then compute **Kendall correlations between phenotype-specific 3D distributions**, Bonferroni-corrected. Example: SCN1A Dravet vs SCN2A autism, R = 0.31.
  - *Gap for us:*
    - one gene family;
    - **patient counts rather than distinct variants** (pseudoreplication, cf. our first paper);
    - a smoothed-count correlation rather than a pairwise distance test;
    - no position or domain control (alignment position is shared by construction);
    - no network proximity.
- **Brünger T, Pérez-Palma E, Montanucci L, et al.** Conserved patterns across ion channels correlate with variant pathogenicity and clinical phenotypes. *Brain* 146(3):923–934, 2023. doi:10.1093/brain/awac305. PMID 36036558 **[V, full text grepped]**. Pore-axis distance is associated with seizure onset and ID (SCN, GRIN genes), and SCN LoF variants lie near the extracellular pore while GoF variants lie on the intracellular side.
- **Brunklaus A, et al. 2022 *Brain*** and **Heyne HO, et al. 2020 *Sci Transl Med*** (paralog/GoF–LoF prediction for channels) are already in project notes.
- **Wen J, Zeng S, Bonzel C-L, et al. (Cai T).** Phenotypic prediction of missense variants via deep contrastive learning (PheMART). *Nat Biomed Eng*, published online 14 Apr 2026. doi:10.1038/s41551-026-01636-4. PMID 41981312. PMCID PMC13229611 **[V, full text partly read]**.
  - *What they did:* Contrastive learning maps missense variants (protein-language-model embeddings of wild-type vs mutant; **no structure input**) and 4,179 ClinVar phenotypes (UMLS CUIs) into a shared space.
  - Within-gene allelic examples are qualitative only: TUBB3 torsion dystonia vs cortical dysplasia, and CASR. Metrics are computed across all variants, not within genes.
  - *Gap:* No formal within-gene test and no structure. It is a strong **ML competitor for Q1**: a reviewer may ask whether structure beats a PheMART-like sequence model.
- **Di Donato N, et al.** Molecular genotype-phenotype correlation in ACTB- and ACTG1-related non-muscle actinopathies. *Am J Hum Genet* 113(2):324–341, 2026. doi:10.1016/j.ajhg.2025.12.007. PMID 41529692 **[V]**. "At least eight distinct clinical disorders": polymerisation-dysregulating missense → neuronal migration defects, degradation-enhancing variants → milder. A recent single-gene allelic-disorder paper; candidate positive control.
- **Lemmers RJLF, et al.** SMCHD1 mutation spectrum for FSHD2 and BAMS reveals disease-specific localisation of variants in the ATPase domain. *J Med Genet* 56(10):693–700, 2019. doi:10.1136/jmedgenet-2019-106168. PMID 31243061 **[V]**.
- **Piticchio SG, Hosseini N, Grigelioniene G, Orellana L.** A structure-aware framework for genomic variant interpretation in genetic skeletal disorders. *bioRxiv* 2026. doi:10.64898/2026.03.15.711892 **[V, abstract read from PDF; preprint]**. 674 skeletal-disorder genes with AF2, multimers and ClinVar variants. Mostly descriptive: coverage, multimers and an interface "structural equivalence" principle. No within-gene disease-label test.
- **Guzmán-Vega FJ, Cardona-Londoño KJ, González-Álvarez AC, et al. (Arold ST).** VarLand: a pipeline to map the structural landscape of missense variants at the proteome scale. *J Biol Chem* 302(2):111071, 2025. doi:10.1016/j.jbc.2025.111071. PMID 41419200 **[V, full text read on PMC]**. 29 structural features (AlphaFold, DSSP, FoldX, M-CSA) with enrichment vs gnomAD by MONDO disease category. **It does not compare diseases within a gene and computes no inter-variant distances.**
- **Gao M, Zhou H, Skolnick J.** Insights into disease-associated mutations in the human proteome through protein structural analysis. *Structure* 23(7):1362–1369, 2015. doi:10.1016/j.str.2015.03.028. PMID 26027735 **[V]**. Disease vs neutral only.
- **David A, Razali R, Wass MN, Sternberg MJE.** Protein–protein interaction sites are hot spots for disease-associated nonsynonymous SNPs. *Hum Mutat* 33(2):359–363, 2012. doi:10.1002/humu.21656. PMID 22072597 **[V title]**.

### 1e. Structure-annotation tools named in the brief (none relates structure to *which* disease)

| Tool | Citation | What it does | Relation to "which disease" |
|---|---|---|---|
| MISCAST | Iqbal S, Hoksza D, Pérez-Palma E, et al. *Nucleic Acids Res* 48(W1):W132–W139, 2020. doi:10.1093/nar/gkaa361. PMID 32402084 **[V title]** | Per-variant structural annotation | Pathogenic vs population only |
| Missense3D | Ittisoponpisan S, Islam SA, Khanna T, et al. Can predicted protein 3D structures provide reliable insights into whether missense variants are disease associated? *J Mol Biol* 431(11):2197–2212, 2019. doi:10.1016/j.jmb.2019.04.009. PMID 30995449 **[V title]**. Also: Missense3D-DB (Khanna T, et al. *Hum Genet* 140:805–812, 2021. doi:10.1007/s00439-020-02246-z. PMID 33502607); Missense3D-PPI (Pennica C, et al. *J Mol Biol* 435:168060, 2023. doi:10.1016/j.jmb.2023.168060. PMID 37356905); Missense3D-TM (Hanna G, et al. *J Mol Biol* 436:168374, 2024. PMID 38182301); Tsitsa I, et al. *J Mol Biol* 2026:170066. doi:10.1016/j.jmb.2026.170066. PMID 42869221 **[V title, all]** | Structural damage flags | Damaging vs not |
| PhyreRisk | Ofoegbu TC, David A, Kelley LA, et al. *J Mol Biol* 431(13):2460–2466, 2019. doi:10.1016/j.jmb.2019.04.043. PMID 31075275 **[V title]** | Maps variants to 3D | Visual only |
| VarMap | Stephenson JD, Laskowski RA, Nightingale A, Hurles ME, Thornton JM. *Bioinformatics* 35(22):4854–4856, 2019. doi:10.1093/bioinformatics/btz482. PMID 31192369 **[V title]** | Genomic → protein/structure mapping | — |
| VarSite | Laskowski RA, Stephenson JD, Sillitoe I, Orengo CA, Thornton JM. *Protein Sci* 29(1):111–119, 2020. doi:10.1002/pro.3746. PMID 31606900 **[V title]** | Disease variants by structural context, with disease grouping | Descriptive disease-class enrichment |
| ProtVar | Stephenson JD, et al. *Nucleic Acids Res* 52(W1):W140–W147, 2024. doi:10.1093/nar/gkae413. PMID 38769064 **[V title]** | Includes predicted interfaces and pockets | Useful for an interface covariate |
| mutation3D | See §1c | Cancer clustering | — |
| "POSE" | **Not found.** I could not identify a published structural-variant tool called POSE. Ask the requester for the reference; do not cite. | | |

### 1f. Mode-of-inheritance predictors (Q3 baselines)

- **Quinodoz M, Royer-Bertrand B, Cisarova K, et al.** DOMINO: using machine learning to predict genes associated with dominant disorders. *Am J Hum Genet* 101(4):623–629, 2017. doi:10.1016/j.ajhg.2017.09.001. PMID 28985496 **[V title]**. Gene level.
- **Danzi MC, Dohrn MF, Fazal S, et al.** Deep structured learning for variant prioritization in Mendelian diseases (MAVERICK). *Nat Commun* 14:4167, 2023. doi:10.1038/s41467-023-39306-7. PMID 37443090 **[V]**. Variant-level benign / dominant-pathogenic / recessive-pathogenic classifier.
- **Petrazzini BO, Balick DJ, Forrest IS, et al.** Ensemble and consensus approaches to prediction of recessive inheritance for missense variants in human disease (MOI-Pred, ConMOI). *Cell Rep Methods* 4(12):100914, 2024. doi:10.1016/j.crmeth.2024.100914. PMID 39657681 **[V]**.
- **Implication for Q3.** MAVERICK and MOI-Pred already predict AD vs AR at the variant level from sequence and annotation features. Q3 must show that structural context (pLDDT region, RWR, burial) separates AD from AR **within genes** beyond these predictors' scores, or at least beyond AlphaMissense. If MAVERICK scores were trained on ClinVar inheritance, beware circularity (§4c).

### 1g. Network / graph proximity of residues (Q2's RWR component)

The previous notes (§2c of `paper/notes/literature.md`) already cover these methods **[V]**: Chennubhotla & Bahar 2006 (Markov propagation and commute times on residue networks), del Sol 2006, Atilgan 2004, Amor 2016, the Dokholyan "Ohm" method, Ponzoni & Bahar 2018/Rhapsody 2020 (ENM features for pathogenicity), Cowen 2017 (review) and Smedley 2014 (ExomeWalker).

**What I searched for:**
- PubMed: "(residue interaction network OR residue contact network OR protein structure network) AND (random walk OR diffusion OR network propagation OR shortest path) AND (disease mutations OR missense variants) AND phenotype" (18 hits, none relevant); "network proximity AND missense AND residue".
- Europe PMC.
- Two web searches.

I found **no study that uses RWR / personalised PageRank / diffusion similarity on residue contact graphs to compare *pairs of variants* and relate that similarity to whether they cause the same disease.** The closest is the allostery and communication literature (Chennubhotla & Bahar; Verkhivker 2019 BBA, PMID 30339916 **[V title]**, a cancer hotspot communication case study). Absence of evidence after this search is not proof, so do a final bioRxiv check before submission.

**Reviewer risk to pre-empt.** On a Cα 8 Å graph with α = 0.3, the RWR profile similarity of two seeds is dominated by their few-hop overlap. It is therefore approximately a monotone function of graph distance, and hence of 3D distance. Q2 must report RWR proximity **conditional on** both Euclidean Cα distance and shortest-path hop distance, and on sequence separation; otherwise "RWR" will be read as a relabelled distance. A useful variant: compare the *region-mass fingerprints* (Core/Semi-core/Disorder) of the two seeds, which is **not** a distance, alongside profile overlap.

---

## 2. Positive controls (Task 2)

**Selection criteria.**
- Published missense-variant disease separation with a structural or regional basis.
- Single AlphaFold model (UniProt length ≤ 2,700 aa).
- Many ClinVar P/LP missense variants.

OMIM IDs and inheritance are from HPO annotations **[V]**; lengths are from UniProt **[V]**.

"Basis" statements are **[V]** when I read them in an abstract, otherwise **[UNVERIFIED detail]**. I did not count ClinVar P/LP missense variants per condition; do this in the pipeline before freezing the list.

### 2a. Tier 1: strong, structure-explicit, ≤ 2,700 aa (recommended core set)

| # | Gene (UniProt, aa) | Diseases (OMIM, inheritance) | Structural / regional basis | Key citations | Q3? |
|---|---|---|---|---|---|
| 1 | **LMNA** (P02545, 664) | FPLD2 151660 AD; DCM1A 115200 AD; EDMD2 181350 AD; EDMD3 616516 AR; HGPS 176670; MADA 248370 AR; CMT2B1 605588 AR | FPLD2 hotspot R482W/Q in the Ig-fold vs muscle-disease R453W have "different" molecular consequences **[V]** (Krimm 2002). Striated-muscle variants span the rod and Ig-fold. HGPS is the c.1824C>T cryptic splice: **exclude** it. | Shackleton S, et al. *Nat Genet* 24:153–156, 2000. doi:10.1038/72807. PMID 10655060 **[V title; reused]**. Krimm I, et al. *Structure* 10(6):811–823, 2002. doi:10.1016/S0969-2126(02)00777-3. PMID 12057196 **[V]**. Raffaele Di Barletta M, et al. Different mutations in the LMNA gene cause autosomal dominant and autosomal recessive EDMD. *Am J Hum Genet* 66:1407–1412, 2000. doi:10.1086/302869. PMID 10739764 **[V title]**. Novelli G, et al. MADA. *Am J Hum Genet* 71:426–431, 2002. doi:10.1086/341908. PMID 12075506 **[V title]**. | Yes (EDMD2 vs EDMD3; MADA R527H) |
| 2 | **PTPN11** (Q06124, 593) | NS1 163950 AD; LEOPARD/NSML 151100 AD; metachondromatosis 156250 AD; JMML 607785 | NS GoF at the N-SH2/PTP interface; NSML variants are catalytically impaired with DN effects; metachondromatosis is LoF | Tartaglia M, et al. *Am J Hum Genet* 70:1555–1563, 2002. doi:10.1086/340847. PMID 11992261 **[V title; reused]**. Kontaridis MI, et al. PTPN11 (Shp2) mutations in LEOPARD syndrome have dominant negative, not activating, effects. *J Biol Chem* 281:6785–6792, 2006. doi:10.1074/jbc.M513068200. PMID 16377799 **[V title]**. | — |
| 3 | **FGFR3** (P22607, 806) | ACH 100800 (G380R, TM); TD1 187600 (extracellular cysteine creation, e.g. R248C, Y373C); TD2 187601 (K650E, kinase); HCH 146000 (N540K, kinase); SADDAN 616482 (K650M); Muenke 602849 (P250R); CATSHL 610474 (LoF; AD/AR) | Classic residue/domain → disease; GoF strength ranking | Rousseau F, et al. *Nature* 371:252–254, 1994. doi:10.1038/371252a0. PMID 8078586 **[V title]**. Tavormina PL, et al. Thanatophoric dysplasia (types I and II) caused by distinct mutations in FGFR3. *Nat Genet* 9:321–328, 1995. doi:10.1038/ng0395-321. PMID 7773297 **[V title]**. Bellus GA, et al. A recurrent mutation in the tyrosine kinase domain of FGFR3 causes hypochondroplasia. *Nat Genet* 10:357–359, 1995. doi:10.1038/ng0795-357. PMID 7670477 **[V title]**. | CATSHL (AD/AR) |
| 4 | **RET** (P07949, 1114) | MEN2A 171400 (extracellular cysteines, C634 etc.); MEN2B 162300 (M918T, kinase); FMTC 155240; Hirschsprung (LoF; OMIM 142623 is not mapped to RET in HPO, verify) | GoF cysteine-rich domain vs kinase-domain activating vs dispersed LoF | Mulligan LM, et al. *Nature* 363:458–460, 1993. doi:10.1038/363458a0. PMID 8099202 **[V title]**. Hofstra RM, et al. *Nature* 367:375–376, 1994. doi:10.1038/367375a0. PMID 7906866 **[V title]**. Romeo G, et al. Point mutations affecting the tyrosine kinase domain of RET in Hirschsprung's disease. *Nature* 367:377–378, 1994. doi:10.1038/367377a0. PMID 8114938 **[V title]**. Eng C, et al. *JAMA* 276:1575–1579, 1996. PMID 8918855 **[V title]**. | — |
| 5 | **KCNQ1** (P51787, 676) | LQT1 192500 AD; JLNS1 220400 AR; SQT2 609621 AD; familial AF3 607554 AD | Transmembrane vs C-terminal location and DN vs haploinsufficiency predict LQT1 severity **[V]**. GoF (S140G in S1) → AF. Recessive JLN alleles: KCNQ1 is a homotetramer with KCNE1, so DN variants sit at subunit interfaces **[UNVERIFIED detail]**. | Moss AJ, et al. *Circulation* 115:2481–2489, 2007. doi:10.1161/CIRCULATIONAHA.106.665406. PMID 17470695 **[V]**. Neyroud N, et al. *Nat Genet* 15:186–189, 1997. doi:10.1038/ng0297-186. PMID 9020846 **[V title]**. Chouabe C, et al. Properties of KvLQT1 K⁺ channel mutations in Romano-Ward and Jervell and Lange-Nielsen inherited cardiac arrhythmias. *EMBO J* 16:5472–5479, 1997. doi:10.1093/emboj/16.17.5472. PMID 9312006 **[V title]**. Chen YH, et al. KCNQ1 gain-of-function mutation in familial atrial fibrillation. *Science* 299:251–254, 2003. doi:10.1126/science.1077771. PMID 12522251 **[V title]**. Bellocq C, et al. *Circulation* 109:2394–2397, 2004. doi:10.1161/01.CIR.0000130409.72142.FE. PMID 15159330 **[V title]**. | **Yes** (RWS vs JLNS) |
| 6 | **SCN4A** (P35499, 1836) | HyperPP 170500; PMC 168300; PAM 608390; HypoPP2 613345 (all AD GoF; HypoPP = S4 arginine gating-pore current); CMS16 614198 AR; congenital myopathy 22A/22B 620351/620369 AR (LoF) | VSD S4 arginines → gating pore (HypoPP) vs inactivation-gate GoF (PMC/HyperPP) vs biallelic LoF | Sokolov S, et al. Gating pore current in an inherited ion channelopathy. *Nature* 446:76–78, 2007. doi:10.1038/nature05598. PMID 17330043 **[V title]**. Zaharieva IT, et al. Loss-of-function mutations in SCN4A cause severe foetal hypokinesia or 'classical' congenital myopathy. *Brain* 139:674–691, 2016. doi:10.1093/brain/awv352. PMID 26700687 **[V title]**. Tsujino A, et al. *PNAS* 100:7377–7382, 2003. doi:10.1073/pnas.1230273100. PMID 12766226 **[V title]**. | **Yes** (AD GoF vs AR LoF) |
| 7 | **CLCN1** (P35523, 988) | Thomsen 160800 AD; Becker 255700 AR; congenital myopathy 255300 AR | DN (dimeric channel) vs recessive LoF. Many alleles behave dominantly in some families and recessively in others **[UNVERIFIED detail]**, so this is a **mixed / weak** control. | Koch MC, et al. The skeletal muscle chloride channel in dominant and recessive human myotonia. *Science* 257:797–800, 1992. doi:10.1126/science.1379744. PMID 1379744 **[V title]**. Pusch M. Myotonia caused by mutations in the muscle chloride channel gene CLCN1. *Hum Mutat* 19:423–434, 2002. doi:10.1002/humu.10063. PMID 11933197 **[V title]**. Suetterlin K, et al. 223 families with myotonia. *Brain* 145:607–620, 2022. doi:10.1093/brain/awab344. PMID 34529042 **[V title]**. | **Yes** |
| 8 | **GJB2** (P29033, 226) | DFNB1A 220290 AR; DFNA3A 601544 AD; KID 148210 AD; Vohwinkel 124500 AD; PPK+deafness 148350 AD; Bart-Pumphrey 149200 AD; hystrix-like ichthyosis+deafness 602540 AD | Dominant syndromic (skin) missense variants concentrate in the N-terminus/first extracellular loop (e.g. D50N KID, D66H Vohwinkel). Recessive deafness alleles are dispersed LoF **[UNVERIFIED as a quantified clustering claim; per reviews]**. Connexin hexamer/dodecamer: **interfaces matter**. | Kelsell DP, et al. *Nature* 387:80–83, 1997. doi:10.1038/387080a0. PMID 9139825 **[V title]**. Richard G, et al. Missense mutations in GJB2 … cause KID syndrome. *Am J Hum Genet* 70:1341–1348, 2002. doi:10.1086/339986. PMID 11912510 **[V title]**. Maestrini E, et al. D66H causes Vohwinkel's syndrome. *Hum Mol Genet* 8:1237–1243, 1999. doi:10.1093/hmg/8.7.1237. PMID 10369869 **[V title]**. Iossa S, et al. GJB2 gene mutations in syndromic skin diseases with sensorineural hearing loss. *Curr Genomics* 12:475–785, 2011. doi:10.2174/138920211797904098. PMID 22547955 **[V title]**. Srinivas M, et al. Human diseases associated with connexin mutations. *BBA Biomembr* 1860:192–201, 2018. doi:10.1016/j.bbamem.2017.04.024. PMID 28457858 **[V title]**. | **Yes** (the canonical AD/AR case) |
| 9 | **CASR** (P41180, 1078) | FHH1 145980 AD (LoF); NSHPT 239200 AD/AR; ADH1 601198 AD (GoF) | 3D modelling: VFTD missense mutations cluster within 10 Å of Ca²⁺-binding sites. The same residues (173, 221) give activating vs inactivating changes **[V]**. NSHPT is mostly homozygous/compound-heterozygous truncating **[V]**. | Hannan FM, et al. Identification of 70 calcium-sensing receptor mutations … clustering of extracellular domain mutations at calcium-binding sites. *Hum Mol Genet* 21:2768–2778, 2012. doi:10.1093/hmg/dds105. PMID 22422767 **[V]**. | FHH vs NSHPT is **zygosity, not mechanism** (null-type Q3 control). FHH vs ADH is a GoF/LoF positive. Note that same-residue opposite effects make this a *hard* case for 3D distance. |
| 10 | **GCK** (P35557, 465) | MODY2 125851 AD (inactivating); PNDM1 606176 AR (homozygous inactivating); HHF3 602485 AD (activating) | "The majority of activating mutations cluster in … the allosteric activator site"; inactivating mutations are dispersed **[V]** | Osbak KK, et al. *Hum Mutat* 30:1512–1526, 2009. doi:10.1002/humu.21110. PMID 19790256 **[V]**. | HHF3 vs MODY2: **strong positive** (GoF pocket vs dispersed LoF). MODY2 vs PNDM1: **null-type Q3 control** (same LoF class; dosage). |
| 11 | **KCNJ11** (Q14654, 390) and **ABCC8** (Q09428, 1581) | KCNJ11: PNDM2 618856 AD (activating), TNDM3 610582, HHF2 601820 AR (LoF). ABCC8: PNDM3 618857 AD/AR, TNDM2 610374, HHF1 256450 AD/AR | Activating mutations at ATP-binding / Kir6.2–SUR1 interfaces vs LoF trafficking/gating **[UNVERIFIED detail]**. **Octameric KATP: monomer models miss interfaces.** | Gloyn AL, et al. Activating mutations in … Kir6.2 and permanent neonatal diabetes. *N Engl J Med* 350:1838–1849, 2004. doi:10.1056/NEJMoa032922. PMID 15115830 **[V title]**. Pipatpolkai T, Usher S, Stansfeld PJ, Ashcroft FM. New insights into KATP channel gene mutations and neonatal diabetes mellitus. *Nat Rev Endocrinol* 16:378–393, 2020. doi:10.1038/s41574-020-0351-y. PMID 32376986 **[V title]**. | **Yes** |
| 12 | **SCN5A** (Q14524, 2016) | LQT3 603830 AD (GoF); Brugada1 601144 AD (LoF); DCM1E 601154; SSS1 608567 AR; PFHB1A 113900 | LQTS: S3–S5+S6 and DIII/DIV linker are high-probability regions; BrS: transmembrane segments only **[V]** | Kapplinger JD, et al. *Circ Cardiovasc Genet* 8:582–595, 2015. doi:10.1161/CIRCGENETICS.114.000831. PMID 25904541 **[V]**. Brünger 2025 (§1d). | Partly (SSS AR) |
| 13 | **SCN2A** (Q99250, 2005) | DEE11 613721 AD; BFNIS3 607745 AD; EA9 618924 AD | Early-infantile epilepsy = GoF; later-onset epilepsy / ID / ASD = LoF and truncating **[V]** | Wolff M, et al. *Brain* 140:1316–1336, 2017. doi:10.1093/brain/awx054. PMID 28379373 **[V]**. Brünger 2023 (§1d). | — |
| 14 | **SCN1A** (P35498, 2009) | Dravet 607208; GEFS+2 604403; FHM3 609634; DEE6B 619317 (all AD) | FHM3 GoF missense vs Dravet LoF. Dravet vs GEFS+ is partly a severity spectrum (see §3). | Dichgans M, et al. Mutation in … SCN1A in familial hemiplegic migraine. *Lancet* 366:371–377, 2005. doi:10.1016/S0140-6736(05)66786-4. PMID 16054936 **[V title]**. Brünger 2025. | — |
| 15 | **SCN9A** (Q15858, 1988) | Erythermalgia 133020 AD; PEPD 167400 AD; CIP 243000 AR; HSAN2D 201300 AR | GoF (activation shift = IEM; impaired fast inactivation = PEPD) vs biallelic LoF. "Allelic variants underlie distinct channel defects and phenotypes." | Yang Y, et al. *J Med Genet* 41:171–174, 2004. doi:10.1136/jmg.2003.012153. PMID 14985375 **[V title]**. Fertleman CR, et al. *Neuron* 52:767–774, 2006. doi:10.1016/j.neuron.2006.10.006. PMID 17145499 **[V title]**. Cox JJ, et al. *Nature* 444:894–898, 2006. doi:10.1038/nature05413. PMID 17167479 **[V title]**. | **Yes** |
| 16 | **CACNA1A** (O00555, 2506) | FHM1 141500 AD (GoF missense); EA2 108500 AD (LoF, many PTVs); DEE42 617106; SCA6 183086 (CAG expansion: **exclude**) | GoF vs LoF; EA2 is truncating-rich, so restrict to missense | Ophoff RA, et al. *Cell* 87:543–552, 1996. doi:10.1016/S0092-8674(00)81373-2. PMID 8898206 **[V title]**. | — |
| 17 | **TP63** (Q9H3D4, 680) | EEC3 604292; AEC/Hay-Wells 106260; RHS 129400; ADULT 103285; LMS 603543; SHFM4 605289 (all AD) | EEC = DNA-binding domain missense; AEC = SAM-domain missense **[V title]** (van Bokhoven 2001; McGrath 2001). Badonyi 2025: RHS (SAM, LoF-like) vs SHFM4 (surface, non-LoF) **[V]**. | van Bokhoven H, et al. *Am J Hum Genet* 69:481–492, 2001. doi:10.1086/323123. PMID 11462173 **[V title]**. McGrath JA, et al. Hay-Wells syndrome is caused by heterozygous missense mutations in the SAM domain of p63. *Hum Mol Genet* 10:221–229, 2001. doi:10.1093/hmg/10.3.221. PMID 11159940 **[V title]**. | Tetramer (interfaces) |
| 18 | **SMCHD1** (A6NHR9, 2005) | FSHD2 158901 (digenic); BAMS 603457 AD | BAMS missense variants are restricted to the extended ATPase domain, while FSHD2 variants are anywhere: "largely disease-specific localisation" **[V]** | Lemmers RJLF, et al. *J Med Genet* 56:693–700, 2019. PMID 31243061 **[V]**. Gordon CT, et al. *Nat Genet* 49:249–255, 2017. doi:10.1038/ng.3765. PMID 28067911 **[V title]**. Shaw ND, et al. *Nat Genet* 49:238–248, 2017. doi:10.1038/ng.3743. PMID 28067909 **[V title]** (reports overlap: some FSHD2-type variants in BAMS). | — |
| 19 | **TRPV4** (Q9HBA0, 871) | Neuromuscular: CMT2C 606071, SPSMA 181405, dHMN8/CDSMA 600175. Skeletal: metatropic 156530, SMDK 184252, brachyolmia 113500, SEDM 184095, parastremmatic 168400. Arthropathy 606835. (All AD) | Neuropathy variants in the ankyrin-repeat domain ("alterations in the ankyrin domain", Auer-Grumbach 2010 title) vs skeletal variants more widespread **[UNVERIFIED as a quantified claim]**. Both classes are GoF. | Deng HX, et al. *Nat Genet* 42:165–169, 2010. doi:10.1038/ng.509. PMID 20037587 **[V]**. Auer-Grumbach M, et al. *Nat Genet* 42:160–164, 2010. doi:10.1038/ng.508. PMID 20037588 **[V title]**. Landouré G, et al. *Nat Genet* 42:170–174, 2010. doi:10.1038/ng.512. PMID 20037586 **[V title]**. Nilius B, Voets T. The puzzle of TRPV4 channelopathies. *EMBO Rep* 14:152–163, 2013. doi:10.1038/embor.2012.219. PMID 23306656 **[V title]**. | — (tetramer) |
| 20 | **DNM2** (P50570, 870) | CNM1 160150 AD; CMT DIB 606482 AD; LCCS5 615368 AR | CMT mutations in the PH domain (title) vs CNM in the middle and PH domains **[UNVERIFIED detail]** | Züchner S, et al. Mutations in the pleckstrin homology domain of dynamin 2 cause dominant intermediate CMT. *Nat Genet* 37:289–294, 2005. doi:10.1038/ng1514. PMID 15731758 **[V title]**. Bitoun M, et al. *Nat Genet* 37:1207–1209, 2005. doi:10.1038/ng1657. PMID 16227997 **[V title]**. | Small AR |
| 21 | **TUBB3** (Q13509, 450) | CFEOM3A 600638 AD; complex cortical dysplasia CDCBM1 614039 AD | Distinct residue sets perturb kinesin interactions (CFEOM) vs other effects (malformations) **[V title only]**. Also a PheMART example. | Tischfield MA, et al. Human TUBB3 mutations perturb microtubule dynamics, kinesin interactions, and axon guidance. *Cell* 140:74–87, 2010. doi:10.1016/j.cell.2009.12.011. PMID 20074521 **[V title]**. | — (heterodimer) |
| 22 | **STAT1** (P42224, 750) | CMC/IMD31C 614162 AD (GoF); IMD31A 614892 AD (LoF, DN, MSMD); MSMD AR 613796 (LoF) | "All of these [GoF] mutations affect the coiled-coil domain"; LoF alleles elsewhere **[V]** | Liu L, et al. *J Exp Med* 208:1635–1648, 2011. doi:10.1084/jem.20110958. PMID 21727188 **[V]**. | **Yes** |
| 23 | **STAT3** (P40763, 770) | HIES 147060 AD (DN, DNA-binding/SH2); infantile multisystem autoimmunity 615952 AD (GoF) | DN LoF vs GoF | Milner JD, et al. *Blood* 125:591–599, 2015. doi:10.1182/blood-2014-09-602763. PMID 25359994 **[V title]**. | — |
| 24 | **TGFBI** (Q15582, 683) | Lattice I 122200 (R124C); Avellino 607541 (R124H); Groenouw I 121900 (R555W); Reis-Bücklers 608470 (R124L); Thiel-Behnke 602082 (R555Q); lattice IIIA 608471 (all AD) | Residue-specific: different substitutions at R124/R555 → different dystrophies | Munier FL, et al. Kerato-epithelin mutations in four 5q31-linked corneal dystrophies. *Nat Genet* 15:247–251, 1997. doi:10.1038/ng0397-247. PMID 9054935 **[V title]**. | **A key "same residue, different disease" stress test:** 3D distance is 0, so only substitution-aware features (AlphaMissense) can separate. |
| 25 | **RHO** (P08100, 348) | RP4 613731 AD/AR; CSNB1 610445 AD (G90D, A292E, T94I constitutive activation); fundus albipunctatus 136880 | Class-based RP mechanisms (misfolding vs other) vs constitutively active CSNB alleles near the Schiff-base region **[V title / UNVERIFIED detail]** | Rao VR, et al. Rhodopsin mutation G90D and a molecular mechanism for congenital night blindness. *Nature* 367:639–642, 1994. doi:10.1038/367639a0. PMID 8107847 **[V title]**. Athanasiou D, et al. *Prog Retin Eye Res* 62:1–23, 2018. doi:10.1016/j.preteyeres.2017.10.002. PMID 29042326 **[V title]**. | Yes (adRP vs arRP) |
| 26 | **PRPH2** (P23942, 346) | RP7 608133 AD/AR; patterned MD 169150; vitelliform MD3 608161; CACD2 613105; fundus albipunctatus 136880 | D2-loop / tetramer variants; RP vs macular; reported distinct spatial patterns (Sun 2026, preprint) | Boon CJ, et al. The spectrum of retinal dystrophies caused by mutations in the peripherin/RDS gene. *Prog Retin Eye Res* 27:213–235, 2008. doi:10.1016/j.preteyeres.2008.01.002. PMID 18328765 **[V title]**. Heath Jeffery RC, et al. *Invest Ophthalmol Vis Sci* 65(5):22, 2024. doi:10.1167/iovs.65.5.22. PMID 38743414 **[V title]**. Sun YJ, et al. medRxiv 2026 (in project notes) **[V abstract]**. | Partly. Known for intrafamilial variability, so it is **mixed** (see §3). |
| 27 | **PIEZO1** (Q92508, 2521; under the 2,700 limit) | DHS1 194380 AD (GoF); lymphatic dysplasia/LMPH3 616843 AR (LoF) | GoF clustered in anchor, inner helix and C-terminal domains **[V]** (Rosato 2026) vs biallelic LoF | Zarychanski R, et al. *Blood* 120:1908–1915, 2012. doi:10.1182/blood-2012-04-422253. PMID 22529292 **[V title]**. Fotiou E, et al. *Nat Commun* 6:8085, 2015. doi:10.1038/ncomms9085. PMID 26333996 **[V title]**. Rosato BE, et al. *Blood* 148:229–244, 2026. doi:10.1182/blood.2025031632. PMID 42013021 **[V]**. | **Yes** (trimer) |
| 28 | **ATP1A3** (P13637, 1013) | AHC2 614820; DYT12/RDP 128235; CAPOS 601338 (E818K only); DEE99 619606 (all AD) | Distinct residue sets per syndrome; CAPOS is a single recurrent variant | Heinzen EL, et al. *Nat Genet* 44:1030–1034, 2012. doi:10.1038/ng.2358. PMID 22842232 **[V title]**. de Carvalho Aguiar P, et al. *Neuron* 43:169–175, 2004. doi:10.1016/j.neuron.2004.06.028. PMID 15260953 **[V title]**. Demos MK, et al. *Orphanet J Rare Dis* 9:15, 2014. doi:10.1186/1750-1172-9-15. PMID 24468074 **[V title]**. | — |
| 29 | **AARS1** (P49588, 968) | CMT2N 613287 AD; DEE29 616339 AR; HDLS2 619661 AD; TTD8 619691 AR | DEE29: ATP-binding / acceptor recognition, LoF. CMT2N: anticodon-binding domain and putative dimer interface, DN **[V via Badonyi 2025]** | Latour P, et al. *Am J Hum Genet* 86:77–82, 2010. doi:10.1016/j.ajhg.2009.12.005. PMID 20045102 **[V title]**. Simons C, et al. *Am J Hum Genet* 96:675–681, 2015. doi:10.1016/j.ajhg.2015.02.012. PMID 25817015 **[V title]**. | **Yes** |
| 30 | **CLCN7** (P51798, 805) | OPTA2 166600 AD (DN); OPTB4 611490 AR; HOD 618541 AD | Dominant variants more clustered (EDC 1.36 vs 1.03) **[V via Badonyi 2025]** | Cleiren E, et al. *Hum Mol Genet* 10:2861–2867, 2001. doi:10.1093/hmg/10.25.2861. PMID 11741829 **[V title]**. | **Yes** (phenocopy: same disease name family, so a good DN test) |
| 31 | **FLNA** (P21333, 2647; X-linked) | PVNH 300049 (LoF); OPD1 311300, OPD2 304120, FMD 305620, MNS 309350 (GoF; localised missense) | OPD-spectrum mutations "all clustered into four regions: the actin-binding domain and rod repeats 3, 10 and 14/15", GoF, vs LoF in PVNH **[V]** | Robertson SP, et al. *Nat Genet* 33:487–491, 2003. doi:10.1038/ng1119. PMID 12612583 **[V]**. | X-linked; exclude from Q3 |
| 32 | **MYH7** (P12883, 1935) | HCM1 192600; DCM1S 613426; distal myopathy 1 (Laing) 160500; MSM AD 608358 / AR 255160 | Cardiomyopathy variants in head/converter/mesa vs skeletal myopathies in the rod **[UNVERIFIED for the rod claim]**. HCM converter/mesa enrichment **[V]**. | Homburger 2016 (§1c) **[V]**. | Yes (MSM AR) |

**Additional candidates from Badonyi 2025 Table 4 [V]:** KRAS (CFC2 vs somatic myeloma: drop somatic labels), BRAF, MTOR (Smith-Kingsmore 616638 vs FCD), CAV3 (LQT9 vs rippling muscle disease), ACTA1 (AR vs AD congenital myopathy), ADCY5, LHCGR (hypogonadism vs male precocious puberty: a clean LoF/GoF pair), NPR2, TWNK, POLR3B, SCN4A. **LHCGR and CAV3 are attractive small additions.**

### 2b. Classic allelic genes we cannot model with a single AlphaFold model (> 2,700 aa)

These need fragment handling or should be excluded; say so explicitly.

- **RYR1** (5,038 aa): MHS1 145600 / CCD 117000 (AD/AR) / multiminicore 255320 (AR) / King-Denborough 619542. Robinson R, et al. Mutations in RYR1 in malignant hyperthermia and central core disease. *Hum Mutat* 27:977–989, 2006. doi:10.1002/humu.20356. PMID 16917943 **[V title]**; Hernández Trapero 2025 **[V]**. Note: MHS1 is "POLYGENIC" in HPO `genes_to_disease`.
- **COL7A1** (2,944 aa): DDEB 131750 (glycine substitutions, DN) vs RDEB 226600. Varki R, et al. *J Med Genet* 44:181–192, 2007. doi:10.1136/jmg.2006.045302. PMID 16971478 **[V title]**.
- **FBN1, ITPR1, NF1, KMT2D**: see project notes (§3).

### 2c. Recommended frozen sets

- **Q1/Q2 core:** LMNA, PTPN11, FGFR3, RET, SCN4A, KCNQ1, GJB2, CASR (FHH vs ADH), GCK (HHF3 vs MODY2), SCN5A, SCN2A, SCN9A, TP63, SMCHD1, TRPV4, STAT1, FLNA, PIEZO1, AARS1, ATP1A3.
- **Stress tests:** TGFBI (same residue, different disease); CASR residues 173/221 (same residue, opposite mechanism).
- **Q3 positives (AD vs AR with different mechanism):** GJB2, KCNQ1, SCN4A, SCN9A, PIEZO1, STAT1, AARS1, CLCN7, KCNJ11/ABCC8, LMNA (EDMD2 vs EDMD3/MADA), RHO.
- **Q3 nulls (AD vs AR differ by zygosity/dosage, same LoF class):** GCK MODY2 vs PNDM1; CASR FHH1 vs NSHPT.

---

## 3. Negative / null controls (Task 3)

Each entry: gene, disease labels, and why the labels are a clinical spectrum, modifier effect or dosage effect rather than a variant-location effect.

| Gene | Labels | Why the labels are not location-driven | Citation |
|---|---|---|---|
| **LMNA (within striated muscle)** | EDMD2 181350 vs DCM1A 115200 vs LGMD1B / L-CMD 613205 | In 53 AD-EDMD patients, "no clear correlation between the phenotype and type or localization of the mutations", with "marked inter- and intra-familial variability", ranging from full EDMD to cardiac-only. **Use LMNA as both a positive (FPLD2 vs striated) and a null (EDMD2 vs DCM1A) control within one protein.** | Bonne G, et al. *Ann Neurol* 48(2):170–180, 2000. PMID 10939567 **[V]** |
| **VCP** | IBMPFD1 167320 vs FTD/ALS6 613954 vs CMT2Y 616687 | 231 individuals, 36 families, 15 mutations: no correlation between mutation and the incidence of myopathy, PDB or FTD, apart from one R159C exception. "Large interfamilial and intrafamilial variation." | Al-Obeidi E, et al. *Clin Genet* 93:119–125, 2018. doi:10.1111/cge.13095. PMID 28692196 **[V]** |
| **FGFR2 (Crouzon vs Pfeiffer)** | Crouzon 123500 vs Pfeiffer 101600 (vs Jackson-Weiss 123150) | "Identical mutations in unrelated individuals giving different phenotypes." **Apert 101200 (S252W/P253R linker) remains a positive**, so FGFR2 also gives within-gene positive and null contrasts. | Rutland P, et al. *Nat Genet* 9:173–176, 1995. doi:10.1038/ng0295-173. PMID 7719345 **[V]**; Wilkie AO, et al. Apert syndrome results from localized mutations of FGFR2 and is allelic with Crouzon syndrome. *Nat Genet* 9:165–172, 1995. doi:10.1038/ng0295-165. PMID 7719344 **[V title]** |
| **PRNP** | FFI 600072 vs CJD 123400 (D178N) | The disease is "determined by a DNA polymorphism" (codon 129 in cis), not by variant position. GSS 137440 is partly position-driven. | Goldfarb LG, et al. *Science* 258:806–808, 1992. doi:10.1126/science.1439789. PMID 1439789 **[V title]** |
| **GCK** | MODY2 125851 vs PNDM1 606176 | Same inactivating alleles in heterozygous vs homozygous state **[V]** | Osbak 2009 **[V]** |
| **CASR** | FHH1 145980 vs NSHPT 239200 | NSHPT is largely homozygous or compound-heterozygous for FHH-type LoF alleles **[V]** | Hannan 2012 **[V]** |
| **ABCA4** (AR only) | STGD1 248200 vs CRD3 604116 vs RP19 601718 | The phenotype is set by the *combination* of two alleles (residual activity), so a single-variant label is ill-defined **[V title; standard model, detail UNVERIFIED]** | Cremers FPM, et al. *Prog Retin Eye Res* 79:100861, 2020. doi:10.1016/j.preteyeres.2020.100861. PMID 32278709 **[V title]** |
| **SCN1A (Dravet vs GEFS+2)** | 607208 vs 604403 | Severity spectrum with overlapping missense alleles. **FHM3 is a positive.** **[UNVERIFIED as a cited claim; find a review before use]** | — |
| **GNAS** | PHP1a vs PPHP | Parent-of-origin imprinting (from project notes) **[UNVERIFIED]** | — |
| **PRPH2** | RP7 vs macular dystrophies | Marked intrafamilial phenotypic variability for identical variants is reported in reviews **[UNVERIFIED detail]**, yet a preprint reports spatial separation. Treat as **exploratory**, not as a positive or null. | Boon 2008 **[V title]** |
| **Single-disease null cohorts (HPO-level)** | STXBP1, POGZ, EHMT1, TRAF7, U2AF2, ANKRD11 | See project notes. Useful for false-positive calibration, though not multi-disease. | `paper/notes/literature.md` §3 |

---

## 4. Pitfalls and recommended statistics (Task 4)

### 4a. Disease-label quality (ClinVar)

- **ClinVar structure.**
  - Submitter-level condition names are mapped to MedGen concepts.
  - Many submissions carry "not provided", "not specified" or generic labels (e.g. "Cardiomyopathy", "Inborn genetic diseases", "<Gene>-related disorder").
  - The same disease appears under submitter-specific names.

  **Recommendations.**
  - Use the SCV level (submission) rather than the aggregate RCV/VCV condition list.
  - Drop generic and "not provided" conditions.
  - Map MedGen → OMIM/MONDO.
  - Require ≥ 1 star.
  - Keep variants with **exactly one specific condition after collapsing synonyms**.

  Cite: Landrum MJ, et al. ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Res* 46(D1):D1062–D1067, 2018. doi:10.1093/nar/gkx1153. PMID 29165669 **[V title]**. Landrum MJ, et al. ClinVar: improvements to accessing data. *Nucleic Acids Res* 48(D1):D835–D844, 2020. doi:10.1093/nar/gkz972. PMID 31777943 **[V title]**. Landrum MJ, et al. ClinVar: updates to support classifications of both germline and somatic variants. *Nucleic Acids Res* 53(D1):D1313–D1321, 2025. doi:10.1093/nar/gkae1090. PMID 39578691 **[V title]**. Landrum MJ, et al. *Nucleic Acids Res* 44(D1):D862–D868, 2016. doi:10.1093/nar/gkv1222. PMID 26582918 **[V title]**.
- **Classification discordance.** Yang S, Lincoln SE, Kobayashi Y, et al. Sources of discordance among germ-line variant classifications in ClinVar. *Genet Med* 19:1118–1126, 2017. doi:10.1038/gim.2017.60. PMID 28569743 **[V title]**.
- **Lumping vs splitting (label hierarchy).** OMIM splits phenotypic series (e.g. DCM 1A vs EDMD2), while ClinGen and MONDO often lump ("LMNA-related disorder"). The analysis is sensitive to this choice: pre-specify the label level and run a sensitivity analysis at a lumped level. Cite: Thaxton C, Goldstein J, DiStefano M, et al. Lumping versus splitting: how to approach defining a disease to enable accurate genomic curation. *Cell Genom* 2(5):100131, 2022. doi:10.1016/j.xgen.2022.100131. PMID 35754516 **[V title]**.
- **Ontology harmonisation.**
  - Vasilevsky NA, Matentzoglu NA, Toro S, et al. Mondo: unifying diseases for the world, by the world. *medRxiv* 2022. doi:10.1101/2022.04.13.22273750 **[V via Europe PMC preprint record]**.
  - Amberger JS, Bocchini CA, Scott AF, Hamosh A. OMIM.org: leveraging knowledge across phenotype-gene relationships. *Nucleic Acids Res* 47(D1):D1038–D1043, 2019. doi:10.1093/nar/gky1151. PMID 30445645 **[V title]**.
  - Use MONDO's `exactMatch` / `equivalentTo` to OMIM, and test robustness with OMIM-only labels (as in Badonyi 2025).
  - Use HPO semantic similarity between disease pairs to grade "how different" two labels are. Badonyi used Lin similarity on HPO (2024-02-08) via the R packages ontologyIndex / ontologySimilarity **[V]**. This lets near-synonymous pairs (EDMD2/DCM1A) be down-weighted.
- **Phenopacket-store labels** are curator-assigned per individual and may reflect the paper's title disease. Combine them with ClinVar only after deduplicating variants. Treat agreement between sources as a quality filter.
- **Label leakage from position.** Some OMIM entries are *defined* by a variant or region (CAPOS = ATP1A3 E818K; HGPS = c.1824C>T; Apert = S252W/P253R). Such labels trivially "separate" by position. Flag them, and keep them in a separate stratum or remove them from the Q1 headline test.

### 4b. Ascertainment, recurrence and founder variants

- **Unit of analysis.** Use distinct variants, never individuals (see our first paper). Hotspot genes (FGFR3 G380R, RET C634/M918T, GJB2 c.35delG, a PTV) contribute one data point per residue. Consider collapsing to **distinct residues** for distance analyses, since several substitutions at one residue give zero-distance "pairs".
- **Founder and population-specific variants.** These inflate counts in single labels and give apparent clusters. Manrai AK, Funke BH, Rehm HL, et al. Genetic misdiagnoses and the potential for health disparities. *N Engl J Med* 375:655–665, 2016. doi:10.1056/NEJMsa1507092. PMID 27532831 **[V title]** (misclassification arising from ancestry-specific variants).
- **Research-focus bias.** Multi-phenotype genes are over-represented among well-studied genes (Badonyi 2025 notes this **[V]**). Domain-targeted sequencing (e.g. RET hotspot exons, FGFR3 hotspot assays) means some labels are only ever tested in certain exons. Mitigation: a per-gene permutation that preserves each label's variant count and a sequence-window structure.
- **Circularity.**
  - AlphaMissense and other baselines may have been tuned on ClinVar.
  - ClinVar P/LP calls may themselves rest on "located in a known hotspot / domain" (ACMG PM1) or "same residue as known pathogenic" (PS1/PM5), creating artificial clustering of pathogenic variants by label.
  - Cite: Grimm DG, Azencott CA, Aicheler F, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat* 36:513–523, 2015. doi:10.1002/humu.22768. PMID 25684150 **[V title]**.
  - Sensitivity check: restrict to variants with ≥ 2 independent submitters, or to variants also present in phenopacket-store with published phenotypes.

### 4c. AlphaFold caveats

- **Monomer models miss interfaces.** This matters most for DN and GoF alleles, which are enriched at interfaces (Gerasimavicius 2022 **[V]**). Many positive-control genes are obligate oligomers: connexin hexamers (GJB2), K⁺ channel tetramers (KCNQ1, TRPV4), KATP octamers (KCNJ11+ABCC8), TP63/STAT dimers or tetramers, PIEZO1 trimers, CLCN dimers.

  **Recommendation.** Add an *interface* covariate from experimental complexes (PDB), predicted complexes and homo-oligomer models, and rerun the RWR on the **assembly graph** for a subset of genes. Resources:
  - Evans R, O'Neill M, Pritzel A, et al. Protein complex prediction with AlphaFold-Multimer. *bioRxiv* 2021. doi:10.1101/2021.10.04.463034 **[V Crossref]**.
  - Abramson J, Adler J, Dunger J, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493–500, 2024. doi:10.1038/s41586-024-07487-w. PMID 38718835 **[V title]**.
  - Burke DF, Bryant P, Barrio-Hernandez I, et al. Towards a structurally resolved human protein interaction network. *Nat Struct Mol Biol* 30:216–225, 2023. doi:10.1038/s41594-022-00910-8. PMID 36690744 **[V title]**.
  - Humphreys IR, et al. Computed structures of core eukaryotic protein complexes. *Science* 374:eabm4805, 2021. doi:10.1126/science.abm4805. PMID 34762488 **[V title]**.
  - Zhang J, et al. Predicting protein-protein interactions in the human proteome. *Science* 390:eadt1630, 2025. doi:10.1126/science.adt1630. PMID 40997207 **[V title]**.
  - Schweke H, et al. An atlas of protein homo-oligomerization across domains of life. *Cell* 187:999–1010, 2024. doi:10.1016/j.cell.2024.01.022. PMID 38325366 **[V title]**.
  - Interactome INSIDER (Meyer 2018) and ProtVar (Stephenson 2024) for predicted interfaces.
- **Disordered regions and motifs.** pLDDT-defined "Disorder" is not inert: motif–domain interfaces in IDRs carry many pathogenic variants (Hubrich 2026 *NSMB*, **[V]** in project notes; Ruff & Pappu 2021). Pre-specify that pLDDT-region features test *monomer* context only.
- **Single static conformation.** Allosteric or state-dependent GoF (channels, GCK, PTPN11 open/closed) may not be captured. Homburger 2016 used pre- and post-power-stroke models **[V]**.
- **Length limit.** Proteins > 2,700 aa are fragments only (AlphaFold DB; project notes). This excludes RYR1, COL7A1, FBN1, ITPR1, NF1, KMT2D. Report this exclusion as a limitation. Varadi M, et al. AlphaFold Protein Structure Database in 2024. *Nucleic Acids Res* 52(D1):D368–D375, 2024. doi:10.1093/nar/gkad1011. PMID 37933859 **[V title]**.
- **Isoform and numbering.** Map ClinVar HGVS to the UniProt canonical sequence used by AlphaFold, as Badonyi did with VEP **[V]**. MANE and UniProt can differ (NF1, ITPR1 examples in project notes).

### 4d. Recommended statistics

**General principles.**
- **The unit is the distinct variant**, or the distinct residue for distance-based tests.
- **Variant-level permutation within gene**, stratified to preserve label sizes. Phipson B, Smyth GK 2010 (in the manuscript) for permutation p-values.
- Combine genes with a **gene-stratified** statistic: sum of per-gene standardised statistics, or Fisher/Stouffer combination with a permutation reference. Report per-gene results alongside.

**Q1 (structural context differs by disease, beyond position/domain/AM).**
- Within gene, fit a multinomial or pairwise logistic model of disease label on structural features after residualising each feature on spline(position), domain and AlphaMissense. The first paper used carrier-weighted least squares on position; here use variant-weighted.
- Alternatively use **PERMANOVA** on the feature vector, with permutations restricted within sequence blocks. Anderson MJ. A new method for non-parametric multivariate analysis of variance. *Austral Ecol* 26:32–46, 2001. doi:10.1111/j.1442-9993.2001.01070.pp.x **[V Crossref]**.
- Compare cross-validated AUC of nested models (position+domain+AM vs + structure) with **gene-grouped CV** to avoid leakage.

**Q2 (pairwise proximity).**
1. **Mantel / partial Mantel**: same-label indicator matrix vs 3D distance (or 1 − RWR similarity), partialling out sequence distance.
   - Mantel N. The detection of disease clustering and a generalized regression approach. *Cancer Res* 27:209–220, 1967. PMID 6018555 **[V title]**.
   - Smouse PE, Long JC, Sokal RR. Multiple regression and correlation extensions of the Mantel test of matrix correspondence. *Syst Zool* 35(4):627–632, 1986. doi:10.2307/2413122 **[V Crossref]**.

   **Caveat.** Partial Mantel tests are anti-conservative when both matrices are spatially autocorrelated, which they are here: sequence and 3D distance are strongly correlated.
   - Guillot G, Rousset F. Dismantling the Mantel tests. *Methods Ecol Evol* 4(4):336–344, 2013. doi:10.1111/2041-210x.12018 **[V Crossref]**.
   - Legendre P, Fortin M-J, Borcard D. Should the Mantel test be used in spatial analysis? *Methods Ecol Evol* 6:1239–1247, 2015. doi:10.1111/2041-210X.12425 **[V Crossref]**.
   - Legendre P, Fortin M-J. Comparison of the Mantel test and alternative approaches for detecting complex multivariate relationships in the spatial analysis of genetic data. *Mol Ecol Resour* 10:831–844, 2010. doi:10.1111/j.1755-0998.2010.02866.x. PMID 21565094 **[V title]**.

   **Preferred:** a **label permutation that preserves sequence structure**. Options:
   - permute disease labels among variants within sliding sequence windows or within domains;
   - circular-shift labels along the sequence;
   - Moran spectral randomisation.

   The null then retains 1D clustering, and the test asks whether 3D/RWR adds beyond it.
2. **Two-label spatial point-pattern tests.**
   - The case–control difference of K-functions, D(t) = K₁₁(t) − K₂₂(t), or cross-K between labels. Diggle PJ, Chetwynd AG. Second-order analysis of spatial clustering for inhomogeneous populations. *Biometrics* 47:1155–1163, 1991. PMID 1742435 **[V title]**. Diggle PJ, Gómez-Rubio V, Brown PE, Chetwynd AG, Gooding S. Second-order analysis of inhomogeneous spatial point processes using case-control data. *Biometrics* 63:550–557, 2007. doi:10.1111/j.1541-0420.2006.00683.x. PMID 17688507 **[V title]**.
   - Sivley 2018's univariate Ripley K with position permutation is the protein-structure precedent **[V]**.
   - The 3DNT sphere test with label permutation (Nason 2026) and PSCAN / POINT / Jin 2022 are case–control analogues that can be repurposed with disease A vs disease B as "case/control".
3. **Pair-level regression with a nested baseline**: logit P(same label) ~ s(sequence distance) + s(3D distance) + s(RWR similarity) + gene effects. Inference by the structured permutation above, because pairs are not independent.
4. **Leave-one-variant-out nearest-label classifier**, following PathProx (Sivley 2018 RTEL1): predict a held-out variant's disease from proximity-weighted votes. Compare sequence-, 3D- and RWR-based kernels by gene-grouped AUC. This is interpretable and robust.

**Q3 (AD vs AR within mixed genes).**
- Same structure as Q1 with a binary label, plus AlphaMissense and (if allowed) MAVERICK / MOI-Pred / LoGoFunc as baselines.
- Include an mLOF-style group statistic (ΔΔG + EDC, Badonyi 2025) as a comparator.
- Calibrate on the GCK and CASR null pairs: no difference expected, because AD and AR alleles are mechanistically equivalent.

**Multiple testing and power.**
- Pre-register the gene set and label level, and apply BH across genes.
- Report the number of distinct variants per label (e.g. require ≥ 5 per label and ≥ 2 labels).
- Run a power simulation with the first paper's design-effect logic.

---

## 5. Venues (Task 5)

Word-limit details are reused from `paper/notes/literature.md` where verified there **[V]**; otherwise **[UNVERIFIED]**.

| Outcome | Venue | Why |
|---|---|---|
| Positive or clearly incremental | **Am J Hum Genet** | Published Wang-adjacent work (Guo 2013), Sivley 2018, Turner-type analyses and GPSEA. High bar on novelty beyond Badonyi 2025 **[UNVERIFIED limits]**. |
| Positive or clearly incremental | **Genome Biology** / **Genome Medicine** | Published Brünger 2025 and LoGoFunc **[UNVERIFIED limits]**. |
| Positive or clearly incremental | **Nature Communications** | Published Gerasimavicius 2022 and Badonyi 2025. A positive multi-gene result with a new network-proximity angle would fit, but it is a stretch if the result is mixed. |
| Positive or clearly incremental | **Cell Genomics** | Published Ng 2025 3D-CLUMP and Thaxton 2022 **[UNVERIFIED limits]**. |
| Mixed / methods-forward | **Human Mutation** (Wiley/Hindawi, OA) | Published Hernández Trapero 2025, Robinson 2006 and Osbak 2009. Natural home for variant-interpretation methods with controls **[UNVERIFIED limits]**. |
| Mixed / methods-forward | **HGG Advances** | Phenopacket-store venue **[limits UNVERIFIED; author page 403 previously]**. |
| Mixed / methods-forward | **PLOS Computational Biology** | No word limit for research articles; abstract ≤ 300 words **[V, project notes]**. |
| Mixed / methods-forward | **Bioinformatics** | Original Paper ~5,000 words / 7 pages **[V, project notes]**. |
| Mixed / methods-forward | **NAR Genomics and Bioinformatics** | No stated limit for standard papers **[V, project notes]**. |
| Mixed / methods-forward | **Protein Science** or **Structure** | If the emphasis shifts to structural mechanism (interfaces, oligomers) **[UNVERIFIED limits]**. |
| Mixed / methods-forward | **European Journal of Human Genetics** | Clinical-genetics audience **[UNVERIFIED]**. |
| Null | **PLOS ONE**, **BMC Research Notes** (accepts null results; ≤ 2,000-word research note **[V, project notes]**), or a Bioinformatics Application Note framed around the tool | |

---

## 6. Search log (for reproducibility)

- **PubMed (E-utilities) queries included:**
  - allelic disorders / same gene different diseases × protein structure × clustering (1,107 hits, screened top-ranked);
  - edgetic / 3D interactome titles;
  - residue network × random walk / diffusion × variant × phenotype (18 and 2 hits; none relevant);
  - mode-of-inheritance prediction;
  - gene-specific author/title queries for every control.
- **Europe PMC:** "same gene" AND "distinct phenotypes" AND "protein structure" AND missense (48 hits; screened); "phenotype-specific" AND "spatial clustering" (13 hits). Full-text XML read for Brünger 2023/2025, Badonyi 2025, Gerasimavicius 2022 and Ng 2025.
- **Web search (5 queries)** surfaced Nason 2026 (medRxiv), VarLand 2025, PheMART 2026 and Piticchio 2026 (bioRxiv); each was then read.
- **Not found:** a tool named "POSE"; any multi-gene pairwise same- vs different-disease distance test; any residue-level RWR variant-proximity study.
