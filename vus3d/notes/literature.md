# Literature notes: 3D label propagation for missense VUS classification (vus3d)

Compiled 2026-11-03 for the proposed tool. For a missense VUS, it propagates known ClinVar P/LP and B/LB labels from *other* residues of the same protein over the AlphaFold Cα contact graph (8 Å; RWR seeded at the query residue, or label propagation). The result is a gene-specific 3D "hotspot" score, which is then combined with AlphaMissense and calibrated as ACMG/AMP PM1 evidence. Evaluation uses a temporal holdout (train on ClinVar Jan 2023; test on missense first classified P/LP or B/LB after that, up to ClinVar Nov 2026). Comparators are AlphaMissense, 1D label density (MutScore-like) and 3D-radius label density. Same-residue labels are excluded because they are PS1/PM5 evidence.

**How sources were checked.** Bibliographic data come from PubMed E-utilities (esearch/esummary/efetch), PMC full text (efetch), Europe PMC REST, Crossref, and the bioRxiv/medRxiv and arXiv APIs.

**Marking.**
- **[V]**: record verified, and the claim was checked in the abstract or full text.
- **[V title]**: record verified, but the specific claim was not checked.
- **[V preprint]**: record verified on bioRxiv/medRxiv; not peer reviewed.
- **[UNVERIFIED]**: not checked against a primary source.

Citations already verified in `paper/notes/literature.md` and `allelic/notes/literature.md` are reused here and marked "(project notes)".

---

## 0. Bottom line

**Novelty verdict.** The general idea is **largely done**. The specific package is not.

*What has been done:*
- **Using other variants' ClinVar labels in the 3D neighbourhood as a predictor** has been done and published in a high-profile journal. DeMAG (Luppino 2023, *Nat Commun*) built a "partners score": the sum of +1/−1 labels of AlphaFold2 residues within 11 Å Cα plus co-evolving partners. This score was the top feature (median importance ≈ 60%). DeMAG was gene-specific in 59 ACMG SF genes and was tested on ClinVar submissions after Dec 2017.
- **3D proximity to P vs B labels** was done earlier: PathProx (Sivley 2018 RTEL1; now in VUStruct 2026).
- **Gene-specific 3D structural "penetrance density"** comes from Kroncke 2020 and O'Neill 2023 (SCN5A, KCNQ1, KCNH2) and Dauda 2026 (RYR2).
- **1D gene-specific positional label scores** with a temporal test exist: MutScore (Quinodoz 2022), which trained on ClinVar Nov 2020 and tested on entries added through Sep 2021.
- **Data-driven hotspots mapped to PM1 or other ACMG strengths** exist: Walsh 2019, Waring 2021, CardioClassifier 2018, and HCSeeker 2025. HCSeeker is a KDE hotspot/coldspot database for PM1 and benign criteria in 889 genes.
- **Genome-wide 3D neighbourhood tests on ClinVar P/LP vs B/LB**: Nason/Finucane 2026 (3DNT) found 872/1,121 proteins with a significant neighbourhood.

*What I did not find* (in PubMed/Europe PMC/medRxiv/bioRxiv searches up to 2026-11-03):
1. Graph diffusion (RWR or label propagation) of P **and** B labels on AlphaFold contact graphs, tested against a matched 3D-radius and a matched 1D-window label density.
2. A test of **increment over AlphaMissense** under a **post-AlphaMissense-era temporal holdout**.
3. **Calibration of the incremental evidence to PM1 strength** with the ClinGen local-posterior method, respecting the SVI rule that PP3 + PM1 sum to at most Strong.
4. Reporting **within-gene** (gene-stratified) rather than pooled metrics.

*"Already done" risk: HIGH* if the paper is framed as "3D label propagation improves VUS classification". Reviewers will cite DeMAG, PathProx, MutScore, Kroncke and 3DNT. *Differentiation must rest on:*
- (i) diffusion vs radius vs 1D, compared head to head;
- (ii) increment over AlphaMissense, conditional on AlphaMissense;
- (iii) a temporal holdout with leakage stratification (PM1-coded and AM-era labels);
- (iv) calibrated PM1 points with coverage and false-positive rates per gene;
- (v) a pre-registered, honest null if the 3D increment over 1D is ≈ 0. That outcome is plausible, given the PhenoBridge findings that AlphaFold context added nothing beyond sequence position.

**Expected effect size.** Prior label-density and clustering features gave **+0.01 to +0.015 pooled AUROC over REVEL** (MutScore) and **no AUROC gain** over REVEL (DeMAG: 0.96 vs 0.96). DeMAG's gain was mainly in specificity: 0.82 vs 0.62 at default thresholds. Pooled AUROC is inflated by gene base rates: a gene-identity null reaches 0.921 (Arnoult 2026 preprint). Over AlphaMissense I expect **ΔAUROC within gene ≈ 0.00–0.02**. Any gain should concentrate in AD/GoF genes with dense, clustered P labels. The increment of 3D over 1D is probably smaller still. I found no published head-to-head of label propagation vs AlphaMissense.

---

## 1. Closest prior work: what each did, and the gap

### 1a. Label-using 3D neighbourhood predictors (closest)

**Luppino F, Adzhubei IA, Cassa CA, Toth-Petroczy A.** DeMAG predicts the effects of variants in clinically actionable genes by integrating structural and evolutionary epistatic features. *Nat Commun* 14:2230, 2023. doi:10.1038/s41467-023-37661-z. PMID 37076482. PMC10115847 **[V, full text read]**. **The single closest paper.**

*What they did:*
- A supervised classifier for 59 ACMG SF v2.0 genes, extended to 316 genes.
- **The "partners score" uses the clinical labels of partner residues.** Partners are residues within **< 11 Å Cα in AlphaFold2 models** (cutoffs 4–11 Å were tested and 11 Å was chosen) plus co-evolving pairs. Benign positions score −1, pathogenic +1, mixed or unannotated 0, and the scores are summed per residue.
- The partners score annotated 53% of previously unlabelled positions in the TP53 example.
- **The partners score was the most important feature (median importance ≈ 60%).** Removing the epistatic and structural features lowered CV specificity from 85% to 69% and sensitivity from 87% to 82%.
- Gene-specific MAF thresholds were used to add benign variants.

*Evaluation:*
- **Temporal:** ClinVar variants submitted **after Dec 2017** (852 P, 433 B). DeMAG AUC 0.96 vs REVEL 0.96 and VEST4 0.95. Specificity 0.82 vs REVEL 0.62; MCC 0.76 vs 0.69.
- Common-variant misclassification in the Estonian Biobank: 0.43 vs REVEL 0.51.

*Gap vs vus3d:*
- One-hop radius sum, no diffusion.
- No 1D-window comparator.
- No AlphaMissense; the paper predates AlphaMissense.
- No ACMG strength calibration.
- Pooled AUC; per-gene AUCs were shown only in CV.
- No explicit same-residue exclusion is described. I did not check the methods for this; same-position labels may feed the score **[UNVERIFIED]**.

**Sivley RM, Sheehan JH, Kropski JA, et al.** Three-dimensional spatial analysis of missense variants in RTEL1 identifies pathogenic variants in patients with Familial Interstitial Pneumonia. *BMC Bioinformatics* 19:18, 2018. doi:10.1186/s12859-018-2010-z. PMID 29361909 **[V]** (project notes).
- A "pathogenic proximity" (PathProx) score: 3D proximity to ClinVar P vs gnomAD neutral variants.
- Leave-one-out AUC 0.85 (7 P vs 29 neutral). 5/6 segregating VUS were predicted pathogenic.
- Single gene, homology model, no calibration.

**Moth CW, Sheehan JH, Al Mamun A, Sivley RM, et al.; Capra JA, Meiler J.** VUStruct: a compute pipeline for high throughput and personalized structural biology. *PLoS Comput Biol* 22(5):e1014183, 2026. doi:10.1371/journal.pcbi.1014183. PMID 42081565 (preprint doi:10.1101/2024.08.06.606224) **[V]**.
- Productionises PathProx ("pathogenicity prediction through spatial clustering analysis") for Undiagnosed Diseases Network cases (> 175).
- The VUStruct documentation says PathProx needs ≥ 3 pathogenic variants and that leave-one-out AUC is often weak **[V web page meilerlab.org/VUStruct, not peer reviewed]**.
- No benchmark against AlphaMissense and no calibration.

**Kroncke BM, Smith DK, Zuo Y, Glazer AM, Roden DM, Blume JD.** A Bayesian method to estimate variant-induced disease penetrance. *PLoS Genet* 16(6):e1008862, 2020. doi:10.1371/journal.pgen.1008862. PMID 32569262 **[V, full text grepped]**.
- For SCN5A–Brugada, a "**penetrance density**": the penetrance of 3D-neighbouring variants averaged into the prior (see S1 Text).
- Penetrance density and peak current carry orthogonal information.
- "The relatively small improvement in R² when all predictors are included suggests most information contained in the sequence-based predictive features is recapitulated by both peak current and penetrance density."
- The authors note that "penetrance density is based on the observed cases and controls", which is a circularity they handle by design.

The same group's follow-ups:
- **O'Neill MJ, et al.** Continuous Bayesian variant interpretation accounts for incomplete penetrance among Mendelian cardiac channelopathies. *Genet Med* 25(3):100355, 2023. doi:10.1016/j.gim.2022.12.002. PMID 36496179 **[V]**. KCNQ1, KCNH2 and SCN5A; heterogeneity within "hot spot" regions.
- **Dauda K, Yamauchi K, Ku M, et al.; Kroncke BM.** Structural evaluation of RYR2-CPVT missense variants and continuous Bayesian estimates of their penetrance. *Circ Genom Precis Med* 19(2):e005457, 2026. doi:10.1161/CIRCGEN.125.005457. PMID 41889371 **[V]**. Spearman with observed penetrance: Bayesian prior 0.19 vs AlphaMissense 0.18, REVEL 0.16, ClinVar 0.083. This is **a direct example of a tiny increment over AlphaMissense.**
- Kroncke BM, et al. Exploiting ion channel structure to assess rare variant pathogenicity. *Heart Rhythm* 15(6):890–894, 2018. doi:10.1016/j.hrthm.2018.01.021. PMID 29330128 **[V title]**.

*Gap:* single genes and penetrance (cases and controls) rather than ClinVar P/B. No diffusion, no PM1 calibration.

### 1b. 1D clustering / label-density predictors (the 1D comparator family)

**Quinodoz M, Peter VG, Cisarova K, et al.** Analysis of missense variants in the human genome reveals widespread gene-specific clustering and improves prediction of pathogenicity. *Am J Hum Genet* 109(3):457–470, 2022. doi:10.1016/j.ajhg.2022.01.006. PMID 35120630 **[V, full text read]** (project notes).
- *Positional score:* a per-transcript random forest on position, for transcripts with ≥ 10 PLP and ≥ 1 BLB. Trained on ClinVar of 2020-11-21 plus gnomAD variants with AF > max-AF-PLP (BS1-like), adding 354,888 BLB variants.
- P clustering was found in 387/840 genes, B clustering in 172/840.
- *MutScore:* a random forest of the positional score plus 17 other features. The positional score was the top feature.
- **Temporal test set 1:** ClinVar entries present on 2021-09-19 but not in the 2020-11-21 training data (1,867 P / 459 B, as stated in the methods). AUC 0.937 vs REVEL 0.924 and VEST4 0.919.
- **Test set 2:** HGMD DM entries added from 2017 onwards vs gnomAD (14,327 / 13,248). AUC 0.915 vs REVEL 0.904.
- Trained tools lost more AUC on new vs old entries (−0.017) than untrained tools (−0.0024), p = 5.8×10⁻⁵. This is an explicit temporal-drift result.
- To avoid type-2 circularity in subsets, B variants were restricted to genes that also carry P variants.
- *Gap:* 1D only, no AlphaMissense, pooled AUC, no ACMG calibration.

**Pérez-Palma E, May P, Iqbal S, et al.; Lal D.** Identification of pathogenic variant enriched regions across genes and gene families. *Genome Res* 30(1):62–71, 2020. doi:10.1101/gr.252601.119. PMID 31871067 **[V]**.
- Gene-family alignment burden analysis found 465 PERs in 1,252 genes.
- ClinVar missense variants inside PERs have OR = 106 for P vs B.
- This is a **region annotation, not a calibrated likelihood ratio.** The OR is in-sample and pooled, so it overstates usable PM1 strength.

**Samocha KE, Kosmicki JA, Karczewski KJ, et al.; Daly MJ.** Regional missense constraint improves variant deleteriousness prediction. *bioRxiv* 2017. doi:10.1101/148353 **[V preprint, record]** (MPC).
- Successor: **Chao KR, et al.** The landscape of regional missense mutational intolerance quantified from 125,748 exomes. *bioRxiv* 2024. doi:10.1101/2024.04.11.588920. PMID 38645134 **[V title]**.
- Population constraint, not labels. A label-free 1D baseline.

**Wiel L, et al.** MetaDome: pathogenicity analysis of genetic variants through aggregation of homologous human protein domains. *Hum Mutat* 40(8):1030–1038, 2019. doi:10.1002/humu.23798. PMID 31116477 **[V title]**.
- Update: **Wiel L, et al.** MetaDome 2027. *Nucleic Acids Res* 2026. doi:10.1093/nar/gkag1102. PMID 42963945 **[V title]**.
- Related: Wiel L, et al. De novo mutation hotspots in homologous protein domains identify function-altering mutations in neurodevelopmental disorders. *Am J Hum Genet* 2023. doi:10.1016/j.ajhg.2022.12.001. PMID 36563679 **[V title]**.

**Zhang X, Theotokis PI, Li N, et al.; Whiffin N, Ware JS.** Genetic constraint at single amino acid resolution in protein domains improves missense variant prioritisation and gene discovery. *Genome Med* 16:88, 2024. doi:10.1186/s13073-024-01358-9. PMID 38992748 **[V]** (Homologous Missense Constraint score).

**Turner TN, et al.** 2015 *HMG* (CLUMP, 1D) and **Lal D, et al.** 2020 *Genome Med* (gene-family) are **[V]** (project notes).

### 1c. Data-driven hotspots mapped to PM1 / ACMG strength (the calibration precedent)

**Walsh R, Mazzarotto F, Whiffin N, et al.; Ware JS.** Quantitative approaches to variant classification increase the yield and precision of genetic testing in Mendelian diseases: the case of hypertrophic cardiomyopathy. *Genome Med* 11:5, 2019. doi:10.1186/s13073-019-0616-z. PMID 30696458 **[V]**.
- Case–control etiological fraction (EF) and clustering identify regions enriched in cases (≤ 6,179 HCM cases).
- "Variant location was a better discriminator of pathogenicity than generic computational functional prediction algorithms."
- Non-truncating classes with EF ≥ 0.95 were found in 5 genes, giving 14–20% more actionable cases.
- They propose ACMG adaptations, i.e. an EF → evidence mapping.

**Waring A, Harper A, Salatino S, et al.; Farrall M.** Data-driven modelling of mutational hotspots and in silico predictors in hypertrophic cardiomyopathy. *J Med Genet* 58(8):556–564, 2021. doi:10.1136/jmedgenet-2020-106922. PMID 32732227 **[V]**.
- BIN-test / ClusterBurden, and a GAM "hotspot+" model that combines burden, clustering and in-silico predictors.
- It "effectively integrate[s] ACMG criteria PM1 and PP3 to yield strong or moderate evidence of pathogenicity for 31.8% of examined VUS".
- **This is the direct precedent for a joint PM1+PP3 model.** It is 1D, case–control, and single-disease.

**Whiffin N, Walsh R, Govind R, et al.** CardioClassifier: disease- and gene-specific computational decision support for clinical genome interpretation. *Genet Med* 20(10):1246–1254, 2018. doi:10.1038/gim.2017.258. PMID 29369293 **[V]**. Gene- and disease-specific automation of 14 computational ACMG rules, including data-derived hotspot rules.

**Yuan X, Xia X, Li J, Zhao G.** HCSeeker: a classification tool for human genetic variant hot and cold spots designed for PM1 and benign criteria in the ACMG-AMP guideline. *Genet Med* 27(12):101591, 2025. doi:10.1016/j.gim.2025.101591. PMID 41017653 **[V abstract]**.
- KDE + EM on variant positions found 988 hot spots and 682 cold spots across 889 genes, with a public database.
- **The abstract reports no LR calibration and no temporal validation.** I did not check the full text **[UNVERIFIED beyond abstract]**.
- **This is the most direct "automated PM1" competitor (1D).** Its cold-spot idea maps onto our B-label propagation.

**Liu S, Feng X, Wu Y, Bu F.** Calibration and refinement of ACMG/AMP criteria for variant classification with BayesQuantify. *J Med Genet* 2025. doi:10.1136/jmg-2025-110863. PMID 40973199 **[V]**.
- An R package that calibrates categorical or continuous evidence to ACMG strengths with LR+ / local LR+ and a bootstrap 95% lower bound.
- They replicated the PP3/BP4 thresholds; PM2 reached Supporting or Moderate depending on the prior.
- **A ready tool for the PM1 calibration step.**

**ClinGen VCEP PM1 definitions (expert hotspot rules = a mandatory baseline, and a leakage source)**
- **PTEN:** Mester JL, et al. Gene-specific criteria for PTEN variant curation. *Hum Mutat* 39:1581–1592, 2018. doi:10.1002/humu.23636. PMID 30311380 **[V, full text]**. PM1 = residues in catalytic motifs 90–94, 123–130 and 166–168 (NP_000305.3).
- **TP53:** Fortuno C, et al. Specifications of the ACMG/AMP variant interpretation guidelines for germline TP53 variants. *Hum Mutat* 42:223–236, 2021. doi:10.1002/humu.24152. PMID 33300245 **[V, full text]**.
  - The "functional domain without benign variation" arm of PM1 is not applicable.
  - PM1 applies at hotspot codons 175, 245, 248, 249, 273 and 282, plus somatic hotspots with ≥ 10 occurrences in cancerhotspots.org.
  - **PS1 cannot be used where PM1 has been applied.**
  - Update: Fortuno C, Frone MN, et al. *Genome Med* 17:128, 2025. doi:10.1186/s13073-025-01536-3. PMID 41126324 **[V]**. Likelihood-ratio-based strength modifications.
- **MYH7:** Kelly MA, et al. *Genet Med* 20:351–359, 2018. doi:10.1038/gim.2017.218. PMID 29300372 **[V title]**. The PM1 region definition (head domain residues) is **[UNVERIFIED]**.
- **RASopathy:** Gelb BD, et al. *Genet Med* 20:1334–1345, 2018. doi:10.1038/gim.2018.3. PMID 29493581 **[V title]**. Update: Wilcox EH, et al. *Genet Med Open* 3:103430, 2025. doi:10.1016/j.gimo.2025.103430. PMID 40496714 **[V]**.
- **RYR1 (MHS):** Johnston JJ, et al. *Genet Med* 23:1288–1295, 2021. doi:10.1038/s41436-021-01125-w. PMID 33767344 **[V title]**.

### 1d. 3D clustering without per-variant label transfer (context; mostly project notes)

- **Nason E, Gerges S, Satterstrom FK, et al.; Daly MJ, Finucane HK.** Systematic identification of disease-associated 3D neighborhoods in protein structures. *medRxiv* 2026. doi:10.64898/2026.05.29.26354366 **[V preprint]** (project notes).
  - The 3DNT test uses 15 Å Fisher spheres with label permutation.
  - ClinVar P/LP vs B/LB gave significant clusters in 872 genes.
  - **It is a discovery test, not a calibrated variant classifier.** This is the newest and highest-profile competitor; position vus3d as "turning 3DNT-type neighbourhoods into calibrated per-variant evidence".
- **Sivley 2018 AJHG** **[V]** (project notes). Only 12% of proteins with significant 3D clustering were found by 1D analyses. **This is the strongest prior argument that 3D could add over 1D.**
- **Iqbal 2023 Brain** (Essential3D) **[V]**: missense at Essential3D sites is 8-fold enriched in NDD patients. **Brünger 2023 Brain** **[V]**: 3D spatial distance scoring of pathogenic vs population variants in 30 ion channels. Both are in the project notes.
- **Cancer 3D hotspot tools:**
  - HotMAPS: Tokheim 2016 *Cancer Res* **[V title]** (project notes).
  - mutation3D: Meyer MJ, et al. *Hum Mutat* 37(5):447–456, 2016. doi:10.1002/humu.22963. PMID 26841357 **[V title]**.
  - 3D Hotspots: Gao J, et al. 3D clusters of somatic mutations in cancer reveal numerous rare mutations as functional targets. *Genome Med* 9:4, 2017. doi:10.1186/s13073-016-0393-x. PMID 28115009 **[V title]**.
  - HotSpot3D: Niu B, et al. *Nat Genet* 48:827–837, 2016. doi:10.1038/ng.3586. PMID 27294619 **[V title]**.
  - CLUMPS: Kamburov 2015 **[V]** (project notes).
  - These are somatic recurrence tools, not germline P/B labels.
- **Population-constraint 3D:**
  - Hicks 2019 *PNAS* **[V]** (project notes).
  - MTR3D: Silk M, et al. *Nucleic Acids Res* 49:W438–W445, 2021. doi:10.1093/nar/gkab428. PMID 34050760 **[V title]**.
  - MTR3D-AF2: Kovacs AS, et al. *Protein Sci* 33:e5112, 2024. doi:10.1002/pro.5112. PMID 39031445 **[V title]**.
  - COSMIS: Li B, Roden DM, Capra JA. The 3D mutational constraint on amino acid sites in the human proteome. *Nat Commun* 13:3273, 2022. doi:10.1038/s41467-022-30936-x. PMID 35672414 **[V title]**.
  - These are label-free 3D baselines worth including as a "3D without labels" control.
- **Per-variant structural annotators:** Missense3D (Ittisoponpisan 2019, *J Mol Biol* 431:2197–2212, PMID 30995449), MISCAST (Iqbal 2020, *NAR* 48:W132, PMID 32402084) and VarSite (Laskowski RA, et al. *Protein Sci* 29:111–119, 2020. doi:10.1002/pro.3746. PMID 31606900) are all **[V title]**. They describe structural damage and do not transfer labels.
- **Proteins with complexes and functional constraints:** Li B, Jin B, Capra JA, Bush WS. Integration of protein structure and population-scale DNA sequence data for disease gene discovery and variant interpretation. *Annu Rev Biomed Data Sci* 2022. doi:10.1146/annurev-biodatasci-122220-112147. PMID 35508071 **[V title]**. A review; useful for the introduction.

### 1e. Structure-based supervised or deep predictors (not label propagation)

- **Gao H, Hamp T, Ede J, et al.; Farh KK.** The landscape of tolerated genetic variation in humans and primates. *Science* 380(6648):eabn8153, 2023. doi:10.1126/science.abn8197 (Crossref: article number eabn8153, DOI suffix abn8197). PMID 37262156 **[V title]**. PrimateAI-3D: a voxelised 3D CNN trained on primate common variants, **with no ClinVar labels**.
- **Cheng J, Novati G, Pan J, et al.** AlphaMissense. *Science* 381:eadg7492, 2023. doi:10.1126/science.adg7492. PMID 37733863 **[V]** (project notes). Trained on population frequency "without explicitly training on" clinical data. Its classification thresholds were set on ClinVar **[UNVERIFIED detail]**.
- **Schmidt A, Röner S, Mai K, et al.; Kircher M, Ludwig KU.** Predicting the pathogenicity of missense variants using features derived from AlphaFold2 (AlphScore). *Bioinformatics* 39(5):btad280, 2023. doi:10.1093/bioinformatics/btad280. PMID 37084271 **[V]**. Alone it is worse than CADD/REVEL; added to them it gives modest gains.
- **Zhao H, Du H, Zhao S, et al.** SIGMA leverages protein structural information to predict the pathogenicity of missense variants. *Cell Rep Methods* 4(1):100687, 2024. doi:10.1016/j.crmeth.2023.100687. PMID 38211594 **[V]**. AF2 residue-graph features; AUC 0.933, and SIGMA+ (combined) 0.966. The test is pooled.
- **Danner M, Begemann M, Elbracht M, Kurth I, Krause J.** Utilizing protein structure graph embeddings to predict the pathogenicity of missense variants. *NAR Genom Bioinform* 7(3):lqaf097, 2025. doi:10.1093/nargab/lqaf097. PMID 40708850 **[V]**. Graph-autoencoder embeddings of ESMFold structures; they enhance CADD. This is a **GNN on residue graphs, but the labels train a classifier; they are not propagated within the protein.**
- **Banerjee A, et al.; Bahar I.** Accurate identification and mechanistic evaluation of pathogenic missense variants with Rhapsody-2. *PNAS* 122:e2418100122, 2025. doi:10.1073/pnas.2418100122. PMID 40314982 **[V title]**. Rhapsody (2020) and ENM dynamics are in the project notes.
- **Stein D, et al.** LoGoFunc. *Genome Med* 15:103, 2023. doi:10.1186/s13073-023-01261-9. PMID 38037155 **[V title]**.
- **Munro D, Singh M.** DeMaSk: a deep mutational scanning substitution matrix and its use for variant impact prediction. *Bioinformatics* 36:5322–5329, 2021 (PubMed year 2021). doi:10.1093/bioinformatics/btaa1030. PMID 33325500 **[V title]**.
- **Wu Y, Li R, Sun S, Weile J, Roth FP.** Improved pathogenicity prediction for rare human missense variants (VARITY). *Am J Hum Genet* 108:1891–1906, 2021. doi:10.1016/j.ajhg.2021.08.012. PMID 34551312 **[V title]**.
- **Cheng J, Lawrence C, Niepert M.** VEGN: Variant Effect Prediction with Graph Neural Networks. *arXiv* 2106.13642, 2021 **[V arXiv record]**. The graph is a gene–variant graph, not a residue contact graph.
- **Gene-specific structure plus AlphaMissense for PP3 (ClinGen-adjacent):**
  - Ramadane-Morchadi L, et al. ACMG/AMP interpretation of BRCA1 missense variants: structure-informed scores add evidence strength granularity to the PP3/BP4 computational evidence. *Am J Hum Genet* 112(5):993–1002, 2025. doi:10.1016/j.ajhg.2024.12.011. PMID 40233743 **[V]**. Likelihood ratios come from the BRIDGES case–control data; AlphaMissense outperforms ΔΔG and BayesDel; BP4 depends on RSA ≤ 60%.
  - Rotenberg N, et al. Integration of protein stability and AlphaMissense scores improves bioinformatic impact prediction for p53 missense and in-frame amino acid deletion variants. *Am J Hum Genet* 112(5):1003–1014, 2025. doi:10.1016/j.ajhg.2025.01.012. PMID 40233742 **[V]**.
  - Together these are the **template for "LR-calibrated structure + AlphaMissense" in a ClinGen setting**, but they use PP3, not PM1, and do not transfer labels.
- **Structure-based PP3/BP4 and PM1-like evidence (POLE/POLD1):** Arnaudi M, et al.; Papaleo E. Interpreting the effects of DNA polymerase variants at the structural level using MAVISp and molecular dynamics simulations. *bioRxiv* 2025. doi:10.1101/2025.08.04.653638 **[V preprint, abstract]**. The search snippet says it applies "structure-based PM1 supporting evidence" **[UNVERIFIED in full text]**.

### 1f. Searched, nothing found

- "label propagation" + variant pathogenicity + protein structure: 0 relevant hits on PubMed or Europe PMC; the only hit, metPropagate 2020, is gene-level.
- GNN trained on ClinVar labels over AlphaFold residue graphs with label transfer: none found. The closest are SIGMA, Danner 2025 and VEGN, which are not label transfer.
- "POSTMAN", "PathoGN": no such methods found.
- A ClinGen SVI document dedicated to PM1 calibration: none found as of 2026-11. The Pejaver 2022 cap on PP3+PM1 is the only SVI-level PM1 rule I could verify (§3e).
- GitHub projects (Missense-Structural-Explorer, SPHARE3D) implement 3D hotspot/PM1 ideas but are not peer reviewed. Cite only if needed.

---

## 2. Expected effect size (what prior work implies)

| Source | Comparison | Gain | Notes |
|---|---|---|---|
| MutScore (Quinodoz 2022) **[V]** | MutScore vs REVEL, temporal ClinVar set | **+0.013 AUROC** (0.937 vs 0.924) | Pooled; before AlphaMissense; includes 17 non-positional features |
| MutScore **[V]** | HGMD post-2017 vs gnomAD | **+0.011** (0.915 vs 0.904) | Pooled |
| DeMAG (Luppino 2023) **[V]** | vs REVEL, ClinVar after Dec 2017 | **AUC 0.96 vs 0.96 (≈0)**; specificity +0.20 at default thresholds; MCC +0.07 | The gain is in calibration/threshold behaviour, not ranking |
| DeMAG CV ablation **[V]** | Without partners and structure features | Specificity 85→69%, sensitivity 87→82% | In-sample CV; label leakage is likely to inflate it |
| AlphScore (Schmidt 2023) **[V]** | Added to CADD/REVEL | "Performance increased" (modest) | 3D features without labels |
| SIGMA (Zhao 2024) **[V]** | SIGMA+ combined | 0.933 → 0.966 | Pooled; possible leakage |
| Dauda 2026, RYR2 **[V]** | Bayesian prior (structure + clinical) vs AlphaMissense | Spearman **0.19 vs 0.18** | Single gene; penetrance outcome |
| Kroncke 2020, SCN5A **[V]** | Penetrance density + function + in-silico tools | "Relatively small improvement in R² when all predictors are included" | |
| Bergquist 2025 **[V]** | AlphaMissense/ESM1b/VARITY vs older tools | "At best only modest improvements" in the evidence/FPR trade-off | Calibration benchmark |
| Arnoult 2026 **[V preprint]** | Gene-identity null (gene P fraction) | **AUROC 0.921**; 4 predictors exceed it by only **0.036–0.044** | Within-gene evaluation **inverts** the ranking: AlphaMissense rises to 1st. Leave-one-out correction is worth only 0.002 |
| Isakov 2024 **[V preprint]** | ClinVar variants created after 2020 | "Mild yet significant decrease" for tools trained on older data | Temporal drift |

**Interpretation for vus3d.**
- (1) The documented pooled gains of label-density features over strong meta-predictors are ~0.01 AUROC. These gains came before AlphaMissense and use pooled metrics that partly reward gene base rates.
- (2) AlphaMissense already encodes much 3D context. Over AlphaMissense, a realistic **within-gene ΔAUROC is 0–0.02**, larger only in AD/GoF genes with many clustered P labels (e.g. RASopathy genes, MYH7, ion channels, TP53 DBD).
- (3) **No paper I found reports that clustering adds nothing once AlphaMissense is included**, but no paper tested it either. The Dauda 2026 RYR2 result (0.19 vs 0.18) and Bergquist's "modest" are the nearest evidence for a small increment.
- (4) The 3D-over-1D increment is the novel and riskiest claim. Sivley 2018 suggests 3D-only clusters exist (88% of 3D-clustered proteins were not found in 1D). The PhenoBridge results suggest AlphaFold context often reduces to sequence position. **Power the study for ΔAUROC ≈ 0.005–0.01 between 3D and 1D label scores**, and pre-register equivalence bounds.
- (5) The plausible clinical value is in **evidence yield and specificity**: more VUS reaching calibrated PM1_Supporting or Moderate, and fewer false PM1 calls in B-dense regions (the cold-spot side, cf. HCSeeker and DeMAG's specificity gain). Ranking gains are expected to be small.

---

## 3. Pitfalls

### 3a. Circularity and leakage
- **Grimm DG, Azencott CA, Aicheler F, et al.** The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat* 36(5):513–523, 2015. doi:10.1002/humu.22768. PMID 25684150 **[V title]** (project notes).
  - Type 1: the same variants in training and test.
  - Type 2: same-gene label structure.
  - **Label propagation exploits type-2 structure by design**, so the comparators must also use same-gene labels (MutScore-like 1D, DeMAG-like radius, and the gene-prior null).
- **PM1-coded labels.** VCEP and lab classifications in the test window often *used* PM1 hotspot rules derived from the same label clustering: PTEN catalytic motifs, TP53 hotspot codons, RASopathy and MYH7 regions. A hotspot score will "rediscover" curators' rules.
  - Mitigation: (i) stratify or exclude test variants whose ClinGen Evidence Repository record lists PM1 (erepo.clinicalgenome.org exposes per-criterion codes for VCEP-curated variants) **[UNVERIFIED API detail]**; (ii) report results separately for genes with and without a VCEP PM1 rule; (iii) include the **expert VCEP PM1 rule itself as a baseline**.
- **PP3 from AlphaMissense after 2023.** AlphaMissense was released Sep 2023. ClinGen calibrated it for PP3 in Bergquist (bioRxiv Sep 2024; *Genet Med* 2025). Test-window labels may therefore *contain* AlphaMissense evidence, which inflates AlphaMissense and its combinations in the temporal test. **Stratify the test set by first-classification date** (2023-01 → 2024-09 vs after), and by ≥ 2-star / VCEP status.
- **PS1/PM5 overlap.** Exclude same-residue labels when scoring, as planned. TP53 rules forbid PS1 together with PM1 (Fortuno 2021 **[V]**). Also run a sensitivity analysis excluding labels within ±1–2 residues in sequence, because adjacent-residue labels act as near-PM5.
- **PM1 + PP3 double counting.** Pejaver 2022 **[V full text]**: "we recommend that laboratories limit the sum of the evidence strength of PP3 and PM1 to strong … the same as limiting the sum of PP3 and PM1 to 4 points". Pejaver also notes that meta-predictors partly capture PM1/PS1/PM5 data. **Therefore calibrate the label score *conditional on* AlphaMissense** (the incremental LR within AlphaMissense bins), or calibrate the joint score as combined PP3+PM1 points capped at 4. Do not stack a separately calibrated PM1 on top of AlphaMissense-PP3.
- **Duplicated and reclassified entries.**
  - Collapse to one record per protein change: several nucleotide changes can give the same amino-acid change.
  - Use the variation ID's *first* P/LP or B/LB date from ClinVar archive VCF/XML snapshots, not "last evaluated".
  - Drop variants that were already P/B before the cutoff and changed later. Or analyse them separately: reclassification is ~15.5% at 5 years (Liang 2026, below) and P↔B flips are rare (Bonetti 2024).
  - Exclude conflicting interpretations.
- **Transcript and model mapping.** Map variants to the UniProt canonical isoform used by AlphaFold DB, and drop mismatches. Proteins > 2,700 aa come only as fragments (project notes). Low-pLDDT regions give spurious graph edges, and disordered regions are where AlphaMissense is weakest (project notes).

### 3b. Label imbalance and benign scarcity
- Many genes have **only P (or only B)** labels. In those genes a within-gene AUROC is undefined, and any pooled gain is a pure gene-prior effect (Arnoult 2026).
  - MutScore required ≥ 10 P and ≥ 1 B per transcript, and restricted B to genes with P to curb type-2 circularity **[V]**.
  - The propagation needs both labels to express a cold-spot signal. Report coverage: the number of genes with ≥ k P and ≥ k B before the cutoff.
- **gnomAD as a benign proxy.** Precedents:
  - Pejaver 2022 used gnomAD rare variants in Mendelian genes for the prior (0.0441) and for post-hoc checks **[V]**.
  - MutScore added gnomAD variants with AF > max-AF-PLP **[V]**.
  - DeMAG used gene-specific MAF thresholds **[V]**.
  - Risks: proxy benigns are common and concentrate in tolerant or disordered regions, so they are "easy" (Radjasandirane 2025). They also include recessive carriers and hypomorphs.
  - *Recommendation:* ClinVar B/LB only for the primary test. gnomAD-proxy B (AF-matched, excluding BA1/BS1-trivial variants) for training augmentation and as a sensitivity analysis. Never use allele-frequency features.
- **Test-set gene concentration.** Post-2023 labels are dominated by a few genes with VCEP or MAVE bulk curation (BRCA1/2, TP53, PTEN, MSH2, KCNH2, etc.), and those labels come from PS3/BS3 functional data. Report per-gene results and leave-one-gene-out influence.

### 3c. Evaluation: pooled vs per gene
- Pooled AUROC is inflated by gene base rates. The gene-identity null reaches 0.921, and within-gene evaluation changes which predictor ranks first (Arnoult 2026 **[V preprint]**).
- AlphaMissense's own paper reports gene-balanced evaluations **[UNVERIFIED detail]**.
- *Recommendation:* the primary metric is a **gene-stratified AUROC**. Count concordant P–B pairs only within the same gene, weighting genes by their number of P×B pairs, with a gene-level (cluster) bootstrap or a van Elteren-type stratified test. Report the unweighted mean per-gene AUROC over genes with ≥ 5 P and ≥ 5 B test labels as secondary. Always report the gene-prior null.

### 3d. Temporal design
Precedents for date-split ClinVar tests:
- Pejaver 2022: a ClinVar 2019 calibration set and a ClinVar 2020 set of variants added after Dec 2019 **[V]**. Bergquist 2025 reused this design **[V]**.
- MutScore (Nov 2020 → Sep 2021) **[V]**; DeMAG (after Dec 2017) **[V]**.
- Isakov 2024 (stratified by creation date) **[V preprint]**; Radjasandirane 2025 (ClinVar after 2021-05) **[V]**.
- RENOVO (Favalli 2021 *AJHG* 108:682–695. doi:10.1016/j.ajhg.2021.03.010. PMID 33761318 **[V]**): tested on variants whose interpretation changed. Bonetti E, et al. *J Transl Med* 22:713, 2024. doi:10.1186/s12967-024-05508-w. PMID 39085881 **[V]**: 16 quarterly ClinVar snapshots, 2020–2024; 82.6% of 10,196 meaningfully reclassified VUS were correctly predicted in 2020.

Reclassification dynamics:
- **Liang Y, Zhang G, Babb L, Rehm HL, et al.** Temporal analysis of clinical variant reclassification (2012–2025). *Genome Med* 2026. doi:10.1186/s13073-026-01805-9 **[V Crossref record; content from search snippet only: cumulative first-reclassification 2.5% at 1 y, 15.5% at 5 y — UNVERIFIED]**.
- **Harrison SM, Rehm HL.** Is 'likely pathogenic' really 90% likely? Reclassification data in ClinVar. *Genome Med* 11:72, 2019. doi:10.1186/s13073-019-0688-9. PMID 31752965 **[V title]**.
- **Sharo AG, et al.** ClinVar and HGMD genomic variant classification accuracy has improved over time, as measured by implied disease burden. *Genome Med* 15:51, 2023. doi:10.1186/s13073-023-01199-y. PMID 37443081 **[V title]**.
- **Fowler DM, Rehm HL.** Will variants of uncertain significance still exist in 2030? *Am J Hum Genet* 111:5–10, 2024. doi:10.1016/j.ajhg.2023.11.005. PMID 38086381 **[V title]**.

Other benchmarking cautions:
- **Livesey BJ, Marsh JA.** Variant effect predictor correlation with functional assays is reflective of clinical classification performance. *Genome Biol* 26:104, 2025. doi:10.1186/s13059-025-03575-w. PMID 40264194 **[V]**. DMS and clinical benchmarks agree "especially for predictors that have not been directly trained on human clinical variants". Earlier work: Livesey & Marsh *Mol Syst Biol* 16:e9380, 2020 (doi:10.15252/msb.20199380, PMID 32627955) and 19:e11474, 2023 (doi:10.15252/msb.202211474, PMID 37310135) **[V title]**.
  - **Suggestion:** a DMS/MAVE secondary validation. Use MAVE labels as P/B proxies in genes with dense pre-2023 ClinVar labels. This tests whether label propagation tracks function, not curation habits.
- **Gunning AC, et al.** Assessing performance of pathogenicity predictors using clinically relevant variant datasets. *J Med Genet* 58(8):547–555, 2021. doi:10.1136/jmedgenet-2020-107003. PMID 32843488 **[V]**. Performance is "substantially lower" on a clinically representative set: REVEL AUC 0.82.
- **Tabet DR, et al.** Benchmarking computational variant effect predictors by their ability to infer human traits. *Genome Biol* 25:172, 2024. doi:10.1186/s13059-024-03314-7. PMID 38951922 **[V title]**. A label-free benchmark.
- **Radjasandirane R, et al.** Insights for variant clinical interpretation based on a benchmark of 65 variant effect predictors. *Genomics* 117:111036, 2025. doi:10.1016/j.ygeno.2025.111036. PMID 40127826 (preprint doi:10.1101/2024.07.08.602580) **[V abstract]**. Most ClinVar variants are "easy to predict". The search snippet says nearby same-domain pathogenic variants make variants easy **[UNVERIFIED in full text]**.
- **Isakov O, Fluss R, Marek-Yagel D, Ben-Shachar S.** The impact of clinical and molecular variant properties on calibration and performance of variant effect prediction tools. *bioRxiv* 2024. doi:10.1101/2024.10.02.614907 **[V preprint]**.
- **Arnoult D, Harrizi S, Nait Irahal I, Mostafa K.** Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. *bioRxiv* 2026. doi:10.64898/2026.08.21.746180 **[V preprint, abstract]**. Releases `genenull`.

### 3e. Calibration (ClinGen SVI), with the numbers verified
- **Tavtigian SV, Greenblatt MS, Harrison SM, et al.** Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med* 20(9):1054–1060, 2018. doi:10.1038/gim.2017.210. PMID 29300386 **[V, full text]**.
  - With prior 0.10 and exponent X = 2, odds of pathogenicity are **Very Strong 350, Strong 18.7, Moderate 4.3 (≈ 4.33), Supporting 2.08**.
  - Benign odds: Supporting **0.48** and Strong **0.053**; Moderate is 0.23 by reciprocity **[derived]**.
  - **The values in the brief (2.08 / 4.33 / 18.7) are confirmed as Tavtigian's.**
- **Tavtigian SV, Harrison SM, Boucher KM, Biesecker LG.** Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Hum Mutat* 41(10):1734–1737, 2020. doi:10.1002/humu.24088. PMID 32720330 **[V]**. Points are proportional to log(odds). Supporting = 1, Moderate = 2, Strong = 4, Very Strong = 8; the 4-point cap is stated in Pejaver 2022.
- **Pejaver V, Byrne AB, Feng BJ, et al.** Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet* 109(12):2163–2177, 2022. doi:10.1016/j.ajhg.2022.10.013. PMID 36413997 **[V, full text]**.
  - Prior = **0.0441** (estimated from gnomAD); c = 1,124.
  - **LR+ thresholds: Supporting 2.406, Moderate 5.790, Strong 33.53, Very Strong 1,124.** Posterior thresholds: 0.0999 / 0.2108 / 0.6073 / 0.9811. BP4 uses the reciprocals: 0.4156 / 0.1727 / 0.0298.
  - Method: the local posterior in an adaptive ε-window; each threshold is the smallest score at which the **one-sided 95% bootstrap lower bound** of local LR+ stays ≥ the threshold (10,000 bootstraps).
  - **Note:** with Pejaver's prior the thresholds are *higher* than Tavtigian's 2.08 / 4.33 / 18.7. Pick one framework and state the prior. For a gene-specific PM1, consider gene-specific priors, and report both.
  - "Limit the sum of PP3 and PM1 to strong" **[V]**.
- **Bergquist T, Stenton SL, Nadeau EAW, et al.; Pejaver V.** Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med* 27(6):101402, 2025. doi:10.1016/j.gim.2025.101402. PMID 40084623 **[V, full text]**.
  - AlphaMissense point bins (−3 … +4): ≤ 0.070 (−3); 0.071–0.099 (−2); 0.100–0.169 (−1); 0.170–0.791 (0); **0.792–0.905 (+1, PP3_Supporting); 0.906–0.971 (+2, Moderate); 0.972–0.989 (+3); ≥ 0.990 (+4, Strong)**.
  - The column-to-point mapping is inferred from the table layout **[V numbers; mapping inferred]**.
  - **Use these as the AlphaMissense-alone comparator for evidence yield.**
- **Stenton SL, et al.** Assessment of the evidence yield for the calibrated PP3/BP4 computational recommendations. *Genet Med* 26(11):101213, 2024. doi:10.1016/j.gim.2024.101213. PMID 39030733 **[V]**. Per-proband yield: a median of 1 PP3_Strong and 3–5 PP3_Moderate. A template for reporting PM1 yield.
- **Brnich SE, et al.** Recommendations for application of the functional evidence PS3/BS3 criterion using the ACMG/AMP sequence variant interpretation framework. *Genome Med* 12:3, 2019. doi:10.1186/s13073-019-0690-2. PMID 31892348 **[V title]**. The OddsPath approach used for PS3/BS3.
- **Richards S, et al.** ACMG/AMP standards and guidelines. *Genet Med* 17(5):405–424, 2015. doi:10.1038/gim.2015.30. PMID 25741868 **[V title]**. The PM1 definition: "mutational hot spot and/or critical and well-established functional domain … without benign variation".
- **Tooling:** BayesQuantify (Liu 2025 **[V]**) implements LR+ / local LR+ calibration with a bootstrap lower bound.
- **Small-n per gene.** Per-gene LR estimation is unstable. Use a hierarchical or pooled calibration (calibrate on score quantiles across genes, with gene-level random effects) and report per-gene intervals. This mirrors the Pejaver/ClinGen caution about gene-specific applicability (Bergquist 2025: "gene- or disease-specific applicability and variability").

---

## 4. Recommended primary metric and comparison

**Primary hypothesis (pre-registered).** In the temporal test set (first P/LP or B/LB after 2023-01; ClinVar Nov 2026; ≥ 1-star, non-conflicting, unique protein change; same-residue labels excluded), adding the 3D-RWR label score to AlphaMissense improves within-gene discrimination over adding a matched 1D-window label score to AlphaMissense.

- **Primary metric:** the **gene-stratified AUROC** (within-gene P–B pair concordance pooled across genes, weighted by P×B pairs). Report **Δ = AUROC(AM + RWR-label) − AUROC(AM + 1D-label)** with a gene-cluster bootstrap 95% CI.
  - Fit all models only on pre-cutoff labels: a logistic stack of logit(AM) plus the label score.
  - Pre-specify an equivalence margin (e.g. ±0.005) so that a null result is interpretable.
- **Key secondary comparisons:**
  - AM + RWR vs AM alone (does the label score add anything?);
  - AM + RWR vs AM + 3D-radius (DeMAG-like, 8 and 11 Å) (does diffusion beat radius?);
  - AM + RWR vs the **expert VCEP PM1 rule** in genes that have one;
  - the gene-prior null (`genenull`) and label-free 3D constraint (COSMIS/MTR3D) as floors.
- **Co-primary clinical metric:** the calibrated incremental evidence.
  - Calibrate the RWR score *conditional on AM bins* on pre-cutoff data (Pejaver local posterior, prior 0.0441, one-sided 95% bootstrap bound). Then verify the realised LR+ per strength bin on the temporal set, as Pejaver did with the ClinVar 2020 set.
  - Report: the fraction of test P/LP reaching PM1_Supporting or Moderate; the false-PM1 rate among test B/LB; the net reclassifiable VUS under the combined PP3 + PM1 ≤ 4-point cap.
- **Stratifications:**
  - Date (before vs after Sep 2024, when AlphaMissense was calibrated for PP3);
  - VCEP or PM1-coded vs not;
  - AD/GoF vs AR/LoF genes;
  - pLDDT ≥ 70 vs < 70 at the query residue;
  - label density before the cutoff (bins of P and B counts within the graph neighbourhood);
  - genes with vs without both P and B labels.
- **Secondary validation free of curation circularity:** MAVE/DMS-derived functional classes in genes with pre-2023 ClinVar labels (Livesey & Marsh 2025 rationale).

---

## 5. Venues

| Outcome | Venue | Rationale |
|---|---|---|
| Positive calibrated PM1 tool with ClinGen-style validation | **Genetics in Medicine** / **Genet Med Open** | Published Pejaver-line calibrations (Bergquist 2025, Stenton 2024), HCSeeker (PM1) and VCEP specs |
| Positive, methodologically strong, multi-gene | **American Journal of Human Genetics** | Published MutScore, Pejaver 2022, BRCA1/TP53 structure+AM LR papers (2025). High bar after DeMAG and 3DNT |
| Positive or "3D ≈ 1D" benchmark with rigorous circularity analysis | **Genome Medicine** / **Genome Biology** | Published Walsh 2019, Livesey & Marsh 2025, temporal reclassification studies, LoGoFunc |
| Null or modest-increment benchmark, or a methods/software paper | **HGG Advances**, **Human Mutation**, **Bioinformatics**, **NAR Genomics & Bioinformatics**, **PLoS Comput Biol** | Accept careful negative or benchmark results; Human Mutation publishes VCEP specs and Bayesian frameworks |
| Preprint first | medRxiv (clinical) or bioRxiv (methods) | Competitors (3DNT, Arnoult 2026) are active preprints; establish priority |

---

## 6. Search log (for reproducibility)

- **PubMed:**
  - `PM1 AND (hotspot OR calibration) AND ACMG` (15 hits; HCSeeker, Waring found)
  - `PM1 AND (AlphaFold OR "protein structure" OR three-dimensional)` (10 hits; case reports only)
  - `ClinVar AND "graph neural network" AND (AlphaFold OR "protein structure") AND missense` (0 hits)
  - `label propagation AND (protein structure OR residue) AND variant pathogenicity` (1 irrelevant hit)
  - Author/title verification queries for every PMID above.
- **Europe PMC:**
  - `TITLE_ABS:("label propagation") AND variant AND pathogenic*` (1 hit: metPropagate)
  - hotspot + AlphaFold/structure + ACMG/PM1/classification (7 hits)
  - AlphaMissense + cluster/hotspot/spatial + ClinVar + P&B (2 hits; gene-specific)
  - temporal/time-split/prospective + ClinVar + VEP/AlphaMissense (4 hits)
  - VUStruct, DeMAG.
- **Web:** 3D clustering + ClinVar + AlphaFold neighbourhood score; PM1 + AlphaFold + calibrated LR; GNN + residue graph + ClinVar; ClinGen SVI PM1 calibration; PathProx; temporal holdout + AlphaMissense; MutScore AUC; "clustering adds little" + AlphaMissense; germline hotspot PM1 3D preprint.
- **Full text read or grepped:** Pejaver 2022 (PMC9748256), Bergquist 2025 (PMC12208618), Tavtigian 2018 (PMC6336098), MutScore (PMC8948164), DeMAG (PMC10115847), Kroncke 2020 (PMC7347235), PTEN VCEP (PMC6329583), TP53 VCEP (PMC8374922).
- **Not fetched (HTTP 429):** Arnoult 2026 full text; only the abstract was read.
