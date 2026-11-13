"""Assemble paper.html: number {{key}} citations by first appearance and append the reference list."""
import re
from pathlib import Path
H = Path(__file__).resolve().parent
REFS = {
 "richards": "Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. <i>Genet Med</i> 17(5):405&ndash;424, 2015. doi:10.1038/gim.2015.30",
 "fowler": "Fowler DM, Rehm HL. Will variants of uncertain significance still exist in 2030? <i>Am J Hum Genet</i> 111:5&ndash;10, 2024. doi:10.1016/j.ajhg.2023.11.005",
 "pejaver": "Pejaver V, Byrne AB, Feng BJ, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. <i>Am J Hum Genet</i> 109(12):2163&ndash;2177, 2022. doi:10.1016/j.ajhg.2022.10.013",
 "bergquist": "Bergquist T, Stenton SL, Nadeau EAW, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. <i>Genet Med</i> 27(6):101402, 2025. doi:10.1016/j.gim.2025.101402",
 "alphamissense": "Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. <i>Science</i> 381:eadg7492, 2023. doi:10.1126/science.adg7492",
 "sivley18": "Sivley RM, Dou X, Meiler J, Bush WS, Capra JA. Comprehensive analysis of constraint on the spatial distribution of missense variants in human protein structures. <i>Am J Hum Genet</i> 102(3):415&ndash;426, 2018. doi:10.1016/j.ajhg.2018.01.017",
 "demag": "Luppino F, Adzhubei IA, Cassa CA, Toth-Petroczy A. DeMAG predicts the effects of variants in clinically actionable genes by integrating structural and evolutionary epistatic features. <i>Nat Commun</i> 14:2230, 2023. doi:10.1038/s41467-023-37661-z",
 "vustruct": "Moth CW, Sheehan JH, Al Mamun A, et al. VUStruct: a compute pipeline for high throughput and personalized structural biology. <i>PLoS Comput Biol</i> 22(5):e1014183, 2026. doi:10.1371/journal.pcbi.1014183",
 "mutscore": "Quinodoz M, Peter VG, Cisarova K, et al. Analysis of missense variants in the human genome reveals widespread gene-specific clustering and improves prediction of pathogenicity. <i>Am J Hum Genet</i> 109(3):457&ndash;470, 2022. doi:10.1016/j.ajhg.2022.01.006",
 "hcseeker": "Yuan X, Xia X, Li J, Zhao G. HCSeeker: a classification tool for human genetic variant hot and cold spots designed for PM1 and benign criteria in the ACMG-AMP guideline. <i>Genet Med</i> 27(12):101591, 2025. doi:10.1016/j.gim.2025.101591",
 "kroncke": "Kroncke BM, Smith DK, Zuo Y, Glazer AM, Roden DM, Blume JD. A Bayesian method to estimate variant-induced disease penetrance. <i>PLoS Genet</i> 16(6):e1008862, 2020. doi:10.1371/journal.pgen.1008862",
 "walsh": "Walsh R, Mazzarotto F, Whiffin N, et al. Quantitative approaches to variant classification increase the yield and precision of genetic testing in Mendelian diseases: the case of hypertrophic cardiomyopathy. <i>Genome Med</i> 11:5, 2019. doi:10.1186/s13073-019-0616-z",
 "nason": "Nason E, Gerges S, Satterstrom FK, et al. Systematic identification of disease-associated 3D neighborhoods in protein structures. <i>medRxiv</i> 2026 (preprint). doi:10.64898/2026.05.29.26354366",
 "grimm": "Grimm DG, Azencott CA, Aicheler F, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. <i>Hum Mutat</i> 36(5):513&ndash;523, 2015. doi:10.1002/humu.22768",
 "arnoult": "Arnoult D, Harrizi S, Nait Irahal I, Mostafa K. Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity predictors. <i>bioRxiv</i> 2026 (preprint). doi:10.64898/2026.08.21.746180",
 "clinvar": "Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. <i>Nucleic Acids Res</i> 46(D1):D1062&ndash;D1067, 2018. doi:10.1093/nar/gkx1153",
 "mane": "Morales J, Pujar S, Loveland JE, et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. <i>Nature</i> 604:310&ndash;315, 2022. doi:10.1038/s41586-022-04558-8",
 "vep": "McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. <i>Genome Biol</i> 17:122, 2016. doi:10.1186/s13059-016-0974-4",
 "afdb": "Varadi M, Bertoni D, Magana P, et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. <i>Nucleic Acids Res</i> 52(D1):D368&ndash;D375, 2024. doi:10.1093/nar/gkad1011",
 "kohler": "K&ouml;hler S, Bauer S, Horn D, Robinson PN. Walking the interactome for prioritization of candidate disease genes. <i>Am J Hum Genet</i> 82(4):949&ndash;958, 2008. doi:10.1016/j.ajhg.2008.02.013",
 "tavtigian": "Tavtigian SV, Greenblatt MS, Harrison SM, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. <i>Genet Med</i> 20(9):1054&ndash;1060, 2018. doi:10.1038/gim.2017.210",
 "pten": "Mester JL, Ghosh R, Pesaran T, et al. Gene-specific criteria for PTEN variant curation: recommendations from the ClinGen PTEN Expert Panel. <i>Hum Mutat</i> 39:1581&ndash;1592, 2018. doi:10.1002/humu.23636",
 "tp53": "Fortuno C, Lee K, Olivier M, et al. Specifications of the ACMG/AMP variant interpretation guidelines for germline TP53 variants. <i>Hum Mutat</i> 42:223&ndash;236, 2021. doi:10.1002/humu.24152",
 "rasopathy": "Gelb BD, Cav&eacute; H, Dillon MW, et al. ClinGen's RASopathy Expert Panel consensus methods for variant interpretation. <i>Genet Med</i> 20:1334&ndash;1345, 2018. doi:10.1038/gim.2018.3",
 "dauda": "Dauda K, Yamauchi K, Ku M, et al. Structural evaluation of RYR2-CPVT missense variants and continuous Bayesian estimates of their penetrance. <i>Circ Genom Precis Med</i> 19(2):e005457, 2026. doi:10.1161/CIRCGEN.125.005457",
}
body = (H / "body.html").read_text()
order = []
def num(m):
    keys = re.findall(r"\{\{(\w+)\}\}", m.group(0))
    ns = []
    for k in keys:
        if k == "REFERENCES":
            return m.group(0)
        if k not in REFS:
            raise KeyError(k)
        if k not in order:
            order.append(k)
        ns.append(order.index(k) + 1)
    ns = sorted(set(ns)); parts = []; i = 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
            j += 1
        parts.append(f"{ns[i]}&ndash;{ns[j]}" if j - i >= 2 else ",".join(map(str, ns[i:j + 1])))
        i = j + 1
    return "[" + ",".join(parts) + "]"
body = re.sub(r"(?:\{\{\w+\}\})+", num, body)
unused = set(REFS) - set(order)
assert not unused, unused
refs = "<ol class=\"refs\">\n" + "\n".join(f"<li>{REFS[k]}.</li>" for k in order) + "\n</ol>"
body = body.replace("{{REFERENCES}}", refs)
(H / "paper.html").write_text((H / "_head.html").read_text() + "\n" + body)
print(len(order), "references")
