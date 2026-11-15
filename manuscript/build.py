"""Assemble paper.html from body.html + style.css.

- {{key}} citations become superscript numbers in order of first use (AMA style).
- {{REFERENCES}} becomes the numbered reference list; every reference must be cited.
- Word counts (abstract; Introduction to Conclusions) are printed for journal limits.
All DOIs below were checked against Crossref on 2026-11-15 (books 'kish' and 'donner' have no DOI).
"""
import html
import re
from pathlib import Path

H = Path(__file__).resolve().parent
REFS = {
    "wakap": "Nguengang Wakap S, Lambert DM, Olry A, et al. Estimating cumulative point prevalence of rare diseases: analysis of the Orphanet database. <i>Eur J Hum Genet</i>. 2020;28(2):165-173. doi:10.1038/s41431-019-0508-0",
    "richards": "Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. <i>Genet Med</i>. 2015;17(5):405-424. doi:10.1038/gim.2015.30",
    "jacobsen": "Jacobsen JOB, Baudis M, Baynam GS, et al. The GA4GH Phenopacket schema defines a computable representation of clinical data. <i>Nat Biotechnol</i>. 2022;40(6):817-820. doi:10.1038/s41587-022-01357-4",
    "hpo": "K&ouml;hler S, Gargano M, Matentzoglu N, et al. The Human Phenotype Ontology in 2021. <i>Nucleic Acids Res</i>. 2021;49(D1):D1207-D1217. doi:10.1093/nar/gkaa1043",
    "danis": "Danis D, Bamshad MJ, Bridges Y, et al. A corpus of GA4GH phenopackets: case-level phenotyping for genomic diagnostics and discovery. <i>HGG Adv</i>. 2025;6(1):100371. doi:10.1016/j.xhgg.2024.100371",
    "gpsea": "Rekerle L, Danis D, Rehburg F, et al. GA4GH phenopacket-driven characterization of genotype-phenotype correlations in Mendelian disorders. <i>Am J Hum Genet</i>. 2026;113(1):57-70. doi:10.1016/j.ajhg.2025.12.001",
    "jumper": "Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. <i>Nature</i>. 2021;596(7873):583-589. doi:10.1038/s41586-021-03819-2",
    "varadi": "Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. <i>Nucleic Acids Res</i>. 2022;50(D1):D439-D444. doi:10.1093/nar/gkab1061",
    "alphamissense": "Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. <i>Science</i>. 2023;381(6664):eadg7492. doi:10.1126/science.adg7492",
    "sivley": "Sivley RM, Dou X, Meiler J, Bush WS, Capra JA. Comprehensive analysis of constraint on the spatial distribution of missense variants in human protein structures. <i>Am J Hum Genet</i>. 2018;102(3):415-426. doi:10.1016/j.ajhg.2018.01.017",
    "badonyi": "Badonyi M, Marsh JA. Prevalence of loss-of-function, gain-of-function and dominant-negative mechanisms across genetic disease phenotypes. <i>Nat Commun</i>. 2025;16:8392. doi:10.1038/s41467-025-63234-3",
    "brunger": "Br&uuml;nger T, Ivaniuk A, P&eacute;rez-Palma E, et al. Conserved missense variant pathogenicity and correlated phenotypes across paralogous genes. <i>Genome Biol</i>. 2025;26:197. doi:10.1186/s13059-025-03663-x",
    "wright": "Wright CF, West B, Tuke M, et al. Assessing the pathogenicity, penetrance, and expressivity of putative disease-causing variants in a population setting. <i>Am J Hum Genet</i>. 2019;104(2):275-286. doi:10.1016/j.ajhg.2018.12.015",
    "hurlbert": "Hurlbert SH. Pseudoreplication and the design of ecological field experiments. <i>Ecol Monogr</i>. 1984;54(2):187-211. doi:10.2307/1942661",
    "lazic": "Lazic SE. The problem of pseudoreplication in neuroscientific studies: is it affecting your analysis? <i>BMC Neurosci</i>. 2010;11:5. doi:10.1186/1471-2202-11-5",
    "kish": "Kish L. <i>Survey Sampling</i>. Wiley; 1965",
    "donner": "Donner A, Klar N. <i>Design and Analysis of Cluster Randomization Trials in Health Research</i>. Arnold; 2000",
    "vep": "McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. <i>Genome Biol</i>. 2016;17(1):122. doi:10.1186/s13059-016-0974-4",
    "mane": "Morales J, Pujar S, Loveland JE, et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. <i>Nature</i>. 2022;604(7905):310-315. doi:10.1038/s41586-022-04558-8",
    "tunya": "Tunyasuvunakool K, Adler J, Wu Z, et al. Highly accurate protein structure prediction for the human proteome. <i>Nature</i>. 2021;596(7873):590-596. doi:10.1038/s41586-021-03828-1",
    "chennubhotla": "Chennubhotla C, Bahar I. Markov propagation of allosteric effects in biomolecular systems: application to GroEL-GroES. <i>Mol Syst Biol</i>. 2006;2:36. doi:10.1038/msb4100075",
    "tong": "Tong H, Faloutsos C, Pan JY. Fast random walk with restart and its applications. In: <i>Proceedings of the Sixth International Conference on Data Mining (ICDM'06)</i>. IEEE; 2006:613-622. doi:10.1109/ICDM.2006.70",
    "cowen": "Cowen L, Ideker T, Raphael BJ, Sharan R. Network propagation: a universal amplifier of genetic associations. <i>Nat Rev Genet</i>. 2017;18(9):551-562. doi:10.1038/nrg.2017.38",
    "winkler": "Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE. Permutation inference for the general linear model. <i>Neuroimage</i>. 2014;92:381-397. doi:10.1016/j.neuroimage.2014.01.060",
    "phipson": "Phipson B, Smyth GK. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. <i>Stat Appl Genet Mol Biol</i>. 2010;9:Article39. doi:10.2202/1544-6115.1585",
    "bh": "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. <i>J R Stat Soc Series B Stat Methodol</i>. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x",
    "storey": "Storey JD, Tibshirani R. Statistical significance for genomewide studies. <i>Proc Natl Acad Sci U S A</i>. 2003;100(16):9440-9445. doi:10.1073/pnas.1530509100",
    "shackleton": "Shackleton S, Lloyd DJ, Jackson SN, et al. LMNA, encoding lamin A/C, is mutated in partial lipodystrophy. <i>Nat Genet</i>. 2000;24(2):153-156. doi:10.1038/72807",
    "basson": "Basson CT, Huang T, Lin RC, et al. Different TBX5 interactions in heart and limb defined by Holt-Oram syndrome mutations. <i>Proc Natl Acad Sci U S A</i>. 1999;96(6):2919-2924. doi:10.1073/pnas.96.6.2919",
    "vandevondel": "Van de Vondel L, De Winter J, Beijer D, et al. De novo and dominantly inherited SPTAN1 mutations cause spastic paraplegia and cerebellar ataxia. <i>Mov Disord</i>. 2022;37(6):1175-1186. doi:10.1002/mds.28959",
    "rabin": "Rabin R, Radmanesh A, Glass IA, et al. Genotype-phenotype correlation at codon 1740 of SETD2. <i>Am J Med Genet A</i>. 2020;182(9):2037-2048. doi:10.1002/ajmg.a.61724",
    "boutboul": "Boutboul D, Kuehn HS, Van de Wyngaert Z, et al. Dominant-negative IKZF1 mutations cause a T, B, and myeloid cell combined immunodeficiency. <i>J Clin Invest</i>. 2018;128(7):3071-3087. doi:10.1172/JCI98164",
    "xian": "Xian J, Parthasarathy S, Ruggiero SM, et al. Assessing the landscape of STXBP1-related disorders in 534 individuals. <i>Brain</i>. 2022;145(5):1668-1683. doi:10.1093/brain/awab327",
    "sirmaci": "Sirmaci A, Spiliopoulos M, Brancati F, et al. Mutations in ANKRD11 cause KBG syndrome, characterized by intellectual disability, skeletal malformations, and macrodontia. <i>Am J Hum Genet</i>. 2011;89(2):289-294. doi:10.1016/j.ajhg.2011.06.007",
    "low": "Low K, Ashraf T, Canham N, et al. Clinical and genetic aspects of KBG syndrome. <i>Am J Med Genet A</i>. 2016;170(11):2835-2846. doi:10.1002/ajmg.a.37842",
    "martinez": "Martinez-Cayuelas E, Blanco-Kelly F, Lopez-Grondona F, et al. Clinical description, molecular delineation and genotype-phenotype correlation in 340 patients with KBG syndrome: addition of 67 new patients. <i>J Med Genet</i>. 2023;60(7):644-654. doi:10.1136/jmg-2022-108632",
    "breslow": "Breslow NE, Clayton DG. Approximate inference in generalized linear mixed models. <i>J Am Stat Assoc</i>. 1993;88(421):9-25. doi:10.1080/01621459.1993.10594284",
    "ruff": "Ruff KM, Pappu RV. AlphaFold and implications for intrinsically disordered proteins. <i>J Mol Biol</i>. 2021;433(20):167208. doi:10.1016/j.jmb.2021.167208",
    "hubrich": "Hubrich D, Alvarado Valverde J, Lee CY, et al. Variant characterization in the intrinsically disordered human proteome. <i>Nat Struct Mol Biol</i>. 2026;33(8):1183-1193. doi:10.1038/s41594-026-01846-z",
    "firth": "Firth D. Bias reduction of maximum likelihood estimates. <i>Biometrika</i>. 1993;80(1):27-38. doi:10.1093/biomet/80.1.27",
    "deboer": "de Boer E, Ockeloen CW, Kampen RA, et al. Missense variants in ANKRD11 cause KBG syndrome by impairment of stability or transcriptional activity of the encoded protein. <i>Genet Med</i>. 2022;24(10):2051-2064. doi:10.1016/j.gim.2022.06.007",
    "lindeboom": "Lindeboom RGH, Supek F, Lehner B. The rules and impact of nonsense-mediated mRNA decay in human cancers. <i>Nat Genet</i>. 2016;48(10):1112-1118. doi:10.1038/ng.3664",
    "walz": "Walz K, Cohen D, Neilsen PM, et al. Characterization of ANKRD11 mutations in humans and mice related to KBG syndrome. <i>Hum Genet</i>. 2015;134(2):181-190. doi:10.1007/s00439-014-1509-2",
}


def words(fragment):
    text = re.sub(r"<(table|figure|div class=\"box\").*?</\1>", " ", fragment, flags=re.S)
    text = html.unescape(re.sub(r"<[^>]+>", " ", text))
    return len(re.findall(r"[A-Za-z0-9][\w.,%'–−-]*", text))


def main():
    body = (H / "body.html").read_text()
    order = []

    def cite(m):
        keys = re.findall(r"\{\{(\w+)\}\}", m.group(0))
        if any(k.isupper() or k.startswith("WC_") for k in keys):
            return m.group(0)
        nums = []
        for k in keys:
            if k not in REFS:
                raise KeyError(k)
            if k not in order:
                order.append(k)
            nums.append(order.index(k) + 1)
        nums = sorted(set(nums))
        parts, i = [], 0
        while i < len(nums):
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            parts.append(f"{nums[i]}&ndash;{nums[j]}" if j - i >= 2 else ",".join(map(str, nums[i:j + 1])))
            i = j + 1
        return f'[{",".join(parts)}]'

    # LaTeX/Vancouver bracket style: citation before punctuation, preceded by a space
    body = re.sub(r"([.,;:])((?:\{\{[a-z]\w*\}\})+)", r" \2\1", body)
    body = re.sub(r"(?<=[^\s>(\[}])((?:\{\{[a-z]\w*\}\})+)", r" \1", body)
    body = re.sub(r"(?:\{\{\w+\}\})+", cite, body)
    unused = set(REFS) - set(order)
    assert not unused, f"uncited references: {unused}"
    refs = '<ol class="refs">\n' + "\n".join(f"<li>{REFS[k]}.</li>" for k in order) + "\n</ol>"
    body = body.replace("{{REFERENCES}}", refs)
    abstract = body[body.index('<div class="abstract"'):body.index('<div class="kw">')]
    main_text = body[body.index('<h2 class="num">Introduction'):body.index("<h2>Data and code availability</h2>")]
    head = ('<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            "<title>Variant recurrence in literature-derived Mendelian cohorts</title>\n"
            f"<style>\n{(H / 'style.css').read_text()}\n</style>\n</head>\n")
    (H / "paper.html").write_text(head + body)
    print(f"{len(order)} references; abstract {words(abstract)} words; main text {words(main_text)} words")


if __name__ == "__main__":
    main()
