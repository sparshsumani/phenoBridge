"""Shared paths and helpers for the allelic-heterogeneity data pipeline (phenoBridge/allelic).

All scripts are run as `python3 -I aNN_*.py` (isolated mode: the script directory is not put on
sys.path automatically, so it is added here explicitly, together with ../multigene for reuse).

Public helpers for the analysis stage:
    ca_coords(uniprot)                -> (pos, xyz[L,3], plddt) of the AlphaFold model used
    contact_graph(uniprot)            -> (pos, W, A): column-stochastic W = A D^-1 of the 8 A C-alpha graph
    rwr_profile(uniprot, residue, alpha=0.3) -> pandas.Series pos -> stationary RWR probability
                                         (exact sparse solve p = alpha (I-(1-alpha)W)^-1 e_s, normalised to 1;
                                          identical to mechanistic_mvp.rwr / multigene.structfeat)
"""
import os, sys, json, hashlib, time, subprocess
from pathlib import Path
from functools import lru_cache

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
MULTIGENE = REPO / "multigene"
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

DATA = Path(os.environ.get("PHENOBRIDGE_DATA", "/home/user/data_ext"))
RAW = DATA / "raw"
CACHE = DATA / "cache"
DERIVED_PPKT = DATA / "derived"            # multigene / phenopacket-store outputs (input here)
INTER_PPKT = DATA / "intermediate"
OUT = DATA / "derived_allelic"
INTER = OUT / "intermediate"
RFDIR = OUT / "residue_features"
GRAPHDIR = OUT / "contact_graphs"
AF_TAR = RAW / "alphafold/UP000005640_9606_HUMAN_v6.tar"
AF_DIR = RAW / "alphafold/models"
AF_LISTING = INTER_PPKT / "af_tar_listing.txt"   # produced by multigene/s05 (tar -tf)
for p in (OUT, INTER, RFDIR, GRAPHDIR, AF_DIR):
    p.mkdir(parents=True, exist_ok=True)

# ---- external inputs (downloaded by a01_download.py) ----
URLS = {
    "clinvar": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
    "mondo": "https://purl.obolibrary.org/obo/mondo.obo",
    "hpoa": "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa",
    "mane_summary": "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz",
    "mane_protein": "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.refseq_protein.faa.gz",
    "uniprot_human": ("https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=%28organism_id%3A9606%29+AND+"
                      "%28reviewed%3Atrue%29&fields=accession%2Cid%2Cgene_primary%2Cgene_synonym%2Cxref_hgnc%2C"
                      "xref_mane-select%2Clength%2Csequence%2Csequence_version"),
    "alphamissense": "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz",
    "mim2gene": "https://omim.org/static/omim/data/mim2gene.txt",
}
FILES = {
    "clinvar": RAW / "clinvar/variant_summary.txt.gz",
    "mondo": RAW / "mondo/mondo.obo",
    "hpoa": RAW / "hpoa/phenotype.hpoa",
    "mane_summary": RAW / "mane/MANE.GRCh38.v1.5.summary.txt.gz",
    "mane_protein": RAW / "mane/MANE.GRCh38.v1.5.refseq_protein.faa.gz",
    "uniprot_human": RAW / "uniprot/uniprot_human_reviewed.tsv",
    "alphamissense": RAW / "alphamissense_allelic/am_subset_allelic.tsv.gz",  # filtered stream, see a09
    "mim2gene": RAW / "omim/mim2gene.txt",
}

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U"}
AA1TO3 = {v: k for k, v in AA3.items()}

MIN_VARIANTS_PER_DISEASE = 5
MAX_AF_SINGLE_FRAGMENT_LEN = 2700


def mg_import(name):
    """Import a module from ../multigene (which itself does `from common import ...` and must see
    multigene/common.py, not this file).  The multigene `common` is swapped into sys.modules only
    while the module is imported; afterwards the module keeps its own references."""
    import importlib, importlib.util
    if name in sys.modules:
        return sys.modules[name]
    mine = sys.modules.get("common")
    mg = sys.modules.get("_multigene_common")
    if mg is None:
        spec = importlib.util.spec_from_file_location("_multigene_common", MULTIGENE / "common.py")
        mg = importlib.util.module_from_spec(spec); spec.loader.exec_module(mg)
        sys.modules["_multigene_common"] = mg
    sys.modules["common"] = mg
    if str(MULTIGENE) not in sys.path:
        sys.path.append(str(MULTIGENE))
    try:
        mod = importlib.import_module(name)
    finally:
        if mine is not None:
            sys.modules["common"] = mine
        else:
            sys.modules.pop("common", None)
    return mod


def save_json(obj, path):
    Path(path).write_text(json.dumps(obj, indent=1, default=str))


def load_json(path):
    return json.loads(Path(path).read_text())


def af_model_path(acc):
    return AF_DIR / f"AF-{acc}-F1-model_v6.pdb.gz"


def read_fasta_gz(path):
    import gzip
    out, cur, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur:
                    out[cur] = "".join(buf)
                cur = line[1:].split()[0]; buf = []
            else:
                buf.append(line)
    if cur:
        out[cur] = "".join(buf)
    return out


# ---------------- analysis helpers ----------------
@lru_cache(maxsize=64)
def _residue_table(uniprot):
    import pandas as pd
    return pd.read_parquet(RFDIR / f"{uniprot}.parquet")


def ca_coords(uniprot):
    """C-alpha coordinates of the AlphaFold v6 model used for the features.
    Returns (pos[int, L], xyz[float, L x 3], plddt[float, L]) in residue order."""
    import numpy as np
    df = _residue_table(uniprot)
    return df.pos.to_numpy(int), df[["ca_x", "ca_y", "ca_z"]].to_numpy(float), df.plddt.to_numpy(float)


@lru_cache(maxsize=64)
def contact_graph(uniprot):
    """Rebuild the 8 A C-alpha contact graph (all residues, undirected, unweighted) from the saved
    edge list contact_graphs/<UNIPROT>.parquet (columns i, j = residue numbers, i < j, dist_A).
    Returns (pos, W, A) with W = A D^-1, D = max(1, degree) exactly as in mechanistic_mvp.rwr."""
    import numpy as np, pandas as pd, scipy.sparse as sp
    pos = _residue_table(uniprot).pos.to_numpy(int)
    e = pd.read_parquet(GRAPHDIR / f"{uniprot}.parquet")
    idx = {p: k for k, p in enumerate(pos)}
    r = e.i.map(idx).to_numpy(); c = e.j.map(idx).to_numpy()
    n = len(pos)
    A = sp.csr_matrix((np.ones(2 * len(r)), (np.r_[r, c], np.r_[c, r])), shape=(n, n))
    deg = np.maximum(1, np.asarray(A.sum(axis=0)).ravel())
    return pos, (A @ sp.diags(1.0 / deg)).tocsc(), A


@lru_cache(maxsize=16)
def _lu(uniprot, alpha):
    import scipy.sparse as sp
    from scipy.sparse.linalg import splu
    pos, W, _ = contact_graph(uniprot)
    return splu((sp.identity(len(pos), format="csc") - (1 - alpha) * W).tocsc())


def rwr_profile(uniprot, residue, alpha=0.3):
    """Full stationary RWR vector seeded at `residue` (AlphaFold/UniProt numbering), exact solve.
    Returns pandas.Series indexed by residue number, summing to 1."""
    import numpy as np, pandas as pd
    pos, W, _ = contact_graph(uniprot)
    k = np.where(pos == int(residue))[0]
    if not len(k):
        raise KeyError(f"residue {residue} not in model of {uniprot}")
    e = np.zeros(len(pos)); e[k[0]] = 1.0
    p = alpha * _lu(uniprot, float(alpha)).solve(e)
    return pd.Series(p / p.sum(), index=pos, name=f"rwr_{uniprot}_{residue}_a{alpha}")
