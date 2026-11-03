"""Shared paths and helpers for the temporal ClinVar "3D label propagation" VUS dataset (phenoBridge/vus3d).

Re-uses ../allelic/common.py (loaded as module `_allelic_common`; its constants and helpers are
re-exported here) and, through it, ../multigene.  All scripts run as `python3 -I vNN_*.py`.

Outputs go to $PHENOBRIDGE_DATA/derived_vus3d (default /home/user/data_ext/derived_vus3d).

Public helpers for the analysis stage (same semantics as allelic/common.py, but reading
derived_vus3d/residue_features and derived_vus3d/contact_graphs):
    ca_coords(uniprot)                -> (pos, xyz[L,3], plddt) of the AlphaFold v6 model used
    contact_graph(uniprot)            -> (pos, W, A): W = A D^-1 (D = max(1, degree)) of the 8 A C-alpha graph
    rwr_profile(uniprot, residue, alpha=0.3) -> pandas.Series pos -> stationary RWR probability
                                         (exact sparse solve p = alpha (I-(1-alpha)W)^-1 e_s, normalised to 1;
                                          identical to mechanistic_mvp.rwr / multigene.structfeat)
    rwr_matrix_solve(uniprot, b, alpha=0.3) -> alpha (I-(1-alpha)W)^-1 b for an arbitrary seed vector b
                                         (multi-seed propagation of labelled residues in one solve)
    residue_table(uniprot)            -> per-residue feature table (pandas.DataFrame)
"""
import os, sys, importlib, importlib.util
from pathlib import Path
from functools import lru_cache

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
ALLELIC = REPO / "allelic"
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))


def _load_allelic_common():
    if "_allelic_common" in sys.modules:
        return sys.modules["_allelic_common"]
    spec = importlib.util.spec_from_file_location("_allelic_common", ALLELIC / "common.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["_allelic_common"] = m
    spec.loader.exec_module(m)
    # allelic/common.py puts allelic/ at sys.path[0]; keep vus3d first so `import common` is ours
    while str(ALLELIC) in sys.path:
        sys.path.remove(str(ALLELIC))
    sys.path.append(str(ALLELIC))
    return m


AC = _load_allelic_common()
# re-exports
DATA, RAW, CACHE = AC.DATA, AC.RAW, AC.CACHE
AF_TAR, AF_DIR, AF_LISTING = AC.AF_TAR, AC.AF_DIR, AC.AF_LISTING
FILES, URLS, AA3, AA1TO3 = AC.FILES, AC.URLS, AC.AA3, AC.AA1TO3
MAX_AF_SINGLE_FRAGMENT_LEN = AC.MAX_AF_SINGLE_FRAGMENT_LEN
save_json, load_json, read_fasta_gz, af_model_path, mg_import = AC.save_json, AC.load_json, AC.read_fasta_gz, AC.af_model_path, AC.mg_import
ALLELIC_OUT, ALLELIC_RFDIR, ALLELIC_GRAPHDIR = AC.OUT, AC.RFDIR, AC.GRAPHDIR

OUT = DATA / "derived_vus3d"
INTER = OUT / "intermediate"
RFDIR = OUT / "residue_features"
GRAPHDIR = OUT / "contact_graphs"
for p in (OUT, INTER, RFDIR, GRAPHDIR):
    p.mkdir(parents=True, exist_ok=True)

# ---- ClinVar releases ----
CLINVAR_ARCHIVE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2023/variant_summary_2023-01.txt.gz"
RELEASES = {
    "2023": dict(file=RAW / "clinvar_2023_01/variant_summary_2023-01.txt.gz", url=CLINVAR_ARCHIVE_URL, label="2023-01"),
    "2026": dict(file=FILES["clinvar"], url=URLS["clinvar"], label="2026-11-01"),
}
CUTOFF_DATE = "2023-01-31"          # end of the train release month (release file dated 2023-01-05)
AM_SUBSET = RAW / "alphamissense_vus3d/am_subset_vus3d.tsv.gz"
AM_FULL_TMP = RAW / "alphamissense_vus3d/_full_tmp/AlphaMissense_aa_substitutions.tsv.gz"
MIN_TRAIN_RESIDUES = 5               # protein gets graph + residue features if >=1 TEST and >=5 TRAIN-labelled residues
THRESHOLDS = (5, 10, 20)


def al_import(name):
    """Import an allelic aNN module (which does `from common import ...` and must see allelic/common.py)."""
    full = "_allelic_" + name
    if full in sys.modules:
        return sys.modules[full]
    mine = sys.modules.get("common")
    sys.modules["common"] = AC
    try:
        spec = importlib.util.spec_from_file_location(full, ALLELIC / f"{name}.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[full] = mod
        spec.loader.exec_module(mod)
    finally:
        if mine is not None:
            sys.modules["common"] = mine
        else:
            sys.modules.pop("common", None)
    return mod


# ---------------- analysis helpers ----------------
@lru_cache(maxsize=128)
def residue_table(uniprot):
    import pandas as pd
    return pd.read_parquet(RFDIR / f"{uniprot}.parquet")


def ca_coords(uniprot):
    """C-alpha coordinates of the AlphaFold v6 model used for the features: (pos, xyz[L,3], plddt)."""
    df = residue_table(uniprot)
    return df.pos.to_numpy(int), df[["ca_x", "ca_y", "ca_z"]].to_numpy(float), df.plddt.to_numpy(float)


@lru_cache(maxsize=128)
def contact_graph(uniprot):
    """8 A C-alpha contact graph from contact_graphs/<UNIPROT>.parquet (i<j residue numbers, dist_A).
    Returns (pos, W, A) with W = A D^-1, D = max(1, degree) exactly as in mechanistic_mvp.rwr."""
    import numpy as np, pandas as pd, scipy.sparse as sp
    pos = residue_table(uniprot).pos.to_numpy(int)
    e = pd.read_parquet(GRAPHDIR / f"{uniprot}.parquet")
    idx = {p: k for k, p in enumerate(pos)}
    r = e.i.map(idx).to_numpy(); c = e.j.map(idx).to_numpy()
    n = len(pos)
    A = sp.csr_matrix((np.ones(2 * len(r)), (np.r_[r, c], np.r_[c, r])), shape=(n, n))
    deg = np.maximum(1, np.asarray(A.sum(axis=0)).ravel())
    return pos, (A @ sp.diags(1.0 / deg)).tocsc(), A


@lru_cache(maxsize=32)
def _lu(uniprot, alpha):
    import scipy.sparse as sp
    from scipy.sparse.linalg import splu
    pos, W, _ = contact_graph(uniprot)
    return splu((sp.identity(len(pos), format="csc") - (1 - alpha) * W).tocsc())


def rwr_matrix_solve(uniprot, b, alpha=0.3):
    """alpha * (I - (1-alpha) W)^-1 b  (b: length-L vector in residue order of residue_table)."""
    import numpy as np
    return alpha * _lu(uniprot, float(alpha)).solve(np.asarray(b, dtype=float))


def rwr_profile(uniprot, residue, alpha=0.3):
    """Full stationary RWR vector seeded at `residue` (AlphaFold/UniProt numbering), exact solve.
    Returns pandas.Series indexed by residue number, summing to 1."""
    import numpy as np, pandas as pd
    pos, W, _ = contact_graph(uniprot)
    k = np.where(pos == int(residue))[0]
    if not len(k):
        raise KeyError(f"residue {residue} not in model of {uniprot}")
    e = np.zeros(len(pos)); e[k[0]] = 1.0
    p = rwr_matrix_solve(uniprot, e, alpha)
    return pd.Series(p / p.sum(), index=pos, name=f"rwr_{uniprot}_{residue}_a{alpha}")
