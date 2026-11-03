"""Vectorised, numerically equivalent re-implementation of Bio.PDB.HSExposureCA / HSExposureCB
(Hamelryck 2005) as used in allelic/a09 (radius 13 A, offset 0).  Biopython's version is a pure-Python
O(L^2) double loop that takes ~75 % of a09's per-protein time; this one gives the same counts
(verified against the Biopython values stored in derived_allelic/residue_features, see v07 --verify-hse).

Biopython logic reproduced:
  * polypeptides: CaPPBuilder (standard amino acids with CA; consecutive residues connected when CA-CA < 4.3 A)
  * HSE-CA: pseudo-CB direction b = norm(norm(ca_i - ca_{i-1}) + norm(ca_i - ca_{i+1})); undefined (NaN) for the
    first/last residue of each polypeptide
  * HSE-CB: real CA->CB vector; Gly: pseudo CB = N rotated by -120 deg about the CA->C axis (Bio rotaxis);
    NaN when atoms are missing
  * every other CA (all polypeptides, self excluded) with |d| < radius is counted 'up' when
    arccos(clip(d.b / (|d||b|), -1, 1)) < pi/2, else 'down'
"""
import numpy as np
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.vectors import rotaxis


def _residues(model):
    """Residues of all CaPPBuilder polypeptides (in order) and peptide boundaries."""
    from Bio.PDB.Polypeptide import CaPPBuilder
    ppl = CaPPBuilder().build_peptides(model)
    res, start, end = [], [], []
    for pp in ppl:
        for i, r in enumerate(pp):
            res.append(r); start.append(i == 0); end.append(i == len(pp) - 1)
    return res, np.array(start), np.array(end)


def _counts(ca, b, valid, radius):
    up = np.full(len(ca), np.nan); dn = np.full(len(ca), np.nan)
    for i in np.where(valid)[0]:
        d = ca - ca[i]
        nrm = np.sqrt((d * d).sum(1))
        m = nrm < radius
        m[i] = False
        dd = d[m]; nd = nrm[m]
        c = (dd @ b[i]) / (nd * np.sqrt(b[i] @ b[i]))
        ang = np.arccos(np.clip(c, -1, 1))
        u = int((ang < np.pi / 2).sum())
        up[i] = u; dn[i] = int(m.sum()) - u
    return up, dn


def hse(model, radius=13.0):
    """Returns {resnum: (ca_up, ca_down, cb_up, cb_down)} (NaN where Biopython sets no value)."""
    res, start, end = _residues(model)
    ca = np.array([r["CA"].get_vector().get_array() for r in res], dtype=float)
    n = len(res)
    # HSE-CA
    b_ca = np.zeros((n, 3)); v_ca = ~(start | end)
    d1 = ca[1:-1] - ca[:-2]; d3 = ca[1:-1] - ca[2:]
    d1 /= np.linalg.norm(d1, axis=1)[:, None]; d3 /= np.linalg.norm(d3, axis=1)[:, None]
    bb = d1 + d3; bb /= np.linalg.norm(bb, axis=1)[:, None]
    b_ca[1:-1] = bb
    # HSE-CB
    b_cb = np.zeros((n, 3)); v_cb = np.zeros(n, bool)
    for k, r in enumerate(res):
        if r.get_resname() == "GLY":
            try:
                n_v = r["N"].get_vector(); c_v = r["C"].get_vector(); ca_v = r["CA"].get_vector()
            except Exception:
                continue
            n_v = n_v - ca_v; c_v = c_v - ca_v
            rot = rotaxis(-np.pi * 120.0 / 180.0, c_v)
            b_cb[k] = n_v.left_multiply(rot).get_array(); v_cb[k] = True
        elif r.has_id("CB") and r.has_id("CA"):
            b_cb[k] = (r["CB"].get_vector() - r["CA"].get_vector()).get_array(); v_cb[k] = True
    au, ad = _counts(ca, b_ca, v_ca, radius)
    cu, cd = _counts(ca, b_cb, v_cb, radius)
    return {int(r.id[1]): (au[k], ad[k], cu[k], cd[k]) for k, r in enumerate(res)}
