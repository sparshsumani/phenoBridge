"""Shared paths and cached REST helpers for the KBG follow-up analysis."""
import hashlib
import json
import sys
import time
import types
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
CACHE = HERE / "cache"
OUT = HERE / "out"
FIG = HERE / "fig"
for d in (CACHE, OUT, FIG):
    d.mkdir(exist_ok=True)

COHORT_JSON = REPO / "ANKRD11_KBGS_individuals.json"
PDB = REPO / "AF-Q6UB99-F1-model_v6.pdb"
PPKT_DIR = Path("/home/user/data_ext/ppkt/monarch-initiative-phenopacket-store-3f36198/"
                "notebooks/ANKRD11/phenopackets")
MANE_ENST = "ENST00000301030"   # MANE Select, = NM_013275.6
CDS_LEN = 7992                  # 2663 aa * 3 + stop
PROT_LEN = 2663


def import_mvp():
    """Import the repo's mechanistic_mvp with gradio stubbed out."""
    sys.modules.setdefault("gradio", types.ModuleType("gradio"))
    sys.path.insert(0, str(REPO))
    import mechanistic_mvp as M  # noqa
    return M


def _key(*parts):
    return hashlib.sha1(json.dumps(parts, sort_keys=True).encode()).hexdigest()[:16]


def cached_get(url, params=None, name=None):
    f = CACHE / f"{name or 'get'}_{_key(url, params)}.json"
    if f.exists():
        return json.loads(f.read_text())
    for attempt in range(5):
        r = requests.get(url, params=params, headers={"Accept": "application/json"}, timeout=120)
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(2 + 3 * attempt); continue
        r.raise_for_status()
        f.write_text(r.text)
        return r.json()
    r.raise_for_status()


def cached_post(url, payload, params=None, name=None):
    f = CACHE / f"{name or 'post'}_{_key(url, payload, params)}.json"
    if f.exists():
        return json.loads(f.read_text())
    for attempt in range(5):
        r = requests.post(url, json=payload, params=params,
                          headers={"Accept": "application/json", "Content-Type": "application/json"},
                          timeout=300)
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(2 + 3 * attempt); continue
        r.raise_for_status()
        f.write_text(r.text)
        return r.json()
    r.raise_for_status()
