"""Shared paths and helpers for the multigene data pipeline."""
import os, json, time, hashlib, sys, types
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = Path(os.environ.get("PHENOBRIDGE_DATA", "/home/user/data_ext"))
PPKT_ROOT = DATA / "ppkt/monarch-initiative-phenopacket-store-3f36198/notebooks"
RAW = DATA / "raw"
CACHE = DATA / "cache"
DERIVED = DATA / "derived"
INTER = DATA / "intermediate"
AF_TAR = RAW / "alphafold/UP000005640_9606_HUMAN_v6.tar"
AF_DIR = RAW / "alphafold/models"
for p in (RAW, CACHE, DERIVED, INTER, AF_DIR):
    p.mkdir(parents=True, exist_ok=True)

MIN_INDIVIDUALS = 20


def import_mvp():
    """Import the repo's mechanistic_mvp with gradio stubbed (UI only)."""
    sys.modules.setdefault("gradio", types.ModuleType("gradio"))
    sys.path.insert(0, str(REPO))
    import mechanistic_mvp as M
    return M


def cached_request(kind, key, fetch, ext="json"):
    """Disk cache for REST calls: CACHE/kind/<sha1(key)>.json. fetch() -> python obj."""
    d = CACHE / kind
    d.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha1(key.encode()).hexdigest()
    f = d / f"{h}.{ext}"
    if f.exists():
        return json.loads(f.read_text())
    obj = fetch()
    if obj is not None:
        f.write_text(json.dumps(obj))
    return obj


def http_json(url, method="GET", payload=None, headers=None, tries=6):
    import requests
    hdr = {"Accept": "application/json", "Content-Type": "application/json"}
    hdr.update(headers or {})
    for i in range(tries):
        try:
            if method == "GET":
                r = requests.get(url, headers=hdr, timeout=120)
            else:
                r = requests.post(url, headers=hdr, data=json.dumps(payload), timeout=300)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503, 504):
                wait = float(r.headers.get("Retry-After", 2 ** i))
                time.sleep(wait)
                continue
            if r.status_code in (400, 404):
                return {"__error__": r.status_code, "text": r.text[:2000]}
            time.sleep(2 ** i)
        except Exception as e:  # network hiccup
            time.sleep(2 ** i)
    raise RuntimeError(f"failed {method} {url}")
