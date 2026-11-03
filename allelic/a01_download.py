#!/usr/bin/env python3
"""Step 1a: download (once) every external flat file used by the allelic pipeline and record
its version (HTTP Last-Modified / ETag, size, sha256, in-file version header where present).

Files go to $DATA/raw/<source>/ ; a provenance record is written to
$DATA/derived_allelic/intermediate/downloads.json.  Existing files are not re-downloaded
(delete a file to refresh it).  AlphaMissense is streamed + filtered later (a09), not stored whole.
"""
import hashlib, json, re, subprocess, sys, gzip
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import URLS, FILES, INTER, save_json


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def head(url):
    r = subprocess.run(["curl", "-sSIL", "--max-time", "60", url], capture_output=True, text=True)
    hd = {}
    for line in r.stdout.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            hd[k.strip().lower()] = v.strip()
    return {k: hd.get(k) for k in ("last-modified", "etag", "content-length")}


def main():
    prov_f = INTER / "downloads.json"
    prov = json.loads(prov_f.read_text()) if prov_f.exists() else {}
    for key in ("clinvar", "mondo", "hpoa", "mane_summary", "mane_protein", "uniprot_human", "mim2gene"):
        url, out = URLS[key], FILES[key]
        out.parent.mkdir(parents=True, exist_ok=True)
        if not out.exists():
            hd = head(url) if key != "uniprot_human" else {}
            print("downloading", key, url, flush=True)
            subprocess.run(["curl", "-sSL", "--http1.1", "--retry", "5", "--fail", "-o", str(out) + ".tmp", url], check=True)
            Path(str(out) + ".tmp").rename(out)
            rec = dict(url=url, file=str(out), http=hd)
        else:
            rec = prov.get(key, dict(url=url, file=str(out), http=None))
        rec["size"] = out.stat().st_size
        rec["sha256"] = sha256(out)
        # in-file version strings
        if key == "mondo":
            with open(out) as fh:
                for line in fh:
                    if line.startswith("data-version:"):
                        rec["version"] = line.split(":", 1)[1].strip(); break
                    if line.startswith("[Term]"):
                        break
        elif key == "hpoa":
            with open(out) as fh:
                for line in fh:
                    if not line.startswith("#"):
                        break
                    if "version" in line.lower() or "date" in line.lower():
                        rec.setdefault("header", []).append(line.strip())
        elif key == "mim2gene":
            with open(out) as fh:
                rec["version"] = [l.strip() for l in fh.readlines()[:3] if "Generated" in l]
        elif key == "clinvar":
            rec["version"] = "variant_summary.txt.gz, Last-Modified " + str((rec.get("http") or {}).get("last-modified"))
        prov[key] = rec
        print(key, json.dumps({k: v for k, v in rec.items() if k != "sha256"}), flush=True)
    save_json(prov, prov_f)


if __name__ == "__main__":
    main()
