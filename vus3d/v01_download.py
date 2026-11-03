#!/usr/bin/env python3
"""Step 1: download (once) the two ClinVar variant_summary releases and the full AlphaMissense table,
and record provenance (exact file name, URL, HTTP Last-Modified / Content-Length, local size, sha256).

  * ClinVar 2023-01 (train snapshot): archive/2023/variant_summary_2023-01.txt.gz (the archive
    directory lists monthly files per year; 2023-01 is the release of 2023-01-05)
  * ClinVar current (2026-11-01, test snapshot): raw/clinvar/variant_summary.txt.gz, already downloaded by
    allelic/a01 (re-downloaded only if absent)
  * AlphaMissense_aa_substitutions.tsv.gz: streamed with resumable curl (-C -, retries; the proxy can cut
    long transfers) into raw/alphamissense_vus3d/_full_tmp/; v06 filters it and deletes the full file.
Output intermediate/downloads.json
"""
import hashlib, json, subprocess, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RELEASES, INTER, URLS, AM_FULL_TMP, AM_SUBSET, save_json, load_json


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


def fetch(url, out, tries=40):
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(out) + ".part")
    for i in range(tries):
        r = subprocess.run(["curl", "-sS", "--http1.1", "-C", "-", "-o", str(tmp), url])
        if r.returncode == 0:
            break
        print("download interrupted (curl exit", r.returncode, ") -> resume", flush=True)
        time.sleep(3)
    else:
        raise SystemExit(f"download failed: {url}")
    subprocess.run(["gzip", "-t", str(tmp)], check=True)
    tmp.rename(out)


def main():
    pf = INTER / "downloads.json"
    prov = load_json(pf) if pf.exists() else {}
    for rel, d in RELEASES.items():
        key = f"clinvar_{rel}"
        out = d["file"]
        if not out.exists():
            print("downloading", d["url"], flush=True)
            fetch(d["url"], out)
        rec = prov.get(key) or {}
        if rec.get("size") != out.stat().st_size:
            rec = dict(url=d["url"], file=str(out), file_name=out.name, release=d["label"], http=head(d["url"]),
                       http_checked_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                       size=out.stat().st_size, sha256=sha256(out))
            cl = (rec["http"] or {}).get("content-length")
            rec["server_size_equals_local"] = (int(cl) == rec["size"]) if cl and cl.isdigit() else None
        prov[key] = rec
        print(key, json.dumps({k: v for k, v in rec.items() if k != "sha256"}), flush=True)
    # AlphaMissense full table (only while the filtered subset does not exist yet)
    rec = prov.get("alphamissense") or dict(url=URLS["alphamissense"])
    if not AM_SUBSET.exists():
        if not AM_FULL_TMP.exists():
            print("downloading AlphaMissense", flush=True)
            fetch(URLS["alphamissense"], AM_FULL_TMP)
        rec.update(http=head(URLS["alphamissense"]), full_file=str(AM_FULL_TMP), full_size=AM_FULL_TMP.stat().st_size)
    prov["alphamissense"] = rec
    save_json(prov, pf)


if __name__ == "__main__":
    main()
