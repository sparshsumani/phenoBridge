#!/usr/bin/env python3
"""Fetch single AlphaFold PDB members from the 27 GB swissprot_pdb_v6.tar on the EBI FTP with
HTTP Range requests (binary search over the alphabetically ordered tar), for accessions that the
human proteome tar does not contain (e.g. GNAS P63092: the proteome set models only Q5JWF2).
Usage: python3 fetch_af_from_swissprot_tar.py P63092 [ACC ...]   -> RAW/alphafold/models/
"""
import re, sys, requests
from common import AF_DIR

URL = "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/swissprot_pdb_v6.tar"
CHUNK = 1 << 20
HDR = re.compile(rb"(AF-([A-Z0-9]+)-F(\d+)-model_v6\.pdb\.gz)\x00")


def get(a, b):
    r = requests.get(URL, headers={"Range": f"bytes={a}-{b}"}, timeout=120)
    assert r.status_code == 206, r.status_code
    return r.content


def headers_in(off):
    """[(abs_offset, name, acc, frag, size)] for tar headers fully inside the chunk at off."""
    buf = get(off, off + CHUNK - 1)
    out = []
    for m in HDR.finditer(buf):
        s = m.start()
        if (off + s) % 512 or s + 512 > len(buf) or buf[s + 257:s + 262] != b"ustar":
            continue
        size = int(buf[s + 124:s + 136].strip(b"\x00 ").decode(), 8)
        out.append((off + s, m.group(1).decode(), m.group(2).decode(), int(m.group(3)), size))
    return out


def total_size():
    return int(requests.head(URL, timeout=60).headers["Content-Length"])


def fetch(acc, frag=1):
    target = (acc, frag)
    lo, hi = 0, total_size() - CHUNK
    while hi - lo > CHUNK:
        mid = (lo + hi) // 2 // 512 * 512
        hs = headers_in(mid)
        if not hs:
            mid += CHUNK; hs = headers_in(mid)
        for h in hs:
            if (h[2], h[3]) == target:
                return save(h)
        first = (hs[0][2], hs[0][3])
        if target < first:
            hi = mid
        else:
            lo = mid
    for off in range(lo, hi + 2 * CHUNK, CHUNK // 2):
        for h in headers_in(off):
            if (h[2], h[3]) == target:
                return save(h)
    raise SystemExit(f"{acc} not found")


def save(h):
    off, name, acc, frag, size = h
    data = get(off + 512, off + 512 + size - 1)
    (AF_DIR / name).write_bytes(data)
    print("saved", AF_DIR / name, size, "bytes")
    return AF_DIR / name


if __name__ == "__main__":
    for a in sys.argv[1:]:
        if not (AF_DIR / f"AF-{a}-F1-model_v6.pdb.gz").exists():
            fetch(a)
