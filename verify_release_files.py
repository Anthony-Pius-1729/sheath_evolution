#!/usr/bin/env python3
"""Check versioned project file hashes after downloading the release source archive."""
from pathlib import Path
import csv
import hashlib

ROOT = Path(__file__).resolve().parent
with (ROOT / 'RELEASE_FILES_SHA256.csv').open() as f:
    rows = list(csv.DictReader(f))
for row in rows:
    p = ROOT / row['path']
    if not p.is_file() or p.stat().st_size != int(row['bytes']):
        raise SystemExit('Missing file or size mismatch: ' + row['path'])
    if hashlib.sha256(p.read_bytes()).hexdigest() != row['sha256']:
        raise SystemExit('Checksum mismatch: ' + row['path'])
print(f'Verified {len(rows)} versioned project files; no scientific calculations performed.')
