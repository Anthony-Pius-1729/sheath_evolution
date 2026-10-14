#!/usr/bin/env python3
"""Record hashes after intentional package edits; does not perform reconstruction."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    excluded = {'package_manifest.json', 'verification_report.json'}
    files = {p for p in HERE.rglob('*') if p.is_file() and p.name not in excluded
             and '__pycache__' not in p.parts and 'reruns' not in p.parts}
    audit = json.loads((ROOT / 'structural_provenance_audit.json').read_text())
    files.update(ROOT / p for p in audit['source_sha256'])
    files.update(ROOT / p for p in ['structural_provenance_audit.json', 'FROZEN_ANALYSIS.md',
        'sheath_pilot/qc_output/R_afp.csv', 'sheath_pilot/qc_output/R_diffocin.csv',
        'results/hypothesis_test_primary.json', 'results/sensitivity.json'])
    manifest = dict(scope='Reconstruction evidence and required retained inputs; refreshed only after deliberate packaging edits.',
        original_frozen_sha256=hashlib.sha256((HERE / 'FROZEN_ANALYSIS_original.md').read_bytes()).hexdigest(),
        files={str(p.relative_to(ROOT)): dict(bytes=p.stat().st_size,
            sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(files)})
    (HERE / 'package_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f"Recorded {len(files)} package/input hashes.")


if __name__ == '__main__':
    main()
