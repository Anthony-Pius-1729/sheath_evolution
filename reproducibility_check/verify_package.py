#!/usr/bin/env python3
"""Verify supplied reconstruction evidence without rerunning structural calculations."""
from pathlib import Path
import csv
import hashlib
import json
import math
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FAMILIES = [('afp', 'F1', 'AFP Afp2', 84), ('afp', 'F2', 'AFP Afp3', 82),
            ('diffocin', 'F1', 'Diffocin', 58)]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify_manifest():
    manifest = json.loads((HERE / 'package_manifest.json').read_text())
    for rel, expected in manifest['files'].items():
        p = ROOT / rel
        assert p.is_file() and digest(p) == expected['sha256'], rel
        assert p.stat().st_size == expected['bytes'], rel
    return manifest


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    manifest = verify_manifest()
    audit = json.loads((ROOT / 'structural_provenance_audit.json').read_text())
    for rel, h in audit['source_sha256'].items():
        assert digest(ROOT / rel) == h
        if rel.startswith('sheath_pilot/'):
            assert digest(HERE / Path(rel).name) == h
    original = (HERE / 'FROZEN_ANALYSIS_original.md').read_bytes()
    assert digest(HERE / 'FROZEN_ANALYSIS_original.md') == manifest['original_frozen_sha256']
    assert (ROOT / 'FROZEN_ANALYSIS.md').read_bytes().startswith(original)
    comparisons = []
    initial_summary = {(r['system'], r['family']): r for r in read(HERE / 'reconstruction_comparison.csv')}
    for system, family, protein, count in FAMILIES:
        def eligible(path):
            records = [r for r in read(path) if r['family'] == family and r['in_domain'] == 'True']
            values = {r['seqres_pos']: r for r in records}
            assert len(values) == len(records)
            return values
        old = eligible(ROOT / f'sheath_pilot/qc_output/R_{system}.csv')
        raw = eligible(HERE / f'qc_output/R_{system}.csv')
        fixed = eligible(HERE / f'qc_output_zeromode/R_{system}_{family}.csv')
        assert len(old) == count and set(old) == set(raw) == set(fixed)
        def differences(other):
            return {col: max(abs(float(old[p][col]) - float(other[p][col])) for p in old)
                    for col in ['R', 'rsa', 'wcn', 'gnm']}
        saved, corrected = differences(raw), differences(fixed)
        for p in old:
            for col in ['canon', 'n_union', 'n_inter', 'n_ext', 'n_con', 'resid_ext', 'resid_con']:
                assert old[p][col] == raw[p][col] == fixed[p][col]
        assert all(math.isfinite(v) and v < (1.7e-13 if k == 'gnm' else 1e-12)
                   for k, v in saved.items())
        assert all(math.isfinite(v) and v < (1.7e-13 if k == 'gnm' else 1e-12)
                   for k, v in corrected.items())
        first = initial_summary[(system, family)]
        assert int(first['gnm_n_mismatched']) == count and float(first['gnm_max_abs_diff']) > 1e6
        assert int(first['sites_shared']) == count and int(first['sites_only_retained']) == int(first['sites_only_reconstructed']) == 0
        assert all(float(first[f'{k}_max_abs_diff']) < 1e-12 for k in ['R', 'rsa', 'wcn'])
        comparisons.append(dict(protein=protein, eligible_sites=count,
            identical_eligible_sets=True, saved_reconstruction_max_absolute_difference=saved,
            initial_GNM_discrepancy_reported_in_summary=float(first['gnm_max_abs_diff']),
            explicit_zero_mode_max_absolute_difference=corrected))
    reference = json.loads((HERE / 'reconstructed_references.json').read_text())
    assert len(reference['proteins']) == 3
    verbose = json.loads((HERE / 'recovered_identifiers.json').read_text())
    for protein in reference['proteins']:
        selections = set()
        for record in verbose[protein['system']]['reference_subunits_selected']:
            m = re.match(r"Subunit\(key='([^']+)', family='([^']+)', ring=(-?\d+), azim=(-?\d+),", record)
            assert m
            selections.add((m[1], m[2], int(m[3]), int(m[4])))
        for state in ['extended', 'contracted']:
            selected = protein[state]
            assert (selected['key'], protein['family'], selected['ring'], selected['azimuth']) in selections
    coordinates = json.loads((HERE / 'coordinate_inputs.json').read_text())['inputs']
    assert {r['accession'] for r in coordinates} == {'6RBN', '6RC8', '8V3X', '8V3Y'}
    for coordinate in coordinates:
        assert digest(ROOT / coordinate['path']) == coordinate['sha256']
    result = dict(status='passed', scope='Saved reconstruction evidence for the three fitted proteins only; no structural calculations rerun.',
        eligible_sites=224, fitted_sites=223, comparisons=comparisons,
        original_producer_hashes_verified=True, original_frozen_specification_preserved=True,
        coordinate_files_verified=4, historical_coordinate_byte_identity_verified=False,
        historical_subunit_execution_history_recovered=False, package_files_checked=len(manifest['files']))
    (HERE / 'verification_report.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
