#!/usr/bin/env python3
"""Audit saved residue measurements and joins without recomputing structures or fits."""
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent
STRUCTURES = [
    ('r2pyocin', 'F1', 'R2 pyocin', 257, 143),
    ('afp', 'F1', 'AFP Afp2', 184, 84),
    ('afp', 'F2', 'AFP Afp3', 190, 82),
    ('salcis', 'F1', 'SalCIS', 168, 68),
    ('diffocin', 'F1', 'Diffocin', 128, 58),
]
PRIMARY = [
    ('afp_Afp2', 'AFP Afp2', 'orthology/afp_sf1_afp2_like.afa', 84, 27, 51),
    ('afp_Afp3', 'AFP Afp3', 'orthology/afp_sf1_afp3_like.afa', 81, 31, 64),
    ('diffocin', 'Diffocin', 'orthology/diffocin_primary.afa', 58, 33, 50),
]


def rows(path):
    with (ROOT / path).open() as f:
        return list(csv.DictReader(f))


def equal(a, b):
    if a == b:
        return True
    try:
        return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12)
    except ValueError:
        return False


def audit():
    inputs = {'supplementary_data/TableS1_residue_mapping.csv',
              'results/hypothesis_test_primary.json'}
    mapping = rows('supplementary_data/TableS1_residue_mapping.csv')
    assert len(mapping) == 927
    structural = []
    by_protein = {}
    for system, family, protein, count, eligible in STRUCTURES:
        path = f'sheath_pilot/qc_output/R_{system}.csv'
        inputs.add(path)
        original = [r for r in rows(path) if r['family'] == family]
        exported = [r for r in mapping if r['protein'] == protein]
        assert len(original) == len(exported) == count, protein
        assert [{k: v for k, v in r.items() if k != 'protein'} for r in exported] == original
        lookup = {r['canon']: r for r in original}
        assert len(lookup) == count, protein
        by_protein[protein] = lookup
        for r in original:
            union, inter, ext, con = (int(r[k]) for k in ('n_union', 'n_inter', 'n_ext', 'n_con'))
            assert union == ext + con - inter and 0 <= inter <= min(ext, con)
            assert (r['in_domain'] == 'True') == (union >= 3)
            assert int(r['seqres_pos']) > 0
            if union >= 3:
                assert equal(r['R'], str(1 - inter / union))
                for k in ('rsa', 'wcn', 'gnm'):
                    assert math.isfinite(float(r[k]))
            else:
                assert r['R'] == ''
        assert sum(r['in_domain'] == 'True' for r in original) == eligible
        structural.append(dict(protein=protein, recorded_positions=count, eligible_positions=eligible))

    models = {r['dataset']: r for r in json.loads((ROOT / 'results/hypothesis_test_primary.json').read_text())}
    fitted = []
    for dataset, protein, alignment, count, nseq, variable in PRIMARY:
        path = f'results/sitetable_primary_{dataset}.csv'
        rate_path = f'rates/{dataset}.res'
        inputs.update([path, alignment, rate_path])
        table = rows(path)
        assert len(table) == count == models[dataset]['sites_fitted']
        seqs = {}
        name = None
        for line in (ROOT / alignment).read_text().splitlines():
            if line.startswith('>'):
                name = line[1:].split()[0]
                assert name not in seqs
                seqs[name] = ''
            elif name:
                seqs[name] += line.strip()
        assert len(seqs) == nseq and len({len(s) for s in seqs.values()}) == 1
        positions = {pos: col for pos, col in enumerate(
            (col for col, aa in enumerate(seqs['REF']) if aa != '-'), 1)}
        rate_rows = {}
        for line in (ROOT / rate_path).read_text().splitlines():
            fields = line.split()
            if fields and fields[0].isdigit():
                assert fields[0] not in rate_rows
                rate_rows[fields[0]] = fields
        seen = set()
        varying = 0
        for r in table:
            assert r['canon'] not in seen
            seen.add(r['canon'])
            original = by_protein[protein][r['canon']]
            assert original['in_domain'] == 'True'
            assert all(equal(v, r[k]) for k, v in original.items()), (protein, r['canon'])
            pos = int(r['seqres_pos'])
            column = [s[positions[pos]] for s in seqs.values()]
            non_gap = sum(aa != '-' for aa in column)
            assert int(r['n_seq']) == nseq and int(r['n_aa']) == non_gap
            assert equal(r['occupancy'], str(non_gap / nseq)) and float(r['occupancy']) >= .8
            assert r['aa'] == seqs['REF'][positions[pos]]
            saved_rate = rate_rows[str(pos)]
            assert saved_rate[1] == r['aa'] and equal(saved_rate[2], r['rate'])
            # The interval has no internal whitespace in the retained Rate4Site output.
            assert equal(saved_rate[4], r['rate_std'])
            assert saved_rate[5] == f'{non_gap}/{nseq}'
            varying += len(set(column) - {'-'}) >= 2
        assert varying == variable
        fitted.append(dict(dataset=dataset, fitted_positions=count, sequences=nseq, variable_positions=variable))

    report = dict(
        scope='Saved-file integrity and correspondence audit; no coordinates, contacts, rates or models recomputed.',
        status='passed', structural=structural, recorded_positions=927,
        eligible_positions=435, primary=fitted, primary_fitted_positions=223,
        verified=['complete retained mapping export', 'contact-count arithmetic and eligibility',
                  'remodeling values from retained counts', 'structural-to-fitted site joins',
                  'alignment reference, occupancy and observed variability', 'saved Rate4Site scores and posterior SDs'],
        not_verified=['original atomic-coordinate bytes', 'historical selected subunit keys',
                      'individual contact partner identities', 'spatial dependence coverage'],
        source_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sorted(inputs)},
    )
    (ROOT / 'RETAINED_DATA_VALIDATION.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    audit()
