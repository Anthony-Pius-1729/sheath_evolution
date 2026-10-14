#!/usr/bin/env python3
"""Optional bounded reconstruction from packaged coordinates; never fit rate models."""
from __future__ import annotations
import argparse
import json
import platform
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gnm-mode', choices=['explicit', 'historical'], default='explicit',
                        help='Explicit zero-mode exclusion (default), or original pseudoinverse.')
    args = parser.parse_args()
    from verify_package import verify_manifest
    verify_manifest()
    import numpy as np
    import pandas as pd
    import Bio
    import structural_qc as sq
    if args.gnm_mode == 'explicit':
        from covariates_zeromode_fix import monomer_gnm_fluctuation
        sq.cov.monomer_gnm_fluctuation = monomer_gnm_fluctuation

    runs = HERE / 'reruns'
    runs.mkdir(exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='reconstruction_', dir=runs))
    sq.STRUCT = HERE / 'structures'
    sq.OUT = output
    assert sq.OUT.resolve() != (ROOT / 'sheath_pilot/qc_output').resolve()
    selection_calls = []
    pick = sq.pick_reference

    def record_pick(subunits, family):
        sub = pick(subunits, family)
        if sub is not None:
            selection_calls.append(dict(family=family, key=sub.key, ring=int(sub.ring), azimuth=int(sub.azim)))
        return sub

    sq.pick_reference = record_pick
    report, comparisons, references = [], [], []
    for system, label, extended, contracted in [
        ('afp', 'Serratia entomophila AFP', '6RBN', '6RC8'),
        ('diffocin', 'Clostridioides difficile diffocin', '8V3X', '8V3Y'),
    ]:
        selection_calls.clear()
        new = sq.analyse_system(system, label, extended, contracted, report)
        if new is None:
            raise RuntimeError(f'No reconstructed table for {system}')
        old = pd.read_csv(ROOT / f'sheath_pilot/qc_output/R_{system}.csv')
        for family in sorted(old.family.unique()):
            a = old[(old.family == family) & old.in_domain].set_index('seqres_pos')
            b = new[(new.family == family) & new.in_domain].set_index('seqres_pos')
            same = set(a.index) == set(b.index)
            shared = sorted(set(a.index) & set(b.index))
            delta = {col: float(np.max(np.abs(a.loc[shared, col] - b.loc[shared, col])))
                     for col in ['R', 'rsa', 'wcn', 'gnm']}
            comparisons.append(dict(system=system, family=family, sites=len(a),
                identical_eligible_sets=same, max_absolute_difference=delta,
                agreement_within_1e_10=same and all(v <= 1e-10 for v in delta.values())))
            final_calls = [v for v in selection_calls if v['family'] == family][-2:]
            assert len(final_calls) == 2
            references.append(dict(system=system, family=family, extended=final_calls[0], contracted=final_calls[1]))
    (output / 'producer_report.txt').write_text(chr(10).join(report))
    result = dict(gnm_mode=args.gnm_mode, output=str(output), comparisons=comparisons,
        reconstructed_references=references, historical_execution_history_recovered=False,
        retained_files_written=False, python=platform.python_version(),
        numpy=np.__version__, pandas=pd.__version__, biopython=Bio.__version__,
        passed=all(r['agreement_within_1e_10'] for r in comparisons))
    (output / 'run_report.json').write_text(json.dumps(result, indent=2) + chr(10))
    print(json.dumps(result, indent=2))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
