#!/usr/bin/env python3
"""Bounded structural reproducibility check (reviewer concern 1).

Question: can the documented procedure, run against the publicly deposited
structures, reconstruct the retained residue-level structural measurements that
the regression consumes?

Design constraints, taken from the review:
  * The producer is used UNMODIFIED. Its sources are copied here and their
    SHA-256 checked against manuscript_latex/structural_provenance_audit.json
    before anything runs, so the check cannot silently test different code.
  * Output goes to this directory only. The retained
    sheath_pilot/qc_output/R_*.csv is opened read-only and never rewritten, so
    reconstructed values can never substitute for the historical ones.
  * Scope is bounded to the proteins the regression actually uses: AFP Afp2 and
    Afp3 (6RBN/6RC8) and diffocin (8V3X/8V3Y). The excluded systems contribute
    no regression value and are not re-run.

Outcome A: reconstruction reproduces the retained values -> record the recovered
identifiers and the agreement statistics.
Outcome B: it does not -> report the discrepancy. The retained measurements and
the existing disclosure stand; nothing is overwritten either way.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
ROOT = HERE.parent
RETAINED = ROOT / "sheath_pilot" / "qc_output"
AUDIT = ROOT / "manuscript_latex" / "structural_provenance_audit.json"

SYSTEMS = {
    "afp": ("Serratia entomophila AFP", "6RBN", "6RC8"),
    "diffocin": ("Clostridioides difficile diffocin", "8V3X", "8V3Y"),
}
TOL = 1e-6


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify_sources() -> dict:
    want = json.loads(AUDIT.read_text())["source_sha256"]
    got, bad = {}, []
    for rel, expect in want.items():
        if not rel.startswith("sheath_pilot/"):
            continue
        local = HERE / Path(rel).name
        if not local.exists():
            bad.append(f"{rel}: missing from the check directory")
            continue
        h = sha256(local)
        got[rel] = h
        if h != expect:
            bad.append(f"{rel}: sha256 {h[:16]} != archived {expect[:16]}")
    if bad:
        sys.exit("ABORT: producer sources do not match the archive:\n  " + "\n  ".join(bad))
    return got


def main() -> int:
    checksums = verify_sources()
    print(f"producer sources verified against the archive ({len(checksums)} files)")

    sys.path.insert(0, str(HERE))
    import structural_qc as sq  # noqa: E402  - imported after the checksum gate

    sq.STRUCT = HERE / "structures"
    sq.OUT = HERE / "qc_output"
    sq.OUT.mkdir(parents=True, exist_ok=True)
    assert sq.OUT.resolve() != RETAINED.resolve(), "output would overwrite retained data"

    report: list[str] = []
    rows = []
    for key, (label, pe, pc) in SYSTEMS.items():
        print(f"reconstructing {key} ({pe}/{pc}) ...", flush=True)
        df = sq.analyse_system(key, label, pe, pc, report)
        if df is None:
            print(f"  {key}: producer returned no frame")
            continue
        df.to_csv(sq.OUT / f"R_{key}.csv", index=False)
        rows.append((key, df))
    (sq.OUT / "reconstruction_report.txt").write_text("\n".join(report))

    cmp_rows = []
    for key, new in rows:
        old = pd.read_csv(RETAINED / f"R_{key}.csv")
        for fam in sorted(set(old.family) & set(new.family)):
            o = old[(old.family == fam) & old.in_domain].set_index("seqres_pos")
            n = new[(new.family == fam) & new.in_domain].set_index("seqres_pos")
            shared = sorted(set(o.index) & set(n.index))
            rec = {"system": key, "family": fam,
                   "sites_retained": len(o), "sites_reconstructed": len(n),
                   "sites_shared": len(shared),
                   "sites_only_retained": len(set(o.index) - set(n.index)),
                   "sites_only_reconstructed": len(set(n.index) - set(o.index))}
            for col in ("R", "rsa", "wcn", "gnm"):
                if col not in o.columns or col not in n.columns:
                    continue
                a = o.loc[shared, col].astype(float)
                b = n.loc[shared, col].astype(float)
                d = (a - b).abs()
                rec[f"{col}_max_abs_diff"] = float(d.max()) if len(d) else None
                rec[f"{col}_n_mismatched"] = int((d > TOL).sum()) if len(d) else None
                rec[f"{col}_spearman"] = (float(a.corr(b, method="spearman"))
                                          if len(d) > 2 else None)
            cmp_rows.append(rec)
    cmp = pd.DataFrame(cmp_rows)
    cmp.to_csv(HERE / "reconstruction_comparison.csv", index=False)
    print()
    print(cmp.to_string(index=False))

    cols = [c for c in ("R", "rsa", "wcn", "gnm") if f"{c}_n_mismatched" in cmp.columns]
    exact = bool(len(cmp)) and all(
        (cmp[f"{c}_n_mismatched"].fillna(0) == 0).all() for c in cols)
    same_sites = bool(len(cmp)) and (cmp.sites_only_retained.sum() == 0
                                     and cmp.sites_only_reconstructed.sum() == 0)
    verdict = ("A: reconstruction reproduces the retained measurements"
               if exact and same_sites else
               "B: reconstruction does NOT reproduce the retained measurements")
    json.dump({"verdict": verdict, "exact_value_agreement": exact,
               "identical_site_sets": same_sites,
               "scope": "AFP (6RBN/6RC8) and diffocin (8V3X/8V3Y) - the three "
                        "proteins used in the regression",
               "tolerance": TOL, "producer_sha256": checksums,
               "retained_never_written": True, "comparison": cmp_rows},
              open(HERE / "reconstruction_verdict.json", "w"), indent=1)
    print(f"\nVERDICT {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
