#!/usr/bin/env python3
"""Regenerate the categorical pilot's verdict artifacts under the NATIVE rule.

WHY THIS SCRIPT EXISTS
----------------------
`pilot.py` counted every mapped residue, including residues whose interface
classification rests on contacts to symmetry-generated neighbours. For the
*Vibrio* extended state (5MXN, a single deposited ring) every inter-ring contact
must be generated from the deposited helical operator, so that run reported
VipB as PASS and "1 / 3 families passing" -- a property of the symmetry
operator, not of the deposition.

The governing rule for this project is NATIVE COORDINATES ONLY, under which the
pilot fails in 0 of 3 families. That verdict existed only in README prose; the
machine-readable artifacts still said PASS. This script regenerates
`output/summary.csv` and `output/qc_report.txt` under the native rule so the
record is self-consistent, and leaves the original run intact in
`output/archive_symmetry_inclusive_2026-10-10/`.

WHAT IS AND IS NOT RECOMPUTED
-----------------------------
Recomputed here, from the residue tables alone:
  - per-class residue counts, smallest interface class, per-family verdict
  - the post-hoc three-way collapse (still post hoc, still not the verdict)

Carried over verbatim from the archived run (structural metadata, independent of
the counting rule): cross-state identity, modelled/mapped residue counts,
substitutions, unmatched internal residues, ring counts, lattice failures.

NOT regenerated, because it requires the coordinates: the per-state structural
diagnostics, the contact-cutoff sweep and the per-class covariate means. Those
remain in the archived report and the regenerated report says where they stand.

SCOPE
-----
This is the historical CATEGORICAL pilot (four-way persistence classes,
pre-registered minimum class size 20). It has nothing to do with the frozen
continuous remodeling metric R_i in `qc_output/` (see `QC_RESULTS.md`). This
script touches `output/` only.

NO SCIENTIFIC RULE IS CHANGED HERE. The pre-registered threshold, the contact
definition, the class definitions, the engineering shell and the native-coordinate
rule are all as previously specified; only their consistent application is fixed.

Usage:  python3 regenerate_native_summary.py [--dry-run]
Needs pandas only.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCRIPT_VERSION = "regenerate_native_summary.py v1.0 (2026-10-12)"

HERE = Path(__file__).parent
OUT_DIR = HERE / "output"
ARCHIVE = OUT_DIR / "archive_symmetry_inclusive_2026-10-10"

# ---- rule constants, mirrored from pilot.py (NOT redefined here) ------------
CONTACT_CUTOFF = 5.0          # A, heavy atom -- pre-registered primary
MIN_CLASS_SIZE = 20           # pre-registered, never revised after seeing counts
ENGINEERING_SHELL = 8.0       # A, exclusion shell around internal indels
INTERFACE_CLASSES = ["persistent_same_partner", "persistent_changed_partner",
                     "extended_only", "contracted_only"]
ALL_CLASSES = INTERFACE_CLASSES + ["non_interface"]

# ---- the two filtering rules, stated once -----------------------------------
# NATIVE:      drop residues whose classification depends on a symmetry-generated
#              contact (uses_symmetry_generated is True). Such a residue leaves
#              the table entirely; it is not reclassified.
# ENGINEERING: drop residues within ENGINEERING_SHELL of an internal indel
#              (engineering_excluded is True). pilot.py already applied this to
#              its counts; it is applied here identically, to every family.
NATIVE_COL = "uses_symmetry_generated"
ENGINEERING_COL = "engineering_excluded"

# ---- independently verified expectations ------------------------------------
# Verified 2026-10-12 directly from the residue CSVs, before this script was
# written. Every block below must reproduce exactly or the script aborts.
# keys: (system, entity) -> {class: n}
EXPECT_NATIVE_ENGINEERING = {          # the regenerated canonical counts
    ("algoriphagus", "P1_Putative phage tail sheath protein FI_692aa"):
        {"persistent_same_partner": 86, "persistent_changed_partner": 2,
         "extended_only": 45, "contracted_only": 24, "non_interface": 354},
    ("vibrio", "P1_Type VI secretion protein_473aa"):
        {"persistent_same_partner": 73, "persistent_changed_partner": 17,
         "extended_only": 4, "contracted_only": 50, "non_interface": 178},
    ("vibrio", "P2_Type VI secretion protein_159aa"):
        {"persistent_same_partner": 56, "persistent_changed_partner": 7,
         "extended_only": 4, "contracted_only": 3, "non_interface": 25},
}
EXPECT_NATIVE_ONLY = {                 # native rule, engineering shell NOT applied
    ("algoriphagus", "P1_Putative phage tail sheath protein FI_692aa"):
        {"persistent_same_partner": 86, "persistent_changed_partner": 2,
         "extended_only": 45, "contracted_only": 24, "non_interface": 354},
    ("vibrio", "P1_Type VI secretion protein_473aa"):
        {"persistent_same_partner": 73, "persistent_changed_partner": 17,
         "extended_only": 4, "contracted_only": 50, "non_interface": 178},
    ("vibrio", "P2_Type VI secretion protein_159aa"):   # README prose table
        {"persistent_same_partner": 57, "persistent_changed_partner": 7,
         "extended_only": 4, "contracted_only": 4, "non_interface": 25},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def residue_path(system: str, entity: str) -> Path:
    """Reproduce pilot.py's output filename rule exactly."""
    safe = "".join(ch if ch.isalnum() else "_" for ch in entity)[:40]
    return OUT_DIR / f"residues_{system}_{safe}.csv"


def counts(df: pd.DataFrame) -> dict:
    vc = df["persistence_class"].value_counts().to_dict()
    return {cls: int(vc.get(cls, 0)) for cls in ALL_CLASSES}


def smallest(c: dict) -> int:
    return min(c[x] for x in INTERFACE_CLASSES)


def verdict(c: dict) -> str:
    return "PASS" if smallest(c) >= MIN_CLASS_SIZE else "FAIL"


def check(label: str, got: dict, want: dict, key) -> None:
    if got != want:
        sys.exit(f"ABORT: {label} mismatch for {key}\n  computed: {got}\n  expected: {want}\n"
                 "The residue tables are not the ones these expectations were "
                 "verified against. Nothing was written.")


def atomic_write(path: Path, text: str) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="run every check and print the report, write nothing")
    args = ap.parse_args()

    archived_summary = ARCHIVE / "summary.csv"
    if not archived_summary.exists():
        sys.exit(f"ABORT: archived run not found at {archived_summary}. The original "
                 "symmetry-inclusive outputs must be preserved before regenerating.")
    arch = pd.read_csv(archived_summary)

    rows, blocks = [], []
    inputs = [(archived_summary, sha256(archived_summary))]

    for _, a in arch.iterrows():
        system, entity = a["system"], a["entity"]
        key = (system, entity)
        rp = residue_path(system, entity)
        if not rp.exists():
            sys.exit(f"ABORT: residue table missing for {key}: {rp}")
        inputs.append((rp, sha256(rp)))
        df = pd.read_csv(rp)

        for col in ("persistence_class", NATIVE_COL, ENGINEERING_COL):
            if col not in df.columns:
                sys.exit(f"ABORT: {rp.name} lacks required column {col!r}")

        native_mask = ~df[NATIVE_COL].astype(bool)
        eng_mask = ~df[ENGINEERING_COL].astype(bool)

        c_all_eng = counts(df[eng_mask])                 # what pilot.py counted
        c_native = counts(df[native_mask])               # native, no eng shell
        c_final = counts(df[native_mask & eng_mask])     # native + eng shell

        # ---- assertions ----------------------------------------------------
        # 1. the archived run is reproducible from these residue tables, so the
        #    tables and the archived verdict are the same run
        archived_counts = {cls: int(a[f"n_{cls}"]) for cls in ALL_CLASSES}
        check("archived (all-contacts + engineering)", c_all_eng, archived_counts, key)
        # 2. the native rule reproduces the independently verified counts
        check("native-only", c_native, EXPECT_NATIVE_ONLY[key], key)
        check("native + engineering", c_final, EXPECT_NATIVE_ENGINEERING[key], key)
        # 3. the native rule can never pass a family the archived rule failed
        assert smallest(c_final) <= smallest(c_all_eng), (
            "native filtering increased the smallest class -- impossible")

        n_dropped_sym = int(df[NATIVE_COL].astype(bool).sum())
        n_dropped_eng = int(df[ENGINEERING_COL].astype(bool).sum())

        rows.append({
            # ---- carried over verbatim from the archived run (rule-independent)
            "system": system, "entity": entity, "identity": a["identity"],
            "modelled_ext": a["modelled_ext"], "modelled_con": a["modelled_con"],
            "mapped": a["mapped"], "substitutions": a["substitutions"],
            "unmatched_internal": a["unmatched_internal"],
            "rings_ext": a["rings_ext"], "rings_con": a["rings_con"],
            "lattice_failures_ext": a["lattice_failures_ext"],
            "lattice_failures_con": a["lattice_failures_con"],
            # ---- recomputed under the native rule
            "counting_rule": "native_coordinates_only+engineering_shell",
            "engineering_excluded": n_dropped_eng,
            "symmetry_dependent_excluded": n_dropped_sym,
            "symmetry_generated_used": False,
            **{f"n_{x}": c_final[x] for x in ALL_CLASSES},
            "smallest_interface_class": smallest(c_final),
            "verdict": verdict(c_final),
            # ---- the superseded numbers, kept inline for auditability
            "archived_smallest_interface_class": int(a["smallest_interface_class"]),
            "archived_verdict": a["verdict"],
        })

        persistent = c_final["persistent_same_partner"] + c_final["persistent_changed_partner"]
        state_specific = c_final["extended_only"] + c_final["contracted_only"]
        blocks.append("\n".join([
            "",
            f"  --- {system} / {entity} ---",
            f"      residue table         : {rp.name}",
            f"      mapped residues       : {int(a['mapped'])}",
            f"      excluded, symmetry-dependent : {n_dropped_sym}",
            f"      excluded, engineering shell  : {n_dropped_eng} "
            f"(within {ENGINEERING_SHELL} A of an internal indel)",
            "",
            f"      persistence class      {'archived':>10s} {'native':>10s} "
            f"{'NATIVE+ENG':>12s}",
            *[f"        {cls:30s} {c_all_eng[cls]:10d} {c_native[cls]:10d} "
              f"{c_final[cls]:12d}" for cls in ALL_CLASSES],
            "",
            f"      smallest interface class = {smallest(c_final)} "
            f"(threshold {MIN_CLASS_SIZE}) -> {verdict(c_final)}"
            f"      [archived rule gave {int(a['smallest_interface_class'])} "
            f"-> {a['verdict']}]",
            f"      [POST HOC] collapsed three-way scheme: persistent={persistent}, "
            f"state_specific={state_specific}, "
            f"non_interface={c_final['non_interface']} "
            f"-> smallest={min(persistent, state_specific)}",
            "      (post hoc, chosen after seeing the counts -- NOT the verdict)",
        ]))

    sdf = pd.DataFrame(rows)
    n_pass = int((sdf["verdict"] == "PASS").sum())

    head = [
        "CONTRACTILE SHEATH INTERFACE-PERSISTENCE STRUCTURAL PILOT",
        "VERDICT ARTIFACTS REGENERATED UNDER THE NATIVE-COORDINATE RULE",
        "",
        f"generated   : {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"script      : {SCRIPT_VERSION}",
        "supersedes  : output/archive_symmetry_inclusive_2026-10-10/ "
        "(original pilot.py run, 2026-10-10 23:20)",
        "",
        "AUTHORITATIVE. The archived run counted residues whose interface class",
        "depends on symmetry-generated contacts and therefore reported 1/3 passing.",
        "",
        "Counting rule applied here, to every family identically:",
        "  (1) NATIVE      -- a residue whose classification depends on a contact to a",
        "                     symmetry-generated neighbour is dropped from the table.",
        "                     It is not reclassified.",
        f"  (2) ENGINEERING -- a residue within {ENGINEERING_SHELL} A of an internal indel is",
        "                     dropped (as pilot.py already did for its own counts).",
        "",
        "Pre-registered parameters, unchanged:",
        f"  contact cutoff          {CONTACT_CUTOFF} A heavy atom",
        f"  minimum class size      {MIN_CLASS_SIZE} residues",
        f"  engineering shell       {ENGINEERING_SHELL} A",
        "  partner identity        (canonical entity, screw index n, C6 index k)",
        "",
        "No evolutionary rates are estimated in this pilot. No scientific rule was",
        "changed by this regeneration; only its consistent application.",
        "",
        "Scope: this is the historical CATEGORICAL pilot. The frozen continuous",
        "remodeling metric R_i lives in qc_output/ (see QC_RESULTS.md) and is not",
        "affected by, and does not depend on, anything in this directory.",
        "",
        "Inputs (sha256):",
        *[f"  {h}  {p.relative_to(OUT_DIR.parent)}" for p, h in inputs],
        "",
        "Recomputed here      : class counts, smallest class, verdict, post-hoc collapse.",
        "Carried over verbatim: identity, modelled/mapped counts, substitutions,",
        "                       unmatched internal, ring counts, lattice failures",
        "                       (structural metadata, independent of counting rule).",
        "NOT regenerated      : per-state structural diagnostics, contact-cutoff sweep,",
        "                       per-class covariate means. These need the coordinates;",
        "                       they remain in the archived report and were computed",
        "                       there WITHOUT the native filter. Do not quote them as",
        "                       native-rule numbers.",
        "",
        "=" * 78,
        "PER-FAMILY COUNTS",
        "=" * 78,
        "",
    ]
    tail = [
        "",
        "=" * 78,
        "OVERALL (native-coordinate rule)",
        "=" * 78,
        sdf[["system", "entity", "mapped"] + [f"n_{x}" for x in INTERFACE_CLASSES]
            + ["smallest_interface_class", "verdict"]].to_string(index=False),
        "",
        f"  protein families passing: {n_pass} / {len(sdf)}",
        "",
        "  The pre-registered four-way gate is not met in any family. The decisive",
        "  class is persistent_changed_partner for Algoriphagus (2) and the",
        "  state-specific classes for Vibrio: the lattice keeps its connectivity",
        "  through contraction, so subunits slide and tilt against the same (n, k)",
        "  neighbours rather than swapping partners. The class is well defined and",
        "  simply not populated -- a property of the mechanism, not of the data.",
    ]
    text = "\n".join(head + blocks + tail) + "\n"

    if args.dry_run:
        print(text)
        print("[dry run] nothing written")
        return 0

    sdf.to_csv(OUT_DIR / "summary.csv.tmp", index=False)
    os.replace(OUT_DIR / "summary.csv.tmp", OUT_DIR / "summary.csv")
    atomic_write(OUT_DIR / "qc_report.txt", text)
    print(text)
    print(f"wrote {OUT_DIR / 'summary.csv'}")
    print(f"wrote {OUT_DIR / 'qc_report.txt'}")
    print(f"families passing: {n_pass} / {len(sdf)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
