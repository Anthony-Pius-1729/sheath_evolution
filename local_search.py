#!/usr/bin/env python3
"""Step 1b + steps 2-3 for the UniProt route: local blastp, then the frozen gates.

Searches the structural reference against the locally built candidate database
from `uniprot_retrieve.py`, then applies the SAME pre-registered gates as the
NCBI route by importing them from `filter_hits.py`:

    identity window      30-90% to the structural reference   (ID_LO, ID_HI)
    reference coverage   >= 80% of the reference length        (MIN_REF_COV)
    length compatibility 0.75-1.30x the reference              (LEN_LO, LEN_HI)
    dereplication        one representative per species, best bitscore

There is no second definition of any gate in this file -- the constants are
imported, so a change to the frozen protocol can only be made in one place.
Identity and coverage are computed exactly as `filter_hits.parse` computes them
from BLAST XML: identity = nident/align_len, reference coverage =
(qend-qstart+1)/reference_length, best HSP per subject chosen by bitscore.

Output goes to `filtered/{family}_filtered.csv` in the schema the synteny stage
already reads, plus `source` and `accession_source` columns so a row's route is
never ambiguous. `accession` carries the RefSeq accession where UniProt provides
one, because that is what `synteny.py` resolves through NCBI Identical Protein
Groups; rows without one keep their UniProt accession and are labelled
accordingly -- synteny will record them as having no genomic location, which is
its documented behaviour and not a silent loss.

Usage: python3 local_search.py --families r2pyocin
Needs: blastp + makeblastdb on PATH, pandas.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.resolve()))
import filter_hits as fh  # noqa: E402  -- the single definition of the frozen gates

CAND = HERE / "retrieval_uniprot"
OUT = HERE / "filtered"

EVALUE = fh_evalue = 1e-5         # same E-value ceiling as the frozen BLAST step
MAX_TARGETS = 20000
OUTFMT = ("6 qseqid sseqid pident nident length mismatch gapopen qstart qend "
          "sstart send slen evalue bitscore")
COLS = OUTFMT.split()[1:]


def iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def tool_version(tool: str) -> str:
    try:
        out = subprocess.run([tool, "-version"], capture_output=True, text=True,
                             timeout=60).stdout.strip().splitlines()
        return out[0] if out else "?"
    except Exception as exc:
        sys.exit(f"ABORT: {tool} is not usable ({exc}). Install BLAST+ first.")


def species_key(organism: str) -> str:
    """Identical to filter_hits: strip strain detail to a species-level key."""
    return " ".join(organism.replace("Candidatus ", "").split()[:2])


def write_atomic(path: Path, text: str) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--families", default=None, help="comma-separated subset")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    refs = json.loads((HERE / "reference_sequences.json").read_text())
    families = ([f.strip() for f in args.families.split(",")] if args.families
                else list(refs))

    blastp_v, makedb_v = tool_version("blastp"), tool_version("makeblastdb")
    print(f"{blastp_v} / {makedb_v}")

    for fam in families:
        if fam not in refs:
            sys.exit(f"ABORT: unknown family {fam!r}")
        cand_fa = CAND / f"{fam}_candidates.fasta"
        cand_tsv = CAND / f"{fam}_candidates.tsv"
        for p in (cand_fa, cand_tsv):
            if not p.exists():
                sys.exit(f"ABORT: {p} missing. Run uniprot_retrieve.py first.")

        ref_seq = refs[fam]["seq"]
        ref_len = len(ref_seq)
        ref_uniprot = (refs[fam].get("uniprot") or [None])[0]
        meta = pd.read_csv(cand_tsv, sep="\t", dtype=str).fillna("")
        print(f"\n=== {fam}: {len(meta)} candidates, reference {ref_len} aa ===")

        work = Path(tempfile.mkdtemp(prefix=f"blast_{fam}_"))
        try:
            subprocess.run(["makeblastdb", "-in", str(cand_fa), "-dbtype", "prot",
                            "-out", str(work / "db"), "-title", f"{fam}_candidates"],
                           check=True, capture_output=True, text=True)
            qpath = work / "query.fasta"
            qpath.write_text(f">REFERENCE_{fam}_{ref_uniprot}\n{ref_seq}\n")
            res = subprocess.run(
                ["blastp", "-query", str(qpath), "-db", str(work / "db"),
                 "-evalue", str(EVALUE), "-max_target_seqs", str(MAX_TARGETS),
                 "-seg", "no", "-outfmt", OUTFMT],
                check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            sys.exit(f"ABORT: {exc.cmd[0]} failed: {exc.stderr[:400]}")
        finally:
            db_files = sorted(p.name for p in work.glob("db.*"))

        rows = [dict(zip(COLS, line.split("\t")))
                for line in res.stdout.strip().split("\n") if line.strip()]
        if not rows:
            sys.exit(f"ABORT: {fam}: blastp returned no hits at E <= {EVALUE}")
        hsp = pd.DataFrame(rows)
        for c in ("nident", "length", "qstart", "qend", "slen"):
            hsp[c] = hsp[c].astype(int)
        for c in ("evalue", "bitscore", "pident"):
            hsp[c] = hsp[c].astype(float)

        # best HSP per subject, by bitscore -- as filter_hits.parse does per Hit
        hsp = hsp.sort_values("bitscore", ascending=False).drop_duplicates("sseqid")
        hsp["pident"] = hsp["nident"] / hsp["length"].clip(lower=1)
        hsp["ref_cov"] = (hsp["qend"] - hsp["qstart"] + 1) / ref_len
        hsp["len_ratio"] = hsp["slen"] / ref_len

        m = meta.rename(columns={"Entry": "uniprot_accession",
                                 "Organism": "organism",
                                 "Protein names": "title",
                                 "Organism (ID)": "taxid",
                                 "RefSeq": "refseq", "EMBL": "embl",
                                 "Proteomes": "proteome", "Reviewed": "reviewed",
                                 "Pfam": "pfam", "InterPro": "interpro",
                                 "Taxonomic lineage": "lineage"})
        keep_cols = [c for c in ("uniprot_accession", "organism", "title", "taxid",
                                 "refseq", "embl", "proteome", "reviewed",
                                 "pfam", "interpro", "lineage")
                     if c in m.columns]
        df = hsp.merge(m[keep_cols], left_on="sseqid",
                       right_on="uniprot_accession", how="left")
        missing = df["uniprot_accession"].isna().sum()
        if missing:
            sys.exit(f"ABORT: {missing} hits have no metadata row; the candidate "
                     "FASTA and TSV are out of step.")

        df["species"] = df["organism"].fillna("?").map(species_key)
        df["refseq_first"] = (df["refseq"].fillna("").str.split(";").str[0]
                              .str.split("[").str[0].str.strip())
        df["accession"] = df.apply(
            lambda r: r["refseq_first"] or r["uniprot_accession"], axis=1)
        df["accession_source"] = df["refseq_first"].map(
            lambda v: "refseq" if v else "uniprot")
        df["slen"] = df["slen"].astype(int)
        df["source"] = "uniprotkb_pfam_family_scan+local_blastp"

        # the full hit table, pre-gate, so the gates' effect is auditable
        df.to_csv(CAND / f"{fam}_blastp_hits.csv", index=False)

        n_hits = len(df)
        is_ref = df["uniprot_accession"] == ref_uniprot
        ref_row = df[is_ref]
        gate = ((df.pident >= fh.ID_LO) & (df.pident <= fh.ID_HI)
                & (df.ref_cov >= fh.MIN_REF_COV)
                & (df.len_ratio >= fh.LEN_LO) & (df.len_ratio <= fh.LEN_HI))
        passed = df[gate].copy()
        n_gate = len(passed)
        final = (passed.sort_values("bitscore", ascending=False)
                 .drop_duplicates("species").reset_index(drop=True))
        final["family"] = fam

        out_cols = ["accession", "accession_source", "uniprot_accession", "title",
                    "pfam", "interpro", "lineage",
                    "organism", "species", "taxid", "slen", "bitscore", "pident",
                    "ref_cov", "len_ratio", "evalue", "refseq", "embl", "proteome",
                    "reviewed", "family", "source"]
        final[[c for c in out_cols if c in final.columns]].to_csv(
            OUT / f"{fam}_filtered.csv", index=False)

        prov = {
            "family": fam, "searched_at": iso(),
            "route": "uniprotkb_pfam_family_scan+local_blastp",
            "reference": {"uniprot": ref_uniprot, "length": ref_len},
            "candidate_db": {
                "fasta": cand_fa.name, "sha256": sha256(cand_fa),
                "sequences": int(len(meta)), "files": db_files},
            "tools": {"blastp": blastp_v, "makeblastdb": makedb_v},
            "search_params": {"evalue": EVALUE, "max_target_seqs": MAX_TARGETS,
                              "seg": "no", "outfmt": OUTFMT,
                              "matrix": "BLOSUM62 (blastp default)"},
            "frozen_gates": {"ID_LO": fh.ID_LO, "ID_HI": fh.ID_HI,
                             "MIN_REF_COV": fh.MIN_REF_COV,
                             "LEN_LO": fh.LEN_LO, "LEN_HI": fh.LEN_HI,
                             "imported_from": "filter_hits.py"},
            "funnel": {"candidates": int(len(meta)),
                       "blastp_hits_evalue": n_hits,
                       "pass_frozen_gates": n_gate,
                       "species_dereplicated": int(len(final))},
            "reference_recovery": {
                "in_candidate_set": bool(is_ref.any()),
                "self_hit_pident": (round(float(ref_row.pident.iloc[0]), 4)
                                    if len(ref_row) else None),
                "self_hit_ref_cov": (round(float(ref_row.ref_cov.iloc[0]), 4)
                                     if len(ref_row) else None),
                "passes_gates": bool(gate[is_ref].any()) if is_ref.any() else None,
                "note": "the reference is expected to be retrieved and then "
                        "excluded by the 90% identity ceiling; that is the gate "
                        "working, not a retrieval failure",
            },
            "identity_distribution": {
                k: round(float(v), 4) for k, v in
                final.pident.describe(percentiles=[.1, .25, .5, .75, .9]).items()},
            "coverage_distribution": {
                k: round(float(v), 4) for k, v in
                final.ref_cov.describe(percentiles=[.1, .25, .5, .75, .9]).items()},
            "genomic_provenance": {
                "accession_source_counts":
                    final.accession_source.value_counts().to_dict(),
                "with_embl_xref": int((final.embl.fillna("") != "").sum()),
                "with_proteome": int((final.proteome.fillna("") != "").sum()),
                "distinct_species": int(final.species.nunique()),
                "distinct_genera": int(final.species.str.split().str[0].nunique())},
        }
        write_atomic(CAND / f"{fam}_search_provenance.json",
                     json.dumps(prov, indent=2))
        shutil.rmtree(work, ignore_errors=True)

        print(f"  candidates {len(meta)} -> blastp hits {n_hits} "
              f"-> gates {n_gate} -> species {len(final)}")
        print(f"  identity  median {final.pident.median():.3f} "
              f"range {final.pident.min():.3f}-{final.pident.max():.3f}")
        print(f"  coverage  median {final.ref_cov.median():.3f} "
              f"min {final.ref_cov.min():.3f}")
        print(f"  genera {final.species.str.split().str[0].nunique()}, "
              f"RefSeq-resolvable {int((final.accession_source == 'refseq').sum())}"
              f"/{len(final)}")
        print(f"  reference self-hit: pident "
              f"{prov['reference_recovery']['self_hit_pident']}, passes gates "
              f"{prov['reference_recovery']['passes_gates']}")
        print(f"  wrote {OUT / (fam + '_filtered.csv')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
