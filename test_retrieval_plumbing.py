#!/usr/bin/env python3
"""End-to-end test of the retrieval PLUMBING only. No protocol content.

`blast_retrieve.py`'s submit / persist / resume / backoff / terminal-state
machinery can be observed on the real refseq_protein run, but the
download -> parse -> validate -> persist leg only executes when a search
actually returns, and refseq_protein searches have been sitting in NCBI's queue
for over an hour. This harness exercises that leg against a small database so it
is verified before the real submissions are made.

Two modes:

  --offline (default)  builds a synthetic BlastOutput fixture and drives the
                       real `validate()` / `atomic_write_bytes()` / re-read path
                       with it, including corrupted bodies that must be
                       rejected. Deterministic, no network. The fixture is
                       SYNTHETIC -- it proves the parse/persist/validation
                       logic, not that NCBI returned anything.

  --live               submits ONE query to `swissprot` with blast_retrieve's own
                       `submit()`, polls with its `poll()`, fetches with its
                       `fetch()`, then runs the same checks on the real body.
                       Needs NCBI's queue to actually return; as of 2026-10-12 a
                       swissprot blastp sat WAITING for 22 min / 37 polls, so
                       this mode may legitimately time out with the plumbing
                       still correct.

WHAT IT IS NOT: this is not part of the frozen retrieval protocol and produces
no scientific output. `swissprot` is used because it answers in minutes; the
frozen protocol's database is refseq_protein and is not changed by this file.
Nothing is written inside `retrieval/`; everything lands in a temp directory.

Usage: python3 test_retrieval_plumbing.py [--live] [--max-minutes 20]
Exit 0 = the plumbing works end to end. Non-zero = it does not; read the output.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import sys
import tempfile
import time
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import blast_retrieve as br  # noqa: E402

TEST_DB = "swissprot"


def fixture(n_hits: int = 60) -> str:
    """A synthetic but structurally real blastp BlastOutput document."""
    hits = []
    for i in range(1, n_hits + 1):
        hits.append(f"""    <Hit>
      <Hit_num>{i}</Hit_num>
      <Hit_id>ref|WP_{100000000 + i}.1|</Hit_id>
      <Hit_def>phage tail sheath family protein [synthetic organism {i}]</Hit_def>
      <Hit_accession>WP_{100000000 + i}</Hit_accession>
      <Hit_len>386</Hit_len>
      <Hit_hsps>
        <Hsp>
          <Hsp_num>1</Hsp_num>
          <Hsp_bit-score>{700 - i}</Hsp_bit-score>
          <Hsp_evalue>{1e-180 * i:.3e}</Hsp_evalue>
          <Hsp_identity>{386 - i}</Hsp_identity>
          <Hsp_align-len>386</Hsp_align-len>
          <Hsp_qseq>{'A' * 386}</Hsp_qseq>
          <Hsp_hseq>{'A' * 386}</Hsp_hseq>
          <Hsp_midline>{'+' * 386}</Hsp_midline>
        </Hsp>
      </Hit_hsps>
    </Hit>""")
    return f"""<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>{br.PROGRAM}</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.17.0+</BlastOutput_version>
  <BlastOutput_db>{br.DATABASE}</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>synthetic fixture query</BlastOutput_query-def>
  <BlastOutput_query-len>386</BlastOutput_query-len>
  <BlastOutput_iterations>
   <Iteration>
    <Iteration_iter-num>1</Iteration_iter-num>
    <Iteration_hits>
{chr(10).join(hits)}
    </Iteration_hits>
   </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--live", action="store_true",
                    help=f"submit a real query to {TEST_DB} instead of using the fixture")
    ap.add_argument("--max-minutes", type=float, default=20.0)
    args = ap.parse_args()

    tmp = Path(tempfile.mkdtemp(prefix="retrieval_plumbing_"))
    # redirect every side effect away from the real run
    br.OUT = tmp
    br.LOG_PATH = tmp / "test.log"
    if args.live:
        br.DATABASE = TEST_DB        # test-only; validate() checks against this
        br.POLL_MIN = 30
        br.POLL_MAX = 60

    refs = json.loads((Path(__file__).parent / "reference_sequences.json").read_text())
    fam = "r2pyocin"
    seq = refs[fam]["seq"]
    print(f"mode       : {'live' if args.live else 'offline (synthetic fixture)'}")
    print(f"temp dir   : {tmp}")
    print(f"database   : {br.DATABASE}"
          + ("  (test only -- protocol uses refseq_protein)" if args.live else ""))
    print(f"query      : {fam}, {len(seq)} aa\n")

    checks: list[tuple[str, bool, str]] = []

    def record(label: str, ok: bool, detail: str = "") -> None:
        checks.append((label, ok, detail))
        print(f"  [{'ok ' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""))

    try:
        if args.live:
            rid = br.submit(fam, seq)
            record("submit returns an RID", bool(rid), rid)

            deadline = time.monotonic() + args.max_minutes * 60
            state, detail, waits = "waiting", "", 0
            while time.monotonic() < deadline:
                state, detail = br.poll(rid)
                print(f"    poll {waits + 1}: {state} ({detail})", flush=True)
                if state != "waiting":
                    break
                waits += 1
                time.sleep(br.POLL_MIN)
            record("poll reaches a non-waiting state", state != "waiting",
                   f"{state}: {detail}")
            if state != "ready":
                record("search returned hits", False, f"state={state}")
                raise SystemExit(1)
            xml = br.fetch(rid)
            record("fetch returns a body", len(xml) > 0, f"{len(xml.encode())} bytes")
        else:
            rid = "FIXTURE00000"
            xml = fixture()
            record("fixture built", len(xml.encode()) > br.MIN_RESULT_BYTES,
                   f"{len(xml.encode())} bytes, synthetic")

        n_hits, qdef = br.validate(xml, rid)
        record("validate accepts the result", n_hits > 0,
               f"{n_hits} hits, query-def {qdef[:40]!r}")

        # deliberately unusable bodies must all be rejected
        bads = [("truncated body", xml[:len(xml) // 3]),
                ("non-XML body", "<html>maintenance</html>" * 600),
                ("undersized body", xml[:200]),
                ("wrong root element", xml.replace("BlastOutput>", "Something>")
                                          .replace("<BlastOutput ", "<Something ")),
                ("wrong program", xml.replace(f"<BlastOutput_program>{br.PROGRAM}<",
                                              "<BlastOutput_program>blastn<")),
                ("wrong database", xml.replace(f"<BlastOutput_db>{br.DATABASE}<",
                                               "<BlastOutput_db>nr<")),
                # padded so it clears the size floor and reaches the Hit check
                ("hitless body", re.sub(r"<Hit>.*</Hit>",
                                        "<!--" + "x" * 11000 + "-->",
                                        xml, flags=re.S))]
        for label, bad in bads:
            try:
                br.validate(bad, rid)
                record(f"validate rejects a {label}", False, "it accepted it")
            except br.Malformed as exc:
                record(f"validate rejects a {label}", True, str(exc)[:60])

        data = xml.encode()
        path = tmp / f"{fam}.xml"
        br.atomic_write_bytes(path, data)
        on_disk = path.read_bytes()
        record("file persisted", path.exists(), f"{path.name}, {len(on_disk)} bytes")
        record("bytes on disk match fetched bytes",
               hashlib.sha256(on_disk).hexdigest() == hashlib.sha256(data).hexdigest())
        record("no .tmp left behind", not list(tmp.glob("*.tmp")))
        n2, _ = br.validate(on_disk.decode(), rid)
        record("file on disk re-validates", n2 == n_hits, f"{n2} hits")
        record("result clears the size floor", len(on_disk) >= br.MIN_RESULT_BYTES,
               f"{len(on_disk)} >= {br.MIN_RESULT_BYTES}")
    finally:
        failed = [c for c in checks if not c[1]]
        print(f"\n{len(checks) - len(failed)}/{len(checks)} checks passed")
        if not failed:
            shutil.rmtree(tmp, ignore_errors=True)
            print("temp dir removed")
        else:
            print(f"temp dir kept for inspection: {tmp}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
