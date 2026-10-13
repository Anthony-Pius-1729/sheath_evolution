#!/usr/bin/env python3
"""Step 1 of the frozen retrieval protocol: broad homolog retrieval.

One search per structural reference, each against refseq_protein. RefSeq is
chosen over nr deliberately: WP_ accessions resolve to genomic locations through
NCBI's Identical Protein Groups, which is what the synteny stage needs. Searching
nr would return many sequences with no usable genomic context.

Nothing is filtered here beyond BLAST's own E-value. The 30-90% identity window
and every other gate is applied downstream, so the raw retrieval stays auditable.
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "retrieval"
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

EXPECT = 1e-5
HITLIST = 1000
SUBMIT_GAP = 12      # NCBI asks for >=10 s between submissions
POLL_GAP = 60        # and >=60 s between polls of one job


def post(params):
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(BLAST, data=data,
                                 headers={"User-Agent": "sheath-evo-pilot/1.0"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=180) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as exc:
            print(f"      request failed ({exc}); retry {attempt + 1}/5", flush=True)
            time.sleep(20)
    raise RuntimeError("BLAST request failed five times")


def submit(name, seq):
    txt = post({"CMD": "Put", "PROGRAM": "blastp", "DATABASE": "refseq_protein",
                "QUERY": seq, "HITLIST_SIZE": HITLIST, "EXPECT": EXPECT,
                "FILTER": "F"})
    m = re.search(r"^\s*RID = (\S+)", txt, re.M)
    if not m:
        raise RuntimeError(f"no RID returned for {name}")
    print(f"  submitted {name}: RID={m.group(1)}", flush=True)
    return m.group(1)


def ready(rid):
    txt = post({"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"})
    if "Status=READY" in txt:
        return "ready" if "ThereAreHits=yes" in txt else "nohits"
    if "Status=FAILED" in txt or "Status=UNKNOWN" in txt:
        return "failed"
    return "waiting"


def fetch_xml(rid):
    return post({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML",
                 "ALIGNMENTS": HITLIST, "DESCRIPTIONS": HITLIST})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    refs = json.load(open(HERE / "reference_sequences.json"))

    jobs = {}
    for name, rec in refs.items():
        done = OUT / f"{name}.xml"
        if done.exists() and done.stat().st_size > 10000:
            print(f"  {name}: already retrieved, skipping", flush=True)
            continue
        jobs[name] = submit(name, rec["seq"])
        time.sleep(SUBMIT_GAP)

    pending = dict(jobs)
    while pending:
        for name, rid in list(pending.items()):
            state = ready(rid)
            print(f"    {name} ({rid}): {state}", flush=True)
            if state == "ready":
                xml = fetch_xml(rid)
                (OUT / f"{name}.xml").write_text(xml)
                print(f"  {name}: saved {len(xml)} bytes", flush=True)
                pending.pop(name)
            elif state in ("failed", "nohits"):
                print(f"  {name}: {state.upper()}", flush=True)
                pending.pop(name)
            time.sleep(5)
        if pending:
            time.sleep(POLL_GAP)

    print("retrieval complete", flush=True)


if __name__ == "__main__":
    main()
