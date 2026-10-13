#!/usr/bin/env python3
"""Step 1, alternative route: candidate retrieval from UniProtKB.

WHY THIS EXISTS
---------------
The frozen route (`blast_retrieve.py`, NCBI BLAST URL API) is intact and remains
the fallback. It is not the blocker's fault: on 2026-10-12 a single blastp of the
386 aa R2 pyocin reference against refseq_protein sat at Status=WAITING for 131
min across 17 polls, and a control of the same query against swissprot -- a
database small enough to answer in under a minute under normal load -- also
stayed WAITING for 22 min across 37 polls. The endpoint is queueing or
throttling, and no amount of client-side correctness fixes that.

This route retrieves the same kind of candidates without a remote queue, and is
more reproducible for it: the UniProt release, the exact query and the local
search parameters are all recorded with the output, so the candidate set can be
reconstructed exactly.

WHAT IT DOES
------------
  1. reads the reference's own UniProt accession from reference_sequences.json;
  2. reads that entry's Pfam signatures FROM UNIPROT -- the family definition is
     taken from the reference, never hardcoded;
  3. streams every bacterial UniProtKB entry carrying any of those signatures,
     with sequence, organism, taxid and the RefSeq/EMBL/proteome cross-references
     that the synteny stage needs to resolve a genomic locus;
  4. writes the candidate set, its metadata and a provenance record.

NO BIOLOGICAL SELECTION HAPPENS HERE. This is the retrieval stage, and it is
deliberately broad: no identity window, no coverage or length gate, no
dereplication. The frozen gates (30-90% identity, >=80% reference coverage,
0.75-1.30x length, one representative per species) are applied downstream by
`local_search.py`, which imports them from `filter_hits.py` so there is one
definition of each. The hypothesis, the gates and the orthology rules are
unchanged by this file.

KNOWN LIMITATION, STATED PLAINLY
--------------------------------
A profile-family query has different recall from a BLAST search. A true ortholog
whose UniProt entry carries none of the reference's Pfam signatures will not be
retrieved here, where BLAST would have found it on sequence similarity alone.
The reverse is also true -- this route reaches entries BLAST ranked below its
hit-list cutoff. Treat the two as different retrieval routes with different
recall, not as substitutes, and report which one produced a given dataset.

Usage:
  python3 uniprot_retrieve.py --families r2pyocin
  python3 uniprot_retrieve.py                      # all five
Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "retrieval_uniprot"
UNIPROT = "https://rest.uniprot.org"
PAGE = 500
PAUSE = 0.34
TIMEOUT = 180
TRIES = 4
USER_AGENT = "sheath-evo-pilot/2.0"

TAXON_BACTERIA = 2          # retrieval scope: Bacteria
FIELDS = ["accession", "id", "reviewed", "protein_name", "organism_name",
          "organism_id", "length", "sequence", "xref_refseq", "xref_embl",
          "xref_proteomes", "lineage",
          # domain architecture: needed to tell orthologous sheaths from other
          # contractile-system sheath families that share a Pfam hit
          "xref_pfam", "xref_interpro"]


def iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def fetch(url: str) -> tuple[str, dict]:
    last = None
    for attempt in range(TRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as fh:
                # normalise header case -- a plain dict() of the message would
                # make these lookups case-sensitive
                return (fh.read().decode("utf-8"),
                        {k.lower(): v for k, v in fh.headers.items()})
        except Exception as exc:
            last = f"{type(exc).__name__}: {exc}"
            if attempt + 1 < TRIES:
                time.sleep(3 + 4 * attempt)
    raise RuntimeError(f"UniProt request failed after {TRIES} tries: {last}\n  {url}")


def entry(accession: str) -> dict:
    txt, _ = fetch(f"{UNIPROT}/uniprotkb/{accession}.json")
    return json.loads(txt)


def pfam_signatures(rec: dict) -> list[str]:
    return sorted({x["id"] for x in rec.get("uniProtKBCrossReferences", [])
                   if x["database"] == "Pfam"})


def stream_search(query: str) -> tuple[list[dict], dict]:
    """Paginate a UniProtKB search, following the RFC-5988 'next' link."""
    url = (f"{UNIPROT}/uniprotkb/search?"
           + urllib.parse.urlencode({"query": query, "format": "tsv",
                                     "fields": ",".join(FIELDS), "size": PAGE}))
    rows, header, meta, page = [], None, {}, 0
    while url:
        txt, hdrs = fetch(url)
        page += 1
        if not meta:
            meta = {"release": hdrs.get("x-uniprot-release"),
                    "release_date": hdrs.get("x-uniprot-release-date"),
                    "total_results": hdrs.get("x-total-results")}
        lines = txt.rstrip("\n").split("\n")
        if header is None:
            header = lines[0].split("\t")
        body = lines[1:] if lines and lines[0].startswith(header[0]) else lines
        for line in body:
            if line.strip():
                rows.append(dict(zip(header, line.split("\t"))))
        # the Link header's URL contains commas (the fields list), so it cannot
        # be split on commas -- match the bracketed URL directly
        m = re.search(r'<([^>]+)>\s*;\s*rel="next"', hdrs.get("link", ""))
        url = m.group(1) if m else None
        print(f"    page {page}: {len(rows)} records", flush=True)
        if url:
            time.sleep(PAUSE)
    meta["pages"] = page
    return rows, meta


def write_atomic(path: Path, text: str) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--families", default=None, help="comma-separated subset")
    ap.add_argument("--taxon", type=int, default=TAXON_BACTERIA)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    refs = json.loads((HERE / "reference_sequences.json").read_text())
    families = ([f.strip() for f in args.families.split(",")] if args.families
                else list(refs))
    for fam in families:
        if fam not in refs:
            sys.exit(f"ABORT: unknown family {fam!r}. Known: {', '.join(refs)}")

    for fam in families:
        rec = refs[fam]
        accs = rec.get("uniprot") or []
        if not accs:
            sys.exit(f"ABORT: {fam} has no UniProt accession in "
                     "reference_sequences.json; the family definition must come "
                     "from the reference entry, not from this script.")
        acc = accs[0]
        print(f"\n=== {fam}: reference {acc} ===", flush=True)

        ref_entry = entry(acc)
        ref_seq = ref_entry["sequence"]["value"]
        # The UniProt entry is used for ONE purpose: to define the Pfam family
        # query. The blastp query sequence is always the PDB seqres in
        # reference_sequences.json, so the structural R_i mapping is unaffected
        # by any difference here. A same-length, >=98%-identical entry is the
        # same gene in a different isolate and is accepted for the signature
        # role, with every substitution recorded in the provenance file. Any
        # indel, or identity below the floor, aborts: that is a different gene.
        VARIANT_ID_FLOOR = 0.98
        variant: list[dict] = []
        if ref_seq != rec["seq"]:
            if len(ref_seq) != len(rec["seq"]):
                sys.exit(f"ABORT: {fam}: UniProt {acc} length differs from "
                         f"reference_sequences.json ({len(ref_seq)} vs "
                         f"{len(rec['seq'])} aa). Resolve before retrieving.")
            variant = [{"pos": i + 1, "pdb": a, "uniprot": b}
                       for i, (a, b) in enumerate(zip(rec["seq"], ref_seq))
                       if a != b]
            ident = 1 - len(variant) / len(ref_seq)
            if ident < VARIANT_ID_FLOOR:
                sys.exit(f"ABORT: {fam}: UniProt {acc} is only {ident:.3f} "
                         f"identical to reference_sequences.json "
                         f"({len(variant)} substitutions). Resolve before "
                         f"retrieving.")
            print(f"  reference is a STRAIN VARIANT of {acc}: "
                  f"{len(variant)} substitutions, identity {ident:.4f} "
                  f"(accepted for Pfam-signature use; blastp query remains the "
                  f"PDB seqres). Substitutions: "
                  + ", ".join(f"{d['pdb']}{d['pos']}{d['uniprot']}"
                              for d in variant))
        else:
            print(f"  reference sequence matches reference_sequences.json "
                  f"({len(ref_seq)} aa)")

        sigs = pfam_signatures(ref_entry)
        if not sigs:
            sys.exit(f"ABORT: {fam}: no Pfam signature on {acc}; this route "
                     "cannot define a candidate set for it.")
        query = ("(" + " OR ".join(f"xref:pfam-{s}" for s in sigs) + ")"
                 f" AND taxonomy_id:{args.taxon}")
        print(f"  Pfam signatures on {acc}: {', '.join(sigs)}")
        print(f"  query: {query}", flush=True)

        rows, meta = stream_search(query)
        print(f"  retrieved {len(rows)} records "
              f"(UniProt {meta.get('release')}, header total "
              f"{meta.get('total_results')})")

        seen, fasta, kept = set(), [], []
        for r in rows:
            a, s = r.get("Entry", ""), (r.get("Sequence") or "").strip()
            if not a or not s or a in seen:
                continue
            seen.add(a)
            kept.append(r)
            fasta.append(f">{a} {r.get('Organism', '')}\n"
                         + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)))

        fasta_path = OUT / f"{fam}_candidates.fasta"
        meta_path = OUT / f"{fam}_candidates.tsv"
        prov_path = OUT / f"{fam}_provenance.json"

        write_atomic(fasta_path, "\n".join(fasta) + "\n")
        cols = list(kept[0].keys()) if kept else FIELDS
        write_atomic(meta_path, "\t".join(cols) + "\n"
                     + "\n".join("\t".join((r.get(c) or "").replace("\t", " ")
                                           for c in cols) for r in kept) + "\n")

        with_refseq = sum(1 for r in kept if (r.get("RefSeq") or "").strip())
        with_embl = sum(1 for r in kept if (r.get("EMBL") or "").strip())
        ref_present = acc in seen
        prov = {
            "family": fam, "retrieved_at": iso(),
            "route": "uniprotkb_pfam_family_scan",
            "reference": {"uniprot": acc, "length": len(ref_seq),
                          "pdb_vs_uniprot_substitutions": variant,
                          "pdb_uniprot_identity":
                              round(1 - len(variant) / len(ref_seq), 5),
                          "organism": ref_entry["organism"]["scientificName"],
                          "taxid": ref_entry["organism"]["taxonId"],
                          "pdb": rec.get("pdb"),
                          "pfam_signatures": sigs,
                          "sequence_sha256":
                              hashlib.sha256(ref_seq.encode()).hexdigest()},
            "query": query, "taxon_scope": args.taxon,
            "uniprot_release": meta.get("release"),
            "uniprot_release_date": meta.get("release_date"),
            "server_total_results": meta.get("total_results"),
            "pages_fetched": meta.get("pages"),
            "records_written": len(kept),
            "reference_in_candidate_set": ref_present,
            "with_refseq_xref": with_refseq,
            "with_embl_xref": with_embl,
            "fields": FIELDS,
            "files": {
                p.name: {"bytes": p.stat().st_size,
                         "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                for p in (fasta_path, meta_path)},
            "note": "retrieval only -- no identity, coverage, length or "
                    "dereplication filter applied here",
        }
        write_atomic(prov_path, json.dumps(prov, indent=2))

        print(f"  wrote {fasta_path.name}, {meta_path.name}, {prov_path.name}")
        print(f"  reference {acc} in candidate set : {ref_present}")
        print(f"  RefSeq cross-reference           : {with_refseq}/{len(kept)}")
        print(f"  EMBL cross-reference             : {with_embl}/{len(kept)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
