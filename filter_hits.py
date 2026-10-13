#!/usr/bin/env python3
"""Steps 2-3 of the frozen protocol: identity window and architecture gates.

Pre-registered, applied identically to every family:
  initial identity window   30-90% to the structural reference
  reference coverage        >= 80% of the reference length aligned
  length compatibility      subject length within 0.75-1.30x the reference

The identity window is a RETRIEVAL window, not an inclusion rule. The hard gates
downstream are alignment- and synteny-based; this stage only removes fragments,
extreme length variants and sequences too close or too far to be informative.

One representative per organism is kept (best bitscore), both to dereplicate the
heavy strain-level redundancy in RefSeq and to keep the synteny stage tractable.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd

HERE = Path(__file__).parent
RET = HERE / "retrieval"
OUT = HERE / "filtered"

ID_LO, ID_HI = 0.30, 0.90
MIN_REF_COV = 0.80
LEN_LO, LEN_HI = 0.75, 1.30

ORG = re.compile(r"\[([^\[\]]+)\]\s*$")


def parse(path, ref_len):
    root = ET.parse(path).getroot()
    rows = []
    for hit in root.iter("Hit"):
        acc = hit.findtext("Hit_accession", "")
        hdef = hit.findtext("Hit_def", "")
        hlen = int(hit.findtext("Hit_len", "0") or 0)
        best = None
        for hsp in hit.iter("Hsp"):
            bits = float(hsp.findtext("Hsp_bit-score", "0"))
            if best is None or bits > best[0]:
                ident = int(hsp.findtext("Hsp_identity", "0"))
                alen = int(hsp.findtext("Hsp_align-len", "1"))
                qfrom = int(hsp.findtext("Hsp_query-from", "0"))
                qto = int(hsp.findtext("Hsp_query-to", "0"))
                best = (bits, ident / max(alen, 1), (qto - qfrom + 1) / ref_len,
                        float(hsp.findtext("Hsp_evalue", "1")))
        if best is None:
            continue
        m = ORG.search(hdef)
        org = m.group(1) if m else "?"
        # strip strain detail to a species-level key for dereplication
        sp = " ".join(org.replace("Candidatus ", "").split()[:2])
        rows.append({"accession": acc, "title": hdef[:110], "organism": org,
                     "species": sp, "slen": hlen, "bitscore": best[0],
                     "pident": best[1], "ref_cov": best[2], "evalue": best[3]})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    refs = json.load(open(HERE / "reference_sequences.json"))
    summary = []

    for name, rec in refs.items():
        path = RET / f"{name}.xml"
        if not path.exists():
            print(f"{name}: no retrieval file yet")
            continue
        ref_len = len(rec["seq"])
        df = parse(path, ref_len)
        n0 = len(df)

        df["len_ratio"] = df["slen"] / ref_len
        keep = df[(df.pident >= ID_LO) & (df.pident <= ID_HI)
                  & (df.ref_cov >= MIN_REF_COV)
                  & (df.len_ratio >= LEN_LO) & (df.len_ratio <= LEN_HI)].copy()
        n1 = len(keep)

        keep = (keep.sort_values("bitscore", ascending=False)
                    .drop_duplicates("species").reset_index(drop=True))
        keep["family"] = name
        keep.to_csv(OUT / f"{name}_filtered.csv", index=False)

        summary.append({
            "family": name, "ref_len": ref_len, "raw_hits": n0,
            "pass_gates": n1, "species_dereplicated": len(keep),
            "pident_median": round(keep.pident.median(), 3) if len(keep) else None,
            "pident_min": round(keep.pident.min(), 3) if len(keep) else None,
            "pident_max": round(keep.pident.max(), 3) if len(keep) else None,
            "genera": keep.species.str.split().str[0].nunique() if len(keep) else 0,
        })
        print(f"{name:12s} raw={n0:5d} -> gates={n1:5d} -> species={len(keep):4d} "
              f"(pident {keep.pident.min():.2f}-{keep.pident.max():.2f} "
              f"median {keep.pident.median():.2f})" if len(keep) else
              f"{name:12s} raw={n0} -> NOTHING PASSED")

    if summary:
        s = pd.DataFrame(summary)
        s.to_csv(OUT / "filter_summary.csv", index=False)
        print()
        print(s.to_string(index=False))


if __name__ == "__main__":
    main()
