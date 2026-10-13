#!/usr/bin/env python3
"""Step 4 of the frozen protocol: genomic-neighbourhood verification.

Sequence similarity alone cannot distinguish a positional ortholog from a
paralogous sheath elsewhere in the genome, another contractile system, or a
free-standing prophage tail gene. That matters most exactly where these systems
live -- prophage-like loci carrying several related tail genes.

For each candidate the genomic location is resolved through NCBI Identical
Protein Groups, the flanking +/-15 kb is read as a feature table, and the
products are scored against contractile-machine marker categories.

This stage SCORES locus context; it does not by itself exclude anything.

A fixed annotation count is not a biological definition of orthology. These loci
differ in architecture, gene length and annotation quality across AFP, pyocin,
SalCIS and diffocin, and a genuine ortholog can score zero simply because its
neighbours are annotated "hypothetical protein". The score is therefore carried
forward as evidence and combined with gene-tree placement downstream; the
distribution is reported at the checkpoint before any threshold becomes an
exclusion rule.

  window                +/- 15 kb around the candidate CDS
  marker categories     tube, baseplate, tape measure, spike/PAAR/VgrG,
                        terminator, fibre, ATPase/ClpV, generic phage tail
  n_categories          distinct marker categories found, excluding sheath itself
  synteny_supported     n_categories >= 2  (a LABEL, not a filter)
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
IN = HERE / "filtered"
OUT = HERE / "synteny"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

WINDOW = 15000
MIN_CATEGORIES = 2
MAX_PER_FAMILY = 250
PAUSE = 0.4           # <= 3 requests/s without an API key

MARKERS = {
    "tube": r"tail tube|tube protein|\bhcp\b|tail tube protein",
    "baseplate": r"baseplate|base plate",
    "tape_measure": r"tape measure|tape-measure",
    "spike": r"\bpaar\b|\bvgrg\b|spike|puncturing",
    "terminator": r"tail terminator|terminator protein",
    "fibre": r"tail fiber|tail fibre|receptor.binding",
    "atpase": r"\bclpv\b|aaa family atpase|aaa\+ atpase",
    "phage_tail": r"phage tail|tail protein|tail assembly|tail length",
    "sheath": r"tail sheath|sheath protein",
}
COMPILED = {k: re.compile(v, re.I) for k, v in MARKERS.items()}


def get(url, params, tries=4):
    q = f"{url}?{urllib.parse.urlencode(params)}"
    for a in range(tries):
        try:
            with urllib.request.urlopen(q, timeout=120) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as exc:
            if a == tries - 1:
                return ""
            time.sleep(3 + 3 * a)
    return ""


def ipg_locations(accessions):
    """Resolve protein accessions to one genomic location each (batched)."""
    loc = {}
    for i in range(0, len(accessions), 50):
        chunk = accessions[i:i + 50]
        txt = get(f"{EUTILS}/efetch.fcgi",
                  {"db": "protein", "id": ",".join(chunk),
                   "rettype": "ipg", "retmode": "text"})
        time.sleep(PAUSE)
        for line in txt.splitlines()[1:]:
            f = line.split("\t")
            if len(f) < 11:
                continue
            nuc, start, stop, strand, prot = f[2], f[3], f[4], f[5], f[6]
            if not nuc or not start or prot not in chunk or prot in loc:
                continue
            try:
                loc[prot] = {"nuc": nuc, "start": int(start), "stop": int(stop),
                             "strand": strand, "assembly": f[10] if len(f) > 10 else ""}
            except ValueError:
                continue
    return loc


def neighbourhood(nuc, start, stop):
    lo = max(1, min(start, stop) - WINDOW)
    hi = max(start, stop) + WINDOW
    txt = get(f"{EUTILS}/efetch.fcgi",
              {"db": "nuccore", "id": nuc, "seq_start": lo, "seq_stop": hi,
               "rettype": "ft", "retmode": "text"})
    time.sleep(PAUSE)
    return [m.group(1).strip() for m in re.finditer(r"\n\t{3}product\t(.+)", txt)]


def score(products):
    cats = {k for k, rx in COMPILED.items()
            for p in products if rx.search(p)}
    return cats


def select(df, cap):
    """Spread the synteny budget across genera and identity deciles rather than
    taking the top hits, which would all be near neighbours of the reference."""
    if len(df) <= cap:
        return df
    df = df.copy()
    df["genus"] = df.species.str.split().str[0]
    df["bin"] = pd.qcut(df.pident, 10, labels=False, duplicates="drop")
    out, seen = [], defaultdict(int)
    for _, row in df.sample(frac=1.0, random_state=0).iterrows():
        key = (row["genus"], row["bin"])
        if seen[key] < max(1, cap // (df.genus.nunique() or 1)):
            out.append(row)
            seen[key] += 1
        if len(out) >= cap:
            break
    if len(out) < cap:
        rest = df.drop([r.name for r in out], errors="ignore")
        out += [r for _, r in rest.head(cap - len(out)).iterrows()]
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    files = sorted(IN.glob("*_filtered.csv"))
    if not files:
        print("no filtered input; run filter_hits.py first")
        return

    for path in files:
        fam = path.stem.replace("_filtered", "")
        df = pd.read_csv(path)
        sel = select(df, MAX_PER_FAMILY).reset_index(drop=True)
        print(f"\n{fam}: {len(df)} candidates -> {len(sel)} sent to synteny", flush=True)

        locs = ipg_locations(sel.accession.tolist())
        print(f"  genomic location resolved for {len(locs)}/{len(sel)}", flush=True)

        rows = []
        for n, (_, r) in enumerate(sel.iterrows(), 1):
            acc = r.accession
            L = locs.get(acc)
            if not L:
                rows.append({**r.to_dict(), "nuc": None, "n_categories": 0,
                             "categories": "", "synteny_supported": False,
                             "reason": "no genomic location"})
                continue
            prods = neighbourhood(L["nuc"], L["start"], L["stop"])
            cats = score(prods)
            other = cats - {"sheath"}
            rows.append({**r.to_dict(), "nuc": L["nuc"], "assembly": L["assembly"],
                         "n_neighbour_cds": len(prods),
                         "n_categories": len(other),
                         "categories": ";".join(sorted(other)),
                         "synteny_supported": len(other) >= MIN_CATEGORIES,
                         "reason": "" if len(other) >= MIN_CATEGORIES
                                   else "sparse/uninformative annotation"})
            if n % 25 == 0:
                print(f"    {n}/{len(sel)} neighbourhoods read", flush=True)

        out = pd.DataFrame(rows)
        out.to_csv(OUT / f"{fam}_synteny.csv", index=False)
        dist = out.n_categories.value_counts().sort_index().to_dict()
        print(f"  {fam}: n_categories distribution {dist}  "
              f"(>=2 in {int(out.synteny_supported.sum())}/{len(out)})", flush=True)


if __name__ == "__main__":
    main()
