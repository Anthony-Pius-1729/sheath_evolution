#!/usr/bin/env python3
"""Ordered gene-neighbourhood extraction and comparison — positional orthology.

The frozen synteny stage asks WHICH marker categories occur within +/-15 kb. That
question is answered identically by every clade in the R2 retrieval set (363/363
resolved loci score >=2 categories), so it cannot tell a positional ortholog from
any other contractile-tail gene. This module asks the next question instead:
in WHAT ORDER, at what spacing, and on which strand.

It adds no scoring to the frozen protocol. `synteny.py` and its MARKERS are
imported unchanged and used only to label genes; the ordered comparison defined
here is a new, separate discriminator reported alongside, never substituted for,
the frozen category score.

COMPARISON DEFINED BEFORE ANY RESULT WAS LOOKED AT
--------------------------------------------------
For a locus:
  * take every CDS whose midpoint lies within +/-WINDOW of the anchor sheath CDS;
  * label each with its frozen marker category, or `other` / `hypothetical`;
  * orient the neighbourhood so the anchor reads on the + strand, so that two
    loci are comparable regardless of which strand each was deposited on;
  * index genes by SLOT: signed rank relative to the anchor (anchor = 0,
    next gene downstream = +1, upstream = -1).

Two loci are then compared by:
  * `slot_agreement`  -- over slots -N..+N, the fraction where both loci have a
                         gene and their categories match. Order-sensitive; the
                         frozen category score is not.
  * `lcs_categories`  -- length of the longest common subsequence of the two
                         ordered category strings, normalised by the shorter.
                         Tolerates insertions/deletions that slot agreement
                         punishes.
  * `landmark_shared` -- flanking non-tail genes (products outside every marker
                         category) shared between the two loci, by product text.
                         Landmarks come from the reference locus as retrieved;
                         none are hardcoded.

No threshold on these quantities is defined here, and none is used to include or
exclude a sequence. They are reported per locus and per pair.

Usage: imported by positional_orthology_r2pyocin.py
"""

from __future__ import annotations

import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.resolve()))
import synteny as syn_mod  # noqa: E402 -- frozen MARKERS/COMPILED/WINDOW/PAUSE

EUTILS = syn_mod.EUTILS
WINDOW = syn_mod.WINDOW          # +/-15 kb, frozen
SLOTS = 8                        # slots either side of the anchor to compare
PAUSE = syn_mod.PAUSE


def _get(url: str, params: dict, tries: int = 4) -> str:
    q = f"{url}?{urllib.parse.urlencode(params)}"
    for a in range(tries):
        try:
            with urllib.request.urlopen(q, timeout=120) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception:
            if a == tries - 1:
                return ""
            time.sleep(3 + 3 * a)
    return ""


def feature_table(nuc: str, lo: int, hi: int) -> list[dict]:
    """Ordered CDS list from an NCBI feature table, with strand and product."""
    txt = _get(f"{EUTILS}/efetch.fcgi",
               {"db": "nuccore", "id": nuc, "seq_start": lo, "seq_stop": hi,
                "rettype": "ft", "retmode": "text"})
    time.sleep(PAUSE)
    feats, cur = [], None
    for line in txt.splitlines():
        m = re.match(r"^(<?\d+)\t(>?\d+)\t(\w+)", line)
        if m:
            if cur and cur["key"] == "CDS":
                feats.append(cur)
            a = int(m.group(1).lstrip("<>"))
            b = int(m.group(2).lstrip("<>"))
            cur = {"key": m.group(3), "start": min(a, b), "stop": max(a, b),
                   "strand": "+" if a <= b else "-", "product": "",
                   "protein_id": "", "gene": "", "locus_tag": ""}
            continue
        q = re.match(r"^\t\t\t(\w+)\t(.*)$", line)
        if q and cur:
            k, v = q.group(1), q.group(2)
            if k in ("product", "gene", "locus_tag"):
                cur[k] = cur[k] or v
            elif k == "protein_id":
                cur["protein_id"] = v.replace("ref|", "").strip("|")
    if cur and cur["key"] == "CDS":
        feats.append(cur)
    feats.sort(key=lambda f: (f["start"] + f["stop"]) / 2)
    return feats


def categorise(product: str) -> str:
    for cat, rx in syn_mod.COMPILED.items():
        if rx.search(product):
            return cat
    if re.search(r"hypothetical|uncharacteri|DUF", product, re.I):
        return "hypothetical"
    return "other"


def neighbourhood(nuc: str, start: int, stop: int, anchor_protein: str | None,
                  window: int = WINDOW) -> dict | None:
    """Ordered, strand-normalised neighbourhood around the anchor CDS."""
    lo = max(1, min(start, stop) - window)
    hi = max(start, stop) + window
    feats = feature_table(nuc, lo, hi)
    if not feats:
        return None
    mid = (min(start, stop) + max(start, stop)) / 2 - lo + 1
    # the anchor is the CDS carrying the query protein id, else the one whose
    # midpoint is closest to the queried coordinates
    idx = None
    if anchor_protein:
        for i, f in enumerate(feats):
            if f["protein_id"].split(".")[0] == anchor_protein.split(".")[0]:
                idx = i
                break
    if idx is None:
        idx = min(range(len(feats)),
                  key=lambda i: abs((feats[i]["start"] + feats[i]["stop"]) / 2 - mid))
    for f in feats:
        f["category"] = categorise(f["product"])
    flip = feats[idx]["strand"] == "-"
    ordered = list(reversed(feats)) if flip else feats
    aidx = len(feats) - 1 - idx if flip else idx
    out = []
    for i, f in enumerate(ordered):
        out.append({**f, "slot": i - aidx,
                    "rel_strand": ("+" if (f["strand"] == "+") != flip else "-")})
    return {"nuc": nuc, "anchor_slot": 0, "flipped": flip, "genes": out,
            "anchor_product": ordered[aidx]["product"],
            "anchor_protein_id": ordered[aidx]["protein_id"]}


def slot_categories(nb: dict, slots: int = SLOTS) -> dict[int, str]:
    return {g["slot"]: g["category"] for g in nb["genes"]
            if -slots <= g["slot"] <= slots}


def slot_agreement(a: dict, b: dict, slots: int = SLOTS) -> tuple[float, int, int]:
    ca, cb = slot_categories(a, slots), slot_categories(b, slots)
    shared = [s for s in range(-slots, slots + 1) if s in ca and s in cb and s != 0]
    if not shared:
        return 0.0, 0, 0
    match = sum(1 for s in shared if ca[s] == cb[s])
    return match / len(shared), match, len(shared)


def lcs_categories(a: dict, b: dict, slots: int = SLOTS) -> float:
    sa = [a_["category"] for a_ in a["genes"] if abs(a_["slot"]) <= slots]
    sb = [b_["category"] for b_ in b["genes"] if abs(b_["slot"]) <= slots]
    if not sa or not sb:
        return 0.0
    dp = [[0] * (len(sb) + 1) for _ in range(len(sa) + 1)]
    for i, x in enumerate(sa, 1):
        for j, y in enumerate(sb, 1):
            dp[i][j] = dp[i - 1][j - 1] + 1 if x == y else max(dp[i - 1][j], dp[i][j - 1])
    return dp[-1][-1] / min(len(sa), len(sb))


def landmarks(nb: dict, slots: int = SLOTS) -> set[str]:
    """Flanking genes that are NOT contractile-tail markers, by product text."""
    out = set()
    for g in nb["genes"]:
        if abs(g["slot"]) > slots or g["category"] in ("other", "hypothetical"):
            if g["category"] in ("other",) and g["product"]:
                out.add(g["product"].lower().strip())
        elif g["category"] == "other" and g["product"]:
            out.add(g["product"].lower().strip())
    return out
