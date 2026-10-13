#!/usr/bin/env python3
"""AFP check 2: find and characterise candidate AFP sheath loci, positionally.

SCIENTIFIC QUESTION
-------------------
Do homologous AFP loci preserve distinguishable Afp2-like and Afp3-like
positions, despite containing several related sheath genes?

WHY SIMILARITY CANNOT MAKE THE ASSIGNMENT
-----------------------------------------
Reference Afp2 (Q6HAD7, 354 aa, slot 0) and Afp3 (Q6HAD6, 451 aa, slot +1) are
only 38.1% identical to one another, and the third tandem sheath gene of the
same locus (WP_169558217.1, 415 aa, slot +2) hits Afp3 at 140 bits and Afp2 at
136 bits — a tie. A best-hit label would therefore interchange the two
evolutionary identities silently in an unknown fraction of loci. This script
NEVER uses a sheath gene's own similarity to assign it.

HOW THE ASSIGNMENT IS MADE INSTEAD, DEFINED BEFORE ANY RESULT WAS SEEN
----------------------------------------------------------------------
  1. resolve each sampled candidate to a genomic locus and take its ordered CDS
     neighbourhood within +/-WINDOW, oriented so the anchor CDS reads + (the
     reference locus is also oriented + at its Afp2 anchor);
  2. blastp every neighbour against all 36 proteins of the reference pADAP
     locus; a neighbour is SHEATH-LIKE if its best reference hit is one of the
     three reference sheath genes (slots 0, +1, +2);
  3. anchor the locus using ONLY NON-SHEATH neighbours that hit a reference
     protein: each such pair gives a candidate_slot - reference_slot offset. The
     locus is anchored if >= MIN_ANCHORS non-sheath anchors agree on one modal
     offset and the anchored pairs are collinear (Spearman rho >= MIN_RHO);
  4. with the offset fixed, reference slots 0 / +1 / +2 map onto candidate
     slots. A sheath-like CDS at the mapped slot is labelled afp2_like /
     afp3_like / third_like BY POSITION ALONE;
  5. anything else is UNASSIGNED, with a recorded reason. Unassigned loci are
     output, not dropped.

Each sheath gene's own best-hit identity is recorded in a separate column so
check 3 can test concordance between position and sequence — but it plays no
part in the label.

Strictly diagnostic. No ortholog set is selected, no alignment QC is run, no
rates are estimated, and nothing frozen is modified.

Usage: python3 afp_loci.py [--n 120] [--window-kb 25]
Outputs: orthology/afp_loci.csv            one row per candidate locus
         orthology/afp_sheath_runs.csv     one row per sheath gene in a run
         orthology/afp_locus_tracks.json   every CDS with its reference mapping
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.resolve()))
import locus_architecture as la  # noqa: E402
import positional_orthology_r2pyocin as po  # noqa: E402
import synteny as syn_mod  # noqa: E402

OUT = HERE / "orthology"
CAND = HERE / "retrieval_uniprot"
# Family-parameterised. The DEFAULTS reproduce the AFP run exactly (three tandem
# reference sheath genes, shared Pfam pool, symmetric prescreen); --ref-track,
# --sheath-slots, --candidates and --prescreen point the same procedure at any
# other family. Nothing about the comparison itself is family-specific.
REF_TRACK_DEFAULT = "afp_reference_track.json"
CACHE_DEFAULT = "afp_locus_cache.json"
SHEATH_SLOTS_DEFAULT = "0:afp2_like,1:afp3_like,2:third_like"
EVALUE = 1e-5
MIN_ANCHORS = 3          # non-sheath flanking anchors needed to fix the offset
MIN_RHO = 0.8            # collinearity of anchored pairs
MIN_ANCHOR_FRAC = 0.6    # share of anchors agreeing on the modal offset


def ref_locus(path: Path) -> tuple[list[dict], dict[str, int], dict[int, str]]:
    track = json.loads(path.read_text())
    track = track["track"] if isinstance(track, dict) else track
    slot_of = {g["protein_id"]: g["slot"] for g in track if g["protein_id"]}
    by_slot = {g["slot"]: g["protein_id"] for g in track if g["protein_id"]}
    return track, slot_of, by_slot


def blast_best(query: dict[str, str], ref: dict[str, str]) -> dict[str, tuple[str, float, float]]:
    work = Path(tempfile.mkdtemp(prefix="afp_"))
    (work / "r.fa").write_text("".join(f">{k}\n{v}\n" for k, v in ref.items()))
    (work / "q.fa").write_text("".join(f">{k}\n{v}\n" for k, v in query.items()))
    subprocess.run(["makeblastdb", "-in", str(work / "r.fa"), "-dbtype", "prot",
                    "-out", str(work / "db")], check=True, capture_output=True)
    r = subprocess.run(["blastp", "-query", str(work / "q.fa"), "-db", str(work / "db"),
                        "-evalue", str(EVALUE), "-max_target_seqs", "10", "-seg", "no",
                        "-outfmt", "6 qseqid sseqid bitscore pident evalue"],
                       check=True, capture_output=True, text=True)
    best: dict[str, tuple[str, float, float]] = {}
    allhits: dict[str, list] = {}
    for line in r.stdout.strip().split("\n"):
        if not line.strip():
            continue
        q, s, b, pid, e = line.split("\t")
        allhits.setdefault(q, []).append((s, float(b), float(pid)))
        if q not in best or float(b) > best[q][1]:
            best[q] = (s, float(b), float(pid))
    return best, allhits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=120, help="candidate loci to sample")
    ap.add_argument("--window-kb", type=float, default=25.0)
    ap.add_argument("--genera", default=None,
                    help="restrict to these genera and dereplicate by SPECIES "
                         "instead of genus (for building one subfamily's set)")
    ap.add_argument("--out-prefix", default="afp")
    ap.add_argument("--ref-track", default=REF_TRACK_DEFAULT)
    ap.add_argument("--cache", default=CACHE_DEFAULT)
    ap.add_argument("--sheath-slots", default=SHEATH_SLOTS_DEFAULT,
                    help="reference slots holding sheath genes, as "
                         "slot:label[,slot:label...]")
    ap.add_argument("--candidates", default="afp_Afp2_candidates.tsv")
    ap.add_argument("--prescreen", default="afp_prescreen.csv",
                    help="CSV with a query-accession column and a bitscore "
                         "column; afp_prescreen.csv (q,bits) or a "
                         "*_blastp_hits.csv (uniprot_accession,bitscore)")
    ap.add_argument("--dedup", choices=("genus", "species"), default=None,
                    help="one representative per genus (default family-wide) "
                         "or per species (default with --genera)")
    ap.add_argument("--min-bits", type=float, default=200.0,
                    help="symmetric screen: bitscore to the best of the three "
                         "reference sheath gene(s)")
    args = ap.parse_args()
    window = int(args.window_kb * 1000)
    OUT.mkdir(parents=True, exist_ok=True)
    REF_SHEATH_SLOTS = {int(k): v for k, v in
                        (x.split(":") for x in args.sheath_slots.split(","))}
    CACHE = OUT / args.cache

    track, ref_slot_of, ref_by_slot = ref_locus(OUT / args.ref_track)
    ref_ids = [g["protein_id"] for g in track if g["protein_id"]]
    sheath_ref_ids = {ref_by_slot[s] for s in REF_SHEATH_SLOTS if s in ref_by_slot}
    print(f"reference locus: {len(ref_ids)} CDS, sheath genes at slots "
          f"{sorted(REF_SHEATH_SLOTS)} = {sorted(sheath_ref_ids)}")

    # ---- candidate sampling -------------------------------------------------
    # Enrichment only, and symmetric over the three reference sheath genes: a
    # candidate qualifies on its bitscore to the BEST of the three, so the screen
    # cannot favour an Afp2 or an Afp3 assignment. It decides where to look, not
    # what anything is called.
    cands = pd.read_csv(CAND / args.candidates, sep="\t", dtype=str).fillna("")
    cands["refseq_first"] = (cands.RefSeq.str.split(";").str[0]
                             .str.split("[").str[0].str.strip())
    cands = cands[cands.refseq_first != ""].copy()
    cands["genus"] = cands.Organism.str.split().str[0]
    sp = OUT / args.prescreen
    screen = pd.read_csv(sp if sp.exists() else CAND / args.prescreen)
    qcol = "q" if "q" in screen.columns else "uniprot_accession"
    bcol = "bits" if "bits" in screen.columns else "bitscore"
    screen = screen.groupby(qcol, as_index=False)[bcol].max()
    keep = set(screen[screen[bcol] >= args.min_bits][qcol])
    cands = cands[cands.Entry.isin(keep)]
    import re as _re
    cands["up"] = cands.Proteomes.map(
        lambda s_: (_re.findall(r"UP\d{9}", s_) or [None])[0])
    multi = set(cands.up.value_counts()[lambda v: v > 1].index) - {None}
    cands["multicopy"] = cands.up.isin(multi)
    print(f"candidate pool after the symmetric screen (bits >= {args.min_bits}): "
          f"{len(cands)} across {cands.genus.nunique()} genera; "
          f"{int(cands.multicopy.sum())} in proteomes with >1 detected copy")
    cands["species"] = cands.Organism.str.split().str[:2].str.join(" ")
    if args.genera:
        want = {g.strip() for g in args.genera.split(",")}
        cands = cands[cands.genus.isin(want)].copy()
        print(f"  restricted to {len(want)} genera: {len(cands)} members, "
              f"{cands.species.nunique()} species")
    dedup_on = args.dedup or ("species" if args.genera else "genus")
    sample = (cands.sort_values("multicopy", ascending=False)
                   .sample(frac=1.0, random_state=0)
                   .sort_values("multicopy", ascending=False)
                   .drop_duplicates(dedup_on).head(args.n))
    if len(sample) < args.n and not args.genera:
        rest = cands[~cands.Entry.isin(sample.Entry)]
        sample = pd.concat([sample, rest.head(args.n - len(sample))])
    print(f"sampled {len(sample)} candidates across {sample.genus.nunique()} genera")

    locs = syn_mod.ipg_locations(sample.refseq_first.tolist())
    print(f"  loci resolved {len(locs)}/{len(sample)}")

    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    nbs = {}
    for n, (_, r) in enumerate(sample.iterrows(), 1):
        a = r.refseq_first
        if a not in locs:
            continue
        key = f"{a}@{window}"
        if key in cache:
            nbs[r.Entry] = cache[key]
            continue
        L = locs[a]
        nb = la.neighbourhood(L["nuc"], L["start"], L["stop"], a, window=window)
        if nb:
            cache[key] = nb
            nbs[r.Entry] = nb
        if n % 20 == 0:
            print(f"    {n}/{len(sample)} neighbourhoods", flush=True)
            CACHE.write_text(json.dumps(cache))
    CACHE.write_text(json.dumps(cache))
    print(f"  neighbourhoods for {len(nbs)} candidates")

    all_ids = list(ref_ids)
    for nb in nbs.values():
        all_ids += [g["protein_id"] for g in nb["genes"] if g["protein_id"]]
    print(f"fetching {len(set(all_ids))} proteins")
    seqs = po.fetch_proteins(all_ids)
    ref_seqs = {i: seqs[i] for i in ref_ids if i in seqs}
    print(f"  reference proteins with sequence {len(ref_seqs)}/{len(ref_ids)}")

    loc_rows, run_rows, tracks = [], [], {}
    org_of = dict(zip(sample.Entry, sample.Organism))
    for acc, nb in nbs.items():
        qs = {g["protein_id"]: seqs[g["protein_id"]] for g in nb["genes"]
              if g["protein_id"] in seqs}
        if not qs:
            continue
        best, allhits = blast_best(qs, ref_seqs)
        genes = []
        for g in nb["genes"]:
            pid = g["protein_id"]
            hit = best.get(pid)
            rslot = ref_slot_of.get(hit[0]) if hit else None
            genes.append({
                "slot": g["slot"], "protein_id": pid, "product": g["product"],
                "rel_strand": g["rel_strand"], "start": g["start"], "stop": g["stop"],
                "ref_hit": hit[0] if hit else None, "ref_slot": rslot,
                "bitscore": hit[1] if hit else None,
                "pident": hit[2] if hit else None,
                "is_sheath_like": bool(hit and hit[0] in sheath_ref_ids),
            })
        tracks[acc] = genes
        df = pd.DataFrame(genes)

        # ---- anchor on NON-sheath flanking homology only -------------------
        anch = df[(~df.is_sheath_like) & df.ref_slot.notna()]
        offsets = (anch.slot - anch.ref_slot).astype(int)
        n_anch = len(anch)
        modal_off, modal_n, rho = None, 0, None
        if n_anch:
            c = Counter(offsets)
            modal_off, modal_n = c.most_common(1)[0]
            if n_anch >= 4 and anch.ref_slot.nunique() > 1:
                rho = float(spearmanr(anch.slot, anch.ref_slot).statistic)
        anchored = (n_anch >= MIN_ANCHORS and modal_n / max(n_anch, 1) >= MIN_ANCHOR_FRAC
                    and (rho is None or rho >= MIN_RHO))

        # ---- sheath runs ---------------------------------------------------
        sl = sorted(df[df.is_sheath_like].slot.tolist())
        runs, cur = [], []
        for s in sl:
            if cur and s == cur[-1] + 1:
                cur.append(s)
            else:
                if cur:
                    runs.append(cur)
                cur = [s]
        if cur:
            runs.append(cur)
        longest = max(runs, key=len) if runs else []

        # ---- positional labels, from the offset alone -----------------------
        expected = ({s + modal_off: lab for s, lab in REF_SHEATH_SLOTS.items()}
                    if anchored else {})
        for s in sl:
            row = df[df.slot == s].iloc[0]
            label = expected.get(s, "unassigned")
            reason = ""
            if not anchored:
                reason = (f"locus not anchored (non-sheath anchors={n_anch}, "
                          f"modal share={modal_n}/{n_anch}, rho={rho})")
            elif label == "unassigned":
                reason = "sheath gene outside the mapped reference slots"
            run_rows.append({
                "candidate": acc, "organism": org_of.get(acc, ""),
                "nuc": nb["nuc"], "protein_id": row.protein_id, "slot": int(s),
                "slot_in_run": (longest.index(s) if s in longest else None),
                "run_length": len(longest),
                "rel_strand": row.rel_strand, "length_aa": None,
                "positional_label": label, "reason": reason,
                "anchored": anchored, "offset": modal_off,
                # recorded for check 3, NOT used for the label
                "seq_best_hit": row.ref_hit,
                "seq_best_hit_label": REF_SHEATH_SLOTS.get(
                    ref_slot_of.get(row.ref_hit, 99), "other"),
                "seq_bitscore": row.bitscore, "seq_pident": row.pident,
            })

        gaps = []
        for a_, b_ in zip(longest, longest[1:]):
            ga = df[df.slot == a_].iloc[0]
            gb = df[df.slot == b_].iloc[0]
            gaps.append(int(gb.start - ga.stop))
        loc_rows.append({
            "candidate": acc, "organism": org_of.get(acc, ""), "nuc": nb["nuc"],
            "n_cds": len(df), "n_sheath_like": len(sl),
            "n_runs": len(runs), "longest_run": len(longest),
            "run_slots": ";".join(map(str, longest)),
            "run_strands": ";".join(df[df.slot.isin(longest)].rel_strand),
            "run_intergenic_bp": ";".join(map(str, gaps)),
            "non_sheath_anchors": n_anch, "modal_offset": modal_off,
            "modal_anchor_share": round(modal_n / n_anch, 3) if n_anch else None,
            "collinearity_rho": None if rho is None else round(rho, 3),
            "anchored": anchored,
            "flanking_ref_slots": ";".join(
                map(str, sorted(anch.ref_slot.dropna().astype(int).unique()))),
        })

    loci = pd.DataFrame(loc_rows).sort_values(["anchored", "longest_run"],
                                              ascending=[False, False])
    runs_t = pd.DataFrame(run_rows)
    for a, s in seqs.items():
        runs_t.loc[runs_t.protein_id == a, "length_aa"] = len(s)
    loci.to_csv(OUT / f"{args.out_prefix}_loci.csv", index=False)
    runs_t.to_csv(OUT / f"{args.out_prefix}_sheath_runs.csv", index=False)
    (OUT / f"{args.out_prefix}_locus_tracks.json").write_text(json.dumps(tracks))

    print(f"\nloci characterised: {len(loci)}")
    print("  anchored on non-sheath flanking homology: "
          f"{int(loci.anchored.sum())} / {len(loci)}")
    print("  sheath-run length distribution: "
          f"{loci.longest_run.value_counts().sort_index().to_dict()}")
    print("  positional labels: "
          f"{runs_t.positional_label.value_counts().to_dict()}")
    if len(runs_t):
        labs = set(REF_SHEATH_SLOTS.values())
        both = runs_t[runs_t.positional_label.isin(labs)]
        agree = int((both.positional_label == both.seq_best_hit_label).sum())
        print(f"  position vs sequence best hit (recorded only): "
              f"{agree}/{len(both)} agree")
    print(f"\nwrote {OUT}/{args.out_prefix}_" + "{loci,sheath_runs,locus_tracks}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
