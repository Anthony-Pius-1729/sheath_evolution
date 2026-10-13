#!/usr/bin/env python3
"""Positional orthology for R2 pyocin, by gene-order homology.

THE PROBLEM THIS ADDRESSES
--------------------------
39 of the 57 same-species paralogs that dereplication discarded fall inside the
candidate ortholog clade, interleaved and not monophyletic. Sequence similarity
chose between those copies by a margin under 20 bits in over half the affected
species. Domain architecture, protein names and the frozen marker-category score
cannot separate them: every resolved locus in the whole retrieval set scores >=2
categories.

WHY PRODUCT TEXT CANNOT BE THE DISCRIMINATOR
--------------------------------------------
The R2 reference locus (NP_249313.1, NC_002516.2 ~681 kb) is annotated almost
entirely as "bacteriophage protein" and "hypothetical protein" -- it carries NOT
ONE of the frozen marker-category names. Scored by the frozen rule the true
reference locus would read 0 categories, "sparse/uninformative annotation". Any
comparison based on product vocabulary would therefore score the reference
against itself as a non-match. Annotation vocabulary is a property of who
deposited the genome, not of the biology.

WHAT IS COMPARED INSTEAD
------------------------
Gene-order conservation by SEQUENCE HOMOLOGY of the neighbours, which is the
standard positional-orthology criterion and is immune to annotation vocabulary:

  1. take the ordered CDS neighbourhood of each locus within +/-WINDOW_KB,
     oriented so the sheath anchor reads on the + strand, slots signed from it;
  2. blastp every neighbour against the reference locus's neighbours;
  3. a query gene is ASSIGNED to the reference slot of its best hit (E <= 1e-3);
  4. report, per locus:
       n_homologous      genes with any reference counterpart
       n_slots_hit       distinct reference slots recovered
       exact_slot_match  genes whose own slot equals its counterpart's slot
       collinearity      Spearman rho between query slot and reference slot
       regulators        homologs of the reference's PrtN / PrtR / PtrB
                         regulatory cassette, named FROM the reference locus
  5. for each (clade member, discarded paralog) pair, report both copies'
     values side by side.

No threshold is applied and nothing is excluded by this script. It asks whether
gene order separates the copies at all; the decision is made in the report.

Usage: python3 positional_orthology_r2pyocin.py [--pairs N] [--window-kb 25]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.resolve()))
import locus_architecture as la  # noqa: E402
import synteny as syn_mod  # noqa: E402

OUT = HERE / "orthology"
CACHE = OUT / "locus_cache.json"
REF_PROTEIN = "NP_249313.1"       # Q9S574 in the PAO1 genome
EVALUE = 1e-3


def fetch_proteins(ids: list[str]) -> dict[str, str]:
    seqs, ids = {}, [i for i in dict.fromkeys(ids) if i]
    for i in range(0, len(ids), 150):
        chunk = ids[i:i + 150]
        q = (f"{syn_mod.EUTILS}/efetch.fcgi?"
             + urllib.parse.urlencode({"db": "protein", "id": ",".join(chunk),
                                       "rettype": "fasta", "retmode": "text"}))
        for attempt in range(4):
            try:
                with urllib.request.urlopen(q, timeout=180) as fh:
                    txt = fh.read().decode("utf-8", "replace")
                break
            except Exception:
                if attempt == 3:
                    txt = ""
                time.sleep(3 + 3 * attempt)
        name = None
        for line in txt.splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
        time.sleep(syn_mod.PAUSE)
        print(f"    proteins {min(i + 150, len(ids))}/{len(ids)}", flush=True)
    return {k: "".join(v) for k, v in seqs.items() if v}


def blast_map(query: dict[str, str], ref: dict[str, str]) -> dict[str, tuple[str, float]]:
    """best reference hit per query protein -> (ref_protein_id, bitscore)"""
    if not query or not ref:
        return {}
    work = Path(tempfile.mkdtemp(prefix="posorth_"))
    (work / "ref.fa").write_text("".join(f">{k}\n{v}\n" for k, v in ref.items()))
    (work / "q.fa").write_text("".join(f">{k}\n{v}\n" for k, v in query.items()))
    subprocess.run(["makeblastdb", "-in", str(work / "ref.fa"), "-dbtype", "prot",
                    "-out", str(work / "db")], check=True, capture_output=True)
    r = subprocess.run(["blastp", "-query", str(work / "q.fa"), "-db", str(work / "db"),
                        "-evalue", str(EVALUE), "-max_target_seqs", "5", "-seg", "no",
                        "-outfmt", "6 qseqid sseqid bitscore evalue"],
                       check=True, capture_output=True, text=True)
    best: dict[str, tuple[str, float]] = {}
    for line in r.stdout.strip().split("\n"):
        if not line.strip():
            continue
        q, s, b, _ = line.split("\t")
        b = float(b)
        if q not in best or b > best[q][1]:
            best[q] = (s, b)
    return best


def score_locus(nb: dict, ref_nb: dict, ref_slot_of: dict[str, int],
                best: dict[str, tuple[str, float]], reg_ids: dict[str, str]) -> dict:
    rows = []
    for g in nb["genes"]:
        pid = g["protein_id"]
        hit = best.get(pid)
        rslot = ref_slot_of.get(hit[0]) if hit else None
        rows.append({"slot": g["slot"], "protein_id": pid, "product": g["product"],
                     "category": g["category"], "rel_strand": g["rel_strand"],
                     "ref_hit": hit[0] if hit else None,
                     "ref_slot": rslot, "bitscore": hit[1] if hit else None})
    df = pd.DataFrame(rows)
    hom = df[df.ref_slot.notna()]
    exact = int((hom.slot == hom.ref_slot).sum())
    rho = (float(spearmanr(hom.slot, hom.ref_slot).statistic)
           if len(hom) >= 4 and hom.ref_slot.nunique() > 1 else None)
    found = {name: (pid in set(hom.ref_hit)) for name, pid in reg_ids.items()}
    return {"table": df,
            "n_genes": len(df), "n_homologous": len(hom),
            "n_slots_hit": int(hom.ref_slot.nunique()),
            "exact_slot_match": exact,
            "collinearity_rho": None if rho is None else round(rho, 3),
            "anchor_homologous": bool(df[df.slot == 0].ref_slot.notna().any()),
            **{f"has_{k}": v for k, v in found.items()},
            "n_regulators": int(sum(found.values()))}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=0, help="limit pairs (0 = all)")
    ap.add_argument("--window-kb", type=float, default=25.0)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    window = int(args.window_kb * 1000)

    pairs = pd.read_csv(OUT / "r2pyocin_L2_inparalog_pairs.csv")
    if args.pairs:
        pairs = pairs.head(args.pairs)
    hits = pd.read_csv(HERE / "retrieval_uniprot" / "r2pyocin_blastp_hits.csv")
    acc_of = dict(zip(hits.uniprot_accession, hits.refseq_first))
    print(f"{len(pairs)} in-paralog pairs; window +/-{args.window_kb:g} kb")

    # ---- reference locus ---------------------------------------------------
    rloc = syn_mod.ipg_locations([REF_PROTEIN])[REF_PROTEIN]
    ref_nb = la.neighbourhood(rloc["nuc"], rloc["start"], rloc["stop"],
                              REF_PROTEIN, window=window)
    ref_ids = [g["protein_id"] for g in ref_nb["genes"] if g["protein_id"]]
    ref_slot_of = {g["protein_id"]: g["slot"] for g in ref_nb["genes"]}
    print(f"reference locus {rloc['nuc']}: {len(ref_nb['genes'])} CDS, "
          f"anchor {ref_nb['anchor_protein_id']}")
    reg_ids = {}
    for g in ref_nb["genes"]:
        p = g["product"].lower()
        for key, tag in (("prtn", "PrtN"), ("prtr", "PrtR"), ("ptrb", "PtrB")):
            if key in p.replace("-", "").replace(" ", ""):
                reg_ids[tag] = g["protein_id"]
    print(f"  regulatory cassette from the reference: "
          + ", ".join(f"{k}={v}@slot{ref_slot_of[v]:+d}" for k, v in reg_ids.items()))

    # ---- resolve every locus ----------------------------------------------
    members = pairs.clade_member.tolist()
    paras = pairs.paralog.tolist()
    want = {u: acc_of.get(u) for u in members + paras}
    missing = [u for u, a in want.items() if not isinstance(a, str) or not a]
    accs = [a for a in want.values() if isinstance(a, str) and a]
    print(f"resolving {len(accs)} loci ({len(missing)} without a RefSeq accession)")
    locs = syn_mod.ipg_locations(accs)
    print(f"  resolved {len(locs)}/{len(accs)}")

    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    nbs = {}
    for n, (u, a) in enumerate(want.items(), 1):
        if not a or a not in locs:
            continue
        key = f"{a}@{window}"
        if key in cache:
            nbs[u] = cache[key]
            continue
        L = locs[a]
        nb = la.neighbourhood(L["nuc"], L["start"], L["stop"], a, window=window)
        if nb:
            cache[key] = nb
            nbs[u] = nb
        if n % 10 == 0:
            print(f"    {n}/{len(want)} neighbourhoods", flush=True)
            CACHE.write_text(json.dumps(cache))
    CACHE.write_text(json.dumps(cache))
    print(f"  neighbourhoods for {len(nbs)}/{len(want)} copies")

    # ---- protein sequences and blast --------------------------------------
    all_ids = list(ref_ids)
    for nb in nbs.values():
        all_ids += [g["protein_id"] for g in nb["genes"] if g["protein_id"]]
    print(f"fetching {len(set(all_ids))} neighbour proteins")
    seqs = fetch_proteins(all_ids)
    ref_seqs = {i: seqs[i] for i in ref_ids if i in seqs}
    print(f"  reference neighbours with sequence: {len(ref_seqs)}/{len(ref_ids)}")

    rows, tables = [], {}
    for u, nb in nbs.items():
        qs = {g["protein_id"]: seqs[g["protein_id"]] for g in nb["genes"]
              if g["protein_id"] in seqs}
        best = blast_map(qs, ref_seqs)
        sc = score_locus(nb, ref_nb, ref_slot_of, best, reg_ids)
        tables[u] = sc.pop("table")
        rows.append({"uniprot": u, "refseq": want[u], "nuc": nb["nuc"],
                     "role": "clade_member" if u in set(members) else "paralog",
                     **sc})
    per_locus = pd.DataFrame(rows)
    per_locus.to_csv(OUT / "r2pyocin_positional_loci.csv", index=False)

    # ---- pair comparison ---------------------------------------------------
    idx = per_locus.set_index("uniprot")
    prows = []
    for _, p in pairs.iterrows():
        m, g = p.clade_member, p.paralog
        if m not in idx.index or g not in idx.index:
            prows.append({"species": p.species, "clade_member": m, "paralog": g,
                          "resolved": False})
            continue
        M, G = idx.loc[m], idx.loc[g]
        prows.append({
            "species": p.species, "clade_member": m, "paralog": g, "resolved": True,
            "member_slots_hit": M.n_slots_hit, "paralog_slots_hit": G.n_slots_hit,
            "member_exact": M.exact_slot_match, "paralog_exact": G.exact_slot_match,
            "member_rho": M.collinearity_rho, "paralog_rho": G.collinearity_rho,
            "member_regulators": M.n_regulators, "paralog_regulators": G.n_regulators,
            "delta_slots_hit": M.n_slots_hit - G.n_slots_hit,
            "delta_exact": M.exact_slot_match - G.exact_slot_match,
            "delta_regulators": M.n_regulators - G.n_regulators,
            "bitscore_margin": p.get("patristic_to_member"),
        })
    pair_tab = pd.DataFrame(prows)
    pair_tab.to_csv(OUT / "r2pyocin_positional_pairs.csv", index=False)

    (OUT / "r2pyocin_positional_tracks.json").write_text(
        json.dumps({u: t.to_dict("records") for u, t in tables.items()}))

    r = pair_tab[pair_tab.resolved]
    print(f"\npairs resolved: {len(r)}/{len(pair_tab)}")
    if len(r):
        print(f"  member has MORE reference slots : {int((r.delta_slots_hit > 0).sum())}")
        print(f"  paralog has more                : {int((r.delta_slots_hit < 0).sum())}")
        print(f"  tied                            : {int((r.delta_slots_hit == 0).sum())}")
        print(f"  member has more regulators      : {int((r.delta_regulators > 0).sum())}")
        print(f"  paralog has more regulators     : {int((r.delta_regulators < 0).sum())}")
        print(f"  median slots hit  member {r.member_slots_hit.median():.1f} "
              f"vs paralog {r.paralog_slots_hit.median():.1f}")
        print(f"  median exact-slot member {r.member_exact.median():.1f} "
              f"vs paralog {r.paralog_exact.median():.1f}")
    print(f"\nwrote {OUT}/r2pyocin_positional_" + "{loci,pairs}.csv, _tracks.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
