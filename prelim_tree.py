#!/usr/bin/env python3
"""Preliminary gene tree for a single family -- retrieval validation only.

This is the Stage A/B QC screen of the frozen protocol, run for one family to
show that the retrieved candidates form a coherent clade around the structural
reference and that locus context is informative within it. It is NOT the
orthology decision and NOT a final phylogeny:

  * FastTree LG+Gamma, as the protocol specifies for QC screens; model selection
    and rigorous ML inference remain open and are deliberately not touched here.
  * No orthology cut is applied. The frozen rule is tree placement COMBINED WITH
    locus context, four-way, and the clade boundary is a decision to be made on
    the evidence -- so this script REPORTS the nested-clade profile from the
    reference to the root and the four-way table at each level, and fixes
    nothing. Any table it prints below a chosen level is provisional.
  * No rates. No Rate4Site.

Inputs : filtered/{fam}_filtered.csv, synteny/{fam}_synteny.csv (optional),
         retrieval_uniprot/{fam}_candidates.fasta, reference_sequences.json
Outputs: evolution_prelim/{fam}_prelim.afa, .nwk, _clade_profile.csv,
         _leaf_table.csv, _summary.json

Usage: python3 prelim_tree.py --family r2pyocin
Needs: pyfamsa, dendropy, pandas, tools/FastTree
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import dendropy
import pandas as pd
from pyfamsa import Aligner, Sequence

HERE = Path(__file__).parent
FILTERED = HERE / "filtered"
SYNTENY = HERE / "synteny"
CAND = HERE / "retrieval_uniprot"
OUT = HERE / "evolution_prelim"
FASTTREE = HERE / "tools" / "FastTree"
REF_LABEL = "REFERENCE"


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif name:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", default="r2pyocin")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    fam = args.family

    OUT.mkdir(parents=True, exist_ok=True)
    refs = json.loads((HERE / "reference_sequences.json").read_text())
    ref_seq = refs[fam]["seq"]
    filt = pd.read_csv(FILTERED / f"{fam}_filtered.csv")
    pool = read_fasta(CAND / f"{fam}_candidates.fasta")

    syn_path = SYNTENY / f"{fam}_synteny.csv"
    syn = pd.read_csv(syn_path) if syn_path.exists() else None
    if syn is not None:
        syn_by_acc = dict(zip(syn.accession, syn.n_categories.fillna(0)))
        sup_by_acc = dict(zip(syn.accession, syn.synteny_supported.fillna(False)))
    else:
        syn_by_acc, sup_by_acc = {}, {}

    records = [Sequence(REF_LABEL.encode(), ref_seq.encode())]
    labels = {REF_LABEL: {"uniprot": (refs[fam].get("uniprot") or [None])[0],
                          "organism": refs[fam].get("organism"), "pident": 1.0,
                          "ref_cov": 1.0, "accession": "REFERENCE",
                          "n_categories": None, "synteny_supported": None,
                          "is_reference": True}}
    missing = 0
    for _, r in filt.iterrows():
        acc = r["uniprot_accession"]
        seq = pool.get(acc)
        if not seq:
            missing += 1
            continue
        records.append(Sequence(acc.encode(), seq.encode()))
        labels[acc] = {"uniprot": acc, "organism": r["organism"],
                       "species": r["species"], "pident": float(r["pident"]),
                       "ref_cov": float(r["ref_cov"]),
                       "accession": r["accession"],
                       "accession_source": r["accession_source"],
                       "n_categories": syn_by_acc.get(r["accession"]),
                       "synteny_supported": sup_by_acc.get(r["accession"]),
                       "is_reference": False}
    if missing:
        sys.exit(f"ABORT: {missing} filtered rows had no sequence in the "
                 "candidate FASTA; the two files are out of step.")
    print(f"{fam}: aligning {len(records)} sequences "
          f"(reference + {len(records) - 1} candidates)", flush=True)

    msa = Aligner(threads=args.threads).align(records)
    afa = OUT / f"{fam}_prelim.afa"
    afa.write_text("".join(
        f">{s.id.decode()}\n{s.sequence.decode()}\n" for s in msa))
    width = len(msa[0].sequence.decode())
    print(f"  MSA {len(msa)} x {width} -> {afa.name}", flush=True)

    nwk = OUT / f"{fam}_prelim.nwk"
    with open(afa) as fi, open(nwk, "w") as fo:
        res = subprocess.run([str(FASTTREE), "-lg", "-gamma", "-quiet"],
                             stdin=fi, stdout=fo, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        sys.exit(f"ABORT: FastTree failed: {res.stderr[:400]}")
    print(f"  tree -> {nwk.name}", flush=True)

    tree = dendropy.Tree.get(path=str(nwk), schema="newick",
                             preserve_underscores=True)
    tree.reroot_at_midpoint(update_bipartitions=False)
    ref_leaf = next((l for l in tree.leaf_node_iter()
                     if l.taxon and l.taxon.label == REF_LABEL), None)
    if ref_leaf is None:
        sys.exit("ABORT: the reference leaf is not in the tree")

    # nested clades from the reference out to the root: the Stage B profile
    rows, node, depth = [], ref_leaf.parent_node, 0
    while node is not None:
        members = [l.taxon.label for l in node.leaf_iter()
                   if l.taxon and l.taxon.label != REF_LABEL]
        if members:
            sub = [labels[m] for m in members if m in labels]
            n_syn = sum(1 for s in sub if s.get("synteny_supported") is True)
            n_scored = sum(1 for s in sub if s.get("n_categories") is not None)
            idents = [s["pident"] for s in sub]
            rows.append({
                "level": depth, "clade_size": len(members),
                "stem_length": round(node.edge_length or 0.0, 5),
                "support": node.label or "",
                "pident_min": round(min(idents), 3),
                "pident_median": round(sorted(idents)[len(idents) // 2], 3),
                "pident_max": round(max(idents), 3),
                "synteny_scored": n_scored,
                "synteny_supported": n_syn,
                "synteny_supported_frac": (round(n_syn / n_scored, 3)
                                           if n_scored else None),
            })
        node, depth = node.parent_node, depth + 1
    profile = pd.DataFrame(rows)
    profile.to_csv(OUT / f"{fam}_clade_profile.csv", index=False)

    leaf_tbl = pd.DataFrame([
        {"label": k, **v} for k, v in labels.items() if k != REF_LABEL])
    # four-way classification at every nested level, all provisional
    clade_sets, node, depth = {}, ref_leaf.parent_node, 0
    while node is not None:
        clade_sets[depth] = {l.taxon.label for l in node.leaf_iter()
                             if l.taxon and l.taxon.label != REF_LABEL}
        node, depth = node.parent_node, depth + 1
    for d, members in clade_sets.items():
        leaf_tbl[f"in_clade_L{d}"] = leaf_tbl.label.isin(members)
    leaf_tbl.to_csv(OUT / f"{fam}_leaf_table.csv", index=False)

    scored = leaf_tbl[leaf_tbl.n_categories.notna()]
    fourway = {}
    for d in sorted(clade_sets):
        if not (0 < len(clade_sets[d]) < len(leaf_tbl)):
            continue
        inc = scored[f"in_clade_L{d}"]
        sup = scored.synteny_supported.fillna(False).astype(bool)
        fourway[f"L{d}"] = {
            "clade_size": int(len(clade_sets[d])),
            "in_clade_synteny_supported": int((inc & sup).sum()),
            "in_clade_not_supported": int((inc & ~sup).sum()),
            "out_clade_synteny_supported": int((~inc & sup).sum()),
            "out_clade_not_supported": int((~inc & ~sup).sum()),
        }

    summary = {
        "family": fam, "built_at": datetime.now(timezone.utc)
        .strftime("%Y-%m-%dT%H:%M:%SZ"),
        "status": "PRELIMINARY QC SCREEN -- not an orthology decision, "
                  "not a final phylogeny, no rates computed",
        "inputs": {"filtered_rows": int(len(filt)),
                   "sequences_aligned": len(records),
                   "synteny_rows": int(len(syn)) if syn is not None else 0},
        "msa": {"tool": "pyfamsa (FAMSA)", "taxa": len(msa), "columns": width},
        "tree": {"tool": "FastTree LG+Gamma", "rooting": "midpoint",
                 "note": "QC screen only; model selection remains open"},
        "reference_placement": {
            "sister_clade_size": int(profile.clade_size.iloc[0]),
            "nested_levels": int(len(profile)),
            "levels_to_half_the_pool": int(
                (profile.clade_size < len(leaf_tbl) / 2).sum()),
        },
        "four_way_by_level_PROVISIONAL": fourway,
        "synteny": {
            "candidates_scored": int(scored.shape[0]),
            "supported": int(scored.synteny_supported.fillna(False).sum()),
            "n_categories_distribution":
                scored.n_categories.value_counts().sort_index()
                .astype(int).to_dict() if len(scored) else {},
        },
    }
    (OUT / f"{fam}_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\n  reference sister clade: {profile.clade_size.iloc[0]} members, "
          f"stem {profile.stem_length.iloc[0]}")
    print(profile.head(12).to_string(index=False))
    print(f"\n  wrote {OUT}/{fam}_" + "{prelim.afa,prelim.nwk,clade_profile.csv,"
          "leaf_table.csv,summary.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
