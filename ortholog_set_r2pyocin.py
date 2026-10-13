#!/usr/bin/env python3
"""R2 pyocin: paralog test on the candidate ortholog clade, then feasibility.

Two things happen here, in order, and the second is only meaningful if the
first comes out clean.

STAGE 1 -- paralog test.
  The frozen gates keep one representative per species by bitscore. 221 of the
  967 species carry more than one gated sheath candidate, so for those species
  the retained sequence is a CHOICE among paralogs, and a wrong choice puts a
  paralog inside the ortholog set. This stage puts the dropped paralogs back:
  the candidate clade's members plus every same-species gated candidate that
  dereplication discarded are realigned and re-treed, and each dropped paralog
  is scored as falling INSIDE the clade (an in-paralog -- the clade is then not
  a clean ortholog set) or OUTSIDE it (dereplication picked the positional
  ortholog and the clade is clean).

STAGE 2 -- evolutionary feasibility at the structural sites.
  For the resulting set: number of orthologs, phylogenetic diversity, occupancy
  and substitution signal at the 143 reference positions that carry a frozen
  R_i value, local alignment reliability around those positions, extreme branch
  lengths, and taxonomic redundancy. The frozen 80% occupancy threshold is
  applied as-is. No rates are estimated and Rate4Site is not run.

The clade boundary is a PARAMETER (--level), not a result: the script reports
what the set looks like at the level it is given. Choosing the level is a
judgement made on the orthology evidence, documented in the report.

Usage: python3 ortholog_set_r2pyocin.py --level 2
Outputs: orthology/r2pyocin_paralog_test.csv, r2pyocin_paralog_test.nwk,
         orthology/r2pyocin_feasibility.json, _site_diagnostics.csv,
         orthology/r2pyocin_ortholog_set.{fasta,afa,nwk}
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

import dendropy
import pandas as pd
from pyfamsa import Aligner, Sequence

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.resolve()))
import filter_hits as fh  # noqa: E402

CAND = HERE / "retrieval_uniprot"
FILTERED = HERE / "filtered"
PRELIM = HERE / "evolution_prelim"
QC = HERE / "sheath_pilot" / "qc_output"
OUT = HERE / "orthology"
FASTTREE = HERE / "tools" / "FastTree"
FAM = "r2pyocin"
REF_LABEL = "REFERENCE"
OCC_PRIMARY = 0.80            # frozen
OCC_SENS = (0.70, 0.90)       # frozen sensitivity
WINDOW = 5                    # +/- columns for local alignment reliability


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def align(records, threads=4):
    return Aligner(threads=threads).align(records)


def write_afa(msa, path: Path):
    path.write_text("".join(f">{s.id.decode()}\n{s.sequence.decode()}\n"
                            for s in msa))


def fasttree(afa: Path, nwk: Path):
    with open(afa) as fi, open(nwk, "w") as fo:
        r = subprocess.run([str(FASTTREE), "-lg", "-gamma", "-quiet"],
                           stdin=fi, stdout=fo, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.exit(f"ABORT: FastTree failed: {r.stderr[:300]}")


def clade_members(tree, label=REF_LABEL, exclude=frozenset()):
    """Smallest clade containing the reference and all non-excluded leaves."""
    ref = next(l for l in tree.leaf_node_iter() if l.taxon.label == label)
    return ref


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, default=2,
                    help="nested-ancestor level defining the candidate clade")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    refs = json.loads((HERE / "reference_sequences.json").read_text())[FAM]
    ref_seq = refs["seq"]
    leaves = pd.read_csv(OUT / f"{FAM}_leaves.csv")
    pool = read_fasta(CAND / f"{FAM}_candidates.fasta")
    hits = pd.read_csv(CAND / f"{FAM}_blastp_hits.csv")

    clade = leaves[leaves.group <= args.level]
    print(f"candidate clade at level {args.level}: {len(clade)} members, "
          f"identity {clade.pident.min():.3f}-{clade.pident.max():.3f}, "
          f"{clade.species.str.split().str[0].nunique()} genera")

    # ---------------- stage 1: put the dropped paralogs back ---------------
    gated = hits[(hits.pident >= fh.ID_LO) & (hits.pident <= fh.ID_HI)
                 & (hits.ref_cov >= fh.MIN_REF_COV)
                 & (hits.len_ratio >= fh.LEN_LO)
                 & (hits.len_ratio <= fh.LEN_HI)].copy()
    gated["species"] = gated.organism.fillna("?").map(
        lambda o: " ".join(o.replace("Candidatus ", "").split()[:2]))
    kept = set(clade.uniprot_accession)
    dropped = gated[gated.species.isin(set(clade.species))
                    & ~gated.uniprot_accession.isin(kept)]
    print(f"same-species paralogs discarded by dereplication: {len(dropped)} "
          f"across {dropped.species.nunique()} species")

    recs = [Sequence(REF_LABEL.encode(), ref_seq.encode())]
    recs += [Sequence(a.encode(), pool[a].encode()) for a in clade.uniprot_accession]
    recs += [Sequence(f"PARALOG_{a}".encode(), pool[a].encode())
             for a in dropped.uniprot_accession]
    msa = align(recs, args.threads)
    afa_p = OUT / f"{FAM}_L{args.level}_paralog_test.afa"
    nwk_p = OUT / f"{FAM}_L{args.level}_paralog_test.nwk"
    write_afa(msa, afa_p)
    fasttree(afa_p, nwk_p)

    t = dendropy.Tree.get(path=str(nwk_p), schema="newick",
                          preserve_underscores=True)
    t.reroot_at_midpoint(update_bipartitions=False)
    ref_leaf = next(l for l in t.leaf_node_iter() if l.taxon.label == REF_LABEL)
    # the clade = smallest ancestor containing the reference and >=90% of the
    # retained members; paralogs falling inside it are in-paralogs
    retained = set(clade.uniprot_accession)
    node, chosen = ref_leaf.parent_node, None
    while node is not None:
        labs = {l.taxon.label for l in node.leaf_iter()}
        if len(labs & retained) >= 0.9 * len(retained):
            chosen = node
            break
        node = node.parent_node
    if chosen is None:
        sys.exit("ABORT: no ancestor recovers 90% of the clade in the paralog test")
    inside = {l.taxon.label for l in chosen.leaf_iter()}
    inpara = sorted(x for x in inside if x.startswith("PARALOG_"))
    outpara = sorted(f"PARALOG_{a}" for a in dropped.uniprot_accession
                     if f"PARALOG_{a}" not in inside)
    # is the duplication recent (paralog is the kept member's sister) or
    # ancient (the in-paralogs form their own subclade)?
    pdm_p = dendropy.PhylogeneticDistanceMatrix.from_tree(t)
    tax = {l.taxon.label: l for l in t.leaf_node_iter()}
    inpara_mono = None
    if inpara:
        mrca = t.mrca(taxon_labels=inpara)
        leaves_under = {l.taxon.label for l in mrca.leaf_iter()}
        inpara_mono = {"mrca_leaves": len(leaves_under),
                       "in_paralogs": len(inpara),
                       "monophyletic": leaves_under == set(inpara)}
    sister_rows = []
    for a in dropped.uniprot_accession:
        lab = f"PARALOG_{a}"
        if lab not in inside or lab not in tax:
            continue
        sp = dropped[dropped.uniprot_accession == a].iloc[0].species
        mates = clade[clade.species == sp].uniprot_accession.tolist()
        if not mates:
            continue
        mate = mates[0]
        d = pdm_p.distance(tax[lab].taxon, tax[mate].taxon)
        sis = tax[lab].parent_node is tax[mate].parent_node
        sister_rows.append({"paralog": a, "species": sp, "clade_member": mate,
                            "patristic_to_member": round(float(d), 4),
                            "is_sister_leaf": bool(sis)})
    sis_tab = pd.DataFrame(sister_rows)
    if len(sis_tab):
        sis_tab.to_csv(OUT / f"{FAM}_L{args.level}_inparalog_pairs.csv", index=False)
        print(f"  in-paralog pairs: {int(sis_tab.is_sister_leaf.sum())}/"
              f"{len(sis_tab)} are sister leaves (recent duplication); "
              f"median patristic distance to their clade member "
              f"{sis_tab.patristic_to_member.median():.3f}")
    if inpara_mono:
        print(f"  in-paralogs monophyletic: {inpara_mono['monophyletic']} "
              f"(their MRCA spans {inpara_mono['mrca_leaves']} leaves)")

    rows = []
    for a in dropped.uniprot_accession:
        g = dropped[dropped.uniprot_accession == a].iloc[0]
        rows.append({"paralog": a, "species": g.species, "pident": round(g.pident, 4),
                     "bitscore": g.bitscore,
                     "falls_inside_clade": f"PARALOG_{a}" in inside})
    ptab = pd.DataFrame(rows)
    ptab.to_csv(OUT / f"{FAM}_L{args.level}_paralog_test.csv", index=False)
    n_in, n_out = int(ptab.falls_inside_clade.sum()), int((~ptab.falls_inside_clade).sum())
    print(f"  clade recovered with {len(inside & retained)}/{len(retained)} members")
    print(f"  dropped paralogs INSIDE the clade : {n_in}  <- in-paralogs")
    print(f"  dropped paralogs OUTSIDE          : {n_out}")

    # ---------------- stage 2: feasibility of the ortholog set -------------
    recs2 = [Sequence(REF_LABEL.encode(), ref_seq.encode())]
    recs2 += [Sequence(a.encode(), pool[a].encode()) for a in clade.uniprot_accession]
    msa2 = align(recs2, args.threads)
    afa = OUT / f"{FAM}_L{args.level}_ortholog_set.afa"
    nwk = OUT / f"{FAM}_L{args.level}_ortholog_set.nwk"
    write_afa(msa2, afa)
    (OUT / f"{FAM}_L{args.level}_ortholog_set.fasta").write_text(
        f">{REF_LABEL}|{refs['uniprot'][0]}|{refs['pdb']}\n{ref_seq}\n"
        + "".join(f">{r.uniprot_accession}|{r.accession}|{r.organism}\n"
                  f"{pool[r.uniprot_accession]}\n" for _, r in clade.iterrows()))
    fasttree(afa, nwk)

    seqs = read_fasta(afa)
    width = len(seqs[REF_LABEL])
    n_seq = len(seqs)
    ref_aln = seqs[REF_LABEL]
    # reference position (1-based) -> alignment column
    pos2col, p = {}, 0
    for i, ch in enumerate(ref_aln):
        if ch != "-":
            p += 1
            pos2col[p] = i

    R = pd.read_csv(QC / f"R_{FAM}.csv")
    sites = R[R.in_domain].copy()
    sites["seqres_pos"] = sites.seqres_pos.astype(int)
    srows = []
    for _, s in sites.iterrows():
        col = pos2col.get(s.seqres_pos)
        if col is None:
            srows.append({"seqres_pos": s.seqres_pos, "R": s.R, "column": None,
                          "occupancy": 0.0, "n_aa": 0, "entropy": 0.0,
                          "window_occupancy": 0.0})
            continue
        column = [v[col] for v in seqs.values()]
        aa = [c for c in column if c != "-"]
        occ = len(aa) / n_seq
        cnt = Counter(aa)
        ent = -sum((c / len(aa)) * math.log2(c / len(aa)) for c in cnt.values()) \
            if aa else 0.0
        lo, hi = max(0, col - WINDOW), min(width, col + WINDOW + 1)
        wocc = sum(1 for j in range(lo, hi)
                   for v in seqs.values() if v[j] != "-") / (n_seq * (hi - lo))
        srows.append({"seqres_pos": s.seqres_pos, "R": round(float(s.R), 4),
                      "column": col, "occupancy": round(occ, 4),
                      "n_aa": len(cnt), "entropy": round(ent, 4),
                      "major_aa_frac": round(max(cnt.values()) / len(aa), 4) if aa else 0,
                      "window_occupancy": round(wocc, 4)})
    sdiag = pd.DataFrame(srows)
    sdiag.to_csv(OUT / f"{FAM}_L{args.level}_site_diagnostics.csv", index=False)

    t2 = dendropy.Tree.get(path=str(nwk), schema="newick", preserve_underscores=True)
    t2.reroot_at_midpoint(update_bipartitions=False)
    tips = [l.edge.length or 0.0 for l in t2.leaf_node_iter()]
    tree_len = sum(n.edge_length or 0.0 for n in t2.preorder_node_iter())
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(t2)
    dists = [pdm.distance(a.taxon, b.taxon)
             for i, a in enumerate(t2.leaf_node_iter())
             for j, b in enumerate(t2.leaf_node_iter()) if j > i]
    q3, q1 = pd.Series(tips).quantile(0.75), pd.Series(tips).quantile(0.25)
    iqr = q3 - q1
    outliers = [l.taxon.label for l in t2.leaf_node_iter()
                if (l.edge.length or 0) > q3 + 3 * iqr]

    genera = Counter(clade.species.str.split().str[0])
    tot = sum(genera.values())
    shannon = -sum((c / tot) * math.log(c / tot) for c in genera.values())

    feas = {
        "family": FAM, "clade_level": args.level,
        "status": "feasibility diagnostics only -- no rates estimated",
        "paralog_test": {
            "dropped_paralogs_tested": int(len(ptab)),
            "inside_clade_in_paralogs": n_in,
            "outside_clade": n_out,
            "clade_members_recovered": f"{len(inside & retained)}/{len(retained)}",
            "in_paralog_topology": inpara_mono,
            "in_paralog_pairs_sister_leaves": (int(sis_tab.is_sister_leaf.sum())
                                               if len(sis_tab) else 0),
            "in_paralog_pairs_tested": int(len(sis_tab)),
            "in_paralog_median_patristic": (round(float(sis_tab.patristic_to_member.median()), 4)
                                            if len(sis_tab) else None),
            "verdict": ("clean -- dereplication kept the clade member in every case"
                        if n_in == 0 else
                        f"{n_in} in-paralogs: the clade is not a clean ortholog set"),
        },
        "n_orthologs": int(len(clade)),
        "phylogenetic_diversity": {
            "tree_length": round(tree_len, 3),
            "mean_patristic_distance": round(float(pd.Series(dists).mean()), 4),
            "median_patristic_distance": round(float(pd.Series(dists).median()), 4),
            "max_patristic_distance": round(float(max(dists)), 4),
            "identity_min": round(float(clade.pident.min()), 4),
            "identity_median": round(float(clade.pident.median()), 4),
        },
        "alignment": {"columns": width, "sequences": n_seq,
                      "columns_occupancy_ge_50pct":
                          int(sum(1 for i in range(width)
                                  if sum(1 for v in seqs.values() if v[i] != "-")
                                  / n_seq >= 0.5)),
                      "reference_columns": len(pos2col)},
        "structural_sites": {
            "n_sites": int(len(sdiag)),
            "mapped_to_a_column": int(sdiag.column.notna().sum()),
            f"occupancy_ge_{int(OCC_PRIMARY*100)}pct":
                int((sdiag.occupancy >= OCC_PRIMARY).sum()),
            f"occupancy_ge_{int(OCC_SENS[0]*100)}pct":
                int((sdiag.occupancy >= OCC_SENS[0]).sum()),
            f"occupancy_ge_{int(OCC_SENS[1]*100)}pct":
                int((sdiag.occupancy >= OCC_SENS[1]).sum()),
            "median_occupancy": round(float(sdiag.occupancy.median()), 4),
            "median_window_occupancy": round(float(sdiag.window_occupancy.median()), 4),
            "sites_window_occupancy_lt_80pct":
                int((sdiag.window_occupancy < 0.80).sum()),
            "invariant_sites": int((sdiag.n_aa <= 1).sum()),
            "sites_major_aa_ge_95pct": int((sdiag.major_aa_frac >= 0.95).sum()),
            "median_distinct_aa": float(sdiag.n_aa.median()),
            "median_entropy": round(float(sdiag.entropy.median()), 4),
            "entropy_quartiles": [round(float(x), 4) for x in
                                  sdiag.entropy.quantile([.25, .5, .75])],
        },
        "branch_lengths": {
            "median_tip": round(float(pd.Series(tips).median()), 5),
            "max_tip": round(float(max(tips)), 5),
            "outlier_tips_gt_q3_plus_3iqr": len(outliers),
            "outlier_labels": outliers[:20],
        },
        "taxonomic_redundancy": {
            "distinct_species": int(clade.species.nunique()),
            "distinct_genera": len(genera),
            "largest_genus": genera.most_common(1)[0][0],
            "largest_genus_frac": round(genera.most_common(1)[0][1] / tot, 3),
            "genus_shannon": round(shannon, 3),
            "effective_genera": round(math.exp(shannon), 2),
        },
    }
    (OUT / f"{FAM}_L{args.level}_feasibility.json").write_text(json.dumps(feas, indent=2))

    print(f"\northolog set: {len(clade)} sequences, {width} columns")
    print(f"  structural sites mapped {feas['structural_sites']['mapped_to_a_column']}"
          f"/{len(sdiag)}; occupancy >=80% in "
          f"{feas['structural_sites']['occupancy_ge_80pct']}")
    print(f"  invariant {feas['structural_sites']['invariant_sites']}, "
          f"major-aa>=95% {feas['structural_sites']['sites_major_aa_ge_95pct']}, "
          f"median distinct aa {feas['structural_sites']['median_distinct_aa']}, "
          f"median entropy {feas['structural_sites']['median_entropy']}")
    print(f"  tree length {feas['phylogenetic_diversity']['tree_length']}, "
          f"mean patristic {feas['phylogenetic_diversity']['mean_patristic_distance']}")
    print(f"  genera {len(genera)}, largest "
          f"{feas['taxonomic_redundancy']['largest_genus']} "
          f"{feas['taxonomic_redundancy']['largest_genus_frac']:.0%}, "
          f"effective {feas['taxonomic_redundancy']['effective_genera']}")
    print(f"  outlier tips {len(outliers)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
