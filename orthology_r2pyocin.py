#!/usr/bin/env python3
"""R2 pyocin orthology validation: clade structure x architecture x locus.

Retrieving homologs is not establishing orthology. This script assembles the
evidence needed to decide where the orthologous group ends, and deliberately
does NOT decide by identity:

  1. the reference-containing clade and each successive SISTER group outwards
     (sister[k] = members gained at ancestor k, i.e. anc[k] - anc[k-1]);
  2. per group: Pfam/InterPro architecture, protein-name vocabulary, length,
     taxonomy, and genomic-neighbourhood marker categories;
  3. duplicated sheath genes -- organisms contributing more than one candidate,
     and candidates sharing a genomic contig;
  4. alignment-width diagnostics, to test whether the 1310-column MSA over a
     386-residue reference reflects family mixing or a few gappy insertions.

`--extend-synteny N` resolves up to N additional loci per group through the
unchanged `synteny.py` helpers, so architecture and locus evidence are available
for the groups the boundary decision actually turns on. That is a diagnostic
extension of coverage, not a change to the frozen scoring rule: WINDOW,
MARKERS and MIN_CATEGORIES are imported, never redefined.

Outputs: orthology/r2pyocin_groups.csv, _leaves.csv, _duplicates.csv,
         _occupancy.csv
Usage:   python3 orthology_r2pyocin.py [--extend-synteny 40]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import dendropy
import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.resolve()))
import synteny as syn_mod  # noqa: E402 -- frozen WINDOW/MARKERS/MIN_CATEGORIES
import filter_hits as fh  # noqa: E402 -- the single definition of the frozen gates

fh_gates = {"ID_LO": fh.ID_LO, "ID_HI": fh.ID_HI, "MIN_REF_COV": fh.MIN_REF_COV,
            "LEN_LO": fh.LEN_LO, "LEN_HI": fh.LEN_HI}

FILTERED = HERE / "filtered"
SYNTENY = HERE / "synteny"
PRELIM = HERE / "evolution_prelim"
OUT = HERE / "orthology"
REF_LABEL = "REFERENCE"
FAM = "r2pyocin"

# vocabulary for reading protein names -- descriptive only, scores nothing
NAME_CLASSES = {
    "sheath_generic": r"tail sheath|sheath protein|phage tail sheath",
    "t6ss": r"type vi|t6ss|vipa|vipb|tssb|tssc",
    "afp": r"antifeeding|\bafp\b",
    "pvc": r"\bpvc\b|photorhabdus virulence",
    "pyocin": r"pyocin",
    "baseplate_or_other": r"baseplate|tube|spike|fiber|fibre",
    "hypothetical": r"hypothetical|uncharacteri",
}


def load_tree():
    t = dendropy.Tree.get(path=str(PRELIM / f"{FAM}_prelim.nwk"), schema="newick",
                          preserve_underscores=True)
    t.reroot_at_midpoint(update_bipartitions=False)
    ref = next((l for l in t.leaf_node_iter()
                if l.taxon and l.taxon.label == REF_LABEL), None)
    if ref is None:
        sys.exit("ABORT: reference leaf absent from the tree")
    return t, ref


def nested_groups(ref_leaf):
    """anc[k] and sister[k] = anc[k] - anc[k-1], outwards from the reference."""
    anc, node = [], ref_leaf.parent_node
    stems, supports = [], []
    while node is not None:
        anc.append({l.taxon.label for l in node.leaf_iter()
                    if l.taxon and l.taxon.label != REF_LABEL})
        stems.append(node.edge_length or 0.0)
        supports.append(node.label or "")
        node = node.parent_node
    sisters, prev = [], set()
    for s in anc:
        sisters.append(s - prev)
        prev = s
    return anc, sisters, stems, supports


def architecture_profile(sub: pd.DataFrame, ref_arch: frozenset) -> dict:
    archs = sub.pfam.fillna("").map(
        lambda s: frozenset(x for x in s.split(";") if x))
    exact = int((archs == ref_arch).sum())
    missing_core = int(archs.map(lambda a: not (a & ref_arch)).sum())
    top = Counter(";".join(sorted(a)) or "none" for a in archs).most_common(2)
    return {"arch_identical_to_reference": exact,
            "arch_identical_frac": round(exact / max(len(sub), 1), 3),
            "arch_no_shared_pfam": missing_core,
            "arch_top": " | ".join(f"{k} ({v})" for k, v in top)}


def name_profile(sub: pd.DataFrame) -> dict:
    t = sub.title.fillna("").str.lower()
    out = {}
    for cls, rx in NAME_CLASSES.items():
        out[f"name_{cls}"] = int(t.str.contains(rx, regex=True).sum())
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extend-synteny", type=int, default=0,
                    help="resolve up to N additional loci per group")
    ap.add_argument("--max-groups-extended", type=int, default=6)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    filt = pd.read_csv(FILTERED / f"{FAM}_filtered.csv")
    syn = pd.read_csv(SYNTENY / f"{FAM}_synteny.csv")
    extra_path = SYNTENY / f"{FAM}_synteny_extra.csv"
    if extra_path.exists():
        syn = pd.concat([syn, pd.read_csv(extra_path)], ignore_index=True)
        syn = syn.drop_duplicates("accession", keep="first")

    refs = json.loads((HERE / "reference_sequences.json").read_text())[FAM]
    ref_arch = frozenset(["PF22671", "PF04984", "PF17482"])

    tree, ref_leaf = load_tree()
    anc, sisters, stems, supports = nested_groups(ref_leaf)

    filt = filt.set_index("uniprot_accession", drop=False)
    group_of = {}
    for k, members in enumerate(sisters):
        for m in members:
            group_of[m] = k
    filt["group"] = filt.uniprot_accession.map(group_of)
    unassigned = int(filt.group.isna().sum())
    if unassigned:
        sys.exit(f"ABORT: {unassigned} filtered rows are not in the tree")
    filt["group"] = filt.group.astype(int)

    # ---- optional targeted synteny extension ------------------------------
    if args.extend_synteny:
        have = set(syn.accession)
        want = []
        for k in sorted(filt.group.unique())[:args.max_groups_extended]:
            sub = filt[(filt.group == k) & (filt.accession_source == "refseq")]
            need = [a for a in sub.accession if a not in have]
            # spread across genera rather than taking the top of the table
            sub2 = sub[sub.accession.isin(need)].copy()
            sub2["genus"] = sub2.species.str.split().str[0]
            pick = (sub2.sample(frac=1.0, random_state=0)
                        .drop_duplicates("genus").head(args.extend_synteny))
            if len(pick) < args.extend_synteny:
                rest = sub2[~sub2.accession.isin(pick.accession)]
                pick = pd.concat([pick, rest.head(args.extend_synteny - len(pick))])
            want += pick.accession.tolist()
        want = list(dict.fromkeys(want))
        print(f"extending synteny for {len(want)} accessions", flush=True)
        if want:
            locs = syn_mod.ipg_locations(want)
            print(f"  locations resolved {len(locs)}/{len(want)}", flush=True)
            rows = []
            for n, acc in enumerate(want, 1):
                L = locs.get(acc)
                if not L:
                    rows.append({"accession": acc, "nuc": None, "n_categories": 0,
                                 "categories": "", "synteny_supported": False,
                                 "n_neighbour_cds": 0,
                                 "reason": "no genomic location"})
                    continue
                prods = syn_mod.neighbourhood(L["nuc"], L["start"], L["stop"])
                cats = syn_mod.score(prods) - {"sheath"}
                n_sheath = sum(1 for p in prods
                               if syn_mod.COMPILED["sheath"].search(p))
                rows.append({"accession": acc, "nuc": L["nuc"],
                             "assembly": L.get("assembly", ""),
                             "n_neighbour_cds": len(prods),
                             "n_sheath_neighbours": n_sheath,
                             "n_categories": len(cats),
                             "categories": ";".join(sorted(cats)),
                             "synteny_supported": len(cats) >= syn_mod.MIN_CATEGORIES,
                             "reason": ""})
                if n % 25 == 0:
                    print(f"    {n}/{len(want)} neighbourhoods", flush=True)
            new = pd.DataFrame(rows)
            if extra_path.exists():
                new = pd.concat([pd.read_csv(extra_path), new], ignore_index=True)
                new = new.drop_duplicates("accession", keep="last")
            new.to_csv(extra_path, index=False)
            syn = pd.concat([pd.read_csv(SYNTENY / f"{FAM}_synteny.csv"), new],
                            ignore_index=True).drop_duplicates("accession")

    syn_cat = dict(zip(syn.accession, syn.n_categories.fillna(0)))
    syn_sup = dict(zip(syn.accession, syn.synteny_supported.fillna(False)))
    syn_nuc = dict(zip(syn.accession, syn.nuc))
    syn_cats = dict(zip(syn.accession, syn.categories.fillna("")))
    filt["n_categories"] = filt.accession.map(syn_cat)
    filt["synteny_supported"] = filt.accession.map(syn_sup)
    filt["nuc"] = filt.accession.map(syn_nuc)
    filt["categories"] = filt.accession.map(syn_cats)
    if "n_sheath_neighbours" in syn.columns:
        filt["n_sheath_neighbours"] = filt.accession.map(
            dict(zip(syn.accession, syn.n_sheath_neighbours)))

    # ---- per-group comparison ---------------------------------------------
    rows = []
    for k in sorted(filt.group.unique()):
        sub = filt[filt.group == k]
        scored = sub[sub.n_categories.notna()]
        rows.append({
            "group": k, "clade_size_cumulative": len(anc[k]),
            "n": len(sub), "stem_to_parent": round(stems[k], 5),
            "support": supports[k],
            "pident_min": round(sub.pident.min(), 3),
            "pident_median": round(sub.pident.median(), 3),
            "pident_max": round(sub.pident.max(), 3),
            "len_median": int(sub.slen.median()),
            "len_min": int(sub.slen.min()), "len_max": int(sub.slen.max()),
            **architecture_profile(sub, ref_arch),
            **name_profile(sub),
            "synteny_scored": len(scored),
            "synteny_supported": int(scored.synteny_supported.fillna(False).sum()),
            "synteny_supported_frac": (
                round(float(scored.synteny_supported.fillna(False).mean()), 3)
                if len(scored) else None),
            "median_n_categories": (float(scored.n_categories.median())
                                    if len(scored) else None),
            "n_sheath_neighbours_median": (
                float(scored.n_sheath_neighbours.median())
                if "n_sheath_neighbours" in scored and
                scored.n_sheath_neighbours.notna().any() else None),
            "members_from_multiparalog_species": None,   # filled in below
            "paralog_margin_median": None,               # filled in below
            "insertion_residues_median": None,           # filled in below
            "genera": sub.species.str.split().str[0].nunique(),
            "top_genera": " | ".join(
                f"{g} ({c})" for g, c in
                Counter(sub.species.str.split().str[0]).most_common(3)),
        })
    groups = pd.DataFrame(rows)
    groups.to_csv(OUT / f"{FAM}_groups.csv", index=False)

    # ---- duplicated sheath genes ------------------------------------------
    # Paralogy has to be read BEFORE dereplication: the frozen rule keeps one
    # representative per species by bitscore, so a genome carrying several
    # prophage sheaths contributes exactly one row here and the duplication is
    # invisible in the filtered table.
    hits = pd.read_csv(HERE / "retrieval_uniprot" / f"{FAM}_blastp_hits.csv")
    gated = hits[(hits.pident >= fh_gates["ID_LO"]) & (hits.pident <= fh_gates["ID_HI"])
                 & (hits.ref_cov >= fh_gates["MIN_REF_COV"])
                 & (hits.len_ratio >= fh_gates["LEN_LO"])
                 & (hits.len_ratio <= fh_gates["LEN_HI"])].copy()
    gated["species"] = gated.organism.fillna("?").map(
        lambda o: " ".join(o.replace("Candidatus ", "").split()[:2]))
    kept = dict(zip(filt.species, filt.uniprot_accession))
    dup_rows = []
    for sp, sub in gated.groupby("species"):
        if len(sub) < 2:
            continue
        sub = sub.sort_values("bitscore", ascending=False)
        top2 = sub.bitscore.iloc[:2].tolist()
        dup_rows.append({
            "basis": "species_paralogs", "value": sp,
            "n_gated_candidates": int(len(sub)),
            "distinct_sequences": int(sub.uniprot_accession.nunique()),
            "kept_representative": kept.get(sp),
            "group_of_kept": int(filt.loc[filt.species == sp, "group"].iloc[0])
            if (filt.species == sp).any() else None,
            "bitscore_top": float(top2[0]),
            "bitscore_runner_up": float(top2[1]),
            "bitscore_margin": round(float(top2[0] - top2[1]), 1),
            "pident_range": f"{sub.pident.min():.3f}-{sub.pident.max():.3f}",
        })
    dups = pd.DataFrame(dup_rows).sort_values("n_gated_candidates", ascending=False)
    dups.to_csv(OUT / f"{FAM}_duplicates.csv", index=False)

    # contig-level: DISTINCT accessions on one contig is real tandem duplication;
    # the same accession repeated is the shared-RefSeq artifact, not a paralog
    contig = (filt[filt.nuc.notna()].groupby("nuc")
              .agg(distinct_accessions=("accession", "nunique"),
                   rows=("accession", "size")).reset_index())
    tandem = contig[contig.distinct_accessions > 1]

    # ---- alignment width diagnostics --------------------------------------
    afa = PRELIM / f"{FAM}_prelim.afa"
    seqs, name = {}, None
    for line in afa.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    width = len(next(iter(seqs.values())))
    occ = [sum(1 for s in seqs.values() if s[i] != "-") / len(seqs)
           for i in range(width)]
    ref_aln = seqs[REF_LABEL]
    ref_cols = [i for i, ch in enumerate(ref_aln) if ch != "-"]
    occ_df = pd.DataFrame({"column": range(width), "occupancy": occ,
                           "reference_has_residue": [c != "-" for c in ref_aln]})
    occ_df.to_csv(OUT / f"{FAM}_occupancy.csv", index=False)

    sparse = {i for i, o in enumerate(occ) if o < 0.05}
    ins = {k: sum(1 for i in sparse if v[i] != "-") for k, v in seqs.items()}
    filt["insertion_residues"] = filt.uniprot_accession.map(ins)
    groups["insertion_residues_median"] = groups.group.map(
        lambda k: float(filt[filt.group == k].insertion_residues.median()))
    groups.to_csv(OUT / f"{FAM}_groups.csv", index=False)
    print(f"reference insertion residues (occupancy <5%): {ins[REF_LABEL]}")

    filt.drop(columns=["pfam", "interpro"], errors="ignore").to_csv(
        OUT / f"{FAM}_leaves.csv", index=False)

    print(f"\nMSA {len(seqs)} x {width}; reference occupies {len(ref_cols)} columns")
    print(f"columns with <1% occupancy: {sum(1 for o in occ if o < 0.01)}  "
          f"<5%: {sum(1 for o in occ if o < 0.05)}  "
          f">=50%: {sum(1 for o in occ if o >= 0.5)}")
    multi_sp = set(dups.value)
    groups["members_from_multiparalog_species"] = groups.group.map(
        lambda k: int(filt[(filt.group == k)].species.isin(multi_sp).sum()))
    groups["multiparalog_frac"] = (groups.members_from_multiparalog_species
                                   / groups.n).round(3)
    margin = dict(zip(dups.value, dups.bitscore_margin))
    groups["paralog_margin_median"] = groups.group.map(
        lambda k: (round(float(pd.Series(
            [margin[s] for s in filt[filt.group == k].species if s in margin]
        ).median()), 1) if any(s in margin for s in filt[filt.group == k].species)
            else None))
    groups.to_csv(OUT / f"{FAM}_groups.csv", index=False)
    print(f"\nparalogy: {len(dups)} species carry >1 gated candidate "
          f"({int(dups.n_gated_candidates.sum() - len(dups))} extra); "
          f"{len(tandem)} contigs carry >1 DISTINCT accession")
    cols = ["group", "n", "clade_size_cumulative", "stem_to_parent", "support",
            "pident_median", "len_median", "arch_identical_frac",
            "synteny_scored", "synteny_supported_frac", "median_n_categories",
            "n_sheath_neighbours_median", "multiparalog_frac",
            "paralog_margin_median", "insertion_residues_median", "genera"]
    print("\n" + groups[cols].to_string(index=False))
    print(f"\nwrote {OUT}/{FAM}_" + "{groups,leaves,duplicates,occupancy}.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
