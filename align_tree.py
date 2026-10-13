#!/usr/bin/env python3
"""Steps 4-10: orthology by tree placement + locus context, then align and STOP.

Orthology is NOT decided by a fixed annotation count. A genuine ortholog can
score zero synteny markers because its neighbours are annotated "hypothetical
protein", and the four systems differ in locus architecture and annotation
quality. So:

  Stage A  preliminary MSA + tree over ALL filtered candidates
  Stage B  locate the reference's clade; report the stem-length profile from the
           reference to the root so deep subfamily splits are visible
  Stage C  classify each candidate by clade membership x synteny support
  Stage D  subsample within the orthologous clade
  Stage E  final MSA + tree, structural-site occupancy

The trees here are QC SCREENS (FastTree, LG+gamma) for judging whether the
sequences form one coherent orthologous clade. They are not the final
phylogenies; model selection and rigorous ML inference are a later decision.

NO EVOLUTIONARY RATES ARE COMPUTED.
"""

from __future__ import annotations

import json
import subprocess
import time
import urllib.parse
import urllib.request
from pathlib import Path

import dendropy
import numpy as np
import pandas as pd
from pyfamsa import Aligner, Sequence
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

HERE = Path(__file__).parent
SYN = HERE / "synteny"
OUT = HERE / "evolution"
QC = HERE / "sheath_pilot" / "qc_output"
FASTTREE = HERE / "tools" / "FastTree"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

TARGET_N, MIN_N = 50, 30
OCC_PRIMARY, OCC_SENS = 0.80, [0.70, 0.90]
PAUSE = 0.4
REF = "REFERENCE"

R_FILE = {"r2pyocin": ("R_r2pyocin.csv", "F1"), "salcis": ("R_salcis.csv", "F1"),
          "diffocin": ("R_diffocin.csv", "F1"), "afp_Afp2": ("R_afp.csv", "F1"),
          "afp_Afp3": ("R_afp.csv", "F2")}


def fetch_seqs(accs):
    seqs = {}
    for i in range(0, len(accs), 100):
        q = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(
            {"db": "protein", "id": ",".join(accs[i:i + 100]),
             "rettype": "fasta", "retmode": "text"})
        try:
            with urllib.request.urlopen(q, timeout=180) as fh:
                txt = fh.read().decode("utf-8", "replace")
        except Exception as exc:
            print(f"    fetch chunk {i} failed: {exc}", flush=True)
            continue
        acc = None
        for line in txt.splitlines():
            if line.startswith(">"):
                acc = line[1:].split()[0]
                seqs[acc] = []
            elif acc:
                seqs[acc].append(line.strip())
        time.sleep(PAUSE)
    return {k: "".join(v) for k, v in seqs.items() if v}


def align(records):
    msa = Aligner(guide_tree="upgma", threads=0).align(
        [Sequence(k.encode(), v.encode()) for k, v in records])
    return [(s.id.decode(), s.sequence.decode()) for s in msa]


def write_fasta(path, records):
    with open(path, "w") as fh:
        for n, s in records:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def fasttree(fin, fout):
    with open(fin) as fi, open(fout, "w") as fo:
        r = subprocess.run([str(FASTTREE), "-lg", "-gamma", "-quiet"],
                           stdin=fi, stdout=fo, stderr=subprocess.PIPE, text=True)
    return r.returncode, (r.stderr or "")[-300:]


def identity_matrix(aln):
    arr = np.frombuffer("".join(s for _, s in aln).encode(), dtype="S1").reshape(
        len(aln), -1)
    gap = arr == b"-"
    n = len(aln)
    M = np.ones((n, n))
    for i in range(n):
        both = ~gap[i] & ~gap
        same = (arr[i] == arr) & both
        cnt = both.sum(axis=1)
        M[i] = np.where(cnt > 0, same.sum(axis=1) / np.maximum(cnt, 1), 0.0)
    M = (M + M.T) / 2
    np.fill_diagonal(M, 1.0)
    return M


def clade_profile(tree_path, labels):
    """Walk from the reference leaf to the root, recording clade size and the
    stem branch subtending each clade. A deep subfamily split shows up as an
    unusually long stem: that is where a set of 'homologs' stops being one
    orthologous group."""
    t = dendropy.Tree.get(path=str(tree_path), schema="newick",
                          preserve_underscores=True)
    try:
        t.reroot_at_midpoint(update_bipartitions=False)
    except Exception:
        pass
    leaf = next((l for l in t.leaf_node_iter()
                 if l.taxon and l.taxon.label.strip("'") == REF), None)
    if leaf is None:
        return [], {}
    prof, node = [], leaf.parent_node
    while node is not None:
        members = {x.taxon.label.strip("'") for x in node.leaf_iter() if x.taxon}
        prof.append({"n": len(members), "stem": node.edge.length or 0.0,
                     "members": members})
        node = node.parent_node
    # terminal branch lengths, for outlier detection
    tips = {l.taxon.label.strip("'"): (l.edge.length or 0.0)
            for l in t.leaf_node_iter() if l.taxon}
    return prof, tips


def choose_clade(prof, syn_map):
    """Largest clade below the first large jump in stem length."""
    if not prof:
        return set()
    stems = np.array([p["stem"] for p in prof])
    best, chosen = None, prof[-1]["members"]
    med = np.median(stems[stems > 0]) if (stems > 0).any() else 0.0
    for i, p in enumerate(prof):
        if p["n"] < 4:
            continue
        if med > 0 and p["stem"] > 4 * med:
            chosen = prof[i]["members"]
            best = i
            break
    if best is None:
        # no clear split: take the largest clade that is still majority-syntenic
        for p in reversed(prof):
            mem = [m for m in p["members"] if m != REF]
            if mem and np.mean([syn_map.get(m, False) for m in mem]) >= 0.5:
                chosen = p["members"]
                break
    return chosen


def subsample(aln, target, prefer):
    if len(aln) <= target:
        return list(range(len(aln)))
    M = identity_matrix(aln)
    D = np.clip(1.0 - M, 0, None)
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    lab = fcluster(Z, t=target, criterion="maxclust")
    chosen = []
    for c in sorted(set(lab)):
        idx = np.where(lab == c)[0]
        pref = [i for i in idx if prefer.get(aln[i][0], False)]
        pool = pref if pref else list(idx)
        sub = D[np.ix_(pool, pool)]
        chosen.append(int(pool[int(np.argmin(sub.mean(axis=1)))]))
    ref = next((i for i, (n, _) in enumerate(aln) if n == REF), None)
    if ref is not None and ref not in chosen:
        chosen.append(ref)
    return sorted(set(chosen))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    refs = json.load(open(HERE / "reference_sequences.json"))
    tables = {f.stem.replace("_synteny", ""): pd.read_csv(f)
              for f in SYN.glob("*_synteny.csv")}
    if not tables:
        print("no synteny output yet")
        return

    # AFP paired-locus statistic -- reported, and used only as preference
    afp_stat = {}
    if "afp_Afp2" in tables and "afp_Afp3" in tables:
        n2 = set(tables["afp_Afp2"].nuc.dropna())
        n3 = set(tables["afp_Afp3"].nuc.dropna())
        shared = n2 & n3
        afp_stat = {"loci_Afp2": len(n2), "loci_Afp3": len(n3),
                    "shared_loci": len(shared),
                    "frac_Afp2_paired": len(shared) / max(len(n2), 1),
                    "frac_Afp3_paired": len(shared) / max(len(n3), 1)}
        for k in ("afp_Afp2", "afp_Afp3"):
            tables[k]["paired_locus"] = tables[k].nuc.isin(shared)
        print(f"\nAFP paired-locus recovery: {len(shared)} loci carry both "
              f"(Afp2 {afp_stat['frac_Afp2_paired']:.1%}, "
              f"Afp3 {afp_stat['frac_Afp3_paired']:.1%})", flush=True)

    rows = []
    for fam, t in sorted(tables.items()):
        ref = refs[fam]
        print(f"\n=== {fam} ({len(t)} candidates)", flush=True)
        seqs = fetch_seqs(t.accession.tolist())
        ok = [(a, s) for a, s in seqs.items()
              if 0.75 <= len(s) / len(ref["seq"]) <= 1.30]
        records = [(REF, ref["seq"])] + ok
        print(f"  sequences usable: {len(ok)}", flush=True)
        if len(ok) < MIN_N:
            rows.append({"family": fam, "status": "INSUFFICIENT",
                         "candidates": len(t), "usable": len(ok)})
            continue

        # Stage A: preliminary alignment + tree
        aln0 = align(records)
        write_fasta(OUT / f"{fam}_prelim_msa.fasta", aln0)
        code0, _ = fasttree(OUT / f"{fam}_prelim_msa.fasta",
                            OUT / f"{fam}_prelim.tree")

        syn_map = dict(zip(t.accession, t.synteny_supported))
        prof, tips = clade_profile(OUT / f"{fam}_prelim.tree",
                                   [n for n, _ in aln0]) if code0 == 0 else ([], {})
        clade = choose_clade(prof, syn_map)
        in_clade = {n for n, _ in aln0 if n in clade or n == REF}

        # Stage C: classify
        cls = []
        for n, _ in aln0:
            if n == REF:
                continue
            cls.append({"accession": n, "in_clade": n in in_clade,
                        "synteny": bool(syn_map.get(n, False))})
        cdf = pd.DataFrame(cls)
        n_orth = int((cdf.in_clade & cdf.synteny).sum())
        n_clade_only = int((cdf.in_clade & ~cdf.synteny).sum())
        n_syn_only = int((~cdf.in_clade & cdf.synteny).sum())
        n_neither = int((~cdf.in_clade & ~cdf.synteny).sum())
        print(f"  clade={len(in_clade)}  orth(clade+syn)={n_orth} "
              f"clade-only={n_clade_only} syn-only(ambiguous paralog)={n_syn_only} "
              f"neither={n_neither}", flush=True)

        # Stage D/E: subsample within the clade, realign, final tree
        keep = [(n, s) for n, s in records if n in in_clade]
        if len(keep) < MIN_N:
            keep = records
            print(f"  clade too small ({len(in_clade)}); keeping all candidates "
                  f"and flagging for manual inspection", flush=True)
        a1 = align(keep)
        idx = subsample(a1, TARGET_N, {a: bool(syn_map.get(a, False)) for a, _ in keep})
        sel = [keep[i] for i in idx]
        aln = align(sel)
        write_fasta(OUT / f"{fam}_msa.fasta", aln)
        code, err = fasttree(OUT / f"{fam}_msa.fasta", OUT / f"{fam}.tree")
        _, tips_f = clade_profile(OUT / f"{fam}.tree", [n for n, _ in aln])

        M = identity_matrix(aln)
        iu = np.triu_indices(len(M), 1)
        med_id = float(np.median(M[iu]))

        arr = np.array([list(s) for _, s in aln])
        refrow = dict(aln)[REF]
        col_of_pos, p = {}, 0
        for c, ch in enumerate(refrow):
            if ch != "-":
                p += 1
                col_of_pos[p] = c
        ref_cov = len(col_of_pos) / len(ref["seq"])
        gapfrac = (arr == "-").mean(axis=0)

        occ = []
        fn, fkey = R_FILE.get(fam, (None, None))
        if fn and (QC / fn).exists():
            rdf = pd.read_csv(QC / fn)
            rdf = rdf[(rdf.family == fkey) & rdf.in_domain & rdf.seqres_pos.notna()]
            for _, r in rdf.iterrows():
                c = col_of_pos.get(int(r.seqres_pos))
                occ.append({"canon": r.canon, "seqres_pos": int(r.seqres_pos),
                            "R": r.R,
                            "occupancy": np.nan if c is None
                            else float((arr[:, c] != "-").mean())})
        odf = pd.DataFrame(occ)
        if len(odf):
            odf.to_csv(OUT / f"{fam}_site_occupancy.csv", index=False)
        o = odf.occupancy.dropna() if len(odf) else pd.Series(dtype=float)

        tv = np.array(list(tips_f.values())) if tips_f else np.array([0.0])
        mad = float(np.median(np.abs(tv - np.median(tv)))) or 1e-9
        outliers = int((tv > np.median(tv) + 5 * mad).sum())

        rows.append({
            "family": fam, "status": "OK", "candidates": len(t),
            "usable": len(ok), "clade": len(in_clade),
            "orth_clade_and_syn": n_orth, "clade_only": n_clade_only,
            "ambiguous_paralog": n_syn_only, "neither": n_neither,
            "final_N": len(aln), "median_identity": round(med_id, 3),
            "ref_coverage": round(ref_cov, 3),
            "sites_in_domain": len(odf),
            "sites_occ_ge_70": int((o >= 0.70).sum()),
            "sites_occ_ge_80": int((o >= OCC_PRIMARY).sum()),
            "sites_occ_ge_90": int((o >= 0.90).sum()),
            "sites_unmapped": int(odf.occupancy.isna().sum()) if len(odf) else 0,
            "aln_len": arr.shape[1],
            "cols_gap_gt50pct": int((gapfrac > 0.5).sum()),
            "longest_gappy_run": int(max((len(list(g)) for k, g in
                                          __import__("itertools").groupby(gapfrac > 0.5)
                                          if k), default=0)),
            "max_tip_branch": round(float(tv.max()), 3),
            "long_branch_outliers": outliers,
            "tree": "ok" if code == 0 else "FAILED",
        })
        print(f"  final N={len(aln)} medianID={med_id:.3f} refcov={ref_cov:.3f} "
              f"sites>=80%occ={int((o >= OCC_PRIMARY).sum())}/{len(odf)}", flush=True)

    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "checkpoint_qc.csv", index=False)
    if afp_stat:
        json.dump(afp_stat, open(OUT / "afp_paired_locus.json", "w"), indent=1)
    print("\n" + rep.to_string(index=False))
    print("\nSTOP: alignments and QC trees built. Rate4Site NOT run.")


if __name__ == "__main__":
    main()
