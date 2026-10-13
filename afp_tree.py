#!/usr/bin/env python3
"""AFP check 3: one joint phylogeny testing the positional assignments.

Builds a SINGLE tree containing reference Afp2, reference Afp3, the third
reference sheath gene, and every sheath-like copy recovered from the diagnostic
loci — assigned and unassigned alike. Separate Afp2 and Afp3 trees are not
built, because that would presuppose the distinction under test.

The positional labels from check 2 are carried onto the tips as ANNOTATION only;
the tree is built from sequence alone and knows nothing about them. Three
statistics, all defined before the tree was inspected:

  nearest-labelled-neighbour concordance
      for each labelled tip, whether its closest labelled tip by patristic
      distance carries the same label. Assumption-light: no clade definition,
      no rooting.
  label MRCA purity
      the smallest clade containing every tip of a label, and what fraction of
      its tips carry that label.
  within- vs between-label distance
      median patristic distance within each label against the median to the
      other labels.

Discordant tips are reported individually, never removed. The third reference
gene is examined for whether it forms its own lineage, sits inside Afp2 or Afp3,
or is unresolved.

Strictly diagnostic. No ortholog set is selected and no rates are estimated.

Usage: python3 afp_tree.py
Outputs: orthology/afp_joint.{afa,nwk}, orthology/afp_tree_tips.csv,
         orthology/afp_check3.json
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

import dendropy
import pandas as pd
from pyfamsa import Aligner, Sequence

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.resolve()))
import positional_orthology_r2pyocin as po  # noqa: E402

OUT = HERE / "orthology"
FASTTREE = HERE / "tools" / "FastTree"
REF = {"WP_010895804.1": "REF_Afp2", "WP_010895805.1": "REF_Afp3",
       "WP_169558217.1": "REF_third"}


def main() -> int:
    runs = pd.read_csv(OUT / "afp_sheath_runs.csv")
    ids = sorted(set(runs.protein_id.dropna()) | set(REF))
    print(f"{len(runs)} sheath copies from {runs.candidate.nunique()} loci; "
          f"{len(ids)} distinct proteins to align")
    seqs = po.fetch_proteins(ids)
    missing = [i for i in ids if i not in seqs]
    if missing:
        print(f"  no sequence for {len(missing)}: {missing[:5]}")

    label_of, org_of = {}, {}
    for _, r in runs.iterrows():
        if r.protein_id in seqs:
            label_of.setdefault(r.protein_id, r.positional_label)
            org_of.setdefault(r.protein_id, r.organism)
    for pid, nm in REF.items():
        label_of[pid] = {"REF_Afp2": "afp2_like", "REF_Afp3": "afp3_like",
                         "REF_third": "third_like"}[nm]

    recs, tip_label = [], {}
    for pid in ids:
        if pid not in seqs:
            continue
        tip = REF.get(pid, pid).replace(".", "_")
        tip_label[tip] = label_of.get(pid, "unassigned")
        recs.append(Sequence(tip.encode(), seqs[pid].encode()))
    msa = Aligner(threads=4).align(recs)
    afa = OUT / "afp_joint.afa"
    afa.write_text("".join(f">{s.id.decode()}\n{s.sequence.decode()}\n" for s in msa))
    nwk = OUT / "afp_joint.nwk"
    with open(afa) as fi, open(nwk, "w") as fo:
        r = subprocess.run([str(FASTTREE), "-lg", "-gamma", "-quiet"],
                           stdin=fi, stdout=fo, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.exit(f"ABORT: FastTree failed: {r.stderr[:300]}")
    print(f"  MSA {len(msa)} x {len(msa[0].sequence.decode())}; tree written")

    t = dendropy.Tree.get(path=str(nwk), schema="newick", preserve_underscores=True)
    t.reroot_at_midpoint(update_bipartitions=False)
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(t)
    tips = {l.taxon.label: l for l in t.leaf_node_iter()}
    labelled = [x for x in tips if tip_label.get(x) in
                ("afp2_like", "afp3_like", "third_like")]

    # ---- nearest labelled neighbour ---------------------------------------
    rows = []
    for a in labelled:
        best, bd = None, None
        for b in labelled:
            if a == b:
                continue
            d = pdm.distance(tips[a].taxon, tips[b].taxon)
            if bd is None or d < bd:
                best, bd = b, d
        rows.append({"tip": a, "label": tip_label[a], "nearest": best,
                     "nearest_label": tip_label[best], "distance": round(float(bd), 4),
                     "concordant": tip_label[a] == tip_label[best],
                     "organism": org_of.get(a.replace("_", ".", 1) if a.startswith("WP") else a, "")})
    nn = pd.DataFrame(rows)

    # ---- MRCA purity -------------------------------------------------------
    purity = {}
    for lab in ("afp2_like", "afp3_like", "third_like"):
        mem = [x for x in labelled if tip_label[x] == lab]
        if len(mem) < 2:
            purity[lab] = {"n": len(mem), "mrca_tips": None, "purity": None}
            continue
        mrca = t.mrca(taxon_labels=mem)
        under = [l.taxon.label for l in mrca.leaf_iter()]
        purity[lab] = {"n": len(mem), "mrca_tips": len(under),
                       "purity": round(sum(1 for x in under
                                           if tip_label.get(x) == lab) / len(under), 3),
                       "monophyletic": set(under) == set(mem)}

    # ---- within vs between -------------------------------------------------
    dist = {}
    for lab in ("afp2_like", "afp3_like", "third_like"):
        mem = [x for x in labelled if tip_label[x] == lab]
        oth = [x for x in labelled if tip_label[x] != lab]
        w = [pdm.distance(tips[a].taxon, tips[b].taxon)
             for i, a in enumerate(mem) for b in mem[i + 1:]]
        bt = [pdm.distance(tips[a].taxon, tips[b].taxon) for a in mem for b in oth]
        dist[lab] = {"n": len(mem),
                     "within_median": round(float(pd.Series(w).median()), 3) if w else None,
                     "between_median": round(float(pd.Series(bt).median()), 3) if bt else None}

    # ---- where does the third reference sit? ------------------------------
    third = tips.get("REF_third")
    sister = None
    if third is not None and third.parent_node is not None:
        sib = [l.taxon.label for l in third.parent_node.leaf_iter()
               if l.taxon.label != "REF_third"]
        sister = {"sister_tips": len(sib),
                  "sister_labels": dict(Counter(tip_label.get(x, "?") for x in sib))}
    d23 = (round(float(pdm.distance(tips["REF_Afp2"].taxon, tips["REF_Afp3"].taxon)), 3)
           if "REF_Afp2" in tips and "REF_Afp3" in tips else None)

    res = {
        "tips_total": len(tips), "tips_labelled": len(labelled),
        "tips_unassigned": len(tips) - len(labelled),
        "label_counts": dict(Counter(tip_label[x] for x in labelled)),
        "nearest_labelled_neighbour_concordance": {
            "overall": round(float(nn.concordant.mean()), 3),
            "by_label": {k: round(float(v), 3)
                         for k, v in nn.groupby("label").concordant.mean().items()},
            "discordant": int((~nn.concordant).sum()),
        },
        "mrca_purity": purity,
        "within_vs_between_patristic": dist,
        "reference_third_gene": sister,
        "reference_afp2_afp3_patristic": d23,
    }
    nn.to_csv(OUT / "afp_tree_tips.csv", index=False)
    (OUT / "afp_check3.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
