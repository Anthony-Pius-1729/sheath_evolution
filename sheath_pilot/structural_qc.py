#!/usr/bin/env python3
"""Structural QC of the frozen remodeling metric across held-out systems.

NO evolutionary rates are computed. No homologs are retrieved. This run answers
only whether the frozen R_i is structurally well behaved:

  1. how many residues meet the pre-registered domain |C_E u C_C| >= 3
  2. the distribution of R_i
  3. whether R_i varies continuously rather than collapsing onto 0/1
  4. stability across contact cutoffs and domain thresholds
  5. whether lattice reconstruction is reliable (and how it compares with the
     deposited helical record, which is reported, never used)
  6. whether Diffocin behaves anomalously
  7. how strongly R_i correlates with monomer RSA, packing and ENM flexibility

Separating this from the evolutionary step is deliberate: it makes it impossible
to tune R_i after seeing which definition best predicts evolutionary rate.
"""

from __future__ import annotations

import sys
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
from scipy.stats import spearmanr

import covariates as cov
import lattice as lat
import remodeling as rem

warnings.filterwarnings("ignore")

HERE = Path(__file__).parent
STRUCT = HERE / "structures"
OUT = HERE / "qc_output"

IDENTITY_MIN = 0.90

HELD_OUT = {
    "pvcpnf":    ("Photorhabdus asymbiotica PVCpnf", "6J0B", "6J0C"),
    "r2pyocin":  ("Pseudomonas aeruginosa R2 pyocin", "3J9Q", "3J9R"),
    "salcis":    ("Salmonella salamae SalCIS",        "9R9H", "9R9I"),
    "cis_sc":    ("Streptomyces coelicolor CIS_Sc",   "8BL4", "8BKY"),
    "diffocin":  ("Clostridioides difficile diffocin", "8V3X", "8V3Y"),
    "afp":       ("Serratia entomophila AFP",         "6RBN", "6RC8"),
}
EXPLORATORY = {
    "algocis":   ("Algoriphagus machipongonensis eCIS", "7AE0", "7AEK"),
    "vibrio":    ("Vibrio cholerae T6SS",               "5MXN", "3J9G"),
    "plpvc1":    ("Photorhabdus luminescens PlPVC1",    "9QGL", "9QGP"),
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

ALIGNER = Align.PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.open_gap_score = -11
ALIGNER.extend_gap_score = -1


def fetch(pid):
    STRUCT.mkdir(parents=True, exist_ok=True)
    p = STRUCT / f"{pid}.cif"
    if not p.exists():
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid}.cif", p)
    return p


def load(pid):
    p = fetch(pid)
    return MMCIFParser(QUIET=True).get_structure(pid, str(p))[0], MMCIF2Dict(str(p))


def chain_entities(cif):
    eids = cif.get("_entity.id", [])
    descs = cif.get("_entity.pdbx_description", [])
    out = {}
    for i, eid in enumerate(cif.get("_entity_poly.entity_id", [])):
        desc = descs[eids.index(eid)] if eid in eids else f"entity_{eid}"
        seq = cif["_entity_poly.pdbx_seq_one_letter_code_can"][i].replace("\n", "")
        for ch in cif["_entity_poly.pdbx_strand_id"][i].split(","):
            out[ch.strip()] = {"eid": eid, "desc": desc, "seq": seq}
    return out


def modelled(chain):
    return [r for r in chain if is_aa(r, standard=True) and "CA" in r]


def seq_of(residues):
    return "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in residues)


def identity(a, b):
    if not a or not b:
        return 0.0
    aln = ALIGNER.align(a, b)[0]
    n = sum(1 for (e0, e1), (c0, c1) in zip(aln.aligned[0], aln.aligned[1])
            for o in range(e1 - e0) if a[e0 + o] == b[c0 + o])
    return n / min(len(a), len(b))


def canonical_maps(res_e, res_c):
    """Map author residue numbers in each state onto a shared canonical index
    (the alignment column), so partner residues are directly comparable."""
    aln = ALIGNER.align(seq_of(res_e), seq_of(res_c))[0]
    me, mc = {}, {}
    k = 0
    for (e0, e1), (c0, c1) in zip(aln.aligned[0], aln.aligned[1]):
        for o in range(e1 - e0):
            me[res_e[e0 + o].get_id()[1]] = k
            mc[res_c[c0 + o].get_id()[1]] = k
            k += 1
    return me, mc


def pick_reference(subunits, family):
    """The subunit of `family` with the most complete surrounding shell: it must
    sit in an interior ring, otherwise its contact set is truncated at the end of
    the deposited segment and R_i is measured against missing neighbours."""
    cands = [s for s in subunits if s.family == family]
    if not cands:
        return None
    rings = sorted({s.ring for s in subunits})
    interior = [s for s in cands if rings[0] < s.ring < rings[-1]] or cands
    cents = np.array([s.centroid for s in subunits])
    return max(interior, key=lambda s: int(np.sum(
        np.linalg.norm(cents - s.centroid, axis=1) < 90.0)))


def analyse_system(key, label, pdb_e, pdb_c, report):
    report.append("")
    report.append("=" * 78)
    report.append(f"{label}   extended {pdb_e} / contracted {pdb_c}")
    report.append("=" * 78)

    model_e, cif_e = load(pdb_e)
    model_c, cif_c = load(pdb_c)
    ent_e, ent_c = chain_entities(model_e and cif_e), chain_entities(cif_c)

    chains_e = {ch.id: modelled(ch) for ch in model_e if len(modelled(ch)) >= 20}
    chains_c = {ch.id: modelled(ch) for ch in model_c if len(modelled(ch)) >= 20}

    # ---- cross-state family matching, by sequence ------------------------
    rep_e, rep_c = {}, {}
    for cid, res in chains_e.items():
        e = ent_e.get(cid)
        if e:
            rep_e.setdefault(e["eid"], (cid, res, e["desc"]))
    for cid, res in chains_c.items():
        e = ent_c.get(cid)
        if e:
            rep_c.setdefault(e["eid"], (cid, res, e["desc"]))

    fam_e, fam_c, families = {}, {}, []
    used = set()
    for eid_e, (cid_e, res_e, desc) in sorted(rep_e.items()):
        best, score = None, 0.0
        for eid_c, (cid_c, res_c, _d) in rep_c.items():
            if eid_c in used:
                continue
            s = identity(seq_of(res_e), seq_of(res_c))
            if s > score:
                best, score = eid_c, s
        if best is None or score < IDENTITY_MIN:
            continue
        used.add(best)
        fkey = f"F{len(families) + 1}"
        families.append({"key": fkey, "desc": desc, "identity": score,
                         "eid_e": eid_e, "eid_c": best,
                         "res_e": res_e, "res_c": rep_c[best][1],
                         "cid_e": cid_e, "cid_c": rep_c[best][0]})
        for cid in chains_e:
            if ent_e.get(cid, {}).get("eid") == eid_e:
                fam_e[cid] = fkey
        for cid in chains_c:
            if ent_c.get(cid, {}).get("eid") == best:
                fam_c[cid] = fkey

    if not families:
        report.append("  HALTED: no cross-state matched families")
        return None

    for f in families:
        report.append(f"  {f['key']}: {f['desc'][:46]:48s} identity={f['identity']:.3f}  "
                      f"modelled {len(f['res_e'])} ext / {len(f['res_c'])} con")

    # ---- canonical residue indices ---------------------------------------
    canon = {}
    for f in families:
        canon[f["key"]] = canonical_maps(f["res_e"], f["res_c"])

    # ---- build subunits, expanding only where the deposition needs it -----
    def build(model, cif, fammap, state):
        nch = sum(1 for ch in model if fammap.get(ch.id))
        try:
            olig = int((cif.get("_pdbx_struct_assembly.oligomeric_count") or ["0"])[0])
        except ValueError:
            olig = 0
        need = olig > sum(1 for ch in model if ch.id in fammap) and len(
            lat.assembly_operators(cif)) > 1
        subs, expanded = lat.build_subunits(model, cif, fammap, need)
        which = 0 if state == "ext" else 1
        for s in subs:
            base = s.key.split("#")[0]
            m = canon[s.family][which]
            s.canon = [m.get(r) for r in s.resids]
        return subs, expanded

    subs_e, exp_e = build(model_e, cif_e, fam_e, "ext")
    subs_c, exp_c = build(model_c, cif_c, fam_c, "con")
    report.append(f"  subunits: {len(subs_e)} extended{' (EXPANDED)' if exp_e else ''}, "
                  f"{len(subs_c)} contracted{' (EXPANDED)' if exp_c else ''}")
    if exp_e != exp_c:
        report.append("  WARNING asymmetric expansion: one state generated, the other "
                      "native. This is the defect that invalidated Vibrio.")

    anchor = families[0]["key"]
    try:
        latt_e = lat.derive(subs_e, anchor)
        latt_c = lat.derive(subs_c, anchor)
    except lat.LatticeError as exc:
        report.append(f"  HALTED: {exc}")
        return None

    ref_e = pick_reference(subs_e, anchor)
    ref_c = pick_reference(subs_c, anchor)
    try:
        hint = lat.orient_to(latt_e, ref_e, ref_c)
        latt_c = lat.derive(subs_c, anchor, axis_hint=hint)
        ref_c = pick_reference(subs_c, anchor)
    except lat.LatticeError as exc:
        report.append(f"  HALTED during axis orientation: {exc}")
        return None

    # symmetric ring completion from the measured lattice
    add_e = lat.complete_rings(latt_e)
    add_c = lat.complete_rings(latt_c)
    subs_e, subs_c = latt_e.subunits, latt_c.subunits
    if add_e or add_c:
        report.append(f"  ring completion   : +{add_e} extended, +{add_c} contracted "
                      f"subunits from the coordinate-derived screw "
                      f"(repeat={latt_e.repeat} ring(s))")

    clash_e = lat.count_clashes(subs_e) if exp_e else 0
    clash_c = lat.count_clashes(subs_c) if exp_c else 0

    report.append(f"  lattice extended  : rings={latt_e.n_rings} sizes={latt_e.ring_sizes} "
                  f"C{latt_e.cn} rise={latt_e.rise:.1f} A twist={latt_e.twist:.1f} deg")
    report.append(f"  lattice contracted: rings={latt_c.n_rings} sizes={latt_c.ring_sizes} "
                  f"C{latt_c.cn} rise={latt_c.rise:.1f} A twist={latt_c.twist:.1f} deg")
    ratio = abs(latt_c.rise / latt_e.rise) if latt_e.rise else float("nan")
    report.append(f"  contraction ratio : {ratio:.2f}"
                  + ("   <-- ANOMALOUS, expected ~0.42-0.45" if ratio > 0.55 else ""))
    if exp_e or exp_c:
        report.append(f"  expansion clashes : {clash_e} extended, {clash_c} contracted "
                      f"({'OK' if max(clash_e, clash_c) == 0 else 'FAIL - operators suspect'})")

    report.append("  metadata QC (coordinates vs deposited helical record):")
    for state, L, cif in (("ext", latt_e, cif_e), ("con", latt_c, cif_c)):
        for name, measured, dep, verdict in lat.metadata_qc(L, cif):
            flag = "" if verdict in ("PASS", "n/a") else "   <-- METADATA ERROR"
            report.append(f"      {state} {name:6s} coords={measured:>12s} "
                          f"deposited={dep:>14s}  {verdict}{flag}")

    # ---- contact sets and R_i --------------------------------------------
    rows_all = []
    for f in families:
        fkey = f["key"]
        rref_e = pick_reference(subs_e, fkey)
        rref_c = pick_reference(subs_c, fkey)
        if rref_e is None or rref_c is None:
            continue
        cs_e, drop_e, gen_e = rem.contact_sets(rref_e, subs_e)
        cs_c, drop_c, gen_c = rem.contact_sets(rref_c, subs_c)

        # A reference in a terminal ring has no ring-1 neighbour, so its contact
        # set is truncated by where the deposited segment ends rather than by
        # biology. Truncation is symmetric across states here, but it shrinks the
        # domain and biases R upward by removing shared contacts.
        def shell(sub, subs):
            rings = {x.ring for x in subs}
            return (sub.ring - 1 in rings) and (sub.ring + 1 in rings)
        shell_e, shell_c = shell(rref_e, subs_e), shell(rref_c, subs_c)

        prim = rem.remodeling_index(cs_e, cs_c)
        df = pd.DataFrame(prim)
        df["family"] = fkey
        df["system"] = key

        # monomer covariates, computed in both states, primary = mean
        inv_e = {v: k for k, v in canon[fkey][0].items()}
        inv_c = {v: k for k, v in canon[fkey][1].items()}
        chain_e = next(ch for ch in model_e if ch.id == f["cid_e"])
        chain_c = next(ch for ch in model_c if ch.id == f["cid_c"])
        res_e_l, res_c_l = modelled(chain_e), modelled(chain_c)
        cvs = {}
        for tag, chn, resl in (("e", chain_e, res_e_l), ("c", chain_c, res_c_l)):
            cvs[f"rsa_{tag}"] = cov.monomer_rsa(chn, resl)
            cvs[f"wcn_{tag}"] = cov.monomer_wcn(resl)
            cvs[f"gnm_{tag}"] = cov.monomer_gnm_fluctuation(resl)
        for nm in ("rsa", "wcn", "gnm"):
            ve = [cvs[f"{nm}_e"].get((" ", inv_e.get(c), " "), np.nan) for c in df["canon"]]
            vc = [cvs[f"{nm}_c"].get((" ", inv_c.get(c), " "), np.nan) for c in df["canon"]]
            df[nm] = np.nanmean(np.vstack([ve, vc]), axis=0)

        # Carry author residue numbers and the position in the entity SEQRES into
        # the output. The evolutionary stage needs each structural site expressed
        # in reference-sequence coordinates to locate its column in a homolog
        # alignment; a canonical index alone cannot be placed there.
        df["resid_ext"] = [inv_e.get(c) for c in df["canon"]]
        df["resid_con"] = [inv_c.get(c) for c in df["canon"]]
        seqres = ent_e.get(f["cid_e"], {}).get("seq", "")
        if seqres:
            a2 = ALIGNER.align(seq_of(f["res_e"]), seqres)[0]
            mod2seq = {}
            for (m0, m1), (s0, s1) in zip(a2.aligned[0], a2.aligned[1]):
                for o in range(m1 - m0):
                    mod2seq[f["res_e"][m0 + o].get_id()[1]] = s0 + o + 1
            df["seqres_pos"] = [mod2seq.get(r) if r is not None else None
                                for r in df["resid_ext"]]
        else:
            df["seqres_pos"] = np.nan

        dom = df[df["in_domain"]]
        summ = rem.distribution_summary(dom["R"])
        report.append(f"\n  --- {fkey} ({f['desc'][:40]}) ---")
        report.append(f"      neighbour shell complete  : "
                      f"{'yes' if shell_e else 'NO (terminal ring)'} ext / "
                      f"{'yes' if shell_c else 'NO (terminal ring)'} con")
        report.append(f"      contacts dropped (partner unmatched): "
                      f"{drop_e} ext / {drop_c} con residues")
        report.append(f"      residues with any contact : {len(df)}")
        report.append(f"      in primary domain (>=3)   : {len(dom)}")
        if summ:
            report.append(f"      R: mean={summ['mean']:.3f} median={summ['median']:.3f} "
                          f"sd={summ['sd']:.3f}  q10={summ['q10']:.3f} q90={summ['q90']:.3f}")
            report.append(f"      R distribution: frac<=0.05={summ['frac_zero']:.3f}  "
                          f"frac>=0.95={summ['frac_one']:.3f}  "
                          f"interior={summ['frac_interior']:.3f}"
                          + ("   <-- COLLAPSED to 0/1" if summ['frac_interior'] < 0.25 else ""))

        # sensitivity: cutoff and domain threshold
        line = []
        for cut in rem.SWEEP_CUTOFFS:
            r = rem.remodeling_index(cs_e, cs_c, cutoff=cut)
            v = [x["R"] for x in r if x["in_domain"]]
            line.append(f"{cut:.1f}A:n={len(v)},m={np.mean(v):.2f}" if v else f"{cut:.1f}A:--")
        report.append("      cutoff sweep : " + "  ".join(line))
        line = []
        for mu in rem.SWEEP_MIN_UNION:
            r = rem.remodeling_index(cs_e, cs_c, min_union=mu)
            v = [x["R"] for x in r if x["in_domain"]]
            line.append(f">={mu}:n={len(v)},m={np.mean(v):.2f}" if v else f">={mu}:--")
        report.append("      domain sweep : " + "  ".join(line))

        # robustness: lattice information stripped
        rstrip = rem.remodeling_index(cs_e, cs_c, strip_lattice=True)
        vs = [x["R"] for x in rstrip if x["in_domain"]]
        both = pd.DataFrame(rstrip).set_index("canon")["R"]
        joined = df.set_index("canon")["R"].to_frame("lat").join(both.to_frame("strip")).dropna()
        rho = spearmanr(joined["lat"], joined["strip"]).statistic if len(joined) > 5 else np.nan
        report.append(f"      lattice-stripped: n={len(vs)} mean={np.mean(vs):.3f} "
                      f"  rho(lattice, stripped)={rho:.3f}")

        # covariate correlations
        report.append("      Spearman rho(R, covariate) in domain:")
        for nm in ("rsa", "wcn", "gnm"):
            sub = dom[["R", nm]].dropna()
            if len(sub) > 10:
                rr = spearmanr(sub["R"], sub[nm])
                report.append(f"        {nm.upper():4s} rho={rr.statistic:+.3f} "
                              f"p={rr.pvalue:.2g}  n={len(sub)}")

        # generated-coordinate dependence of the persistent set
        gd = pd.DataFrame(rem.generated_dependence(cs_e, cs_c, gen_e, gen_c))
        if (add_e or add_c) and len(gd):
            both = gd.dropna(subset=["R", "R_native_only"])
            fg = gd["f_generated"].dropna()
            rho_g = (spearmanr(both["R"], both["R_native_only"]).statistic
                     if len(both) > 5 else np.nan)
            report.append(f"      GENERATED-COORDINATE DEPENDENCE")
            report.append(f"        residues whose persistent set uses generated coords: "
                          f"{int((gd['n_inter_generated'] > 0).sum())}/{int((gd['n_inter'] > 0).sum())}")
            report.append(f"        f_generated: mean={fg.mean():.3f} median={fg.median():.3f} "
                          f"frac(f=0)={float((fg == 0).mean()):.3f} "
                          f"frac(f>0.5)={float((fg > 0.5).mean()):.3f}")
            report.append(f"        R with completion   mean={both['R'].mean():.3f} (n={len(both)})")
            report.append(f"        R native-only       mean={both['R_native_only'].mean():.3f}")
            report.append(f"        delta={both['R_native_only'].mean() - both['R'].mean():+.3f}  "
                          f"rho={rho_g:.3f}"
                          + ("   <-- MATERIAL CHANGE" if abs(both['R_native_only'].mean()
                                                             - both['R'].mean()) > 0.10
                             or (rho_g == rho_g and rho_g < 0.8) else "   OK provenance only"))
            df = df.merge(gd[["canon", "f_generated", "R_native_only",
                              "n_inter_generated"]], on="canon", how="left")

        df["identity"] = f["identity"]
        df["contraction_ratio"] = ratio
        rows_all.append(df)

    if not rows_all:
        return None
    out = pd.concat(rows_all, ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / f"R_{key}.csv", index=False)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = [
        "FROZEN REMODELING METRIC -- STRUCTURAL QC",
        "",
        "R_i = 1 - |C_E & C_C| / |C_E u C_C|",
        "C_i^s = {(partner_family, d_ring, d_azimuth, partner_residue)}",
        "",
        f"  primary cutoff      {rem.PRIMARY_CUTOFF} A      sweep {rem.SWEEP_CUTOFFS}",
        f"  primary domain      |union| >= {rem.PRIMARY_MIN_UNION}   sweep {rem.SWEEP_MIN_UNION}",
        "  lattice             coordinate-derived; deposited metadata is QC only",
        "  partner residues    canonical alignment index, both states",
        "",
        "NO evolutionary rates are computed in this run.",
    ]

    frames = []
    for group, systems in (("HELD-OUT", HELD_OUT), ("EXPLORATORY", EXPLORATORY)):
        report.append("")
        report.append("#" * 78)
        report.append(f"# {group}")
        report.append("#" * 78)
        for key, (label, pe, pc) in systems.items():
            try:
                df = analyse_system(key, label, pe, pc, report)
            except Exception as exc:  # noqa: BLE001 - report, never mask
                report.append(f"  ERROR in {key}: {type(exc).__name__}: {exc}")
                continue
            if df is not None:
                df["group"] = group
                frames.append(df)

    if frames:
        allr = pd.concat(frames, ignore_index=True)
        allr.to_csv(OUT / "R_all.csv", index=False)
        dom = allr[allr["in_domain"]]
        report.append("")
        report.append("=" * 78)
        report.append("SUMMARY (primary domain only)")
        report.append("=" * 78)
        g = dom.groupby(["group", "system", "family"])["R"].agg(
            ["count", "mean", "median", "std"]).round(3)
        report.append(g.to_string())

    text = "\n".join(report)
    (OUT / "structural_qc_report.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    sys.exit(main())
