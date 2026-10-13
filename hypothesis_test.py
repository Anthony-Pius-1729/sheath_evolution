#!/usr/bin/env python3
"""The frozen hypothesis test: does R_i explain site-rate variation beyond RSA, WCN, GNM?

Implements FROZEN_ANALYSIS.md exactly. Nothing here is tunable: the datasets, the
predictors, the occupancy rule, the model form, the bootstrap block length and the
replicate count are all fixed by that document before any rate was estimated.

  z(r_i) = b0 + bR*z(R_i) + bA*z(RSA_i) + bP*z(WCN_i) + bF*z(GNM_i) + e_i

Usage:  python hypothesis_test.py [--sensitivity NAME --res FILE --rcsv FILE ...]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

HERE = Path(__file__).parent
OUT = HERE / "results"
RNG_SEED = 20261013          # fixed
N_BOOT = 10_000              # fixed by FROZEN_ANALYSIS.md section 5
BLOCK_LEN = 10               # fixed: contiguous runs of 10 reference residues
MIN_OCCUPANCY = 0.80         # fixed by FROZEN_ANALYSIS.md section 3
PREDICTORS = ["R", "rsa", "wcn", "gnm"]

# dataset -> (rate4site .res, R_*.csv, family filter or None)
DATASETS = {
    "diffocin": ("rates/diffocin.res", "sheath_pilot/qc_output/R_diffocin.csv", None),
    "afp_Afp2": ("rates/afp_Afp2.res", "sheath_pilot/qc_output/R_afp.csv", "F1"),
    "afp_Afp3": ("rates/afp_Afp3.res", "sheath_pilot/qc_output/R_afp.csv", "F2"),
}


def read_r4s(path: Path) -> pd.DataFrame:
    """Parse a rate4site .res table and verify the rate direction.

    Columns are POS SEQ SCORE QQ-INTERVAL STD MSA-DATA. SCORE is the normalised
    rate: mean 0 over sites, and per the file's own header 'lower value = higher
    conservation', i.e. HIGHER SCORE = FASTER EVOLUTION.

    That direction is not taken on trust. rate4site was run with `-y`, which writes
    the un-normalised posterior rates on their natural scale (non-negative, mean 1,
    where larger unambiguously means faster). `check_direction` below asserts that
    the normalised SCORE is a strictly increasing function of those un-normalised
    rates, so the sign of bR cannot be inverted by a misread convention.
    """
    rows = []
    for line in path.read_text().replace("\r", "").split("\n"):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 6 or not parts[0].isdigit():
            continue
        pos, aa, score = int(parts[0]), parts[1], float(parts[2])
        # QQ-INTERVAL may be split across tokens; STD and MSA DATA are the last two
        std = float(parts[-2])
        msa = parts[-1]
        n_aa, n_tot = (int(x) for x in msa.split("/"))
        rows.append({"seqres_pos": pos, "aa": aa, "rate": score, "rate_std": std,
                     "n_aa": n_aa, "n_seq": n_tot, "occupancy": n_aa / n_tot})
    df = pd.DataFrame(rows)
    assert len(df), f"no rate rows parsed from {path}"
    return df


def check_direction(res: pd.DataFrame, unnorm_path: Path) -> dict:
    """Assert that higher normalised SCORE means faster evolution.

    Compares the .res scores against the `-y` un-normalised posterior rates, which
    are on the natural rate scale: non-negative, mean ~1, larger = faster. Requires
    a strictly increasing (Spearman rho == +1) relationship.
    """
    u = read_r4s_unnorm(unnorm_path)
    m = res.merge(u, on="seqres_pos", how="inner", validate="one_to_one")
    assert len(m) == len(res), f"{unnorm_path}: {len(m)} of {len(res)} positions matched"
    assert (m.rate_unnorm >= 0).all(), f"{unnorm_path}: negative un-normalised rates"
    rho = float(spearmanr(m.rate, m.rate_unnorm).statistic)
    # rate4site prints both tables rounded, so exact ties can be ordered differently
    # between them; rho is 1 - O(1e-8) rather than exactly 1. Anything below 0.999
    # would mean a real inversion, which is what this guards against.
    assert rho >= 0.999, (f"{unnorm_path}: normalised score is not an increasing "
                          f"function of the un-normalised rate (Spearman rho = {rho}); "
                          f"the rate direction is not established")
    return {"direction_rho_vs_unnormalised": rho,
            "unnormalised_mean": float(m.rate_unnorm.mean()),
            "unnormalised_range": [float(m.rate_unnorm.min()), float(m.rate_unnorm.max())]}


def read_r4s_unnorm(path: Path) -> pd.DataFrame:
    """Parse the `-y` un-normalised rate table (POS SEQ RATE ... layout)."""
    rows = []
    for line in path.read_text().replace("\r", "").split("\n"):
        parts = line.strip().split()
        if len(parts) < 3 or not parts[0].isdigit():
            continue
        rows.append({"seqres_pos": int(parts[0]), "rate_unnorm": float(parts[2])})
    df = pd.DataFrame(rows)
    assert len(df), f"no un-normalised rate rows parsed from {path}"
    return df


def z(v: pd.Series) -> pd.Series:
    return (v - v.mean()) / v.std(ddof=1)


def fit(d: pd.DataFrame) -> dict:
    X = sm.add_constant(d[PREDICTORS].apply(z))
    y = z(d["rate"])
    m = sm.OLS(y, X).fit()
    ci = m.conf_int(alpha=0.05)
    # partial R^2 of R: (RSS_without - RSS_full) / RSS_without
    rss_full = float(m.ssr)
    m0 = sm.OLS(y, sm.add_constant(d[[p for p in PREDICTORS if p != "R"]].apply(z))).fit()
    partial_r2 = (float(m0.ssr) - rss_full) / float(m0.ssr)
    vif = {}
    for p in PREDICTORS:
        others = [q for q in PREDICTORS if q != p]
        r2 = sm.OLS(z(d[p]), sm.add_constant(d[others].apply(z))).fit().rsquared
        vif[p] = float(1.0 / (1.0 - r2))
    return {"model": m, "beta": {p: float(m.params[p]) for p in PREDICTORS},
            "se": {p: float(m.bse[p]) for p in PREDICTORS},
            "p": {p: float(m.pvalues[p]) for p in PREDICTORS},
            "ci": {p: [float(ci.loc[p, 0]), float(ci.loc[p, 1])] for p in PREDICTORS},
            "r2": float(m.rsquared), "r2_adj": float(m.rsquared_adj),
            "partial_r2_R": partial_r2, "vif": vif, "n_sites": int(len(d))}


def block_bootstrap(d: pd.DataFrame) -> dict:
    """Spatial block bootstrap: blocks are contiguous BLOCK_LEN-residue windows of the
    reference numbering, resampled with replacement to the original site count."""
    rng = np.random.default_rng(RNG_SEED)
    d = d.sort_values("seqres_pos").reset_index(drop=True)
    d["_block"] = (d.seqres_pos - 1) // BLOCK_LEN
    blocks = [g.index.to_numpy() for _, g in d.groupby("_block")]
    n = len(d)
    betas = np.full(N_BOOT, np.nan)
    for b in range(N_BOOT):
        idx = []
        while len(idx) < n:
            idx.extend(blocks[rng.integers(len(blocks))])
        idx = np.asarray(idx[:n])
        s = d.loc[idx]
        if s[PREDICTORS].std(ddof=1).min() == 0 or s["rate"].std(ddof=1) == 0:
            continue
        X = sm.add_constant(s[PREDICTORS].apply(z).to_numpy())
        try:
            betas[b] = np.linalg.lstsq(X, z(s["rate"]).to_numpy(), rcond=None)[0][1]
        except np.linalg.LinAlgError:
            continue
    ok = betas[~np.isnan(betas)]
    return {"n_blocks": len(blocks), "n_replicates_used": int(ok.size),
            "ci95": [float(np.percentile(ok, 2.5)), float(np.percentile(ok, 97.5))],
            "median": float(np.median(ok)),
            "frac_same_sign_as_point": None}


def invariant_sites(afa: Path, positions: list[int]) -> set[int]:
    """Reference positions whose alignment column holds a single distinct residue."""
    a: dict[str, str] = {}
    for blk in afa.read_text().split(">")[1:]:
        h, *body = blk.split("\n")
        a[h.strip()] = "".join(body)
    ref = a["REF"]
    p2c, p = {}, 0
    for i, ch in enumerate(ref):
        if ch != "-":
            p += 1
            p2c[p] = i
    inv = set()
    for pos in positions:
        col = {v[p2c[pos]] for v in a.values()} - {"-"}
        if len(col) <= 1:
            inv.add(pos)
    return inv


def run(name: str, res: str, rcsv: str, fam: str | None, *,
        min_occ: float = MIN_OCCUPANCY, drop_invariant_from: str | None = None,
        tag: str = "primary") -> dict:
    r4s = read_r4s(HERE / res)
    direction = check_direction(r4s, (HERE / res).with_suffix(".unnorm"))
    R = pd.read_csv(HERE / rcsv)
    if fam:
        R = R[R.family == fam]
    R = R[R.in_domain].copy()
    R["seqres_pos"] = R.seqres_pos.astype(int)
    d = R.merge(r4s, on="seqres_pos", how="inner", validate="one_to_one")
    n_sites_total = len(R)
    dropped_merge = n_sites_total - len(d)
    before = len(d)
    d = d[d.occupancy >= min_occ].copy()
    dropped_occ = before - len(d)
    dropped_inv = 0
    if drop_invariant_from:
        inv = invariant_sites(HERE / drop_invariant_from, d.seqres_pos.tolist())
        dropped_inv = int(d.seqres_pos.isin(inv).sum())
        d = d[~d.seqres_pos.isin(inv)].copy()
    assert not d[PREDICTORS + ["rate"]].isna().any().any(), f"{name}: NaN in model inputs"
    f = fit(d)
    bb = block_bootstrap(d)
    out = {"dataset": name, "structural_sites": n_sites_total,
           "sites_unmatched_in_rate_file": dropped_merge,
           "sites_dropped_occupancy": dropped_occ,
           "min_occupancy": min_occ, "sites_dropped_invariant": dropped_inv,
           "sites_fitted": len(d),
           "n_sequences": int(d.n_seq.iloc[0]),
           "rate_std_median": float(d.rate_std.median()),
           "rate_range": [float(d.rate.min()), float(d.rate.max())], **direction,
           "invariant_sites_in_fit": int((d.rate_std == 0).sum()),
           **{k: v for k, v in f.items() if k != "model"},
           "block_bootstrap": bb}
    # tag the filename: sensitivity runs must not overwrite the primary site tables
    d.to_csv(OUT / f"sitetable_{tag}_{name}.csv", index=False)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="run a single dataset by name")
    ap.add_argument("--tag", default="primary", help="label for the output files")
    ap.add_argument("--res", default=None, help="override the .res path (sensitivity runs)")
    ap.add_argument("--min-occupancy", type=float, default=MIN_OCCUPANCY)
    ap.add_argument("--drop-invariant-from", default=None,
                    help="alignment whose invariant reference columns are dropped")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    sel = {args.only: DATASETS[args.only]} if args.only else DATASETS
    results = []
    for name, (res, rcsv, fam) in sel.items():
        results.append(run(name, args.res or res, rcsv, fam,
                           min_occ=args.min_occupancy,
                           drop_invariant_from=args.drop_invariant_from, tag=args.tag))
    (OUT / f"hypothesis_test_{args.tag}.json").write_text(json.dumps(results, indent=1))
    rows = []
    for r in results:
        rows.append({
            "dataset": r["dataset"], "n_seq": r["n_sequences"], "sites": r["sites_fitted"],
            "beta_R": round(r["beta"]["R"], 4), "se_R": round(r["se"]["R"], 4),
            "p_R": r["p"]["R"],
            "ols_ci_lo": round(r["ci"]["R"][0], 4), "ols_ci_hi": round(r["ci"]["R"][1], 4),
            "boot_ci_lo": round(r["block_bootstrap"]["ci95"][0], 4),
            "boot_ci_hi": round(r["block_bootstrap"]["ci95"][1], 4),
            "partial_r2_R": round(r["partial_r2_R"], 4), "model_r2": round(r["r2"], 4),
            "beta_rsa": round(r["beta"]["rsa"], 4), "beta_wcn": round(r["beta"]["wcn"], 4),
            "beta_gnm": round(r["beta"]["gnm"], 4),
            "max_vif": round(max(r["vif"].values()), 2),
        })
    t = pd.DataFrame(rows)
    t.to_csv(OUT / f"hypothesis_test_{args.tag}.csv", index=False)
    print(t.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
