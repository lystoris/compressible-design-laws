#!/usr/bin/env python3
"""verify.py — headline-number reproduction gate.

Runs a checklist of the paper's headline claims against the reproduced canon in
`results/` (and, for the stats-validation decoupling check, the committed audited
golden `tests/fixtures/golden/effdim_grid_t1.csv`). Prints a PASS/FAIL table and
exits non-zero if any check fails.

Note: `verify.py` reads the committed canonical `results/`, so `make verify` on its own
validates the frozen canon. It becomes an end-to-end reproduction check only after
`make reproduce` regenerates `results/` from source (and the slow golden tests in tests/
are where the engines are actually re-run and compared to the audited outputs).

Targets are grounded in docs/reconciliation-notes.md (R-1..R-4): they reflect what
this code release ACTUALLY reproduces, not the paper's original plan placeholders.
Notably:
  - R-3: the reconstructed R5 decoupling grid's eta2_deff (0.647) is checked via a
    dominance test (eta2_deff >> eta2_dnom), not against the paper's 0.73 headline
    eta2 — only the confound-controlled partial correlation (-0.843, matching the
    paper's -0.84) is checked as an exact target for the reconstruction. The
    stats-validation check (against the committed audited golden t1 grid) DOES
    check eta2_deff against 0.748, since that is real audited data, not a
    reconstruction.
  - R-4: Fig6B selection payoff is checked as a pooled-median dominance claim
    (law > black box for N>=100), since the exact per-cell seeds behind the
    original round-07 run are lost.

Run: /usr/bin/python3 verify.py
"""
from __future__ import annotations

import os
import sys

import pandas as pd
from scipy import stats as scipy_stats

from cdl.stats import eta2, partial_spearman

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
GOLDEN = os.path.join(HERE, "tests", "fixtures", "golden")

ROWS = []  # (name, got, expected_str, tol_str, passed)


def check(name, got, expected, tol):
    """Numeric check: PASS if |got - expected| <= tol."""
    passed = abs(got - expected) <= tol
    ROWS.append((name, f"{got:.4f}", f"{expected:.4f}", f"{tol:.4f}", passed))
    return passed


def check_bool(name, got_value, condition_str, passed):
    """Boolean/condition check: PASS if `passed` is True. `got_value` is shown as-is."""
    ROWS.append((name, str(got_value), condition_str, "-", passed))
    return passed


def main():
    all_pass = True

    # --- R1/R3 sweep falloff (results/sweep_grid.csv) ---
    sweep = pd.read_csv(os.path.join(RESULTS, "sweep_grid.csv"))
    med_by_d = sweep.groupby("d")["r2_law"].median()
    for d, expected, tol in [(3, 0.977, 0.05), (9, 0.749, 0.05),
                              (20, 0.436, 0.05), (67, 0.045, 0.05)]:
        got = float(med_by_d.loc[d])
        all_pass &= check(f"sweep median r2_law @ d={d}", got, expected, tol)

    rho_sel, _ = scipy_stats.spearmanr(sweep["r2_law"], 1 - sweep["law_top1_regret"])
    all_pass &= check("sweep selection: spearman(r2_law, 1-law_top1_regret)",
                       float(rho_sel), 0.704, 0.08)

    # --- R3 curated partials (results/summary_t3.json) ---
    with open(os.path.join(RESULTS, "summary_t3.json")) as f:
        import json
        summary_t3 = json.load(f)
    all_pass &= check("curated partial_effd (summary_t3)",
                       summary_t3["partial_effd"], -0.705, 0.08)
    all_pass &= check("curated partial_nomd (summary_t3)",
                       summary_t3["partial_nomd"], 0.601, 0.10)

    # --- R3 ProteinGym partials (results/proteingym_summary_t3.json) ---
    with open(os.path.join(RESULTS, "proteingym_summary_t3.json")) as f:
        pg_summary_t3 = json.load(f)
    all_pass &= check("proteingym partial_effd (proteingym_summary_t3)",
                       pg_summary_t3["partial_effd"], -0.80, 0.05)
    all_pass &= check("proteingym partial_nomd (proteingym_summary_t3)",
                       pg_summary_t3["partial_nomd"], 0.30, 0.08)

    # --- R2 decoupling: stats-validation against the committed audited golden ---
    golden_effdim = pd.read_csv(os.path.join(GOLDEN, "effdim_grid_t1.csv"))
    eta2_deff_g, partial_deff_g = [], []
    for _, sub in golden_effdim.groupby("family"):
        eta2_deff_g.append(eta2(sub.r2_law, sub.d_eff))
        partial_deff_g.append(partial_spearman(sub.r2_law, sub.d_eff, sub.d_nom))
    eta2_deff_g = sum(eta2_deff_g) / len(eta2_deff_g)
    partial_deff_g = sum(partial_deff_g) / len(partial_deff_g)
    all_pass &= check("[stats-validation] golden eta2_deff (per-family-avg)",
                       eta2_deff_g, 0.748, 0.02)
    all_pass &= check("[stats-validation] golden partial_deff (per-family-avg)",
                       partial_deff_g, -0.837, 0.03)

    # --- R2 decoupling: reconstruction against results/effdim_grid_t1.csv ---
    recon_effdim = pd.read_csv(os.path.join(RESULTS, "effdim_grid_t1.csv"))
    eta2_deff_r, eta2_dnom_r, partial_deff_r = [], [], []
    for _, sub in recon_effdim.groupby("family"):
        eta2_deff_r.append(eta2(sub.r2_law, sub.d_eff))
        eta2_dnom_r.append(eta2(sub.r2_law, sub.d_nom))
        partial_deff_r.append(partial_spearman(sub.r2_law, sub.d_eff, sub.d_nom))
    eta2_deff_r = sum(eta2_deff_r) / len(eta2_deff_r)
    eta2_dnom_r = sum(eta2_dnom_r) / len(eta2_dnom_r)
    partial_deff_r = sum(partial_deff_r) / len(partial_deff_r)
    all_pass &= check("[reconstruction] partial_deff (effdim_grid_t1, per-family-avg)",
                       partial_deff_r, -0.843, 0.06)
    dominance = eta2_deff_r > 50 * eta2_dnom_r
    all_pass &= check_bool(
        "[reconstruction] dominance eta2_deff > 50*eta2_dnom",
        f"eta2_deff={eta2_deff_r:.4f}, eta2_dnom={eta2_dnom_r:.4f}",
        "eta2_deff > 50*eta2_dnom",
        dominance,
    )

    # --- R4 panel headline (results/panel_t3.csv) ---
    panel = pd.read_csv(os.path.join(RESULTS, "panel_t3.csv"))
    bc_r2 = float(panel.loc[panel["id"] == "beta-carotene", "r2_law"].iloc[0])
    all_pass &= check("beta-carotene r2_law (panel_t3)", bc_r2, 0.92, 0.04)

    proteingym = pd.read_csv(os.path.join(RESULTS, "proteingym_t3.csv"))
    full_compressors = pd.concat([
        panel.loc[panel["r2_law"] >= 0.7, ["id"]],
        proteingym.loc[proteingym["r2_law"] >= 0.7, ["id"]],
    ])
    n_full = len(full_compressors)
    is_only_bc = n_full == 1 and full_compressors["id"].tolist() == ["beta-carotene"]
    all_pass &= check_bool(
        "full-compressor count (r2_law>=0.7 across 82 landscapes)",
        f"n={n_full}, ids={full_compressors['id'].tolist()}",
        "n==1 and id==beta-carotene",
        is_only_bc,
    )

    # --- R5 simulator (results/panel_t3.csv) ---
    vl_r2 = float(panel.loc[panel["id"] == "vanlent-simulator", "r2_law"].iloc[0])
    all_pass &= check("vanlent-simulator r2_law (panel_t3)", vl_r2, 0.667, 0.03)

    # --- Selection payoff (results/selection_regime.csv) ---
    sel = pd.read_csv(os.path.join(RESULTS, "selection_regime.csv"))
    sel_ge100 = sel[sel["N"] >= 100]
    law_med = float(sel_ge100["law_ov"].median())
    bb_med = float(sel_ge100["bb_ov"].median())
    payoff = law_med > bb_med
    all_pass &= check_bool(
        "selection payoff: pooled median law_ov > bb_ov (N>=100)",
        f"law_ov={law_med:.1f}, bb_ov={bb_med:.1f}",
        "law_ov > bb_ov",
        payoff,
    )

    # --- Print table ---
    headers = ("claim", "got", "expected", "tol", "status")
    widths = [max(len(str(r[i])) for r in ([headers] + ROWS)) for i in range(5)]

    def fmt_row(r):
        status = "PASS" if r[4] is True or r[4] == "PASS" else ("FAIL" if r[4] in (False, "FAIL") else r[4])
        cells = [str(r[0]), str(r[1]), str(r[2]), str(r[3]), status]
        return " | ".join(c.ljust(w) for c, w in zip(cells, widths))

    print(fmt_row(headers))
    print("-+-".join("-" * w for w in widths))
    n_pass = 0
    for r in ROWS:
        status = "PASS" if r[4] else "FAIL"
        n_pass += 1 if r[4] else 0
        print(fmt_row((r[0], r[1], r[2], r[3], status)))

    print()
    print(f"{n_pass}/{len(ROWS)} checks passed")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
