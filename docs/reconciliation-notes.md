# Reconciliation notes — code release vs manuscript

Deltas found while making the code reproduce the paper's reported numbers. Policy: **the paper's
reported numbers are ground truth** — the code is made to reproduce them, and any inconsistency in the
manuscript *text* (as opposed to its reported values) is flagged here for a tracked-changes correction
in `paper/Patterns/Xu_etal_Patterns_Manuscript.docx`.

## R-1. Panel subsample cap: Methods says 1,500; the actual runs used 2,000

**Status:** code aligned to 2,000 (Task 5); **manuscript text correction pending** (Task 17).

- The published `panel_t3.csv` golden and all three trajectory summaries (`summary_t{1,2,3}.json`,
  round-08) were generated with a **max_n = 2,000** subsample cap. Confirmed: every trajectory's
  `panel_*.csv` has `max(N) = 2000`.
- Those summaries give curated `partial_effd` = −0.724 / −0.622 / −0.705 (t1/t2/t3), median
  **−0.705, range [−0.62, −0.72]**, and `partial_nomd` median **+0.55 [+0.51, +0.60]** — which match
  the manuscript's reported values **exactly**. So the paper's headline curated partials are the
  2,000-cap numbers.
- The manuscript **Methods** ("Design encoding and dimensionality") says: *"For datasets larger than
  1,500 measured designs, a random subsample of 1,500 was drawn."* This **1,500 is stale**; the runs
  behind every figure and reported number used 2,000.
- Verification: at `max_n = 2000`, the four subsample-sensitive datasets reproduce the golden
  bit-for-bit — gb1 −0.114, aav 0.214, avgfp 0.102, poelwijk 0.318 (r2_law), with eff_d matching to
  the reported decimal. Datasets with N < 1,500 reproduce identically at any cap ≥ their size.

**Resolution:**
1. `scripts/run_panel.py` default `--max-n` set to **2000** (reproduces golden + reported partials).
2. **Manuscript:** correct Methods "a random subsample of 1,500 was drawn" → **"2,000"** (tracked
   change, author "Jiabao Xu"), and mirror into `manuscript-master.md`. Gated into Task 17.

## R-2. `mva-243` manifest row has no data file

**Status:** flagged; manifest cleanup pending.

- `data/cleaned/manifest.csv` lists `mva-243` (Mukherjee 2022) with status `needs-review`, but there
  is no `mva-243.csv` — it is the optional dataset that never arrived (MANIFEST "MVA-243 table … add
  if supplements arrive"). `run_panel` prints `[err] mva-243: No such file` and skips it, leaving the
  correct **13** curated datasets, so results are unaffected.
- **Resolution:** drop the `mva-243` row from the released `manifest.csv` (it is not part of the
  published panel), so a reproducer sees no error. To be done when the canonical `data/` manifest is
  finalized for the Zenodo upload.
