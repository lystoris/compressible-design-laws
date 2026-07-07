.PHONY: data reproduce figures verify test test-all

# Download and unpack the Zenodo data archive (10.5281/zenodo.21238983) into data/.
# Safe to re-run: no-ops if data/cleaned/manifest.csv is already present.
data:
	/usr/bin/python3 scripts/fetch_data.py

# Full reproduction: results/ -> figure data -> figures -> tables.
# Requires data/ to be populated first (`make data`).
reproduce:
	/usr/bin/python3 scripts/run_panel.py --traj t3
	/usr/bin/python3 scripts/run_proteingym.py --traj t3
	/usr/bin/python3 scripts/run_sweep.py --traj t3
	/usr/bin/python3 scripts/run_decoupling.py --traj t1
	/usr/bin/python3 scripts/run_simulator.py
	/usr/bin/python3 scripts/build_figure_data.py
	Rscript figures/R/make_all.R
	if [ -f tables/make_tables.py ]; then /usr/bin/python3 tables/make_tables.py; fi

# Just the figures (assumes results/ already populated).
figures:
	/usr/bin/python3 scripts/build_figure_data.py
	Rscript figures/R/make_all.R

# Headline-number reproduction gate: checks 17 numbers reported in the paper
# against results/, exits non-zero on any mismatch.
verify:
	/usr/bin/python3 verify.py

# Fast test suite (excludes slow reproduction tests).
test:
	/usr/bin/python3 -m pytest -q -m "not slow"

# Full test suite, including slow reproduction tests.
test-all:
	/usr/bin/python3 -m pytest -q
