# Compressible Design Laws

Reproducibility code for:

> Xu, J., Tu, W. & Xu, J. "Effective dimensionality governs when combinatorial design
> landscapes compress into interpretable laws." *Patterns* (Cell Press), 2026 (in
> submission).

- **Data (Zenodo):** [10.5281/zenodo.21238983](https://doi.org/10.5281/zenodo.21238983)
- **ProteinGym source subset (Zenodo):** [10.5281/zenodo.15293562](https://doi.org/10.5281/zenodo.15293562)
- **Code (GitHub):** [github.com/lystoris/compressible-design-laws](https://github.com/lystoris/compressible-design-laws)

## Overview

Combinatorial design landscapes (e.g. multi-gene metabolic engineering panels, combinatorial
protein DMS assays) are sometimes well described by a low-dimensional, interpretable law and
sometimes are not. This repository provides the analysis engine (`src/cdl/`) and reproduction
scripts (`scripts/`) that measure a design landscape's **compressibility** (how well a sparse,
low-order symbolic law fits it) and its **effective dimensionality** (the number of design
variables that actually carry signal, as opposed to the nominal number of variables), and shows
that effective dimensionality -- not nominal dimensionality -- governs compressibility, across a
panel of curated real datasets, the combinatorial subset of ProteinGym, controlled synthetic
sweeps, and a mechanistic pathway simulator.

## Install

Requires Python >= 3.9.

```bash
pip install -e .
```

For development/tests:

```bash
pip install -e ".[test]"
```

All of the paper's **headline numbers** reproduce in pure Python (`verify.py`; see below) --
no R is required to check the quantitative claims. Rendering the publication-quality **figures**
additionally requires **R (>= 4.2)** with the packages `ggplot2`, `patchwork`, `scales` and
`dplyr` installed.

## Get the data

Data are not stored in this git repository (`data/` is gitignored) -- they are hosted on Zenodo:

```bash
make data
```

This runs `scripts/fetch_data.py`, which downloads and unpacks the Zenodo record
[10.5281/zenodo.21238983](https://doi.org/10.5281/zenodo.21238983) into `data/`. It is safe to
re-run: if `data/cleaned/manifest.csv` is already present it no-ops. If the automatic download
can't run (e.g. no network access), it prints manual download instructions pointing at the same
DOI. See [`data/README.md`](data/README.md) for the full expected layout and per-dataset
provenance.

## Reproduce

```bash
make reproduce   # regenerate results/, figure data, figures, and tables from data/
make verify      # check the paper's 17 headline numbers against results/
```

`make verify` runs `verify.py`, which checks 17 numbers reported in the paper against the
contents of `results/` (and, for the stats-validation decoupling check, an audited golden
fixture) and **exits non-zero if any check fails**. It is the fastest way to confirm the code
release reproduces the paper's quantitative claims.

Other useful targets:

```bash
make figures     # just rebuild figure data + render figures (Fig 2-6, S1-S2)
make test        # fast test suite (excludes slow reproduction tests)
make test-all    # full test suite, including slow reproduction tests
```

## Licence

Code is released under the [MIT licence](LICENSE). Datasets fetched via `make data` carry their
own licences -- see the Zenodo record and `data/cleaned/manifest.csv` for per-dataset provenance
and terms.

## Citation

See [`CITATION.cff`](CITATION.cff).
