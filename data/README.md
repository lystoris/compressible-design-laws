# Data archive — `compressible-design-laws`

This folder is the **data payload** for the code release accompanying Xu, Tu & Xu,
"Effective dimensionality governs when combinatorial design landscapes compress into interpretable
laws" (*Patterns*). It is **not committed to git** (the repo `.gitignore` excludes `data/**`); it is
uploaded to Zenodo and fetched by `scripts/fetch_data.py` / `make data`.

Upload this whole `data/` folder (or a zip of it) as the Zenodo record, then record the minted DOI in
`data/README.md` and the repo's top-level `README.md`.

## Layout

```
data/
├── cleaned/                              # 13 curated design→phenotype datasets, uniform schema
│   ├── manifest.csv                      #   id, paper, doi, organism, target, phenotype_name, N,
│   │                                     #   nominal_d, encoding, group, source, status, notes
│   ├── SCHEMA.md                         #   column contract for the cleaned tables
│   ├── beta-carotene.csv  tryptophan.csv  limonene.csv  isoprenol.csv
│   ├── jervis-RBS.csv  pca-paperA.csv  vanlent-pCA.csv  smanski-nif.csv
│   ├── poelwijk.csv  gb1.csv  avgfp.csv  aav.csv
│   └── vanlent-simulator.csv             #   enumerated 6^7 simulator ground truth
└── proteingym/
    ├── proteingym_ref.csv                # ProteinGym reference; row filter includes_multiple_mutants
    └── DMS_ProteinGym_substitutions/     # the 69 COMBINATORIAL (multi-mutant) DMS assays used
        └── <DMS_id>.csv                  #   columns used: mutant, DMS_score
```

## Provenance

- **Curated cleaned datasets** were reduced to a uniform schema from the originally published sources
  listed (with DOI/accession) in `cleaned/manifest.csv`. Each table ends in a single continuous
  `phenotype` column; replicate designs aggregated to their mean. See `cleaned/SCHEMA.md`.
- **ProteinGym**: only the 69 assays with `includes_multiple_mutants = true` in `proteingym_ref.csv`
  are included here (the combinatorial subset the population analysis uses). The full 217-assay
  substitution benchmark is available from the official ProteinGym Zenodo record
  **10.5281/zenodo.15293562**; this folder is a faithful subset of it, not a re-derivation.
- The code subsamples each assay to at most 1,500 designs at runtime (seed-locked), so the raw assay
  files are kept unmodified to preserve exact reproducibility.

## Regenerating (optional)

The curated cleaned tables are produced from raw sources by the cleaning pipeline
`paper/datasets/clean_all.py` in the authors' working tree; that pipeline is out of scope for this
code release (data is distributed pre-cleaned here). Provenance for every dataset is in
`cleaned/manifest.csv`.
