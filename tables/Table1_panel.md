# Table 1. The compressibility panel

| Dataset | System | Encoding | Nominal d | Effective d | N | Compressibility R2det | Black-box ceiling | Compresses | Domain |
|---|---|---|---|---|---|---|---|---|---|
| beta-carotene (Tu et al., 2026) | Y. lipolytica / β-carotene | 3 integer gene copies | 3 | 1.26 [1.25, 1.26] | 32 | 0.92 [0.92, 0.92] | 0.97 [0.97, 0.97] | yes | positive |
| vanlent-simulator (van Lent et al., 2023) | in silico kinetic / flux | 7 enzyme levels (full 6^7) | 7 | 2.86 [2.85, 2.88] | 2000 | 0.67 [0.63, 0.69] | 0.89 [0.87, 0.90] | partial | simulator |
| jervis-RBS (Jervis et al., 2019) | E. coli / limonene | combinatorial RBS (one-hot) | 14 | 3.38 [3.22, 3.52] | 60 | 0.54 [0.54, 0.54] | 0.58 [0.58, 0.58] | partial | negative |
| vanlent-pCA (van Lent et al., 2025) | S. cerevisiae / p-coumaric acid | 19 genes × promoter strength | 19 | 5.74 [5.71, 5.79] | 403 | 0.49 [0.49, 0.49] | 0.51 [0.51, 0.51] | partial | negative |
| isoprenol (Carruthers et al., 2025) | P. putida / isoprenol | binary multiplex-CRISPRi (≤67) | 121 | 12.10 [12.10, 12.29] | 401 | 0.43 [0.43, 0.43] | 0.45 [0.45, 0.45] | partial | negative |
| limonene (Radivojević et al., 2020) | E. coli / limonene | 9 pathway protein levels | 9 | 5.80 [5.55, 6.17] | 30 | 0.39 [0.39, 0.39] | 0.62 [0.62, 0.62] | no | negative |
| pca-paperA (Moreno-Paz et al., 2024) | S. cerevisiae / p-coumaric acid | 6 categorical pathway factors | 5 | 2.08 [2.07, 2.10] | 91 | 0.38 [0.38, 0.38] | 0.69 [0.69, 0.69] | no | positive |
| poelwijk (Poelwijk et al., 2019) | eqFP611 / brightness | 13 binary positions (complete 2^13) | 13 | 3.27 [3.05, 4.02] | 2000 | 0.32 [0.24, 0.43] | 0.94 [0.93, 0.94] | no | keystone |
| aav (Bryant et al., 2021) | AAV2 capsid / viability | 28-position segment | 28 | 22.98 [22.91, 23.61] | 2000 | 0.22 [0.21, 0.24] | 0.56 [0.53, 0.56] | no | contrast |
| avgfp (Sarkisyan et al., 2016) | avGFP / brightness | ~233 positions | 233 | 84.13 [80.81, 93.18] | 2000 | 0.09 [0.07, 0.10] | 0.32 [0.23, 0.40] | no | positive |
| smanski-nif (Smanski et al., 2014) | E. coli (nif) / nitrogenase | refactored cluster part-slots | 115 | 49.13 [48.68, 52.31] | 80 | 0.04 [0.04, 0.04] | -0.02 [-0.02, -0.02] | no | positive |
| gb1 (Wu et al., 2016) | protein G B1 / fitness | 4 sites × 20 aa | 4 | 3.92 [3.89, 3.93] | 2000 | 0.03 [-0.11, 0.08] | 0.19 [0.13, 0.21] | no | contrast |
| tryptophan (Zhang et al., 2020) | S. cerevisiae / tryptophan | 5 genes × 6 promoters | 5 | 3.38 [3.37, 3.39] | 307 | -0.27 [-0.27, -0.27] | 0.13 [0.13, 0.13] | no | negative |
