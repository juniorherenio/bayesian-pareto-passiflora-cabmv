# Bayesian multi-trait selection with uncertainty-aware Pareto ranking under genotype × season interaction and CABMV pressure in passion fruit (*Passiflora* spp.)

[![DOI](https://doi.org/10.5281/zenodo.20328435)

> **Gonçalves Júnior, D.H., Albuquerque, D.P., Viana, A.P., Cerri, R.** (2026). Bayesian multi-trait selection with uncertainty-aware Pareto ranking under genotype × season interaction and CABMV pressure in passion fruit hybrids. *Scientific Reports*. [DOI pending]

---

## Overview

This repository contains all data and R scripts used in the above manuscript. The study proposes a novel analytical framework for parent selection in *Passiflora* spp. breeding programs under chronic *Cowpea aphid-borne mosaic virus* (CABMV) pressure, integrating:

- Bayesian univariate mixed models with heterogeneous G×E variance structure (MCMCglmm)
- Proportional Odds Model (CLMM) for individual-level disease pressure characterisation
- Bayesian Disease Severity Index for Plant Breeding (BDSIP)
- Uncertainty-aware Pareto frontier optimisation over MCMC posterior samples
- Unsupervised probabilistic clustering (UMAP + GMM)

The framework is designed to be transferable to other pathosystem × crop combinations with longitudinal ordinal disease data and repeated-measures phenotyping at the individual level.

---

## Repository structure

```
bayesian-pareto-passiflora-cabmv/
├── data/
│   ├── raw/
│   │   └── dataset_completo.csv          # Full phenotypic dataset with repeated measures
│   ├── processed/
│   │   ├── pom_indices.csv               # POM-derived individual disease pressure indices
│   │   ├── aacpd_tres_versoes.csv        # Three AUDPC versions per individual × harvest date
│   │   ├── bdsip_parametros.csv          # Raw BDSIP parameters (8 per individual)
│   │   └── bdsip_completo.csv            # BDSIP parameters + standardised scores + composite index
│   └── results/
│       ├── params_GxA_POM.csv            # Bayesian univariate model parameters (7 traits)
│       ├── h2_estacao_GxA_POM.csv        # Narrow-sense heritability by trait × season
│       ├── rg_bivariado.csv              # Genetic correlations from bivariate models
│       ├── ranking_genitores.csv         # Resilience index and ranking (95 individuals)
│       ├── pareto_bayesiano.csv          # Bayesian Pareto frontier results
│       ├── genitores_recomendados.csv    # Recommended parents (consensus Pareto)
│       ├── clusters_genitores.csv        # GMM cluster assignments with UMAP coordinates
│       └── resumo_clusters.csv           # Cluster summary statistics
├── scripts/
│   ├── 00_setup.R                        # R environment setup and package loading
│   ├── 01_data_preparation.R             # Data wrangling and repeated-measures structure
│   ├── 02_audpc_calculation.R            # Partial accumulated AUDPC (three versions)
│   ├── 03_proportional_odds_model.R      # CLMM for ordinal CABMV severity (ordinal package)
│   ├── 04_bdsip_index.R                  # BDSIP: 8-parameter dynamic severity index
│   ├── 05_disease_metric_comparison.R    # DIC comparison: AUDPC vs BDSIP vs POM
│   ├── 06_bayesian_univariate_GxE.R      # Bayesian univariate models with idh G×E (MCMCglmm)
│   ├── 07_genetic_correlations_bivariate.R  # Bivariate Bayesian models for genetic correlations
│   ├── 08_resilience_index.R             # Resilience index computation (IR_commercial, IR_resistance, IR_final)
│   ├── 09_pareto_frontier.R              # Point-estimate and Bayesian Pareto frontier (NSGA-II)
│   ├── 10_umap_gmm_clustering.R          # UMAP dimensionality reduction + GMM clustering
│   └── 11_figures.R                      # All manuscript figures
└── README.md
```

---

## Data description

### `data/raw/dataset_completo.csv`
Full phenotypic dataset. Each row is one individual × harvest date observation (203 observations, 95 individuals, repeated-measures structure).

| Column | Description | Unit |
|---|---|---|
| `individuo` | Individual plant ID (globally unique) | — |
| `bloco` | Block (randomised complete block design, 4 blocks) | — |
| `familia` | Full-sib family (1 = 293×355; 2 = 501×355; 3 = 516×355) | — |
| `cruzamento` | Cross description | — |
| `P1` | Maternal parent | — |
| `P2` | Paternal parent (common tester: 355) | — |
| `data` | Harvest date | date |
| `estacao` | Season (inverno/primavera/verao/outono) | — |
| `ano` | Year | — |
| `estacao_ano` | Season-year (5 categories) | — |
| `n_frutos` | Number of harvested fruits | count |
| `comp_mm` | Fruit length | mm |
| `diam_mm` | Fruit diameter | mm |
| `ind_formato` | Shape index (length/diameter) | dimensionless |
| `massa_fruto_g` | Mean fruit mass | g |
| `massa_polpa_g` | Pulp mass (raw pulp with seeds and aril) | g |
| `rend_polpa_pct` | Pulp yield | % |
| `esp_casca_mm` | Rind thickness (mean of 4 points) | mm |
| `sst_brix` | Total soluble solids | °Brix |
| `aacpd_planta` | Partial accumulated AUDPC — whole plant | — |
| `aacpd_folha` | Partial accumulated AUDPC — young leaves | — |
| `aacpdm` | Mean AUDPC (mean of plant and leaf) | — |
| `parcela` | Plot ID (family × block) | — |

---

### `data/processed/pom_indices.csv`
Individual-level disease pressure indices derived from the Proportional Odds Model (CLMM).

| Column | Description |
|---|---|
| `individuo` | Individual plant ID |
| `P_severo_media` | Mean posterior probability of severe disease P(score ≥ 3) across all assessments |
| `P_severo_pico` | Peak posterior probability of severe disease |
| `P_severo_final` | Final-assessment posterior probability |
| `P_severo_var` | Temporal variance of posterior severity probability |
| `n_aval` | Number of assessments |
| `familia` | Full-sib family |
| `morreu_fusarium` | Death by Fusarium spp. (1 = yes, 0 = no; right-censored in survival analysis) |
| `POM_idx` | Composite POM index (standardised) |

---

### `data/processed/aacpd_tres_versoes.csv`
Three AUDPC versions per individual × harvest date, used for disease metric comparison (Table 4).

| Column | Description |
|---|---|
| `individuo` | Individual plant ID |
| `data_coleta` | Harvest date |
| `aacpd_acum` | Accumulated AUDPC up to harvest date (trapezoidal method) |
| `aacpd_periodo` | Period AUDPC (assessments within the harvest season) |
| `aacpd_jan75` | 75-day windowed AUDPC (fruit development period) |

---

### `data/processed/bdsip_parametros.csv`
Raw BDSIP parameters (8 per individual), without standardisation.

| Column | Description | Unit |
|---|---|---|
| `individuo` | Individual plant ID | — |
| `mu_max` | Severity peak | score (1–4) |
| `t_pico` | Time to peak severity | days |
| `slope_up` | Ascending-phase progression rate | score/day |
| `slope_down` | Descending-phase rate | score/day |
| `t_surv_obs` | Observed survival time to score ≥ 3 | days |
| `evento` | Event indicator (1 = reached score ≥ 3; 0 = censored) | — |
| `area_beta` | Area under curve on Beta scale (Smithson & Verkuilen 2006) | — |
| `estab` | Temporal stability (variance of scores) | — |
| `prop_acima` | Proportion of assessments with score ≥ 3 | — |
| `n_aval` | Number of assessments | — |
| `familia` | Full-sib family | — |

---

### `data/processed/bdsip_completo.csv`
All BDSIP parameters plus standardised scores and composite BDSIP index.

Includes all columns from `bdsip_parametros.csv` plus standardised versions (suffix `_sc`) and the composite `BDSIP` score (mean of five standardised components: `mu_max_sc`, `-t_surv_sc`, `area_beta_sc`, `slope_up_sc`, `estab_sc`).

---

### `data/results/params_GxA_POM.csv`
Bayesian univariate model parameters for the 7 evaluated traits.

| Column | Description |
|---|---|
| `trait` | Trait name |
| `h2_pond` | Observation-weighted narrow-sense heritability |
| `beta_POM` | Regression coefficient on POM covariate |
| `P_beta_neg` | Posterior probability P(β < 0) |
| `DIC` | Deviance Information Criterion |
| `ESS_min` | Minimum effective sample size across variance components |

---

### `data/results/h2_estacao_GxA_POM.csv`
Narrow-sense heritability by trait × season-year.

| Column | Description |
|---|---|
| `trait` | Trait name |
| `estacao` | Season-year |
| `h2_media` | Posterior mean heritability |
| `h2_025` | 2.5th percentile of posterior distribution |
| `h2_975` | 97.5th percentile of posterior distribution |

---

### `data/results/rg_bivariado.csv`
Genetic correlations from bivariate Bayesian models (7 trait pairs).

| Column | Description |
|---|---|
| `par` | Trait pair label |
| `trait1` | First trait |
| `trait2` | Second trait |
| `rg` | Posterior mean genetic correlation |
| `rg_025` | 2.5th percentile |
| `rg_975` | 97.5th percentile |
| `P_pos` | Posterior probability P(rg > 0) |
| `ESS` | Effective sample size of rg chain |
| `interpretacao` | Biological interpretation |

---

### `data/results/ranking_genitores.csv`
Resilience index components and final ranking for all 95 individuals.

| Column | Description |
|---|---|
| `individuo` | Individual plant ID |
| `familia` | Full-sib family |
| `cruzamento` | Cross description |
| `IR_comercial` | Commercial resilience index |
| `IR_res_sc` | Resistance index (standardised) |
| `IR_estabilidade` | Stability component |
| `IR_final` | Final resilience index |
| `rank_final` | Final ranking position |
| `rank_comercial` | Commercial ranking position |
| `rank_resistencia` | Resistance ranking position |
| `P_severo_media` | Mean POM severity probability |
| `estab_media` | Mean stability |
| `var_pred_medio` | Mean prediction variance |

---

### `data/results/pareto_bayesiano.csv`
Bayesian Pareto frontier results (86 individuals with POM data).

| Column | Description |
|---|---|
| `individuo` | Individual plant ID |
| `P_pareto` | Posterior probability of being Pareto-optimal (500 MCMC samples) |
| `contagem` | Count of samples in which individual was on the Pareto frontier |
| `familia` | Full-sib family |
| `cruzamento` | Cross description |
| `IR_comercial` | Commercial index |
| `IR_res_sc` | Resistance index (standardised) |
| `IR_final` | Final resilience index |
| `P_severo_media` | Mean POM severity probability |
| `var_pred_medio` | Mean prediction variance |

---

### `data/results/clusters_genitores.csv`
GMM cluster assignments with UMAP coordinates for all 95 individuals.

| Column | Description |
|---|---|
| `individuo` | Individual plant ID |
| `cluster` | GMM cluster assignment (1, 2, or 3) |
| `entropia` | Posterior entropy of cluster assignment (0 = certain) |
| `umap1` | UMAP dimension 1 |
| `umap2` | UMAP dimension 2 |
| `familia` | Full-sib family |
| `IR_comercial` | Commercial index |
| `IR_res_sc` | Resistance index |
| `IR_final` | Final index |
| `P_severo_media` | Mean POM severity probability |
| `var_pred_medio` | Mean prediction variance |
| `P_pareto` | Posterior Pareto-optimal probability |

---

## Computational requirements

| Software | Version |
|---|---|
| R | 4.5.3 |
| MCMCglmm | — |
| ordinal | — |
| uwot | — |
| mclust | — |
| survival | — |
| betareg | — |
| ggplot2 | — |

**MCMC performance note:** Chains were executed sequentially (not in parallel) to preserve OpenBLAS matrix-level parallelisation. On the test machine (Intel i5-12th gen, 16 GB RAM, Windows 11, OpenBLAS 0.3.33), each univariate chain required approximately 15–20 minutes. Using `makeCluster` parallelisation on Windows causes child processes to revert to the default BLAS, negating the speedup.

---

## Reproducibility

To reproduce all analyses:

1. Clone this repository
2. Open `scripts/00_setup.R` and install required packages
3. Run scripts in numerical order (01 through 11)
4. All intermediate outputs are saved to `data/processed/` and `data/results/`
5. Figures are saved by `scripts/11_figures.R`

---

## License

- **Code** (`/scripts`): [MIT License](LICENSE)
- **Data** (`/data`): [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)

---

## Citation

If you use this code or data, please cite:

> Gonçalves Júnior, D.H., Albuquerque, D.P., Viana, A.P., Cerri, R. (2026). Bayesian multi-trait selection with uncertainty-aware Pareto ranking under genotype × season interaction and CABMV pressure in passion fruit hybrids. *Scientific Reports*. [DOI pending]

---

## Contact

Deurimar Herênio Gonçalves Júnior  
Postdoctoral Researcher — Universidade Federal do Espírito Santo (UFES)  
[your email here]
