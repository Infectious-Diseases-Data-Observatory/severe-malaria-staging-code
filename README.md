# Peripheral Parasite Staging in Severe Falciparum Malaria

## Project overview

This repository contains the R, Python, and Stan code used to investigate the prognostic value of peripheral parasite developmental staging in severe falciparum malaria.

The study uses individual patient data from the AQUAMAT and SEAQUAMAT trials to examine associations between parasite developmental stage and mortality. Hierarchical models were used to estimate mean parasite age and parasite age standard deviation, accounting for the circular 48-hour developmental cycle and sequestration of mature parasites. Plasma PfHRP2 concentrations were used to estimate parasite biomass.

Generalized Optimal Sparse Decision Trees (GOSDT) were used to investigate whether parasite staging improves mortality risk classification alongside conventional clinical and laboratory predictors. Models incorporating staging were compared with models without staging, including through repeated held-out validation.

## Repository structure

The repository is organised into five principal directories:

| Directory | Description |
|---|---|
| `Biomass/` | Estimation and comparison of parasite biomass |
| `Charts/` | Scripts for generating figures and visualisations |
| `Data/` | Data preparation, exploratory analyses, and sequestration weights |
| `GOSDT/` | Mortality classification using optimal sparse decision trees |
| `Staging/` | Parasite age modelling, bimodality screening, and posterior analyses |

Where applicable, `Complete case/` contains analyses restricted to patients with complete data, while `Imputed/` contains analyses using multiply imputed datasets.

### Biomass

The `Biomass/` directory contains scripts for estimating circulating, total, and sequestered parasite biomass using peripheral parasitaemia and plasma PfHRP2 concentrations.

Both `Complete case/code/` and `Imputed/code/` contain the following script:

| Script | Description |
|---|---|
| `PfHRP2_biomass_analysis.R` | Estimates parasite biomass and compares biomass measures between survival outcomes, including sensitivity analyses using alternative PfHRP2 conversion coefficients. |

Generated CSV summaries and figures are excluded from the public repository.

### Charts

The `Charts/` directory contains R Markdown scripts used to generate figures for the manuscript and supplementary material.

Scripts are organised into `Complete case/R Markdown/` and `Imputed/R Markdown/`.

| Script | Description |
|---|---|
| `bar_observed_staging.Rmd` | Compares observed peripheral parasite stage distributions between survivors and patients who died. |
| `bar_observed_vs_predicted.Rmd` | Compares observed peripheral stage proportions with predictions from the hierarchical parasite age model. |
| `heatmap_mu_age.Rmd` | Visualises inferred parasite age distributions across patients in AQUAMAT and SEAQUAMAT. |
| `PfHRP2_vs_parasitaemia.Rmd` | Examines the relationship between plasma PfHRP2 concentration and peripheral parasitaemia. |
| `scatter_mu_sigma.Rmd` | Visualises the relationship between estimated mean parasite age and parasite age standard deviation. |
| `spline_regression.Rmd` | Examines associations between inferred parasite age parameters and mortality using generalised additive models. |
| `stage_proportions_mortality.Rmd` | Examines associations between inferred latent parasite stage proportions and mortality. |

Not all scripts are present in both analysis directories.

Generated figures, rendered Markdown documents, and other outputs are excluded from the public repository.

### Data

The `Data/` directory contains scripts for data preparation, exploratory checks, replication of historical staging findings, and calculation of sequestration weights.

| Script | Description |
|---|---|
| `Heatmap_dat_all.py` | Generates exploratory heatmaps from the original dataset. |
| `Heatmap_dat_all_check.py` | Performs checks associated with the exploratory heatmap analysis. |
| `Linear_regression_dat_all.py` | Examines relationships between selected variables in the original dataset using linear regression. |
| `Linear_regression_check.py` | Performs additional checks of the linear regression analyses. |
| `Replication_staging_mortality.R` | Reassesses historical associations between peripheral parasite developmental stages and mortality in AQUAMAT and SEAQUAMAT. |
| `Sequestration_weights.R` | Calculates sequestration weights for each developmental stage from observed peripheral staging, using tiny rings as the reference stage. |

The sequestration weights were calculated using patients with complete staging observations and a total recorded stage count between 93 and 107, inclusive. The resulting weights were used as fixed inputs to the hierarchical parasite age models.

The original and imputed patient datasets, along with generated results and diagnostic outputs, are excluded from the public repository.

### GOSDT

The `GOSDT/` directory contains R and Python scripts used to construct, analyse, visualise, and validate mortality classification trees.

#### R scripts

The `GOSDT/R/` directory contains:

| Script | Description |
|---|---|
| `Impute_dat_all.R` | Performs multiple imputation of missing clinical and laboratory variables. |
| `Thresholds_for_mortality.R` | Estimates candidate mortality risk thresholds using statistical models. |
| `Mortality_for_thresholds.R` | Examines mortality associated with candidate predictor thresholds. |
| `Original_thresholds_binned.R` | Converts predictors from the complete-case dataset into binary features. |
| `Original_thresholds_binned_no_staging.R` | Prepares binary predictors without staging variables for the complete-case analysis. |
| `Imputed_thresholds_binned.R` | Constructs binary predictors for the multiply imputed datasets. |
| `Original_analysis_script.R` | Analyses decision tree structures and classification results from the complete-case analysis. |
| `Original_analysis_no_staging.R` | Analyses complete-case trees fitted without staging predictors. |
| `Imputed_analysis_script.R` | Analyses decision trees across the multiply imputed datasets. |
| `Imputed_analysis_no_staging.R` | Analyses imputed-data trees fitted without staging predictors. |
| `Original_plotting_script.R` | Produces decision tree visualisations for the complete-case analysis. |
| `Original_plotting_no_staging.R` | Produces complete-case decision tree visualisations without staging predictors. |
| `Imputed_plotting_script.R` | Produces decision tree visualisations for the imputed analyses. |
| `Imputed_plotting_no_staging.R` | Produces imputed-data decision tree visualisations without staging predictors. |

#### Python scripts

The `GOSDT/Python/` directory contains scripts for GOSDT classification and evaluation.

The filenames use `Original` for complete-case analyses and `Imputed` for analyses using multiply imputed datasets.

| Script pattern | Description |
|---|---|
| `Original_gosdt.py`, `Imputed_gosdt.py` | Fit GOSDT models incorporating staging predictors. |
| `*_gosdt_no_staging.py` | Fit GOSDT models excluding staging predictors. |
| `*_gosdt_validation.py` | Evaluate models using a held-out validation partition. |
| `*_gosdt_validation_no_staging.py` | Evaluate held-out models without staging predictors. |
| `*_gosdt_repeated_validation.py` | Perform repeated held-out validation of models incorporating staging. |
| `*_gosdt_repeated_validation_no_staging.py` | Perform repeated held-out validation without staging predictors. |
| `*_gosdt_staging_comparison.py` | Compare classification performance with and without staging predictors. |

The GOSDT analyses used binary candidate predictors derived from clinical, laboratory, and parasite staging variables. Models were fitted separately with and without staging to assess the incremental predictive contribution of parasite staging.

Generated decision trees, patient-level predictions, validation datasets, fitted analysis objects, and performance outputs are excluded from the public repository.

### Staging

The `Staging/` directory contains R and Stan code used to infer parasite age distributions from observed peripheral staging.

#### Stan models

The `Staging/code/Stan/` directory contains:

| File | Description |
|---|---|
| `biomass_hierarchical_unimodal.stan` | Hierarchical unimodal wrapped Student-t model of parasite age, incorporating fixed sequestration weights and a multinomial likelihood for observed stage counts. |
| `biomass_hierarchical_bimodal.stan` | Hierarchical bimodal mixture model of parasite age, incorporating two wrapped Student-t components, mixture weights, and fixed sequestration weights. |

Both models represent the 48-hour parasite developmental cycle and account for the reduced peripheral representation of mature parasite stages.

#### R scripts

The `Staging/Complete case/code/R/` and `Staging/Imputed/code/R/` directories contain corresponding analysis scripts.

| Script | Description |
|---|---|
| `bimodality_screen.R` | Compares unimodal and bimodal parasite age distributions using BIC and classifies patients according to the consistency of evidence for bimodality across three screening settings. |
| `trimodality_screen.R` | Compares bimodal and trimodal parasite age distributions using BIC and assesses the consistency of evidence for trimodality across three screening settings. |
| `hierarchical_main.R` | Fits the hierarchical parasite age models and processes posterior estimates. |
| `export_latent.R` | Exports inferred latent parasite stage proportions and associated model estimates. |
| `traceplots.R` | Generates diagnostic plots for the hierarchical models. |
| `examples.R` | Produces examples of observed staging and inferred parasite age distributions for selected patients. |

The screening procedures use observed staging data, which were not imputed. One copy of each patient's staging observations was retained for the imputed analysis. The resulting parasite age estimates were subsequently attached to the multiply imputed clinical datasets.

The screening procedures perform numerical optimisation for individual patients, while the subsequent hierarchical models are fitted using Stan.

Fitted Stan model objects, posterior summaries, patient-level latent distributions, diagnostic plots, and other generated outputs are excluded from the public repository.

## Software requirements

The project uses R, Python, and Stan. The versions below describe the software environments recorded for this repository.

### R

- **R version:** 4.6.1
- **Operating system:** macOS Golden Gate 27.0.1 (Apple Silicon, aarch64)

Principal R packages:

| Package | Version |
|---|---|
| `data.table` | 1.18.4 |
| `rstan` | 2.32.7 |
| `StanHeaders` | 2.32.10 |
| `posterior` | 1.7.0 |
| `mice` | 3.19.0 |
| `mgcv` | 1.9-4 |
| `ggplot2` | 4.0.3 |
| `scales` | 1.4.0 |
| `knitr` | 1.51 |

Additional R packages may be required by individual scripts.

### Stan

- **Stan version:** 2.32.2
- **R interface:** `rstan` 2.32.7

The hierarchical models were implemented using Stan and fitted through R.

### Python

- **Python version:** 3.13.7

Principal Python packages recorded in the project virtual environment:

| Package | Version |
|---|---|
| `gosdt` | 1.0.4 |
| `scikit-learn` | 1.9.1 |
| `pandas` | 3.0.6 |
| `numpy` | 2.5.3 |
| `scipy` | 1.18.1 |
| `cloudpickle` | 3.1.2 |
| `joblib` | 1.6.0 |

These versions describe the inspected environments and are not necessarily minimum compatible versions. Software versions may have changed during development.

## Data availability

This project uses individual patient data from the AQUAMAT and SEAQUAMAT severe malaria trials.

The underlying patient datasets are not included in this repository because they contain sensitive clinical information and are subject to data access restrictions.

Fitted model objects, patient-level predictions, posterior estimates, and generated outputs have also been excluded where they may contain or reveal restricted information.

The repository provides the analysis code and statistical model definitions, but does not contain all files required for independent reproduction of the results. Access to the underlying data is subject to the relevant data-sharing arrangements.
