# IPO-IRT_Paper

**Supplementary materials for *Interval-Pseudo-Order Item Response Theory Models: Grounding Psychometric Measurement in Representational Theory*.**

This repository contains the R code and supporting materials for reproducing the analyses reported in the manuscript. The proposed interval-pseudo-order item response theory (IPO-IRT) framework distinguishes three latent states—ignorance, partial knowledge, and mastery—and implements two specifications: **IPO-U**, with a common hesitation-band half-width, and **IPO-K**, with a common half-width and normalized item-specific relative sharpness.

## Contents

The repository contains four R scripts supporting the theoretical implementation, simulation study, empirical analyses, and supplementary investigations of IPO-IRT.

- **`IPO-IRT Utilities.R`** — Core implementation of the IPO-IRT framework, including the uniform-threshold (IPO-U) and relative item-sharpness (IPO-K) models. Provides state-probability functions, item characteristic curves, marginal maximum likelihood (MML) and expectation-maximization (EM) estimators, EAP ability scoring, response simulation, information-criterion comparisons, and pairwise InterModel Vigorish (IMV) comparisons.

- **`IPO-IRT Simulation.R`** — Reproduces the Monte Carlo simulation study, including parameter recovery, latent-score recovery, information-criterion model selection, pairwise predictive comparisons, and misspecification analyses involving conventional logistic IRT models. Also includes the supplementary investigation of posterior-predictive versus EAP plug-in IMV protocols. Simulations can be executed in parallel using `future` and `future.apply`. Results are saved by condition, and functions are provided to generate publication figures.

- **`IPO-IRT ICAR CV Protocol.R`** — Implements cross-validated comparisons using the Letter and Number Series and Matrix Reasoning subtests of the International Cognitive Ability Resource (ICAR). Fits IPO-U, IPO-K, Rasch, 2PL, and fixed-guessing 2.5PL models and calculates pairwise IMV under two prediction protocols: posterior-predictive integration and EAP plug-in prediction. Provides the common prediction functions used in the main empirical analysis.

- **`IPO-IRT ICAR Empirical Example.R`** — Reproduces the empirical analyses reported in the manuscript. Includes model estimation with multiple starting values, profile-likelihood analyses of the common half-width, information criteria, likelihood-ratio tests, parameter-identification diagnostics, and five-fold cross-validated IMV comparisons. Generates the empirical results tables and figures.

## Reproducing the simulation study

1. Place `IPO-IRT Utilities.R` and `IPO-IRT Simulation.R` in the same working directory.
2. Set the working directory to the folder containing these scripts.
3. Run the simulation script. Its current executable block calls `run_all_simulations()` with 1,000 replications for the principal conditions, posterior-predictive and EAP plug-in IMV analyses, and the optional prediction-protocol supplement.

**Important:** The script starts the full simulation when sourced because its main `run_all_simulations()` call is not commented out. To inspect functions, conduct a smaller test, or generate figures from saved results without rerunning everything, first comment out that executable call and its immediately following summary calls.

For a smaller test, use a reduced number of replications after loading the functions:

```r
results_test <- run_all_simulations(workers = 2, strategy = "multisession", n_reps = 2,
                                    n_imv_reps = 2, protocol_n_reps = 2, run_protocol_supplement = TRUE,
                                    overwrite = TRUE, results_dir = "IPO_IRT_Simulation_Test")
```

The main simulation varies the calibration sample size (*N* = 250, 500, 1,000, 2,000), the number of items (*J* = 10, 20, 50, 100), and the IPO generating model (IPO-U or IPO-K). Each replication uses an independent test sample for latent-score recovery and predictive comparison. A separate misspecification arm generates data from Rasch, 2PL, and 3PL models. The simulation also examines differences between the two IMV prediction protocols for short tests.

## Reproducing the empirical analyses

The empirical examples use the **International Cognitive Ability Resource (ICAR)** Letter and Number Series and Matrix Reasoning subtests, distributed through the [**Item Response Warehouse (IRW)**][https://itemresponsewarehouse.org/]. The items have six response options; the corresponding chance-response probability is therefore fixed at 1/6. The analyses use the SAPA data and retain respondents with at least one observed response on the corresponding subtest.

### 1. Prepare the repository

Download the R scripts into the same working directory. The scripts expect the utility file to be named `IPO-IRT Utilities.R` and the other scripts to retain the filenames listed above. The empirical scripts obtain the data through `irw::irw_fetch()`. On first use, IRW access may require authentication through [Redivis][https://redivis.com/]. Retrieved datasets are subsequently cached locally to avoid repeated downloads.

### 2. Run the main empirical analysis

Execute:

```r
source("IPO-IRT ICAR Empirical Example.R")
```

This script automatically loads the required functions from `IPO-IRT ICAR CV Protocol.R` and `IPO-IRT Utilities.R`, retrieves the ICAR data, and performs the empirical analyses for both subtests.

The procedure includes fitting IPO-U, IPO-K, Rasch, 2PL, and fixed-guessing 2.5PL models; selecting IPO solutions from multiple starting values; profiling the common hesitation-band half-width; computing AIC, BIC, HQ, and likelihood-ratio statistics; and conducting five-fold cross-validation.

Pairwise IMV comparisons are calculated using both posterior-predictive integration and EAP plug-in prediction. Posterior-predictive IMV is the primary protocol reported in the manuscript.

The analysis uses parallel computation through `future` and `future.apply`. Depending on available hardware, the full analysis can be computationally demanding.

### 3. Retrieve the empirical results

The main script generates the following files:

| File | Description |
|---|---|
| `icar_sapa_cache.rds` | Locally cached ICAR response data |
| `icar_empirical_results.txt` | Detailed estimation results and model diagnostics |
| `icar_empirical_fits.rds` | Saved fitted models, likelihood profiles, and cross-validation summaries |
| `icar_empirical_cv_preds.rds` | Saved fold-level predictions for both IMV protocols |
| `icar_empirical_tables.tsv` | Manuscript Table 1 and Supplementary Table S5 |
| `fig_icar_profile_h.png` / `.pdf` | Profile-likelihood figure for the common half-width |
| `fig_icar_imv.png` / `.pdf` | Cross-validated pairwise IMV figure |

The fitted-model objects and predictions are preserved so that selected analyses and reports can be reproduced without repeating every estimation step.