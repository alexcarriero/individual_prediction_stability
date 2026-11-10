# IndivPredStability

[![arXiv](https://img.shields.io/badge/arXiv-2610.33744-b31b1b.svg)](https://arxiv.org/abs/2610.33744)

Code to reproduce the simulation study in:

> Carriero A\*, van Dam MLG\*, de Hond AAH, Riley RD, Collins GS, Moons KGM, van Smeden M, Erler NS.
> **Estimating Stability of Individual-Specific Predictions from AI-based prediction models using Bootstrap and Cross-validation Resampling.**
> *arXiv preprint* arXiv:2610.33744 [stat.ME] (2026). https://arxiv.org/abs/2610.33744
>
> \* These authors contributed equally to this work.

---

## Overview

This repository contains everything needed to reproduce the results of [*Estimating Stability of Individual-Specific Predictions from AI-based prediction models using Bootstrap and Cross-validation Resampling*](https://arxiv.org/abs/2610.33744): the data-generating mechanisms (DGMs), the model-development and resampling procedures, the simulation scripts used on the HPC, the scripts that calculate the performance measures, and the code for the tables and figures. All data are simulated, so no external data are needed.

### Simulation scenarios

| | |
|---|---|
| **Full factorial design** | 2 DGMs × 2 sample sizes = 4 scenarios (Table 1) |
| **Data-generating mechanisms (DGMs)** | **(1) Linear:** logistic DGM with w₀ = −1.4649 and w₁ = 0.1610: six weak effects (w₁), three strong effects (3w₁) and three null effects, tuned to AUC ≈ 0.80.<br><br>**(2) Non-linear:** a random forest trained on a preparatory dataset (n = 10,000) generated with interactions and higher-order terms (Eq. 3) |
| **Sample sizes** | N = 458 (minimum required sample size from `pmsampsize`) and N = 229 (half of that) |
| **Initial datasets** | 100 per scenario |

### Per initial dataset

| Step | Description |
|---|---|
| **1. Original models** | Fit the two models, elastic net logistic regression (`glmnet`) and random forest (`ranger`), to the initial dataset, and obtain an in-sample original prediction for every individual |
| **2. Reference** | Independent sampling: draw 200 new datasets of size N from the DGM, refit both models to each, and predict for every individual in the initial dataset |
| **3. Resampling methods** | Non-parametric bootstrap, parametric bootstrap, out-of-bag bootstrap, and repeated 10-, 5- and 2-fold cross-validation. Each method gives 200 predictions per individual and model |
| **4. Performance measures** | Per individual: prediction variance, empirical 95% range width, MAPE (vs. the original prediction), tMAPE (vs. the true DGM probability), and 2.5th and 97.5th percentiles. Each resampling method is then compared with the reference (signed difference and absolute relative deviation). Initial models are also assessed with AUC and calibration |

---

## Software

The simulation was run in **R 4.6.0** on the high-performance computing (HPC) cluster at UMC Utrecht (Rocky Linux 8.10, SLURM), using these packages:

- **Simulation:** `glmnet` 5.0, `ranger` 0.18.0, `MASS` 7.3-65, `pROC` 1.19.0.1, `tidyverse` 2.0.0 (`dplyr` 1.2.1, `tidyr` 1.3.2, `purrr` 1.2.2, `ggplot2` 4.0.3)
- **Processing results:** `arrow`, `data.table`, `dplyr`, `here`
- **Tables and figures:** [Quarto](https://quarto.org), `tidyverse`, `kableExtra`, `cowplot`, `patchwork`, `ggh4x`, `RColorBrewer`, `scales`, `mgcv`, `pROC`

The full `sessionInfo()` of the R installation used on the HPC is in [`sessionInfo.txt`](sessionInfo.txt). It was recorded after the simulation finished, from the same R installation and package library.

---

## Repository structure

```
.
├── simulation_code/
│   ├── 0_set_up/
│   │   ├── helpers.R                  # data generation, original model fits, AUC/Brier
│   │   ├── set_up_linear.Rmd          # scenarios 1–2: seeds, initial datasets, original models
│   │   ├── set_up_non_linear.Rmd      # scenarios 3–4: builds the random forest DGM, then as above
│   │   ├── seeds.RDS                  # dataset seeds, scenarios 1–2
│   │   ├── seeds_non_linear.RDS       # dataset seeds, scenarios 3–4
│   │   └── rf_dgm.RDS                 # the non-linear (random forest) DGM used in the paper
│   └── 1_stability_simulation/
│       ├── execute.R                  # runs one method on one initial dataset
│       ├── simulation_helpers.R       # model fitting + all resampling methods (linear DGM)
│       ├── simulation_helpers_nl.R    # same, for the non-linear DGM
│       ├── rf_dgm.RDS                 # copy of the DGM, read by simulation_helpers_nl.R
│       └── shell_*.sh                 # one SLURM array job per method
└── results/
    └── scripts/
        ├── parquet_conversion.R            # raw .RDS output -> Parquet
        ├── calculate_evaluation_criteria.R # per-individual instability measures
        └── figures_and_tables.qmd          # tables and figures in the paper
```

Scenarios are numbered as follows throughout the code:

| Scenario | DGM | N |
|---|---|---|
| 1 | Linear | 458 |
| 2 | Linear | 229 |
| 3 | Non-linear | 458 |
| 4 | Non-linear | 229 |

---

## Reproducing the results

The simulation code was run step by step, with a few settings changed by hand between runs. The steps below say exactly what to change.

### Step 1: Generate the initial datasets and fit the original models

In `simulation_code/0_set_up/`, first create the output folders:

```bash
cd simulation_code/0_set_up
mkdir -p raw_data/scenario_{1..4} raw_data/external_{1,2}
mkdir -p results/models/scenario_{1..4} results/predictions/scenario_{1..4}
mkdir -p results/hyperparameters/{elastic_net,random_forest}/scenario_{1..4}
```

Then run the chunks of `set_up_linear.Rmd` (scenarios 1–2) and `set_up_non_linear.Rmd` (scenarios 3–4) in order. Each one:

1. creates the dataset seeds (`seeds.RDS`, `seeds_non_linear.RDS`; both are included)
2. generates 100 initial datasets per scenario in `raw_data/`
3. fits elastic net and random forest to each initial dataset, saving the models, original predictions and selected hyperparameters in `results/`

`set_up_non_linear.Rmd` also builds the non-linear DGM: it generates a preparatory dataset (n = 10,000), then fits the random forest that is saved as `rf_dgm.RDS`. This step is slow. The DGM used in the paper is included, so you can skip this chunk and use the saved file.

### Step 2: Run the stability simulation

In `simulation_code/1_stability_simulation/`, make the initial datasets available, and create the output and log folders:

```bash
cd simulation_code/1_stability_simulation
cp -r ../0_set_up/raw_data .
mkdir -p results/scenario_{1..4} logs
```

`execute.R` runs **one method on one initial dataset of one scenario**. Before each run, edit two lines in `execute.R`:

- **Scenario:** set `sc` to 1, 2, 3 or 4.
- **Helper file:** use `source("simulation_helpers.R")` for scenarios 1–2, and `source("simulation_helpers_nl.R")` for scenarios 3–4.

In the `shell_*.sh` files, replace the path to `Rscript` and the `--mail-user` address with your own.

Each SLURM script submits an array of 100 tasks, one per initial dataset:

| SLURM script | Method | Output file (in `results/scenario_<sc>/`) | Time limit |
|---|---|---|---|
| `shell_independent.sh` | Independent sampling (reference) | `scenario_<sc>_independent_<i>.RDS` | 25 h |
| `shell_bootstrap.sh` | Non-parametric bootstrap | `scenario_<sc>_bootstrapping_<i>.RDS` | 25 h |
| `shell_bootstrap_par.sh` | Parametric bootstrap | `scenario_<sc>_bootstrapPAR_<i>.RDS` | 25 h |
| `shell_bootstrap_635.sh` | Out-of-bag bootstrap | `scenario_<sc>_bootstrap635_<i>.RDS` | 25 h |
| `shell_cross_validation_10.sh` | Repeated 10-fold CV | `scenario_<sc>_crossvalidation_<i>.RDS` | 60 h |
| `shell_cross_validation_5.sh` | Repeated 5-fold CV | `scenario_<sc>_crossvalidation5_<i>.RDS` | 60 h |
| `shell_cross_validation_2.sh` | Repeated 2-fold CV | `scenario_<sc>_crossvalidation2_<i>.RDS` | 60 h |

For example, to run the non-parametric bootstrap for all 100 initial datasets:

```bash
sbatch shell_bootstrap.sh
```

To run a single initial dataset without SLURM (here dataset 1 with the bootstrap):

```bash
Rscript execute.R 1 bootstrap
```

Run times differ a lot between methods. Relative to the non-parametric bootstrap, repeated 2-fold CV takes about 2×, repeated 5-fold CV and the parametric bootstrap about 3×, and repeated 10-fold CV and the out-of-bag bootstrap about 4–6×.

### Step 3: Convert results to Parquet and calculate the instability measures

Run the scripts in `results/scripts/` from the repository root (they use `here::here()`):

1. **`parquet_conversion.R`** collects the `.RDS` output from Steps 1 and 2 and writes it as Parquet to `results/parquet_v2/`, `results/original_predictions_v2/` and `results/hyperparameters_v2/`. Before running, point the input paths (lines 34–57) to your output:
   - `resampling_source_root` → `simulation_code/1_stability_simulation/results`
   - `original_source_root` → `simulation_code/0_set_up/results/predictions`
   - `hyperparameter_source_root` → `simulation_code/0_set_up/results/hyperparameters`
2. **`calculate_evaluation_criteria.R`** calculates the per-individual measures (variance, 95% range width, MAPE, tMAPE, 2.5th and 97.5th percentiles) and writes them to `results/evaluation_v2/`.

### Step 4: Produce the tables and figures

`results/scripts/figures_and_tables.qmd` reads the Parquet files from Step 3 and produces the tables and figures in the paper. Render it from the repository root:

```bash
quarto render results/scripts/figures_and_tables.qmd
```

It currently covers:

- **Main text:** Table 3, Figures 1–3
- **Supplementary material:** S-Table 1, S-Figures 1–9, 12–21 and 26–29

The code for Figures 4–7 and S-Figures 10–11, 22–25 and 30–31 will be added in a later update.

---

## Citation

```bibtex
@article{carriero2026stability,
  title   = {Estimating Stability of Individual-Specific Predictions from {AI}-based
             prediction models using Bootstrap and Cross-validation Resampling},
  author  = {Carriero, Alex and van Dam, Martijn L. G. and de Hond, Anne A. H. and
             Riley, Richard D. and Collins, Gary S. and Moons, Karel G. M. and
             van Smeden, Maarten and Erler, Nicole S.},
  journal = {arXiv preprint arXiv:2610.33744},
  year    = {2026},
  url     = {https://arxiv.org/abs/2610.33744}
}
```

## Licence

The code in this repository is released under the [MIT License](LICENSE).

## Contact

Alex Carriero, Julius Center for Health Sciences and Primary Care, UMC Utrecht: a.j.carriero@umcutrecht.nl

For questions about the code, please open a [GitHub issue](../../issues).
