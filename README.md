
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Repository overview

This repository contains the code accompanying the paper **“Collective
Outlier Detection and Enumeration with Conformalized Closed Testing.”**

It provides implementations of ACODE (Automatic Conformal Outlier
Detection and Enumeration), closed-testing procedures based on rank
tests, and all simulation and real-data experiments reported in the
paper, including those in the appendices. The repository includes:

- An R package (`nout/`) implementing ACODE and the closed-testing
  procedures.
- Python and R experiment scripts (`experiments/`) used to run the
  simulation studies and real-data analyses, together with shell
  launchers that drive parameter sweeps on either a Slurm cluster or a
  single machine.
- Core Python utilities (`methods/`) for data generation, conformal
  scoring, ACODE orchestration, and interfacing with the R package for
  closed testing via `rpy2`.
- R scripts (`experiments/figures_and_tables/`) that assemble the final
  figures and tables in the paper from the experiment outputs.
- A metadata layer (`experiments/metadata/`) mapping every figure and
  table in the manuscript to the scripts that produce it (see the
  *Reproducibility map* section below).
- Precomputed critical value tables (`tables/`) for commonly used local
  tests (Fisher, Wilcoxon–Mann–Whitney, and Shiraishi tests against
  Lehmann alternatives), used to accelerate small-sample calculations
  where asymptotic approximations may be unreliable.
- A `data/` folder containing the datasets used in the experiments (see
  the *Data* section below).

> **Platform.** This code was developed and tested on Linux. macOS
> should work for the local (non-cluster) workflow, though shell scripts
> assume GNU coreutils and bash 4+ (install via
> `brew install bash coreutils`); the cluster workflow requires Slurm.
> Windows is not supported.

The sections below describe how the repository is organized, how to
install the dependencies, and how to reproduce every figure and table
reported in the paper (including those in the appendices).

# Data

> **Note.** The `data/` folder is not hosted on GitHub due to size and
> possible licensing restrictions. All processed data files are included
> in the supplementary material. Original sources are listed below; see
> Appendix A6.3 of the manuscript for further details.

| Dataset                      | File(s)                                                                                                                         | Description                                                       | Source / Prior usage                                                                                             |
|------------------------------|---------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------|
| **LHC Olympics 2020 (LHCO)** | `events_anomalydetection_v2.features.h5`, `events_anomalydetection_Z_XY_qqq.features.h5`, `events_anomalydetection_Z_XY_qqq.h5` | Community benchmark for anomaly detection in high-energy physics. | Kasieczka et al., 2021.                                                                                          |
| **ALOI**                     | `aloi.arff`                                                                                                                     | 27 variables, 1,508 outliers, 48,026 inliers.                     | Amsterdam Library of Object Images (Geusebroek et al., 2005). Previously used in Bates et al., 2023.             |
| **Covertype**                | `cover.mat`                                                                                                                     | 10 variables, 2,747 outliers, 286,048 inliers.                    | Covertype dataset from UCI / ODDS repository. Previously used in Bates et al., 2023.                             |
| **CreditCard**               | `creditcard.csv`                                                                                                                | 30 variables, 492 outliers, 284,315 inliers.                      | Kaggle Credit Card Fraud Dataset. Used in Bates et al., 2023; and Marandon et al., 2024.                         |
| **Mammography**              | `mammography.mat`                                                                                                               | 6 variables, 260 outliers, 10,923 inliers.                        | Mammography dataset from ODDS repository. Used in Bates et al., 2023; Liang et al., 2024; Marandon et al., 2024. |
| **Pendigits**                | `pendigits.mat`                                                                                                                 | 16 variables, 156 outliers, 6,714 inliers.                        | Pen-Based Recognition of Handwritten Digits dataset (ODDS repository). Used in Bates et al., 2023.               |
| **Shuttle**                  | `shuttle.mat`                                                                                                                   | 9 variables, 3,511 outliers, 45,586 inliers.                      | Statlog Shuttle dataset (ODDS repository). Used in Bates et al., 2023; and Marandon et al., 2024.                |

# Code organization

## R package `nout/` (required)

All closed-testing procedures rely on the R package **`nout`**, included
in this repository under `nout/`. This package is called from both the R
experiments and the Python experiments (via `rpy2`), so installing
`nout` is required for all workflows.

Before running experiments, install the package locally:

``` r
install.packages("devtools")  # if needed
devtools::install("nout/")
```

## Methods (`methods/`)

Core utilities used by the Python experiments: data handling, synthetic
data generators, conformal score construction, the ACODE implementation,
and the interface to the R package for closed testing.

- `methods/models.py` — synthetic data generators used in the simulation
  studies.
- `methods/my_utils.py` — utilities for experiment bookkeeping and
  evaluation.
- `methods/conformal.py` — ACODE implementation. Imports the R package
  `nout` via `rpy2` (`importr("nout")`) and exposes the closed-testing
  routines for use inside the Python code.

## Experiments (`experiments/`)

This is where the simulation and real-data experiments are launched and
where their outputs land. The folder layout is:

    experiments/
    ├── exp_A.py / exp_A.sh                          # Python experiment driver (synthetic, real tabular, LHCO)
    ├── exp_B1_global_testing.R / .sh                # R experiment: global testing (Appendix B1)
    ├── exp_B2_enumeration.R    / .sh                # R experiment: enumeration (Appendix B2)
    ├── exp_C1_time.R                                # R experiment: timing comparison (Appendix C1)
    ├── exp_C2_simes.R                               # R experiment: Simes-permutation comparison (Appendix C2)
    ├── make_figures_tables.sh                       # builds all figures and tables sequentially
    ├── submit_exp_A.sh                              # Grid launcher for exp_A (for Figures 1-4, Figures A2-A10, Figures A13-A22, Table A4)
    ├── submit_exp_B1.sh                             # Grid launcher for exp_B1 (for Figure A11)
    ├── submit_exp_B2.sh                             # Grid launcher for exp_B2 (for Figure A12)
    ├── submit_exp_C1.sh                             # Wrapper for exp_C1 (for Figure A1)
    ├── submit_exp_C2.sh                             # Wrapper for exp_C2 (for Table A2)
    ├── run_all_experiments.sh                       # Meta-launcher: submits every experiment via map_submit.tsv
    ├── figures_and_tables/                          # R scripts that turn experiment outputs into final figures/tables
    │   ├── make_fig*.R, make_tab*.R                 # one make_*.R per paper figure/table group
    │   ├── utils_plotting_*.R                       # helper functions for making figures and tables
    │   ├── figures/                                 # produced figures (.pdf)
    │   ├── tables/                                  # produced tables (.tex)
    │   └── logs/                                    # per-script R output logs
    ├── metadata/                                    # provenance map: paper item ↔ scripts ↔ result directory
    │   ├── map_full.tsv, map_full.md, map_full.tex  # consolidated map (see Reproducibility map section)
    │   ├── map_submit.tsv, map_results.tsv          # extracted from submit_*.sh headers
    │   ├── map_scripts.tsv, map_figures.tsv         # extracted from make_*.R headers
    │   ├── extract_submit_metadata.sh               # rebuilds the submit-side maps
    │   ├── extract_make_metadata.sh                 # rebuilds the make-side maps
    │   └── cross_reference_metadata.sh              # joins everything into map_full.tsv / .md
    ├── logs/                                        # cluster job logs (one subdir per result key)
    └── results/                                     # experiment outputs (one subdir per result key)

All submit launchers accept the same flags for consistency:

- `--cluster` (default) submits jobs via `sbatch`.
- `--local` runs jobs sequentially on the current machine.
- `--dry-run` prints the commands without executing or submitting them.
- `--force` re-runs experiments even if result files already exist.

They also skip any configuration whose output already exists, so
launching the same script a second time only fills in the missing
pieces.

# Reproducibility map (paper item → scripts)

| Paper item | Submit script    | Results directory | Make script   | Result key |
|------------|------------------|-------------------|---------------|------------|
| Figure 1   | submit_exp_A.sh  | results/fig1/     | make_fig1.R   | fig1       |
| Figure 2   | submit_exp_A.sh  | results/fig2/     | make_fig2.R   | fig2       |
| Figure 3   | submit_exp_A.sh  | results/fig3/     | make_fig3.R   | fig3       |
| Figure 4   | submit_exp_A.sh  | results/fig4/     | make_fig4.R   | fig4       |
| Figure A1  | submit_exp_C1.sh | results/figA1/    | make_figA1.R  | figA1      |
| Figure A2  | submit_exp_A.sh  | results/figA2/    | make_figA2.R  | figA2      |
| Figure A3  | submit_exp_A.sh  | results/fig3/     | make_fig3.R   | fig3       |
| Figure A4  | submit_exp_A.sh  | results/fig2/     | make_fig2.R   | fig2       |
| Figure A5  | submit_exp_A.sh  | results/figA5/    | make_fig2.R   | figA5      |
| Figure A6  | submit_exp_A.sh  | results/figA6/    | make_fig3.R   | figA6      |
| Figure A7  | submit_exp_A.sh  | results/fig2/     | make_fig2.R   | fig2       |
| Figure A8  | submit_exp_A.sh  | results/figA6/    | make_fig3.R   | figA6      |
| Figure A9  | submit_exp_A.sh  | results/figA5/    | make_fig2.R   | figA5      |
| Figure A10 | submit_exp_A.sh  | results/figA10/   | make_figA10.R | figA10     |
| Figure A11 | submit_exp_B1.sh | results/figA11/   | make_figA11.R | figA11     |
| Figure A12 | submit_exp_B2.sh | results/figA12/   | make_figA12.R | figA12     |
| Figure A13 | submit_exp_A.sh  | results/fig1/     | make_fig1.R   | fig1       |
| Figure A14 | submit_exp_A.sh  | results/figA14/   | make_figA14.R | figA14     |
| Figure A15 | submit_exp_A.sh  | results/figA15/   | make_figA15.R | figA15     |
| Figure A16 | submit_exp_A.sh  | results/figA15/   | make_figA15.R | figA15     |
| Figure A17 | submit_exp_A.sh  | results/figA15/   | make_figA15.R | figA15     |
| Figure A18 | submit_exp_A.sh  | results/figA15/   | make_figA15.R | figA15     |
| Figure A19 | submit_exp_A.sh  | results/figA15/   | make_figA15.R | figA15     |
| Figure A20 | submit_exp_A.sh  | results/figA15/   | make_figA15.R | figA15     |
| Figure A21 | submit_exp_A.sh  | results/figA21/   | make_figA21.R | figA21     |
| Figure A22 | submit_exp_A.sh  | results/figA21/   | make_figA21.R | figA21     |
| Table A2   | submit_exp_C2.sh | results/tabA2/    | make_tabA2.R  | tabA2      |
| Table A4   | submit_exp_A.sh  | results/fig1/     | make_fig1.R   | fig1       |

``` bash
cd experiments/metadata
./extract_submit_metadata.sh
./extract_make_metadata.sh
./cross_reference_metadata.sh
```

This is also the table to consult when only a subset of figures needs to
be regenerated.

# Reproducibility workflow

## Dependencies

#### Version of primary software used

- R version: **4.4.2**
- Python version: **3.12.3**

#### Libraries and dependencies used by the code

The R package `nout/` (installed locally from this repository) depends
on the following R packages (from its `DESCRIPTION` `Imports:` field):

- `foreach` (1.5.2)
- `hommel` (1.8)
- `multcomp` (1.4-30)
- `stats` (base R)
- `sumSome` (1.1.1)
- `Iso` (0.0-21)
- `mixmodel` (\>= 0.5; installed version: 0.5; not on CRAN — installed
  from GitHub via `devtools::install_github("rohitpatra/mixmodel")`)
- `fitdistrplus` (1.2-6)
- `GoFKernel` (2.1-3)

In addition, the repository uses the following R packages for experiment
orchestration, result aggregation, and figure/table production:

- `tidyverse` (2.0.0), including `ggplot2` (3.5.2), `dplyr` (1.1.4),
  `tidyr` (1.3.1), `readr` (2.1.5), and `tibble`
- `progress` (1.2.3)
- `ggh4x` (0.3.1)
- `scales` (1.4.0)
- `latex2exp` (0.9.6)
- `RColorBrewer` (1.1.3)
- `kableExtra` (1.4.0)
- `combinat` (0.0.8)
- `knitr` (1.50)
- `devtools` (2.4.5) — used to install the local `nout/` package
- `testthat` (3.2.1.1) — only required to run unit tests for `nout/`

The Python experiments and ACODE implementation rely on standard
scientific computing and machine learning libraries:

- `numpy` (1.26.4)
- `pandas` (2.2.3)
- `scipy` (1.15.2)
- `scikit-learn` (1.6.1)
- `statsmodels` (0.14.4)
- `tqdm` (4.67.1)
- `matplotlib` (3.10.1)
- `seaborn` (0.13.2)
- `joblib` (1.4.2)
- `mat73` (0.65) — required for reading HDF5-based `.mat` files; depends
  on the system HDF5 library (`libhdf5` on Linux, available via
  `brew install hdf5` on macOS)

The Python code interfaces with the R implementation of closed testing
via:

- `rpy2` (3.5.11)

## Workflow

#### Precomputed results included

The supplementary materials include the precomputed experiment outputs
under `experiments/results/`, organized by result key:

    experiments/results/
    ├── fig1/    fig2/    fig3/    fig4/
    ├── figA1/   figA2/   figA5/   figA6/   figA10/
    ├── figA11/  figA12/  figA14/  figA15/  figA21/
    └── tabA2/

This means every figure and table in the manuscript can be reproduced
without re-running the experiments — skip Step 1 below and proceed
directly to Step 2. The experiment-execution stage (Step 1) is
documented for completeness and for anyone wishing to regenerate results
from scratch.

**Note.** The `experiments/results/` folder is not hosted on GitHub due
to size.

#### Step 0: Install dependencies

To install all dependencies from scratch:

``` r
# R --- core install
install.packages("devtools")

devtools::install_github("rohitpatra/mixmodel")   # not on CRAN

devtools::install("nout/")                        # local package

install.packages(c(
  "tidyverse", "progress", "ggplot2", "ggh4x", "scales", "latex2exp",
  "RColorBrewer", "kableExtra", "combinat", "knitr"
))
```

``` bash
# Python
pip install numpy pandas scipy scikit-learn statsmodels tqdm rpy2 \
            matplotlib seaborn joblib mat73
```

#### Step 1: Run experiments

Reproducing every figure and table involves two stages: (1) running the
experiments to populate `experiments/results/`, and (2) building the
figures and tables from those results.

To reproduce a **specific** figure or table, look up its row in the
“reproducibility map” table (also saved in
`experiments/metadata/map_full.md`), then run the corresponding submit
script with the result key as the argument.

For example, to regenerate the results behind Figure 2 (which requires
`results/fig2/`):

``` bash
cd experiments
./submit_exp_A.sh fig2 --cluster        # submit ~thousands of jobs via sbatch
# or, for local sequential execution (slow):
./submit_exp_A.sh fig2 --local
# or, to just see what would happen:
./submit_exp_A.sh fig2 --local --dry-run
```

To reproduce **all** experiments in the paper:

``` bash
cd experiments
./run_all_experiments.sh                # cluster mode, asks for confirmation
./run_all_experiments.sh --dry-run -y   # preview the full job list, no submission
```

`run_all_experiments.sh` reads `metadata/map_submit.tsv` and invokes
each submit script with every result key it produces. The script prints
a warning before submitting because the total number of jobs is in the
thousands; this stage is intended for a cluster (see *Computing time*
below).

#### Step 2: Build figures and tables

Once the experiments have completed and `experiments/results/` is
populated, build every figure and table in the paper with:

``` bash
cd experiments
./make_figures_tables.sh
```

This reads `metadata/map_scripts.tsv` and runs each `make_*.R` in turn,
redirecting the per-script R output to `figures_and_tables/logs/`.
Output PDFs land in `figures_and_tables/figures/` and `.tex` tables in
`figures_and_tables/tables/`. A summary at the end reports how many
scripts succeeded and where the outputs are.

To rebuild a **single** figure or table instead, run the corresponding
`make_*.R` directly:

``` bash
cd experiments/figures_and_tables
Rscript make_fig2.R                     # produces Figure 2, Figure A4, Figure A7
Rscript make_tabA2.R                    # produces Table A2
```

## Computing time

Each individual experiment runs in seconds to a few minutes on a
standard laptop. The aggregate computational cost arises from the large
number of configurations and repetitions, not from any single run. On a
Slurm cluster with dozens of parallel cores, the full pipeline completes
in a few hours; on a single laptop, completing it sequentially is
impractical.
