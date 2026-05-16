#### Precomputed results included

The supplementary materials include the precomputed experiment outputs under `experiments/results/`, organized by result key:

```
experiments/results/
├── fig1/    fig2/    fig3/    fig4/
├── figA1/   figA2/   figA5/   figA6/   figA10/
├── figA11/  figA12/  figA14/  figA15/  figA21/
└── tabA2/
```

This means every figure and table in the manuscript can be reproduced without re-running the experiments --- skip Step 1 below and proceed directly to Step 2. The experiment-execution stage (Step 1) is documented for completeness and for anyone wishing to regenerate results from scratch.

**Note.** The `experiments/results/` folder is not hosted on GitHub due to size. 

#### Step 0: Install dependencies

To install all dependencies from scratch:

~~~r
# R --- core install
install.packages("devtools")

devtools::install_github("rohitpatra/mixmodel")   # not on CRAN

devtools::install("nout/")                        # local package

install.packages(c(
  "tidyverse", "progress", "ggplot2", "ggh4x", "scales", "latex2exp",
  "RColorBrewer", "kableExtra", "combinat", "knitr"
))
~~~

~~~bash
# Python
pip install numpy pandas scipy scikit-learn statsmodels tqdm rpy2 \
            matplotlib seaborn joblib mat73
~~~


#### Step 1: Run experiments

Reproducing every figure and table involves two stages: (1) running the experiments to populate `experiments/results/`, and (2) building the figures and tables from those results.

To reproduce a **specific** figure or table, look up its row in the "reproducibility map" table (also saved in `experiments/metadata/map_full.md`), then run the corresponding submit script with the result key as the argument. 

For example, to regenerate the results behind Figure 2 (which requires `results/fig2/`):

~~~bash
cd experiments
./submit_exp_A.sh fig2 --cluster        # submit ~thousands of jobs via sbatch
# or, for local sequential execution (slow):
./submit_exp_A.sh fig2 --local
# or, to just see what would happen:
./submit_exp_A.sh fig2 --local --dry-run
~~~

To reproduce **all** experiments in the paper:

~~~bash
cd experiments
./run_all_experiments.sh                # cluster mode, asks for confirmation
./run_all_experiments.sh --dry-run -y   # preview the full job list, no submission
~~~

`run_all_experiments.sh` reads `metadata/map_submit.tsv` and invokes each submit script with every result key it produces. The script prints a warning before submitting because the total number of jobs is in the thousands; this stage is intended for a cluster (see *Computing time* below).

#### Step 2: Build figures and tables

Once the experiments have completed and `experiments/results/` is populated, build every figure and table in the paper with:

~~~bash
cd experiments
./make_figures_tables.sh
~~~

This reads `figures_and_tables/script_map.tsv` and runs each `make_*.R` in turn, redirecting the per-script R output to `figures_and_tables/logs/`. Output PDFs land in `figures_and_tables/figures/` and `.tex` tables in `figures_and_tables/tables/`. A summary at the end reports how many scripts succeeded and where the outputs are.

To rebuild a **single** figure or table instead, run the corresponding `make_*.R` directly:

~~~bash
cd experiments/figures_and_tables
Rscript make_fig2.R                     # produces Figure 2, Figure A4, Figure A7
Rscript make_tabA2.R                    # produces Table A2
~~~

