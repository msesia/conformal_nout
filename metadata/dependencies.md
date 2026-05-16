#### Version of primary software used

- R version: **4.4.2**
- Python version: **3.12.3**

#### Libraries and dependencies used by the code

The R package `nout/` (installed locally from this repository) depends on the following R packages (from its `DESCRIPTION` `Imports:` field):

- `foreach` (1.5.2)
- `hommel` (1.8)
- `multcomp` (1.4-30)
- `stats` (base R)
- `sumSome` (1.1.1)
- `Iso` (0.0-21)
- `mixmodel` (>= 0.5; installed version: 0.5; not on CRAN — installed from GitHub via `devtools::install_github("rohitpatra/mixmodel")`)
- `fitdistrplus` (1.2-6)
- `GoFKernel` (2.1-3)

In addition, the repository uses the following R packages for experiment orchestration, result aggregation, and figure/table production:

- `tidyverse` (2.0.0), including `ggplot2` (3.5.2), `dplyr` (1.1.4), `tidyr` (1.3.1), `readr` (2.1.5), and `tibble`
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

The Python experiments and ACODE implementation rely on standard scientific computing and machine learning libraries:

- `numpy` (1.26.4)
- `pandas` (2.2.3)
- `scipy` (1.15.2)
- `scikit-learn` (1.6.1)
- `statsmodels` (0.14.4)
- `tqdm` (4.67.1)
- `matplotlib` (3.10.1)
- `seaborn` (0.13.2)
- `joblib` (1.4.2)
- `mat73` (0.65) — required for reading HDF5-based `.mat` files; depends on the system HDF5 library (`libhdf5` on Linux, available via `brew install hdf5` on macOS)

The Python code interfaces with the R implementation of closed testing via:

- `rpy2` (3.5.11)
