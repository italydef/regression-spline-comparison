# Regression spline fitting of piecewise-smooth functions — a comparative study

*A reproducible comparison of 15 regression-spline methods for fitting
piecewise-smooth functions with jumps and cusps, on simulated and real data.*

Code and data to reproduce the simulation study and real-data analyses in:

> U. Amato, A. Antoniadis, I. De Feis, and I. Gijbels.
> *On algorithms for regression spline fitting of piecewise-smooth functions: a comparative study.*
> Manuscript, 2026 (unpublished).

The study compares **15 regression-spline fitting procedures** — EAS, RC, HAS, TF,
FAPS, EAPS, SARS, ASPL, FREEK, BAPS, DMS, BARS, CPR, AMP, SCHA — on six simulated
test functions and on two real datasets (the USD/ISK exchange rate and a well-log
series). There is no build system: "running" means sourcing the R scripts below in
order.

## Requirements

Tested with **R 4.6.0 (2026-04-24)** on macOS (Intel iMac, Ventura 13.7; Apple M2
MacBook Pro, Sequoia 15.7). The code also runs on Intel/Linux; only the native
components (see below) are architecture-specific.

Two methods require external software:

- **MATLAB R2024b** — used only by **RC** and **BAPS** (through the `R.matlab`
  bridge). Set `MATLAB_BIN` at the top of each `run_*.R` to your MATLAB binary; this
  is the only machine-specific path in the code.
- **JAGS** — required by **DMS** (`miscF` → `R2jags` → `rjags`).
  Download: https://sourceforge.net/projects/mcmc-jags/files/

### CRAN packages

Install from CRAN (versions used in the study):

```r
install.packages(c(
  "R.matlab",     # 3.7.0   (RC, BAPS)
  "wavethresh",   # 4.7.3   (noise-scale estimate in the SIMULATION)
  "genlasso",     # 1.6.1   (TF)
  "mgcv",         # 1.9-4   (EAPS)
  "cpr",          # 0.4.1   (CPR)
  "aspline",      # 0.2.0   (ASPL)
  "splines2",     # 0.5.4   (ASPL; must be attached, see note below)
  "assist",       # 3.1.9   (HAS)
  "fields",       # 17.3    (EAS)
  "pso",          # 1.0.4   (EAPS)
  "miscF",        # 0.1-5   (DMS; pulls in rjags -> needs JAGS installed)
  "GeDS"          # 0.3.5   (REAL DATA only: noise-scale for RC/SARS/AMP)
))
```

Note the simulation/real-data split: the **simulation** uses `wavethresh` for the
RC/SARS/AMP noise-scale estimate and does **not** need `GeDS`; the **real-data**
scripts use `GeDS::NGeDS` and do **not** need `wavethresh`. `splines2` must be
*attached* (`library(splines2)`), not merely installed, because ASPL calls
`bSpline()` un-namespaced.

### Bundled packages (in `CODES_TO_BE_INSTALLED/`)

These are provided in this repository and installed from the local tarballs/sources.
See **[`CODES_TO_BE_INSTALLED/PACKAGE_NOTES.md`](CODES_TO_BE_INSTALLED/PACKAGE_NOTES.md)** for full details and
original download links, and **[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)**
for provenance and licenses.

```r
# order matters: SemiPar before AdaptFitOS
install.packages("CODES_TO_BE_INSTALLED/SemiPar_1.0-4.2.tar.gz",        repos = NULL, type = "source")
install.packages("CODES_TO_BE_INSTALLED/AdaptFitOS_0.69.tar.gz",        repos = NULL, type = "source")  # FAPS
install.packages("CODES_TO_BE_INSTALLED/bsml_1.5-1.tar.gz",             repos = NULL, type = "source")  # HAS
install.packages("CODES_TO_BE_INSTALLED/freeknotsplines_1.0.1.tar.gz",  repos = NULL, type = "source")  # FREEK
install.packages("CODES_TO_BE_INSTALLED/sars_1.0.tar.gz",               repos = NULL, type = "source")  # SARS
# SCHACE (SCHA): unzip CODES_TO_BE_INSTALLED/SCHACE-main.zip and install per its README
# BARS: see CODES_TO_BE_INSTALLED/INSTALL_bars.pdf (builds an executable + shared library)
```

### Native code (rebuild per machine — do not reuse binaries)

Compiled artefacts (`*.so`, `*.o`, `*.dll`) are tied to one OS **and** one CPU and are
deliberately **not** committed. Rebuild them on your machine:

- **smoothness** (Fortran, used by **EAS**): the `smoothness-master/` folder
  contains the `.f90` source and `cp_source.R`. Build the shared library in place:
  ```
  cd smoothness-master
  R CMD SHLIB -O3 localrisk.f90 localfit.f90 -o localmethod64.so
  ```
  Full cross-platform instructions (including the gfortran setup) are in
  **[`USING_smoothness-master.md`](USING_smoothness-master.md)**.
  All three driver scripts reference this one folder
  (`SMOOTH_DIR <- file.path(WORKING_DIR, "smoothness-master")`), so the build
  works the same on any architecture.
- **sars** and **BARS** also contain compiled code; `R CMD INSTALL` (sars) and the
  BARS build notes (`CODES_TO_BE_INSTALLED/INSTALL_bars.pdf`) recompile them for your machine.

## Repository layout

```
.
├── METHODS.R               # method wrappers for the SIMULATION
├── METHODS_REAL_DATA.R     # method wrappers for the REAL DATA
├── run_simulation.R        # driver: simulation
├── run_realdata_usdisk.R   # driver: USD/ISK exchange-rate data
├── run_realdata_well.R     # driver: well-log data
├── make_boxplots.R         # figures from simulation results (RMSE / MXDV boxplots)
├── make_heatmap.R          # figures: ranking heatmap
├── make_estimates.R        # figures: example curve-estimate panels
├── functions.R             # true/noisy test-function figure
├── default_knots.r         # default.knots() helper (FAPS, CPR)
├── varcalib.R              # noise-scale calibration (RC, SARS, AMP)
├── AMPScode/               # AMP method (bundled)
├── OKPSPS-main/            # EAPS method (bundled, from AnFreTh/OKPSPS)
├── R_software/             # sourced helpers (bundled)
├── adaptive_psplinesv06/   # MATLAB code for RC / BAPS (bundled)
├── smoothness-master/      # EAS Fortran source (sourced, not installed; build the .so in place)
├── USING_smoothness-master.md  # how to build the smoothness .so on any OS
├── INPUT/
│   ├── Simulation/         # simulation results land here (generated)
│   └── RealData/           # usdisk.csv, well.txt, well.RData
├── OUTPUT/
│   ├── Simulation/         # simulation figures (generated)
│   └── RealData/           # real-data figures (generated)
├── CODES_TO_BE_INSTALLED/  # third-party packages to INSTALL (tarballs, BARS, SCHACE) + install docs
├── THIRD_PARTY_NOTICES.md  # provenance + licenses of bundled code
└── LICENSE                 # MIT (authors' own code only)
```

## Reproducing the results

Run R **from the repository root** — all paths are relative to the working directory.

### Simulation

```r
source("run_simulation.R")
```

Configurable at the top of `run_simulation.R`:

```r
n_functs     <- 6     # number of test functions
n_methods    <- 15    # number of methods (must equal length(METHODS))
n_sim        <- 50    # replications per function
n_points     <- 256   # design points
bdeg         <- 2     # B-spline degree (used inside some method wrappers)
SNR          <- 3     # signal-to-noise ratio, var(f) / sigma^2
estimate_rep <- 1     # replication kept for the estimate panels (Figs 3.4/3.5)
```

This writes the results to `INPUT/Simulation/`. Then produce the figures:

```r
source("make_boxplots.R")
source("make_heatmap.R")
source("make_estimates.R")
```

Figures are written to `OUTPUT/Simulation/`.

To plot true functions f1, f2, f3, f4, f4, f6 plus a noisy sumulation, run functions.R.
The output plots are written to `OUTPUT/Simulation/`.


### Real data

Independent of the simulation:

```r
source("run_realdata_well.R")     # well-log data
source("run_realdata_usdisk.R")   # USD/ISK exchange-rate data
```

Datasets are read from `INPUT/RealData/`; figures and fits are written to
`OUTPUT/RealData/`. (BAPS is unstable on the well-log data and is intentionally
excluded from that figure.)

## Notes and gotchas

- Each method is wrapped in `tryCatch`; a method that fails leaves `NA`s rather than
  stopping the run.
- `set.seed(123)` is set for reproducibility — keep it.
- `run_simulation.R` sources `METHODS.R` (not `METHODS_REAL_DATA.R`); the two method
  files differ in their noise-scale estimator and must not be merged.

## License and attribution

The authors' own code is released under the **MIT License** (`LICENSE`). Bundled
third-party components remain under their original licenses; see
`THIRD_PARTY_NOTICES.md`. If you use this code, please cite the paper above.
