# Package notes

These notes list the packages required by the study and where each comes from.
(This is a Markdown rendering of `PACKAGE_NOTES.docx`, kept alongside it.)

Tested on an iMac 4.2 GHz Intel Core i7 quad-core (Ventura 13.7.8) and a MacBook Pro
Apple M2 (Sequoia 15.7.7), under **R 4.6.0 (2026-04-24)**.

## Install from CRAN

| Package | Version | Used by |
|---|---|---|
| `R.matlab` | 3.7.0 | RC, BAPS |
| `wavethresh` | 4.7.3 | noise scale (simulation) |
| `genlasso` | 1.6.1 | TF |
| `mgcv` | 1.9-4 | EAPS |
| `cpr` | 0.4.1 | CPR |
| `aspline` | 0.2.0 | ASPL |
| `splines2` | 0.5.4 | ASPL (attach with `library(splines2)`) |
| `assist` | 3.1.9 | HAS |
| `fields` | 17.3 | EAS |
| `pso` | 1.0.4 | EAPS |
| `miscF` | 0.1-5 | DMS — depends on `rjags` 4-17, which depends on JAGS |
| `GeDS` | 0.3.5 | real datasets only |

JAGS: https://sourceforge.net/projects/mcmc-jags/files/

## Provided with the code (in this `install/` folder)

- **SCHACE** — https://github.com/ZhaoyingLuLuLu/SCHACE (`SCHACE-main.zip`)
- **AdaptFitOS** — `AdaptFitOS_0.69.tar.gz`, from
  https://cran.r-project.org/src/contrib/Archive/AdaptFitOS/ (requires SemiPar)
- **SemiPar** — `SemiPar_1.0-4.2.tar.gz`, from
  https://cran.r-project.org/src/contrib/Archive/SemiPar/
- **bsml** — `bsml_1.5-1.tar.gz`, from
  https://cran.r-project.org/src/contrib/Archive/bsml/
- **freeknotsplines** — `freeknotsplines_1.0.1.tar.gz`, from
  https://cran.r-project.org/src/contrib/Archive/freeknotsplines/
- **sars** — `sars_1.0.tar.gz`, obtained from the authors
- **BARS** — see https://www.jstatsoft.org/article/view/v026i01 and
  https://www.stat.cmu.edu/~kass/bars/bars.html; for installation on any OS see
  `INSTALL_bars.pdf` (LaTeX source `INSTALL_BARS.tex` provided too)

## Also used, provided with the code (at the repository root)

- **OKPSPS** — https://github.com/AnFreTh/OKPSPS
- **adaptive_psplinesv06** — obtained from the authors (MATLAB)
- **AMPScode** — obtained from the authors
- **R_software**
- **smoothness-master** — obtained from the authors; for installation on any OS see
  `USING_smoothness-master.md`

## Reproducing the study

**Simulation.** Run `run_simulation.R` (settings at the top: `n_functs`,
`n_methods`, `n_sim`, `n_points`, `bdeg`, `SNR`, `estimate_rep`). Results go to
`INPUT/Simulation`. Plots: run `make_boxplots.R`, `make_estimates.R`,
`make_heatmap.R` → `OUTPUT/Simulation`.

**Real data.** Run `run_realdata_well.R` (well-log) and `run_realdata_usdisk.R`
(USD/ISK). Datasets are in `INPUT/RealData`; outputs in `OUTPUT/RealData`.
