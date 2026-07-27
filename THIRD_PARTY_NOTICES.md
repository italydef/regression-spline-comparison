# Third-party notices

This repository bundles third-party code for reproducibility. The MIT license in
`LICENSE` covers **only** the code written by the authors of the accompanying paper
(the `run_*.R`, `make_*.R`, `METHODS*.R`, `functions.R`, `default_knots.r`, and
`varcalib.R` scripts). Everything listed below is the work of others and remains
under its original license and copyright. It is included here solely to make the
simulation and real-data study reproducible; all credit belongs to the original
authors.

If you reuse any of these components outside this study, cite and comply with the
original source rather than this repository.

## Installable packages (`CODES_TO_BE_INSTALLED_OR_USED/`)

| Component | Source | License |
|---|---|---|
| `AdaptFitOS_0.69.tar.gz` | CRAN archive: https://cran.r-project.org/src/contrib/Archive/AdaptFitOS/ | GPL (≥ 2) |
| `SemiPar_1.0-4.2.tar.gz` (AdaptFitOS dependency) | CRAN archive: https://cran.r-project.org/src/contrib/Archive/SemiPar/ | GPL (≥ 2) |
| `bsml_1.5-1.tar.gz` | CRAN archive: https://cran.r-project.org/src/contrib/Archive/bsml/ | GPL (≥ 2) |
| `freeknotsplines_1.0.1.tar.gz` | CRAN archive: https://cran.r-project.org/src/contrib/Archive/freeknotsplines/ | GPL (≥ 2) |
| `sars_1.0.tar.gz` | Obtained from the authors (Zhou & Shen, SARS method) | See package; permission of original authors |
| `SCHACE-main.zip` | https://github.com/ZhaoyingLuLuLu/SCHACE | MIT |
| `BARS/` | https://www.stat.cmu.edu/~kass/bars/bars.html — JSS v26i01: https://www.jstatsoft.org/article/view/v026i01 | See BARS source / JSS article |
| `smoothness-master/`, `smoothness-master.zip` | Obtained from the authors | See package; permission of original authors |

## Bundled helper code (runtime, at repo root)

| Component | Source | Notes |
|---|---|---|
| `OKPSPS-main/` | https://github.com/AnFreTh/OKPSPS | Used for the EAPS method (`fit_gam_optim`). Check the upstream repo for its current license. |
| `smoothness-master/` | Obtained from the authors | Used for EAS (`pilotQV`). Fortran **source** only; the compiled `localmethod64.so`/`.o` are intentionally **not** shipped and must be built per machine in this folder (see `CODES_TO_BE_INSTALLED_OR_USED/USING_smoothness-master.md`). |
| `adaptive_psplinesv06/` | Obtained from the authors | MATLAB code used by RC and BAPS. |
| `AMPScode/` | Obtained from the authors (Fink & Wells, AMPS method) | Used for the AMP method. |
| `R_software/` | Obtained from the authors | Sourced helpers. |

## Method → original reference

Each compared method traces to a published method; see the paper's reference list.
Key attributions: EAS (Jang & Oh 2011), RC (Ruppert & Carroll 2000), HAS (Luo & Wahba
1997), TF (Tibshirani 2014), FAPS (Wiesenfarth et al. 2012), EAPS (Thielmann et al.
2025), SARS (Zhou & Shen 2001), ASPL (Goepp et al. 2025), FREEK (Miyata & Shen 2003),
BAPS (Baladandayuthapani et al. 2005), DMS (Denison et al. 1998), BARS (DiMatteo et al.
2001), CPR (Dewitt 2017), AMP (Fink & Wells 2004), SCHA (Lu et al. 2024).

> If any component here was shared with you privately and you would prefer it not be
> redistributed publicly, remove that folder/tarball from the repository and replace it
> with a download link in `CODES_TO_BE_INSTALLED_OR_USED/PACKAGE_NOTES.md`.
