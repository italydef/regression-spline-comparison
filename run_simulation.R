# =============================================================================
# run_simulation.R
# -----------------------------------------------------------------------------
# Monte Carlo simulation study for the regression-spline comparison paper
# (Amato, Antoniadis, de Feis, Gijbels).
#
# WHAT IT DOES
#   For each of 6 test functions f1..f6 and each of `n_sim` noisy replications,
#   applies all methods in METHODS and records two criteria (RMSE, MXDV).
#   For one chosen replication it also keeps every method's fitted curve, so the
#   per-method estimate panels (Figures 3.4 / 3.5) can be reproduced.
#
# OUTPUT (single file, consumed by the plotting scripts)
#   <RESULTS_DIR>/simulation_results.RData  containing:
#     rmse3d, mxdv3d : arrays [n_functs, n_methods, n_sim] of the two criteria
#     results_df     : tidy long data frame (test_function, sim_id, method,
#                      RMSE, MXDV)  <- used by make_boxplots.R and make_heatmap.R
#     estimates      : list, one data frame per function with columns
#                      x, noisy, true, <one column per method> for replication
#                      `estimate_rep`  <- used by make_estimates.R
#     methods, test_funcs, n_sim, n_points, SNR, estimate_rep : metadata
#
# The three companion scripts read this file; they do NOT re-run the simulation.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------
# The software is distributed as ONE self-contained folder. Put that folder
# wherever you like, then point R's working directory at it before running:
#   - RStudio: Session > Set Working Directory > To Source File Location, or
#   - console: setwd("/path/to/AMPS")
# Everything below is derived from that root, so no absolute paths need editing.
WORKING_DIR  <- getwd()            # root of the software folder
OKPSPS_DIR   <- file.path(WORKING_DIR, "OKPSPS-main")
SMOOTH_DIR <- file.path(WORKING_DIR, "smoothness-master")
R_software_DIR <- file.path(WORKING_DIR  , "R_software")
AMPSCODE_DIR <- file.path(WORKING_DIR, "AMPScode")
MATLAB_SRS   <- file.path(WORKING_DIR, "adaptive_psplinesv06", "srs")

# The ONLY genuinely machine-specific setting: path to YOUR MATLAB executable.
MATLAB_BIN <- "/Applications/MATLAB_R2024b.app/bin/matlab"

# Where the results file is written:
RESULTS_DIR <- file.path(WORKING_DIR, "INPUT", "Simulation")

# Package -> method map.  Tags:
#   METHOD(S): fn()   = called DIRECTLY by that wrapper in METHODS.R
#   (via sourced X)   = reached only through a sourced helper file
# NB: EAS calls no library directly - pilotQV() is sourced from cp_source.R 
#     (which uses fields::Tps); varcalib() is sourced.
library(R.matlab)         # RC, BAPS: MATLAB bridge - setVariable()/evaluate()/getVariable()
library(wavethresh)       # RC, SARS, AMP: wd()/accessD() for the noise (hatsigma) estimate
library(genlasso)         # TF: trendfilter() / cv.trendfilter()
library(mgcv)             # EAPS: gam()
library(cpr)              # CPR: cp() / bsplines()
library(aspline)          # ASPL: aspline()  (adaptive-ridge knot selection)
library(splines2)         # ASPL: bSpline() - MUST be attached; METHODS.R calls it un-namespaced
library(assist)           # HAS: cubic() reproducing-kernel basis
library(fields)           # EAS (via sourced cp_source.R): Tps() / fastTps()
library(pso)              # EAPS (via sourced OKPSPS/ok_gam.R): psoptim() in fit_gam_optim()
library(miscF)            # DMS: curve.polynomial.rjmcmc()  (pulls in R2jags -> rjags -> JAGS)
library(SCHACE)           # SCHA: main.SCHACE()
library(AdaptFitOS)       # FAPS: asp2()
library(bsml)             # HAS: bsml()  (also supplies stdz())
library(BARS)             # BARS: bars()
library(sars)             # SARS: sars()
library(freeknotsplines)  # FREEK: fit.search.numknots() / fitted.freekt()

library(dplyr)            # results assembly (dplyr::bind_rows)

# -----------------------------------------------------------------------------
# 2. Source the method implementations and helper bases
# -----------------------------------------------------------------------------
source(file.path(WORKING_DIR  , "default_knots.r"))
# Each name in METHODS must resolve to a function(yy, x).
source(file.path(WORKING_DIR  , "METHODS.R"))
# Variance / threshold calibration used by some of the methods.
source(file.path(WORKING_DIR  , "varcalib.R"))

source(file.path(R_software_DIR, "mybases.r"))
source(file.path(R_software_DIR, "mypnormalb.r"))

source(file.path(AMPSCODE_DIR, "AMPS.r"))
source(file.path(AMPSCODE_DIR, "VC.r"))
source(file.path(AMPSCODE_DIR, "FDRnorm.r"))
source(file.path(AMPSCODE_DIR, "quantileknots.r"))
source(file.path(AMPSCODE_DIR, "DifMatrix.r"))
source(file.path(AMPSCODE_DIR, "PLSfit.r"))
source(file.path(AMPSCODE_DIR, "plsmin.r"))

old_wd <- getwd()
setwd(OKPSPS_DIR)
source(file.path("data_simulation", "functions.R"))
source(file.path("module", "ok_gam.R"))
setwd(old_wd)

# smoothness-master: compiled Fortran shared object + its R interface
# (localmethod64.so is produced by the smoothness-master build step)
dyn.load(file.path(SMOOTH_DIR, "localmethod64.so"))
source(file.path(SMOOTH_DIR, "cp_source.R"))

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

# Two-valued sign: +1 for x > 0, -1 otherwise (note sgn(0) = -1, unlike base R).
sgn <- function(x) 2 * (x > 0) - 1

# Performance criteria for a single fit.
compute_criteria <- function(pred, true) {
  c(RMSE = sqrt(mean((pred - true)^2)),
    MXDV = max(abs(pred - true)))
}

# -----------------------------------------------------------------------------
# 4. Simulation settings
# -----------------------------------------------------------------------------
n_functs     <- 6     # number of test functions
n_methods    <- 15    # number of methods (must equal length(METHODS))
n_sim        <- 50    # replications per function
n_points     <- 256   # design points
bdeg         <- 2     # B-spline degree (used inside some method wrappers)
SNR          <- 3     # signal-to-noise ratio, var(f) / sigma^2
estimate_rep <- 1     # which replication to keep for the estimate panels (Fig 3.4/3.5)

# Method dispatch labels. These match both METHODS.R and the manuscript list.
METHODS <- c("EAS", "RC", "HAS", "TF", "FAPS", "EAPS", "SARS",
             "ASPL", "FREEK", "BAPS", "DMS", "BARS", "CPR", "AMP", "SCHA")

methods <- METHODS
stopifnot(length(methods) == n_methods)


# -----------------------------------------------------------------------------
# 5. Test functions f1..f6  
# -----------------------------------------------------------------------------
f1 <- function(x) ifelse(x < 0.5, -1.5, 0.25 * sin(x^2 * pi^1.5))
f2 <- function(x) ifelse(x < 0.6, 1 / (0.01  + (x - 0.30)^2),
                                  1 / (0.015 + (x - 0.65)^2))
f3 <- function(x) 100 / exp(abs(10 * x - 5)) + ((10 * x - 5)^5) / 500
f4 <- function(x) {
  y <- sin(15 * x) + 0.3 * x^2
  y[x <  0.4]           <- y[x <  0.4]           + 2
  y[x >= 0.4 & x < 0.6] <- y[x >= 0.4 & x < 0.6] - 2
  y[x >= 0.6 & x < 0.8] <- y[x >= 0.6 & x < 0.8] + 1
  y
}
f5 <- function(x) sin(5.5 * pi * x) -
  4 * sgn(x[59] - x) - 2 * sgn(x[78] - x) + 3 * sgn(x[180] - x) - 1.75
f6 <- function(x) 2 * sin(4 * pi * x) - 6 * abs(x - 0.4)^(3 / 10) - sgn(0.7 - x)

test_funcs <- paste0("f", 1:n_functs)


# -----------------------------------------------------------------------------
# 6. Start the MATLAB server (needed by the MATLAB-based methods)
# -----------------------------------------------------------------------------
setwd(MATLAB_SRS)
Matlab$startServer(matlab = MATLAB_BIN)
matlab <- Matlab()
isOpen <- open(matlab, trials = 30, interval = 2)  # MATLAB takes ~20-30 s to boot
if (isTRUE(isOpen)) message("MATLAB server is open")


# -----------------------------------------------------------------------------
# 7. Main simulation loop
# -----------------------------------------------------------------------------
rmse3d <- array(0, dim = c(n_functs, n_methods, n_sim))
mxdv3d <- array(0, dim = c(n_functs, n_methods, n_sim))
results_list <- list()
estimates    <- list()          # fitted curves for replication `estimate_rep`
counter <- 1L

set.seed(123)

for (index_fun in seq_along(test_funcs)) {
  f_name <- test_funcs[index_fun]
  cat("test function =", f_name, "\n")

  f_fun  <- get(f_name)
  x      <- seq(0, 1, length.out = n_points)
  y_true <- f_fun(x)
  sig    <- sqrt(var(y_true)) / SNR

  # Knot sets consumed by the method wrappers
  knots1 <- default.knots(x, 12)
  kn     <- default.knots(x, n_points / 4)
  kn.var <- default.knots(kn, 30)

  for (sim_id in seq_len(n_sim)) {
    cat("  simulation run =", sim_id, "\n")
    yy <- as.vector(y_true + rnorm(n_points, mean = 0, sd = sig))

    # For the chosen replication, collect every method's fitted curve
    keep_fit <- (sim_id == estimate_rep)
    if (keep_fit) {
      pred_mat <- matrix(NA_real_, nrow = n_points, ncol = n_methods)
      colnames(pred_mat) <- methods
    }

    for (index_meth in seq_along(methods)) {
      m_name <- methods[index_meth]
      cat("    method =", m_name, "\n")

      m_fun  <- get(m_name)
      y_pred <- m_fun(yy, x)
      crit   <- compute_criteria(y_pred, y_true)

      results_list[[counter]] <- data.frame(
        test_function = f_name,
        sim_id        = sim_id,
        method        = m_name,
        RMSE          = unname(crit["RMSE"]),
        MXDV          = unname(crit["MXDV"])
      )
      rmse3d[index_fun, index_meth, sim_id] <- crit["RMSE"]
      mxdv3d[index_fun, index_meth, sim_id] <- crit["MXDV"]
      if (keep_fit) pred_mat[, index_meth] <- y_pred
      counter <- counter + 1L
    }

    # Store the noisy data, the truth and all fits for the estimate panels
    if (keep_fit)
      estimates[[f_name]] <- data.frame(
        x = x, noisy = yy, true = y_true, pred_mat, check.names = FALSE)
  }
}

close(matlab)


# -----------------------------------------------------------------------------
# 8. Assemble and save results
# -----------------------------------------------------------------------------
results_df <- na.omit(dplyr::bind_rows(results_list))
str(results_df)

out_file <- file.path(RESULTS_DIR,
                      paste0("simulation_results_nsim", n_sim, "_SNR", SNR, ".RData"))

save(rmse3d, mxdv3d, results_df, estimates,
     methods, test_funcs, n_sim, n_points, SNR, estimate_rep,
     file = out_file)

message("Saved: ", out_file)
