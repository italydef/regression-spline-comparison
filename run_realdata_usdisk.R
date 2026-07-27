# =============================================================================
# run_realdata.R
# -----------------------------------------------------------------------------
# Real-data counterpart of run_simulation.R.
#
# Applies the SAME 15 methods ONCE
# to a real observed series (the USD/ISK monthly exchange rate, 2004.01-2015.12)
# and keeps only the fitted curves.
#
# OUTPUT
#   <OUTPUT_DIR>/usdisk_fits.RData   : fit_df (index, date, x, observed, <one
#                                      column per method>) + metadata
#   <OUTPUT_DIR>/usdisk_fits.pdf     : one panel per method, observed points
#                                      with the fitted curve overlaid
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
# smoothness-master: compiled Fortran shared object + its R interface
SMOOTH_DIR <- file.path(WORKING_DIR, "smoothness-master")
R_software_DIR <- file.path(WORKING_DIR  , "R_software")
AMPSCODE_DIR <- file.path(WORKING_DIR, "AMPScode")
MATLAB_SRS   <- file.path(WORKING_DIR, "adaptive_psplinesv06", "srs")

# The ONLY genuinely machine-specific setting: path to YOUR MATLAB executable.
MATLAB_BIN <- "/Applications/MATLAB_R2024b.app/bin/matlab"

# Where the real dataset lives and where output is written:
DATA_DIR   <- WORKING_DIR                       # folder holding usdisk.csv

INPUT_DIR<- file.path(WORKING_DIR,"INPUT","RealData")
file_inp <- file.path(INPUT_DIR,"usdisk.csv")

OUTPUT_DIR <- file.path(WORKING_DIR, "OUTPUT", "RealData")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Package -> method map.  Tags:
#   METHOD(S): fn()   = called DIRECTLY by that wrapper in METHODS_REAL_DATA.R
#   (via sourced X)   = reached only through a sourced helper file
# NB: EAS calls no library directly - pilotQV() is sourced from cp_source.R
#     (which itself uses fields::Tps); varcalib() is sourced from varcalib.R.
library(R.matlab)         # RC, BAPS: MATLAB bridge - setVariable()/evaluate()/getVariable()
library(GeDS)             # RC, SARS, AMP: NGeDS()  (also drives the hatsigma noise scale)
library(miscF)            # DMS: curve.polynomial.rjmcmc()  (pulls in R2jags -> rjags -> JAGS)
library(BARS)             # BARS: bars()
library(bsml)             # HAS: bsml()  (also supplies stdz())
library(genlasso)         # TF: trendfilter() / cv.trendfilter()
library(mgcv)             # EAPS: gam()  (primary fit + 'ad' fallback)
library(AdaptFitOS)       # FAPS: asp2()
library(cpr)              # CPR: cp() / bsplines()
library(sars)             # SARS: sars()
library(aspline)          # ASPL: aspline()  (adaptive-ridge knot selection)
library(freeknotsplines)  # FREEK: fit.search.numknots() / fitted.freekt()
library(SCHACE)           # SCHA: main.SCHACE()
library(pso)              # EAPS (via sourced OKPSPS/ok_gam.R): psoptim() in fit_gam_optim()
library(assist)           # HAS: cubic() reproducing-kernel basis
library(splines2)         # ASPL: bSpline() - MUST be attached; METHODS_REAL_DATA.R calls it un-namespaced
library(fields)           # EAS (via sourced cp_source.R): Tps() / fastTps()

# Data wrangling + plotting for this script only
library(ggplot2)          # fit panels + ggsave() to PDF
library(dplyr)            # data wrangling (bind_rows, etc.)
library(tidyr)            # reshape to long form for plotting

# -----------------------------------------------------------------------------
# 2. Source the method implementations and helper bases  (unchanged)
# -----------------------------------------------------------------------------
source(file.path(WORKING_DIR  , "default_knots.r"))
# Each name in METHODS must resolve to a function(yy, x).
source(file.path(WORKING_DIR  , "METHODS_REAL_DATA.R"))
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
# 3. Helper functions still needed by the method wrappers
# -----------------------------------------------------------------------------
# Two-valued sign: +1 for x > 0, -1 otherwise (note sgn(0) = -1, unlike base R).
sgn <- function(x) 2 * (x > 0) - 1

# -----------------------------------------------------------------------------
# 4. Settings and method dispatch (same 15 methods; no SNR / n_sim / noise)
# -----------------------------------------------------------------------------
n_methods <- 15
bdeg      <- 2

METHODS <- c("EAS", "RC", "HAS", "TF", "FAPS", "EAPS", "SARS",
             "ASPL", "FREEK", "BAPS", "DMS", "BARS", "CPR", "AMP", "SCHA")
methods <- METHODS
stopifnot(length(methods) == n_methods)

# -----------------------------------------------------------------------------
# 5. Read the real dataset  (replaces the test-function block)
# -----------------------------------------------------------------------------
isk <- read.csv(file_inp)

# rows 41:184 = 2004.01 ~ 2015.12  (144 monthly observations)
y_obs <- isk$price[41:184]
stopifnot(!anyNA(y_obs))

n_points <- length(y_obs) # 144

# Numeric design on [0,1] - what the spline methods expect (mirrors the sim).
x <- seq(0, 1, length.out = n_points)
yy <- as.vector(y_obs)

# A clean monthly date vector for the plot axis (Jan 2004 -> Dec 2015).
dates <- seq(as.Date("2004-01-01"), by = "month", length.out = n_points)

# Knot sets consumed by the method wrappers - rebuilt for the real length.
knots1 <- default.knots(x, 12)
kn     <- default.knots(x, n_points / 4)
kn.var <- default.knots(kn, 30)

# -----------------------------------------------------------------------------
# 6. Start the MATLAB server (needed by the MATLAB-based methods)
# -----------------------------------------------------------------------------
setwd(MATLAB_SRS)
Matlab$startServer(matlab = MATLAB_BIN)
matlab <- Matlab()
isOpen <- open(matlab, trials = 30, interval = 2)  # MATLAB takes ~20-30 s to boot
if (isTRUE(isOpen)) message("MATLAB server is open")
setwd(old_wd)

# -----------------------------------------------------------------------------
# 7. Fit every method ONCE to the observed series
# -----------------------------------------------------------------------------
set.seed(123)   # keep it reproducible

fit_mat <- matrix(NA_real_, nrow = n_points, ncol = n_methods)
colnames(fit_mat) <- methods

for (index_meth in seq_along(methods)) {
  m_name <- methods[index_meth]
  cat("method =", m_name, "\n")
  m_fun  <- get(m_name)

  # A failing method leaves NAs and does not abort the whole run.
  fit_mat[, index_meth] <- tryCatch(
    as.vector(m_fun(yy, x)),
    error = function(e) {
      message("  ", m_name, " failed: ", conditionMessage(e))
      rep(NA_real_, n_points)
    }
  )
}

close(matlab)

# -----------------------------------------------------------------------------
# 8. Assemble and save the fits
# -----------------------------------------------------------------------------
fit_df <- data.frame(index    = seq_len(n_points),
  date = dates, x = x, observed = yy, fit_mat, check.names = FALSE)
str(fit_df)

out_file <- file.path(OUTPUT_DIR, "usdisk_fits.RData")
save(fit_df, methods, n_points, file_inp, file = out_file)
message("Saved: ", out_file)

# -----------------------------------------------------------------------------
# Observed data plot
# -----------------------------------------------------------------------------
dates <- isk[41:184, 1]

# Convert "04-Jan" -> "2004-Jan"
date_labels <- sub("^([0-9]{2})-([A-Za-z]{3})$", "20\\1-\\2", dates)

plot_df <- data.frame(date = factor(date_labels, levels = date_labels),
           price = yy)

cp <- which(date_labels == "2008-Sep")

data_plot <- ggplot(plot_df, aes(x = date, y = price, group = 1)) +
  geom_point(shape = 16, colour = "black", size = 1.3) +
  geom_line(colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = cp, colour = "grey50", linetype = 2) +
  annotate("label", x = cp, y = 150, label = "Sep 2008", fill = "white") +
  scale_x_discrete(
    breaks = date_labels[seq(1, length(date_labels), by = 12)]) +
  coord_cartesian(ylim = c(40, 160)) + labs(title = "USD/ISK exchange rate",
    x = "Time", y = "USD/ISK") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(file.path(OUTPUT_DIR, "usdisk_data.pdf"), plot = data_plot,
       width = 9, height = 4.5)

# -----------------------------------------------------------------------------
# 9. Plot: one panel per method, observed points + fitted curve
# -----------------------------------------------------------------------------
plot_long <- fit_df %>%
  pivot_longer(cols = all_of(methods),
               names_to = "method", values_to = "fit") %>%
  mutate(method = factor(method, levels = methods))

fit_plot <- ggplot(plot_long, aes(x = date)) +
  geom_vline(xintercept = as.Date("2008-09-01"),
             colour = "grey70", linetype = 2) +           # 2008-09 reference
  geom_point(aes(y = observed), colour = "grey55", size = 0.5) +
  geom_line(aes(y = fit), colour = "firebrick", linewidth = 0.7) +
  facet_wrap(~ method, ncol = 3) +
  labs(title = "USD/ISK exchange rate",
       x = "Time", y = "USD/ISK") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

out_pdf <- file.path(OUTPUT_DIR, "usdisk_fits.pdf")
ggsave(out_pdf, fit_plot, width = 12, height = 14, device = "pdf")
message("Saved: ", out_pdf)
