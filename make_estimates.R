# =============================================================================
# make_estimates.R
# -----------------------------------------------------------------------------
# Reproduces the per-method estimate figures (Figures 3.4 / 3.5): for one
# simulated noisy data set per test function, a panel of the noisy data plus
# one panel per method showing the fitted curve (blue) against the true
# function (red).
#
# INPUT : INPUT/Simulation/simulation_results_nsim<n_sim>_SNR<SNR>.RData  (written by run_simulation.R)
#         uses the `estimates` list, i.e. the fitted curves saved for
#         replication `estimate_rep`.
# OUTPUT: one file per test function,
#         <OUTPUT_DIR>/estimates_f1.pdf ... estimates_f6.pdf
#
# This script only reads results; it does not re-run the simulation.
# =============================================================================

library(ggplot2)     # per-method estimate panels
library(patchwork)   # compose panels (wrap_plots)

# --- Path: point this at the same folder run_simulation.R wrote to -----------
WORKING_DIR  <- getwd() 
INPUT_DIR<- file.path(WORKING_DIR,"INPUT","Simulation")
OUTPUT_DIR <- file.path(WORKING_DIR,"OUTPUT","Simulation")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

n_sim        <- 50    # replications per function
SNR          <- 3 
file_inp <- file.path(INPUT_DIR,
                      paste0("simulation_results_nsim", n_sim, "_SNR", SNR, ".RData"))
load(file_inp)
# provides: estimates, methods, test_funcs, ...

panel_theme <- theme_minimal(base_size = 8) +
  theme(plot.title = element_text(size = 9),
        axis.title = element_blank())


# -----------------------------------------------------------------------------
# Helper: the full panel set (noisy data + one panel per method) for a function
# -----------------------------------------------------------------------------
estimate_figure <- function(d, methods) {

  # Panel 1: the noisy data
  p_noisy <- ggplot(d, aes(x, noisy)) +
    geom_point(shape = 19, size = 0.5, color = "darkgrey") +
    labs(title = "noisy data") + panel_theme

  # One panel per method: true (red) with the fitted estimate (blue) on top
  method_panels <- lapply(methods, function(m) {
    ggplot(d, aes(x)) +
      geom_line(aes(y = true), color = "red", linewidth = 0.4) +
      geom_line(aes(y = .data[[m]]), color = "blue", linewidth = 0.5) +
      labs(title = paste0(m, " estimate")) + panel_theme
  })

  # 16 panels -> 4 x 4 grid (matches the layout used in the paper)
  wrap_plots(c(list(p_noisy), method_panels), ncol = 4)
}


# -----------------------------------------------------------------------------
# One figure per test function
# -----------------------------------------------------------------------------
for (f_name in test_funcs) {
  d   <- estimates[[f_name]]
  fig <- estimate_figure(d, methods)
  out <- file.path(OUTPUT_DIR,
                   paste0("estimates_", f_name, "_nsim", n_sim, "_SNR", SNR, ".pdf"))
  ggsave(out, fig, width = 11, height = 9, device = "pdf")
  message("Saved: ", out)
}

