# =============================================================================
# make_boxplots.R
# -----------------------------------------------------------------------------
# Reproduces the boxplot figure (Figure 3.2): for each test function,
# side-by-side boxplots of RMSE and MXDV across the n_sim runs, one box per
# method.
#
# INPUT : INPUT/Simulation/simulation_results_nsim<n_sim>_SNR<SNR>.RData  (written by run_simulation.R)
# OUTPUT: <OUTPUT_DIR>/boxplots_SNR3.pdf
#
# This script only reads results; it does not re-run the simulation.
# =============================================================================

library(ggplot2)     # boxplots (geom_boxplot) + theming
library(dplyr)       # data wrangling (%>%, mutate, filter)
library(patchwork)   # compose panels (wrap_plots, p_rmse | p_mxdv)

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
# provides: results_df, methods, test_funcs, n_sim, SNR, ...

# Fix the method order on the x-axis to the manuscript order
results_df <- results_df %>%
  mutate(method = factor(method, levels = methods))

method_colors <- scales::hue_pal()(length(methods))
names(method_colors) <- methods

# -----------------------------------------------------------------------------
# Helper: the two boxplot panels (RMSE, MXDV) for one test function
# -----------------------------------------------------------------------------
boxplot_row <- function(df, index_fun) {
  f_name <- paste0("f", index_fun)
  d      <- filter(df, test_function == f_name)

  base <- list(
    geom_boxplot(),
    scale_fill_manual(values = method_colors, guide = "none"),
    theme_minimal(base_size = 10),
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  p_rmse <- ggplot(d, aes(method, RMSE, fill = method)) + base +
    labs(title = paste0("function f", index_fun, ": SNR=", SNR, ", RMSE"),
         x = "Method", y = paste0("RMSE over ", n_sim, " runs"))

  p_mxdv <- ggplot(d, aes(method, MXDV, fill = method)) + base +
    labs(title = paste0("function f", index_fun, ": SNR=", SNR, ", MXDV"),
         x = "Method", y = paste0("MXDV over ", n_sim, " runs"))

  p_rmse | p_mxdv
}

# -----------------------------------------------------------------------------
# Build one row per function and stack them
# -----------------------------------------------------------------------------
rows       <- lapply(seq_along(test_funcs), function(i) boxplot_row(results_df, i))
final_plot <- wrap_plots(rows, ncol = 1)

out_pdf <- file.path(OUTPUT_DIR,
                     paste0("boxplots_nsim", n_sim, "_SNR", SNR, ".pdf"))
ggsave(out_pdf, final_plot, width = 12, height = 18)
message("Saved: ", out_pdf)