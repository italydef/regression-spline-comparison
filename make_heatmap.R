# =============================================================================
# make_heatmap.R
# -----------------------------------------------------------------------------
# Reproduces the comparison heatmap (Figure 3.3): two side-by-side blocks
# (RMSE and MXDV), rows = methods, columns = test functions. Each cell shows
# the mean criterion over the n_sim runs with its within-function rank in
# parentheses; methods are ordered by overall average rank.
#
# INPUT : INPUT/Simulation/simulation_results_nsim<n_sim>_SNR<SNR>.RData  (written by run_simulation.R)
# OUTPUT: <OUTPUT_DIR>/heatmap_ranking.pdf
#
# This script only reads results; it does not re-run the simulation.
# =============================================================================

library(ggplot2)     # heatmap tiles (geom_tile) + theming
library(dplyr)       # ranking / summary (%>%, group_by, summarise)
library(tidyr)       # pivot_longer() to long form
library(patchwork)   # compose panels (p_rmse + p_mxdv)
  
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
# provides: results_df, methods, test_funcs, ...

# -----------------------------------------------------------------------------
# 1. Mean per (function, method) and within-function ranks (lower = better)
# -----------------------------------------------------------------------------
summary_df <- results_df %>%
  group_by(test_function, method) %>%
  summarise(mean_RMSE = mean(RMSE),
            mean_MXDV = mean(MXDV),
            .groups = "drop") %>%
  group_by(test_function) %>%
  mutate(rank_RMSE = rank(mean_RMSE, ties.method = "average"),
         rank_MXDV = rank(mean_MXDV, ties.method = "average")) %>%
  ungroup()


# -----------------------------------------------------------------------------
# 2. Order methods by overall average rank (best method at the top)
# -----------------------------------------------------------------------------
method_order <- summary_df %>%
  group_by(method) %>%
  summarise(avg_rank = mean((rank_RMSE + rank_MXDV) / 2), .groups = "drop") %>%
  arrange(avg_rank) %>%
  pull(method)

# ggplot draws the first factor level at the bottom, so reverse to put the
# best-ranked method at the top of the heatmap.
summary_df <- summary_df %>%
  mutate(method        = factor(method, levels = rev(method_order)),
         test_function = factor(test_function, levels = test_funcs))


# -----------------------------------------------------------------------------
# 3. Long format: one row per (function, method, criterion)
# -----------------------------------------------------------------------------
plot_df <- summary_df %>%
  pivot_longer(cols = c(mean_RMSE, mean_MXDV, rank_RMSE, rank_MXDV),
               names_to = c(".value", "criterion"),
               names_pattern = "(mean|rank)_(.*)")


# -----------------------------------------------------------------------------
# 4. Heatmap for one criterion (its own colour scale)
# -----------------------------------------------------------------------------
make_heatmap <- function(data, crit_name) {
  ggplot(filter(data, criterion == crit_name),
         aes(test_function, method, fill = mean)) +
    geom_tile(color = "lightgray", linewidth = 0.4) +
    geom_text(aes(label = paste0(sprintf("%.3f", mean),
                                 "\n(", round(rank), ")")),
              size = 4.5, fontface = "bold", colour = "black",
              lineheight = 0.9) +
    # begin/end restrict plasma to its lighter half so black text stays readable
    scale_fill_viridis_c(option = "plasma", direction = -1,
                         begin = 0.4, end = 1) +
    labs(title = crit_name, x = "Test function", y = "Method", fill = "Mean") +
    theme_minimal(base_size = 14) +
    theme(panel.grid   = element_blank(),
          axis.text.x  = element_text(angle = 45, hjust = 1, size = 13),
          axis.text.y  = element_text(size = 13),
          axis.title   = element_text(size = 15),
          legend.title = element_text(size = 13),
          legend.text  = element_text(size = 12),
          plot.title   = element_text(face = "bold", size = 18))
}

p_rmse <- make_heatmap(plot_df, "RMSE")
p_mxdv <- make_heatmap(plot_df, "MXDV")

final_plot <- p_rmse + p_mxdv    # separate fill scales, side by side

out_pdf <- file.path(OUTPUT_DIR,
                     paste0("heatmap_ranking_nsim", n_sim, "_SNR", SNR, ".pdf"))
ggsave(out_pdf, final_plot, width = 14, height = 8, device = "pdf")
message("Saved: ", out_pdf)

