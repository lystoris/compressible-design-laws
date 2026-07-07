# Fig S (related to Fig 2) — controlled sweep phase diagram: compressibility falls with design
# dimensionality ACROSS noise, sample size, and encoding (robustness of Fig 2a).
# Run from repo root: Rscript figures/R/figS_sweep_phase.R
source("figures/R/theme_pub.R")

# sweep_grid.csv: this repo's regenerated sweep (scripts/run_sweep.py --traj t3 -> results/sweep_grid.csv)
sweep <- read_csv(file.path("results", "sweep_grid.csv"),
                  show_col_types = FALSE)
sweep <- sweep %>% mutate(comp = pmax(r2_law, 0))   # floor at 0 for display

# one row per robustness factor: median compressibility vs d, by epistasis k, faceted by the factor
row_plot <- function(df, facet_var, lab_fn, tag) {
  agg <- df %>% group_by(d, k, .data[[facet_var]]) %>%
    summarise(r2 = median(comp, na.rm = TRUE), .groups = "drop")
  agg$facet_lab <- lab_fn(agg[[facet_var]])
  agg$facet_lab <- factor(agg$facet_lab,
                          levels = unique(agg$facet_lab[order(agg[[facet_var]])]))
  ggplot(agg, aes(d, r2, colour = factor(k), group = k)) +
    geom_line(linewidth = 0.5) + geom_point(size = 0.9) +
    facet_wrap(~facet_lab, nrow = 1, labeller = label_parsed) +
    scale_x_log10() + ylim(0, 1) +
    scale_colour_brewer(palette = "YlOrRd", name = expression("epistasis " * italic(k))) +
    labs(x = expression("design dimensionality " * italic(d)),
         y = expression("compressibility (" * R[det]^2 * ")"), tag = tag) +
    theme_pub()
}

pa <- row_plot(sweep, "sigma", function(v) paste0("sigma==", v), "a")   # measurement noise
pb <- row_plot(sweep, "N", function(v) paste0("italic(N)==", v), "b")   # sample size
pc <- row_plot(sweep, "encoding", function(v) paste0("'", v, "'"), "c") # encoding

fig <- (pa / pb / pc) + plot_layout(guides = "collect") & theme(legend.position = "right")
save_fig(fig, "FigS1_sweep_phase", width_mm = 140, height_mm = 150)
