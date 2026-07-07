# Fig S (related to Fig 3) — the independent effective-d estimator recovers TRUE effective d.
# Calibration of the RF participation-ratio estimator used in the population panel.
# Build data first: /usr/bin/python3 scripts/build_calibration_data.py
# Run: Rscript figures/R/figS_calibration.R
source("figures/R/theme_pub.R")

cal <- read_csv(file.path(DATA_DIR, "figS_calibration.csv"), show_col_types = FALSE)
fam_lab <- c(additive = "additive", michaelis_menten = "Michaelis-Menten", random_gp = "random GP")
d_nom <- max(cal$d_nom)

p <- ggplot(cal, aes(d_eff_set, eff_d_est, group = d_eff_set)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_hline(yintercept = d_nom, linetype = "dotted", colour = COL_NOM, linewidth = 0.3) +
  geom_boxplot(width = 0.9, outlier.size = 0.3, linewidth = 0.3, fill = "grey92") +
  facet_wrap(~family, nrow = 1, labeller = as_labeller(fam_lab)) +
  annotate("text", x = 2, y = d_nom, vjust = -0.5, hjust = 0, size = 2, parse = TRUE,
           colour = COL_NOM, label = paste0("nominal~italic(d)==", d_nom)) +
  scale_x_continuous(breaks = sort(unique(cal$d_eff_set))) +
  labs(x = expression("true (set) effective " * italic(d)),
       y = expression(atop("estimated effective " * italic(d), "(RF participation ratio)")),
       subtitle = expression("dashed = identity (perfect recovery); dotted = nominal " * italic(d))) +
  theme_pub()

save_fig(p, "FigS2_calibration", width_mm = 140, height_mm = 62)
