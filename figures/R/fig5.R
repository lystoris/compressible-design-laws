# Fig 5 — necessary but not sufficient: low effective d sets the ceiling; clean compressors are rare.
# Run: Rscript figures/R/fig5.R
source("figures/R/theme_pub.R")

pop <- read_csv(file.path(DATA_DIR, "fig3_population.csv"), show_col_types = FALSE)
oof <- read_csv(file.path(DATA_DIR, "betacarotene_oof.csv"), show_col_types = FALSE) %>%
  mutate(split = ifelse(grepl("extrap", split), "extrapolation", "interpolation"))
probe <- read_csv(file.path(DATA_DIR, "fig4_poelwijk_probe.csv"), show_col_types = FALSE)

# (a) among LOW-effective-d landscapes, compressibility spans partial->full; few cross 0.7
low <- pop %>% filter(eff_d <= 6)
frac_full <- mean(low$r2_law >= 0.7)
pa <- ggplot(low, aes(r2_law)) +
  geom_histogram(binwidth = 0.1, fill = COL_EFF, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.7, linetype = "dashed", colour = "grey40") +
  annotate("text",
    x = -0.55, y = Inf, vjust = 1.4, hjust = 0, size = 2.1,
    label = sprintf("only %.0f%% of low-\neff-d landscapes\nreach 0.7", 100 * frac_full)
  ) +
  labs(
    x = expression("compressibility (" * R[det]^2 * ")"), y = expression("landscapes (effective " * italic(d) * " ≤ 6)"),
    title = expression(bold("Low effective ") * bolditalic(d) * bold(" ≠ guaranteed law")), tag = "a"
  ) +
  theme_pub()

# (b) CLEAN compressor: beta-carotene out-of-fold fit (R2det ~0.92), law C + C/B + A*C
r2det <- 1 - sum((oof$observed - oof$predicted)^2) / sum((oof$observed - mean(oof$observed))^2)
pb <- ggplot(oof, aes(observed, predicted, colour = split)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = "dashed") +
  geom_point(size = 1.3, alpha = 0.85) +
  scale_colour_manual(values = c(interpolation = COL_EFF, extrapolation = COL_NOM), name = NULL) +
  annotate("text",
    x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 2.3,
    label = "β-carotene: titre ~ C + C/B + A·C"
  ) +
  annotate("text",
    x = -Inf, y = Inf, hjust = -0.1, vjust = 3.1, size = 2.3, parse = TRUE,
    label = sprintf("R[det]^2 == %.2f", r2det)
  ) +
  labs(
    x = "observed titre", y = "predicted (out-of-fold)",
    title = "A clean compressor", tag = "b"
  ) +
  theme_pub() +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))

# (c) PARTIAL compressor: Poelwijk -- learnable (bb 0.93) but short-law-capped; link refuted
probe <- probe %>% mutate(transform = factor(transform, levels = transform))
# short x labels (full meaning in legend); "interpretable" = binary encoding, linear+interaction terms only
x_labs <- c(
  "raw" = "raw", "log" = "log", "rank-normal" = "rank",
  "interpretable (binary, linear+interaction)" = "interpretable"
)
pc <- ggplot(probe, aes(transform, r2_law)) +
  geom_col(fill = COL_EFF, width = 0.6) +
  geom_hline(aes(yintercept = max(r2_bb)), linetype = "dashed", colour = COL_BB) +
  annotate("text",
    x = 0.6, y = max(probe$r2_bb) + 0.03, hjust = 0, size = 2,
    label = "black-box ceiling 0.93", colour = "grey40"
  ) +
  geom_text(aes(label = sprintf("%.2f", r2_law)), vjust = -0.4, size = 2) +
  scale_x_discrete(labels = x_labs) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = NULL, y = expression("Poelwijk law " * R[det]^2),
    title = "A partial compressor", tag = "c"
  ) +
  theme_pub()

fig <- (pa | pb | pc) + plot_layout(widths = c(1, 1, 1))
save_fig(fig, "Fig5_necessary_not_sufficient", width_mm = 170, height_mm = 60)
