# Fig 3 — controlled foundation: compressibility is dimensionality-driven, predicts selection,
# and is governed by EFFECTIVE (set) d, not nominal d. Run: Rscript figures/R/fig3.R
source("figures/R/theme_pub.R")

# sweep_grid.csv: this repo's regenerated sweep (scripts/run_sweep.py --traj t3 -> results/sweep_grid.csv)
sweep <- read_csv(file.path("results", "sweep_grid.csv"),
  show_col_types = FALSE
)
eff <- read_csv(file.path(DATA_DIR, "fig2_effdim_pooled.csv"), show_col_types = FALSE)
eta <- read_csv(file.path(DATA_DIR, "fig2_eta2.csv"), show_col_types = FALSE)

# compressibility floored at 0 for display (negative R2det = incompressible)
sweep <- sweep %>% mutate(comp = pmax(r2_law, 0))

# (a) compressibility falls with nominal dimensionality (by epistasis order k)
sa <- sweep %>%
  group_by(d, k) %>%
  summarise(r2 = median(comp, na.rm = TRUE), .groups = "drop")
pa <- ggplot(sa, aes(d, r2, colour = factor(k), group = k)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  scale_x_log10() +
  ylim(0, 1) +
  scale_colour_brewer(palette = "YlOrRd", name = expression("epistasis " * italic(k))) +
  labs(
    x = expression("design dimensionality " * italic(d)), y = expression("compressibility (" * R[det]^2 * ", floored)"),
    title = "Compressibility falls with dimensionality", tag = "a"
  ) +
  theme_pub() +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))

# (b) compressibility predicts selection success (audited cluster-robust rho)
sb <- sweep %>% mutate(sel = 1 - law_top1_regret)
pb <- ggplot(sb, aes(comp, sel)) +
  geom_point(colour = COL_EFF, alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", colour = COL_HL, fill = "grey70", linewidth = 0.5) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  annotate("text",
    x = 0.96, y = 0.06, hjust = 1, size = 2.3,
    label = "Spearman rho = 0.71 [0.37, 0.75]\ncluster-robust p = 3.3e-4"
  ) +
  labs(
    x = expression("compressibility (" * R[det]^2 * ", floored)"), y = "selection success (1 − top-1 regret)",
    title = "Compressibility predicts selection", tag = "b"
  ) +
  theme_pub()

# (c) DECOUPLING: compressibility vs effective d SET by decoys, faceted by generator family
# NOTE: fig2_effdim_pooled.csv / fig2_eta2.csv here come from this repo's reconstructed
# decoupling grid (results/effdim_grid_t1.csv), which emits family label "additive" where the
# original round-05 audited grid used "additive_linear" (see docs/reconciliation-notes.md R-3)
# -- both keys are mapped so this labeller works against either source.
fam_lab <- c(additive_linear = "additive", additive = "additive",
             michaelis_menten = "Michaelis-Menten", random_gp = "random GP")
pc <- ggplot(eff, aes(factor(d_eff), r2_law)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, fill = "grey92") +
  stat_summary(fun = median, geom = "line", aes(group = 1), colour = COL_EFF, linewidth = 0.6) +
  facet_wrap(~family, nrow = 1, labeller = as_labeller(fam_lab)) +
  labs(
    x = "effective dimensionality (set by inert decoys)", y = expression(R[det]^2),
    title = expression(bold("Compressibility tracks EFFECTIVE ") * bolditalic(d)), tag = "c"
  ) +
  theme_pub() +
  theme(strip.text = element_text(size = 6, face = "bold"))

# (d) compressibility is flat in nominal d: eta^2(effective) vs eta^2(nominal)
ed <- eta %>%
  tidyr::pivot_longer(c(eta2_eff, eta2_nom), names_to = "axis", values_to = "eta2") %>%
  mutate(axis = recode(axis, eta2_eff = "effective d", eta2_nom = "nominal d"))
pd2 <- ggplot(ed, aes(family, eta2, fill = axis)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_point(position = position_dodge(0.7), size = 0.6, colour = "grey30") +
  scale_fill_manual(values = c("effective d" = COL_EFF, "nominal d" = COL_NOM), name = NULL,
    breaks = c("effective d", "nominal d"),
    labels = expression("effective " * italic(d), "nominal " * italic(d))) +
  scale_x_discrete(labels = fam_lab) +
  labs(
    x = NULL, y = expression("variance explained (" * eta^2 * ")"),
    title = expression(bold("Nominal ") * bolditalic(d) * bold(" explains ~none")), tag = "d"
  ) +
  theme_pub() +
  theme(
    legend.position = c(0.98, 1.1), legend.justification = c(1, 1),
    axis.text.x = element_text(size = 6.8)
  )

fig <- (pa | pb) / (pc | pd2) + plot_layout(heights = c(1, 1))
save_fig(fig, "Fig3_controlled", width_mm = 140, height_mm = 140)
