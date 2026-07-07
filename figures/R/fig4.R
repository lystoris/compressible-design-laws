# Fig 4 (CENTERPIECE) — the population law: compressibility vs INDEPENDENT effective d across ~82 landscapes.
# Run from repo root: Rscript figures/R/fig4.R
source("figures/R/theme_pub.R")

pop <- read_csv(file.path(DATA_DIR, "fig3_population.csv"), show_col_types = FALSE)
part <- read_csv(file.path(DATA_DIR, "fig3_partials.csv"), show_col_types = FALSE)
pseu <- read_csv(file.path(DATA_DIR, "fig3_pseudorep.csv"), show_col_types = FALSE)

lab_ids <- c(
     "beta-carotene", "gb1", "poelwijk", "isoprenol", "avgfp",
     "RL20_AQUAE_Tsuboyama_2023_1GYZ", "CBPA2_HUMAN_Tsuboyama_2023_1O6X"
)
pop <- pop %>% mutate(lab = ifelse(id %in% lab_ids,
     recode(id,
          RL20_AQUAE_Tsuboyama_2023_1GYZ = "RL20 (nom 59)",
          CBPA2_HUMAN_Tsuboyama_2023_1O6X = "CBPA2 (nom 72)",
          `beta-carotene` = "β-carotene", gb1 = "GB1", poelwijk = "Poelwijk",
          isoprenol = "isoprenol", avgfp = "avGFP"
     ), NA_character_
))

# (a) compressibility vs EFFECTIVE d  -- the hero panel
# regimes shown as a continuum (none / partial / full), not a single pass/fail line,
# so the LOESS gradient (the result) leads and no single point carries the argument.
pa <- ggplot(pop, aes(eff_d, r2_law)) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.4, ymax = 0.7, fill = "grey60", alpha = 0.10) +
     annotate("text",
          x = max(pop$eff_d), y = 0.12,
          label = "none", hjust = 1, vjust = 0.5, size = 1.8,
          colour = "grey55", fontface = "italic"
     ) +
     geom_smooth(method = "loess", se = TRUE, colour = COL_EFF, fill = COL_EFF, alpha = 0.18, linewidth = 1.0, span = 1) +
     annotate("text",
          x = max(pop$eff_d), y = -0.42, hjust = 1, vjust = 0.5, size = 2.1, colour = "black",
          parse = TRUE, label = "Spearman~rho == -0.69 ~~ (n == 82)"
     ) +
     geom_point(aes(colour = domain, size = nom_d), alpha = 0.6) +
     geom_text_repel(aes(label = lab),
          size = 2, min.segment.length = 0, box.padding = 0.4,
          max.overlaps = Inf, seg.colour = "grey60", colour = COL_HL, fontface = "bold"
     ) +
     scale_x_log10() +
     scale_colour_manual(values = PAL_DOMAIN, name = "domain") +
     scale_size_continuous(range = c(0.6, 3.6), name = expression("nominal " * italic(d)), trans = "log10", breaks = c(3, 10, 30, 120)) +
     labs(
          x = "effective dimensionality (independent)", y = expression("compressibility (" * R[det]^2 * ")"),
          title = expression(bold("Effective ") * bolditalic(d) * bold(" predicts compressibility")), tag = "a"
     ) +
     theme_pub() +
     theme(legend.position = "top", legend.box = "horizontal")

# (b) compressibility vs NOMINAL d  -- the contrast: no clean trend
pb <- ggplot(pop, aes(nom_d, r2_law)) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.4, ymax = 0.7, fill = "grey60", alpha = 0.10) +
     geom_smooth(method = "loess", se = TRUE, colour = COL_NOM, fill = COL_NOM, alpha = 0.18, linewidth = 1.0, span = 1) +
     annotate("text",
          x = max(pop$nom_d), y = -0.42, hjust = 1, vjust = 0.5, size = 2.1, colour = "black",
          parse = TRUE, label = "Spearman~rho == -0.15 ~~ (n == 82)"
     ) +
     geom_point(aes(colour = domain), size = 1.2, alpha = 0.6) +
     scale_x_log10() +
     scale_colour_manual(values = PAL_DOMAIN, guide = "none") +
     labs(
          x = "nominal dimensionality (knob count)", y = expression(R[det]^2),
          title = expression(bold("Nominal ") * bolditalic(d) * bold(" does not")), tag = "b"
     ) +
     theme_pub()

# (c) partial-correlation forest (both streams, 3 trajectories each)
part <- part %>% mutate(
     axis = factor(axis, levels = c("nominal d", "effective d")),
     grp = paste(stream, axis)
)
pc <- ggplot(part, aes(partial, interaction(axis, stream), colour = axis)) +
     geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
     stat_summary(fun = mean, geom = "point", size = 2) +
     geom_point(alpha = 0.5, size = 1.1, position = position_nudge(y = 0.12)) +
     scale_colour_manual(values = c("effective d" = COL_EFF, "nominal d" = COL_NOM), guide = "none") +
     scale_y_discrete(labels = function(x) {
          lab <- sub("\\..*", "", x)
          parse(text = ifelse(grepl("effective", lab), "'effective'~italic(d)", "'nominal'~italic(d)"))
     }) +
     facet_grid(stream ~ ., scales = "free_y", space = "free_y", switch = "y") +
     labs(
          x = "partial Spearman ρ", y = NULL,
          title = expression(bold("Effective ") * bolditalic(d) * bold(" beats nominal ") * bolditalic(d)), tag = "c"
     ) +
     theme_pub() +
     theme(
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90)
     )

# (d) ProteinGym pseudoreplication-corrected scale (short labels for the narrow right column)
lev_short <- c(
     `all 69 (pooled)` = "pooled", `Tsuboyama within-protocol` = "within-protocol",
     `15 other studies` = "cross-study", `study-collapsed` = "study-collapsed"
)
pseu <- pseu %>%
     mutate(level = recode(level, !!!lev_short)) %>%
     mutate(level = factor(level, levels = rev(unname(lev_short))))
pd <- ggplot(pseu, aes(rho, level)) +
     geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
     geom_segment(aes(x = 0, xend = rho, yend = level), colour = "grey70", linewidth = 0.4) +
     geom_point(aes(size = n), colour = COL_EFF) +
     geom_text(aes(label = sprintf("%.2f (n=%d)", rho, n)), vjust = -1.1, size = 1.9) +
     scale_size_continuous(range = c(1.3, 3.5), guide = "none") +
     scale_x_continuous(limits = c(-1.12, 0.15)) +
     labs(
          x = expression(rho * "(compress, effective " * italic(d) * ")"), y = NULL,
          title = "Robust to pseudoreplication", tag = "d"
     ) +
     theme_pub()

# Layout (D): 2x2 unequal columns (wide a/b left, narrow c/d right). Shared domain/size legend
# collected to the FIGURE TOP; c's colour legend dropped as redundant with its y-axis rows.
# free(pc, "l"): c/d have mismatched left-margin widths (d's level labels are long), which breaks d's
# axis corner; release c's left from alignment (patchwork 1.3+). free(pd) hits a patchwork bug; free(pc) is fine.
# guide_area() puts the collected legend in a dedicated strip ABOVE the panel titles.
panels <- ((pa / pb) | (free(pc, side = "l") / pd)) + plot_layout(widths = c(1.2, 1))
fig <- guide_area() / panels +
     plot_layout(guides = "collect", heights = c(0.08, 1)) &
     theme(legend.position = "top", legend.box = "horizontal")
save_fig(fig, "Fig4_population_law", width_mm = 150, height_mm = 145)
