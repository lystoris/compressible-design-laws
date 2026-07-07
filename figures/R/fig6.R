# Fig 6 — mechanistic simulator: low effective d partially compresses AND the law out-selects the
# black box. Run: Rscript figures/R/fig6.R
source("figures/R/theme_pub.R")

# sim_Ncurve.csv: this repo's regenerated N-curve (scripts/run_simulator.py -> results/sim_Ncurve.csv)
nc  <- read_csv(file.path("results","sim_Ncurve.csv"),
                show_col_types=FALSE)
sel <- read_csv(file.path(DATA_DIR,"fig5_selection.csv"), show_col_types=FALSE)
sens<- read_csv(file.path(DATA_DIR,"fig5_enzyme_sens.csv"), show_col_types=FALSE)

# (a) simulator partially compresses: law R2det vs tuned black box across training size
na <- nc %>% select(N, law=r2_law, `black box`=r2_blackbox) %>%
  tidyr::pivot_longer(-N, names_to="model", values_to="r2")
pa <- ggplot(na, aes(N, r2, colour=model)) +
  geom_line(linewidth=0.6) + geom_point(size=1.2) +
  scale_colour_manual(values=c(law=COL_EFF, `black box`=COL_BB), name=NULL) +
  scale_x_log10() + ylim(0,1) +
  labs(x=expression("training size "*italic(N)), y=expression(R[det]^2),
       title="Simulator partially compresses", tag="a") +
  theme_pub() + theme(legend.position=c(0.98,0.02), legend.justification=c(1,0))

# (b) selection payoff: law recovers more of the top-100 region than the black box, at every N
sb <- sel %>% select(N, law=law_top100, `black box`=bb_top100) %>%
  tidyr::pivot_longer(-N, names_to="model", values_to="overlap")
pb <- ggplot(sb, aes(N, overlap, colour=model)) +
  geom_line(linewidth=0.6) + geom_point(size=1.2) +
  scale_colour_manual(values=c(law=COL_EFF, `black box`=COL_BB), name=NULL) +
  scale_x_log10() +
  labs(x=expression("training size "*italic(N)), y="top-100 designs recovered /100",
       title="The law out-selects the black box", tag="b") +
  theme_pub() + theme(legend.position=c(0.02,0.98), legend.justification=c(0,1))

# (c) effective vs nominal enzymes: 3 of 7 carry signal (EA, EC, EG) -> flux ~ EA*EC/EG
sens <- sens %>% mutate(enzyme=factor(enzyme, levels=enzyme),
                        key=enzyme %in% c("EA","EC","EG"))
pc <- ggplot(sens, aes(enzyme, abs_spearman, fill=key)) +
  geom_col(width=0.65) +
  scale_fill_manual(values=c(`TRUE`=COL_EFF,`FALSE`="grey80"), guide="none") +
  annotate("text", x=5.5, y=max(sens$abs_spearman)*0.9, size=2.2, colour="grey30",
           label="law: flux ∝ EA·EC / EG\n3 of 7 enzymes effective") +
  labs(x="enzyme (design knob)", y="|Spearman| with flux",
       title="Effective ≪ nominal enzymes", tag="c") + theme_pub()

fig <- pa | pb | pc
save_fig(fig, "Fig6_simulator_selection", width_mm=174, height_mm=62)
