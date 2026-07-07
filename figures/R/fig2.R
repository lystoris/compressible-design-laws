# Fig 2 (methodology) — why the ABSOLUTE metric matters: a ratio-to-black-box hides incompressibility
# when the black box itself fails. Run: Rscript figures/R/fig2.R
source("figures/R/theme_pub.R")

# anchor_table.csv: round-03 was not reconstructed for this repo;
# this is the committed, audited round-03-t3 anchor output, copied in verbatim.
anch <- read_csv(file.path(DATA_DIR,"anchor_table.csv"), show_col_types=FALSE)
pop  <- read_csv(file.path(DATA_DIR,"fig3_population.csv"), show_col_types=FALSE)

# (a) per real anchor: absolute R2det vs ratio-to-black-box (Lambda) -- ratio stays high when both fail
aa <- anch %>% select(anchor, `absolute R2det`=r2_law, `ratio to black box`=lambda_clipped) %>%
  tidyr::pivot_longer(-anchor, names_to="metric", values_to="value") %>%
  mutate(anchor=recode(anchor, betacarotene="β-carotene"))
pa <- ggplot(aa, aes(reorder(anchor, value), value, fill=metric)) +
  geom_col(position=position_dodge(0.7), width=0.62) +
  geom_hline(yintercept=0.7, linetype="dashed", colour="grey55") +
  scale_fill_manual(values=c(`absolute R2det`=COL_EFF, `ratio to black box`=COL_NOM), name=NULL,
                     breaks=c("absolute R2det","ratio to black box"),
                     labels=c(expression("absolute "*R[det]^2), "ratio to black box")) +
  labs(x=NULL, y="compressibility metric",
       title="A ratio metric hides incompressibility", tag="a") +
  theme_pub() + theme(legend.position=c(0.02,0.98), legend.justification=c(0,1),
                      axis.text.x=element_text(angle=20, hjust=1))

# (b) population: where the black box is weak (low r2_bb), the ratio is uninformative
pb <- ggplot(pop, aes(r2_bb, r2_law, colour=domain)) +
  geom_abline(slope=1, intercept=0, colour="grey75", linetype="dashed") +
  geom_point(alpha=0.75, size=1.1) +
  scale_colour_manual(values=PAL_DOMAIN, name="domain") +
  annotate("text", x=0.15, y=0.9, hjust=0, size=2,
           label="ratio = law / black box\nmisleads when black box is weak", colour="grey35") +
  labs(x=expression("black-box ceiling ("*R[det]^2*")"),
       y=expression("law compressibility ("*R[det]^2*")"),
       title="Use the absolute metric", tag="b") +
  theme_pub() + theme(legend.position=c(0.98,0.02), legend.justification=c(1,0))

fig <- pa | pb
save_fig(fig, "Fig2_metric", width_mm=174, height_mm=70)
