# Render all figures. Run from repo root: Rscript figures/R/make_all.R
# (Fig 1 is a BioRender schematic; not rendered here.)
# Figure data must be prepped first: /usr/bin/python3 scripts/build_figure_data.py
#   and, for the calibration SI: /usr/bin/python3 scripts/build_calibration_data.py
for (f in c("fig2","fig3","fig4","fig5","fig6","figS_calibration","figS_sweep_phase")) {
  cat("\n==== ", f, " ====\n")
  source(file.path("figures","R", paste0(f,".R")))
}
cat("\nAll figures rendered to figures/ (600 dpi PNG; Figs 2-6 + 2 SI).\n")
