# Shared theme / palette / IO for v2 figures (Cell Systems style, multi-panel)
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(patchwork)
  library(ggrepel); library(scales); library(jsonlite)
})

DATA_DIR <- file.path("figures","data")   # run from repo root: Rscript figures/R/make_all.R
OUT_DIR  <- file.path("figures")

PAL_DOMAIN <- c(metabolic="#1B9E77", protein="#D95F02", simulator="#7570B3", synthetic="#7570B3")
COL_EFF <- "#2C7FB5"; COL_NOM <- "#B0670B"; COL_BB <- "#9A9A9A"; COL_HL <- "#222222"
COMPRESS_LINE <- 0.7

theme_pub <- function(base=7) {
  theme_classic(base_size=base) +
    theme(
      axis.text=element_text(colour="black", size=base-0.5),
      axis.title=element_text(size=base),
      plot.title=element_text(size=base+1, face="bold"),
      plot.subtitle=element_text(size=base-0.5, colour="grey30"),
      plot.tag=element_text(size=base+3, face="bold"),
      legend.key.size=unit(3,"mm"), legend.text=element_text(size=base-1),
      legend.title=element_text(size=base-0.5),
      strip.background=element_blank(), strip.text=element_text(size=base, face="bold"),
      plot.margin=margin(3,3,3,3)
    )
}

# Cell column widths (mm): single 85, 1.5-col 114, double 174
save_fig <- function(plot, name, width_mm, height_mm) {
  ggsave(file.path(OUT_DIR, paste0(name,".png")), plot,
         width=width_mm, height=height_mm, units="mm",
         dpi=600, device="png", bg="white", limitsize=FALSE)
  cat("saved", name, sprintf("(%.0fx%.0f mm)\n", width_mm, height_mm))
}
