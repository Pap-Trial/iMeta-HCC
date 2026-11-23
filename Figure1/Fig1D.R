source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

# Figure 1D - per-patient CopyKAT-derived CNV scores
#
# CNV score definition:
#   per-cell score = mean |CNV| across genomic bins of the CopyKAT CNV matrix,
#   divided by the expected diploid value at the same sequencing depth (LOESS
#   fitted on a per-patient, depth-matched immune-cell reference). A ratio of
#   1.0 therefore represents the diploid baseline.
#
# Normal_cell (diploid reference):
#   250 T/B/plasma/myeloid cells per patient, sampled to match the nFeatureRNA
#   quantile distribution of that patient's hepatocytes 


cnv <- read.csv("Figure1D_perpatient_CNV.csv", stringsAsFactors = FALSE)

cnv_levels <- c(
  "Hepa_Stress",
  "Hepa_Norm",
  "Hepa_Invasive",
  "Hepa_Meta",
  "Hepa_Imm",
  "Normal_cell"
)
cnv$group <- factor(cnv$group, levels = cnv_levels)

cnv_palette <- c(
  Hepa_Stress   = "#1F77B4",
  Hepa_Norm     = "#FF7F0E",
  Hepa_Invasive = "#2CA02C",
  Hepa_Meta     = "#D62728",
  Hepa_Imm      = "#9467BD",
  Normal_cell   = "#8B2500"
)

## patient-level tests vs the diploid reference
stats <- do.call(rbind, lapply(setdiff(cnv_levels, "Normal_cell"), function(g) {
  a <- cnv[cnv$group == g, , drop = FALSE]
  b <- cnv[cnv$group == "Normal_cell", , drop = FALSE]
  m <- merge(
    a[, c("patient", "median_CNV_ratio")],
    b[, c("patient", "median_CNV_ratio")],
    by = "patient",
    suffixes = c(".g", ".n")
  )
  p <- if (nrow(m) >= 6) {
    suppressWarnings(wilcox.test(
      m$median_CNV_ratio.g, m$median_CNV_ratio.n, paired = TRUE
    ))$p.value
  } else {
    suppressWarnings(wilcox.test(
      a$median_CNV_ratio, b$median_CNV_ratio
    ))$p.value
  }
  data.frame(
    group = g,
    n_patients = nrow(a),
    n_paired = nrow(m),
    p = p
  )
}))
stats$p_adj <- p.adjust(stats$p, method = "BH")
stats$label <- ifelse(
  stats$p_adj < 0.001, "***",
  ifelse(stats$p_adj < 0.01, "**",
         ifelse(stats$p_adj < 0.05, "*", "ns"))
)
print(stats)

p <- ggplot(cnv, aes(x = group, y = median_CNV_ratio, fill = group)) +
  geom_hline(
    yintercept = 1, linetype = "dotted", color = "grey40", size = 0.3
  ) +
  geom_boxplot(width = 0.6, outlier.shape = NA, size = 0.3) +
  scale_fill_manual(values = cnv_palette, name = NULL) +
  scale_x_discrete(
    labels = c(
      "Hepa_Stress", "Hepa_Norm", "Hepa_Invasive",
      "Hepa_Meta", "Hepa_Imm", "Normal cell"
    )
  ) +
  scale_y_continuous(limits = c(0.9, 3.05), breaks = c(1, 2, 3)) +
  labs(x = NULL, y = "CNV score") +
  theme_pub(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## brackets: five reference comparisons only
ymax <- max(cnv$median_CNV_ratio, na.rm = TRUE)
for (i in seq_len(nrow(stats))) {
  x1 <- which(cnv_levels == stats$group[i])
  y <- ymax * (1.06 + 0.09 * i)
  p <- p +
    annotate(
      "segment",
      x = x1, xend = 6, y = y, yend = y, size = 0.25
    ) +
    annotate(
      "text",
      x = (x1 + 6) / 2, y = y + 0.02,
      label = stats$label[i], size = 2.2
    )
}
p <- p + coord_cartesian(ylim = c(0.9, ymax * 1.55))

save_panel(p, "Fig1D_CNV_score.pdf", width = 4.2, height = 3.4)
