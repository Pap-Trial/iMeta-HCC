suppressPackageStartupMessages({
  library(ggrepel)
})

source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

obj <- read_annotated()
hepa <- subset(obj, subset = CellType == "Hepatocyte")
Idents(hepa) <- "Group"

run_deg <- function(ident.1) {
  markers <- FindMarkers(
    hepa,
    ident.1 = ident.1,
    ident.2 = "All_N",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  markers$gene <- rownames(markers)
  markers$comparison <- ident.1
  markers
}

deg <- rbind(run_deg("Pos_T"), run_deg("Neg_T"))

deg$threshold <- "ns"
deg$threshold[deg$avg_log2FC > 1 & deg$p_val_adj < 0.05] <- "up"
deg$threshold[deg$avg_log2FC < -1 & deg$p_val_adj < 0.05] <- "down"
deg$threshold <- factor(deg$threshold, levels = c("down", "ns", "up"))

label_df <- deg %>%
  filter(threshold != "ns") %>%
  group_by(comparison) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

comparison_labels <- c(
  Pos_T = "HBV+ tumor vs adjacent",
  Neg_T = "HBV- tumor vs adjacent"
)

p <- ggplot(
  deg,
  aes(x = avg_log2FC, y = -log10(p_val_adj), color = threshold)
) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 1.8,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(down = "#377EB8", ns = "grey80", up = "#E41A1C")
  ) +
  facet_wrap(
    ~comparison,
    nrow = 1,
    labeller = labeller(comparison = comparison_labels)
  ) +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[10]~adjusted~P)
  ) +
  theme_pub(7) +
  theme(legend.position = "top")

save_panel(p, "Fig1G_hepatocyte_DEG_volcano.pdf", width = 7, height = 3.4)
