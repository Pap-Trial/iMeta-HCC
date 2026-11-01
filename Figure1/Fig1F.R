suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

obj <- read_annotated()
hepa <- subset(obj, subset = CellType == "Hepatocyte")
Idents(hepa) <- "Hepa_State"

markers <- FindMarkers(
  hepa,
  ident.1 = "Hepa_Imm",
  ident.2 = NULL,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

genes <- rownames(markers)[which(markers$p_val_adj < 0.05)]
genes <- genes[!is.na(genes)]
if (length(genes) < 5) {
  markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]
  genes <- rownames(markers)[1:min(200, nrow(markers))]
}

ego <- enrichGO(
  gene = genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

go_df <- as.data.frame(ego)
if (nrow(go_df) == 0) {
  stop("No enriched GO terms were found for Hepa_Imm markers.")
}

go_df <- go_df[order(go_df$p.adjust), ]
go_df <- go_df[1:min(10, nrow(go_df)), ]
go_df$Description <- factor(go_df$Description, levels = rev(go_df$Description))

p_go <- ggplot(
  go_df,
  aes(x = -log10(p.adjust), y = Description, fill = -log10(p.adjust))
) +
  geom_col(width = 0.7) +
  scale_fill_gradient(low = "#FDDBC7", high = "#B2182B") +
  labs(x = "-log10(adjusted P)", y = NULL, fill = "-log10(P.adj)") +
  theme_pub(7)

sub <- subset(
  obj,
  subset = CellType %in% c("Hepatocyte", "Macrophage", "Fibroblast")
)
sub$CellType <- factor(
  sub$CellType,
  levels = c("Hepatocyte", "Macrophage", "Fibroblast")
)

p_vln <- VlnPlot(
  sub,
  features = "SPP1",
  group.by = "CellType",
  pt.size = 0,
  cols = celltype_colors[c("Hepatocyte", "Macrophage", "Fibroblast")]
) +
  labs(x = NULL, y = "SPP1 expression") +
  theme_pub(7) +
  NoLegend()

p <- p_go | p_vln

save_panel(p, "Fig1F_Hepa_Imm_GO_SPP1.pdf", width = 7.5, height = 3.4)
