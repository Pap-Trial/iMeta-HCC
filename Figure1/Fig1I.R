suppressPackageStartupMessages({
  library(CellChat)
})

source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

obj <- read_annotated()

sub <- subset(
  obj,
  subset = CellType %in% c("Hepatocyte", "Macrophage", "Fibroblast") &
    !is.na(CellState)
)

group_counts <- table(sub$CellState)
keep_groups <- names(group_counts)[group_counts >= 10]
sub <- subset(sub, subset = CellState %in% keep_groups)
Idents(sub) <- "CellState"

cellchat <- createCellChat(object = sub, group.by = "CellState")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(
  cellchat,
  type = "triMean",
  population.size = TRUE,
  nboot = 100
)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file.path(data_dir, "cellchat_HepF_Mye.rds"))

fib <- subset(obj, subset = CellType == "Fibroblast")
p_umap <- DimPlot(
  fib,
  reduction = "umap",
  group.by = "Fib_State",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.05,
  raster = FALSE,
  cols = fib_colors
) +
  ggtitle("Fibroblast states") +
  theme_pub(7)

save_panel(p_umap, "Fig1I_fibroblast_umap.pdf", width = 3.6, height = 3.2)

pdf(
  file.path(fig_dir, "Fig1I_SPP1_cellchat_network.pdf"),
  width = 5,
  height = 5
)
netVisual_aggregate(
  cellchat,
  signaling = "SPP1",
  layout = "circle",
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE
)
dev.off()

png(
  file.path(fig_dir, "Fig1I_SPP1_cellchat_network.png"),
  width = 5,
  height = 5,
  units = "in",
  res = 600
)
netVisual_aggregate(
  cellchat,
  signaling = "SPP1",
  layout = "circle",
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE
)
dev.off()
