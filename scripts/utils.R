suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size + 1),
      strip.text = element_text(color = "black", size = base_size),
      legend.text = element_text(color = "black", size = base_size - 1),
      legend.title = element_text(color = "black", size = base_size),
      plot.title = element_text(color = "black", size = base_size + 1, hjust = 0.5)
    )
}

save_panel <- function(plot, filename, width, height, dpi = 600) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  ggsave(
    filename = file.path(fig_dir, sub("\\.pdf$", ".png", filename)),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

add_group_metadata <- function(obj) {
  if (!"Group" %in% colnames(obj@meta.data)) {
    is_tumor <- obj$Tissue %in% c("Tumor", "T")
    is_positive <- obj$HBV %in% c("Positive", "Pos")
    obj$Group <- ifelse(
      is_tumor,
      ifelse(is_positive, "Pos_T", "Neg_T"),
      "All_N"
    )
  }
  obj$Group <- factor(obj$Group, levels = c("Pos_T", "Neg_T", "All_N"))
  obj
}

read_annotated <- function() {
  obj <- readRDS(file.path(data_dir, "03_HCC_annotated.rds"))
  add_group_metadata(obj)
}

read_spatial <- function() {
  obj <- readRDS(file.path(data_dir, "04_spatial_annotated.rds"))
  if (!"Group" %in% colnames(obj@meta.data)) {
    obj$Group <- "Pos_T"
  }
  obj
}

pct_table <- function(obj, group.by, split.by = "Group") {
  as.data.frame(table(obj@meta.data[[split.by]], obj@meta.data[[group.by]])) %>%
    setNames(c("Group", "CellType", "Freq")) %>%
    group_by(Group) %>%
    mutate(Percent = Freq / sum(Freq)) %>%
    ungroup()
}

join_layers_compat <- function(obj) {
  if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
    obj <- Seurat::JoinLayers(obj)
  }
  obj
}

build_cellchat <- function(obj, group.by, min.cells = 10, nboot = 100) {
  suppressPackageStartupMessages(library(CellChat))
  cellchat <- createCellChat(object = obj, group.by = group.by)
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    population.size = TRUE,
    nboot = nboot
  )
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat
}
