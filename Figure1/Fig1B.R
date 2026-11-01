source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

obj <- read_annotated()

p <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "CellType",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.05,
  raster = FALSE,
  cols = celltype_colors
) +
  ggtitle("Major cell types") +
  theme_pub(8)

save_panel(p, "Fig1B_major_celltype_umap.pdf", width = 4.2, height = 3.6)
