source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

spatial <- read_spatial()

spatial_colors <- c(hepa_colors, fib_colors, mye_colors, celltype_colors)
spatial_colors <- spatial_colors[!duplicated(names(spatial_colors))]

p_he <- SpatialDimPlot(
  spatial,
  image.alpha = 1,
  pt.size.factor = 0.01,
  alpha = 0,
  crop = FALSE
) +
  ggtitle("H&E")

p_state <- SpatialDimPlot(
  spatial,
  group.by = "SpatialCellType",
  pt.size.factor = 2,
  alpha = 0.8,
  crop = FALSE,
  image.alpha = 0,
  cols = spatial_colors
) +
  ggtitle("Spatial cell states") +
  theme_pub(7)

p <- p_he | p_state

save_panel(p, "Fig2A_spatial_cellstates.pdf", width = 8, height = 4)
