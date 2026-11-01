source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

spatial <- read_spatial()

p_spatial <- SpatialFeaturePlot(
  spatial,
  features = "HBV",
  image.alpha = 0,
  pt.size.factor = 2
) +
  scale_fill_gradientn(
    colours = c(
      "#3B4992FF",
      "lightblue",
      "white",
      "yellow",
      "#CD202CFF",
      "darkred"
    )
  ) +
  ggtitle("HBV signal") +
  theme_pub(7)

p_vln <- VlnPlot(
  spatial,
  features = "HBV",
  group.by = "Hepa_State",
  pt.size = 0,
  cols = hepa_colors
) +
  labs(x = NULL, y = "HBV signal") +
  theme_pub(7) +
  NoLegend()

p <- p_spatial | p_vln

save_panel(p, "Fig2B_HBV_spatial_violin.pdf", width = 8, height = 3.6)
