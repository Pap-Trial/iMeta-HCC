suppressPackageStartupMessages({
  library(STDistance)
  library(ape)
  library(ggtree)
})

source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

spatial <- read_spatial()

coords <- GetTissueCoordinates(spatial)
coords$Newbarcode <- coords$barcode
coords$Sample <- "spatial"
coords$Sampleid <- "spatial"

coords_norm <- normalize_spatial(coords)

meta <- spatial@meta.data
meta$Newbarcode <- rownames(meta)
posi <- merge(coords_norm, meta, by = "Newbarcode", all.x = TRUE)
posi <- posi[!is.na(posi$SpatialCellType), , drop = FALSE]

reference_type <- "Fibroblast"
target_types <- c("Hepa_Stress", "Hepa_Norm", "Hepa_Imm")
target_types <- target_types[target_types %in% posi$SpatialCellType]

distance_results <- calculate_nearest_distances(
  posi,
  reference_type = reference_type,
  target_types = target_types,
  x_col = "pxl_row_in_fullres",
  y_col = "pxl_col_in_fullres",
  id_col = "Newbarcode",
  type_col = "SpatialCellType"
)

p_box <- plot_distance_boxplot(
  distance_results,
  id_col = "Newbarcode",
  show_points = FALSE,
  y_scale = "log10",
  palette = "Dark2"
) +
  labs(x = NULL, y = "Distance (log10)") +
  theme_pub(7)

states <- c(target_types, reference_type)
distance_matrix <- matrix(
  NA_real_,
  nrow = length(states),
  ncol = length(states),
  dimnames = list(states, states)
)
diag(distance_matrix) <- 0

for (i in seq_along(states)) {
  dr <- calculate_nearest_distances(
    posi,
    reference_type = states[i],
    target_types = states[-i],
    x_col = "pxl_row_in_fullres",
    y_col = "pxl_col_in_fullres",
    id_col = "Newbarcode",
    type_col = "SpatialCellType"
  )
  for (j in setdiff(seq_along(states), i)) {
    distance_matrix[i, j] <- mean(dr[[states[j]]], na.rm = TRUE)
  }
}

distance_matrix <- (distance_matrix + t(distance_matrix)) / 2
tree <- ape::nj(as.dist(distance_matrix))

p_tree <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(size = 1.8, align = TRUE) +
  ggtitle("NJ tree") +
  theme_pub(7)

sub <- subset(spatial, subset = CellType %in% c("Hepatocyte", "Fibroblast"))
sub$State <- ifelse(
  sub$CellType == "Hepatocyte",
  as.character(sub$Hepa_State),
  as.character(sub$Fib_State)
)

p_umap <- DimPlot(
  sub,
  reduction = "umap",
  group.by = "State",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.05,
  raster = FALSE,
  cols = c(hepa_colors, fib_colors)
) +
  ggtitle("Hepatocyte and fibroblast states") +
  theme_pub(7)

p <- p_umap / p_tree / p_box + plot_layout(heights = c(1.1, 1, 1))

save_panel(p, "Fig2C_spatial_proximity.pdf", width = 4.5, height = 7)
