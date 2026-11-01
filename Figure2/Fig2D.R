suppressPackageStartupMessages({
  library(STDistance)
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

target_label <- "Hepa_Imm"
target_type <- target_label

distance_results <- calculate_nearest_distances(
  posi,
  reference_type = reference_type,
  target_types = target_type,
  x_col = "pxl_row_in_fullres",
  y_col = "pxl_col_in_fullres",
  id_col = "Newbarcode",
  type_col = "SpatialCellType"
)

distance_results$barcode <- distance_results$Newbarcode
fib_barcodes <- posi$Newbarcode[posi$SpatialCellType == reference_type]
fib <- subset(spatial, cells = fib_barcodes)
cd44_fib <- FetchData(fib, vars = "CD44")
cd44_fib$barcode <- rownames(cd44_fib)

df_cor <- merge(
  distance_results,
  cd44_fib,
  by = "barcode",
  all.x = TRUE
)
df_cor$distance <- df_cor[[target_type]]
df_cor <- df_cor[
  !is.na(df_cor$distance) & !is.na(df_cor$CD44),
  ,
  drop = FALSE
]

cor_test <- cor.test(df_cor$distance, df_cor$CD44, method = "spearman")
cor_label <- paste0(
  "rho = ",
  round(unname(cor_test$estimate), 2),
  "\nP = ",
  format.pval(cor_test$p.value, digits = 2)
)

spatial$ViolinType <- as.character(spatial$SpatialCellType)
violin_levels <- c("Fibroblast", "Hepa_Stress", "Hepa_Norm", "Hepa_Imm")
spatial$ViolinType <- factor(spatial$ViolinType, levels = violin_levels)
sub <- subset(spatial, subset = ViolinType %in% violin_levels)

violin_colors <- c(Fibroblast = "#FFD92F", hepa_colors)
p_vln <- VlnPlot(
  sub,
  features = "CD44",
  group.by = "ViolinType",
  pt.size = 0,
  cols = violin_colors
) +
  labs(x = NULL, y = "CD44 expression") +
  theme_pub(7) +
  NoLegend()

p_cor <- ggplot(df_cor, aes(x = distance, y = CD44)) +
  geom_point(size = 0.7, alpha = 0.5, color = "#4682B4") +
  geom_smooth(method = "lm", se = TRUE, color = "#C94052", size = 0.6) +
  annotate(
    "text",
    x = min(df_cor$distance, na.rm = TRUE),
    y = max(df_cor$CD44, na.rm = TRUE),
    label = cor_label,
    hjust = 0,
    vjust = 1,
    size = 2.2
  ) +
  labs(
    x = paste0("Distance to nearest ", target_label, " (um)"),
    y = "Fibroblast CD44 expression"
  ) +
  theme_pub(7)

p <- p_vln | p_cor

save_panel(p, "Fig2D_CD44_distance.pdf", width = 7, height = 3.2)
