suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(GenomicRanges)
})

source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

obj <- read_annotated()
hepa <- subset(obj, subset = CellType == "Hepatocyte")

if (!"BCL3_regulon" %in% colnames(hepa@meta.data)) {
  scenic_db_dir <- file.path(data_dir, "cisTarget_databases")

  counts <- GetAssayData(hepa, assay = "RNA", slot = "counts")
  scenicOptions <- initializeScenic(
    org = "hgnc",
    dbDir = scenic_db_dir,
    nCores = 4
  )

  genes_kept <- geneFiltering(
    counts,
    scenicOptions = scenicOptions,
    minCountsPerGene = 3 * 0.01 * ncol(counts),
    minSamples = 0.01 * ncol(counts)
  )
  expr_mat <- counts[genes_kept, , drop = FALSE]
  expr_mat <- log2(expr_mat + 1)

  runGenie3(expr_mat, scenicOptions, nParts = 10)
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(
    scenicOptions,
    coexMethod = "top5perTarget"
  )
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expr_mat)
  scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

  regulon_auc <- loadInt(scenicOptions, "aucell_regulonAUC")
  bcl3_row <- grep("BCL3", rownames(regulon_auc))[1]
  if (is.na(bcl3_row)) {
    stop("BCL3 regulon was not found in the SCENIC output.")
  }
  hepa$BCL3_regulon <- as.numeric(regulon_auc[bcl3_row, colnames(hepa)])


}

p_state <- DimPlot(
  hepa,
  reduction = "umap",
  group.by = "Hepa_State",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.05,
  raster = FALSE,
  cols = hepa_colors
) +
  ggtitle("Hepatocyte states") +
  theme_pub(7)

p_regulon <- FeaturePlot(
  hepa,
  features = "BCL3_regulon",
  reduction = "umap",
  pt.size = 0.05,
  raster = FALSE,
  cols = c("lightgrey", "#D62728")
) +
  ggtitle("BCL3 regulon activity") +
  theme_pub(7)

p <- p_state | p_regulon

save_panel(p, "Fig1H_Hepa_BCL3_regulon.pdf", width = 7.5, height = 3.4)
