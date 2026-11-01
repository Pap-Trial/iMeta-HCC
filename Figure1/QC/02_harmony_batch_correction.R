suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
})

data_dir <- "data"
obj <- readRDS(file.path(data_dir, "01_HCC_doublet_decontX.rds"))
DefaultAssay(obj) <- "decontX"

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = 20,
  verbose = FALSE
)

obj$Batch <- interaction(obj$Patient, obj$Tissue, drop = TRUE)

obj <- RunHarmony(
  obj,
  group.by.vars = "Batch",
  dims.use = 1:20,
  reduction.save = "harmony",
  verbose = FALSE
)

obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "umap",
  verbose = FALSE
)
obj <- RunTSNE(
  obj,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "tsne",
  check_duplicates = FALSE
)
obj <- FindNeighbors(
  obj,
  reduction = "harmony",
  dims = 1:20,
  verbose = FALSE
)
obj <- FindClusters(obj, resolution = 1.5, verbose = FALSE)

saveRDS(obj, file.path(data_dir, "02_HCC_harmony.rds"))
