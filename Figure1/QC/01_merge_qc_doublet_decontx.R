suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(decontX)
  library(SingleCellExperiment)
  library(dplyr)
  library(Matrix)
})

data_dir <- "data"
seurat_list <- readRDS(file.path(data_dir, "00_seurat_list_raw.rds"))

hb_genes_all <- c(
  "HBA1", "HBA2", "HBB", "HBD", "HBE1",
  "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"
)

seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.MT"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  hb_genes <- intersect(hb_genes_all, rownames(obj))
  obj[["percent.HB"]] <- PercentageFeatureSet(obj, features = hb_genes)
  obj
})

seurat_list <- lapply(seurat_list, function(obj) {
  subset(
    obj,
    subset = nFeature_RNA > 200 &
      nCount_RNA > 2000 &
      percent.MT < 25 &
      percent.HB < 5
  )
})

run_doubletfinder <- function(obj) {
  obj <- NormalizeData(obj, verbose = FALSE)
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
    npcs = 30,
    verbose = FALSE
  )
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 1.5, verbose = FALSE)

  sweep_res <- paramSweep(obj, PCs = 1:30, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  homotypic_prop <- modelHomotypic(obj$seurat_clusters)
  n_exp <- round(0.075 * ncol(obj))
  n_exp_adj <- round(n_exp * (1 - homotypic_prop))

  obj <- doubletFinder(
    obj,
    PCs = 1:30,
    pN = 0.25,
    pK = pk,
    nExp = n_exp,
    sct = FALSE
  )
  obj <- doubletFinder(
    obj,
    PCs = 1:30,
    pN = 0.25,
    pK = pk,
    nExp = n_exp_adj,
    sct = FALSE,
    reuse.pANN = paste0("pANN_0.25_", pk, "_", n_exp)
  )

  obj$Doublet <- paste0(
    obj@meta.data[[ncol(obj@meta.data)]],
    "_",
    obj@meta.data[[ncol(obj@meta.data) - 1]]
  )
  obj
}

seurat_list <- lapply(seurat_list, run_doubletfinder)
seurat_list <- lapply(seurat_list, function(obj) {
  subset(obj, subset = Doublet == "Singlet_Singlet")
})

merged <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "HCC"
)

if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  merged <- JoinLayers(merged)
}

sce <- as.SingleCellExperiment(merged)
sce <- decontX(sce, batch = sce$orig.ident)

merged$contamination <- sce$decontX_contamination
decontx_assay <- grep("decont", names(assays(sce)), ignore.case = TRUE, value = TRUE)[1]
decontx_counts <- assays(sce)[[decontx_assay]]

if ("CreateAssay5Object" %in% getNamespaceExports("Seurat")) {
  merged[["decontX"]] <- CreateAssay5Object(counts = decontx_counts)
} else {
  merged[["decontX"]] <- CreateAssayObject(counts = round(decontx_counts))
}

DefaultAssay(merged) <- "decontX"
merged <- NormalizeData(
  merged,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

saveRDS(merged, file.path(data_dir, "01_HCC_doublet_decontX.rds"))
