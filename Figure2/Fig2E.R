suppressPackageStartupMessages({
  library(readxl)
  library(pheatmap)
  library(dplyr)
})

data_dir <- ""     # <<< 如路径不同请修改
setwd(data_dir)

source("01_common/config.R")

tma_file <- "Figure2E.xlsx"

tma <- read_excel(tma_file, sheet = 1)
tma$Tissue <- dplyr::recode(
  as.character(tma$Tissue),
  Normal = "Para",
  N = "Para",
  T = "Tumor"
)

vars <- c(
  "Meandistance_TumortoCOL1A1",
  "Percent_SPP1",
  "Percent_COL1A1",
  "Percent_GPC3"
)
if (!all(vars %in% colnames(tma))) {
  stop("Figure2E.xlsx is missing one or more required columns.")
}

tumor <- as.data.frame(tma[tma$Tissue == "Tumor", vars, drop = FALSE])
para <- as.data.frame(tma[tma$Tissue == "Para", vars, drop = FALSE])

if (nrow(tumor) != nrow(para)) {
  stop("Tumor and Para rows are not paired.")
}

tumor_mat <- t(as.matrix(tumor))
para_mat <- t(as.matrix(para))
colnames(tumor_mat) <- paste0("S", seq_len(ncol(tumor_mat)))
colnames(para_mat) <- paste0("S", seq_len(ncol(para_mat)))

dat <- rbind(para_mat, tumor_mat)
rownames(dat) <- c(
  "Distance_Para",
  "SPP1+_Para",
  "COL1A1+_Para",
  "GPC3+_Para",
  "Distance_Tumor",
  "SPP1+_Tumor",
  "COL1A1+_Tumor",
  "GPC3+_Tumor"
)

dat_scaled <- scale(dat)
dat_scaled[is.na(dat_scaled)] <- 0

fig_dir <- "outputs"
dir.create(fig_dir, showWarnings = FALSE)

pheatmap(
  dat_scaled,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("#A6CEE3", "white", "#F28E8B"))(50),
  fontsize_row = 7,
  fontsize_col = 5,
  border_color = "white",
  main = "TMA marker positivity and spatial distance",
  filename = file.path(fig_dir, "Fig2E_TMA_heatmap.pdf"),
  width = 8,
  height = 4
)

pheatmap(
  dat_scaled,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("#A6CEE3", "white", "#F28E8B"))(50),
  fontsize_row = 7,
  fontsize_col = 5,
  border_color = "white",
  main = "TMA marker positivity and spatial distance",
  filename = file.path(fig_dir, "Fig2E_TMA_heatmap.png"),
  width = 8,
  height = 4
)
