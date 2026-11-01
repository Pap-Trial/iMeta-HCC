suppressPackageStartupMessages({
  library(CellChat)
})

source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

obj <- read_annotated()

sub_all <- subset(
  obj,
  subset = CellType %in% c("Hepatocyte", "Macrophage", "Fibroblast") &
    !is.na(CellState)
)
keep_groups <- names(table(sub_all$CellState))[table(sub_all$CellState) >= 10]
sub_all <- subset(sub_all, subset = CellState %in% keep_groups)

states <- unique(sub_all$CellState)
lineage_map <- setNames(
  ifelse(
    grepl("^Hepa", states),
    "Hepatocyte",
    ifelse(
      grepl("^Fib", states),
      "Fibroblast",
      ifelse(grepl("^Macro|^Mono", states), "Macrophage", "Other")
    )
  ),
  states
)

run_cellchat <- function(obj) {
  set.seed(123)
  cellchat <- createCellChat(object = obj, group.by = "CellState")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    population.size = TRUE,
    nboot = 100
  )
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat
}

downsample_states <- function(obj, n = 400, seed = 123) {
  set.seed(seed)
  cells <- unlist(lapply(split(colnames(obj), obj$CellState), function(x) {
    if (length(x) > n) sample(x, n) else x
  }), use.names = FALSE)
  subset(obj, cells = cells)
}

compute_spp1 <- function(obj, design) {
  cellchat <- run_cellchat(obj)
  prob <- cellchat@netP$prob
  if (!"SPP1" %in% dimnames(prob)[[3]]) {
    stop("SPP1 pathway was not found in the CellChat result.")
  }
  prob_mat <- prob[, , "SPP1"]
  target_states <- colnames(prob_mat)
  target_fib <- target_states[lineage_map[target_states] == "Fibroblast"]
  if (length(target_fib) == 0) {
    stop("No fibroblast target states were found.")
  }

  prob_fib <- prob_mat[, target_fib, drop = FALSE]
  denom <- colSums(prob_fib, na.rm = TRUE)
  denom[denom == 0] <- 1
  contrib <- sweep(prob_fib, 2, denom, "/")

  df <- as.data.frame.table(contrib, stringsAsFactors = FALSE)
  names(df) <- c("SourceState", "TargetState", "Contribution")
  df$Source <- lineage_map[df$SourceState]
  df$Source[is.na(df$Source)] <- "Other"
  df$Design <- design
  df
}

sub_pos <- subset(sub_all, subset = Group == "Pos_T")
keep_pos <- names(table(sub_pos$CellState))[table(sub_pos$CellState) >= 10]
sub_pos <- subset(sub_pos, subset = CellState %in% keep_pos)
sub_eq <- downsample_states(sub_all, n = 400, seed = 123)

designs <- list(
  "All cells" = sub_all,
  "Pos_T" = sub_pos,
  "Equal-n (400 cells)" = sub_eq
)

df <- dplyr::bind_rows(lapply(names(designs), function(design) {
  compute_spp1(designs[[design]], design)
}))

df$Design <- factor(df$Design, levels = names(designs))
df$SourceState <- factor(df$SourceState, levels = names(lineage_map))

source_colors <- c(
  Hepatocyte = "#6A3D9A",
  Macrophage = "#E31A1C",
  Fibroblast = "#FFD92F",
  Other = "grey70"
)
df$Source <- factor(df$Source, levels = names(source_colors))

df_lineage <- df %>%
  group_by(Design, Source) %>%
  summarise(Contribution = mean(Contribution, na.rm = TRUE), .groups = "drop")

p_bar <- ggplot(df_lineage, aes(x = Design, y = Contribution, fill = Source)) +
  geom_col(width = 0.7, position = "fill", color = "white", size = 0.2) +
  scale_fill_manual(values = source_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "SPP1 signaling contribution") +
  theme_pub(7)

p_box <- df %>%
  filter(Design == "Pos_T") %>%
  ggplot(aes(x = SourceState, y = Contribution, fill = Source)) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    color = "black",
    size = 0.3,
    alpha = 0.7
  ) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.5) +
  scale_fill_manual(values = source_colors, drop = FALSE) +
  facet_grid(~Source, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Contribution in Pos_T") +
  theme_pub(7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

p <- p_bar / p_box + plot_layout(heights = c(1, 1.2))

save_panel(p, "Fig1J_SPP1_contribution.pdf", width = 8, height = 5)
