source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

obj <- read_annotated()
hepa <- subset(obj, subset = CellType == "Hepatocyte")

p_umap <- DimPlot(
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

props <- hepa@meta.data %>%
  count(Group, Hepa_State) %>%
  group_by(Group) %>%
  mutate(Percent = n / sum(n)) %>%
  ungroup()

props$Group <- factor(props$Group, levels = c("Pos_T", "Neg_T", "All_N"))
props$Hepa_State <- factor(props$Hepa_State, levels = names(hepa_colors))

p_bar <- ggplot(props, aes(x = Group, y = Percent, fill = Hepa_State)) +
  geom_col(width = 0.7, color = "white", size = 0.2) +
  scale_fill_manual(values = hepa_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion") +
  theme_pub(7)

p <- p_umap | p_bar

save_panel(p, "Fig1E_hepatocyte_umap_proportion.pdf", width = 7.2, height = 3.2)
