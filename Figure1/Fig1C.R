source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

obj <- read_annotated()
meta <- obj@meta.data
meta$Sample <- paste(meta$Patient, meta$Tissue, sep = "_")

sample_order <- meta %>%
  distinct(Sample, Patient, Tissue, HBV) %>%
  arrange(Patient, Tissue)

props <- meta %>%
  count(Sample, CellType) %>%
  group_by(Sample) %>%
  mutate(Percent = n / sum(n)) %>%
  ungroup()

props$Sample <- factor(props$Sample, levels = sample_order$Sample)
props$CellType <- factor(props$CellType, levels = names(celltype_colors))

p_bar <- ggplot(props, aes(x = Sample, y = Percent, fill = CellType)) +
  geom_col(width = 0.9, color = NA) +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Cell proportion") +
  theme_pub(6) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

annotation <- sample_order
annotation$Sample <- factor(annotation$Sample, levels = sample_order$Sample)

p_hbv <- ggplot(annotation, aes(x = Sample, y = 1, fill = HBV)) +
  geom_tile() +
  scale_fill_manual(values = hbv_colors, drop = FALSE) +
  theme_void() +
  theme(legend.position = "bottom")

p_tissue <- ggplot(annotation, aes(x = Sample, y = 1, fill = Tissue)) +
  geom_tile() +
  scale_fill_manual(values = tissue_colors, drop = FALSE) +
  theme_void() +
  theme(legend.position = "bottom")

p <- p_bar / p_hbv / p_tissue + plot_layout(heights = c(10, 0.5, 0.5))

save_panel(p, "Fig1C_cell_composition.pdf", width = 8, height = 5)
