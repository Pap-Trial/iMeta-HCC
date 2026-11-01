source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

obj <- read_annotated()

df <- obj@meta.data %>%
  filter(CellType == "Hepatocyte") %>%
  mutate(
    CNV_class = ifelse(Group == "All_N", "Normal", as.character(Hepa_State))
  )

cnv_levels <- c(
  "Normal",
  "Hepa_Norm",
  "Hepa_Stress",
  "Hepa_Invasive",
  "Hepa_Meta",
  "Hepa_Imm"
)
df$CNV_class <- factor(df$CNV_class, levels = cnv_levels)

fill_values <- c(Normal = "grey70", hepa_colors)

p <- ggplot(df, aes(x = CNV_class, y = CNV_score, fill = CNV_class)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.7) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    size = 0.3
  ) +
  scale_fill_manual(values = fill_values, guide = "none") +
  labs(x = NULL, y = "CNV score") +
  theme_pub(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_panel(p, "Fig1D_CNV_score.pdf", width = 4.2, height = 3.4)
