suppressPackageStartupMessages({
  library(readxl)
  library(survival)
  library(survminer)
  library(pROC)
})

data_dir <- ""     # <<< 如路径不同请修改
setwd(data_dir)

tma_file <- 'Figure2G.xlsx'

source("01_common/config.R")
source("01_common/palettes.R")
source("01_common/utils.R")

tma <- readxl::read_excel(tma_file)

d <- tma %>%
  filter(HBV == "Positive", Tissue == "Tumor") %>%
  mutate(
    OS = as.numeric(as.character(OS)),
    OS.Month = as.numeric(OS.Month),
    dist = as.numeric(Meandistance_TumortoCOL1A1),
    spp1 = as.numeric(Percent_SPP1),
    gpc3 = as.numeric(Percent_GPC3),
    col1a1 = as.numeric(Percent_COL1A1)
  ) %>%
  filter(!is.na(dist), !is.na(OS), !is.na(OS.Month))

scatter_panel <- function(df, x, y, xlab, ylab) {
  df <- df[!is.na(df[[x]]) & !is.na(df[[y]]), , drop = FALSE]
  cor_test <- cor.test(df[[x]], df[[y]], method = "pearson")
  p_label <- if (cor_test$p.value < 0.001) {
    "P < 0.001"
  } else {
    paste0("P = ", round(cor_test$p.value, 3))
  }
  cor_label <- paste0(
    "r = ",
    round(unname(cor_test$estimate), 3),
    "\n",
    p_label
  )
  print(round(unname(cor_test$estimate), 3))
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(size = 0.9, alpha = 0.7, color = "#4682B4") +
    geom_smooth(method = "lm", se = TRUE, color = "#C94052", size = 0.6) +
    annotate(
      "text",
      x = min(df[[x]], na.rm = TRUE),
      y = max(df[[y]], na.rm = TRUE),
      label = cor_label,
      hjust = 0,
      vjust = 1,
      size = 2.2
    ) +
    labs(x = xlab, y = ylab) +
    theme_pub(7)
}

p_g1 <- scatter_panel(
  d,
  "spp1",
  "gpc3",
  "Percentage of SPP1",
  "Percentage of GPC3"
)
p_g2 <- scatter_panel(
  d,
  "spp1",
  "col1a1",
  "Percentage of SPP1",
  "Percentage of COL1A1"
)

dist_cutoff <- median(d$dist, na.rm = TRUE)
d$group <- factor(
  ifelse(d$dist <= dist_cutoff, "High", "Low"),
  levels = c("High", "Low")
)

fit <- survfit(Surv(OS.Month, OS) ~ group, data = d)
p_km <- ggsurvplot(
  fit,
  data = d,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = FALSE,
  palette = c("#C94052", "#4F7A5E"),
  size = 0.8,
  xlim = c(0, 60),
  break.time.by = 10,
  legend.title = "",
  legend.labs = c("High", "Low"),
  surv.median.line = "hv",
  xlab = "Time (months)",
  ylab = "Survival probability",
  ggtheme = theme_pub(7)
)$plot

roc_obj <- roc(d$OS, -d$dist, quiet = TRUE)
auc_value <- as.numeric(auc(roc_obj))
p_roc <- ggroc(roc_obj, legacy.axes = TRUE, size = 0.8) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  annotate(
    "text",
    x = 0.55,
    y = 0.15,
    label = paste0("AUC = ", round(auc_value, 3)),
    size = 2.5
  ) +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_pub(7)

p <- (p_g1 | p_g2) / (p_km | p_roc) + plot_layout(heights = c(1, 1.2))

save_panel(p, "Fig2G_TMA_prognosis.pdf", width = 7.2, height = 6)
