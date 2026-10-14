#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S5c (step 2) - meta-cluster proportions according to the CMS of the sample (n = 42 classified)
# Kruskal-Wallis test across CMS1-4, Dunn post-hoc tests with Holm correction (rstatix).
# Inputs : Tables/CMS_pseudobulk_48.csv (01_*.R); Tables/Proportions_metaclusters_per_sample.csv
#          (analyses/06_ecotype_composition_S5ab_S6_S7)
# Outputs: Figures/FigS5c/boxplot_cms_<meta-cluster>.pdf ; Tables/FigS5c_kruskal.csv, FigS5c_dunn.csv
# Usage  : Rscript analyses/07_CMS_pseudobulk_S5c/02_FigS5c_ecotypes_by_CMS.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(ggtext); library(rstatix) })
set.seed(1)   # fixed jitter of the points (no effect on the statistics)

cms  <- read.csv(file.path(tables.dir, "CMS_pseudobulk_48.csv"))
prop <- read.csv(file.path(tables.dir, "Proportions_metaclusters_per_sample.csv"), colClasses = c(cl = "character"))
d0 <- inner_join(prop, cms[, c("sample", "CMS_global")], by = c("Sample" = "sample")) |>
  filter(CMS_global %in% c("CMS1", "CMS2", "CMS3", "CMS4")) |>
  mutate(cms_group = factor(CMS_global, levels = c("CMS1", "CMS2", "CMS3", "CMS4")), Proportion = Proportion / 100)
out_dir <- file.path(figures.dir, "FigS5c"); dir.create(out_dir, showWarnings = FALSE)
PAL <- c(CMS1 = "#DE8D25", CMS2 = "#1E5BA0", CMS3 = "#BC6095", CMS4 = "#3E8E61")

kw_all <- list(); dunn_all <- list()
for (k in meta_order) {
  d <- filter(d0, cl == k)
  kw <- kruskal_test(d, Proportion ~ cms_group)
  kw_all[[k]] <- mutate(kw, metacluster = meta_names[[k]])
  dunn_all[[k]] <- mutate(dunn_test(d, Proportion ~ cms_group, p.adjust.method = "holm"), metacluster = meta_names[[k]])
  n <- count(d, cms_group)
  lab <- if (kw$p < 0.001) "italic(p) < 0.001" else paste0("italic(p) == ", signif(kw$p, 2))
  p <- ggplot(d, aes(cms_group, Proportion, fill = cms_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) + geom_jitter(width = 0.2, size = 2, alpha = 0.3) +
    scale_fill_manual(values = PAL) + scale_x_discrete(labels = setNames(paste0(n$cms_group, "<br>(n=", n$n, ")"), n$cms_group)) +
    annotate("text", x = 2.5, y = max(d$Proportion) * 1.1, label = lab, parse = TRUE, size = 6) +
    labs(title = meta_names[[k]], y = "Proportion") + theme_bw() +
    theme(legend.position = "none", axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_markdown(size = rel(1.8)), axis.text.y = element_text(face = "bold", size = rel(1.6)),
          axis.title.y = element_text(face = "bold", size = rel(1.8)), plot.title = element_text(face = "bold", size = rel(1.6), hjust = 0.5))
  ggsave(file.path(out_dir, paste0("boxplot_cms_", k, ".pdf")), p, width = 5, height = 5.5)
}
write.csv(bind_rows(kw_all), file.path(tables.dir, "FigS5c_kruskal.csv"), row.names = FALSE)
write.csv(bind_rows(dunn_all), file.path(tables.dir, "FigS5c_dunn.csv"), row.names = FALSE)
