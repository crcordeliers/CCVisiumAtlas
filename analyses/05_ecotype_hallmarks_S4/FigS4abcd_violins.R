#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S4a-d - molecular hallmarks of the meta-clusters (ssGSEA scores of the 3,231 meta-spots)
#   S4a  CMS2, CMS3 (CMScaller templates), Pelka epithelial, Guinney epithelial: all meta-clusters;
#        cell cycle, WNT (Guinney) and CNV fraction (fastCNV): the 4 tumour meta-clusters, with Wilcoxon P of
#        tumour-cycling vs each other tumour meta-cluster
#   S4b  CMS1, CMS4, Pelka stromal, Guinney mesenchymal: all meta-clusters
#   S4c  neutrophils (MCP-counter); S4d  T cells and B lineage (MCP-counter)
# Input : Results_Harmony/res_10/seurat_obj_with_pelka.RData (R/Step6; also carries the scores computed in R/Step3)
# Output: Figures/FigS4/Violin_<score>.pdf (one file per score)
# Usage : Rscript analyses/05_ecotype_hallmarks_S4/FigS4abcd_violins.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })
set.seed(1)   # fixed jitter of the points (no effect on the statistics)

load(file.path(harmony.dir, "seurat_obj_with_pelka.RData"))       # -> seurat_obj
seurat_obj$merged_meta_clusters <- factor(seurat_obj$merged_meta_clusters, levels = meta_order)
out_dir <- file.path(figures.dir, "FigS4"); dir.create(out_dir, showWarnings = FALSE)

# score column -> panel title
ALL <- c(CMS2 = "CMS2", CMS3 = "CMS3", Epithelial = "Pelka Epithelial", Epithelial_sig = "Guinney Epithelial",
         CMS1 = "CMS1", CMS4 = "CMS4", Stromal = "Pelka Stromal", Mesenchymal_sig = "Guinney Mesenchymal",
         Neutrophils = "Neutrophils", T_cells = "T cells", B_lineage = "B cells")
TUM <- c(CellCycle = "Cell cycle", WNT = "WNT", CNV_fraction = "CNV fraction")

violin <- function(df, title, ylab) {
  ggplot(df, aes(cluster, value, fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.3, height = 0, size = 1.5, alpha = 0.3) +
    scale_fill_manual(values = meta_colors) +
    labs(y = ylab, title = title) + theme_bw() +
    theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_blank(),
          axis.text = element_text(face = "bold", size = rel(1.6)), axis.title.y = element_text(face = "bold", size = rel(1.8)),
          plot.title = element_text(face = "bold", size = rel(1.6)))
}

for (v in names(ALL)) {
  df <- data.frame(cluster = seurat_obj$merged_meta_clusters, value = seurat_obj@meta.data[[v]])
  ggsave(file.path(out_dir, paste0("Violin_", v, "_metaclusters.pdf")), violin(df, ALL[[v]], "Enrichment score"), width = 6, height = 4)
}

# tumour meta-clusters, with P values (tumour-cycling vs the 3 others)
tum <- seurat_obj$merged_meta_clusters %in% c("0", "1", "7", "8")
for (v in names(TUM)) {
  df <- data.frame(cluster = droplevels(seurat_obj$merged_meta_clusters[tum]), value = seurat_obj@meta.data[[v]][tum])
  ymax <- max(df$value, na.rm = TRUE); ymin <- min(df$value, na.rm = TRUE)
  comps <- list(c("0", "1"), c("0", "7"), c("0", "8"))
  pv <- bind_rows(lapply(comps, function(cmp) {
    p <- wilcox.test(df$value[df$cluster == cmp[1]], df$value[df$cluster == cmp[2]])$p.value
    data.frame(group1 = cmp[1], group2 = cmp[2], p = p,
               p_label = if (p < 1e-4) "P < 0.0001" else if (p < 1e-3) "P < 0.001" else paste0("P = ", signif(p, 3)))
  }))
  pv$y  <- c(ymax * 1.02, ymax * 1.17, ymax * 1.32)
  lv <- levels(df$cluster)
  pv$x1 <- match(pv$group1, lv); pv$x2 <- match(pv$group2, lv); pv$x <- (pv$x1 + pv$x2) / 2
  p <- violin(df, TUM[[v]], if (v == "CNV_fraction") "CNV fraction" else "Enrichment score") +
    scale_y_continuous(limits = c(ymin * 0.2, ymax * 1.4), expand = expansion(mult = c(0.05, 0.10))) +
    geom_segment(data = pv, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.7) +
    geom_text(data = pv, aes(x = x, y = y * 1.05, label = p_label), size = 5, fontface = "bold", inherit.aes = FALSE)
  ggsave(file.path(out_dir, paste0("Violin_", v, "_tumor_metaclusters.pdf")), p, width = 4, height = 4)
  write.csv(pv, file.path(out_dir, paste0("Violin_", v, "_tumor_metaclusters_P.csv")), row.names = FALSE)
}
