#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S2a,b - spatial organisation of the meta-clusters
#   S2a  projection of the meta-clusters on the 48 Visium capture areas (one PDF per sample)
#   S2b  Moran's I of each meta-cluster in each sample (binary indicator "spot in cluster", k = 6 nearest spots,
#        row-standardised weights, spdep::moran.mc with 999 permutations, one-sided); boxplots per meta-cluster in the
#        display order of Fig. 2c, dashed line = median
# Inputs : Processed_seurat/processed_seurat_<Sample>_withmeta.RData (R/Step5)
# Outputs: Figures/FigS2a_projection/<Sample>.pdf, Figures/FigS2b_MoranI_by_metacluster.pdf, Tables/FigS2b_MoranI.csv
# Usage  : Rscript analyses/02_spatial_organisation_S2ab/FigS2ab_projection_and_MoranI.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(purrr); library(spdep); library(ggplot2) })

dir.create(file.path(figures.dir, "FigS2a_projection"), showWarnings = FALSE)

get_xy <- function(seu) {
  xy <- try(GetTissueCoordinates(seu, filter = "tissue"), silent = TRUE)
  if (!inherits(xy, "try-error")) return(data.frame(barcode = rownames(xy), x = xy[[1]], y = xy[[2]]))
  img <- seu@images[[1]]@coordinates
  data.frame(barcode = rownames(img), x = img$pxl_col_in_fullres, y = img$pxl_row_in_fullres)
}

moran_one <- function(seu, cluster_var, cluster_id, k = 6, nsim = 999) {
  coords <- get_xy(seu)
  keep <- intersect(rownames(seu@meta.data), coords$barcode)
  coords <- coords[match(keep, coords$barcode), ]
  x <- as.numeric(seu@meta.data[keep, cluster_var] == cluster_id)
  if (var(x) == 0) return(data.frame(cluster = cluster_id, MoranI = NA, p_value = NA, n_spots = sum(x)))
  lw <- nb2listw(knn2nb(knearneigh(as.matrix(coords[, c("x", "y")]), k = k)), style = "W", zero.policy = TRUE)
  t <- moran.mc(x, listw = lw, nsim = nsim, zero.policy = TRUE, alternative = "greater")
  data.frame(cluster = cluster_id, MoranI = unname(t$statistic), p_value = t$p.value, n_spots = sum(x))
}

set.seed(1)
res <- map_dfr(Sample.names, function(s) {
  load(file.path(processed.dir, paste0("processed_seurat_", s, "_withmeta.RData")))   # -> sobj
  p <- SpatialDimPlot(sobj, group.by = "merged_meta_clusters", pt.size.factor = 2.5, alpha = 1) +
    scale_fill_manual(values = meta_colors) + ggtitle(s) + theme(legend.position = "bottom")
  ggsave(file.path(figures.dir, "FigS2a_projection", paste0(s, ".pdf")), p, width = 6, height = 6)
  cl <- sort(unique(na.omit(as.character(sobj$merged_meta_clusters))))
  map_dfr(cl, ~ moran_one(sobj, "merged_meta_clusters", .x)) |> mutate(sample = s)
})
write.csv(res, file.path(tables.dir, "FigS2b_MoranI.csv"), row.names = FALSE)
message("median Moran's I ", round(median(res$MoranI, na.rm = TRUE), 2), "; P < 0.05 in ",
        sum(res$p_value < 0.05, na.rm = TRUE), "/", sum(!is.na(res$p_value)), " pairs")

res$cluster <- factor(res$cluster, levels = meta_order)
p <- ggplot(res, aes(cluster, MoranI, fill = cluster)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = median(res$MoranI, na.rm = TRUE), linetype = 2) +
  scale_fill_manual(values = meta_colors, guide = "none") +
  labs(x = "Metacluster IDs", y = "Moran's I", title = "Moran's I distribution by metacluster") + theme_classic()
ggsave(file.path(figures.dir, "FigS2b_MoranI_by_metacluster.pdf"), p, width = 6, height = 4)
