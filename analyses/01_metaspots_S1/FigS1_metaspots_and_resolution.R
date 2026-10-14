#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S1 - construction of the integrated meta-spot dataset
#   S1a  number of spots per meta-spot, per sample (boxplots)                      [re-written, see note]
#   S1b  distribution of the number of spots per meta-spot (dashed line = median)  [re-written, see note]
#   S1c  number of meta-clusters as a function of the Seurat resolution (0-2, step 0.01; dashed line = 0.5)
#   S1d  clustree of the meta-clusters for resolutions 0.2-1.2
#   S1e  alluvial plot: patient-specific meta-spots (by sample) -> cohort-wide meta-clusters
# Inputs : Results_Harmony/res_10/seurat_object_res_10.rds and scaled_data_list_10.RData (R/Step3, R/Step5)
# Outputs: Figures/FigS1a..e_*.pdf, Tables/FigS1c_resolution_sweep.csv
# Usage  : Rscript analyses/01_metaspots_S1/FigS1_metaspots_and_resolution.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(clustree); library(ggalluvial) })

seurat_obj <- readRDS(file.path(harmony.dir, "seurat_object_res_10.rds"))
load(file.path(harmony.dir, "scaled_data_list_10.RData"))     # per-sample spot tables with Cluster_res_10

## ---- S1a, S1b: spots per meta-spot ----
sp <- bind_rows(lapply(names(scaled_data_list), function(s)
  data.frame(sample = s, cluster = scaled_data_list[[s]]$Cluster_res_10))) |>
  count(sample, cluster, name = "n_spots")
stopifnot(nrow(sp) == ncol(seurat_obj))
message("meta-spots per sample: median ", median(table(sp$sample)), " (", paste(range(table(sp$sample)), collapse = "-"),
        "); spots per meta-spot: median ", median(sp$n_spots))
p <- ggplot(sp, aes(sample, n_spots)) + geom_boxplot(fill = "grey70", outlier.size = 0.5) +
  labs(x = "Sample", y = "Number of spots per meta-spot", title = "Distribution of the number of spots per metaspot across samples") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
ggsave(file.path(figures.dir, "FigS1a_spots_per_metaspot_by_sample.pdf"), p, width = 12, height = 4)
med <- median(sp$n_spots)
p <- ggplot(sp, aes(n_spots)) + geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#9AAAB5", colour = "white") +
  geom_density(colour = "blue") + geom_vline(xintercept = med, colour = "blue", linetype = 2) +
  annotate("text", x = med, y = 0, label = paste0("Median = ", med), vjust = -0.5, hjust = -0.1, colour = "blue") +
  labs(x = "Number of spots per metaspot", y = "Density") + theme_classic()
ggsave(file.path(figures.dir, "FigS1b_spots_per_metaspot_histogram.pdf"), p, width = 5, height = 4)

## ---- S1c: number of meta-clusters vs resolution (Harmony SNN graph of the meta-spots) ----
graph_name <- "RNA_snn"
df_knee <- bind_rows(lapply(seq(0, 2, by = 0.01), function(r) {
  tmp <- FindClusters(seurat_obj, resolution = r, graph.name = graph_name, verbose = FALSE)
  data.frame(resolution = r, n_clusters = nlevels(Idents(tmp)))
}))
write.csv(df_knee, file.path(tables.dir, "FigS1c_resolution_sweep.csv"), row.names = FALSE)
p <- ggplot(df_knee, aes(resolution, n_clusters)) + geom_point(size = 1.9, alpha = 0.2) +
  geom_smooth(method = "loess", span = 0.5, se = FALSE, linewidth = 1.4) +
  geom_vline(xintercept = res_meta, linetype = 2) +
  labs(x = "Seurat resolution", y = "Number of meta-clusters", title = "Knee plot (metaspots, Harmony)") +
  theme_classic(base_size = 14)
ggsave(file.path(figures.dir, "FigS1c_knee_metaclusters.pdf"), p, width = 6.5, height = 4.2)

## ---- S1d: clustree ----
obj <- seurat_obj
for (r in seq(0.2, 1.2, by = 0.1))
  if (!paste0(graph_name, "_res.", r) %in% colnames(obj@meta.data))
    obj <- FindClusters(obj, resolution = r, graph.name = graph_name, verbose = FALSE)
p <- clustree(obj, prefix = paste0(graph_name, "_res."), prop_filter = 0.02, node_size = "size", edge_arrow = TRUE) +
  theme_classic(base_size = 13) + labs(title = "Cluster transitions across resolutions (clustree)")
ggsave(file.path(figures.dir, "FigS1d_clustree_metaclusters.pdf"), p, width = 10, height = 6)

## ---- S1e: alluvial, samples -> meta-clusters (one line per meta-spot) ----
df_plot <- seurat_obj@meta.data |>
  transmute(sample = factor(sample, levels = sort(unique(sample))),
            meta = factor(merged_meta_clusters, levels = as.character(0:10))) |>
  count(sample, meta, name = "n")
p <- ggplot(df_plot, aes(axis1 = sample, axis2 = meta, y = n)) +
  geom_alluvium(aes(fill = sample), alpha = 0.65, knot.pos = 0.4, width = 1/10, show.legend = FALSE) +
  geom_stratum(width = 1/5, fill = "white", color = "grey50") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.6) +
  theme_void()
ggsave(file.path(figures.dir, "FigS1e_alluvial_samples_to_metaclusters.pdf"), p, width = 9.5, height = 12)
