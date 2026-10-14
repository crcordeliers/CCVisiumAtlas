# Fig. 1b - meta-spot construction (example sample #74_D1: per-sample clusters at resolution 10, UMAP and section) and
# UMAP of the 3,231 meta-spots coloured by meta-cluster. One figure per PDF file.

source("R/00_config.R")

library(Seurat)
library(RColorBrewer)
library(ggplot2)

save_pdf <- function(p, name, w = 9, h = 7) {
  pdf(file.path(figures.dir, name), width = w, height = h)
  print(p)
  dev.off()
}

# Meta-spots of one sample (resolution-10 clusters)
load(file.path(processed.dir, "processed_seurat_Ech_074_D1.RData"))   # -> sobj
num_clusters  <- length(unique(sobj$Cluster_res_10))
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_clusters)
save_pdf(DimPlot(sobj, group.by = "Cluster_res_10", cols = color_palette),
         "Fig1b_Ech_074_D1_metaspots_UMAP.pdf")
save_pdf(SpatialDimPlot(sobj, group.by = "Cluster_res_10", pt.size.factor = 2.5, alpha = 0.8) +
           scale_fill_manual(values = color_palette),
         "Fig1b_Ech_074_D1_metaspots_spatial.pdf")

# UMAP of the 3,231 meta-spots coloured by meta-cluster
seurat_obj <- readRDS(file.path(harmony.dir, "seurat_object_res_10.rds"))
save_pdf(DimPlot(seurat_obj, reduction = "umap", group.by = "merged_meta_clusters", pt.size = 1) +
           scale_color_manual(values = meta_colors, labels = meta_names, breaks = meta_order),
         "Fig1b_metaclusters_UMAP.pdf")

# Meta-clusters on the tissue section
load(file.path(processed.dir, "processed_seurat_Ech_074_D1_withmeta.RData"))   # -> sobj
save_pdf(SpatialDimPlot(sobj, group.by = "merged_meta_clusters", pt.size.factor = 2.5, alpha = 1) +
           scale_fill_manual(values = meta_colors, labels = meta_names, breaks = meta_order) +
           theme(legend.position = "bottom"),
         "Fig1b_Ech_074_D1_metaclusters_spatial.pdf")
