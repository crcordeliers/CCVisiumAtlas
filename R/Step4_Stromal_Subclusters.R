# Step 4 - Stromal sub-clusters
# Meta-spots of the two stromal meta-clusters (MC2 Stroma-A, MC4 Stroma-B; n = 774) are re-clustered at
# resolution 0.5 on the (subsetted) meta-spot graph, giving 5 stromal sub-clusters. Markers of the
# sub-clusters give the stromal part of Supplementary Table S2.

source("R/00_config.R")

library(Seurat)
library(GSVA)
library(dplyr)

seurat_obj <- readRDS(file.path(harmony.dir, "seurat_object_res_10.rds"))
load(file.path(data.dir, "signatures.RData"))

Cluster2_4 <- subset(seurat_obj, subset = merged_meta_clusters %in% stromal_metaclusters)

# FindClusters overwrites RNA_snn_res.0.5 / seurat_clusters with the sub-clusters:
# the meta-cluster stays available in `merged_meta_clusters`.
Cluster2_4 <- FindClusters(Cluster2_4, resolution = res_stromal)
Cluster2_4$stromal_subcluster      <- as.character(Cluster2_4$seurat_clusters)
Cluster2_4$stromal_subcluster_name <- unname(stromal_names[Cluster2_4$stromal_subcluster])

## ---- Markers: Supplementary Table S2 (stromal sub-clusters) ----
seurat.markers <- FindAllMarkers(Cluster2_4, only.pos = TRUE)
seurat.markers$cluster_name <- unname(stromal_names[as.character(seurat.markers$cluster)])
write.csv(seurat.markers, file.path(tables.dir, "Markers_stromal_subclusters_all.csv"), row.names = FALSE)

top20 <- seurat.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup()
write.csv(top20, file.path(tables.dir, "TableS2_stromal_subclusters_top20.csv"), row.names = FALSE)

ClusterListStromalSubcluster <- split(top20$gene, top20$cluster)
save(ClusterListStromalSubcluster, file = file.path(harmony.dir, "ClusterListStromalSubcluster.RData"))

## ---- Signature scores (ssGSEA) on the stromal meta-spots ----
sig_gsva <- as.data.frame(t(gsva(ssgseaParam(as.matrix(Cluster2_4[["RNA"]]$scale.data), signatures, normalize = TRUE))))
Cluster2_4 <- AddMetaData(Cluster2_4, sig_gsva)

saveRDS(Cluster2_4, file = file.path(harmony.dir, "Cluster2_4.rds"))

## ---- Map sub-clusters back to the full object and to spots ----
# (stored as character: assigning a factor into an NA vector would store its integer codes 1..5)
sc <- setNames(rep(NA_character_, ncol(seurat_obj)), colnames(seurat_obj))
sc[colnames(Cluster2_4)] <- Cluster2_4$stromal_subcluster
seurat_obj$stromal_subcluster      <- sc
seurat_obj$stromal_subcluster_name <- unname(stromal_names[sc])
saveRDS(seurat_obj, file = file.path(harmony.dir, "seurat_object_res_10_with5subsets.rds"))

load(file.path(harmony.dir, "scaled_data_list_10.RData"))
sc <- setNames(seurat_obj$stromal_subcluster, colnames(seurat_obj))
for (Sample in names(scaled_data_list)) {
  scaled_data_list[[Sample]]$stromal_subcluster <- unname(sc[scaled_data_list[[Sample]]$metaspot])
}
save(scaled_data_list, file = file.path(harmony.dir, "scaled_data_list_10.RData"))
