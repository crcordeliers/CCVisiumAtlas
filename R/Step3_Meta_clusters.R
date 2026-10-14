# Step 3 - Meta-clusters
# Integration of the 3,231 meta-spots (PCA on all genes -> Harmony on sample -> Louvain, resolution 0.5),
# signature scores, cluster markers (Supplementary Table S2, meta-clusters) and propagation of the
# meta-cluster labels to every spot.
#
# Requirements:
#  * harmony 1.x (harmony >= 2.0 uses a different algorithm and gives different meta-clusters).
#  * Default seeds (RunPCA/RunUMAP seed.use = 42).
#  * `features` passed to RunPCA built as below (column names of the meta-spot table, including the leading
#    "Cluster_res_10" column, which Seurat drops because it is not in scale.data); another list of the same
#    genes changes the irlba starting vector and the result.

source("R/00_config.R")

library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(tibble)
library(GSVA)

stopifnot("harmony 1.x is required (e.g. 1.2.4)" = packageVersion("harmony") < "2.0.0")

data_clin <- readRDS(file.path(data.dir, "data_clin_48_filtered.RDS"))
load(file.path(results_res.dir, "gene_expression_list.RData"))
load(file.path(results_res.dir, "signature_cnv_list.RData"))
load(file.path(data.dir, "signatures.RData"))

res_key <- paste0("res_", res_metaspot)

## ---- Meta-spot matrix and metadata ----
all_samples_mean_exp <- as.data.frame(do.call(rbind, lapply(gene_expression_list, `[[`, res_key)))
all_samples_metadata <- as.data.frame(do.call(rbind, lapply(signature_cnv_list, `[[`, res_key)))
rownames(all_samples_mean_exp) <- all_samples_mean_exp$sample_res
rownames(all_samples_metadata) <- all_samples_metadata$sample_res

meta <- data.frame(metaspot = rownames(all_samples_mean_exp),
                   sample   = sapply(strsplit(rownames(all_samples_mean_exp), "_res"), `[`, 1)) %>%
  left_join(data_clin, by = "sample") %>%
  left_join(all_samples_metadata %>% select(-all_of(paste0("Cluster_res_", res_metaspot))), by = c("metaspot" = "sample_res"))
rownames(meta) <- meta$metaspot

## ---- Seurat object on meta-spots (the averaged scaled matrix is used directly as scale.data) ----
metaspot_exp <- dplyr::select(all_samples_mean_exp, -sample_res)   # keeps the Cluster_res_10 column: see note above
seurat_obj <- CreateSeuratObject(counts = t(metaspot_exp), meta.data = meta)
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj[["RNA"]]$scale.data <- t(metaspot_exp)

seurat_obj <- RunPCA(seurat_obj, features = rownames(t(metaspot_exp)), npcs = 50)   # seed.use = 42 (default)
seurat_obj$batch <- sapply(strsplit(rownames(seurat_obj@meta.data), "_res"), `[`, 1)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "batch")
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:15)             # seed.use = 42 (default)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = res_meta)
seurat_obj <- RunTSNE(seurat_obj)

seurat_obj$merged_meta_clusters <- as.character(seurat_obj$seurat_clusters)
seurat_obj$metacluster_name     <- unname(meta_names[seurat_obj$merged_meta_clusters])

## ---- Markers: Supplementary Table S2 (meta-clusters) ----
# The "data" layer of this object contains the averaged scaled values (not log-normalised counts), so
# avg_log2FC values are not true log fold-changes; the signature definition (ranking + threshold) uses them.
seurat.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
seurat.markers$cluster_name <- unname(meta_names[as.character(seurat.markers$cluster)])
write.csv(seurat.markers, file.path(tables.dir, "Markers_metaclusters_all.csv"), row.names = FALSE)

top20 <- seurat.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup()
write.csv(top20, file.path(tables.dir, "TableS2_metaclusters_top20.csv"), row.names = FALSE)

ClusterList <- split(top20$gene, top20$cluster)
save(ClusterList, file = file.path(harmony.dir, "ClusterList.RData"))

## ---- Signature scores per meta-spot (mean of scaled expression and ssGSEA) ----
expr_mat_scale <- as.matrix(seurat_obj[["RNA"]]$scale.data)
sig_mean <- sapply(signatures, function(genes) {
  genes <- genes[genes %in% rownames(expr_mat_scale)]
  if (length(genes) > 0) colMeans(expr_mat_scale[genes, , drop = FALSE], na.rm = TRUE) else rep(NA, ncol(expr_mat_scale))
})
colnames(sig_mean) <- paste0(colnames(sig_mean), "_mean")
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(sig_mean))

sig_gsva <- as.data.frame(t(gsva(ssgseaParam(expr_mat_scale, signatures, normalize = TRUE))))
seurat_obj <- AddMetaData(seurat_obj, sig_gsva)

saveRDS(seurat_obj, file = file.path(harmony.dir, "seurat_object_res_10.rds"))

## ---- Propagate meta-cluster labels to spots ----
load(file.path(results_res.dir, "scaled_data_list.RData"))
mc <- setNames(seurat_obj$merged_meta_clusters, colnames(seurat_obj))
for (Sample in names(scaled_data_list)) {
  dt <- scaled_data_list[[Sample]]
  dt$metaspot <- metaspot_id(Sample, dt[[paste0("Cluster_res_", res_metaspot)]])
  dt$merged_meta_clusters <- unname(mc[dt$metaspot])
  scaled_data_list[[Sample]] <- dt
}
save(scaled_data_list, file = file.path(harmony.dir, "scaled_data_list_10.RData"))
