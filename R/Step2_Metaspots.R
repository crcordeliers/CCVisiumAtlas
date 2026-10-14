# Step 2 - Meta-spots
# Per-sample clustering at resolution 10; each cluster becomes a "meta-spot" whose expression is the
# mean of the scaled expression (scale.data) of its spots. Also averages signature scores and CNV
# fraction per meta-spot.

source("R/00_config.R")

library(Seurat)
library(dplyr)
library(future.apply)

load(file.path(data.dir, "signatures.RData"))   # -> signatures (66 gene sets)

## ---- 2a. Resolution-10 clustering + spot-level signature scores ----
plan(multisession, workers = n_workers)
future_lapply(Sample.names, function(Sample) {
  load(file.path(allseurat.dir, paste0("seurat_", Sample, ".RData")))   # -> sobj

  compute_signature_scores <- function(expr_mat, suffix = "") {
    sig_scores <- as.data.frame(t(sapply(signatures, function(genes) {
      genes <- genes[genes %in% rownames(expr_mat)]
      if (length(genes) > 0) colMeans(expr_mat[genes, , drop = FALSE], na.rm = TRUE) else rep(NA, ncol(expr_mat))
    })))
    rownames(sig_scores) <- paste0(rownames(sig_scores), suffix)
    t(sig_scores)
  }
  sobj <- AddMetaData(sobj, compute_signature_scores(as.matrix(sobj@assays$Spatial$scale.data), "_scale"))
  sobj <- AddMetaData(sobj, compute_signature_scores(as.matrix(sobj@assays$Spatial$data), ""))

  # Clustering on the per-sample graph built in Step 1
  sobj <- FindClusters(sobj, resolution = res_metaspot, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = 1:15)
  sobj <- RunTSNE(sobj, dims = 1:15)
  sobj@meta.data[[paste0("Cluster_res_", res_metaspot)]] <- Idents(sobj)

  save(sobj, file = file.path(processed.dir, paste0("processed_seurat_", Sample, ".RData")))
}, future.seed = TRUE)
plan(sequential)

## ---- 2b. Meta-spot means (genes, signatures + CNV) ----
cluster_col <- paste0("Cluster_res_", res_metaspot)

plan(multisession, workers = n_workers)
results <- future_lapply(Sample.names, function(Sample) {
  load(file.path(processed.dir, paste0("processed_seurat_", Sample, ".RData")))   # -> sobj
  md <- sobj@meta.data

  gene_scaled_data <- as.data.frame(t(GetAssayData(sobj, assay = "Spatial", layer = "scale.data")))
  gene_scaled_data[[cluster_col]] <- as.numeric(as.character(md[[cluster_col]]))

  feats <- c("CNV_fraction", names(signatures), paste0(names(signatures), "_scale"))
  signature_cnv_data <- md[, intersect(feats, colnames(md)), drop = FALSE]
  signature_cnv_data[[cluster_col]] <- as.numeric(as.character(md[[cluster_col]]))

  gene_mean_exp <- gene_scaled_data %>%
    group_by(across(all_of(cluster_col))) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  gene_mean_exp$sample_res <- metaspot_id(Sample, gene_mean_exp[[cluster_col]])

  signature_cnv_mean <- signature_cnv_data %>%
    group_by(across(all_of(cluster_col))) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  signature_cnv_mean$sample_res <- metaspot_id(Sample, signature_cnv_mean[[cluster_col]])

  spot_data <- cbind(gene_scaled_data, signature_cnv_data[, setdiff(colnames(signature_cnv_data), cluster_col), drop = FALSE])
  spot_data$Sample <- Sample

  list(gene_expression = setNames(list(gene_mean_exp), paste0("res_", res_metaspot)),
       signature_cnv   = setNames(list(signature_cnv_mean), paste0("res_", res_metaspot)),
       spot_data       = spot_data)
}, future.seed = TRUE)
plan(sequential)

gene_expression_list  <- setNames(lapply(results, `[[`, "gene_expression"), Sample.names)
signature_cnv_list    <- setNames(lapply(results, `[[`, "signature_cnv"), Sample.names)
scaled_data_list      <- setNames(lapply(results, `[[`, "spot_data"), Sample.names)

save(gene_expression_list,  file = file.path(results_res.dir, "gene_expression_list.RData"))
save(signature_cnv_list,    file = file.path(results_res.dir, "signature_cnv_list.RData"))
save(scaled_data_list,      file = file.path(results_res.dir, "scaled_data_list.RData"))
