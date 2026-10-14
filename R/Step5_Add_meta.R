# Step 5 - Add meta-cluster / stromal sub-cluster labels to the per-sample (spot-level) Seurat objects
# Output: Processed_seurat/processed_seurat_<Sample>_withmeta.RData (used for spatial plots and Supplementary Fig. S4e)
#         Tables/spot_labels_all48.rds/.csv (sample, barcode, merged_meta_clusters, stromal_subcluster; 193,418 spots)

source("R/00_config.R")

library(Seurat)

load(file.path(harmony.dir, "scaled_data_list_10.RData"))
spot_labels <- list()

for (Sample in Sample.names) {
  load(file.path(processed.dir, paste0("processed_seurat_", Sample, ".RData")))   # -> sobj
  df <- scaled_data_list[[Sample]]
  bc <- intersect(colnames(sobj), rownames(df))

  mc <- sc <- setNames(rep(NA_character_, ncol(sobj)), colnames(sobj))
  mc[bc] <- as.character(df[bc, "merged_meta_clusters"])
  sc[bc] <- as.character(df[bc, "stromal_subcluster"])

  sobj$merged_meta_clusters    <- mc
  sobj$metacluster_name        <- unname(meta_names[mc])
  sobj$stromal_subcluster      <- sc
  sobj$stromal_subcluster_name <- unname(stromal_names[sc])
  # Combined label: stromal sub-cluster name for stromal spots, meta-cluster name otherwise
  sobj$meta_stromal_label <- ifelse(!is.na(sobj$stromal_subcluster_name), sobj$stromal_subcluster_name, sobj$metacluster_name)

  save(sobj, file = file.path(processed.dir, paste0("processed_seurat_", Sample, "_withmeta.RData")))
  spot_labels[[Sample]] <- data.frame(barcode = colnames(sobj), sample = Sample,
                                      merged_meta_clusters = mc, stromal_subcluster = sc)
}

# Spot-level labels used by the deconvolution, YAP-Rev, CellChat and gradient analyses (analyses/)
spot_labels <- do.call(rbind, spot_labels); rownames(spot_labels) <- NULL
saveRDS(spot_labels, file.path(tables.dir, "spot_labels_all48.rds"))
write.csv(spot_labels, file.path(tables.dir, "spot_labels_all48.csv"), row.names = FALSE)
