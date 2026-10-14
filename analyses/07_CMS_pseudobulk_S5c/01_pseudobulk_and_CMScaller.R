#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S5c (step 1) - CMS classification of each Visium sample on a pseudo-bulk
# Raw counts summed over all in-tissue spots of the sample ("global" pseudo-bulk) and over tumour spots only (meta-clusters
# 0, 1, 7, 8; "tumoral", used for comparison); CMScaller v2.0.1, RNAseq = TRUE, FDR = 0.05, seed 42.
# Inputs : <CRCVISIUM_RAW>/<Sample>/outs/filtered_feature_bc_matrix ; Tables/spot_labels_all48.rds (R/Step5)
# Output : Tables/CMS_pseudobulk_48.csv (sample, CMS_global, FDR_global, CMS_tumoral, FDR_tumoral)
# Usage  : Rscript analyses/07_CMS_pseudobulk_S5c/01_pseudobulk_and_CMScaller.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(CMScaller) })

labels <- readRDS(file.path(tables.dir, "spot_labels_all48.rds"))
pb <- lapply(Sample.names, function(s) {
  cnt <- Read10X(file.path(raw.dir, s, "outs", "filtered_feature_bc_matrix"))
  lab <- labels[labels$sample == s, ]
  bc  <- intersect(colnames(cnt), lab$barcode)
  tum <- bc[as.character(lab$merged_meta_clusters[match(bc, lab$barcode)]) %in% c("0", "1", "7", "8")]
  list(global = rowSums(cnt[, bc, drop = FALSE]), tumoral = rowSums(cnt[, tum, drop = FALSE]))
})
names(pb) <- Sample.names
mat <- function(k) { m <- sapply(pb, `[[`, k); m }
cms_g <- CMScaller(mat("global"),  RNAseq = TRUE, rowNames = "symbol", FDR = 0.05, doPlot = FALSE, seed = 42)
cms_t <- CMScaller(mat("tumoral"), RNAseq = TRUE, rowNames = "symbol", FDR = 0.05, doPlot = FALSE, seed = 42)
out <- data.frame(sample = rownames(cms_g), CMS_global = cms_g$prediction, FDR_global = cms_g$FDR,
                  CMS_tumoral = cms_t[rownames(cms_g), "prediction"], FDR_tumoral = cms_t[rownames(cms_g), "FDR"])
write.csv(out, file.path(tables.dir, "CMS_pseudobulk_48.csv"), row.names = FALSE)
print(table(out$CMS_global, useNA = "ifany"))     # CMS1 10, CMS2 14, CMS3 3, CMS4 15, unclassified 6
