# Step 1 - Per-sample Seurat objects from Space Ranger outputs
# Normalisation, scaling, per-sample PCA/SNN graph and CNV fractions.

source("R/00_config.R")

library(Seurat)

CNV_fractions <- readRDS(cnv.file)

for (Sample in Sample.names) {

  seurat_file <- file.path(allseurat.dir, paste0("seurat_", Sample, ".RData"))
  if (file.exists(seurat_file)) next

  sobj <- Load10X_Spatial(file.path(raw.dir, Sample, "outs"), filename = "filtered_feature_bc_matrix.h5")

  # Normalisation + scaling (all genes: ScaleData is called before FindVariableFeatures)
  sobj <- NormalizeData(sobj, assay = "Spatial")
  sobj <- ScaleData(sobj, assay = "Spatial")

  # Per-sample graph (used in Step 2 for the resolution-10 clustering)
  sobj <- FindVariableFeatures(sobj, assay = "Spatial")
  sobj <- RunPCA(sobj, assay = "Spatial", verbose = FALSE)
  sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:30)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution = 0.8)
  sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:30)

  # CNV fraction per spot (fastCNV)
  if (Sample %in% names(CNV_fractions)) {
    cnv_data <- as.data.frame(CNV_fractions[Sample])
    colnames(cnv_data) <- "CNV_fraction"
    sobj <- AddMetaData(sobj, cnv_data)
  }

  save(sobj, file = seurat_file)
  rm(sobj); gc()
}
