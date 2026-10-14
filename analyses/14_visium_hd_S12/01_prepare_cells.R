#!/usr/bin/env Rscript
# =====================================================================================================================
# Visium HD, step 1 - per-cell table from the segmented Space Ranger output (one sample per run)
#
# Sections processed with Space Ranger 4.0.1 using a custom segmentation: nuclei detected on the CZI image with
# InstanSeg (QuPath extension) and expanded by 5 µm. Builds, for each segmented cell: barcode, polygon centroid
# (full-resolution pixels and µm), area, UMI, number of genes, Space Ranger graph-based cluster.
# Samples: Visium_HD_D1_1856_074_D1 (= Visium Ech_074_D1, relapse), Visium_HD2_A1_1341_395_C1 (Ech_395_C1, relapse),
#          Visium_HD_A1_1898_2177_B1 (Ech_21-77_B1, no relapse), Visium_HD2-D1-1656-295B1 (Ech_295_B1, no relapse).
# Input : <HD_DIR>/data/<SAMPLE>/outs (filtered_feature_cell_matrix.h5, cell_segmentations.geojson, spatial/)
# Output: <HD_DIR>/per_sample/<SAMPLE>/objects/cells_meta_<SAMPLE>.rds, counts_<SAMPLE>.rds
# Usage : SAMPLE=Visium_HD_D1_1856_074_D1 Rscript analyses/14_visium_hd_S12/01_prepare_cells.R
# =====================================================================================================================
suppressPackageStartupMessages({ library(Seurat); library(jsonlite); library(Matrix) })

sample_id <- Sys.getenv("SAMPLE"); stopifnot(nzchar(sample_id))
hd_dir <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd")
sroot  <- file.path(hd_dir, "per_sample", sample_id)
outs   <- file.path(hd_dir, "data", sample_id, "outs")
dir.create(file.path(sroot, "objects"), showWarnings = FALSE, recursive = TRUE)

um_per_px <- fromJSON(file.path(outs, "spatial", "scalefactors_json.json"))$microns_per_pixel
cnt <- Read10X_h5(file.path(outs, "filtered_feature_cell_matrix.h5"))
feats <- fromJSON(file.path(outs, "cell_segmentations.geojson"), simplifyVector = FALSE)$features

# centroid and area of each polygon (signed-area formula)
poly_stats <- function(f) {
  ring <- f$geometry$coordinates[[1]]
  x <- vapply(ring, `[[`, numeric(1), 1); y <- vapply(ring, `[[`, numeric(1), 2)
  n <- length(x); x0 <- x[-n]; y0 <- y[-n]; x1 <- x[-1]; y1 <- y[-1]
  cr <- x0 * y1 - x1 * y0; a <- sum(cr) / 2
  if (abs(a) < 1e-9) return(c(f$properties$cell_id, mean(x), mean(y), 0))
  c(f$properties$cell_id, sum((x0 + x1) * cr) / (6 * a), sum((y0 + y1) * cr) / (6 * a), abs(a))
}
st <- do.call(rbind, lapply(feats, poly_stats))
cells <- data.frame(cell_id = as.integer(st[, 1]), x_px = st[, 2], y_px = st[, 3], area_um2 = st[, 4] * um_per_px^2)
cells$barcode <- sprintf("cellid_%09d-1", cells$cell_id)
stopifnot(all(colnames(cnt) %in% cells$barcode))
cells <- cells[match(colnames(cnt), cells$barcode), ]
cells$x_um <- cells$x_px * um_per_px; cells$y_um <- cells$y_px * um_per_px
cells$nUMI <- colSums(cnt); cells$nGene <- colSums(cnt > 0)
gcl <- read.csv(file.path(outs, "analysis", "clustering", "gene_expression_graphclust", "clusters.csv"))
cells$graphclust <- gcl$Cluster[match(cells$barcode, gcl$Barcode)]

saveRDS(cells, file.path(sroot, "objects", paste0("cells_meta_", sample_id, ".rds")))
saveRDS(cnt, file.path(sroot, "objects", paste0("counts_", sample_id, ".rds")))
message(sample_id, ": ", ncol(cnt), " cells, median area ", round(median(cells$area_um2), 1), " µm2")
