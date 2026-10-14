#!/usr/bin/env Rscript
# =====================================================================================================================
# Visium HD, step 2 - RCTD (spacexr, doublet mode) with the 17-category Pelka reference, segmented cells
# UMI_min = 50 (the RCTD default of 100 would remove ~22% of the cells; median ~220 UMI per cell).
# Input : objects from 01_prepare_cells.R ; <HD_DIR>/reference/pelka_reference_17groups.rds (00_*.R)
# Output: <HD_DIR>/per_sample/<SAMPLE>/objects/rctd17_doublet_all_<SAMPLE>.rds
# Usage : SAMPLE=<id> Rscript analyses/14_visium_hd_S12/02_rctd17_doublet.R   (SLURM, 32 cores, ~30 min)
# =====================================================================================================================
suppressPackageStartupMessages({ library(spacexr); library(SpatialExperiment); library(SingleCellExperiment); library(Matrix) })

sample_id <- Sys.getenv("SAMPLE"); stopifnot(nzchar(sample_id))
hd_dir <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd")
sroot  <- file.path(hd_dir, "per_sample", sample_id)
max_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

cells <- readRDS(file.path(sroot, "objects", paste0("cells_meta_", sample_id, ".rds")))
cnt   <- readRDS(file.path(sroot, "objects", paste0("counts_", sample_id, ".rds")))
stopifnot(identical(colnames(cnt), cells$barcode))
coords <- as.matrix(cells[, c("x_um", "y_um")]); colnames(coords) <- c("x", "y"); rownames(coords) <- cells$barcode

reference_se <- readRDS(file.path(hd_dir, "reference", "pelka_reference_17groups.rds"))
spe <- SpatialExperiment(assay = list(counts = cnt), spatialCoords = coords)
rctd_data <- createRctd(spe, reference_se, cell_type_col = "cell_type", UMI_min = 50)
result <- runRctd(rctd_data, rctd_mode = "doublet", max_cores = max_cores)
saveRDS(result, file.path(sroot, "objects", paste0("rctd17_doublet_all_", sample_id, ".rds")))
print(table(colData(result)$spot_class))
