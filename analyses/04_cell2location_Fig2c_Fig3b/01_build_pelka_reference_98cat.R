#!/usr/bin/env Rscript
# =====================================================================================================================
# Cell2location reference, step 1: Pelka et al. 2021 single-cell reference at the finest level (98 categories)
# Paper: Methods "Cell2location deconvolution"; used for Fig. 2c, Fig. 3b (and S8h, matCAF = cS26 + cS27 + cS28).
#
# Input  : GSE178341_crc10x_full_c295v4_submit.h5 (raw counts, GEO) and crc10x_tSNE_cl_global.tsv (Broad Single Cell
#          Portal SCP1162; column ClusterFull = 98 categories, "Tumor " prefix for malignant epithelial states).
# Output : <C2L_DIR>/reference/  matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz (all 370,115 cells)
#                                reference_cellIDs_98cat.tsv (down-sampled to max 1,000 cells per category, seed 42)
# Usage  : Rscript analyses/04_cell2location_Fig2c_Fig3b/01_build_pelka_reference_98cat.R
# Env    : R 4.5 (Seurat 5, rhdf5); run on a cluster node (~60 GB RAM).
# =====================================================================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

pelka_dir <- Sys.getenv("CRCVISIUM_PELKA", "~/data/pelka2021_raw")     # .h5 + broad_scp/crc10x_tSNE_cl_global.tsv
c2l_dir   <- Sys.getenv("CRCVISIUM_C2L",   "~/data/cell2location")
ref_dir   <- file.path(c2l_dir, "reference"); dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

counts <- Read10X_h5(file.path(pelka_dir, "GSE178341_crc10x_full_c295v4_submit.h5"))
tsne_file <- file.path(pelka_dir, "broad_scp", "crc10x_tSNE_cl_global.tsv")
tsne <- read.delim(tsne_file, skip = 2, header = FALSE, stringsAsFactors = FALSE,
                   col.names = strsplit(readLines(tsne_file, n = 1), "\t")[[1]])
rownames(tsne) <- tsne$NAME
stopifnot(all(colnames(counts) %in% rownames(tsne)))
tsne <- tsne[colnames(counts), ]

# "/" is not allowed in category names by spacexr (e.g. "Stem/TA-like"): replaced by "-"
cell_type <- factor(gsub("/", "-", tsne$ClusterFull, fixed = TRUE))
stopifnot(nlevels(cell_type) == 98, !anyNA(cell_type))

# full matrix (genes x cells)
writeMM(counts, file.path(ref_dir, "matrix.mtx"))
writeLines(colnames(counts), file.path(ref_dir, "barcodes.tsv"))
writeLines(rownames(counts), file.path(ref_dir, "features.tsv"))
system(paste(c("gzip -f", file.path(ref_dir, c("matrix.mtx", "barcodes.tsv", "features.tsv"))), collapse = " "))

# down-sampling: max 1,000 cells per category, seed 42 (72,382 cells; same subset as the RCTD reference)
obj <- CreateSeuratObject(counts = counts, meta.data = data.frame(cell_type = cell_type, row.names = colnames(counts)))
set.seed(42)
Idents(obj) <- "cell_type"
obj_ds <- subset(obj, downsample = 1000)
cellids <- data.frame(cellID = colnames(obj_ds), cell_type = as.character(obj_ds$cell_type))
write.table(cellids, file.path(ref_dir, "reference_cellIDs_98cat.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
message(nrow(cellids), " cells, ", length(unique(cellids$cell_type)), " categories")
