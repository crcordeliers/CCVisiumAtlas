#!/usr/bin/env Rscript
# =====================================================================================================================
# Pelka et al. 2021 references for RCTD (Visium HD cell annotation, Supplementary Fig. S12)
#   17 broad categories : GEO clMidwayPr for non-epithelial cells; epithelial cells split into malignant (Broad SCP
#                         ClusterMidway "EpiT") and normal ("Epi"/"EpiN") epithelium.
#   98 fine categories  : Broad SCP ClusterFull ("Tumor " prefix for malignant states).
# Both are down-sampled to max 1,000 cells per category (seed 42) and saved as SingleCellExperiment objects.
# Inputs : <PELKA_DIR>/GSE178341_crc10x_full_c295v4_submit.h5, GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz,
#          broad_scp/crc10x_tSNE_cl_global.tsv
# Outputs: <HD_DIR>/reference/pelka_reference_17groups.rds, pelka_reference_98cats.rds, cell_type_17groups.tsv.gz
# Usage  : Rscript analyses/14_visium_hd_S12/00_build_rctd_references.R
# =====================================================================================================================
suppressPackageStartupMessages({ library(Seurat); library(SingleCellExperiment); library(Matrix) })

pelka_dir <- Sys.getenv("CRCVISIUM_PELKA", "~/data/pelka2021_raw")
hd_dir    <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd")
ref_dir   <- file.path(hd_dir, "reference"); dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

counts <- Read10X_h5(file.path(pelka_dir, "GSE178341_crc10x_full_c295v4_submit.h5"))
clust  <- read.csv(file.path(pelka_dir, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"), stringsAsFactors = FALSE)
colnames(clust)[1] <- "cellID"                       # first column is mislabelled "sampleID" in GEO
rownames(clust) <- clust$cellID
clust <- clust[colnames(counts), ]
tsne_file <- file.path(pelka_dir, "broad_scp", "crc10x_tSNE_cl_global.tsv")
tsne <- read.delim(tsne_file, skip = 2, header = FALSE, stringsAsFactors = FALSE,
                   col.names = strsplit(readLines(tsne_file, n = 1), "\t")[[1]])
rownames(tsne) <- tsne$NAME

# ---- 17 categories ----
map_direct <- c(B = "B", Plasma = "Plasma", TCD4 = "CD4 T", TCD8 = "CD8 T", Macro = "Macrophage-Monocyte",
                Mono = "Macrophage-Monocyte", Tgd = "T-NK-ILC-Other", TZBTB16 = "T-NK-ILC-Other", ILC = "T-NK-ILC-Other",
                NK = "NK", Endo = "Endothelial", DC = "Dendritic cell", Fibro = "Fibroblast-CAF", Mast = "Mast",
                Granulo = "Granulocyte", Peri = "Pericyte", SmoothMuscle = "Smooth muscle", Schwann = "Schwann")
is_epi <- clust$clTopLevel == "Epi"
epi_midway <- tsne[clust$cellID[is_epi], "ClusterMidway"]
stopifnot(all(epi_midway %in% c("EpiT", "Epi", "EpiN")))
grp17 <- rep(NA_character_, nrow(clust))
grp17[!is_epi] <- unname(map_direct[clust$clMidwayPr[!is_epi]])
grp17[is_epi] <- ifelse(epi_midway == "EpiT", "Malignant epithelial", "Normal epithelial")
grp17 <- factor(grp17); stopifnot(nlevels(grp17) == 17, !anyNA(grp17))
write.table(data.frame(cellID = colnames(counts), cell_type = as.character(grp17)),
            gzfile(file.path(ref_dir, "cell_type_17groups.tsv.gz")), sep = "\t", row.names = FALSE, quote = FALSE)

# ---- 98 categories ("/" not allowed by spacexr -> "-") ----
grp98 <- factor(gsub("/", "-", tsne[colnames(counts), "ClusterFull"], fixed = TRUE))
stopifnot(nlevels(grp98) == 98)

downsample_sce <- function(labels) {
  obj <- CreateSeuratObject(counts = counts, meta.data = data.frame(cell_type = labels, row.names = colnames(counts)))
  set.seed(42)
  Idents(obj) <- "cell_type"
  obj <- subset(obj, downsample = 1000)
  SingleCellExperiment(assays = list(counts = LayerData(obj, assay = "RNA", layer = "counts")),
                       colData = S4Vectors::DataFrame(cell_type = obj$cell_type))
}
saveRDS(downsample_sce(grp17), file.path(ref_dir, "pelka_reference_17groups.rds"))
saveRDS(downsample_sce(grp98), file.path(ref_dir, "pelka_reference_98cats.rds"))
