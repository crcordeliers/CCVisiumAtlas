# Step 6 - Pelka et al. 2021 cell-state signatures (ssGSEA) on meta-spots (used in Figs 2 and 4)

source("R/00_config.R")

library(Seurat)
library(GSVA)

seurat_obj <- readRDS(file.path(harmony.dir, "seurat_object_res_10_with5subsets.rds"))



.clean_genes <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  x <- trimws(x)
  x <- x[!grepl("^DESC\\b|^DESC:", x, ignore.case = TRUE)]
  x <- gsub("[^A-Za-z0-9_-]", "", x)
  x <- unique(x)
  x[nzchar(x)]
}

read_signatures_file <- function(path) {
  df <- read.table(path, sep="\t", header=TRUE, check.names=FALSE,
                   quote="", comment.char="", fill=TRUE, stringsAsFactors=FALSE)
  sigs <- lapply(df, .clean_genes)
  names(sigs) <- colnames(df)
  sigs[vapply(sigs, length, 1L) > 0]
}


compute_ssgsea_df <- function(seurat_obj, signatures) {
  mat <- as.matrix(seurat_obj@assays$RNA$scale.data)
  param <- ssgseaParam(mat, signatures, normalize = TRUE)
  gs <- gsva(param)
  as.data.frame(t(gs))
}

added_cols_by_level <- list()


Low   = file.path(data.dir, "Pelka2021_ClusterLow.txt")
sigs1 <- read_signatures_file(Low) 
ss1 <- compute_ssgsea_df(seurat_obj, sigs1)

Top = file.path(data.dir, "Pelka2021_ClusterTop.txt")
sigs2 <- read_signatures_file(Top)
ss2 <- compute_ssgsea_df(seurat_obj, sigs2)

# Add to metadata
seurat_obj <- AddMetaData(seurat_obj, ss1)
added_cols_by_level[[Low]] <- colnames(ss1)

seurat_obj <- AddMetaData(seurat_obj, ss2)
added_cols_by_level[[Top]] <- colnames(ss2)

save(seurat_obj, file = file.path(harmony.dir, "seurat_obj_with_pelka.RData"))
