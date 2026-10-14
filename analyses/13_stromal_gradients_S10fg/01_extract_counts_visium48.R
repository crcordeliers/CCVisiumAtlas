#!/usr/bin/env Rscript
# =====================================================================================================================
# Raw counts per spot of the genes used for the stromal gradients (Supplementary Fig. S10f,g), 48 Visium samples.
# Genes: MSigDB HALLMARK_TGF_BETA_SIGNALING and NABA_COLLAGENS (genesets_MSigDB.tsv) + EPCAM, KRT20 + total UMI.
# Input : All_Seurat/seurat_<Sample>.RData (R/Step1_Object_creation.R; raw counts of the Spatial assay)
# Output: <GRADIENT_DIR>/<Sample>_genesets.csv.gz (barcode, sample, nCount, genes); <GRADIENT_DIR>/visium48_EPCAM_KRT20.csv.gz
# Usage : Rscript analyses/13_stromal_gradients_S10fg/01_extract_counts_visium48.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

gradient_dir <- Sys.getenv("CRCVISIUM_GRADIENT", file.path(work.dir, "gradients"))
dir.create(gradient_dir, recursive = TRUE, showWarnings = FALSE)
gs <- read.delim("analyses/13_stromal_gradients_S10fg/genesets_MSigDB.tsv", header = FALSE)
G <- unique(unlist(strsplit(gs$V2, ",")))

ep <- list()
for (smp in Sample.names) {
  e <- new.env(); load(file.path(allseurat.dir, paste0("seurat_", smp, ".RData")), envir = e)
  so <- get(Filter(function(n) inherits(get(n, e), "Seurat"), ls(e))[1], e); rm(e)
  m <- GetAssayData(so, assay = if ("Spatial" %in% Assays(so)) "Spatial" else DefaultAssay(so), layer = "counts")
  stopifnot(all(abs(m@x - round(m@x)) < 1e-8))
  g <- intersect(G, rownames(m))
  write.csv(data.frame(barcode = colnames(m), sample = smp, nCount = colSums(m), as.matrix(t(m[g, , drop = FALSE])), check.names = FALSE),
            gzfile(file.path(gradient_dir, paste0(smp, "_genesets.csv.gz"))), row.names = FALSE)
  ep[[smp]] <- data.frame(sample = smp, barcode = colnames(m), nCount = colSums(m),
                          EPCAM = as.vector(m["EPCAM", ]), KRT20 = as.vector(m["KRT20", ]))
  rm(so, m); gc()
}
write.csv(do.call(rbind, ep), gzfile(file.path(gradient_dir, "visium48_EPCAM_KRT20.csv.gz")), row.names = FALSE)
