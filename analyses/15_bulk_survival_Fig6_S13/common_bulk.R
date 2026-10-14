# =====================================================================================================================
# Shared by 01_Fig6_survival.R and 02_FigS13_forest_by_arm_MMR.R: bulk cohorts, signatures, ssGSEA scores.
# Signatures: top-20 marker genes of each meta-cluster and stromal sub-cluster (data/reference_TableS2_signatures.csv,
# = Supplementary Table S2, in rank order) and YAP-Rev (AHNAK, ANXA1, CAPN2). ssGSEA (GSVA, normalize = TRUE) on the
# union of the signature genes only, in the order of Table S2 (the gene order breaks ssGSEA ties).
# Creates PETACC.tot, IDEA.tot, PRODIGE.tot with z-scored scores (Metacluster_*, StromalSubcluster_*, YAPRev) and raw
# scores (*_raw), plus meta_names / stromal_names / colors_all for the displayed signatures.
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(GSVA); library(survival); library(dplyr); library(broom); library(ggplot2); library(ggpubr); library(ggtext)
  library(survminer)
})
out_dir <- file.path(figures.dir, "bulk"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
surv_data <- readRDS(Sys.getenv("CRCVISIUM_BULK", file.path(work.dir, "bulk_cohorts.rds")))
PETACC.tot  <- surv_data$PETACC
IDEA.tot    <- surv_data$IDEA
PRODIGE.tot <- surv_data$PRODIGE
PETACC.tot$etat_TTR <- PETACC.tot$recurr_cens
PETACC.tot$TTR_time <- PETACC.tot$recurr

# ---- signatures (Supplementary Table S2) ----
ts2 <- read.csv(file.path(data.dir, "reference_TableS2_signatures.csv"), stringsAsFactors = FALSE)
gene_list <- function(lvl, nm) { x <- ts2[ts2$level == lvl & ts2$signature == nm, ]; x$gene[order(x$rank)] }
mc_codes <- as.character(0:10); sc_codes <- as.character(0:4)          # numeric order of the cluster codes
ClusterListMeta <- setNames(lapply(mc_codes, function(k) gene_list("metacluster", meta_names[[k]])), paste0("Metacluster_", mc_codes))
ClusterListStromalSubcluster <- setNames(lapply(sc_codes, function(k) gene_list("stromal_subcluster", stromal_names[[k]])),
                                         paste0("StromalSubcluster_", sc_codes))
sign <- c(ClusterListMeta, ClusterListStromalSubcluster, list(YAPRev = c("AHNAK", "ANXA1", "CAPN2")))

# displayed signatures: tumour and stromal meta-clusters, stromal sub-clusters, YAP-Rev
meta_names_fp <- setNames(meta_names[c("0", "1", "7", "8", "2", "4")], paste0("Metacluster_", c("0", "1", "7", "8", "2", "4")))
stromal_names_fp <- setNames(stromal_names[stromal_order], paste0("StromalSubcluster_", stromal_order))
colors_all <- c(setNames(meta_colors[c("0", "1", "7", "8", "2", "4")], names(meta_names_fp)),
                setNames(stromal_colors[stromal_order], names(stromal_names_fp)), YAPRev = "#333333")
meta_names <- meta_names_fp; stromal_names <- stromal_names_fp

# ---- ssGSEA scores: z-scored (Fig. 6a-d) and raw (*_raw, Fig. 6e-g) ----
build_ssgsea <- function(dataset_name) {
  dataset <- get(paste0(dataset_name, ".tot"))
  expr <- t(as.matrix(dataset[, intersect(unique(unlist(sign, use.names = FALSE)), colnames(dataset)), drop = FALSE]))
  storage.mode(expr) <- "double"
  expr <- expr[rowSums(is.na(expr)) < ncol(expr), , drop = FALSE]
  if (anyNA(expr)) {                                  # row-median imputation
    med <- apply(expr, 1, function(x) median(x, na.rm = TRUE))
    idx <- which(is.na(expr), arr.ind = TRUE)
    expr[idx] <- med[idx[, 1]]
  }
  sign_mapped <- lapply(sign, function(g) intersect(g, rownames(expr)))
  sign_mapped <- sign_mapped[vapply(sign_mapped, length, integer(1)) >= 2]
  gs_scores <- gsva(ssgseaParam(exprData = expr, geneSets = sign_mapped, normalize = length(sign_mapped) > 1))
  scores_raw <- as.data.frame(t(gs_scores))
  scores_scaled <- as.data.frame(scale(scores_raw))
  colnames(scores_raw) <- paste0(colnames(scores_raw), "_raw")
  stopifnot(identical(rownames(scores_scaled), rownames(dataset)))
  dataset <- cbind(scores_scaled, scores_raw, dataset)
  colnames(dataset) <- make.names(colnames(dataset), unique = TRUE)
  assign(paste0(dataset_name, ".tot"), dataset, inherits = TRUE)
}
for (dataset_name in c("PETACC", "IDEA", "PRODIGE")) build_ssgsea(dataset_name)

