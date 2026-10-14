#!/usr/bin/env Rscript
# =====================================================================================================================
# Visium HD, step 4 - per-cell annotation table and cell-type neighbourhood enrichment (Methods "Visium HD single-cell")
#
# Annotation: type17 = step-2 first_type (non-rejected cells); type98 = step-3 fine state ("unresolved" -> NA).
# Neighbourhood enrichment: pairs of cells whose centroids are <= RADIUS µm apart (default 30); observed number of pairs
# per type pair vs N_PERM random permutations of the labels (positions fixed) -> z-score and log2(observed/expected).
# Computed on all non-rejected cells and on RCTD singlets only, at both levels.
# Inputs : objects of steps 1-3
# Outputs: <HD_DIR>/per_sample/<SAMPLE>/objects/cells_annotated_<SAMPLE>.rds ; tables/nhood_enrichment_r<RADIUS>.csv
# Usage  : SAMPLE=<id> Rscript analyses/14_visium_hd_S12/04_annotate_and_proximity.R [RADIUS]
# =====================================================================================================================
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Matrix)
  library(RANN)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
RADIUS <- if (length(args) >= 1) as.numeric(args[1]) else 30
N_PERM <- 1000
MIN_CELLS <- 50

sample_id <- Sys.getenv("SAMPLE"); stopifnot(nzchar(sample_id))
short <- sub("^Visium_HD2?[-_]", "", sample_id)
hd_dir <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd")
sroot  <- file.path(hd_dir, "per_sample", sample_id)
std_block <- c(Visium_HD_D1_1856_074_D1 = "Ech_074_D1", Visium_HD_A1_1898_2177_B1 = "Ech_21-77_B1",
               Visium_HD2_A1_1341_395_C1 = "Ech_395_C1", "Visium_HD2-D1-1656-295B1" = "Ech_295_B1")[sample_id]
for (d in c("tables", "tmp")) dir.create(file.path(sroot, d), showWarnings = FALSE, recursive = TRUE)
tabdir <- file.path(sroot, "tables")

cells <- readRDS(file.path(sroot, "objects", paste0("cells_meta_", sample_id, ".rds")))
lvl1 <- readRDS(file.path(sroot, "objects", paste0("rctd17_doublet_all_", sample_id, ".rds")))
cd1 <- as.data.frame(colData(lvl1))
cells$spot_class <- cd1$spot_class[match(cells$barcode, rownames(cd1))]
cells$type17 <- as.character(cd1$first_type[match(cells$barcode, rownames(cd1))])
cells$spot_class[is.na(cells$spot_class)] <- "below_UMI_min"
cells$type17[cells$spot_class %in% c("reject", "below_UMI_min")] <- NA

sub_files <- list.files(file.path(sroot, "objects", "rctd98_within_family"), "^subtypes_.*\\.rds$", full.names = TRUE)
if (length(sub_files)) {
  subs <- do.call(rbind, lapply(sub_files, readRDS))
  cells$type98 <- subs$subtype98[match(cells$barcode, subs$barcode)]
  cells$type98[grepl("\\(unresolved\\)$", cells$type98)] <- NA
} else cells$type98 <- NA
saveRDS(cells, file.path(sroot, "objects", paste0("cells_annotated_", sample_id, ".rds")))
write.csv(as.data.frame(table(type17 = cells$type17, useNA = "ifany")), file.path(tabdir, "counts_type17.csv"), row.names = FALSE)
write.csv(as.data.frame(table(type98 = cells$type98, useNA = "ifany")), file.path(tabdir, "counts_type98.csv"), row.names = FALSE)

nhood_enrichment <- function(xy, lab, radius, n_perm) {
  n <- nrow(xy)
  nn <- nn2(xy, xy, k = 60, searchtype = "radius", radius = radius)
  i <- rep(seq_len(n), ncol(nn$nn.idx)); j <- as.vector(nn$nn.idx)
  ok <- j > 0 & j != i
  A <- sparseMatrix(i = i[ok], j = j[ok], x = 1, dims = c(n, n))
  lab <- factor(lab)
  L <- sparseMatrix(i = seq_len(n), j = as.integer(lab), x = 1, dims = c(n, nlevels(lab)))
  obs <- as.matrix(t(L) %*% A %*% L)
  set.seed(42)
  sims <- array(NA_real_, c(nlevels(lab), nlevels(lab), n_perm))
  for (p in seq_len(n_perm)) {
    Lp <- L[sample(n), ]
    sims[, , p] <- as.matrix(t(Lp) %*% A %*% Lp)
  }
  mu <- apply(sims, c(1, 2), mean); sdv <- apply(sims, c(1, 2), sd)
  z <- (obs - mu) / sdv
  dimnames(obs) <- dimnames(z) <- dimnames(mu) <- list(levels(lab), levels(lab))
  list(obs = obs, exp = mu, z = z, log2oe = log2((obs + 1) / (mu + 1)))
}

to_long <- function(r, level, subset) {
  d <- expand.grid(type_a = rownames(r$z), type_b = colnames(r$z), stringsAsFactors = FALSE)
  d$observed <- r$obs[cbind(d$type_a, d$type_b)]
  d$expected <- r$exp[cbind(d$type_a, d$type_b)]
  d$z <- r$z[cbind(d$type_a, d$type_b)]
  d$log2_obs_exp <- r$log2oe[cbind(d$type_a, d$type_b)]
  d$level <- level; d$cells <- subset; d$radius_um <- RADIUS
  d
}

run_one <- function(level, lab_col, subset_name, keep) {
  lab <- cells[[lab_col]]
  keep <- keep & !is.na(lab)
  big <- names(which(table(lab[keep]) >= MIN_CELLS))
  keep <- keep & lab %in% big
  message("level ", level, " / ", subset_name, ": ", sum(keep), " cells, ", length(big), " types")
  if (length(big) < 2) return(NULL)
  r <- nhood_enrichment(as.matrix(cells[keep, c("x_um", "y_um")]), lab[keep], RADIUS, N_PERM)
  to_long(r, level, subset_name)
}

res <- list(
  run_one("17", "type17", "all_nonreject", rep(TRUE, nrow(cells))),
  run_one("17", "type17", "singlet_only", cells$spot_class == "singlet"),
  run_one("98", "type98", "all_nonreject", rep(TRUE, nrow(cells))),
  run_one("98", "type98", "singlet_only", cells$spot_class == "singlet")
)
out <- do.call(rbind, res)
write.csv(out, file.path(tabdir, paste0("nhood_enrichment_r", RADIUS, ".csv")), row.names = FALSE)
message("OK")
