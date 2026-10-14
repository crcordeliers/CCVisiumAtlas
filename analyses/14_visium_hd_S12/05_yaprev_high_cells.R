#!/usr/bin/env Rscript
# =====================================================================================================================
# Visium HD, step 5 - YAP-Rev high tumour cells and their spatial clusters (Methods "Visium HD single-cell analysis")
#
# YAP-Rev (AHNAK + ANXA1 + CAPN2) scored in tumour cells only (step-2 first_type "Malignant epithelial", non-rejected):
#  1) neighbourhood = the K nearest tumour cells (cell included) within MAX_D µm (K = 15);
#  2) neighbourhood rate = sum of YAP-Rev UMI / sum of total UMI (pooled counts, robust to low depth);
#  3) null = N_PERM permutations of the cells (YAP-Rev/total pairs kept together, positions fixed) giving the per-cell
#     null mean and SD -> z-score, one-sided normal P, Benjamini-Hochberg FDR;
#  4) YAP-Rev high = FDR < 0.05; clusters = connected YAP-Rev high cells (link <= 30 µm) of >= MIN_AMAS cells (20).
# Inputs : objects of steps 1-2
# Outputs: <HD_DIR>/per_sample/<SAMPLE>/objects/tumor_cells_yaprev_k<K>_min<MIN_AMAS>_<SAMPLE>.rds ; tables/yaprev_amas_*.csv
# Usage  : SAMPLE=<id> Rscript analyses/14_visium_hd_S12/05_yaprev_high_cells.R 15 20 1000
# =====================================================================================================================
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Matrix)
  library(RANN)
  library(igraph)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1) as.integer(args[1]) else 15          # neighbourhood size (cell included)
MIN_AMAS <- if (length(args) >= 2) as.integer(args[2]) else 20   # minimum cluster size
N_PERM <- if (length(args) >= 3) as.integer(args[3]) else 1000
MAX_D <- 50; FDR_MAX <- 0.05; LINK_D <- 30
TAG <- paste0("k", K, "_min", MIN_AMAS)
GENES <- c("AHNAK", "ANXA1", "CAPN2")

sample_id <- Sys.getenv("SAMPLE"); stopifnot(nzchar(sample_id))
short <- sub("^Visium_HD2?[-_]", "", sample_id)
hd_dir <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd")
sroot  <- file.path(hd_dir, "per_sample", sample_id)
std_block <- c(Visium_HD_D1_1856_074_D1 = "Ech_074_D1", Visium_HD_A1_1898_2177_B1 = "Ech_21-77_B1",
               Visium_HD2_A1_1341_395_C1 = "Ech_395_C1", "Visium_HD2-D1-1656-295B1" = "Ech_295_B1")[sample_id]
for (d in c("tables", "tmp")) dir.create(file.path(sroot, d), showWarnings = FALSE, recursive = TRUE)
tabdir <- file.path(sroot, "tables")

cells <- readRDS(file.path(sroot, "objects", paste0("cells_meta_", sample_id, ".rds")))
cnt <- readRDS(file.path(sroot, "objects", paste0("counts_", sample_id, ".rds")))
lvl1 <- readRDS(file.path(sroot, "objects", paste0("rctd17_doublet_all_", sample_id, ".rds")))
cd1 <- as.data.frame(colData(lvl1))
cells$spot_class <- cd1$spot_class[match(cells$barcode, rownames(cd1))]
cells$type17 <- as.character(cd1$first_type[match(cells$barcode, rownames(cd1))])
cells$type17[is.na(cells$spot_class) | cells$spot_class == "reject"] <- NA

tum <- cells[which(cells$type17 == "Malignant epithelial"), ]
ct <- cnt[, tum$barcode]
tot <- colSums(ct)
n <- nrow(tum)
message("tumour cells: ", n)
stopifnot(all(GENES %in% rownames(ct)))

# ---- tumour kNN neighbourhood ----
xy <- as.matrix(tum[, c("x_um", "y_um")])
nn <- nn2(xy, xy, k = K)
ok <- nn$nn.dists <= MAX_D
W <- sparseMatrix(i = rep(seq_len(n), K)[as.vector(ok)], j = as.vector(nn$nn.idx)[as.vector(ok)],
                  x = 1, dims = c(n, n))
tot_nb <- as.vector(W %*% tot)

hotspots <- function(yap, n_perm, seed = 42) {
  # null by permutation of cells (positions fixed); as the smallest empirical p (1/(n_perm+1)) cannot pass BH over
  # ~43,000 tests, permutations are used to estimate the null mean and SD PER CELL -> z-score, one-sided normal p, BH FDR
  obs <- as.vector(W %*% yap) / tot_nb
  set.seed(seed)
  s1 <- numeric(n); s2 <- numeric(n)
  for (p in seq_len(n_perm)) {
    o <- sample(n)
    rp <- as.vector(W %*% yap[o]) / as.vector(W %*% tot[o])
    s1 <- s1 + rp; s2 <- s2 + rp^2
  }
  mu <- s1 / n_perm
  sdv <- sqrt(pmax(s2 / n_perm - mu^2, 1e-12))
  z <- (obs - mu) / sdv
  pval <- pnorm(z, lower.tail = FALSE)
  list(rate = obs, z = z, p = pval, fdr = p.adjust(pval, "BH"))
}

find_amas <- function(high) {
  idx <- which(high)
  if (length(idx) < MIN_AMAS) return(rep(NA_integer_, n))
  nb <- nn2(xy[idx, , drop = FALSE], xy[idx, , drop = FALSE], k = min(30, length(idx)),
            searchtype = "radius", radius = LINK_D)
  e <- cbind(rep(seq_along(idx), ncol(nb$nn.idx)), as.vector(nb$nn.idx))
  e <- e[e[, 2] > 0 & e[, 1] != e[, 2], , drop = FALSE]
  g <- graph_from_edgelist(e, directed = FALSE)
  g <- add_vertices(g, max(0, length(idx) - vcount(g)))
  comp <- components(g)$membership
  sz <- table(comp)
  keep <- as.integer(names(sz)[sz >= MIN_AMAS])
  out <- rep(NA_integer_, n)
  cm <- comp[seq_along(idx)]
  out[idx[cm %in% keep]] <- match(cm[cm %in% keep], keep)
  out
}

# ---- YAP-Rev ----
yap <- colSums(ct[GENES, , drop = FALSE])
h <- hotspots(yap, N_PERM)
tum$yap_umi <- yap; tum$total_umi <- tot
tum$yap_rate_nb <- h$rate * 1e4   # YAP-Rev UMI per 10,000 UMI (neighbourhood)
tum$yap_z <- h$z; tum$yap_p <- h$p; tum$yap_fdr <- h$fdr
tum$yap_high <- h$fdr < FDR_MAX
tum$amas <- find_amas(tum$yap_high)
n_high <- sum(tum$yap_high); n_in_amas <- sum(!is.na(tum$amas))
message("YAP-Rev high: ", n_high, " cells; in clusters (>= ", MIN_AMAS, "): ", n_in_amas,
        "; clusters: ", length(unique(na.omit(tum$amas))))

# contribution of each gene, inside / outside clusters (UMI per 10,000 UMI)
per_gene <- sapply(GENES, function(g) {
  x <- ct[g, ]
  c(in_amas = 1e4 * sum(x[!is.na(tum$amas)]) / sum(tot[!is.na(tum$amas)]),
    out_amas = 1e4 * sum(x[is.na(tum$amas)]) / sum(tot[is.na(tum$amas)]))
})
per_gene <- rbind(per_gene, fold = per_gene["in_amas", ] / per_gene["out_amas", ])
print(round(per_gene, 2))
write.csv(round(per_gene, 3), file.path(tabdir, paste0("yaprev_amas_per_gene_", TAG, ".csv")))

amas_tab <- do.call(rbind, lapply(split(seq_len(n), tum$amas), function(i) data.frame(
  amas = tum$amas[i[1]], n_cells = length(i), x_um = mean(tum$x_um[i]), y_um = mean(tum$y_um[i]),
  yap_per10k = 1e4 * sum(yap[i]) / sum(tot[i]),
  AHNAK_per10k = 1e4 * sum(ct["AHNAK", i]) / sum(tot[i]),
  ANXA1_per10k = 1e4 * sum(ct["ANXA1", i]) / sum(tot[i]),
  CAPN2_per10k = 1e4 * sum(ct["CAPN2", i]) / sum(tot[i]))))
if (!is.null(amas_tab)) {
  amas_tab <- amas_tab[order(-amas_tab$n_cells), ]
  write.csv(amas_tab, file.path(tabdir, paste0("yaprev_amas_list_", TAG, ".csv")), row.names = FALSE)
  print(head(amas_tab, 15))
}

saveRDS(tum, file.path(sroot, "objects", paste0("tumor_cells_yaprev_", TAG, "_", sample_id, ".rds")))

message("OK")
