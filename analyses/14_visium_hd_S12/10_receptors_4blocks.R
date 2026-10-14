#!/usr/bin/env Rscript
# =====================================================================================================================
# Collagen / matrix receptors in YAP-Rev high vs other tumour cells, 4 Visium HD sections (data of Supplementary Fig. S12f,g)
#  - YAP-Rev high: k15, FDR < 0.05 (05_*.R); positions core / edge / bud (06_*.R); RCTD singlets only.
#  - Ratio YAP-Rev high / other tumour cells, stratified by position (each stratum weighted by its number of YAP-Rev
#    high cells); 95% CI by spatial block bootstrap (200 µm blocks, 500 resamples).
#  - Contamination: fraction of CAF and myeloid transcripts estimated in each group from markers (CAF: COL1A1/COL1A2/
#    COL3A1/DCN/LUM; myeloid: C1QA/C1QB/TYROBP) and reference profiles (singlet CAFs and macrophages/monocytes of the
#    same section) -> ratio expected if the difference came only from contamination.
#  - Tumour specificity: expression in singlet CAFs / in tumour cells.
# Inputs : <HD_DIR>/per_sample/<sample>/objects for the 4 sections (steps 01-06)
# Output : <HD_DIR>/figures_S12/F9_receptors_YAPRevhigh_vs_other_4blocks.csv (genes detected in < 2% of the cells of
#          both groups flagged low_detection)
# Usage  : Rscript analyses/14_visium_hd_S12/10_receptors_4blocks.R
# =====================================================================================================================
suppressPackageStartupMessages({ library(Matrix); library(ggplot2) })
hd_dir <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd"); set.seed(31)
out <- file.path(hd_dir, "figures_S12"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
BL <- data.frame(hd = c("Visium_HD_D1_1856_074_D1", "Visium_HD2_A1_1341_395_C1", "Visium_HD_A1_1898_2177_B1", "Visium_HD2-D1-1656-295B1"),
                 lab = c("074_D1 (R)", "395_C1 (R)", "21-77_B1 (NR)", "295_B1 (NR)"), dir = NA, stringsAsFactors = FALSE)
BL$dir <- file.path(hd_dir, "per_sample", BL$hd, "objects")
COL6SET <- Sys.getenv("RECSET", "") == "COL6"
REC <- if (COL6SET) c("CSPG4", "ANTXR1", "ANTXR2", "ITGA1", "ITGA2", "ITGA3", "ITGA10", "ITGAV", "ITGB1", "ITGB3", "CD44", "SDC4") else
       c("CD44", "SDC1", "SDC4", "DDR1", "DDR2", "ITGA1", "ITGA2", "ITGA3", "ITGA5", "ITGA6", "ITGA10", "ITGA11", "ITGAV",
         "ITGB1", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB8", "CD47")
PCm <- c("RGS5", "PDGFRB", "NOTCH3", "KCNJ8")   # pericyte markers (contamination check for NG2/CSPG4)
PFX <- if (COL6SET) "F10_COL6receptors" else "F9_receptors"
CAFm <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"); MYm <- c("C1QA", "C1QB", "TYROBP")
res <- list()
for (i in seq_len(nrow(BL))) {
  s <- BL$hd[i]; ca <- readRDS(file.path(BL$dir[i], paste0("cells_annotated_", s, ".rds")))
  cnt <- readRDS(file.path(BL$dir[i], paste0("counts_", s, ".rds")))
  pos <- readRDS(file.path(BL$dir[i], paste0("tumor_cells_position_", s, ".rds")))
  k <- readRDS(file.path(BL$dir[i], paste0("tumor_cells_yaprev_k15_min20_", s, ".rds")))
  sg <- ca$spot_class[match(pos$barcode, ca$barcode)] == "singlet"
  pos <- pos[sg, ]; tb <- pos$barcode; yh <- tb %in% k$barcode[k$yap_high]; ps <- as.character(pos$position)
  g <- intersect(c(REC, CAFm, MYm, PCm), rownames(cnt)); M <- t(as.matrix(cnt[g, tb])); tot <- colSums(cnt[, tb])
  caf <- ca$barcode[which(ca$spot_class == "singlet" & ca$type17 == "Fibroblast-CAF")]
  my  <- ca$barcode[which(ca$spot_class == "singlet" & ca$type17 == "Macrophage-Monocyte")]
  pe  <- ca$barcode[which(ca$spot_class == "singlet" & ca$type17 == "Pericyte")]
  prof <- function(b, gg) 1e4 * rowSums(cnt[gg, b, drop = FALSE]) / sum(cnt[, b])
  pc <- prof(caf, g); pm <- prof(my, g); pp <- prof(pe, g)
  rate <- function(idx, gg) 1e4 * colSums(M[idx, gg, drop = FALSE]) / sum(tot[idx])
  strat_ratio <- function(idx) {     # ratio stratified by position, weighted by n YAP-Rev high
    num <- 0; w <- 0
    for (p in c("core", "edge", "bud")) { a <- idx[yh[idx] & ps[idx] == p]; b <- idx[!yh[idx] & ps[idx] == p]
      if (length(a) >= 10 && length(b) >= 10) { num <- num + length(a) * log((rate(a, g) + 1e-3) / (rate(b, g) + 1e-3)); w <- w + length(a) } }
    exp(num / w) }
  obs <- strat_ratio(seq_along(tb))
  blk <- paste(floor(pos$x_um / 200), floor(pos$y_um / 200)); ub <- unique(blk); bi <- split(seq_along(blk), blk)
  bs <- replicate(500, strat_ratio(unlist(bi[sample(ub, length(ub), TRUE)], use.names = FALSE)))
  # expected contamination (not stratified; order of magnitude), YAP-Rev high vs other
  a <- which(yh); b <- which(!yh); ra <- rate(a, g); rb <- rate(b, g)
  alpha <- function(r) c(caf = mean((r[CAFm] - rb[CAFm]) / (pc[CAFm] - rb[CAFm])), my = mean((r[MYm] - rb[MYm]) / (pm[MYm] - rb[MYm])),
                          pe = mean((r[PCm] - rb[PCm]) / (pp[PCm] - rb[PCm])))
  al <- alpha(ra)
  expd <- (rb + al["caf"] * (pc - rb) + al["my"] * (pm - rb) + max(0, al["pe"]) * (pp - rb)) / rb
  r <- data.frame(block = BL$lab[i], gene = intersect(REC, g), ratio = obs[intersect(REC, g)],
                  lo = apply(bs[intersect(REC, g), , drop = FALSE], 1, quantile, .025, na.rm = TRUE),
                  hi = apply(bs[intersect(REC, g), , drop = FALSE], 1, quantile, .975, na.rm = TRUE),
                  expected_contamination = expd[intersect(REC, g)],
                  pct_YAPhigh = 100 * colMeans(M[a, intersect(REC, g), drop = FALSE] > 0), pct_other = 100 * colMeans(M[b, intersect(REC, g), drop = FALSE] > 0),
                  cp10k_YAPhigh = ra[intersect(REC, g)], cp10k_other = rb[intersect(REC, g)],
                  CAF_over_tumor = (pc[intersect(REC, g)] + 1e-3) / (rb[intersect(REC, g)] + 1e-3),
                  Pericyte_over_tumor = (pp[intersect(REC, g)] + 1e-3) / (rb[intersect(REC, g)] + 1e-3), alpha_pericyte = al["pe"],
                  n_YAPhigh = length(a), n_other = length(b), alpha_CAF = al["caf"], alpha_myeloid = al["my"], row.names = NULL)
  res[[s]] <- r; message(s, " : n YAP-high singlets = ", length(a), " ; alpha CAF = ", round(al["caf"], 3), " myeloid = ", round(al["my"], 3))
}
R <- do.call(rbind, res)
R$low_detection <- pmax(R$pct_YAPhigh, R$pct_other) < 2   # detected in < 2% of the cells of both groups: ratio not interpretable
write.csv(R, file.path(out, paste0(PFX, "_YAPRevhigh_vs_other_4blocks.csv")), row.names = FALSE)
print(R[R$block == "074_D1 (R)", c("gene", "ratio", "lo", "hi", "expected_contamination", "pct_YAPhigh", "pct_other", "CAF_over_tumor")], digits = 3, row.names = FALSE)
W <- reshape(R[, c("block", "gene", "ratio")], idvar = "gene", timevar = "block", direction = "wide"); print(W, digits = 2, row.names = FALSE)

