#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S12c - EPCAM in YAP-Rev high tumour cells relative to other tumour cells at the same position
# (all positions, gland core, gland edge, bud), Visium HD 074_D1. EPCAM UMI per 10,000 UMI pooled per group; ratio
# YAP-Rev high / other; 95% CI by spatial block bootstrap (200 µm blocks, 1,000 resamples). All non-rejected tumour
# cells (RCTD first type for doublets).
# Inputs : <HD_DIR>/per_sample/Visium_HD_D1_1856_074_D1/objects (steps 01-06)
# Usage  : Rscript analyses/14_visium_hd_S12/08_FigS12c_EPCAM.R
# =====================================================================================================================
suppressPackageStartupMessages({ library(Matrix); library(ggplot2) })
hd_dir <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd"); set.seed(29)
S <- "Visium_HD_D1_1856_074_D1"; d074 <- file.path(hd_dir, "per_sample", S, "objects")
out <- file.path(hd_dir, "figures_S12"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
cnt <- readRDS(file.path(d074, paste0("counts_", S, ".rds"))); pos <- readRDS(file.path(d074, paste0("tumor_cells_position_", S, ".rds")))
k <- readRDS(file.path(d074, paste0("tumor_cells_yaprev_k15_min20_", S, ".rds"))); ca <- readRDS(file.path(d074, paste0("cells_annotated_", S, ".rds")))
tb <- pos$barcode; yh <- tb %in% k$barcode[k$yap_high]; ps <- as.character(pos$position)
sg <- ca$spot_class[match(tb, ca$barcode)] == "singlet"; tot <- colSums(cnt[, tb])
G <- c("EPCAM"); U <- sapply(G, function(g) cnt[g, tb]); U <- matrix(U, ncol = 1, dimnames = list(NULL, G))
blk <- paste(floor(pos$x_um / 200), floor(pos$y_um / 200)); ub <- unique(blk); bi <- split(seq_along(blk), blk)
POS <- c(all = "All positions", core = "Gland core", edge = "Gland edge", bud = "Bud (1-4 cells)")
rat <- function(idx, sel) { sapply(names(POS), function(p) { s <- sel[idx] & (if (p == "all") TRUE else ps[idx] == p)
  a <- idx[s & yh[idx]]; b <- idx[s & !yh[idx]]
  vapply(G, function(g) (sum(U[a, g]) / sum(tot[a])) / (sum(U[b, g]) / sum(tot[b])), numeric(1)) }) }
for (cl in c("singlets", "all_cells")) {
  sel <- if (cl == "singlets") sg else rep(TRUE, length(tb))
  obs <- matrix(rat(seq_along(tb), sel), nrow = 1, dimnames = list(G, names(POS)))
  bs <- array(replicate(1000, rat(unlist(bi[sample(ub, length(ub), TRUE)], use.names = FALSE), sel)), dim = c(1, length(POS), 1000), dimnames = list(G, names(POS), NULL))
  if (cl == "singlets") next   # singlet subset resampled first (fixed random stream), not reported
  d <- expand.grid(marker = G, position = names(POS), stringsAsFactors = FALSE)
  d$ratio <- obs[cbind(d$marker, d$position)]
  d$lo <- apply(bs, c(1, 2), quantile, .025, na.rm = TRUE)[cbind(d$marker, d$position)]
  d$hi <- apply(bs, c(1, 2), quantile, .975, na.rm = TRUE)[cbind(d$marker, d$position)]
  d$n_YAPRev_high <- sapply(d$position, function(p) sum(sel & yh & (if (p == "all") TRUE else ps == p)))
  d$n_other <- sapply(d$position, function(p) sum(sel & !yh & (if (p == "all") TRUE else ps == p)))
  write.csv(d, file.path(out, paste0("FigS12c_EPCAM_YAPRevhigh_vs_other_", cl, "_074_D1.csv")), row.names = FALSE); print(d, digits = 3)
  d$position <- factor(POS[d$position], levels = POS); d$marker <- factor(d$marker, levels = G)
  g <- ggplot(d, aes(position, ratio, fill = marker)) + geom_col(width = .6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, linewidth = .4) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = .3) +
    scale_fill_manual(values = c(EPCAM = "#2C7FB8"), guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    labs(x = NULL, y = "EPCAM in YAP-Rev-high tumor cells\n(relative to other tumor cells)",
         title = "YAP-Rev-high cells, 074_D1 (all cells)") +
    theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9), legend.position = "bottom")
  f <- paste0("FigS12c_EPCAM_YAPRevhigh_vs_other_", cl, "_074_D1")
  ggsave(file.path(out, paste0(f, ".pdf")), g, width = 3.6, height = 3.4, device = cairo_pdf)
}
