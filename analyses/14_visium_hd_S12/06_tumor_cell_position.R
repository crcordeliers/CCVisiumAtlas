#!/usr/bin/env Rscript
# =====================================================================================================================
# Visium HD, step 6 - position of each tumour cell (gland core / gland edge / bud) and YAP-Rev by position
# (Supplementary Fig. S12b; Methods "Visium HD single-cell analysis")
#
# Tumour cells (step-2 first_type "Malignant epithelial", non-rejected) are linked when <= LINK_D µm apart (20):
#  bud  = connected group of 1-4 cells; edge = larger group and >= 20% non-tumour cells among assigned neighbours
#  within 25 µm; core = otherwise. Also: distance to the nearest assigned non-tumour cell.
# YAP-Rev = AHNAK + ANXA1 + CAPN2 UMI per 10,000 UMI, pooled per class. Controls: CAF / myeloid / epithelial markers per
# class; RCTD singlets only; 200 expression-matched random 3-gene sets; spatial block bootstrap (200 µm) for 95% CIs;
# mixture model giving the increase expected from contamination by neighbouring CAF and myeloid cells only.
# Inputs : objects of steps 1-2
# Outputs: <HD_DIR>/per_sample/<SAMPLE>/objects/tumor_cells_position_<SAMPLE>.rds ; tables/yaprev_position_*.csv
#          (incl. yaprev_position_mixing_model.csv used in S12b)
# Usage  : SAMPLE=<id> Rscript analyses/14_visium_hd_S12/06_tumor_cell_position.R
# =====================================================================================================================
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Matrix)
  library(RANN)
  library(igraph)
  library(ggplot2)
})
LINK_D <- 20; BUD_MAX <- 4; NB_R <- 25; EDGE_FRAC <- 0.2
N_RAND <- 200; N_BOOT <- 1000; BLOCK <- 200
GENES <- c("AHNAK", "ANXA1", "CAPN2")
CTRL <- list(CAF = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
             Myeloid = c("C1QA", "C1QB", "TYROBP", "LYZ", "CD74"),
             Epithelial = c("EPCAM", "KRT8", "KRT18", "CEACAM5", "CDH17"))

sample_id <- Sys.getenv("SAMPLE"); stopifnot(nzchar(sample_id))
short <- sub("^Visium_HD2?[-_]", "", sample_id)
hd_dir <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd")
sroot  <- file.path(hd_dir, "per_sample", sample_id)
tabdir <- file.path(sroot, "tables")
dir.create(tabdir, showWarnings = FALSE, recursive = TRUE)

cells <- readRDS(file.path(sroot, "objects", paste0("cells_meta_", sample_id, ".rds")))
cnt <- readRDS(file.path(sroot, "objects", paste0("counts_", sample_id, ".rds")))
cd <- as.data.frame(colData(readRDS(file.path(sroot, "objects", paste0("rctd17_doublet_all_", sample_id, ".rds")))))
cells$cls <- cd$spot_class[match(cells$barcode, rownames(cd))]
cells$t17 <- as.character(cd$first_type[match(cells$barcode, rownames(cd))])
cells$t17[is.na(cells$cls) | cells$cls == "reject"] <- NA
assigned <- !is.na(cells$t17)
isT <- assigned & cells$t17 == "Malignant epithelial"
xy_all <- as.matrix(cells[, c("x_um", "y_um")])

# ---- connected groups of tumour cells ----
ti <- which(isT); xt <- xy_all[ti, ]
nb <- nn2(xt, xt, k = 30, searchtype = "radius", radius = LINK_D)
e <- cbind(rep(seq_along(ti), ncol(nb$nn.idx)), as.vector(nb$nn.idx))
e <- e[e[, 2] > 0 & e[, 1] < e[, 2], , drop = FALSE]
g <- make_empty_graph(n = length(ti), directed = FALSE) + edges(as.vector(t(e)))
comp <- components(g); csize <- comp$csize[comp$membership]

# ---- fraction of non-tumour neighbours (among assigned cells) ----
ai <- which(assigned); xa <- xy_all[ai, ]
nba <- nn2(xa, xt, k = 40, searchtype = "radius", radius = NB_R)
isT_a <- isT[ai]
frac_nt <- sapply(seq_along(ti), function(r) {
  j <- nba$nn.idx[r, ]; j <- j[j > 0]; j <- j[ai[j] != ti[r]]
  if (!length(j)) NA else mean(!isT_a[j])
})
# distance to the nearest assigned non-tumour cell
nti <- which(assigned & !isT)
dnt <- nn2(xy_all[nti, ], xt, k = 1)$nn.dists[, 1]

tum <- data.frame(barcode = cells$barcode[ti], x_um = xt[, 1], y_um = xt[, 2], cls = cells$cls[ti],
                  comp_size = csize, frac_nontumour = frac_nt, dist_nontumour = dnt)
tum$position <- ifelse(tum$comp_size <= BUD_MAX, "bud",
                       ifelse(!is.na(tum$frac_nontumour) & tum$frac_nontumour >= EDGE_FRAC, "edge", "core"))
tum$position <- factor(tum$position, levels = c("core", "edge", "bud"))
print(table(tum$position))

ct <- cnt[, tum$barcode]; tot <- colSums(ct)
tum$total_umi <- tot
tum$yap_umi <- colSums(ct[GENES, ])
rate <- function(num, idx) 1e4 * sum(num[idx]) / sum(tot[idx])

# ---- table per class: YAP-Rev, genes, controls ----
ctrl_genes <- intersect(unlist(CTRL), rownames(ct))
rows <- c("YAP-Rev (3 genes)", GENES, ctrl_genes)
mk_tab <- function(sel) {
  sapply(levels(tum$position), function(p) {
    idx <- sel & tum$position == p
    c(n_cells = sum(idx), rate(tum$yap_umi, idx), sapply(c(GENES, ctrl_genes), function(gn) rate(ct[gn, ], idx)))
  })
}
tab_all <- mk_tab(rep(TRUE, nrow(tum))); rownames(tab_all) <- c("n_cells", rows)
tab_sing <- mk_tab(tum$cls == "singlet"); rownames(tab_sing) <- c("n_cells", rows)
fold <- function(tb) cbind(tb, edge_vs_core = tb[, "edge"] / tb[, "core"], bud_vs_core = tb[, "bud"] / tb[, "core"])
out_all <- fold(tab_all); out_sing <- fold(tab_sing)
out_all[1, 4:5] <- NA; out_sing[1, 4:5] <- NA
message("all tumour cells (UMI per 10k)"); print(round(out_all, 2))
message("singlets only"); print(round(out_sing, 2))
write.csv(round(out_all, 3), file.path(tabdir, "yaprev_position_rates_all.csv"))
write.csv(round(out_sing, 3), file.path(tabdir, "yaprev_position_rates_singlets.csv"))

# ---- spatial block bootstrap ----
blk <- paste(floor(tum$x_um / BLOCK), floor(tum$y_um / BLOCK))
ub <- unique(blk); bidx <- split(seq_len(nrow(tum)), blk)
set.seed(42)
boot <- t(replicate(N_BOOT, {
  i <- unlist(bidx[sample(ub, replace = TRUE)], use.names = FALSE)
  p <- tum$position[i]; y <- tum$yap_umi[i]; tt <- tot[i]
  r <- tapply(seq_along(i), p, function(k) 1e4 * sum(y[k]) / sum(tt[k]))
  c(edge_vs_core = r[["edge"]] / r[["core"]], bud_vs_core = r[["bud"]] / r[["core"]])
}))
ci <- apply(boot, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
message("95% CI, block bootstrap (", BLOCK, " um): edge/core ", paste(round(ci[, 1], 2), collapse = "-"),
        "; bud/core ", paste(round(ci[, 2], 2), collapse = "-"))

# ---- control: random 3-gene sets matched on expression ----
gm <- rowMeans(ct)
cand <- names(gm)[gm > 0 & !(names(gm) %in% c(GENES, ctrl_genes))]
brk <- unique(quantile(gm[cand], seq(0, 1, length.out = 21)))
bin <- cut(gm, brk, include.lowest = TRUE); names(bin) <- names(gm)
idx_c <- tum$position == "core"; idx_e <- tum$position == "edge"; idx_b <- tum$position == "bud"
set.seed(7)
rand <- t(sapply(seq_len(N_RAND), function(r) {
  gs <- sapply(GENES, function(g) { pool <- cand[bin[cand] == bin[g]]; sample(pool, 1) })
  y <- colSums(ct[gs, , drop = FALSE])
  c(edge_vs_core = rate(y, idx_e) / rate(y, idx_c), bud_vs_core = rate(y, idx_b) / rate(y, idx_c))
}))
obs <- c(edge_vs_core = out_all["YAP-Rev (3 genes)", "edge_vs_core"], bud_vs_core = out_all["YAP-Rev (3 genes)", "bud_vs_core"])
p_rand <- sapply(names(obs), function(k) (sum(rand[, k] >= obs[k]) + 1) / (N_RAND + 1))
message("random gene sets: edge/core YAP-Rev = ", round(obs[1], 2), " vs random median ",
        round(median(rand[, 1]), 2), " (p = ", signif(p_rand[1], 3), "); bud/core YAP-Rev = ", round(obs[2], 2),
        " vs ", round(median(rand[, 2]), 2), " (p = ", signif(p_rand[2], 3), ")")
write.csv(data.frame(rand), file.path(tabdir, "yaprev_position_random_genesets.csv"), row.names = FALSE)

# ---- mixture model: increase expected from contamination alone ----
# Reference profiles = RCTD-singlet CAF (Fibroblast-CAF) and myeloid (Macrophage-Monocyte) cells of the section.
# Fraction of contaminating transcripts per class estimated from markers ((class - core) / (reference - core)); expected
# YAP-Rev = core + a_caf * (CAF - core) + a_mye * (myeloid - core).
ref_rate <- function(type, genes) {
  b <- cells$barcode[which(cells$cls == "singlet" & cells$t17 == type)]
  m <- cnt[intersect(genes, rownames(cnt)), b, drop = FALSE]
  1e4 * sum(m) / sum(cnt[, b])
}
caf_m <- intersect(CTRL$CAF, rownames(ct)); mye_m <- c("C1QA", "C1QB", "TYROBP")
mix <- do.call(rbind, lapply(list(all = rep(TRUE, nrow(tum)), singlets = tum$cls == "singlet"), function(sel) {
  cr <- function(genes, p) rate(colSums(ct[genes, , drop = FALSE]), sel & tum$position == p)
  out <- NULL
  for (p in c("edge", "bud")) {
    a_caf <- (cr(caf_m, p) - cr(caf_m, "core")) / (ref_rate("Fibroblast-CAF", caf_m) - cr(caf_m, "core"))
    a_mye <- (cr(mye_m, p) - cr(mye_m, "core")) / (ref_rate("Macrophage-Monocyte", mye_m) - cr(mye_m, "core"))
    for (gs in list(GENES, "AHNAK", "ANXA1", "CAPN2")) {
      core <- cr(gs, "core")
      expct <- core + a_caf * (ref_rate("Fibroblast-CAF", gs) - core) + a_mye * (ref_rate("Macrophage-Monocyte", gs) - core)
      out <- rbind(out, data.frame(position = p, genes = paste(gs, collapse = "+"), alpha_CAF = a_caf, alpha_myeloid = a_mye,
                                   core = core, observed = cr(gs, p), expected_contamination_only = expct,
                                   obs_vs_core = cr(gs, p) / core, expected_vs_core = expct / core))
    }
  }
  out
}))
mix$cells <- rep(c("all", "singlets"), each = nrow(mix) / 2)
message("mixture model (UMI per 10k)"); print(mix, digits = 3)
message("singlet CAF reference: YAP-Rev ", round(ref_rate("Fibroblast-CAF", GENES), 2), "; CAF markers ",
        round(ref_rate("Fibroblast-CAF", caf_m), 1), " | myeloid: YAP-Rev ", round(ref_rate("Macrophage-Monocyte", GENES), 2),
        "; C1QA/B+TYROBP ", round(ref_rate("Macrophage-Monocyte", mye_m), 1))
write.csv(mix, file.path(tabdir, "yaprev_position_mixing_model.csv"), row.names = FALSE)

# ---- strict control: genes matched on tumour expression AND on CAF/tumour enrichment ----
bcaf <- cells$barcode[which(cells$cls == "singlet" & cells$t17 == "Fibroblast-CAF")]
gcaf <- rowSums(cnt[, bcaf]) / sum(cnt[, bcaf])
gtum <- rowSums(ct[, idx_c]) / sum(tot[idx_c])
ok2 <- names(gtum)[gtum > 0 & !(names(gtum) %in% c(GENES, ctrl_genes))]
lr <- log2((gcaf[ok2] + 1e-6) / (gtum[ok2] + 1e-6))   # CAF vs tumour-core enrichment
brk_t <- unique(quantile(gtum[ok2], seq(0, 1, length.out = 11)))
brk_l <- unique(quantile(lr, seq(0, 1, length.out = 11)))
key <- paste(cut(gtum[ok2], brk_t, include.lowest = TRUE), cut(lr, brk_l, include.lowest = TRUE)); names(key) <- ok2
lr_y <- log2((gcaf[GENES] + 1e-6) / (gtum[GENES] + 1e-6))
key_y <- paste(cut(pmin(pmax(gtum[GENES], min(brk_t)), max(brk_t)), brk_t, include.lowest = TRUE),
               cut(pmin(pmax(lr_y, min(brk_l)), max(brk_l)), brk_l, include.lowest = TRUE))
message("YAP-Rev bins (tumour expression x CAF enrichment): ", paste(GENES, key_y, sep = "=", collapse = " | "),
        "; pool sizes: ", paste(sapply(key_y, function(k) sum(key == k)), collapse = "/"))
set.seed(11)
rand2 <- t(sapply(seq_len(N_RAND), function(r) {
  gs <- sapply(key_y, function(k) { pool <- ok2[key == k]; if (length(pool)) sample(pool, 1) else NA })
  if (anyNA(gs)) return(c(edge_vs_core = NA, bud_vs_core = NA))
  y <- colSums(ct[gs, , drop = FALSE])
  c(edge_vs_core = rate(y, idx_e) / rate(y, idx_c), bud_vs_core = rate(y, idx_b) / rate(y, idx_c))
}))
p_rand2 <- sapply(names(obs), function(k) (sum(rand2[, k] >= obs[k], na.rm = TRUE) + 1) / (sum(!is.na(rand2[, k])) + 1))
message("strict control (matched tumour x CAF): edge/core median ", round(median(rand2[, 1], na.rm = TRUE), 2),
        " (p = ", signif(p_rand2[1], 3), "); bud/core median ", round(median(rand2[, 2], na.rm = TRUE), 2),
        " (p = ", signif(p_rand2[2], 3), ")")
write.csv(data.frame(rand2), file.path(tabdir, "yaprev_position_random_genesets_tumour_x_CAF_matched.csv"), row.names = FALSE)

# ---- continuous gradient: distance to the nearest non-tumour cell ----
tum$dist_bin <- cut(tum$dist_nontumour, c(0, 10, 20, 30, 50, 80, Inf), right = FALSE)
grad <- do.call(rbind, lapply(levels(tum$dist_bin), function(b) {
  idx <- tum$dist_bin == b & !is.na(tum$dist_bin)
  data.frame(dist_bin = b, n_cells = sum(idx), YAPRev = rate(tum$yap_umi, idx),
             CAF = rate(colSums(ct[intersect(CTRL$CAF, rownames(ct)), ]), idx),
             Myeloid = rate(colSums(ct[intersect(CTRL$Myeloid, rownames(ct)), ]), idx),
             Epithelial = rate(colSums(ct[intersect(CTRL$Epithelial, rownames(ct)), ]), idx))
}))
print(grad)
write.csv(grad, file.path(tabdir, "yaprev_distance_gradient.csv"), row.names = FALSE)
saveRDS(tum, file.path(sroot, "objects", paste0("tumor_cells_position_", sample_id, ".rds")))

message("OK")
