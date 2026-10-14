#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S12a, S12b, S12d - Visium HD section 074_D1 (segmented cells)
#   S12a  map of cell types (RCTD) with YAP-Rev high tumour cells in red          -> FigS12a_celltype_map_YAPRev_074_D1.pdf
#   S12b  YAP-Rev expression of tumour cells at gland edges and in buds relative to gland-core cells; 95% CI by spatial
#         block bootstrap (200 µm blocks, 1,000 resamples); grey = level expected from contamination by neighbouring
#         CAF/myeloid cells (mixture model of 06_tumor_cell_position.R)          -> FigS12b_YAPRev_by_tumor_position_074_D1.pdf
#   S12d  % stromal/immune neighbours (<= 30 µm) of YAP-Rev high clusters vs 2,000 size-matched random groups of tumour
#         cells, cluster-level test                                               -> FigS12d_stromal_neighbors_YAPRev_clusters_074_D1.pdf
# Inputs : <HD_DIR>/per_sample/Visium_HD_D1_1856_074_D1/{objects,tables} (steps 01-06)
# Outputs: <HD_DIR>/figures_S12/ (one PDF per panel)
# Usage  : Rscript analyses/14_visium_hd_S12/07_FigS12a_b_d.R
# =====================================================================================================================
suppressPackageStartupMessages({ library(Matrix); library(RANN); library(ggplot2) })
hd_dir <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd")
S <- "Visium_HD_D1_1856_074_D1"
d074 <- file.path(hd_dir, "per_sample", S, "objects"); tab074 <- file.path(hd_dir, "per_sample", S, "tables")
out <- file.path(hd_dir, "figures_S12"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(23)
YAP <- c("AHNAK", "ANXA1", "CAPN2"); RED <- "#D7301F"; POSL <- c(core = "Gland core", edge = "Gland edge", bud = "Bud (1-4 cells)")
th <- theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9))

ca  <- readRDS(file.path(d074, paste0("cells_annotated_", S, ".rds")))
cnt <- readRDS(file.path(d074, paste0("counts_", S, ".rds")))
pos <- readRDS(file.path(d074, paste0("tumor_cells_position_", S, ".rds")))
k15 <- readRDS(file.path(d074, paste0("tumor_cells_yaprev_k15_min20_", S, ".rds")))
mix <- read.csv(file.path(tab074, "yaprev_position_mixing_model.csv"))

# ---------------- S12b: ratios by position + spatial block bootstrap ----------------
tb <- pos$barcode; tot <- colSums(cnt[, tb])
G <- list(`YAP-Rev` = YAP, EPCAM = "EPCAM", CEACAM5 = "CEACAM5")
U <- sapply(G, function(g) colSums(cnt[g, tb, drop = FALSE]))
blk <- paste(floor(pos$x_um / 200), floor(pos$y_um / 200)); ub <- unique(blk); bi <- split(seq_along(blk), blk)
ratio <- function(idx) { p <- pos$position[idx]
  r <- sapply(names(G), function(g) sapply(c("core", "edge", "bud"), function(q) sum(U[idx[p == q], g]) / sum(tot[idx[p == q]])))
  sweep(r, 2, r["core", ], "/") }
obs <- ratio(seq_along(tb))
boot <- replicate(1000, ratio(unlist(bi[sample(ub, length(ub), replace = TRUE)], use.names = FALSE)))
R <- do.call(rbind, lapply(names(G), function(g) data.frame(marker = g, position = c("core", "edge", "bud"), ratio = obs[, g],
  lo = apply(boot[, g, ], 1, quantile, .025), hi = apply(boot[, g, ], 1, quantile, .975),
  cp10k = 1e4 * sapply(c("core", "edge", "bud"), function(q) sum(U[pos$position == q, g]) / sum(tot[pos$position == q])))))
R$n_cells <- as.vector(table(pos$position)[R$position])
write.csv(R, file.path(out, "FigS12b_position_ratios_074_D1.csv"), row.names = FALSE); print(R, digits = 3)

m <- mix[mix$genes == "AHNAK+ANXA1+CAPN2" & mix$cells == "all", ]
d1 <- rbind(data.frame(position = R$position[R$marker == "YAP-Rev"], what = "Observed", v = R$ratio[R$marker == "YAP-Rev"],
                       lo = R$lo[R$marker == "YAP-Rev"], hi = R$hi[R$marker == "YAP-Rev"]),
            data.frame(position = c("core", "edge", "bud"), what = "Expected from contamination only",
                       v = c(1, m$expected_vs_core[m$position == "edge"], m$expected_vs_core[m$position == "bud"]), lo = NA, hi = NA))
d1$position <- factor(POSL[d1$position], levels = POSL); d1$what <- factor(d1$what, levels = c("Observed", "Expected from contamination only"))
g1 <- ggplot(d1, aes(position, v, fill = what)) + geom_col(position = position_dodge(.75), width = .7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.75), width = .2, linewidth = .4, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = .3) +
  scale_fill_manual(values = c("Observed" = RED, "Expected from contamination only" = "grey75"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = NULL, y = "YAP-Rev expression\n(relative to gland core)", title = "YAP-Rev in tumor cells by position (074_D1)") + th + theme(legend.position = "bottom")
ggsave(file.path(out, "FigS12b_YAPRev_by_tumor_position_074_D1.pdf"), g1, width = 3.8, height = 3.6, device = cairo_pdf)


# ---------------- S12d: stromal neighbours, cluster level ----------------
st <- ca[which(ca$spot_class == "singlet" & !is.na(ca$type17) & !ca$type17 %in% c("Malignant epithelial", "Normal epithelial")), ]
ep <- ca[which(ca$spot_class == "singlet" & ca$type17 %in% c("Malignant epithelial", "Normal epithelial")), ]
inam <- !is.na(k15$amas[match(tb, k15$barcode)]); amb <- tb[inam]
ep2 <- ep[!ep$barcode %in% amb, ]                         # cells of the clusters excluded from neighbours
xy <- as.matrix(pos[, c("x_um", "y_um")])
ns <- rowSums(nn2(as.matrix(st[, c("x_um", "y_um")]), xy, k = 60, searchtype = "radius", radius = 30)$nn.idx > 0)
ne <- rowSums(nn2(as.matrix(ep2[, c("x_um", "y_um")]), xy, k = 60, searchtype = "radius", radius = 30)$nn.idx > 0) - (tb %in% ep2$barcode)
grp_pct <- function(i) 100 * sum(ns[i]) / max(1, sum(ns[i]) + sum(ne[i]))
am_id <- k15$amas[match(tb, k15$barcode)]
obs_am <- tapply(seq_along(tb), am_id, grp_pct); sz <- as.vector(table(am_id[inam]))
message("pooled % stromal neighbours: clusters ", round(grp_pct(which(inam)), 1), " vs other ", round(grp_pct(which(!inam)), 1))
# random groups: the n nearest tumour cells (outside clusters) of a random tumour cell; sizes drawn from the cluster sizes
oth <- which(!inam); kn <- nn2(xy[oth, ], xy[oth, ], k = max(sz))$nn.idx
rand_grp <- function(n) { c0 <- sample(length(oth), 1); grp_pct(oth[kn[c0, 1:n]]) }
null <- replicate(2000, rand_grp(sample(sz, 1)))
null_med <- replicate(10000, median(sample(null, length(obs_am), replace = TRUE)))
p3 <- (sum(null_med >= median(obs_am)) + 1) / (length(null_med) + 1)
d3 <- rbind(data.frame(g = "Random tumor-cell\ngroups (size-matched)", v = null), data.frame(g = "YAP-Rev-high\nclusters", v = as.vector(obs_am)))
d3$g <- factor(d3$g, levels = c("Random tumor-cell\ngroups (size-matched)", "YAP-Rev-high\nclusters"))
top <- max(d3$v) * 1.08
g3 <- ggplot(d3, aes(g, v, fill = g)) + geom_boxplot(width = .5, outlier.shape = NA, alpha = .75, colour = "grey25") +
  geom_jitter(data = d3[d3$g == levels(d3$g)[2], ], width = .12, height = 0, size = 1.6, shape = 21, fill = RED, colour = "black", stroke = .3) +
  annotate("segment", x = 1, xend = 2, y = top, yend = top, linewidth = .4) + annotate("segment", x = c(1, 2), xend = c(1, 2), y = top, yend = top * .97, linewidth = .4) +
  annotate("text", x = 1.5, y = top * 1.035, label = paste0("cluster-level P ", ifelse(p3 < 1e-3, "< 0.001", paste0("= ", signif(p3, 2)))), size = 3, vjust = 0) +
  scale_fill_manual(values = c("grey75", RED), guide = "none") + scale_y_continuous(limits = c(0, top * 1.1), expand = expansion(mult = c(.02, 0))) +
  labs(x = NULL, y = "Stromal/immune neighbors (%)", title = "Stromal contact of YAP-Rev-high clusters (074_D1)") + th
ggsave(file.path(out, "FigS12d_stromal_neighbors_YAPRev_clusters_074_D1.pdf"), g3, width = 3.4, height = 3.6, device = cairo_pdf)
message("S12d: median clusters ", round(median(obs_am), 1), " vs random groups ", round(median(null), 1), " P = ", signif(p3, 3))

# ---------------- S12a: cell-type map ----------------
cc <- ca[which(!is.na(ca$spot_class) & ca$spot_class != "reject" & !is.na(ca$type17)), ]
yh <- k15$barcode[k15$yap_high]
MAT <- c("cS26 (Myofibro)", "cS27 (CXCL14+ CAF)", "cS28 (GREM1+ CAF)")
cc$cat <- with(cc, ifelse(barcode %in% yh, "YAP-Rev-high tumor cells",
  ifelse(type17 == "Malignant epithelial", "Tumor cells",
  ifelse(type17 == "Fibroblast-CAF" & !is.na(type98) & type98 %in% MAT, "matCAF",
  ifelse(type17 == "Fibroblast-CAF", "Other fibroblasts/CAF",
  ifelse(type17 %in% c("Pericyte", "Smooth muscle"), "Pericytes/smooth muscle",
  ifelse(type17 == "Endothelial", "Endothelial",
  ifelse(type17 %in% c("Macrophage-Monocyte", "Dendritic cell", "Granulocyte", "Mast"), "Myeloid",
  ifelse(type17 %in% c("B", "Plasma", "CD4 T", "CD8 T", "NK", "T-NK-ILC-Other"), "Lymphoid/plasma",
  ifelse(type17 == "Normal epithelial", "Normal epithelium", "Other"))))))))))
CAT <- c("Tumor cells" = "#2171B5", "matCAF" = "#006D2C", "Other fibroblasts/CAF" = "#41AB5D", "Pericytes/smooth muscle" = "#A1D99B",
         "Endothelial" = "#807DBA", "Myeloid" = "#FD8D3C", "Lymphoid/plasma" = "#E7BA00", "Normal epithelium" = "#969696", "Other" = "#D9D9D9",
         "YAP-Rev-high tumor cells" = RED)
cc$cat <- factor(cc$cat, levels = names(CAT))
bg <- cc[cc$cat != "YAP-Rev-high tumor cells", ]; fg <- cc[cc$cat == "YAP-Rev-high tumor cells", ]
H <- diff(range(cc$y_um)) / diff(range(cc$x_um))
mk5 <- function(d1, d2, ps, title, xl = NULL, yl = NULL) {
  g <- ggplot() + geom_point(data = d1, aes(x_um, -y_um, colour = cat), size = ps, stroke = 0, shape = 16) +
    geom_point(data = d2, aes(x_um, -y_um, colour = cat), size = ps * 2.5, stroke = 0, shape = 16) +
    scale_colour_manual(values = CAT, name = NULL, drop = FALSE) + guides(colour = guide_legend(override.aes = list(size = 2.5))) +
    labs(title = title) + theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 9, hjust = .5), legend.position = "right", plot.background = element_rect(fill = "white", colour = NA))
  if (is.null(xl)) g + coord_equal() else g + coord_equal(xlim = xl, ylim = -rev(yl), expand = FALSE) +
    annotate("segment", x = xl[2] - 250, xend = xl[2] - 50, y = -yl[2] + 40, yend = -yl[2] + 40, linewidth = .8) +
    annotate("text", x = xl[2] - 150, y = -yl[2] + 80, label = "200 µm", size = 2.6) }
g5 <- mk5(bg, fg, .45, "Visium HD 074_D1: cell types (RCTD) and YAP-Rev-high tumor cells") +
  annotate("segment", x = max(cc$x_um) - 1100, xend = max(cc$x_um) - 100, y = -max(cc$y_um) - 150, yend = -max(cc$y_um) - 150, linewidth = .8) +
  annotate("text", x = max(cc$x_um) - 600, y = -max(cc$y_um) - 380, label = "1 mm", size = 2.8)
ggsave(file.path(out, "FigS12a_celltype_map_YAPRev_074_D1.pdf"), g5, width = 11, height = 9.5 * H + .3, device = cairo_pdf)
