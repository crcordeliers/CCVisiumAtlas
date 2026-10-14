#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S12e - number of matCAFs within 100 µm of YAP-Rev high cluster cells vs other tumour cells (074_D1)
#
# matCAF = cS26 myofibroblasts + cS27 CXCL14+ CAFs + cS28 GREM1+ CAFs (step-3 fine states, non-rejected cells).
# For each tumour cell: number of matCAFs (and of all stromal/immune cells) with centroid <= R µm (R = 50, 100, 200).
# Cluster-level test: for each YAP-Rev high cluster (k15_min20), median count around its cells; statistic = median over
# clusters; null = 10,000 draws of the same number of random tumour cells (one-sided P).
# Inputs : <HD_DIR>/per_sample/Visium_HD_D1_1856_074_D1/objects (steps 01-06)
# Outputs: <HD_DIR>/figures_S12/FigS12e_matCAF_within_100um_074_D1.pdf, matCAF_density_around_yaprev_074_D1.csv
# Usage  : Rscript analyses/14_visium_hd_S12/09_FigS12e_matCAF.R
# =====================================================================================================================
suppressPackageStartupMessages({ library(RANN); library(ggplot2) })
hd_dir <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd")
S <- "Visium_HD_D1_1856_074_D1"; objd <- file.path(hd_dir, "per_sample", S, "objects")
out <- file.path(hd_dir, "figures_S12"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(16)
ca  <- readRDS(file.path(objd, paste0("cells_annotated_", S, ".rds")))
pos <- readRDS(file.path(objd, paste0("tumor_cells_position_", S, ".rds")))
k15 <- readRDS(file.path(objd, paste0("tumor_cells_yaprev_k15_min20_", S, ".rds")))
MATT <- c("cS26 (Myofibro)", "cS27 (CXCL14+ CAF)", "cS28 (GREM1+ CAF)")
ok <- which(!is.na(ca$spot_class) & ca$spot_class != "reject")
mc <- ca[ok, ][!is.na(ca$type98[ok]) & ca$type98[ok] %in% MATT, ]
st <- ca[ok, ][!is.na(ca$type17[ok]) & !ca$type17[ok] %in% c("Malignant epithelial", "Normal epithelial"), ]
xy <- as.matrix(pos[, c("x_um", "y_um")])
cnt <- function(ref, R) rowSums(nn2(as.matrix(ref[, c("x_um", "y_um")]), xy, k = min(nrow(ref), 2000), searchtype = "radius", radius = R)$nn.idx > 0)
RS <- c(50, 100, 200)
D <- sapply(RS, function(R) cnt(mc, R)); colnames(D) <- paste0("matCAF_", RS, "um")
DS <- sapply(RS, function(R) cnt(st, R)); colnames(DS) <- paste0("stroma_", RS, "um")
A <- cbind(data.frame(barcode = pos$barcode, amas_k15 = k15$amas[match(pos$barcode, k15$barcode)]), D, DS)

res <- list(); NULLN <- 10000
for (R in RS) {
  v <- A[[paste0("matCAF_", R, "um")]]; vs <- A[[paste0("stroma_", R, "um")]]; inG <- !is.na(A$amas_k15)
  am <- tapply(v[inG], A$amas_k15[inG], median); obs <- median(am)
  nul <- replicate(NULLN, median(sample(v, length(am))))
  res[[as.character(R)]] <- data.frame(group = "amas_k15", radius_um = R, n_cells = sum(inG), n_clusters = length(am),
    median_matCAF_group = median(v[inG]), median_matCAF_other_tumour = median(v[!inG]),
    pct_matCAF_among_stroma_group = 100 * sum(v[inG]) / sum(vs[inG]), pct_matCAF_among_stroma_other = 100 * sum(v[!inG]) / sum(vs[!inG]),
    median_over_clusters = obs, null_median = median(nul), pct_clusters_above_tumour_median = 100 * mean(am > median(v)),
    P_cluster_level = (sum(nul >= obs) + 1) / (NULLN + 1))
}
T <- do.call(rbind, res); rownames(T) <- NULL
write.csv(T, file.path(out, "matCAF_density_around_yaprev_074_D1.csv"), row.names = FALSE)
print(T, digits = 3)

# ---- figure: matCAF within 100 um ----
lab <- ifelse(!is.na(A$amas_k15), "YAP-Rev-high\ncluster cells", "Other tumor\ncells")
P <- data.frame(g = factor(lab, levels = c("Other tumor\ncells", "YAP-Rev-high\ncluster cells")), v = A$matCAF_100um)
am <- tapply(A$matCAF_100um[!is.na(A$amas_k15)], A$amas_k15[!is.na(A$amas_k15)], median)
pp <- T$P_cluster_level[T$group == "amas_k15" & T$radius_um == 100]
# boxplots without outliers (most other tumour cells have zero matCAF); dots = median per cluster (unit of the test)
set.seed(16)
pts <- data.frame(g = factor("YAP-Rev-high\ncluster cells", levels = levels(P$g)), v = am)
top <- max(quantile(P$v[P$g == levels(P$g)[2]], .75) + 1.5 * IQR(P$v[P$g == levels(P$g)[2]]), max(am)) * 1.08
g <- ggplot(P, aes(g, v, fill = g)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = .75, colour = "grey25", coef = 1.5) +
  geom_jitter(data = pts, aes(g, v), inherit.aes = FALSE, width = .12, height = 0, size = 1.6, shape = 21, fill = "#E6AB02", colour = "black", stroke = .3) +
  annotate("segment", x = 1, xend = 2, y = top, yend = top, linewidth = .4) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = top, yend = top * .97, linewidth = .4) +
  annotate("text", x = 1.5, y = top * 1.035, label = paste0("cluster-level P ", ifelse(pp < 1e-3, "< 0.001", paste0("= ", signif(pp, 2)))), size = 3, vjust = 0) +
  scale_fill_manual(values = c("grey75", "#E6AB02"), guide = "none") +
  scale_y_continuous(limits = c(0, top * 1.1), expand = expansion(mult = c(0.02, 0))) +
  labs(x = NULL, y = "matCAF within 100 \u00b5m (cells)", title = "matCAF around YAP-Rev-high tumor cells (074_D1)") +
  theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9))
f1 <- "FigS12e_matCAF_within_100um_074_D1"
ggsave(file.path(out, paste0(f1, ".pdf")), g, width = 3.4, height = 3.6, device = cairo_pdf)
r <- T[T$radius_um == 100, ]

