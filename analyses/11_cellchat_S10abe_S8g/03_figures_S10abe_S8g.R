#!/usr/bin/env Rscript
# =====================================================================================================================
# Cell-cell communication figures (spatial CellChat, patient level) - Supplementary Fig. S10a, S10b, S10e and S8g
#
#   Section A  S10e: collagen VI communication, Stroma-A1..B3 -> YAP-Rev high
#   Section B  S10a: communication strength, YAP-Rev high <-> stromal sub-clusters (both directions pooled)
#   Section C  S10b: ligand -> receptor diagrams, Stroma-A2 / Stroma-B2 <-> YAP-Rev high (both directions pooled)
#   Section D  S8g: YAP-Rev high spots per patient (% of all spots), Relapse vs No Relapse
#
# Communication strength between two groups in a patient = sum of the probabilities of the interactions significant in
# that patient (CellChat permutation P < 0.05); groups with < 10 spots are excluded. Relapse vs No Relapse: two-sided
# Mann-Whitney (exact if one group <= 8 and no ties, as scipy.stats.mannwhitneyu(method = "auto")); BH within each
# YAP-Rev high - stromal sub-cluster pair at the ligand-receptor level.
# Inputs : <CELLCHAT_DIR>/tables/{LR_prob_per_sample_nonzero.csv.gz, n_spots_per_group_per_sample.csv}
#          (02_extract_per_patient_LR.R); data/data_clin_48_filtered.RDS
# Output : <CELLCHAT_DIR>/figures/
# Usage  : Rscript analyses/11_cellchat_S10abe_S8g/03_figures_S10abe_S8g.R
# =====================================================================================================================
source("R/00_config.R")
# ---------------------------------------------------------------- configuration
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggforce); library(patchwork); library(RColorBrewer)
})

cellchat_dir <- Sys.getenv("CRCVISIUM_CELLCHAT", file.path(work.dir, "CellChat"))
base_dir <- file.path(cellchat_dir, "tables")
out_dir  <- file.path(cellchat_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

lr_file    <- file.path(base_dir, "LR_prob_per_sample_nonzero.csv.gz")
nspot_file <- file.path(base_dir, "n_spots_per_group_per_sample.csv")

YAP  <- "Regenerative"                   # CellChat label of YAP-Rev high spots (stromal spots, YAPRev_raw >= 0.40)
MIN_SPOTS <- 10                          # a group is considered in a patient only with >= 10 spots
COL6 <- c("COL6A1", "COL6A2", "COL6A3")
# immune-related CellChat pathways
IMM  <- c("MIF", "CXCL", "CCL", "MHC-I", "MHC-II", "COMPLEMENT", "ICAM", "VCAM", "ANNEXIN", "SPP1", "CD86", "IL16", "CD6",
          "CD45", "SELPLG", "BAFF", "LT", "TNF", "IFN-II", "CSF", "IL1", "IL6")
LAB  <- c("Stroma A 1" = "Stroma-A1", "Stroma A 2" = "Stroma-A2", "Stroma B 1" = "Stroma-B1",
          "Stroma B 2" = "Stroma-B2", "Stroma B 3" = "Stroma-B3")
NODE <- c("Stroma A 2" = "#B0EACD", "Stroma B 2" = "#F4A9A8")   # official stromal subcluster colours
GRP_COL <- c("No Relapse" = "#9197FF", "Relapse" = "#FC7171")
SPECTRAL_R <- rev(brewer.pal(11, "Spectral"))
PT <- ggplot2::.pt

# partner label -> CellChat groups it pools (canonical metacluster order; Metacluster_8 = Tumor-fibroblasts)
META <- list("Tumor-cycling" = "Tumor-cycling", "Tumor-intermediate" = "Tumor-intermediate",
             "Tumor-infiltrated" = "Tumor-infiltrated", "Tumor-fibroblasts" = "Tumor-fibroblasts",
             "Stroma-A" = c("Stroma A 1", "Stroma A 2"), "Stroma-B" = c("Stroma B 1", "Stroma B 2", "Stroma B 3"),
             "Necrosis" = "Necrosis", "Muscularis propria" = "Muscularis propria", "Lymphoid aggregates" = "Lymphoid aggregates",
             "Nerves - arteries" = "Nerves - arteries", "Normal colon epithelium" = "Normal colon epithelium")
SUBS <- setNames(as.list(names(LAB)), LAB)

# ---- data loading ----
d    <- fread(lr_file)[pval < 0.05]     # interactions significant within each sample
n    <- fread(nspot_file)
clin <- as.data.table(readRDS(file.path(data.dir, "data_clin_48_filtered.RDS")))
grp  <- setNames(ifelse(clin$recurr_cens == 1, "Relapse", "No Relapse"), clin$sample)
NSP  <- setNames(n$n_spots, paste(n$sample, n$group, sep = "|"))
nsp  <- function(s, g) { v <- NSP[paste(s, g, sep = "|")]; ifelse(is.na(v), 0, v) }
valid <- function(s, g) nsp(s, g) >= MIN_SPOTS
samples <- names(grp)[valid(names(grp), YAP)]          # patients with >= 10 YAP-Rev high spots (12 Relapse, 34 No Relapse)

# two-sided Mann-Whitney, as scipy.stats.mannwhitneyu(method='auto'): exact when one sample has <= 8 values and there
# are no ties, otherwise normal approximation with continuity correction
mw_p <- function(a, b) {
  exact <- min(length(a), length(b)) <= 8 && !anyDuplicated(c(a, b))
  suppressWarnings(wilcox.test(a, b, exact = exact, correct = TRUE)$p.value)
}
fmt_p <- function(p) ifelse(p < 0.001, "P < 0.001", paste("P =", formatC(p, digits = 2, format = "g")))

# per sample / sender / receiver: strength (sum of prob), number of interactions, immune and collagen VI strength
agg <- d[, .(strength = sum(prob), n_int = .N, immune = sum(prob[pathway_name %in% IMM]),
             col6 = sum(prob[ligand %in% COL6]), COL6A1 = sum(prob[ligand == "COL6A1"]),
             COL6A2 = sum(prob[ligand == "COL6A2"]), COL6A3 = sum(prob[ligand == "COL6A3"])), by = .(sample, source, target)]
agg[, k := paste(sample, source, target, sep = "|")]
getv <- function(col, s, a, b) { v <- agg[[col]][match(paste(s, a, b, sep = "|"), agg$k)]; ifelse(is.na(v), 0, v) }

# communication between YAP-Rev high and a partner (one or several CellChat groups; several = spot-weighted mean over the
# groups with >= 10 spots in the patient); direction "to" = YAP -> partner, "from" = partner -> YAP, "both" = sum
partner_value <- function(col, s, groups, direction) {
  g <- groups[valid(s, groups)]
  if (!length(g)) return(NA_real_)
  w <- nsp(s, g); w <- w / sum(w)
  v <- switch(direction, to = getv(col, s, YAP, g), from = getv(col, s, g, YAP),
              both = getv(col, s, YAP, g) + getv(col, s, g, YAP))
  sum(w * v)
}
DIR_LAB <- c(to = "YAP-Rev high → %s", from = "%s → YAP-Rev high", both = "YAP-Rev high ↔ %s")
DIR_TAG <- c(to = "YAPRevHigh_to_%s", from = "%s_to_YAPRevHigh", both = "YAPRevHigh_%s_bidirectional")

# ---- grouped boxplot, No Relapse vs Relapse per partner (style of the Collagen VI boxplot) ----
# V: data.table(partner, sample, value); returns the test table
grouped_boxplot <- function(V, parts, ylab, file, width, height, rotate = FALSE) {
  V <- V[!is.na(value)]
  V[, group := factor(grp[sample], levels = names(GRP_COL))]
  V[, i := match(partner, parts) - 1L][, xp := i + (as.integer(group) - 1.5) * 0.36]
  set.seed(3); V[, xj := xp + runif(.N, -0.07, 0.07)]
  st <- V[, .(n_R = sum(group == "Relapse"), n_NR = sum(group == "No Relapse"),
              median_R = median(value[group == "Relapse"]), median_NR = median(value[group == "No Relapse"]),
              mean_R = mean(value[group == "Relapse"]), mean_NR = mean(value[group == "No Relapse"]),
              P = if (sum(group == "Relapse") >= 3 && sum(group == "No Relapse") >= 3)
                    mw_p(value[group == "Relapse"], value[group == "No Relapse"]) else NA_real_,
              top = max(value)), by = .(partner, i)]
  setorder(st, i)
  h <- max(V$value) * 0.03                                   # additive offset for the P brackets
  br <- st[!is.na(P)]
  ticks <- sprintf("%s\n%d NR - %d R", st$partner, st$n_NR, st$n_R)
  p <- ggplot(V) +
    geom_boxplot(aes(xp, value, group = xp, fill = group), width = 0.3, outlier.shape = NA, alpha = 0.55,
                 linewidth = 0.3, colour = "grey25") +
    geom_point(aes(xj, value, fill = group), shape = 21, size = 1.2, stroke = 0.2, colour = "black") +
    geom_segment(data = br, aes(x = i - 0.18, xend = i + 0.18, y = top + 2 * h, yend = top + 2 * h), linewidth = 0.25) +
    geom_segment(data = br, aes(x = i - 0.18, xend = i - 0.18, y = top + h, yend = top + 2 * h), linewidth = 0.25) +
    geom_segment(data = br, aes(x = i + 0.18, xend = i + 0.18, y = top + h, yend = top + 2 * h), linewidth = 0.25) +
    geom_text(data = br, aes(i, top + 2.4 * h, label = fmt_p(P)), vjust = 0, size = 7.5 / PT) +
    scale_fill_manual(values = GRP_COL, name = NULL) +
    scale_x_continuous(breaks = seq_along(parts) - 1, labels = ticks, limits = c(-0.6, length(parts) - 0.4), expand = c(0, 0)) +
    labs(x = NULL, y = ylab) + theme_classic(base_size = 9) +
    theme(axis.text.x = if (rotate) element_text(angle = 45, hjust = 1, vjust = 1, size = 7.5, colour = "black")
                        else element_text(size = 8.5, colour = "black"),
          axis.text.y = element_text(colour = "black"), legend.position = "right", legend.justification = "top")
  ggsave(file.path(out_dir, file), p, width = width, height = height, device = cairo_pdf)
  st[, c("i", "top") := NULL][]
}

# ---------------------------------------------------------------- Section A: collagen VI, Stroma-A1..B2
cat("\n===== A. Collagen VI YAP-Rev high <-> stromal subclusters =====\n")
# Per patient, sum of the probabilities of all significant
# interactions with a COL6A1/A2/A3 ligand; patients with >= 10 spots in both groups
colA <- rbindlist(lapply("from", function(dr) {
  V <- rbindlist(lapply(names(SUBS), function(pl) data.table(partner = pl, sample = samples,
                 value = sapply(samples, function(s) partner_value("col6", s, SUBS[[pl]], dr)))))
  fn <- sprintf("Boxplot_CollagenVI_%s_RelapseVsNoRelapse.pdf", sprintf(DIR_TAG[dr], "substroma"))
  st <- grouped_boxplot(V, names(SUBS), paste0("Collagen VI communication\n", sprintf(DIR_LAB[dr], "stromal subcluster")),
                        fn, 6.8, 3.6)
  st[, direction := sprintf(DIR_LAB[dr], "partner")]
}))
print(colA[, .(direction, partner, n_R, n_NR, P = signif(P, 3))])
fwrite(colA, file.path(out_dir, "CollagenVI_substroma_to_YAPRevHigh_tests.csv"))

# ---------------------------------------------------------------- Section B: every metacluster / stromal subcluster
cat("\n===== B. YAP-Rev high <-> metaclusters / stromal subclusters =====\n")
METRICS <- list(CommStrength = c("strength", "Communication strength"))
SETS <- list(substroma = list(parts = SUBS, word = "stromal subcluster", width = 6.8, rotate = FALSE))
for (metric in names(METRICS)) {
  col <- METRICS[[metric]][1]
  res <- rbindlist(lapply(names(SETS), function(sn) { S <- SETS[[sn]]
    rbindlist(lapply("both", function(dr) {
      V <- rbindlist(lapply(names(S$parts), function(pl) data.table(partner = pl, sample = samples,
                     value = sapply(samples, function(s) partner_value(col, s, S$parts[[pl]], dr)))))
      fn <- sprintf("Boxplot_%s_%s_RelapseVsNoRelapse.pdf", metric, sprintf(DIR_TAG[dr], sn))
      st <- grouped_boxplot(V, names(S$parts), paste0(METRICS[[metric]][2], "\n", sprintf(DIR_LAB[dr], S$word)),
                            fn, S$width, 3.9, S$rotate)
      cbind(data.table(partner_set = sn, direction = sprintf(DIR_LAB[dr], "partner")), st)
    }))
  }))
  fwrite(res, file.path(out_dir, sprintf("%s_YAPRevHigh_partners_tests.csv", metric)))
  cat(metric, ": P < 0.05 ->", res[P < 0.05, paste(partner, sub("partner", "X", direction), signif(P, 2))], sep = "\n  ")
}

# ---------------------------------------------------------------- Section C: two-column diagrams, both directions pooled
cat("\n===== C. Two-column LR diagrams, Stroma-A2 / Stroma-B2 <-> YAP-Rev high (pooled) =====\n")
ann <- unique(d[, .(interaction_name, pathway_name, ligand, receptor)], by = "interaction_name")
# per patient and LR pair: prob(YAP -> stroma) + prob(stroma -> YAP); Mann-Whitney Relapse vs No Relapse (tested if
# detected in >= 3 patients); BH-FDR within each stromal subcluster
lr_pooled <- function(st) {
  ss <- samples[valid(samples, st)]
  x <- d[((source == YAP & target == st) | (source == st & target == YAP)) & sample %in% ss]
  M <- dcast(x[, .(v = sum(prob)), by = .(sample, interaction_name)], sample ~ interaction_name, value.var = "v", fill = 0)
  M <- merge(data.table(sample = ss), M, by = "sample", all.x = TRUE); M[is.na(M)] <- 0
  g <- grp[M$sample]
  res <- rbindlist(lapply(setdiff(names(M), "sample"), function(f) {
    a <- M[[f]][g == "Relapse"]; b <- M[[f]][g == "No Relapse"]
    data.table(stroma = LAB[st], feature = f, n_R = length(a), n_NR = length(b), mean_R = mean(a), mean_NR = mean(b),
               det_R = mean(a > 0), det_NR = mean(b > 0), log2FC = log2((mean(a) + 1e-6) / (mean(b) + 1e-6)),
               P = if (sum(a > 0) + sum(b > 0) >= 3) mw_p(a, b) else NA_real_)
  }))
  res <- merge(res, ann, by.x = "feature", by.y = "interaction_name", sort = FALSE)
  res[!is.na(P), FDR := p.adjust(P, "BH")][]
}
# pairs higher in Relapse at nominal P < 0.05, top 25 by P (ties included), ordered by pathway / ligand / receptor
select_pairs <- function(T, maxn = 25) {
  x <- T[log2FC > 0 & P < 0.05][order(P)]
  if (nrow(x) > maxn) x <- x[P <= x$P[maxn] + 1e-12]
  x[order(pathway_name, ligand, receptor, method = "radix")]
}
# one diagram; col = value mapped to colour / width / opacity (scaled to [vmin, vmax]); bg = No Relapse / Relapse frame
draw_twocol <- function(x, st, H, col = "mean_R", vmin = NULL, vmax = NULL, head = "Higher in Relapse (nominal P < 0.05)",
                        bg = NULL, panel_in = 4.5) {
  if (is.null(vmin)) { vmin <- min(x[[col]]); vmax <- max(x[[col]]) }
  ligs <- unique(x$ligand); recs <- unique(x$receptor)
  yl <- setNames(H - 1 - (seq_along(ligs) - 1) * (H - 1) / max(length(ligs) - 1, 1), ligs)
  yr <- setNames(H - 1 - (seq_along(recs) - 1) * (H - 1) / max(length(recs) - 1, 1), recs)
  x <- copy(x)[, q := (get(col) - vmin) / (vmax - vmin + 1e-12)][order(q)]
  x[, id := factor(seq_len(.N))]
  bz <- x[, .(bx = c(0.35, 1, 1, 1.65), by = c(yl[ligand], yl[ligand], yr[receptor], yr[receptor]), q = q), by = id]
  nf <- sum(x$FDR < 0.05, na.rm = TRUE); mixed <- nf > 0 && nf < nrow(x)
  # label box extents from label length (7.5 pt italic ~ 0.62 em per character), in data units of a 6.8-unit-wide panel
  char <- 0.62 * 7.5 / 72 / (panel_in / 6.8)
  xmin <- 0.2 - char * max(nchar(ligs)); xmax <- 1.8 + char * max(nchar(gsub("_", "+", recs)))
  top <- if (is.null(bg)) H + 2.7 else H + 3.4
  box <- function(x0, x1) data.table(x = c(x0, x1, x1, x0), y = c(-0.8, -0.8, H - 0.5, H - 0.5))
  p <- ggplot() + coord_cartesian(xlim = c(-2.4, 4.4), ylim = c(-1.2, top), expand = FALSE) + theme_void()
  if (!is.null(bg)) {
    fr <- data.table(x = c(-2.3, 4.3, 4.3, -2.3), y = c(-1.15, -1.15, H + 3.3, H + 3.3))
    p <- p + geom_shape(data = fr, aes(x, y), fill = GRP_COL[bg], alpha = 0.07, radius = unit(3, "mm")) +
      geom_shape(data = fr, aes(x, y), fill = NA, colour = GRP_COL[bg], linewidth = 5 / 2.13, radius = unit(3, "mm"))
  }
  for (b in list(box(xmin - 0.2, 0.55), box(1.45, xmax + 0.2)))
    p <- p + geom_shape(data = b, aes(x, y), fill = "white", colour = NODE[st], linewidth = 4 / 2.13, radius = unit(2.5, "mm")) +
      geom_shape(data = b, aes(x, y), fill = NA, colour = "black", linewidth = 1.2 / 2.13, radius = unit(2.5, "mm"))
  p <- p +
    geom_bezier(data = bz, aes(bx, by, group = id, colour = q, linewidth = (0.5 + 3.5 * q) / 2.13, alpha = 0.3 + 0.7 * q),
                lineend = "butt") +
    scale_colour_gradientn(colours = SPECTRAL_R, limits = c(0, 1), breaks = c(0, 1),
                           labels = if (col == "mean_R" && is.null(bg)) c("min", "max") else formatC(c(vmin, vmax), digits = 2, format = "g"),
                           name = if (is.null(bg)) "Comm. prob. (Relapse, both directions), scaled within each panel"
                                  else sprintf("%s: mean comm. prob. (both directions), common scale NR / R", LAB[st]),
                           guide = guide_colourbar(title.position = "bottom", title.hjust = 0.5, barwidth = unit(4.5, "cm"),
                                                   barheight = unit(0.25, "cm"), frame.colour = "black")) +
    scale_linewidth_identity() + scale_alpha_identity() +
    geom_segment(data = data.table(y = yl), aes(x = 0.28, xend = 0.42, y = y, yend = y), colour = "#555555", linewidth = 4 / 2.13, lineend = "round") +
    geom_segment(data = data.table(y = yr), aes(x = 1.58, xend = 1.72, y = y, yend = y), colour = "#555555", linewidth = 4 / 2.13, lineend = "round") +
    geom_text(data = data.table(y = yl, l = ligs), aes(0.2, y, label = l), hjust = 1, size = 7.5 / PT, fontface = "italic") +
    geom_text(data = data.table(y = yr, l = gsub("_", "+", recs)), aes(1.8, y, label = l), hjust = 0, size = 7.5 / PT, fontface = "italic") +
    annotate("text", x = c((xmin - 0.2 + 0.55) / 2, (1.45 + xmax + 0.2) / 2), y = H - 0.2, label = c("Ligands", "Receptors"),
             vjust = 0, size = 10 / PT, fontface = "bold") +
    annotate("segment", x = 0.55, xend = 1.45, y = H + 0.75, yend = H + 0.75, linewidth = 2 / 2.13,
             arrow = arrow(ends = "both", type = "closed", length = unit(1.6, "mm"))) +
    annotate("text", x = 1, y = H + 1.15, vjust = 0, size = 8 / PT, lineheight = 0.95,
             label = sprintf("%s ↔ YAP-Rev high (both directions pooled)\n%s\n%s", LAB[st], head,
                             if (nf == nrow(x)) "all pairs FDR < 0.05" else if (nf == 0) "no pair at FDR < 0.05"
                             else sprintf("* FDR < 0.05 (%d/%d)", nf, nrow(x)))) +
    theme(legend.position = "bottom", legend.title = element_text(size = 7.5), legend.text = element_text(size = 7))
  if (mixed) {                                             # star on the curve near the ligand end (cubic Bezier, t = 0.25)
    t <- 0.25; s <- x[FDR < 0.05]
    s[, `:=`(sx = (1 - t)^3 * 0.35 + 3 * (1 - t)^2 * t + 3 * (1 - t) * t^2 + t^3 * 1.65,
             sy = ((1 - t)^3 + 3 * (1 - t)^2 * t) * yl[ligand] + (3 * (1 - t) * t^2 + t^3) * yr[receptor])]
    p <- p + geom_text(data = s, aes(sx, sy), label = "*", size = 9 / PT, fontface = "bold")
  }
  p
}
STS <- c("Stroma A 2", "Stroma B 2")
Ts <- lapply(setNames(STS, STS), lr_pooled)
fwrite(rbindlist(Ts), file.path(out_dir, "LR_tests_YAPRevHigh_StromaA2_B2_bidirectional_pooled.csv"))
xs <- lapply(Ts, select_pairs)
H <- max(sapply(xs, function(x) max(uniqueN(x$ligand), uniqueN(x$receptor))))
for (st in STS) cat(sprintf("%s <-> YAP-Rev high: %d pairs shown, %d at FDR < 0.05; %d pairs higher in Relapse at P < 0.05\n",
                            LAB[st], nrow(xs[[st]]), sum(xs[[st]]$FDR < 0.05), Ts[[st]][P < 0.05 & log2FC > 0, .N]))
# C1: Relapse, colour = mean pooled probability in Relapse, scaled within each panel
p <- wrap_plots(lapply(STS, function(st) draw_twocol(xs[[st]], st, H)), nrow = 1) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "TwoColumn_LR_StromaA2_StromaB2_YAPRevHigh_bidirectional_pooled_Relapse.pdf"), p,
       width = 9, height = 0.28 * H + 2.5, device = cairo_pdf)

# ---------------------------------------------------------------- Section D: abundance of YAP-Rev high spots
cat("\n===== D. Number / % of YAP-Rev high spots =====\n")
ny  <- setNames(sapply(names(grp), function(s) nsp(s, YAP)), names(grp))
tot <- n[, .(tot = sum(n_spots)), by = sample][match(names(grp), sample), tot]
D <- data.frame(sample = names(grp), group = unname(grp[names(grp)]), pct = 100 * ny / tot)
a <- D$pct[D$group == "Relapse"]; b <- D$pct[D$group == "No Relapse"]; pD <- mw_p(a, b)
source("R/boxplot_relapse.R")
g <- boxplot_relapse(D, "pct", pD, "YAP-Rev high spots", "% of all spots per patient", sprintf("Patients (n = %d)", nrow(D)))
ggsave(file.path(out_dir, "Boxplot_Pct_YAPRevHigh_spots_RelapseVsNoRelapse.pdf"), g, width = 5, height = 5)
abD <- data.table(metric = "pct_of_all_spots", n_R = length(a), n_NR = length(b), median_R = median(a), median_NR = median(b), P = pD)
print(abD); fwrite(abD, file.path(out_dir, "N_YAPRevHigh_spots_tests.csv"))

