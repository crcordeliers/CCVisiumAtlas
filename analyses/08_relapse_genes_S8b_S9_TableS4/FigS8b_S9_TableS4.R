#!/usr/bin/env Rscript
# =====================================================================================================================
# Expression of the recurrence-associated genes (Fig. 4a) by ecotype and relapse status
#   S8b      16 up-regulated genes in the 4 tumour and 2 stromal meta-clusters, No relapse (light) vs Relapse
#            (dark); Wilcoxon P per meta-cluster
#   Table S4 pairwise comparisons between the 6 meta-clusters, separately in patients with and without relapse, for the
#            genes and for the REC / iREC / YAP / RSC signatures (Wilcoxon, BH within gene x relapse group)
#   S9a-d    pMMR patients: scaled expression of the up-regulated genes in Stroma-A1, A2, B1, B2 meta-spots by relapse
#   S9e      pMMR patients: each gene by stromal sub-cluster and relapse status (Wilcoxon P per sub-cluster)
# Inputs : seurat_obj_with_pelka.RData (meta-spots, R/Step6), Cluster2_4.rds (stromal meta-spots, R/Step4),
#          data/signatures_fig3.RData
# Outputs: Figures/FigS8b/, Figures/FigS9/, Tables/TableS4_*.csv
# Usage  : Rscript analyses/08_relapse_genes_S8b_S9_TableS4/FigS8b_S9_TableS4.R
# recurr_cens = 1 means relapse.
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(GSVA); library(dplyr); library(tidyr); library(purrr); library(ggplot2) })

out8 <- file.path(figures.dir, "FigS8b"); out9 <- file.path(figures.dir, "FigS9")
for (d in c(out8, out9)) dir.create(d, showWarnings = FALSE)
UP <- c("TMEM159", "APLP2", "CTXN1", "CAPN2", "DDX3Y", "TBL1X", "RASAL2", "NT5DC2", "CORO1C", "BOK", "AHNAK",
        "LRP10", "ANXA1", "AZIN2", "PPFIBP1", "NFIC")
GENES <- c(UP, "HLA-A")   # HLA-A: Table S4 and Wilcoxon P tables only
p_lab <- function(p) case_when(p < 1e-4 ~ "italic(P) < '0.0001'", p < 1e-3 ~ "italic(P) < '0.001'",
                               !is.na(p) ~ paste0("italic(P) == '", formatC(p, format = "f", digits = 3), "'"), TRUE ~ "italic(P) == NA")

# boxplots of one variable by cluster x relapse (light = no relapse, dark = relapse), Wilcoxon P per cluster
box_by_cluster <- function(df, var, light, dark, ylab, xlabels = NULL) {
  v <- df[[var]]; q <- quantile(v, c(.25, .75), na.rm = TRUE); iqr <- diff(q)
  y_min <- q[1] - 1.5 * iqr; y_max <- q[2] + 1.5 * iqr; rng <- y_max - y_min
  pv <- df |> group_by(cluster) |>
    summarise(pval = tryCatch(wilcox.test(.data[[var]] ~ recurr_cens)$p.value, error = function(e) NA), .groups = "drop") |>
    mutate(p_label = p_lab(pval), y_pos = y_max + 0.07 * rng)
  fills <- setNames(c(light, dark), c(paste0(names(light), "_0"), paste0(names(dark), "_1")))
  p <- ggplot(df, aes(cluster, .data[[var]], fill = interaction(cluster, recurr_cens, sep = "_"))) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = fills) +
    { if (!is.null(xlabels)) scale_x_discrete(labels = xlabels) } +
    geom_text(data = pv, aes(cluster, y_pos, label = p_label), parse = TRUE, size = 7, fontface = "bold", inherit.aes = FALSE) +
    scale_y_continuous(limits = c(y_min - 0.05 * rng, y_max + 0.15 * rng), expand = expansion(mult = c(0.02, 0.05))) +
    labs(title = var, x = "", y = ylab) + theme_bw() +
    theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          axis.text = element_text(face = "bold", size = rel(1.6), color = "black"), axis.title.y = element_text(face = "bold", size = rel(1.8)),
          plot.title = element_text(face = "bold", size = rel(1.6), hjust = 0.5), legend.position = "none")
  list(plot = p, p = pv)
}

## ---- S8b: genes by tumour / stromal meta-cluster and relapse ----
load(file.path(harmony.dir, "seurat_obj_with_pelka.RData"))           # -> seurat_obj
ts <- c("0", "1", "7", "8", "2", "4")
light_mc <- c("0" = "#FF6B6B", "1" = "#69A4DF", "7" = "#FFA7D2", "8" = "#B7B7B7", "2" = "#6ECD6C", "4" = "#FFA984")
dark_mc  <- c("0" = "#9F1718", "1" = "#2D577C", "7" = "#D32092", "8" = "#696969", "2" = "#397737", "4" = "#B05601")
df_mc <- FetchData(seurat_obj, vars = c("recurr_cens", "merged_meta_clusters", GENES)) |>
  mutate(cluster = factor(merged_meta_clusters, levels = ts), recurr_cens = factor(recurr_cens)) |>
  filter(!is.na(cluster))
pv8 <- list()
for (g in GENES) {
  r <- box_by_cluster(df_mc, g, light_mc, dark_mc, "Scaled Expression")
  if (g %in% UP) ggsave(file.path(out8, paste0("FigS8b_", g, "_by_metacluster_relapse.pdf")), r$plot, width = 10, height = 5)
  pv8[[g]] <- mutate(r$p, gene = g)
}
write.csv(bind_rows(pv8), file.path(tables.dir, "FigS8b_Wilcoxon_P.csv"), row.names = FALSE)

## ---- Table S4: pairwise comparisons between meta-clusters, within each relapse group ----
load(file.path(data.dir, "signatures_fig3.RData"))                    # -> signatures (REC, iREC, YAP, RSC)
sig <- as.data.frame(t(gsva(ssgseaParam(as.matrix(seurat_obj@assays$RNA$scale.data), signatures, normalize = TRUE))))
sig_labels <- c(Ogden_REC = "Ogden_REC", Ogden_iREC = "Ogden_iREC", Gregorieff_YAP_DOWN = "Gregorieff_YAP_DOWN",
                Gregorieff_YAP_UP = "Gregorieff_YAP_UP", YAP = "Serra_YAP", RSC = "Ayyaz_RSC")
sig <- sig[, names(sig_labels)]; colnames(sig) <- sig_labels
df_s4 <- cbind(df_mc, sig[rownames(df_mc), ])
relapse_labels <- c("1" = "relapse", "0" = "no relapse")
pairs <- combn(ts, 2, simplify = FALSE)

pairwise <- function(vars) {
  res <- map_dfr(vars, function(v) map_dfr(names(relapse_labels), function(rc) map_dfr(pairs, function(pr) {
    d <- df_s4[df_s4$recurr_cens == rc, ]
    x <- d[[v]][d$cluster == pr[1]]; y <- d[[v]][d$cluster == pr[2]]
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    tibble(variable = v, relapse_status = relapse_labels[[rc]], comparison = paste0("cluster_", pr[1], "_vs_", pr[2]),
           n_cluster_1 = length(x), n_cluster_2 = length(y), median_cluster_1 = median(x), median_cluster_2 = median(y),
           p_value = if (length(x) >= 2 && length(y) >= 2) wilcox.test(x, y, exact = FALSE)$p.value else NA_real_)
  })))
  res |> group_by(variable, relapse_status) |> mutate(p_adjusted = p.adjust(p_value, "BH")) |> ungroup()
}
for (set in list(list("genes", GENES), list("signatures", unname(sig_labels)))) {
  long <- pairwise(set[[2]])
  write.csv(long, file.path(tables.dir, paste0("TableS4_", set[[1]], "_pairwise_long.csv")), row.names = FALSE)
  for (col in c("p_value", "p_adjusted"))
    write.csv(pivot_wider(long[, c("variable", "relapse_status", "comparison", col)], names_from = comparison, values_from = all_of(col)),
              file.path(tables.dir, paste0("TableS4_", set[[1]], "_", col, ".csv")), row.names = FALSE)
}

## ---- S9: pMMR patients, stromal sub-clusters ----
cl24 <- readRDS(file.path(harmony.dir, "Cluster2_4.rds"))
df_sc <- FetchData(cl24, vars = c("recurr_cens", "mmr_status", "stromal_subcluster", GENES)) |>
  filter(mmr_status == "pMMR") |>
  mutate(cluster = factor(stromal_subcluster, levels = stromal_order), recurr_cens = factor(recurr_cens))

# S9a-d: all up-regulated genes in one sub-cluster
for (k in c("0", "2", "1", "3")) {
  dk <- filter(df_sc, cluster == k)
  n <- table(dk$recurr_cens)
  dl <- pivot_longer(dk[, c("recurr_cens", UP)], all_of(UP), names_to = "gene", values_to = "value")
  p <- ggplot(dl, aes(gene, value, fill = recurr_cens)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = relapse_colors, labels = c("0" = paste0("No Relapse n = ", n[["0"]]), "1" = paste0("Relapse n = ", n[["1"]]))) +
    scale_y_continuous(limits = c(-0.5, 1.5)) +
    labs(title = "pMMR", x = "", y = paste0("Scaled Expression - ", stromal_names[[k]]), fill = "") + theme_bw() +
    theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = rel(1.6), color = "black"),
          axis.text.y = element_text(face = "bold", size = rel(1.6), color = "black"), axis.title.y = element_text(face = "bold", size = rel(1.8)),
          plot.title = element_text(face = "bold", size = rel(1.6), hjust = 0.5), legend.text = element_text(face = "bold", size = rel(1.4)),
          legend.position = "top")
  ggsave(file.path(out9, paste0("FigS9_", stromal_names[[k]], "_pMMR.pdf")), p, width = 10, height = 5)
}
# S9e: each gene across sub-clusters
light_sc <- c("0" = "#B9D9F9", "2" = "#BAF4D7", "1" = "#FFD6BA", "3" = "#FFC0C0", "4" = "#FFEEAB")
dark_sc  <- c("0" = "#5A89B0", "2" = "#4FA57F", "1" = "#BD7F3E", "3" = "#CA5F5D", "4" = "#B19C07")
pv9 <- list()
for (g in GENES) {
  r <- box_by_cluster(df_sc, g, light_sc, dark_sc, "Scaled Expression - pMMR",
                      xlabels = c("0" = "A1", "2" = "A2", "1" = "B1", "3" = "B2", "4" = "B3"))
  if (g %in% UP) ggsave(file.path(out9, paste0("FigS9e_", g, "_by_subcluster_relapse_pMMR.pdf")), r$plot, width = 10, height = 5)
  pv9[[g]] <- mutate(r$p, gene = g)
}
write.csv(bind_rows(pv9), file.path(tables.dir, "FigS9e_Wilcoxon_P.csv"), row.names = FALSE)
