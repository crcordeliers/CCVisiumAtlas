#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S10c,d - genes correlated with ANXA1 in the stromal meta-spots (stroma-A + stroma-B), by relapse
#   S10c  for each relapse group, the 50 genes most correlated with ANXA1 (Pearson, all genes of the meta-spot matrix);
#         heatmap of the correlation coefficients in both groups for the genes of the relapse top 50 (genes specific to
#         the relapse group + genes common to both), ranked by decreasing correlation in the relapse group
#   S10d  scatter plots ANXA1 vs COL6A1 / COL6A2 / COL6A3, relapse (red) vs no relapse (blue), Pearson R and P per group
# Input : seurat_obj_with_pelka.RData (meta-spots, R/Step6); "data" layer = averaged scaled expression of the meta-spots
# Output: Figures/FigS10c_heatmap_ANXA1_correlated_genes.pdf, Figures/FigS10d_scatter_ANXA1_<gene>.pdf, Tables/FigS10c_*.csv
# Usage : Rscript analyses/12_ANXA1_correlations_S10cd/FigS10cd_ANXA1_correlations.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(ComplexHeatmap); library(circlize) })

load(file.path(harmony.dir, "seurat_obj_with_pelka.RData"))          # -> seurat_obj
stroma <- subset(seurat_obj, subset = merged_meta_clusters %in% stromal_metaclusters)
gene_A <- "ANXA1"; top_n <- 50

cor_with_A <- function(group) {
  expr <- GetAssayData(subset(stroma, subset = recurr_cens == group), layer = "data")
  a <- expr[gene_A, ]
  stopifnot(mean(a > 0) >= 0.10, sd(a) > 0)
  r <- apply(expr, 1, function(x) cor(x, a, method = "pearson"))
  sort(r[!names(r) %in% gene_A & !is.na(r)], decreasing = TRUE)
}
cors <- list("0" = cor_with_A("0"), "1" = cor_with_A("1"))
top <- lapply(cors, head, top_n)
common <- intersect(names(top[["0"]]), names(top[["1"]]))
genes_hm <- c(setdiff(names(top[["1"]]), names(top[["0"]])), common)     # relapse top 50 (specific + common)

corr_matrix <- cbind(R0 = cors[["0"]][genes_hm], R1 = cors[["1"]][genes_hm])
corr_matrix <- corr_matrix[order(corr_matrix[, "R1"], decreasing = TRUE), , drop = FALSE]
write.csv(data.frame(gene = rownames(corr_matrix), r_no_relapse = corr_matrix[, 1], r_relapse = corr_matrix[, 2]),
          file.path(tables.dir, "FigS10c_ANXA1_correlated_genes.csv"), row.names = FALSE)

col_ha <- HeatmapAnnotation(df = data.frame(RelapseStatus = factor(c("R0", "R1")), row.names = colnames(corr_matrix)),
                            col = list(RelapseStatus = c(R0 = "#9197FF", R1 = "#FC7171")),
                            show_legend = FALSE, show_annotation_name = FALSE)
HM <- Heatmap(corr_matrix, name = "Pearson r", col = colorRamp2(c(0, 0.4, 0.8), c("blue", "white", "red")),
              cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = TRUE,
              column_names_side = "top", column_labels = c("No relapse", "Relapse"), column_names_rot = 0,
              column_names_gp = gpar(fontface = "bold"), column_names_centered = TRUE, top_annotation = col_ha)
pdf(file.path(figures.dir, "FigS10c_heatmap_ANXA1_correlated_genes.pdf"), width = 10, height = 10)
print(HM)
dev.off()

## ---- S10d: ANXA1 vs collagen VI genes ----
data_exp <- FetchData(stroma, vars = c("recurr_cens", gene_A, "COL6A1", "COL6A2", "COL6A3"))
for (g in c("COL6A1", "COL6A2", "COL6A3")) {
  y_min <- min(data_exp$ANXA1); y_off <- 0.06 * diff(range(data_exp$ANXA1))
  st <- data_exp |> group_by(recurr_cens) |>
    summarise(r = cor(ANXA1, .data[[g]]), p = cor.test(ANXA1, .data[[g]])$p.value, .groups = "drop") |>
    mutate(label = paste0("R = ", round(r, 2), ", p = ", ifelse(p < 1e-4, "<0.0001", signif(p, 2))),
           x = min(data_exp[[g]]), y = ifelse(recurr_cens == 1, y_min + 2 * y_off, y_min + y_off))
  p <- ggplot(data_exp, aes(.data[[g]], ANXA1, color = factor(recurr_cens))) +
    geom_point() + geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    scale_color_manual(values = relapse_colors, labels = c("0" = "No relapse", "1" = "Relapse")) +
    geom_text(data = st, aes(x, y, label = label, color = factor(recurr_cens)), fontface = "bold", hjust = 0, vjust = 0,
              size = 5, show.legend = FALSE) +
    labs(x = g, y = "ANXA1", title = "Stroma") + theme_bw() +
    theme(legend.position = "none", axis.text = element_text(face = "bold", size = rel(1.2)),
          axis.title = element_text(face = "bold", size = rel(1.2)), plot.title = element_text(face = "bold", size = rel(1.5)))
  ggsave(file.path(figures.dir, paste0("FigS10d_scatter_ANXA1_", g, ".pdf")), p, width = 4, height = 4)
}
