# Fig. 2a and Supplementary Fig. S4e - spatial ecotypes (meta-clusters)
#   Fig. 2a : heatmap of signature scores (ssGSEA) across the 3,231 meta-spots, split by meta-cluster
#   S4e     : neighbourhood of stroma-A / stroma-B spots (alluvial plot, first crown = 6 nearest spots) and share of
#             stroma-A among tumour-stroma contacts (77.7%), tested per patient against the share expected from the
#             abundance of stroma-A (paired Wilcoxon test; observed > expected in 42/48 samples)
# (The Pelka et al. epithelial / stromal violin plots of S4a,b are in analyses/05_ecotype_hallmarks_S4.)

source("R/00_config.R")

library(Seurat)
library(GSVA)
library(ComplexHeatmap)
library(circlize)
library(Cairo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)
library(RANN)
library(ggalluvial)

## ======================= Fig. 2a : signature heatmap =======================
seurat_obj <- readRDS(file.path(harmony.dir, "seurat_object_res_10.rds"))
load(file.path(data.dir, "signatures.RData"))

sig <- as.data.frame(t(gsva(ssgseaParam(as.matrix(seurat_obj[["RNA"]]$scale.data), signatures, normalize = TRUE))))

signatures_heatmap <- c("CMS1", "CMS2", "CMS3", "CMS4",
                        "T_cells", "NK_cells", "Cytotoxic_lymphocytes", "B_lineage", "Myeloid_dendritic_cells",
                        "Monocytic_lineage", "Neutrophils",
                        "Epithelial_sig", "Endothelial_cells", "Mesenchymal_sig", "Fibroblasts", "Intestinal_stem_cells",
                        "EMT", "TGFB", "WNT", "MYC", "SRC", "JAK_STAT", "Caspases", "Proteasome", "CellCycle",
                        "integrinB3", "VEGF", "Glycolysis", "Fatty_acids")

mc_fac <- factor(seurat_obj$merged_meta_clusters, levels = meta_order)
annotation_heatmap <- columnAnnotation(
  meta_clusters = mc_fac,
  col = list(meta_clusters = meta_colors),
  annotation_legend_param = list(meta_clusters = list(title = "Meta-clusters", at = meta_order, labels = meta_names[meta_order]))
)

Heatmap_metaClusters <- Heatmap(
  t(as.matrix(sig[colnames(seurat_obj), signatures_heatmap])),
  top_annotation    = annotation_heatmap,
  col               = colorRamp2(c(-0.5, 0, 0.5), c("navy", "white", "red")),
  cluster_rows      = FALSE,
  cluster_columns   = TRUE,
  show_row_names    = TRUE,
  show_column_names = FALSE,
  column_split      = mc_fac,
  heatmap_legend_param = list(title = "")
)
CairoPDF(file.path(figures.dir, "Fig2a_Heatmap_metaClusters.pdf"), width = 30, height = 10)
print(Heatmap_metaClusters)
dev.off()

## ======================= S4e : spatial neighbourhood =======================
cluster_var <- "merged_meta_clusters"
lev_mc      <- as.character(0:10)
nn_cols     <- paste0("nn_", lev_mc)
neighbors_rds <- file.path(tables.dir, "neighbors_meta_counts_ALL_rings.rds")

# Visium spot coordinates
get_xy <- function(seu) {
  xy <- try(Seurat::GetTissueCoordinates(seu, filter = "tissue"), silent = TRUE)
  if (!inherits(xy, "try-error")) return(tibble(barcode = rownames(xy), x = xy[[1]], y = xy[[2]]))
  img <- seu@images[[1]]@coordinates
  tibble(barcode = rownames(img), x = img$pxl_col_in_fullres, y = img$pxl_row_in_fullres)
}

# For each spot, meta-cluster counts among its neighbours: first crown (6 nearest), second crown (next 12), both (18)
neighbor_meta_counts_by_rings <- function(seu, cluster_var, levels_all) {
  coords     <- get_xy(seu)
  cells_keep <- intersect(rownames(seu@meta.data), coords$barcode)
  coords     <- coords[match(cells_keep, coords$barcode), , drop = FALSE]
  labels_chr <- as.character(seu@meta.data[cells_keep, cluster_var, drop = TRUE])
  labels     <- factor(labels_chr, levels = levels_all)

  # 18 nearest spots (exact Euclidean kNN; first column of nn2 is the spot itself).
  nn_mat <- RANN::nn2(as.matrix(coords[, c("x", "y")]), k = 19)$nn.idx[, -1, drop = FALSE]

  count_ring <- function(idx_cols, ring_label) {
    counts_mat <- t(apply(nn_mat[, idx_cols, drop = FALSE], 1, function(nn_idx) as.integer(table(labels[nn_idx]))))
    colnames(counts_mat) <- paste0("nn_", levels_all)
    tibble(ring = ring_label, barcode = coords$barcode, meta_cluster = labels_chr) |>
      bind_cols(as_tibble(counts_mat))
  }
  bind_rows(count_ring(1:6, "ring1"), count_ring(7:18, "ring2_only"), count_ring(1:18, "ring1_2"))
}

if (!file.exists(neighbors_rds)) {
  nn_all <- map_dfr(Sample.names, function(sname) {
    load(file.path(processed.dir, paste0("processed_seurat_", sname, "_withmeta.RData")))   # -> sobj
    neighbor_meta_counts_by_rings(sobj, cluster_var, lev_mc) |> mutate(sample = sname, .before = 1)
  })
  saveRDS(nn_all, neighbors_rds)
  write_csv(nn_all, sub("\\.rds$", ".csv", neighbors_rds))
}
nn_all <- readRDS(neighbors_rds)

ring_sel <- "ring1"   # first spatial crown (6 closest spots)

## ---- Neighbour flows of stroma-A (MC2) and stroma-B (MC4) spots ----
flows <- nn_all |>
  filter(ring == ring_sel, meta_cluster %in% stromal_metaclusters) |>
  select(meta_cluster, all_of(nn_cols)) |>
  pivot_longer(all_of(nn_cols), names_to = "neighbor", values_to = "n") |>
  group_by(meta_cluster, neighbor) |>
  summarise(value = sum(n, na.rm = TRUE), .groups = "drop") |>
  filter(value > 0) |>
  mutate(center      = meta_names[meta_cluster],
         target_code = sub("^nn_", "", neighbor),
         neighbor    = meta_names[target_code])
write_csv(flows, file.path(tables.dir, "FigS4e_stroma_neighbors_ring1.csv"))

df_alluv <- flows |> transmute(center, neighbor = factor(neighbor, levels = meta_names), value)

sankey <- ggplot(df_alluv, aes(axis1 = center, axis2 = neighbor, y = value)) +
  geom_alluvium(aes(fill = neighbor), width = 0.3, alpha = 0.55, color = NA) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.3, color = "black") +
  scale_fill_manual(values = meta_colors_named, name = "Meta-clusters") +
  theme_void() +
  theme(legend.position = "right",
        legend.title    = element_text(size = 12, face = "bold"),
        legend.text     = element_text(size = 10))
pdf(file.path(figures.dir, "FigS4e_alluvial_stroma_neighbors.pdf"), width = 10, height = 6)
print(sankey)
dev.off()

## ---- Tumour contacts: share of stroma-A vs stroma-B (77.7% vs 22.3%, 3.48-fold) ----
# Contacts counted from tumour spots (meta-clusters 0, 1, 7, 8) towards stromal spots in the first crown.
tumor_codes <- c("0", "1", "7", "8")
tumor_contacts <- nn_all |>
  filter(ring == ring_sel, meta_cluster %in% tumor_codes) |>
  summarise(stroma_A = sum(nn_2), stroma_B = sum(nn_4))
share_A <- tumor_contacts$stroma_A / (tumor_contacts$stroma_A + tumor_contacts$stroma_B)

# Null hypothesis: stroma-A contacts proportional to the share of stroma-A among stromal spots
n_spots <- nn_all |> filter(ring == ring_sel, meta_cluster %in% stromal_metaclusters) |> count(meta_cluster)
p0 <- n_spots$n[n_spots$meta_cluster == "2"] / sum(n_spots$n)
bt <- binom.test(tumor_contacts$stroma_A, tumor_contacts$stroma_A + tumor_contacts$stroma_B, p = p0)

stats_2D <- data.frame(
  contacts_stromaA = tumor_contacts$stroma_A, contacts_stromaB = tumor_contacts$stroma_B,
  pct_stromaA = round(100 * share_A, 1), pct_stromaB = round(100 * (1 - share_A), 1),
  fold = round(tumor_contacts$stroma_A / tumor_contacts$stroma_B, 2),
  expected_pct_stromaA = round(100 * p0, 1), binom_p = bt$p.value)
print(stats_2D)
write_csv(stats_2D, file.path(tables.dir, "FigS4e_tumor_contacts_stromaA_vs_B.csv"))

## ---- Per patient: observed vs expected share of stroma-A among tumour-stroma contacts (legend of S4e) ----
pp <- nn_all |>
  filter(ring == ring_sel, meta_cluster %in% tumor_codes) |>
  group_by(sample) |> summarise(A = sum(nn_2), B = sum(nn_4), .groups = "drop")
sp <- nn_all |>
  filter(ring == ring_sel, meta_cluster %in% stromal_metaclusters) |>
  count(sample, meta_cluster) |>
  pivot_wider(names_from = meta_cluster, values_from = n, values_fill = 0)
pp <- left_join(pp, sp, by = "sample") |>
  mutate(pct_A = 100 * A / (A + B), expected_pct_A = 100 * `2` / (`2` + `4`))
wt <- wilcox.test(pp$pct_A, pp$expected_pct_A, paired = TRUE)
message("observed > expected in ", sum(pp$pct_A > pp$expected_pct_A), "/", nrow(pp), " samples; paired Wilcoxon P = ", signif(wt$p.value, 3))
write_csv(pp, file.path(tables.dir, "FigS4e_tumor_contacts_per_patient.csv"))
