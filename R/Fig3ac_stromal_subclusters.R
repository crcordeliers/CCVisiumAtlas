# Fig. 3a and Fig. 3c - stromal sub-clusters
#   Fig. 3a : UMAPs of the meta-spots (stromal meta-clusters highlighted), of the 774 stromal meta-spots coloured by
#             meta-cluster and by sub-cluster, and the same on one section (#74_B1)
#   Fig. 3c : average ssGSEA scores of CAF / immune signatures per stromal sub-cluster (dot plot; version of the revised
#             manuscript: Buechler NF, Chen matCAF / myCAF / iCAF, Smillie IAF and MCP-counter populations)
# (Fig. 3b, cell2location, is in analyses/04_cell2location_Fig2c_Fig3b; sub-cluster proportions by MMR status / relapse,
#  Supplementary Fig. S6, are in analyses/06_ecotype_composition_S5ab_S6_S7.)

source("R/00_config.R")

library(Seurat)
library(ggplot2)
library(GSVA)
library(tidyverse)
library(readxl)

output_dir <- figures.dir

## ---- Fig. 3a: stromal meta-clusters and sub-clusters (UMAP, section #74_B1) ----
stroma_only_palette <- setNames(rep("gray", length(meta_colors)), names(meta_colors))
stroma_only_palette[stromal_metaclusters] <- meta_colors[stromal_metaclusters]

seurat_obj <- readRDS(file.path(harmony.dir, "seurat_object_res_10_with5subsets.rds"))
ggsave(file.path(output_dir, "Fig3a_UMAP_all_metaspots_stroma.pdf"),
       DimPlot(seurat_obj, group.by = "merged_meta_clusters") + scale_color_manual(values = stroma_only_palette) + NoLegend(),
       width = 7, height = 7)

Cluster2_4 <- readRDS(file.path(harmony.dir, "Cluster2_4.rds"))
ggsave(file.path(output_dir, "Fig3a_UMAP_stromal_metaclusters.pdf"),        # meta-cluster (RNA_snn_res.0.5 holds the sub-clusters here)
       DimPlot(Cluster2_4, group.by = "merged_meta_clusters") + scale_color_manual(values = stroma_only_palette) + NoLegend(),
       width = 7, height = 7)
ggsave(file.path(output_dir, "Fig3a_UMAP_stromal_subclusters.pdf"),
       DimPlot(Cluster2_4, group.by = "stromal_subcluster") +
         scale_color_manual(values = stromal_colors, labels = stromal_names, breaks = stromal_order),
       width = 8, height = 7)

load(file.path(processed.dir, "processed_seurat_Ech_074_B1_withmeta.RData"))   # -> sobj
ggsave(file.path(output_dir, "Fig3a_Ech_074_B1_spatial_stromal_metaclusters.pdf"),
       SpatialDimPlot(sobj, group.by = "merged_meta_clusters", pt.size.factor = 2.5, alpha = 1) +
         scale_fill_manual(values = stroma_only_palette) + NoLegend(),
       width = 7, height = 7)
ggsave(file.path(output_dir, "Fig3a_Ech_074_B1_spatial_stromal_subclusters.pdf"),
       SpatialDimPlot(sobj, group.by = "stromal_subcluster", pt.size.factor = 2.5) +
         scale_fill_manual(values = stromal_colors, labels = stromal_names, breaks = stromal_order, na.value = "gray90"),
       width = 8, height = 7)


## ---- Fig. 3c: CAF / immune signatures per stromal sub-cluster ----
load(file.path(harmony.dir, "seurat_obj_with_pelka.RData"))   # -> seurat_obj (Step 6), carries the MCP-counter scores
seurat_obj_strom <- subset(seurat_obj, subset = merged_meta_clusters %in% stromal_metaclusters)

C3_vs_C8_Turley <- read_excel(file.path(data.dir, "C3_vs_C8_Turley_Nature2021.xlsx"))
signatures_fibro <- list(
  Chen_iCAF   = c("CCL2", "CXCL12", "CXCL14"),                       # Chen et al., Clin Transl Med 2023
  Chen_myCAF  = c("RGS5", "MYH11", "ACTA2"),
  Chen_matCAF = c("COL10A1", "CTHRC1", "POSTN"),
  Buechler_NF = C3_vs_C8_Turley %>% filter(avg_logFC > 0 & pct.1 > 0.1) %>% pull(Gene),   # Buechler et al., Nature 2021
  Smillie_IAF = c("PLAU", "CHI3L1", "MMP3", "IL1R1", "IL13RA2", "TNFSF11", "MMP10", "OSMR", "IL11", "STRA6",
                  "FAP", "WNT2", "TWIST1", "IL24")                     # Smillie et al., Cell 2019
)
sig_fibro <- as.data.frame(t(gsva(ssgseaParam(as.matrix(seurat_obj_strom@assays$RNA$scale.data), signatures_fibro, normalize = TRUE))))
for (sig_name in colnames(sig_fibro)) seurat_obj_strom[[sig_name]] <- sig_fibro[[sig_name]]

select <- c("Buechler_NF", "Chen_matCAF", "Monocytic_lineage", "T_cells", "NK_cells", "Cytotoxic_lymphocytes", "CD8_T_cells",
            "B_lineage", "Myeloid_dendritic_cells", "Chen_myCAF", "Chen_iCAF", "Smillie_IAF", "Neutrophils", "Endothelial_cells")
rename_map <- c(Monocytic_lineage = "Becht_Monocytic", T_cells = "Becht_T_cells", NK_cells = "Becht_NK_cells",
                Cytotoxic_lymphocytes = "Becht_CTLs", CD8_T_cells = "Becht_CD8_T_cells", B_lineage = "Becht_B_cells",
                Myeloid_dendritic_cells = "Becht_Myeloid_DCs", Neutrophils = "Becht_Neutrophils", Endothelial_cells = "Becht_ECs")

custom_order <- rev(stromal_order)                                 # Stroma-B3 ... Stroma-A1 (bottom to top)
seurat_obj_strom$ident_ordered <- factor(seurat_obj_strom$stromal_subcluster, levels = custom_order,
                                         labels = unname(stromal_names[custom_order]))
Idents(seurat_obj_strom) <- seurat_obj_strom$ident_ordered

select_colours <- c("#0000FF", "#3939FF", "#7171FF", "#AAAAFF", "#E3E3FF", "#F2E3E3", "#D8AAAA", "#BF7171", "#A53939", "#8B0000")
dp_data <- DotPlot(seurat_obj_strom, features = select)$data    # avg.exp.scaled = average score scaled across sub-clusters

DP_fixed <- ggplot(dp_data, aes(y = id, x = features.plot)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, colour = "black", stroke = 0.4, alpha = 0.8) +
  scale_fill_gradientn(colours = select_colours) +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(rename_map), rename_map[x], x)) +
  scale_size(range = c(3, 10)) +
  theme_minimal(base_size = 9) +
  theme(axis.title = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = rel(1.5), angle = 45, hjust = 1), axis.text.y = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.2))) +
  labs(size = "Percent Expressed", fill = "Average Expression")
ggsave(file.path(output_dir, "Fig3c_dotplot_stromal_subclusters.pdf"), DP_fixed, width = 9, height = 5)
