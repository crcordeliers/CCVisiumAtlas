# Fig. 4a and Supplementary Fig. S8a, S8c, S8f - recurrence-associated genes in the stromal compartment
#   Fig. 4a : volcano plot, stromal meta-spots (stroma-A + stroma-B, n = 774) of patients with (n = 247 meta-spots) vs
#             without (n = 527) relapse; Wilcoxon (FindMarkers); |log2FC| > 0.6, adjusted P < 0.01, pct > 0.5 in the
#             group where the gene is higher. The meta-spot "data" layer holds averaged scaled values, so avg_log2FC
#             is not a true log fold-change.
#   S8a     : scaled expression of the up-regulated genes in all meta-spots, by relapse status
#   S8c, S8f: REC / iREC (Ogden), YAP (Gregorieff, Serra) and RSC (Ayyaz) signatures in all vs stromal meta-spots
# (S8b, per tumour / stromal ecotype, and S9, per stromal sub-cluster, are in analyses/08_relapse_genes_S8b_S9_TableS4.)

source("R/00_config.R")

library(ggplot2)
library(Seurat)
library(ggrepel)
library(tidyverse)

## ---- Fig. 4a: volcano plot ----
seurat_obj <- readRDS(file.path(harmony.dir, "seurat_object_res_10.rds"))

seurat_obj_strom <- subset(
  seurat_obj, subset = merged_meta_clusters %in% stromal_metaclusters)

Idents(object = seurat_obj_strom) <- "recurr_cens"

# R1 versus R0 
metacluster_strom_R1vsR0_DEA_markers <- FindMarkers(
  seurat_obj_strom,
  ident.1 = "1",
  ident.2 = "0",
  verbose = TRUE)

# Visualize 
metacluster_strom_R1vsR0_DEA_markers <- metacluster_strom_R1vsR0_DEA_markers %>%
  rownames_to_column(var = "gene") %>%
  mutate(significance = case_when(
    p_val_adj < 0.01 & avg_log2FC > 0.6  & pct.1>0.5~ "Upregulated",
    p_val_adj < 0.01 & avg_log2FC < -0.6  & pct.2>0.5~ "Downregulated",
    TRUE ~ "Not significant"
  ))
save(metacluster_strom_R1vsR0_DEA_markers, file = file.path(tables.dir, "metacluster_strom_R1vsR0_DEA_markers_filtre_pct.RData"))
# Table underlying the volcano plot (all genes)
write.csv(metacluster_strom_R1vsR0_DEA_markers[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "significance")],
          file.path(tables.dir, "Fig4a_volcano_stroma_R1vsR0_DEA.csv"), row.names = FALSE)

# Select top genes up and down
top_genes_up <- metacluster_strom_R1vsR0_DEA_markers %>%
  filter(significance == "Upregulated" & pct.1>0.5) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50)

top_genes_down <- metacluster_strom_R1vsR0_DEA_markers %>%
  filter(significance == "Downregulated" & pct.2>0.5) %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 50)

top_genes <- bind_rows(top_genes_up, top_genes_down)

### filter
metacluster_strom_R1vsR0_DEA_markers_filtre=metacluster_strom_R1vsR0_DEA_markers %>%
  filter(avg_log2FC > -5&avg_log2FC < 5)

### ---- Fig. 4a: volcano plot ---- 
volcano <- ggplot(metacluster_strom_R1vsR0_DEA_markers_filtre, 
                  aes(x = avg_log2FC, y = -log10(p_val), color = significance)) +
  geom_point(alpha = 0.9, size = 3.5) +  
  
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 7,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "gray60"
  ) +
  
  geom_hline(yintercept = -log10(5e-7), linetype = "dashed", color = "darkgray", size = 0.6) +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "darkgray", size = 0.6) +
  
  scale_color_manual(
    values = c(
      "Upregulated" = "#D73027",
      "Downregulated" = "#4575B4",
      "Not significant" = "lightgray"
    )
  ) +
  
  labs(
    title = "Stromal compartment R1 vs R0",
    x = "Log2 Fold Change",
    y = "-Log10 p-Value"
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    legend.title = element_blank(),
    legend.position = "none",  
    plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5),
    axis.title = element_text(size = rel(1.2), face = "bold"),
    axis.text = element_text(size = rel(1.2)),
    panel.grid.minor = element_blank()
  )

pdf(file.path(figures.dir, "Fig4a_volcano.pdf"), width = 9, height = 8)
print(volcano)
dev.off()


## ---- S8a: up-regulated genes, all meta-spots ----

seurat_obj_subset <- readRDS(file.path(harmony.dir, "seurat_object_res_10.rds"))
DEG=c("TMEM159","TFAP2A","KRT5","APLP2","CTXN1","CAPN2","DDX3Y","TBL1X","RASAL2","NT5DC2","CORO1C", 
      "BOK","AHNAK","LRP10","ANXA1","NSG1","AZIN2","PPFIBP1","NFIC","HLA-A")
data_exp <- FetchData(seurat_obj_subset, vars = c("sample", "CNV_fraction","recurr_cens", "mmr_status", "merged_meta_clusters", DEG))
df_plot_sig=data_exp

df_plot_sig$recurr_cens <- as.factor(df_plot_sig$recurr_cens)

top_genes=c("TMEM159","APLP2","CTXN1","CAPN2","DDX3Y","TBL1X","RASAL2","NT5DC2","CORO1C","BOK","AHNAK",  
            "LRP10","ANXA1","AZIN2","PPFIBP1","NFIC","HLA-A")

select_gene=setdiff(top_genes,"HLA-A")

df_long <- df_plot_sig %>%
  select(recurr_cens, all_of(select_gene)) %>%
  pivot_longer(cols = all_of(select_gene), names_to = "gene", values_to = "value")


ymin <- -0.5
ymax <- 1.5  

counts <- df_plot_sig %>%
  ungroup() %>%
  mutate(rc = as.character(recurr_cens)) %>%
  count(rc, name = "n") %>%
  tidyr::complete(rc = c("0","1"), fill = list(n = 0))

n_non <- counts$n[counts$rc == "0"]
n_rec <- counts$n[counts$rc == "1"]

# plot all metaspots
Combined_Boxplot<-ggplot(df_long, aes(x = gene, y = value, fill = recurr_cens)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("0" = "#9197FF", "1" = "#FC7171"), 
                    labels = c("0" = paste0("No Relapse n = ",n_non),
                               "1" = paste0("Relapse n = ",n_rec))) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  theme_bw() +
  labs(title="" , x = "", y = "Scaled Expression - All meta-spots", fill = "Recurrence Status") +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = rel(1.6), color="black"),
        axis.text.y = element_text(face = "bold", size = rel(1.6), color="black"),
        axis.title.y = element_text(face = "bold", size = rel(1.8)),
        plot.title = element_text(face = "bold", size = rel(1.6), hjust = 0.5),
        legend.title=element_text(face = "bold", size = rel(1.4)),
        legend.text=element_text(face = "bold", size = rel(1.4)),
        legend.position="top")
ggsave(file.path(figures.dir, "FigS8a_DEG_boxplot_all_metaspots.pdf"), plot = Combined_Boxplot, width = 12, height = 7)



## ---- S8c, S8f: signatures, all vs stromal meta-spots ----

# signatures
load(file.path(data.dir, "signatures_fig3.RData"))

library(GSVA)

seurat_obj <- readRDS(file.path(harmony.dir, "seurat_object_res_10.rds"))

Sig_all_signatures <- ssgseaParam(as.matrix(seurat_obj@assays$RNA$scale.data), signatures, normalize = TRUE)
signat <- gsva(Sig_all_signatures)
sig <- as.data.frame(t(signat))

data_exp <- FetchData(seurat_obj, vars = c("sample", "CNV_fraction","recurr_cens", "mmr_status", "merged_meta_clusters"))
df_plot_sig=merge(sig, data_exp, by=0)
df_plot_sig$cluster = df_plot_sig$merged_meta_clusters


df_recoded <- df_plot_sig %>%
  mutate(cluster_recoded = "All_clusters") %>%
  bind_rows(
    df_plot_sig %>%
      filter(cluster %in% c(2, 4)) %>%
      mutate(cluster_recoded = as.character(cluster))
  )


df_recoded$cluster_recoded=factor(df_recoded$cluster_recoded, levels=c("All_clusters", "2","4"))
df_recoded$recurr_cens <- factor(df_recoded$recurr_cens, levels = c("0", "1"))


colnames(df_recoded)
output_dir = file.path(figures.dir, "FigS8cf_signatures_R0_R1")
dir.create(output_dir, showWarnings = FALSE)

df_recoded <- df_recoded %>%
  mutate(cluster_recoded_bis = if_else(cluster_recoded %in% c("2", "4"),
                                       "pooled_stroma",
                                       cluster_recoded))

# names as stored in signatures_fig3.RData; display labels give the references
select_sig <- c("RSC", "YAP", "Gregorieff_YAP_UP", "Gregorieff_YAP_DOWN", "Ogden_REC", "Ogden_iREC")
sig_labels <- c(RSC = "Ayyaz_RSC", YAP = "Serra_YAP", Gregorieff_YAP_UP = "Gregorieff_YAP_UP",
                Gregorieff_YAP_DOWN = "Gregorieff_YAP_DOWN", Ogden_REC = "Ogden_REC", Ogden_iREC = "Ogden_iREC")

for (signature_id in select_sig) {
  
  gene_values <- df_recoded[[signature_id]]
  Q1 <- quantile(gene_values, 0.25, na.rm = TRUE)
  Q3 <- quantile(gene_values, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  y_min <- Q1 - 1.5 * IQR_val
  y_max <- Q3 + 1.5 * IQR_val
  y_range <- y_max - y_min
  y_limits <- c(y_min - 0.05 * y_range, y_max + 0.15 * y_range) 
  
  pval_df <- df_recoded %>%
    group_by(cluster_recoded_bis) %>%
    summarise(
      pval = tryCatch(
        wilcox.test(.data[[signature_id]] ~ recurr_cens)$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(
      pval = as.numeric(pval),
      p_label = case_when(
        pval < 1e-4 ~ "italic(P) < '0.0001'",
        pval < 1e-3 ~ "italic(P) < '0.001'",
        !is.na(pval) ~ paste0("italic(P) == '", formatC(pval, format = "f", digits = 3), "'"),
        TRUE ~ "italic(P) == NA"
      ),
      y_pos = y_max + (0.07 * y_range) 
    )
  
  Combined_Boxplot <- ggplot(
    df_recoded,
    aes(
      x = factor(cluster_recoded_bis),
      y = .data[[signature_id]],
      fill = recurr_cens)
  ) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 0.6) +
    scale_fill_manual(values = c("0" = "#9197FF", "1" = "#FC7171"), 
                      labels = c("0" = "No Relapse", "1" = "Relapse"))+
    scale_x_discrete(
      labels = c(
        "All_clusters" = "All metaspots",
        "pooled_stroma" = "Stromal metaspots"))+
    
    labs(title = sig_labels[[signature_id]],
         x = "",
         y = "Enrichment score",
         fill="") +
    theme_bw() +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.02, 0.05))) +
    
    geom_text(
      data = pval_df,
      aes(
        x = factor(cluster_recoded_bis),
        y = y_pos,
        label = p_label
      ),
      parse = TRUE,
      size = 7, #5
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    
    theme(
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(face = "bold",size = rel(1.6), color = "black"),
      axis.text.y = element_text(face = "bold", size = rel(1.6), color = "black"),
      axis.title.y = element_text(face = "bold", size = rel(1.8)),
      plot.title = element_text(face = "bold", size = rel(1.8) ),
      legend.title = element_text(face = "bold", size = rel(1.4)),
      legend.text = element_text(face = "bold", size = rel(1.4)),
      legend.position = "top"
    )
  
  
  pdf_file <- file.path(output_dir, paste0("FigS8cf_", sig_labels[[signature_id]], "_R0_R1_all_vs_stroma.pdf"))
  ggsave(pdf_file, plot = Combined_Boxplot, width = 5, height = 5)
    
}
