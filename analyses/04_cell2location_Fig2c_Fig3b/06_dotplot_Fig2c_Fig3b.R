# =====================================================================================================================
# Fig. 2c and Fig. 3b - Cell2location (98 Pelka cell states) composition of meta-clusters and stromal sub-clusters
#
# For each meta-cluster / stromal sub-cluster and each of the 98 categories:
#   mean_proportion = mean proportion over all spots of the group;
#   pct_present     = % of spots of the group with proportion > 2% (fixed absolute threshold).
# Categories kept if mean proportion >= 1% in at least one group (meta-cluster or sub-cluster). Dot size = pct_present,
# dot colour = z-score of mean_proportion across groups (column-wise), so colours show relative enrichment.
# Inputs : <C2L_DIR>/tables/c2l_98cat_proportions_<Sample>.csv (05_extract_proportions.py),
#          Tables/spot_labels_all48.rds (R/Step5_Add_meta.R), data/pelka_hierarchy_87_to_17.csv (Pelka cell-state families)
# Outputs: Figures/Fig2c_metacluster_C2L_dotplot.pdf, Figures/Fig3b_stromal_subcluster_C2L_dotplot.pdf (+ tables)
# Usage  : Rscript analyses/04_cell2location_Fig2c_Fig3b/06_dotplot_Fig2c_Fig3b.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

C98DIR    <- file.path(Sys.getenv("CRCVISIUM_C2L", "~/data/cell2location"), "tables")
HIER_FILE <- file.path(data.dir, "pelka_hierarchy_87_to_17.csv")

meta_all <- readRDS(file.path(tables.dir, "spot_labels_all48.rds"))

join_list <- lapply(Sample.names, function(s) {
  c2l_df <- read.csv(file.path(C98DIR, paste0("c2l_98cat_proportions_", s, ".csv")), row.names = 1, check.names = FALSE)
  c2l_df$barcode <- rownames(c2l_df)
  merge(meta_all[meta_all$sample == s, ], c2l_df, by = "barcode")
})
final <- bind_rows(join_list)
rm(join_list); gc()
cell_types98 <- setdiff(colnames(final), c("barcode", "sample", "merged_meta_clusters", "stromal_subcluster"))
message(nrow(final), " spots, ", length(cell_types98), " categories")

# ---- 98 categories -> cell family (17 groups), from the Pelka et al. hierarchy ----
hier <- read.csv(HIER_FILE, stringsAsFactors = FALSE)
hier$sub_clean <- gsub("/", "-", hier$cl295v11SubFull, fixed = TRUE)

fam17_of_row <- function(top, mid) {
  if (top == "Epi") return("Normal epithelial")
  if (top == "Mast") return("Mast")
  if (top == "B") return("B")
  if (top == "Plasma") return("Plasma")
  if (top == "Myeloid") {
    if (mid == "DC") return("Dendritic cell")
    if (mid == "Granulo") return("Granulocyte")
    if (mid %in% c("Macro", "Mono")) return("Macrophage-Monocyte")
  }
  if (top == "Strom") {
    if (mid == "Endo") return("Endothelial")
    if (mid == "Fibro") return("Fibroblast-CAF")
    if (mid == "Peri") return("Pericyte")
    if (mid == "Schwann") return("Schwann")
    if (mid == "SmoothMuscle") return("Smooth muscle")
  }
  if (top == "TNKILC") {
    if (mid == "TCD4") return("CD4 T")
    if (mid == "TCD8") return("CD8 T")
    if (mid == "NK") return("NK")
    if (mid %in% c("ILC", "Tgd", "TZBTB16")) return("T-NK-ILC-Other")
  }
  NA_character_
}
hier$family17 <- mapply(fam17_of_row, hier$clTopLevel, hier$clMidwayPr)
fam_lookup <- setNames(hier$family17, hier$sub_clean)
family98 <- fam_lookup[cell_types98]
family98[grepl("^Tumor cE", cell_types98)] <- "Malignant epithelial"
names(family98) <- cell_types98
stopifnot(!anyNA(family98))

family_order <- c("Malignant epithelial", "Normal epithelial", "Fibroblast-CAF", "Endothelial", "Pericyte",
                   "Smooth muscle", "Schwann", "CD8 T", "NK", "B", "Plasma", "Dendritic cell", "CD4 T",
                   "T-NK-ILC-Other", "Macrophage-Monocyte", "Granulocyte", "Mast")

# ---- aggregation per group ----
THRESH <- 0.02
agg_group <- function(data, group_col) {
  data %>%
    filter(!is.na(.data[[group_col]])) %>%
    pivot_longer(all_of(cell_types98), names_to = "cell_type", values_to = "proportion") %>%
    group_by(.data[[group_col]], cell_type) %>%
    summarise(mean_proportion = mean(proportion, na.rm = TRUE),
              pct_present = 100 * mean(proportion > THRESH, na.rm = TRUE),
              n_spots = n(), .groups = "drop") %>%
    rename(group = 1)
}
agg_mc <- agg_group(final, "merged_meta_clusters")
agg_st <- agg_group(final, "stromal_subcluster")

# ---- columns kept: mean proportion >= 1% in at least one group ----
KEEP_THRESH <- 0.01
max_mc <- agg_mc %>% group_by(cell_type) %>% summarise(max_mc = max(mean_proportion), .groups = "drop")
max_st <- agg_st %>% group_by(cell_type) %>% summarise(max_st = max(mean_proportion), .groups = "drop")
notable <- full_join(max_mc, max_st, by = "cell_type") %>%
  mutate(max_st = ifelse(is.na(max_st), 0, max_st), max_mc = ifelse(is.na(max_mc), 0, max_mc))

kept_types <- notable %>% filter(max_mc >= KEEP_THRESH | max_st >= KEEP_THRESH) %>% pull(cell_type)
message("Categories kept: ", length(kept_types), " / ", length(cell_types98))

build_col_order <- function(order_by) {
  base <- data.frame(cell_type = kept_types) %>%
    mutate(family = family98[cell_type]) %>%
    mutate(family = factor(family, levels = family_order)) %>%
    left_join(notable, by = "cell_type")
  if (order_by == "numeric") {
    base <- base %>% mutate(num_key = as.numeric(gsub(".*?([0-9]+).*", "\\1", cell_type))) %>%
      arrange(family, num_key)
  } else {
    base <- base %>% mutate(rank_val = pmax(max_mc, max_st)) %>%
      arrange(family, desc(rank_val))
  }
  col_levels <- base$cell_type
  fam_seq <- as.character(base$family)
  gap_positions <- which(fam_seq[-1] != fam_seq[-length(fam_seq)]) + 0.5
  fam_blocks <- base %>% mutate(pos = row_number()) %>%
    group_by(family) %>% summarise(xmid = mean(pos), .groups = "drop")
  list(col_levels = col_levels, gap_positions = gap_positions, fam_blocks = fam_blocks)
}
ord_numeric <- build_col_order("numeric")

# names and display order from R/00_config.R (meta-cluster 8 = Tumor-fibroblasts)
SC_NAMES <- stromal_names
MC_NAMES <- meta_names
MC_ORDER <- rev(c("0", "1", "7", "8", "2", "4", "3", "5", "6", "9", "10"))
SC_ORDER <- rev(c("0", "2", "1", "3", "4"))

make_dotplot <- function(agg, group_label, scale_by, title, out_basename, group_names = NULL, group_order = NULL, h_base = 1, h_per_row = 0.8, ord = ord_numeric) {
  col_levels <- ord$col_levels; gap_positions <- ord$gap_positions; fam_blocks <- ord$fam_blocks
  d <- agg %>% filter(cell_type %in% kept_types) %>%
    mutate(cell_type = factor(cell_type, levels = col_levels))
  group_levels <- if (!is.null(group_order)) group_order else sort(unique(as.integer(as.character(d$group))))
  d$group <- factor(as.character(d$group), levels = as.character(group_levels))
  if (!is.null(group_names)) levels(d$group) <- group_names[levels(d$group)]

  if (scale_by == "column") {
    d <- d %>% group_by(cell_type) %>% mutate(z = as.numeric(scale(mean_proportion))) %>% ungroup()
  } else {
    d <- d %>% group_by(group) %>% mutate(z = as.numeric(scale(mean_proportion))) %>% ungroup()
  }
  d$z[is.na(d$z)] <- 0
  zmax <- max(abs(d$z), na.rm = TRUE)
  if (!is.finite(zmax) || zmax == 0) zmax <- 1

  p <- ggplot(d, aes(x = cell_type, y = group)) +
    geom_vline(xintercept = gap_positions, color = "grey75", linewidth = 0.3) +
    geom_point(aes(size = pct_present, fill = z), shape = 21, color = "black", stroke = 0.3) +
    geom_text(data = fam_blocks, aes(x = xmid, y = Inf, label = family), inherit.aes = FALSE,
              angle = 45, hjust = 0, vjust = 0, size = 2.6, color = "grey30") +
    coord_cartesian(clip = "off") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0,
                          limits = c(-zmax, zmax), name = paste0("Z-score\n(", scale_by, ")")) +
    scale_size_continuous(name = paste0("% spots\n(prop.>", THRESH * 100, "%)"),
                           range = c(0.3, 7), limits = c(0, 100)) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(t = 70, r = 60, b = 5.5, l = 5.5),
          legend.box.spacing = unit(1.6, "cm")) +
    labs(x = NULL, y = NULL)

  w <- max(14, length(kept_types) * 0.22)
  h <- h_base + h_per_row * length(group_levels)
  ggsave(file.path(figures.dir, paste0(out_basename, ".pdf")), p, width = w, height = h, bg = "white")
  write.csv(d, file.path(tables.dir, paste0(out_basename, ".csv")), row.names = FALSE)
}

make_dotplot(agg_mc, "merged_meta_clusters", "column", "", "Fig2c_metacluster_C2L_dotplot",
             group_names = MC_NAMES, group_order = MC_ORDER, h_base = 0.3, h_per_row = 0.55, ord = ord_numeric)
make_dotplot(agg_st, "stromal_subcluster", "column", "", "Fig3b_stromal_subcluster_C2L_dotplot",
             group_names = SC_NAMES, group_order = SC_ORDER, h_base = 1, h_per_row = 0.8, ord = ord_numeric)
