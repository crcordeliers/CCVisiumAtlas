#!/usr/bin/env Rscript
# =============================================================================
# benchmark_clustering.R
# =============================================================================
# Runs the two R-based alternative spatial clustering methods (BayesSpace,
# BASS) and then benchmarks all five clustering methods -- Louvain (reference,
# assumed already computed), BayesSpace, BASS, CellCharter, GraphST (the last
# two produced by the companion script `benchmark_clustering.py`) -- on
# multi-sample 10x Visium data, across five complementary axes:
#   1. Spatial coherence (Moran's I per cluster)
#   2. Agreement with independent pathology ground truth (ARI, NMI)
#   3. Cross-patient generalizability (patient-mixing entropy)
#   4. Cluster purity vs. pathology
#   5. Pairwise concordance between the five clustering methods themselves
#
# Package versions used in the reference run:
#   BayesSpace 1.20.2, harmony 2.0.5, R 4.5.2 (BayesSpace + BASS)
#   BASS       1.1.0.17
#   mclust     6.1.3, aricode 1.1.0, spdep >= 1.3
#
# Usage:
#   Rscript benchmark_clustering.R --step bayesspace
#   Rscript benchmark_clustering.R --step bass
#   Rscript benchmark_clustering.R --step metrics
#
# Each Visium sample is expected as a standard SpaceRanger output directory
# (`<DATA_ROOT>/<sample>/outs/`).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

# -----------------------------------------------------------------------
# Configuration -- EDIT THESE PATHS for your own environment
# -----------------------------------------------------------------------
# Run from the repository root. Paths can be set with the same environment variables as R/00_config.R.
DATA_ROOT <- Sys.getenv("CRCVISIUM_RAW", "/path/to/visium_spaceranger_outputs")   # one subdir per sample
OUT_DIR <- file.path(Sys.getenv("CRCVISIUM_WORKDIR", "."), "clustering_benchmark")
FIG_DIR <- file.path(OUT_DIR, "figures")
SAMPLES_FILE <- Sys.getenv("CRCVISIUM_SAMPLES", "data/samples_48.txt")             # one sample name per line
SEED <- 42                                            # BayesSpace
BASS_SEED <- 0                                        # BASS
N_CLUSTERS <- 11                                      # matched to the Louvain reference

# Per-spot cluster label files expected as input to the metrics step (one
# CSV per method, columns: sample, barcode, <method>_cluster). Louvain and
# pathology ground truth are external inputs to this script (not computed
# here); CellCharter/GraphST are produced by benchmark_clustering.py.
LOUVAIN_LABELS_CSV <- file.path(OUT_DIR, "louvain_labels_all48.csv")                          # sample, barcode, merged_meta_clusters
PATHOLOGY_CSV <- file.path(OUT_DIR, "pathology_annotation_48samples.csv")                    # sample, barcode, pathology_full, pathology_collapsed
BAYESSPACE_LABELS_CSV <- file.path(OUT_DIR, "spot_bayesspace_cluster_all48.csv")
BASS_LABELS_CSV <- file.path(OUT_DIR, "spot_bass_cluster_all48.csv")
CELLCHARTER_LABELS_CSV <- file.path(OUT_DIR, "spot_cellcharter_cluster_all48.csv")
GRAPHST_LABELS_CSV <- file.path(OUT_DIR, "spot_graphst_cluster_harmony_all48.csv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Space Ranger >= 2 writes tissue_positions.csv (with header); 1.x writes
# tissue_positions_list.csv (no header, same columns)
read_positions <- function(spatial_dir) {
  f <- file.path(spatial_dir, "tissue_positions.csv")
  if (file.exists(f)) return(read.csv(f))
  read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), header = FALSE,
           col.names = c("barcode", "in_tissue", "array_row", "array_col",
                         "pxl_row_in_fullres", "pxl_col_in_fullres"))
}

load_samples <- function() readLines(SAMPLES_FILE)

# =============================================================================
# 1. BayesSpace
# =============================================================================
run_bayesspace <- function() {
  suppressPackageStartupMessages({
    library(BayesSpace)
    library(SingleCellExperiment)
    library(harmony)
    library(Matrix)
    library(Seurat)
  })
  samples <- load_samples()
  set.seed(SEED)  # k-means init (Harmony, BayesSpace) + MCMC

  sce_list <- lapply(samples, function(s) {
    counts <- Read10X_h5(file.path(DATA_ROOT, s, "outs", "filtered_feature_bc_matrix.h5"))
    coords <- read_positions(file.path(DATA_ROOT, s, "outs", "spatial"))
    coords <- coords[match(colnames(counts), coords$barcode), ]
    SingleCellExperiment(
      assays = list(counts = counts),
      colData = DataFrame(sample = s, barcode = colnames(counts),
                           array_row = coords$array_row, array_col = coords$array_col)  # BayesSpace >= 1.11 names
    )
  })
  sce_combined <- do.call(cbind, sce_list)

  # Joint analysis: all 48 samples share the same array grid, so without an
  # offset .find_neighbors() would link spots of different patients at the same
  # position. Each sample is shifted by 200 array rows (even shift keeps the
  # hexagonal row/col parity).
  idx <- match(sce_combined$sample, samples) - 1L
  sce_combined$array_row <- sce_combined$array_row + 200L * idx

  # log-normalize + PCA (n.PCs=50), then Harmony batch correction (batch = sample)
  sce_combined <- spatialPreprocess(sce_combined, n.PCs = 50, platform = "Visium")
  sce_combined <- RunHarmony(sce_combined, "sample", reduction.save = "HARMONY")

  # spatial clustering on the Harmony-corrected embedding
  sce_combined <- spatialCluster(sce_combined, use.dimred = "HARMONY", q = N_CLUSTERS,
                                  nrep = 10000, platform = "Visium")

  out <- as.data.frame(colData(sce_combined)) %>%
    select(sample, barcode, bayesspace_cluster = spatial.cluster)
  write_csv(out, BAYESSPACE_LABELS_CSV)
  message("Wrote: ", BAYESSPACE_LABELS_CSV)
}

# =============================================================================
# 2. BASS
# =============================================================================
run_bass <- function(C = 20, R = N_CLUSTERS) {
  suppressPackageStartupMessages({
    library(BASS)
    library(Seurat)
    library(Matrix)
  })
  samples <- load_samples()
  set.seed(BASS_SEED)  # mclust init + MCMC

  cntm <- list(); xym <- list()
  for (s in samples) {
    counts <- Read10X_h5(file.path(DATA_ROOT, s, "outs", "filtered_feature_bc_matrix.h5"))
    coords <- read_positions(file.path(DATA_ROOT, s, "outs", "spatial"))
    coords <- coords[match(colnames(counts), coords$barcode), ]
    cntm[[s]] <- counts
    xym[[s]] <- as.matrix(coords[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
    rownames(xym[[s]]) <- coords$barcode  # BASS checks names against colnames(counts)
  }

  # BASS object: C = cell/domain types (package vignette default), R = spatial
  # domains (matched to the Louvain reference cluster count)
  BASS_obj <- createBASSObject(cntm, xym, C = C, R = R,
                                beta_method = "SW", init_method = "mclust")
  BASS_obj <- BASS.preprocess(BASS_obj, doLogNormalize = TRUE, geneSelect = "sparkx",
                               nSE = 3000, doPCA = TRUE, scaleFeature = FALSE, nPC = 20)
  # MCMC at package defaults (burnin=2000, nsample=5000)
  BASS_obj <- BASS.run(BASS_obj)
  BASS_obj <- BASS.postprocess(BASS_obj)

  zlabels_list <- BASS_obj@results$z
  out <- bind_rows(lapply(seq_along(samples), function(i) {
    data.frame(sample = samples[i], barcode = colnames(cntm[[samples[i]]]),
               bass_cluster = zlabels_list[[i]])
  }))
  write_csv(out, BASS_LABELS_CSV)
  message("Wrote: ", BASS_LABELS_CSV)
}

# =============================================================================
# 3. Comparison metrics (all five methods)
# =============================================================================
run_metrics <- function() {
  suppressPackageStartupMessages({
    library(mclust)   # adjustedRandIndex
    library(aricode)  # NMI
    library(spdep)    # Moran's I
  })

  METHOD_LEVELS <- c("Louvain", "BayesSpace", "BASS", "CellCharter", "GraphST")

  louvain <- read_csv(LOUVAIN_LABELS_CSV, show_col_types = FALSE) %>%
    rename(cluster = merged_meta_clusters) %>% mutate(method = "Louvain")
  bayesspace <- read_csv(BAYESSPACE_LABELS_CSV, show_col_types = FALSE) %>%
    rename(cluster = bayesspace_cluster) %>% mutate(method = "BayesSpace")
  bass <- read_csv(BASS_LABELS_CSV, show_col_types = FALSE) %>%
    rename(cluster = bass_cluster) %>% mutate(method = "BASS")
  cellcharter <- read_csv(CELLCHARTER_LABELS_CSV, show_col_types = FALSE) %>%
    rename(cluster = cellcharter_cluster) %>% mutate(method = "CellCharter") %>%
    mutate(barcode = sub("^(.*-1)-[0-9]+$", "\\1", barcode))  # drop anndata concat suffix (-<batch>)
  graphst <- read_csv(GRAPHST_LABELS_CSV, show_col_types = FALSE) %>%
    rename(cluster = graphst_cluster_harmony) %>% mutate(method = "GraphST")

  all_labels <- bind_rows(louvain, bayesspace, bass, cellcharter, graphst) %>%
    select(sample, barcode, cluster, method)
  pathology <- read_csv(PATHOLOGY_CSV, show_col_types = FALSE)

  # ---- wide table: one column per method, for spot-level joins ----
  wide <- all_labels %>%
    pivot_wider(names_from = method, values_from = cluster) %>%
    inner_join(pathology, by = c("sample", "barcode"))

  samples <- unique(wide$sample)

  # spot array coordinates for the spatial weights (Moran's I)
  positions <- bind_rows(lapply(samples, function(s)
    read_positions(file.path(DATA_ROOT, s, "outs", "spatial")) %>%
      transmute(sample = s, barcode, row = array_row, col = array_col)))
  wide <- wide %>% inner_join(positions, by = c("sample", "barcode"))
  message("Spots in metrics table: ", nrow(wide), "; missing labels per method: ",
          paste(sapply(METHOD_LEVELS, function(m) sum(is.na(wide[[m]]))), collapse = "/"))

  # -------------------------------------------------------------------
  # 1. Moran's I per cluster per sample (spatial coherence)
  # -------------------------------------------------------------------
  moransi_list <- list()
  for (m in METHOD_LEVELS) {
    for (s in samples) {
      ws <- wide[wide$sample == s, ]
      coords <- as.matrix(ws[, c("row", "col")])  # expects row/col columns joined upstream
      nb <- knearneigh(coords, k = 6) %>% knn2nb()
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
      for (cl in unique(na.omit(ws[[m]]))) {
        ind <- as.numeric(ws[[m]] == cl)
        mt <- tryCatch(moran.mc(ind, lw, nsim = 999, zero.policy = TRUE),
                        error = function(e) NULL)
        if (!is.null(mt)) {
          moransi_list[[length(moransi_list) + 1]] <- data.frame(
            method = m, sample = s, cluster = cl, MoranI = mt$statistic,
            p_value = mt$p.value, n_spots = sum(ind))
        }
      }
    }
  }
  moransi <- bind_rows(moransi_list)
  write_csv(moransi, file.path(OUT_DIR, "moransi_5methods.csv"))

  # -------------------------------------------------------------------
  # 2. ARI / NMI vs. pathology ground truth, per sample
  # -------------------------------------------------------------------
  ari_nmi_path <- bind_rows(lapply(METHOD_LEVELS, function(m) {
    bind_rows(lapply(samples, function(s) {
      ws <- wide[wide$sample == s, ]
      data.frame(sample = s, method = m,
                 ARI = adjustedRandIndex(ws[[m]], ws$pathology_full),
                 NMI = NMI(factor(ws[[m]]), factor(ws$pathology_full), variant = "sqrt"))  # labels as factors (aricode); NMI normalised by sqrt(H(U) H(V))
    }))
  }))
  write_csv(ari_nmi_path, file.path(OUT_DIR, "ari_nmi_vs_pathology_5methods.csv"))

  # -------------------------------------------------------------------
  # 3. Patient-mixing entropy per cluster
  # -------------------------------------------------------------------
  entropy_norm <- function(counts) {
    p <- counts / sum(counts)
    h <- -sum(p * log(p + 1e-12))
    h / log(length(counts))  # normalized to [0,1] by the max possible entropy
  }
  mixing <- all_labels %>%
    group_by(method, cluster) %>%
    summarise(n_samples = n_distinct(sample), n_spots = n(),
              shannon_entropy_norm = entropy_norm(table(sample)), .groups = "drop")
  write_csv(mixing, file.path(OUT_DIR, "patient_mixing_5methods.csv"))

  # -------------------------------------------------------------------
  # 4. Cluster purity vs. pathology (weighted mean dominant-category fraction)
  # -------------------------------------------------------------------
  purity <- bind_rows(lapply(METHOD_LEVELS, function(m) {
    wide %>% group_by(cluster = .data[[m]]) %>%
      summarise(n_spots = n(),
                purity = max(table(pathology_full)) / n(), .groups = "drop") %>%
      mutate(method = m)
  }))
  write_csv(purity, file.path(OUT_DIR, "purity_per_cluster_5methods.csv"))

  # -------------------------------------------------------------------
  # 5. Pairwise ARI / NMI between the five methods
  # -------------------------------------------------------------------
  pairs <- combn(METHOD_LEVELS, 2, simplify = FALSE)
  pairwise <- bind_rows(lapply(pairs, function(pr) {
    bind_rows(lapply(samples, function(s) {
      ws <- wide[wide$sample == s, ]
      data.frame(sample = s, method_1 = pr[1], method_2 = pr[2],
                 ARI = adjustedRandIndex(ws[[pr[1]]], ws[[pr[2]]]),
                 NMI = NMI(factor(ws[[pr[1]]]), factor(ws[[pr[2]]]), variant = "sqrt"))
    }))
  }))
  write_csv(pairwise, file.path(OUT_DIR, "ari_nmi_pairwise_5methods_per_sample.csv"))

  pairwise_summary <- pairwise %>% group_by(method_1, method_2) %>%
    summarise(median_ARI = median(ARI, na.rm = TRUE),
              median_NMI = median(NMI, na.rm = TRUE), .groups = "drop")
  write_csv(pairwise_summary, file.path(OUT_DIR, "ari_nmi_pairwise_5methods_summary.csv"))

  message("All metrics computed and written to ", OUT_DIR)
  invisible(list(moransi = moransi, ari_nmi_path = ari_nmi_path, mixing = mixing,
                  purity = purity, pairwise = pairwise_summary))
}

# =============================================================================
# 4. Figures (benchmark panel + pairwise ARI/NMI heatmap)
# =============================================================================
run_figures <- function() {
  suppressPackageStartupMessages({ library(patchwork); library(ggnewscale) })
  METHOD_LEVELS <- c("Louvain", "BayesSpace", "BASS", "CellCharter", "GraphST")
  METHOD_COLORS <- c("Louvain" = "#4E79A7", "BayesSpace" = "#F28E2B", "BASS" = "#59A14F",
                      "CellCharter" = "#E15759", "GraphST" = "#B07AA1")

  moransi <- read_csv(file.path(OUT_DIR, "moransi_5methods.csv"), show_col_types = FALSE)
  mixing <- read_csv(file.path(OUT_DIR, "patient_mixing_5methods.csv"), show_col_types = FALSE)
  ari_nmi_path <- read_csv(file.path(OUT_DIR, "ari_nmi_vs_pathology_5methods.csv"), show_col_types = FALSE)
  purity <- read_csv(file.path(OUT_DIR, "purity_per_cluster_5methods.csv"), show_col_types = FALSE)
  summ <- read_csv(file.path(OUT_DIR, "ari_nmi_pairwise_5methods_summary.csv"), show_col_types = FALSE)

  for (df in list(moransi, mixing, ari_nmi_path, purity)) df$method <- factor(df$method, levels = METHOD_LEVELS)

  panel <- function(data, yvar, title) {
    ggplot(data, aes(x = method, y = .data[[yvar]], fill = method, color = method)) +
      geom_boxplot(alpha = 0.35, outlier.shape = NA, color = "grey30", width = 0.65) +
      geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
      scale_fill_manual(values = METHOD_COLORS, guide = "none") +
      scale_color_manual(values = METHOD_COLORS, guide = "none") +
      labs(title = title, x = NULL, y = NULL) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
  }
  combined <- (panel(moransi, "MoranI", "Moran's I") |
               panel(mixing, "shannon_entropy_norm", "Patient mixing\n(entropy)") |
               panel(ari_nmi_path, "ARI", "ARI\n(vs pathology)") |
               panel(ari_nmi_path, "NMI", "NMI\n(vs pathology)") |
               panel(purity, "purity", "Cluster purity\n(vs pathology)")) +
    plot_annotation(title = "All metrics at a glance (5 methods)")
  ggsave(file.path(FIG_DIR, "FigS2c_all_metrics_panel_5methods.pdf"), combined, width = 18, height = 6, bg = "white")

  # pairwise ARI (lower triangle) / NMI (upper triangle) heatmap, own min-max scale each
  mat_ari <- matrix(NA, 5, 5, dimnames = list(METHOD_LEVELS, METHOD_LEVELS))
  mat_nmi <- mat_ari
  for (i in seq_len(nrow(summ))) {
    mat_ari[summ$method_1[i], summ$method_2[i]] <- summ$median_ARI[i]
    mat_ari[summ$method_2[i], summ$method_1[i]] <- summ$median_ARI[i]
    mat_nmi[summ$method_1[i], summ$method_2[i]] <- summ$median_NMI[i]
    mat_nmi[summ$method_2[i], summ$method_1[i]] <- summ$median_NMI[i]
  }
  idx <- expand.grid(i = seq_along(METHOD_LEVELS), j = seq_along(METHOD_LEVELS))
  idx$method_1 <- METHOD_LEVELS[idx$i]; idx$method_2 <- METHOD_LEVELS[idx$j]
  idx$value <- NA; idx$metric <- NA
  for (k in seq_len(nrow(idx))) {
    if (idx$i[k] == idx$j[k]) next
    if (idx$i[k] > idx$j[k]) { idx$value[k] <- mat_ari[idx$method_1[k], idx$method_2[k]]; idx$metric[k] <- "ARI" }
    else { idx$value[k] <- mat_nmi[idx$method_1[k], idx$method_2[k]]; idx$metric[k] <- "NMI" }
  }
  idx$method_1 <- factor(idx$method_1, levels = rev(METHOD_LEVELS))
  idx$method_2 <- factor(idx$method_2, levels = METHOD_LEVELS)
  idx_ari <- idx[idx$metric == "ARI" & !is.na(idx$metric), ]
  idx_nmi <- idx[idx$metric == "NMI" & !is.na(idx$metric), ]

  p_heat <- ggplot() +
    geom_tile(data = idx_ari, aes(x = method_2, y = method_1, fill = value), color = "white", linewidth = 1) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = range(idx_ari$value, na.rm = TRUE), name = "ARI") +
    new_scale_fill() +
    geom_tile(data = idx_nmi, aes(x = method_2, y = method_1, fill = value), color = "white", linewidth = 1) +
    scale_fill_gradient(low = "#FFF5EB", high = "#A63603", limits = range(idx_nmi$value, na.rm = TRUE), name = "NMI") +
    geom_text(data = idx, aes(x = method_2, y = method_1, label = ifelse(is.na(value), "", sprintf("%.3f", value))),
              size = 3.8, fontface = "bold") +
    labs(x = NULL, y = NULL) +   # Supplementary Fig. S2d
    coord_fixed() + theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
          axis.text.y = element_text(size = 14, face = "bold"), panel.grid = element_blank())
  ggsave(file.path(FIG_DIR, "FigS2d_ari_nmi_pairwise_heatmap.pdf"), p_heat, width = 8, height = 6.5, bg = "white")

  message("Figures written to ", FIG_DIR)
}

# =============================================================================
# CLI
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
step_idx <- which(args == "--step")
step <- if (length(step_idx)) args[step_idx + 1] else NA

if (is.na(step)) {
  message("Usage: Rscript benchmark_clustering.R --step {bayesspace|bass|metrics|figures}")
} else if (step == "bayesspace") {
  run_bayesspace()
} else if (step == "bass") {
  run_bass()
} else if (step == "metrics") {
  run_metrics()
} else if (step == "figures") {
  run_figures()
} else {
  stop("Unknown --step: ", step)
}
