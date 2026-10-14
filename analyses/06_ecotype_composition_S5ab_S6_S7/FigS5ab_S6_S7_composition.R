#!/usr/bin/env Rscript
# =====================================================================================================================
# Ecotype composition of the 48 tumours (proportion of spots of each meta-cluster / stromal sub-cluster per patient)
#   S5a  number of samples containing each meta-cluster + distribution of its proportion in those samples
#   S5b  proportion of each meta-cluster, dMMR vs pMMR (Wilcoxon)
#   S6a  number of samples containing each stromal sub-cluster (proportion > 1% of the stromal spots) + distribution
#   S6b  proportion of each stromal sub-cluster (among stromal spots), dMMR vs pMMR (Wilcoxon)
#   S6c  same, Relapse vs No relapse (Wilcoxon)
#   S7a,b proportions of the 4 tumour and 2 stromal meta-clusters (recomputed among tumour + stromal spots), per patient,
#        all patients (a) and pMMR (b); patients ordered by total tumour proportion
#   S7c,d pooled tumour (MC 0, 1, 7, 8) and stromal (MC 2, 4) proportions, Relapse vs No relapse (Wilcoxon), all (c) / pMMR (d)
# Input : Results_Harmony/res_10/scaled_data_list_10.RData (spot-level meta-cluster and sub-cluster labels, R/Step3-4);
#         data/data_clin_48_filtered.RDS (MMR status, relapse)
# Output: Figures/FigS5ab_S6_S7/ (one PDF per panel / ecotype) + proportion tables in Tables/
# Usage : Rscript analyses/06_ecotype_composition_S5ab_S6_S7/FigS5ab_S6_S7_composition.R
# Patients: 11 dMMR, 37 pMMR (24 No relapse, 13 Relapse).
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(ggtext) })
set.seed(1)   # fixed jitter of the points (no effect on the statistics)

load(file.path(harmony.dir, "scaled_data_list_10.RData"))
clin <- readRDS(file.path(data.dir, "data_clin_48_filtered.RDS"))
out_dir <- file.path(figures.dir, "FigS5ab_S6_S7"); dir.create(out_dir, showWarnings = FALSE)

spots <- bind_rows(lapply(names(scaled_data_list), function(s)
  data.frame(Sample = s, merged_meta_clusters = as.character(scaled_data_list[[s]]$merged_meta_clusters),
             stromal_subcluster = as.character(scaled_data_list[[s]]$stromal_subcluster)))) |>
  left_join(clin[, c("sample", "recurr_cens", "mmr_status")], by = c("Sample" = "sample")) |>
  mutate(recurr = factor(recurr_cens, levels = c(0, 1), labels = c("No Relapse", "Relapse")))

proportions <- function(d, var) {
  d |> filter(!is.na(.data[[var]])) |>
    count(Sample, recurr, mmr_status, cl = .data[[var]], name = "count") |>
    group_by(Sample) |> mutate(Proportion = 100 * count / sum(count)) |> ungroup() |>
    complete(nesting(Sample, recurr, mmr_status), cl, fill = list(count = 0, Proportion = 0))
}
pmeta <- proportions(spots, "merged_meta_clusters")
pstro <- proportions(spots, "stromal_subcluster")            # among stromal spots only
write.csv(pmeta, file.path(tables.dir, "Proportions_metaclusters_per_sample.csv"), row.names = FALSE)
write.csv(pstro, file.path(tables.dir, "Proportions_stromal_subclusters_per_sample.csv"), row.names = FALSE)

p_lab <- function(p) if (p < 1e-4) "italic(P) < '0.0001'" else if (p < 1e-3) "italic(P) < '0.001'" else
  paste0("italic(P) == '", formatC(p, format = "f", digits = 3), "'")

# boxplot of one proportion between two groups (MMR: coral / grey; relapse: blue / red), Wilcoxon P
two_group_box <- function(d, group, title, file) {
  lev  <- if (group == "mmr_status") c("dMMR", "pMMR") else c("No Relapse", "Relapse")
  fill <- if (group == "mmr_status") c(dMMR = "coral", pMMR = "grey") else
    c("No Relapse" = unname(relapse_colors["0"]), "Relapse" = unname(relapse_colors["1"]))
  d$g <- factor(d[[group]], levels = lev)
  n <- table(d$g)
  labs_x <- setNames(paste0("**", lev, "**<br><span style='font-weight:400; font-size:15pt'>(n = ", n[lev], ")</span>"), lev)
  pval <- wilcox.test(Proportion ~ g, data = d)$p.value
  p <- ggplot(d, aes(g, Proportion, fill = g)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) + geom_jitter(width = 0.2, size = 2, alpha = 0.3) +
    scale_fill_manual(values = fill) + scale_x_discrete(labels = labs_x) +
    annotate("text", x = 1.5, y = max(d$Proportion) * 1.1, label = p_lab(pval), parse = TRUE, size = 6) +
    labs(title = title, y = "Proportion") + theme_bw() +
    theme(legend.position = "none", axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_markdown(size = rel(1.8), lineheight = 1.05, margin = margin(t = 4)),
          axis.text.y = element_text(face = "bold", size = rel(1.6)), axis.title.y = element_text(face = "bold", size = rel(1.8)),
          plot.title = element_text(face = "bold", size = rel(1.6), hjust = 0.5))
  ggsave(file.path(out_dir, file), p, width = 5, height = 5)
  data.frame(panel = file, P = pval)
}

# number of samples containing each cluster (bars) + per-sample proportions in those samples (dots, sorted)
prevalence <- function(d, order, colors, threshold, prefix) {
  cnt <- d |> filter(Proportion > threshold) |> count(cl) |> right_join(data.frame(cl = order), by = "cl") |>
    mutate(n = coalesce(n, 0L), cl = factor(cl, levels = order))
  p <- ggplot(cnt, aes(cl, n, fill = cl)) + geom_col(width = 0.6) + scale_fill_manual(values = colors) +
    labs(x = "", y = "Number of samples") + theme_bw() +
    theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_blank())
  ggsave(file.path(out_dir, paste0(prefix, "_number_of_samples.pdf")), p, width = 8, height = 3)
  y_max <- max(d$Proportion); max_n <- max(cnt$n)
  for (k in order) {
    dk <- d |> filter(cl == k, Proportion > 0) |> arrange(Proportion) |>
      mutate(x_pos = seq_len(n()) + (max_n - n()) / 2)
    p <- ggplot(dk, aes(x_pos, Proportion)) + geom_point(color = colors[[k]], size = 2) +
      xlim(0, max_n) + ylim(0, y_max) + labs(x = "", y = "Proportion") + theme_bw() +
      theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_blank())
    ggsave(file.path(out_dir, paste0(prefix, "_distribution_", k, ".pdf")), p, width = 2, height = 4)
  }
}

stats <- list()
## ---- S5a, S5b: meta-clusters ----
prevalence(pmeta, meta_order, meta_colors, 0, "FigS5a")
for (k in meta_order) stats[[length(stats) + 1]] <- two_group_box(filter(pmeta, cl == k), "mmr_status", meta_names[[k]],
                                                                  paste0("FigS5b_MMR_", k, ".pdf"))
## ---- S6a-c: stromal sub-clusters ----
prevalence(pstro, stromal_order, stromal_colors, 1, "FigS6a")
for (k in stromal_order) {
  stats[[length(stats) + 1]] <- two_group_box(filter(pstro, cl == k), "mmr_status", stromal_names[[k]], paste0("FigS6b_MMR_", k, ".pdf"))
  stats[[length(stats) + 1]] <- two_group_box(filter(pstro, cl == k), "recurr", stromal_names[[k]], paste0("FigS6c_relapse_", k, ".pdf"))
}

## ---- S7: tumour and stromal meta-clusters only ----
ts_codes <- c("0", "1", "7", "8", "2", "4")
pts <- pmeta |> filter(cl %in% ts_codes) |> group_by(Sample) |> mutate(prop_ts = 100 * count / sum(count)) |> ungroup() |>
  mutate(cl = factor(cl, levels = ts_codes))
for (sub in c("all", "pMMR")) {
  d <- if (sub == "all") pts else filter(pts, mmr_status == "pMMR")
  ord <- d |> filter(cl %in% c("0", "1", "7", "8")) |> group_by(Sample) |> summarise(t = sum(prop_ts)) |> arrange(t) |> pull(Sample)
  d$Sample <- factor(d$Sample, levels = ord)
  p <- ggplot(d, aes(Sample, prop_ts, fill = cl)) + geom_col(position = "stack") +
    scale_fill_manual(values = meta_colors) + facet_grid(. ~ recurr, scales = "free_x", space = "free_x") +
    labs(y = "Proportion") + theme_bw() +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size = 16), axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), panel.grid = element_blank(),
          strip.text.x = element_text(size = 14, face = "bold"))
  ggsave(file.path(out_dir, paste0("FigS7", if (sub == "all") "a" else "b", "_proportions_tumor_stroma_", sub, ".pdf")), p, width = 8, height = 6)
  pooled <- d |> mutate(pool = ifelse(cl %in% c("2", "4"), "Stroma", "Tumor")) |>
    group_by(Sample, recurr, pool) |> summarise(Proportion = sum(prop_ts), .groups = "drop")
  for (pl in c("Tumor", "Stroma"))
    stats[[length(stats) + 1]] <- two_group_box(filter(pooled, pool == pl), "recurr", pl,
                                                paste0("FigS7", if (sub == "all") "c" else "d", "_", pl, "_", sub, ".pdf"))
}
write.csv(bind_rows(stats), file.path(tables.dir, "FigS5ab_S6_S7_Wilcoxon_P.csv"), row.names = FALSE)
