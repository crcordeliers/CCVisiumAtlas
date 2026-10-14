#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S8i (left) - YAP-Rev high spots are closer to the invasive front in relapsing tumours
#
# For each of the 31 patients with an annotated front: delta = median(distance to front of YAP-Rev high spots) -
# median(distance of YAP-Rev low spots) (YAPRev_raw >= 0.40). Negative delta = YAP-Rev high spots enriched near the front.
# Point size = number of YAP-Rev high spots. No relapse (n = 23) vs Relapse (n = 8): two-sided Mann-Whitney test.
# Inputs : <YAPREV_DIR>/merged_spots/*.csv, <YAPREV_DIR>/spot_distances/*.csv (03_distance_to_invasive_front.py)
# Output : Figures/FigS8i_delta_distance_to_front.pdf, Tables/FigS8i_delta_distance_to_front.csv
# Usage  : Rscript analyses/10_YAPRev_spots_S8g-j/05_FigS8i_distance_to_front.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
set.seed(1)   # fixed jitter of the points (no effect on the statistics)

yaprev_dir <- Sys.getenv("CRCVISIUM_YAPREV", "~/data/yaprev_spots")
dist_dir   <- Sys.getenv("CRCVISIUM_FRONT_OUT", file.path(yaprev_dir, "spot_distances"))
clin <- readRDS(file.path(data.dir, "data_clin_48_filtered.RDS"))
THR <- 0.40

front_samples <- sub("\\.csv$", "", list.files(dist_dir, pattern = "\\.csv$"))
front_df <- bind_rows(lapply(front_samples, function(s) {
  m  <- read.csv(file.path(yaprev_dir, "merged_spots", paste0(s, "_merged.csv")))
  dd <- read.csv(file.path(dist_dir, paste0(s, ".csv")), header = FALSE, col.names = c("sample", "YAPRev", "distance_um"))
  stopifnot(nrow(m) == nrow(dd))
  data.frame(sample = s, YAPRev_raw = m$YAPRev_raw, distance_um = dd$distance_um)
}))
front_df$recurr <- clin$recurr_cens[match(front_df$sample, clin$sample)]

pp <- front_df %>% group_by(sample, recurr) %>%
  summarise(n_high = sum(YAPRev_raw >= THR),
            median_dist_high = if (n_high > 0) median(distance_um[YAPRev_raw >= THR]) else NA_real_,
            median_dist_low  = median(distance_um[YAPRev_raw < THR]), .groups = "drop") %>%
  mutate(delta_highminuslow = median_dist_high - median_dist_low) %>%
  filter(!is.na(delta_highminuslow)) %>%
  mutate(group = factor(ifelse(recurr == 1, "Relapse", "No relapse"), levels = c("No relapse", "Relapse")))
write.csv(pp, file.path(tables.dir, "FigS8i_delta_distance_to_front.csv"), row.names = FALSE)

wt <- wilcox.test(delta_highminuslow ~ group, data = pp)
p_lab <- if (wt$p.value < 1e-3) "italic(P) < '0.001'" else paste0("italic(P) == '", formatC(wt$p.value, format = "f", digits = 3), "'")
fill <- c("No relapse" = unname(relapse_colors["0"]), "Relapse" = unname(relapse_colors["1"]))
y_range <- diff(range(pp$delta_highminuslow))

p <- ggplot(pp, aes(x = group, y = delta_highminuslow)) +
  geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(size = n_high, color = group), width = 0.12, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  annotate("text", x = 1.5, y = max(pp$delta_highminuslow) + 0.08 * y_range, label = p_lab, parse = TRUE, size = 7, fontface = "bold") +
  scale_color_manual(values = fill, guide = "none") +
  scale_size_continuous(name = "n high spots", breaks = c(25, 250, 1500)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = paste0("YAP-Rev >= ", THR), x = "", y = "Delta distance to front, um") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size = rel(1.6), color = "black"),
        axis.title.y = element_text(face = "bold", size = rel(1.8)), plot.title = element_text(face = "bold", size = rel(1.8)),
        legend.title = element_text(face = "bold", size = rel(1.2)), legend.text = element_text(face = "bold", size = rel(1.2)),
        legend.position = "top")
ggsave(file.path(figures.dir, "FigS8i_delta_distance_to_front.pdf"), p, width = 5, height = 5)
message("n = ", sum(pp$group == "No relapse"), " / ", sum(pp$group == "Relapse"), "; Mann-Whitney P = ", signif(wt$p.value, 3),
        "; median delta (um) NR / R: ", paste(round(tapply(pp$delta_highminuslow, pp$group, median)), collapse = " / "))
