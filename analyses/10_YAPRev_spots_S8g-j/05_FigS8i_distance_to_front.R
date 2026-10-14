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
  mutate(group = factor(ifelse(recurr == 1, "Relapse", "No Relapse"), levels = c("No Relapse", "Relapse")))
write.csv(pp, file.path(tables.dir, "FigS8i_delta_distance_to_front.csv"), row.names = FALSE)

wt <- wilcox.test(delta_highminuslow ~ group, data = pp)
source("R/boxplot_relapse.R")
p <- boxplot_relapse(pp, "delta_highminuslow", wt$p.value, "Invasive front", "\u0394 distance to front (\u00b5m)",
                     sprintf("Patients (n = %d)", nrow(pp)), size = "n_high", size_breaks = c(25, 250, 1500), hline = 0)
ggsave(file.path(figures.dir, "FigS8i_delta_distance_to_front.pdf"), p, width = 5, height = 5, device = cairo_pdf)
message("n = ", sum(pp$group == "No Relapse"), " / ", sum(pp$group == "Relapse"), "; Mann-Whitney P = ", signif(wt$p.value, 3),
        "; median delta (um) No Relapse / Relapse: ", paste(round(tapply(pp$delta_highminuslow, pp$group, median)), collapse = " / "))
