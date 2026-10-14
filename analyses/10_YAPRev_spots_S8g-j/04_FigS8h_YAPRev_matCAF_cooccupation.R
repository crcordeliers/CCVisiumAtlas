#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S8h - patient-level % of spots that are both YAP-Rev high and matCAF high, Relapse vs No relapse
#
# YAP-Rev high : YAPRev_raw >= 0.40 (01/02_*.R).
# matCAF       : sum of cell2location proportions cS26 (Myofibro) + cS27 (CXCL14+ CAF) + cS28 (GREM1+ CAF);
#                matCAF high = top 5% of the pooled distribution (threshold 0.4229, i.e. >= 42% of the spot).
# Co-occupation: the same spot is YAP-Rev high and matCAF high; expressed as % of all spots of the patient.
# Groups (34 No relapse / 14 Relapse) compared with a two-sided Mann-Whitney test.
# Inputs : <YAPREV_DIR>/merged_spots/*.csv ; <C2L_DIR>/tables/c2l_98cat_proportions_<sample>.csv ; data/data_clin_48_filtered.RDS
# Output : Figures/FigS8h_YAPRev_matCAF_cooccupation.pdf, Tables/FigS8h_YAPRev_matCAF_cooccupation.csv
# Usage  : Rscript analyses/10_YAPRev_spots_S8g-j/04_FigS8h_YAPRev_matCAF_cooccupation.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

yaprev_dir <- Sys.getenv("CRCVISIUM_YAPREV", "~/data/yaprev_spots")
c2l_tables <- file.path(Sys.getenv("CRCVISIUM_C2L", "~/data/cell2location"), "tables")
clin <- readRDS(file.path(data.dir, "data_clin_48_filtered.RDS"))
YAP_THR <- 0.40
MATCAF  <- c("cS26 (Myofibro)", "cS27 (CXCL14+ CAF)", "cS28 (GREM1+ CAF)")

d <- bind_rows(lapply(Sample.names, function(s) {
  m <- read.csv(file.path(yaprev_dir, "merged_spots", paste0(s, "_merged.csv")))
  p <- read.csv(file.path(c2l_tables, paste0("c2l_98cat_proportions_", s, ".csv")), check.names = FALSE)
  colnames(p)[1] <- "barcode"
  inner_join(data.frame(barcode = m$barcode, sample = s, YAPRev_raw = m$YAPRev_raw),
             data.frame(barcode = p$barcode, matCAF = rowSums(p[, MATCAF])), by = "barcode")
}))
matcaf_thr <- quantile(d$matCAF, 0.95, names = FALSE)
message("matCAF top-5% threshold: ", round(matcaf_thr, 4))

pp <- d %>% group_by(sample) %>%
  summarise(n_spots = n(), n_co_occupied = sum(YAPRev_raw >= YAP_THR & matCAF >= matcaf_thr), .groups = "drop") %>%
  mutate(pct_of_all_spots = 100 * n_co_occupied / n_spots,
         recurr_cens = clin$recurr_cens[match(sample, clin$sample)],
         group = factor(ifelse(recurr_cens == 1, "Relapse", "No Relapse"), levels = c("No Relapse", "Relapse")))
write.csv(pp, file.path(tables.dir, "FigS8h_YAPRev_matCAF_cooccupation.csv"), row.names = FALSE)

pval <- wilcox.test(pct_of_all_spots ~ group, data = pp)$p.value
source("R/boxplot_relapse.R")
p <- boxplot_relapse(pp, "pct_of_all_spots", pval, "YAP-Rev high & matCAF high", "% of spots per patient",
                     sprintf("Patients (n = %d)", nrow(pp)))
ggsave(file.path(figures.dir, "FigS8h_YAPRev_matCAF_cooccupation.pdf"), p, width = 5, height = 5)
message("Mann-Whitney P = ", signif(pval, 3))
