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
set.seed(1)   # fixed jitter of the points (no effect on the statistics)

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
         group = factor(ifelse(recurr_cens == 1, "Relapse", "No relapse"), levels = c("No relapse", "Relapse")))
write.csv(pp, file.path(tables.dir, "FigS8h_YAPRev_matCAF_cooccupation.csv"), row.names = FALSE)

pval <- wilcox.test(pct_of_all_spots ~ group, data = pp)$p.value
p_lab <- if (pval < 1e-4) "italic(P) < '0.0001'" else if (pval < 1e-3) "italic(P) < '0.001'" else
  paste0("italic(P) == '", formatC(pval, format = "f", digits = 3), "'")
fill <- c("No relapse" = unname(relapse_colors["0"]), "Relapse" = unname(relapse_colors["1"]))

p <- ggplot(pp, aes(x = "matCAF\n(cS26+27+28)", y = pct_of_all_spots, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 0.6) +
  geom_point(aes(color = group), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 1.5, alpha = 0.6) +
  annotate("text", x = 1, y = max(pp$pct_of_all_spots) * 1.12, label = p_lab, parse = TRUE, size = 7, fontface = "bold") +
  scale_fill_manual(values = fill, name = "") + scale_color_manual(values = fill, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(title = "matCAF (cS26+cS27+cS28)", x = "", y = "% spots co-occupied") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size = rel(1.6), color = "black"),
        axis.title.y = element_text(face = "bold", size = rel(1.8)), plot.title = element_text(face = "bold", size = rel(1.8)),
        legend.text = element_text(face = "bold", size = rel(1.4)), legend.position = "top")
ggsave(file.path(figures.dir, "FigS8h_YAPRev_matCAF_cooccupation.pdf"), p, width = 5, height = 5)
message("Mann-Whitney P = ", signif(pval, 3))
