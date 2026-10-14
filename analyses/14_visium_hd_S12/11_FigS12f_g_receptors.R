#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S12f,g - collagen receptors in YAP-Rev high tumour cells (from 10_receptors_4blocks.R)
#  S12f: 074_D1, receptors detected in >= 2% of cells; bars = ratio YAP-Rev high / other tumour cells at the same position
#        (singlets), 95% CI (spatial block bootstrap); purple = CI excluding 1
#  S12g: ITGA2 and SDC4 in the 4 Visium HD sections (red = relapse, blue = no relapse; transparent = not significant)
# Usage : Rscript analyses/14_visium_hd_S12/11_FigS12f_g_receptors.R
# =====================================================================================================================
suppressPackageStartupMessages(library(ggplot2))
out <- file.path(Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd"), "figures_S12")
R <- read.csv(file.path(out, "F9_receptors_YAPRevhigh_vs_other_4blocks.csv"))
R$sig <- R$lo > 1 | R$hi < 1
th <- theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9))

d <- R[R$block == "074_D1 (R)" & pmax(R$pct_YAPhigh, R$pct_other) >= 2 & R$gene != "CD47", ]
d <- d[order(-d$ratio), ]; d$gene <- factor(d$gene, levels = d$gene)
g <- ggplot(d, aes(gene, ratio, fill = sig)) + geom_col(width = .65) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .2, linewidth = .4) + geom_hline(yintercept = 1, linetype = 2, linewidth = .3) +
  scale_fill_manual(values = c(`TRUE` = "#6A51A3", `FALSE` = "grey75"), labels = c(`TRUE` = "95% CI excluding 1", `FALSE` = "NS"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = NULL, y = "Expression in YAP-Rev-high tumor cells\n(relative to other tumor cells)", title = "Collagen receptors, 074_D1 (singlets)") +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(file.path(out, "FigS12f_collagen_receptors_YAPRevhigh_074_D1.pdf"), g, width = 3.8, height = 3.6, device = cairo_pdf)
n074 <- d$n_YAPhigh[1]; no074 <- d$n_other[1]

e <- R[R$gene %in% c("ITGA2", "SDC4"), ]
e$block <- factor(e$block, levels = c("074_D1 (R)", "395_C1 (R)", "21-77_B1 (NR)", "295_B1 (NR)"))
e$lab <- factor(sub(" \\(.*", "", e$block), levels = c("074_D1", "395_C1", "21-77_B1", "295_B1"))
e$grp <- ifelse(grepl("\\(R\\)", e$block), "Relapse", "No relapse")
g2 <- ggplot(e, aes(lab, ratio, fill = grp, alpha = sig)) + geom_col(width = .65) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .2, linewidth = .4, alpha = 1) + geom_hline(yintercept = 1, linetype = 2, linewidth = .3) +
  facet_wrap(~gene, nrow = 1) + scale_fill_manual(values = c(Relapse = "#FC7171", `No relapse` = "#9197FF"), name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .35), guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = NULL, y = "Expression in YAP-Rev-high tumor cells\n(relative to other tumor cells)", title = "ITGA2 and SDC4 in the four Visium HD samples (singlets)") +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(size = 9, face = "italic"))
ggsave(file.path(out, "FigS12g_ITGA2_SDC4_YAPRevhigh_4samples.pdf"), g2, width = 4.2, height = 3.6, device = cairo_pdf)
nn <- unique(e[, c("lab", "n_YAPhigh")]); nn <- nn[order(nn$lab), ]
print(d[, c("gene", "ratio", "lo", "hi", "expected_contamination")], digits = 3)
