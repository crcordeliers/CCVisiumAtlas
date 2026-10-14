# Patient-level boxplot, No Relapse vs Relapse (Supplementary Fig. S8g-i; same style as S8c, S8f)
# df: one row per patient with columns `group` ("No Relapse" / "Relapse") and the variable `y`; pval: P of the test.
boxplot_relapse <- function(df, y, pval, title, ylab, xlab, size = NULL, size_breaks = waiver(), hline = NULL) {
  df$group <- factor(df$group, levels = c("No Relapse", "Relapse"))
  p_lab <- if (pval < 1e-4) "italic(P) < '0.0001'" else if (pval < 1e-3) "italic(P) < '0.001'" else
    paste0("italic(P) == '", formatC(pval, format = "f", digits = 3), "'")
  rng <- range(df[[y]], na.rm = TRUE)
  pos <- ggplot2::position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8, seed = 1)
  pts <- if (is.null(size)) {
    ggplot2::geom_point(position = pos, size = 1.4, alpha = 0.5, colour = "grey20", show.legend = FALSE)
  } else {
    ggplot2::geom_point(ggplot2::aes(size = .data[[size]]), position = pos, alpha = 0.5, colour = "grey20")
  }
  g <- ggplot2::ggplot(df, ggplot2::aes(x = xlab, y = .data[[y]], fill = group))
  if (!is.null(hline)) g <- g + ggplot2::geom_hline(yintercept = hline, linetype = "dashed", colour = "grey60")
  g <- g +
    ggplot2::geom_boxplot(outlier.shape = NA, position = ggplot2::position_dodge(width = 0.8), alpha = 0.6) +
    pts +
    ggplot2::annotate("text", x = 1, y = rng[2] + 0.08 * diff(rng), label = p_lab, parse = TRUE, size = 7, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c("No Relapse" = "#9197FF", "Relapse" = "#FC7171"), name = "") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.06))) +
    ggplot2::labs(title = title, x = "", y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.6), color = "black"),
                   axis.text.y = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.6), color = "black"),
                   axis.title.y = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.8)),
                   plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.8)),
                   legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.4)),
                   legend.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.4)),
                   legend.position = "top")
  if (!is.null(size)) g <- g + ggplot2::scale_size_continuous(name = "n high spots", breaks = size_breaks, range = c(1, 5)) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1), size = ggplot2::guide_legend(order = 2)) +
    ggplot2::theme(legend.box = "vertical", legend.spacing.y = ggplot2::unit(0, "pt"))
  g
}
