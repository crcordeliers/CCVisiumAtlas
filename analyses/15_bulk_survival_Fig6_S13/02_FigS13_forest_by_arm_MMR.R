#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S13 - PETACC8 forest plots (HR per 1-SD increase of each signature, univariate Cox, time to
# recurrence) stratified by treatment arm (a: FOLFOX circles / FOLFOX + cetuximab squares) and by MMR status
# (b: pMMR circles / dMMR squares). Unadjusted models (strata too small for stable multivariable models).
# Input : <BULK_RDS> (see common_bulk.R)
# Output: <WORKDIR>/Figures/bulk/FigS13_PETACC8_forest_by_arm_and_MMR.pdf (+ tables)
# Usage : CRCVISIUM_BULK=<rds> Rscript analyses/15_bulk_survival_Fig6_S13/02_FigS13_forest_by_arm_MMR.R
# =====================================================================================================================
source("analyses/15_bulk_survival_Fig6_S13/common_bulk.R")
suppressPackageStartupMessages(library(patchwork))
PETACC.tot$mmr_status <- factor(PETACC.tot$mmr_status, levels = c("pMMR", "dMMR"))
stopifnot("name_trt" %in% colnames(PETACC.tot))

# =========================================================
# Helpers
# =========================================================
label_map <- c(meta_names, stromal_names, "YAPRev" = "YAP-Rev")
meta_order_labels    <- unname(meta_names)
stromal_order_labels <- unname(stromal_names)
level_labels <- c(meta_order_labels, stromal_order_labels, "YAP-Rev")

get_var_cat_list <- function(df) {
  meta_vars    <- intersect(names(meta_names), colnames(df))
  stromal_vars <- intersect(names(stromal_names), colnames(df))
  yap_var      <- intersect("YAPRev", colnames(df))
  var_list <- c(meta_vars, stromal_vars, yap_var)
  cat_list <- c(
    rep("Metaclusters",        length(meta_vars)),
    rep("Stromal subclusters", length(stromal_vars)),
    if (length(yap_var)) "YAPRev" else character(0)
  )
  list(var_list = var_list, cat_list = cat_list)
}

# Linear HR axis with limits/breaks fixed per panel; CIs whose upper bound exceeds the last break are truncated
# and flagged with ">" (see build_forest_panel).

p_label_col <- function(res) {
  res %>% dplyr::mutate(
    p_label = dplyr::case_when(
      p.value < 1e-4 ~ 'italic(P) < "0.0001"',
      p.value < 1e-3 ~ 'italic(P) < "0.001"',
      TRUE           ~ paste0("italic(P) == ", formatC(p.value, format = "f", digits = 3))
    )
  )
}

base_theme <- function() {
  list(
    ggpubr::theme_pubclean(),
    theme(
      axis.title.y    = element_blank(),
      axis.text.y     = ggtext::element_markdown(size = 13, face = "bold"),
      panel.spacing.y = unit(1.1, "lines"),
      strip.text.y    = element_text(face = "bold"),
      legend.position = "top",
      legend.title    = element_text(face = "bold"),
      plot.margin     = margin(6, 25, 6, 6),
      plot.caption    = element_text(size = 8, color = "grey30", hjust = 0)
    )
  )
}

# Generic forest-plot panel with explicit breaks (y_breaks).
build_forest_panel <- function(res, group_col, dodge, shape_values, shape_labels,
                                legend_name, title, hide_y = FALSE, y_breaks) {
  max_tick <- max(y_breaks)
  min_tick <- min(y_breaks)

  res <- res %>%
    mutate(
      truncated   = conf.high > max_tick,
      ymax_plot   = pmin(conf.high, max_tick),
      arrow_label = ifelse(truncated, ">", "")
    )
  has_trunc <- any(res$truncated)

  # margins proportional to the displayed range (room for the truncation arrow and the P values)
  tick_range  <- max_tick - min_tick
  margin_unit <- tick_range * 0.15
  arrow_y     <- max_tick + 0.3 * margin_unit
  pval_y      <- max_tick + 1.3 * margin_unit
  upper_limit <- max_tick + 1.8 * margin_unit
  lower_limit <- min_tick - 0.3 * margin_unit

  p <- ggplot(res, aes(x = label)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
    geom_pointrange(
      aes(y = estimate, ymin = conf.low, ymax = ymax_plot,
          group = .data[[group_col]], shape = .data[[group_col]], fill = color_hex),
      position = dodge, color = "black", linewidth = 0.25, stroke = 0.25, size = 0.9
    ) +
    (if (has_trunc) geom_text(
      aes(x = label, y = arrow_y, label = arrow_label, group = .data[[group_col]]),
      position = dodge, hjust = 0, size = 3.6, fontface = "bold", color = "black"
    ) else NULL) +
    geom_text(
      aes(x = label, y = pval_y, label = p_label, group = .data[[group_col]]),
      position = dodge, hjust = 1, size = 3.6, color = "black", parse = TRUE
    ) +
    coord_flip(clip = "off") +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    scale_fill_identity() +
    scale_shape_manual(values = shape_values, name = legend_name, labels = shape_labels) +
    scale_y_continuous(
      limits = c(lower_limit, upper_limit),
      breaks = y_breaks,
      labels = scales::label_number(accuracy = 0.01, drop0trailing = TRUE),
      expand = expansion(mult = c(0, 0))
    ) +
    base_theme() +
    labs(
      y = "Hazard Ratio (TTR), unadjusted", x = NULL, title = title,
      caption = if (has_trunc) paste0("> : 95% CI truncated at HR = ", round(max_tick, 2)) else NULL
    )

  if (hide_y) p <- p + theme(axis.text.y = element_blank())
  p
}

# =========================================================
# S13a - stratified by treatment arm (FOLFOX vs FOLFOX + cetuximab), unadjusted
# =========================================================
arm_levels <- c("Folfox", "Folfox+Cetux")

build_cox_by_arm <- function(df) {
  vc <- get_var_cat_list(df)
  var_list <- vc$var_list; cat_list <- vc$cat_list
  time_var <- "recurr"; status_var <- "recurr_cens"

  res <- lapply(seq_along(var_list), function(i) {
    v <- var_list[i]
    per_arm <- lapply(arm_levels, function(arm) {
      df_arm <- df[df$name_trt == arm, , drop = FALSE]
      keep <- stats::complete.cases(df_arm[, c(time_var, status_var, v), drop = FALSE])
      df_a <- df_arm[keep, , drop = FALSE]
      fm   <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, status_var, v))
      fit  <- tryCatch(coxph(fm, data = df_a), error = function(e) NULL)
      out  <- if (!is.null(fit)) broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) else NULL
      if (!is.null(out)) out <- dplyr::filter(out, .data$term == v)
      if (!is.null(out)) out <- dplyr::mutate(out, arm = arm, n_arm = nrow(df_a))
      out
    })
    bind_rows(per_arm) |>
      mutate(var = v, category = cat_list[i],
             label = ifelse(v %in% names(label_map), unname(label_map[[v]]), v),
             color_hex = unname(colors_all[v]))
  })
  bind_rows(res)
}

res_arm <- build_cox_by_arm(PETACC.tot) |>
  filter(is.finite(estimate), is.finite(conf.low), is.finite(conf.high)) |>
  mutate(
    label    = factor(label, levels = rev(level_labels)),
    category = factor(category, levels = c("Metaclusters", "Stromal subclusters", "YAPRev")),
    arm      = factor(arm, levels = arm_levels)
  ) |>
  p_label_col()

n_folfox <- sum(PETACC.tot$name_trt == "Folfox", na.rm = TRUE)
n_cetux  <- sum(PETACC.tot$name_trt == "Folfox+Cetux", na.rm = TRUE)

dodge_arm <- position_dodge2(width = 0.65, preserve = "single", reverse = TRUE)

p_arm <- build_forest_panel(
  res_arm, group_col = "arm", dodge = dodge_arm,
  shape_values = c("Folfox" = 21, "Folfox+Cetux" = 22),
  shape_labels = c("Folfox" = paste0("Folfox (n=", n_folfox, ")"),
                    "Folfox+Cetux" = paste0("Folfox+Cetux (n=", n_cetux, ")")),
  legend_name = "Treatment arm",
  title = "Stratified by treatment arm",
  y_breaks = c(0.75, 1, 1.25, 1.5, 1.75)
)

# =========================================================
# S13b - stratified by MMR status, unadjusted
# =========================================================
mmr_levels <- c("pMMR", "dMMR")

build_cox_by_mmr <- function(df) {
  vc <- get_var_cat_list(df)
  var_list <- vc$var_list; cat_list <- vc$cat_list
  time_var <- "recurr"; status_var <- "recurr_cens"

  res <- lapply(seq_along(var_list), function(i) {
    v <- var_list[i]
    per_mmr <- lapply(mmr_levels, function(m) {
      df_m <- df[!is.na(df$mmr_status) & df$mmr_status == m, , drop = FALSE]
      keep <- stats::complete.cases(df_m[, c(time_var, status_var, v), drop = FALSE])
      df_k <- df_m[keep, , drop = FALSE]
      fm   <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, status_var, v))
      fit  <- tryCatch(coxph(fm, data = df_k), error = function(e) NULL)
      out  <- if (!is.null(fit)) broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) else NULL
      if (!is.null(out)) out <- dplyr::filter(out, .data$term == v)
      if (!is.null(out)) out <- dplyr::mutate(out, mmr = m, n_mmr = nrow(df_k))
      out
    })
    bind_rows(per_mmr) |>
      mutate(var = v, category = cat_list[i],
             label = ifelse(v %in% names(label_map), unname(label_map[[v]]), v),
             color_hex = unname(colors_all[v]))
  })
  bind_rows(res)
}

res_mmr <- build_cox_by_mmr(PETACC.tot) |>
  filter(is.finite(estimate), is.finite(conf.low), is.finite(conf.high)) |>
  mutate(
    label    = factor(label, levels = rev(level_labels)),
    category = factor(category, levels = c("Metaclusters", "Stromal subclusters", "YAPRev")),
    mmr      = factor(mmr, levels = mmr_levels)
  ) |>
  p_label_col()

n_pmmr <- sum(PETACC.tot$mmr_status == "pMMR", na.rm = TRUE)
n_dmmr <- sum(PETACC.tot$mmr_status == "dMMR", na.rm = TRUE)

dodge_mmr <- position_dodge2(width = 0.65, preserve = "single", reverse = TRUE)

p_mmr <- build_forest_panel(
  res_mmr, group_col = "mmr", dodge = dodge_mmr,
  shape_values = c("pMMR" = 21, "dMMR" = 22),
  shape_labels = c("pMMR" = paste0("pMMR (n=", n_pmmr, ")"),
                    "dMMR" = paste0("dMMR (n=", n_dmmr, ")")),
  legend_name = "MMR status",
  title = "Stratified by MMR status",
  hide_y = TRUE,   # right panel: signature names not repeated
  y_breaks = c(0.5, 1, 1.5, 2, 2.5)
)

# =========================================================
# both panels side by side
# =========================================================
p_side_by_side <- p_arm + p_mmr + patchwork::plot_layout(ncol = 2, widths = c(1.35, 1))

pdf(file.path(out_dir, "FigS13_PETACC8_forest_by_arm_and_MMR.pdf"), width = 14, height = 8)
print(p_side_by_side)
dev.off()
write.csv(res_arm, file.path(out_dir, "FigS13a_PETACC8_Cox_by_treatment_arm.csv"), row.names = FALSE)
write.csv(res_mmr, file.path(out_dir, "FigS13b_PETACC8_Cox_by_MMR.csv"), row.names = FALSE)
