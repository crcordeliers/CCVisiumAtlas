#!/usr/bin/env Rscript
# =====================================================================================================================
# Fig. 6 - prognostic value of the spatially derived signatures in three bulk RNA-seq cohorts (time to recurrence)
#   6a-c  forest plots, HR per 1-SD increase of each signature (ssGSEA, z-scored within cohort), univariate (circles)
#         and multivariable (squares) Cox models; PETACC8 (n = 1,733), IDEA-France (n = 1,248), PRODIGE13 (n = 529)
#   text  multivariable Cox PETACC8: Stroma-A2, YAP-Rev, both, and Stroma-A2 x YAP-Rev interaction
#   6d    predicted HR along Stroma-A2 / YAP-Rev with the other program fixed at Q1 or Q3 (PETACC8)
#   6e-g  Kaplan-Meier curves, 4 groups Low/High Stroma-A2 x Low/High YAP-Rev, PETACC8-derived medians of raw ssGSEA
#         scores applied unchanged to IDEA-France and PRODIGE13; log-rank test and unadjusted Cox vs Low/Low
# Adjustment: grade, WHO performance status, baseline high-risk features (pT4 and/or pN2; pT4 only in PRODIGE13),
# MMR status, bowel complication at diagnosis (occlusion or perforation).
# Signatures: top-20 marker genes of each meta-cluster and stromal sub-cluster (data/reference_TableS2_signatures.csv,
# = Supplementary Table S2) and YAP-Rev (AHNAK, ANXA1, CAPN2). ssGSEA (GSVA) on the union of the signature genes only
# (the ranking is restricted to the signature genes, in the order of Table S2 - the order breaks ties).
# Input : <BULK_RDS> (00_build_bulk_dataset.R): list(PETACC, IDEA, PRODIGE) of expression + clinical tables
# Output: <WORKDIR>/Figures/bulk/
# Usage : CRCVISIUM_BULK=bulk_cohorts.rds Rscript analyses/15_bulk_survival_Fig6_S13/01_Fig6_survival.R
# Env   : R >= 4.3, GSVA, survival, survminer, broom, ggpubr, ggtext
# =====================================================================================================================
source("analyses/15_bulk_survival_Fig6_S13/common_bulk.R")

# =====================================================================================================================
# Fig. 6a-c - forest plots
# =====================================================================================================================
adjustment_vars_for <- function(df, dataset_name) {
  vars <- c("recgrade", "recOMS", "recRisk", "mmr_status", "bowel_complication")
  vars <- vars[vars %in% colnames(df)]
  vars <- vars[vapply(vars, function(v) any(!is.na(df[[v]])), logical(1))]
  unique(vars)
}

cox_formula <- function(time, status, pred, adjust) {
  rhs <- c(pred, adjust)
  as.formula(paste0("Surv(", time, ", ", status, ") ~ ", paste(rhs, collapse = " + ")))
}

# =========================================================
# Forest plots Unadjusted vs Adjusted (Metaclusters / Stromal / YAP-Rev) -------
# =========================================================

label_map <- c(meta_names, stromal_names, "YAPRev" = "YAP-Rev")

build_cox_unadj_adj <- function(cohort) {
  df <- get(paste0(cohort, ".tot"))
  cohort_label <- dplyr::recode(
    cohort,
    PETACC = "PETACC8",
    IDEA   = "IDEA-France",
    PRODIGE = "PRODIGE13"
  )
  
  meta_vars    <- intersect(names(meta_names), colnames(df))
  stromal_vars <- intersect(names(stromal_names), colnames(df))
  yap_var      <- intersect("YAPRev", colnames(df))
  
  var_list <- c(meta_vars, stromal_vars, yap_var)
  if (!length(var_list)) return(tibble())
  
  cat_list <- c(
    rep("Metaclusters",        length(meta_vars)),
    rep("Stromal subclusters", length(stromal_vars)),
    if (length(yap_var)) "YAPRev" else character(0)
  )
  
  time_var   <- if (cohort == "PETACC") "recurr"    else "TTR_time"
  status_var <- if (cohort == "PETACC") "recurr_cens" else "etat_TTR"
  adjust_vars <- adjustment_vars_for(df, cohort)
  
  res <- lapply(seq_along(var_list), function(i) {
    v <- var_list[i]
    
    # Unadjusted
    keep_u <- stats::complete.cases(df[, c(time_var, status_var, v), drop = FALSE])
    df_u   <- df[keep_u, , drop = FALSE]
    fm_u   <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, status_var, v))
    fit_u  <- tryCatch(coxph(fm_u, data = df_u), error = function(e) NULL)
    out_u  <- if (!is.null(fit_u)) broom::tidy(fit_u, exponentiate = TRUE, conf.int = TRUE) else NULL
    if (!is.null(out_u)) out_u <- dplyr::filter(out_u, .data$term == v)
    
    # Adjusted
    rhs_a  <- c(v, adjust_vars)
    keep_a <- stats::complete.cases(df[, c(time_var, status_var, rhs_a), drop = FALSE])
    df_a   <- df[keep_a, , drop = FALSE]
    fm_a   <- cox_formula(time = time_var, status = status_var, pred = v, adjust = adjust_vars)
    fit_a  <- tryCatch(coxph(fm_a, data = df_a), error = function(e) NULL)
    out_a  <- if (!is.null(fit_a)) broom::tidy(fit_a, exponentiate = TRUE, conf.int = TRUE) else NULL
    if (!is.null(out_a)) out_a <- dplyr::filter(out_a, .data$term == v)
    
    bind_rows(
      if (!is.null(out_u)) mutate(out_u, model = "Unadjusted") else NULL,
      if (!is.null(out_a)) mutate(out_a, model = "Adjusted")   else NULL
    ) |>
      mutate(
        var       = v,
        cohort    = cohort_label,
        category  = cat_list[i],
        label     = ifelse(v %in% names(label_map), unname(label_map[[v]]), v),
        color_hex = unname(colors_all[v])
      )
  })
  
  bind_rows(res)
}

res_both <- bind_rows(
  build_cox_unadj_adj("PETACC"),
  build_cox_unadj_adj("IDEA"),
  build_cox_unadj_adj("PRODIGE")
)

# remove rows with infinite or missing estimates / CIs
res_both <- res_both |>
  filter(
    is.finite(estimate),
    is.finite(conf.low),
    is.finite(conf.high)
  )


meta_order_labels    <- unname(meta_names)
stromal_order_labels <- unname(stromal_names)
level_labels <- c(meta_order_labels, stromal_order_labels, "YAP-Rev")

res_both <- res_both |>
  mutate(
    label    = factor(.data$label, levels = rev(level_labels)),
    category = factor(.data$category, levels = c("Metaclusters", "Stromal subclusters", "YAPRev")),
    model    = factor(.data$model, levels = c("Unadjusted", "Adjusted"))
  )

plot_one_cohort <- function(cohort_label) {
  dfc <- dplyr::filter(res_both, .data$cohort == cohort_label)
  if (!nrow(dfc)) return(invisible(NULL))
  
  n_cohort <- switch(
    cohort_label,
    "PETACC8"     = nrow(PETACC.tot),
    "IDEA-France" = nrow(IDEA.tot),
    "PRODIGE13"   = nrow(PRODIGE.tot)
  )
  
  dodge <- position_dodge2(width = 0.65, preserve = "single", reverse = TRUE)
  
  dfc <- dfc %>%
    dplyr::mutate(
      p_label = dplyr::case_when(
        p.value < 1e-4 ~ 'italic(P) < "0.0001"',
        p.value < 1e-3 ~ 'italic(P) < "0.001"',
        TRUE           ~ paste0("italic(P) == ", formatC(p.value, format = "f", digits = 3))
      )
    )
  
  ## ---- HR axis: ticks within the CI range, space on the right for the P values ----
  hr_min <- min(dfc$conf.low,  na.rm = TRUE)
  hr_max <- max(dfc$conf.high, na.rm = TRUE)
  span   <- hr_max - hr_min
  
  # ~8 ticks with a round step
  raw_step <- span / 8
  step <- dplyr::case_when(
    raw_step <= 0.10 ~ 0.10,
    raw_step <= 0.25 ~ 0.25,
    TRUE             ~ 0.50
  )
  
  # first tick just above the lowest CI bound
  min_tick <- ceiling(hr_min / step) * step
  
  # last tick just below the highest CI bound
  max_tick <- floor(hr_max / step) * step
  
  # at least one tick < 1 and one > 1
  if (max_tick <= 1) {
    max_tick <- 1 + step
  }
  if (min_tick >= 1) {
    min_tick <- 1 - step
  }
  
  # ticks only between min_tick and max_tick
  breaks_hr <- seq(min_tick, max_tick, by = step)
  
  # axis limits: slightly below min_tick, and above max_tick to leave room for the P values
  lower_limit <- hr_min - 0.5 * step
  upper_limit <- max_tick + 1.7 * step
  
  # P values placed in the area without ticks on the right
  pval_y <- max_tick + 1.4 * step
  ## --------------------------------------------------------------------
  
  
  ggplot(dfc, aes(x = label)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
    
    geom_pointrange(
      aes(
        y = estimate, ymin = conf.low, ymax = conf.high,
        group = model, shape = model, fill = color_hex
      ),
      position = dodge,
      color = "black",
      linewidth = 0.25,
      stroke = 0.25,
      size = 1
    ) +
    
    # P values right-aligned
    geom_text(
      aes(
        x = label,
        y = pval_y,
        label = p_label,
        group = model
      ),
      position = dodge,
      hjust = 1,
      size  = 4.5,
      color = "black",
      parse = TRUE
    ) +
    
    coord_flip(clip = "off") +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    scale_fill_identity() +
    scale_shape_manual(values = c("Unadjusted" = 21, "Adjusted" = 22),
                       name = "Model") +
    scale_y_continuous(
      limits = c(lower_limit, upper_limit),
      breaks = breaks_hr,
      expand = expansion(mult = c(0, 0))
    ) +
    ggpubr::theme_pubclean() +
    theme(
      axis.title.y    = element_blank(),
      axis.text.y     = ggtext::element_markdown(size = 15, face = "bold"),
      panel.spacing.y = unit(1.1, "lines"),
      strip.text.y    = element_text(face = "bold"),
      legend.position = "top",
      legend.title    = element_text(face = "bold"),
      plot.margin     = margin(6, 25, 6, 6)
    ) +
    labs(
      y = "Hazard Ratio (TTR)", x = NULL
    )
}

p_petacc  <- plot_one_cohort("PETACC8")
p_idea    <- plot_one_cohort("IDEA-France")
p_prodige <- plot_one_cohort("PRODIGE13")

if (!is.null(p_petacc)) {
  pdf(file.path(out_dir, "Fig6a_forest_PETACC8.pdf"), width = 8, height = 8)
  print(p_petacc)
  dev.off()
}
if (!is.null(p_idea)) {
  pdf(file.path(out_dir, "Fig6b_forest_IDEA.pdf"), width = 8, height = 8)
  print(p_idea)
  dev.off()
}
if (!is.null(p_prodige)) {
  pdf(file.path(out_dir, "Fig6c_forest_PRODIGE13.pdf"), width = 8, height = 8)
  print(p_prodige)
  dev.off()
}

write.csv(res_both, file.path(out_dir, "Fig6abc_forest_Cox.csv"), row.names = FALSE)


# =====================================================================================================================
# Multivariable Cox PETACC8 with Stroma-A2 x YAP-Rev interaction (numbers given in the Results)
# =====================================================================================================================
safe_adjustments <- function(df) {
  vars <- c("recgrade", "recOMS", "recRisk", "mmr_status", "bowel_complication")
  vars[vars %in% names(df) & vapply(vars, function(v) any(!is.na(df[[v]])), logical(1))]
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) "p < 0.0001" else if (p < 1e-3) "p < 0.001" else sprintf("p = %.3f", p)
}

summ_line <- function(tidy_row) {
  sprintf(
    "HR = %.2f, 95%% CI %.2f–%.2f, %s",
    tidy_row$estimate, tidy_row$conf.low, tidy_row$conf.high, fmt_p(tidy_row$p.value)
  )
}

required_petacc <- c("recurr", "recurr_cens", "StromalSubcluster_2", "YAPRev")
stopifnot(all(required_petacc %in% names(PETACC.tot)))

required_idea_prodige <- c("TTR_time", "etat_TTR", "StromalSubcluster_2", "YAPRev")
stopifnot(all(required_idea_prodige %in% names(IDEA.tot)))
stopifnot(all(required_idea_prodige %in% names(PRODIGE.tot)))

petacc   <- PETACC.tot
adj_vars <- safe_adjustments(petacc)
pred_vars <- c("StromalSubcluster_2", "YAPRev")

cc_train <- complete.cases(petacc[, c("recurr", "recurr_cens", pred_vars, adj_vars), drop = FALSE])
train <- petacc[cc_train, , drop = FALSE]

means <- sapply(train[, pred_vars, drop = FALSE], mean, na.rm = TRUE)
sds   <- sapply(train[, pred_vars, drop = FALSE],  sd,   na.rm = TRUE)

train_std <- train
for (v in pred_vars) train_std[[v]] <- (train[[v]] - means[[v]]) / sds[[v]]

fml_train <- as.formula(
  paste0(
    "Surv(time = recurr, event = recurr_cens) ~ ",
    "StromalSubcluster_2 * YAPRev",
    if (length(adj_vars)) paste0(" + ", paste(adj_vars, collapse = " + ")) else ""
  )
)

fit <- coxph(fml_train, data = train_std, x = TRUE, y = TRUE)
tid <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

row_s2  <- tid %>% filter(term == "StromalSubcluster_2")         %>% slice(1)
row_yap <- tid %>% filter(term == "YAPRev")                      %>% slice(1)
row_int <- tid %>% filter(term == "StromalSubcluster_2:YAPRev")  %>% slice(1)


cat("\n--- Multivariable Cox on PETACC8 (standardised predictors; text of the Results) ---\n")
cat(sprintf("Stroma-A2 (per +1 SD): %s\n",  summ_line(row_s2)))
cat(sprintf("YAP-Rev  (per +1 SD): %s\n",  summ_line(row_yap)))
cat(sprintf("Interaction (Stroma-A2 x YAP-Rev): %s\n", summ_line(row_int)))
# model with both signatures, without interaction
fit_both <- coxph(update(fml_train, . ~ . - StromalSubcluster_2:YAPRev), data = train_std)
print(broom::tidy(fit_both, exponentiate = TRUE, conf.int = TRUE)[1:2, ])

# =====================================================================================================================
# Fig. 6d - interaction curves (PETACC8)
# =====================================================================================================================
# Adjustment of the Fig. 6d model: grade, WHO performance status, baseline high-risk features, MMR status
PETACC_6d <- PETACC.tot
PETACC_6d$bowel_complication <- NA
# -----------------------------
# 1) PETACC8 interaction model
# -----------------------------
fit_interaction_petacc <- function(df) {
  time_var   <- "recurr"
  status_var <- "recurr_cens"
  
  adj <- safe_adjustments(df)
  needed <- c(time_var, status_var, "StromalSubcluster_2", "YAPRev", adj)
  
  cc  <- stats::complete.cases(df[, needed, drop = FALSE])
  dfc <- df[cc, , drop = FALSE]
  
  # z-scored scores used as such (no re-scaling)
  rhs <- paste(c("StromalSubcluster_2 * YAPRev", adj), collapse = " + ")
  fml <- as.formula(
    paste0("Surv(", time_var, ", ", status_var, ") ~ ", rhs)
  )
  
  fit <- coxph(fml, data = dfc, x = TRUE, y = TRUE)
  
  list(
    fit      = fit,
    data_z   = dfc,
    time_var = time_var,
    status_var = status_var
  )
}

fit_PETACC <- fit_interaction_petacc(PETACC_6d)

# -----------------------------
# 2) Effect curves: effect = "YAP" -> HR along YAP-Rev with Stroma-A2 fixed at Q1/Q3;
#                   effect = "StromaA2" -> HR along Stroma-A2 with YAP-Rev fixed at Q1/Q3
# -----------------------------
build_effect_curves_petacc <- function(fit_obj, effect = c("YAP", "StromaA2")) {
  effect <- match.arg(effect)
  
  fit    <- fit_obj$fit
  df_z   <- fit_obj$data_z
  
  # x-axis variable
  x_var <- if (effect == "YAP") "YAPRev" else "StromalSubcluster_2"
  # interacting variable (Low/High strata)
  z_var <- if (effect == "YAP") "StromalSubcluster_2" else "YAPRev"
  
  # x grid over the observed range of the score
  x_range <- range(df_z[[x_var]], na.rm = TRUE)
  x_seq   <- seq(x_range[1], x_range[2], length.out = 200)
  
  # Q1 / Q3 of the other program
  z_q      <- stats::quantile(df_z[[z_var]], probs = c(0.25, 0.75), na.rm = TRUE)
  z_levels <- as.numeric(z_q)
  z_labels <- if (effect == "YAP") {
    c("Low Stroma-A2 (Q1)", "High Stroma-A2 (Q3)")
  } else {
    c("Low YAP-Rev (Q1)", "High YAP-Rev (Q3)")
  }
  
  grid <- expand.grid(
    x = x_seq,
    z = z_levels
  )
  
  # template newdata
  newdata <- df_z[rep(1, nrow(grid)), , drop = FALSE]
  
  # adjustment covariates fixed at reference level / mean
  adj_vars <- safe_adjustments(df_z)
  for (v in adj_vars) {
    if (is.factor(df_z[[v]])) {
      newdata[[v]] <- factor(levels(df_z[[v]])[1], levels = levels(df_z[[v]]))
    } else {
      newdata[[v]] <- mean(df_z[[v]], na.rm = TRUE)
    }
  }
  
  # fill the grid
  newdata[[x_var]] <- grid$x
  newdata[[z_var]] <- grid$z
  
  # reference patient: median of both programs
  refdata <- newdata[1, , drop = FALSE]
  refdata[[x_var]] <- stats::median(df_z[[x_var]], na.rm = TRUE)
  refdata[[z_var]] <- stats::median(df_z[[z_var]], na.rm = TRUE)
  
  ref_lp <- predict(fit, newdata = refdata, type = "lp", se.fit = FALSE)
  pred   <- predict(fit, newdata = newdata, type = "lp", se.fit = TRUE)
  
  lp_centered <- pred$fit - as.numeric(ref_lp)
  se          <- pred$se.fit
  
  HR      <- exp(lp_centered)
  HR_low  <- exp(lp_centered - 1.96 * se)
  HR_high <- exp(lp_centered + 1.96 * se)
  
  out <- data.frame(
    x      = grid$x,
    z      = factor(grid$z, levels = z_levels, labels = z_labels),
    HR     = HR,
    HR_low = HR_low,
    HR_high = HR_high
  )
  
  out$effect <- if (effect == "YAP") "YAP-Rev" else "Stroma-A2"
  out
}

# -----------------------------
# 3) PETACC8 plots
# -----------------------------
curves_PETACC_YAP <- build_effect_curves_petacc(fit_PETACC, effect = "YAP")
curves_PETACC_S2  <- build_effect_curves_petacc(fit_PETACC, effect = "StromaA2")

plot_interaction_petacc <- function(curves_df, effect_label) {
  
  x_lab <- if (effect_label == "YAP-Rev") {
    "YAP-Rev signature score (z-scaled)"
  } else {
    "Stroma-A2 signature score (z-scaled)"
  }
  
  ggplot(curves_df, aes(x = x, y = HR, color = z, fill = z)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
    geom_ribbon(aes(ymin = HR_low, ymax = HR_high),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(
      name   = "",
      values = c(
        "Low Stroma-A2 (Q1)"  = "#1B9E77",
        "High Stroma-A2 (Q3)" = "#D95F02",
        "Low YAP-Rev (Q1)"    = "#3333cc",
        "High YAP-Rev (Q3)"   = "#cc3399"
      )
    ) +
    scale_fill_manual(
      name   = "",
      values = c(
        "Low Stroma-A2 (Q1)"  = "#1B9E77",
        "High Stroma-A2 (Q3)" = "#D95F02",
        "Low YAP-Rev (Q1)"    = "#3333cc",
        "High YAP-Rev (Q3)"   = "#cc3399"
      )
    ) +
    labs(
      x        = x_lab,
      y        = "Relative hazard (HR)"
    ) +
    theme_pubclean(base_size = 15) +
    theme(
      legend.position  = "top",
      legend.title     = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 9)
    )
}

p_petacc_yap <- plot_interaction_petacc(curves_PETACC_YAP, effect_label = "YAP-Rev")
p_petacc_s2  <- plot_interaction_petacc(curves_PETACC_S2,  effect_label = "Stroma-A2")

# -----------------------------
# 4) one PDF per panel
# -----------------------------
pdf(file.path(out_dir, "Fig6d_interaction_YAPRev.pdf"), width = 4.5, height = 5)
print(p_petacc_yap)
dev.off()

pdf(file.path(out_dir, "Fig6d_interaction_StromaA2.pdf"), width = 4.5, height = 5)
print(p_petacc_s2)
dev.off()

# =====================================================================================================================
# Fig. 6e-g - Kaplan-Meier curves, 4 groups (PETACC8 medians of raw scores)
# =====================================================================================================================
suppressPackageStartupMessages(library(survminer))
# HR / CI / P table
hr_table_from_cox <- function(fit){
  s <- summary(fit)
  data.frame(
    level = rownames(s$coefficients),
    HR    = s$conf.int[,"exp(coef)"],
    LCL   = s$conf.int[,"lower .95"],
    UCL   = s$conf.int[,"upper .95"],
    p     = s$coefficients[,"Pr(>|z|)"],
    row.names = NULL, check.names = FALSE
  )
}

# 4 groups
make_groups_from_cuts <- function(df, s2_var, yap_var, s2_cut, yap_cut) {
  S2_hi  <- df[[s2_var]]  > s2_cut
  YAP_hi <- df[[yap_var]] > yap_cut
  
  grp <- dplyr::case_when(
    !S2_hi & !YAP_hi ~ "Low Stroma-A2 / Low YAP-Rev",
    S2_hi & !YAP_hi ~ "High Stroma-A2 / Low YAP-Rev",
    !S2_hi &  YAP_hi ~ "Low Stroma-A2 / High YAP-Rev",
    S2_hi &  YAP_hi ~ "High Stroma-A2 / High YAP-Rev",
    .default        = NA_character_
  )
  
  factor(
    grp,
    levels = c(
      "Low Stroma-A2 / Low YAP-Rev",
      "High Stroma-A2 / Low YAP-Rev",
      "Low Stroma-A2 / High YAP-Rev",
      "High Stroma-A2 / High YAP-Rev"
    )
  )
}

# palette in the order of the levels above
km_palette <- c(
  "#115566",  # Low/Low
  "#ad33ff",  # High Stroma-A2 / Low YAP-Rev
"#00e6ac",  # Low Stroma-A2 / High YAP-Rev
"#cc0066"   # High/High
)

# ---- one cohort, one set of cut-offs ----
plot_km_4groups_generic <- function(df,
                                    cohort_label,
                                    s2_var, yap_var,
                                    s2_cut, yap_cut,
                                    risk_table = TRUE,
                                    annotation_mode = c("both", "unadjusted", "adjusted")) {
  
  annotation_mode <- match.arg(annotation_mode)
  
  # complete cases
  keep <- complete.cases(df[, c("TTR_time", "etat_TTR", s2_var, yap_var), drop = FALSE])
  dd <- df[keep, , drop = FALSE]
  if (!nrow(dd)) return(NULL)
  
  dd$group4 <- make_groups_from_cuts(dd, s2_var, yap_var, s2_cut, yap_cut)
  
  fit_km <- survfit(Surv(TTR_time, etat_TTR) ~ group4, data = dd)
  cols   <- km_palette[seq_along(levels(dd$group4))]
  
  # log-rank P
  km_p <- surv_pvalue(fit_km, data = dd)$pval
  
  # fixed decimal format
  km_p_label <- dplyr::case_when(
    is.na(km_p)   ~ 'italic(P) == "NA"',
    km_p < 1e-4   ~ 'italic(P) < "0.0001"',
    km_p < 0.001  ~ 'italic(P) < "0.001"',
    TRUE          ~ paste0('italic(P) == "', formatC(km_p, format = "f", digits = 3), '"')
  )
  
  # P value added manually below
  p_km <- ggsurvplot(
    fit_km, data = dd,
    risk.table        = risk_table,
    risk.table.height = 0.3,
    risk.table.fontsize = 5,
    pval              = FALSE,
    pval.method       = FALSE,
    conf.int          = FALSE,
    legend.title      = "",
    legend.labs       = levels(dd$group4),
    legend            = c(0.2, 0.25),
    xlab              = "Years",
    ylab              = "Time-to-recurrence",
    xlim              = c(0, 10),
    break.x.by        = 2,
    palette           = cols,
    ggtheme           = theme_pubclean(base_size = 14),
    ylim              = c(0, 1)
  )
  
  # ---- group names in the risk table replaced by coloured bars ----
  
  rt <- p_km$table
  
  # empty y-axis labels
  rt <- rt +
    scale_y_discrete(labels = rep("", length(levels(dd$group4))))
  
  # colours reversed: the risk-table y axis lists the strata in reverse order (the unreversed version gave
  # inverted colour bars in the number-at-risk table)
  df_segments <- data.frame(
    strata = levels(dd$group4),
    col    = rev(cols)
  )
  
  # one segment per group
  rt <- rt +
    geom_segment(
      data        = df_segments,
      inherit.aes = FALSE,
      aes(
        x    = -0.5,
        xend = -0.3,
        y    = strata,
        yend = strata,
        colour = col
      ),
      linewidth = 3
    ) +
    scale_colour_identity()
  
  p_km$table <- rt
  
  
  
  # title
  p_km$plot <- p_km$plot +
    labs(title = cohort_label) +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 2),
      legend.text  = element_text(size = 14),
      legend.key.size = unit(1.2, "lines")
    )
  
  # log-rank P
  p_km$plot <- p_km$plot +
    annotate(
      "text",
      x     = 8,
      y     = 0.40,
      label = km_p_label,
      parse = TRUE,
      size  = 5
    )
  
  # unadjusted Cox
  fit_ua <- coxph(Surv(TTR_time, etat_TTR) ~ group4, data = dd)
  tb_ua  <- hr_table_from_cox(fit_ua)
  
  # adjusted Cox
  adj <- safe_adjustments(dd)
  if (length(adj)) {
    fmla_adj <- as.formula(
      paste0("Surv(TTR_time, etat_TTR) ~ group4 + ", paste(adj, collapse = " + "))
    )
    fit_ad <- coxph(fmla_adj, data = dd)
    tb_ad  <- hr_table_from_cox(fit_ad)
  } else {
    tb_ad <- NULL
  }
  
  pretty <- function(x) sprintf("%.2f", x)
  
  # P formatting in the Cox annotations
  fmt_p_vec <- function(p) {
    dplyr::case_when(
      is.na(p)  ~ "P = NA",
      p < 1e-4  ~ "P < 0.0001",
      p < 1e-3  ~ "P < 0.001",
      TRUE      ~ paste0("P = ", formatC(p, format = "f", digits = 3))
    )
  }
  
  ann <- NULL
  
  # Unadjusted
  if (annotation_mode %in% c("both", "unadjusted")) {
    lab_ua <- tb_ua %>%
      dplyr::filter(grepl("^group4", level)) %>%
      dplyr::mutate(
        level = sub("^group4", "", level),
        p_txt = fmt_p_vec(p),
        line  = paste0(
          level, ": HR ", pretty(HR),
          " (", pretty(LCL), "–", pretty(UCL), "), ",
          p_txt
        )
      ) %>%
      dplyr::pull(line)
    
    ann <- paste0(
      "Cox (unadjusted) vs Low Stroma-A2 / Low YAP-Rev:\n",
      paste(lab_ua, collapse = "\n")
    )
  }
  
  # Adjusted
  if (annotation_mode %in% c("both", "adjusted") && !is.null(tb_ad)) {
    lab_ad <- tb_ad %>%
      dplyr::filter(grepl("^group4", level)) %>%
      dplyr::mutate(
        level = sub("^group4", "", level),
        p_txt = fmt_p_vec(p),
        line  = paste0(
          level, ": HR ", pretty(HR),
          " (", pretty(LCL), "–", pretty(UCL), "), ",
          p_txt
        )
      ) %>%
      dplyr::pull(line)
    
    ann <- paste0(
      if (!is.null(ann)) paste0(ann, "\n\n") else "",
      "Cox (adjusted) vs Low Stroma-A2 / Low YAP-Rev:\n",
      paste(lab_ad, collapse = "\n")
    )
  }
  
  if (!is.null(ann)) {
    p_km$plot <- p_km$plot +
      annotate(
        "label",
        x      = Inf,
        y      = 0.05,
        hjust  = 1.05,
        vjust  = 0,
        label  = ann,
        size   = 4.5,
        color  = "black"
      )
  }
  
  list(
    gg       = p_km,
    hr_unadj = tb_ua,
    hr_adj   = tb_ad,
    n        = nrow(dd)
  )
}

# ---- non-standardised ssGSEA; PETACC8 medians applied unchanged to the 3 cohorts ----

# cut-offs: PETACC8 medians of the raw scores
s2_cut_petacc_raw  <- median(PETACC.tot$StromalSubcluster_2_raw, na.rm = TRUE)
yap_cut_petacc_raw <- median(PETACC.tot$YAPRev_raw,               na.rm = TRUE)


for (ck in list(list("PETACC", PETACC.tot, "PETACC8", "Fig6e"), list("IDEA", IDEA.tot, "IDEA-France", "Fig6f"),
                list("PRODIGE", PRODIGE.tot, "PRODIGE13", "Fig6g"))) {
  km <- plot_km_4groups_generic(df = ck[[2]], cohort_label = ck[[3]], s2_var = "StromalSubcluster_2_raw", yap_var = "YAPRev_raw",
                                s2_cut = s2_cut_petacc_raw, yap_cut = yap_cut_petacc_raw, annotation_mode = "unadjusted")
  print(km$hr_unadj)
  write.csv(km$hr_unadj, file.path(out_dir, paste0(ck[[4]], "_KM4_", ck[[1]], "_Cox_unadjusted.csv")), row.names = FALSE)
  pdf(file.path(out_dir, paste0(ck[[4]], "_KM4_", ck[[1]], ".pdf")), width = 10.2, height = 7)
  print(arrange_ggsurvplots(list(km$gg), ncol = 1, nrow = 1))
  dev.off()
}
