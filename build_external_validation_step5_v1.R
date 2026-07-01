args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
step4_dir <- file.path(project_dir, "manifests", "step 4")
step5_dir <- file.path(project_dir, "manifests", "step 5")
fig_dir <- file.path(step5_dir, "figures")
dir.create(step5_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("build_external_validation_step5_v1.R\n")
cat("Project directory:", project_dir, "\n")
cat("Step 4 directory:", step4_dir, "\n")
cat("Step 5 directory:", step5_dir, "\n")
cat("Figure directory:", fig_dir, "\n")
cat("R version:", R.version.string, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

required <- c("data.table", "ggplot2")
for (pkg in required) {
  cat("Package", pkg, "available:", requireNamespace(pkg, quietly = TRUE), "\n")
}
if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for code-generated figures")
library(data.table)
library(ggplot2)

classifier_rank <- data.table(
  table3_rank = 1:8,
  classifier_id = c(
    "IFNG_18",
    "PAM50",
    "ONCOTYPE_DX_21",
    "CYTOLYTIC_ACTIVITY_2",
    "MAMMAPRINT_70",
    "GGI_128_PROBE",
    "ENDOPREDICT_11",
    "GENE76"
  ),
  classifier_display = c(
    "A. IFN-gamma 18",
    "B. PAM50",
    "C. OncotypeDX",
    "D. Cytolytic activity",
    "E. MammaPrint",
    "F. GGI",
    "G. EndoPredict",
    "H. Wang 76-gene / GENE76"
  )
)

add_classifier_display <- function(dt) {
  dt <- copy(dt)
  dt[classifier_rank, on = "classifier_id", `:=`(
    table3_rank = i.table3_rank,
    classifier_display = i.classifier_display
  )]
  dt[is.na(classifier_display), `:=`(
    table3_rank = 999L,
    classifier_display = classifier_id
  )]
  dt[, classifier_axis := factor(classifier_display, levels = rev(classifier_rank$classifier_display))]
  dt[, classifier_facet := factor(classifier_display, levels = classifier_rank$classifier_display)]
  dt
}

set.seed(20260630)

score_file <- file.path(step4_dir, "classifier_score_matrix_v1.csv")
score_qc_file <- file.path(step4_dir, "classifier_score_qc_v1.csv")
prep_file <- file.path(step4_dir, "dataset_preprocessing_decisions_v1.csv")

cat("\nInput file check:\n")
for (f in c(score_file, score_qc_file, prep_file)) cat(file.exists(f), f, "\n")
if (!file.exists(score_file)) stop("Missing Step 4 score matrix: ", score_file)
if (!file.exists(score_qc_file)) stop("Missing Step 4 score QC: ", score_qc_file)
if (!file.exists(prep_file)) stop("Missing Step 4 preprocessing decisions: ", prep_file)

scores <- fread(score_file, encoding = "UTF-8", showProgress = FALSE)
score_qc <- fread(score_qc_file, encoding = "UTF-8", showProgress = FALSE)
prep <- fread(prep_file, encoding = "UTF-8", showProgress = FALSE)

cat("\nLoaded rows:\n")
cat("scores:", nrow(scores), "\n")
cat("score_qc:", nrow(score_qc), "\n")
cat("prep:", nrow(prep), "\n")

scores[, endpoint_binary_pcr1_rd0 := suppressWarnings(as.integer(endpoint_binary_pcr1_rd0))]
scores[, score_value := suppressWarnings(as.numeric(score_value))]
scores[, coverage_pct := suppressWarnings(as.numeric(coverage_pct))]

valid_scores <- scores[
  !is.na(endpoint_binary_pcr1_rd0) &
    endpoint_binary_pcr1_rd0 %in% c(0L, 1L) &
    is.finite(score_value) &
    score_status != "not_scored_insufficient_coverage"
]

cat("Valid score rows for performance analysis:", nrow(valid_scores), "\n")

data_profile <- scores[, .(
  n_score_rows = .N,
  n_finite_scores = sum(is.finite(score_value)),
  n_missing_scores = sum(!is.finite(score_value)),
  n_samples = uniqueN(gsm_sample_id),
  n_pcr = uniqueN(gsm_sample_id[endpoint_binary_pcr1_rd0 == 1L]),
  n_rd = uniqueN(gsm_sample_id[endpoint_binary_pcr1_rd0 == 0L]),
  pcr_prevalence = round(uniqueN(gsm_sample_id[endpoint_binary_pcr1_rd0 == 1L]) / uniqueN(gsm_sample_id), 4),
  score_status_values = paste(sort(unique(score_status)), collapse = ";")
), by = .(dataset_accession, matrix_label, platform, analysis_set)]
fwrite(data_profile, file.path(step5_dir, "data_profile_step5_v1.csv"))

safe_auc <- function(y, score) {
  ok <- y %in% c(0, 1) & is.finite(score)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

average_precision <- function(y, score) {
  ok <- y %in% c(0, 1) & is.finite(score)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  n_pos <- sum(y == 1L)
  if (n_pos == 0L || length(y) == 0L) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  precision <- tp / (tp + fp)
  mean(precision[y == 1L])
}

auc_boot_ci <- function(y, score, sign = 1, b = 500L) {
  ok <- y %in% c(0, 1) & is.finite(score)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok]) * sign
  if (length(y) < 10L || sum(y == 1L) < 3L || sum(y == 0L) < 3L) {
    return(c(NA_real_, NA_real_, 0L))
  }
  vals <- numeric(0)
  n <- length(y)
  for (i in seq_len(b)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(y[idx])) < 2L) next
    vals <- c(vals, safe_auc(y[idx], score[idx]))
  }
  vals <- vals[is.finite(vals)]
  if (length(vals) < 50L) return(c(NA_real_, NA_real_, length(vals)))
  c(as.numeric(stats::quantile(vals, 0.025, na.rm = TRUE)),
    as.numeric(stats::quantile(vals, 0.975, na.rm = TRUE)),
    length(vals))
}

wilcox_auc_p <- function(y, score) {
  ok <- y %in% c(0, 1) & is.finite(score)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  if (sum(y == 1L) < 2L || sum(y == 0L) < 2L) return(NA_real_)
  tryCatch({
    stats::wilcox.test(score[y == 1L], score[y == 0L], exact = FALSE)$p.value
  }, error = function(e) NA_real_)
}

logit <- function(p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  log(p / (1 - p))
}

fit_glm_safe <- function(dt) {
  dt <- dt[is.finite(score_value) & endpoint_binary_pcr1_rd0 %in% c(0L, 1L)]
  if (nrow(dt) < 10L || length(unique(dt$endpoint_binary_pcr1_rd0)) < 2L) {
    return(list(
      fit = NULL, intercept = NA_real_, beta = NA_real_, status = "not_fit_insufficient_data"
    ))
  }
  fit <- tryCatch(
    suppressWarnings(stats::glm(endpoint_binary_pcr1_rd0 ~ score_value, data = dt, family = stats::binomial())),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(fit = NULL, intercept = NA_real_, beta = NA_real_, status = paste0("glm_error: ", conditionMessage(fit))))
  }
  cf <- stats::coef(fit)
  status <- if (any(!is.finite(cf))) "fit_nonfinite_coefficient" else "fit_ok"
  list(fit = fit, intercept = unname(cf[1]), beta = unname(cf[2]), status = status)
}

safe_cal_intercept <- function(y, pred) {
  ok <- y %in% c(0, 1) & is.finite(pred)
  y <- as.integer(y[ok])
  pred <- pmin(pmax(pred[ok], 1e-6), 1 - 1e-6)
  if (length(y) < 10L || length(unique(y)) < 2L) return(NA_real_)
  out <- tryCatch(
    suppressWarnings(stats::glm(y ~ 1 + offset(logit(pred)), family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(out)) return(NA_real_)
  unname(stats::coef(out)[1])
}

safe_cal_slope <- function(y, pred) {
  ok <- y %in% c(0, 1) & is.finite(pred)
  y <- as.integer(y[ok])
  lp <- logit(pred[ok])
  if (length(y) < 10L || length(unique(y)) < 2L || stats::sd(lp) == 0) {
    return(c(NA_real_, NA_real_))
  }
  out <- tryCatch(
    suppressWarnings(stats::glm(y ~ lp, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(out)) return(c(NA_real_, NA_real_))
  cf <- stats::coef(out)
  c(unname(cf[1]), unname(cf[2]))
}

bin_calibration <- function(dt, n_bins = 5L) {
  dt <- dt[is.finite(predicted_probability) & endpoint_binary_pcr1_rd0 %in% c(0L, 1L)]
  if (nrow(dt) < 10L || length(unique(dt$endpoint_binary_pcr1_rd0)) < 2L) return(data.table())
  # Quantile bins; fall back to equal-width if duplicated breaks collapse.
  probs <- seq(0, 1, length.out = n_bins + 1L)
  br <- unique(as.numeric(stats::quantile(dt$predicted_probability, probs, na.rm = TRUE)))
  if (length(br) < 3L) {
    br <- unique(seq(min(dt$predicted_probability), max(dt$predicted_probability), length.out = n_bins + 1L))
  }
  if (length(br) < 3L) return(data.table())
  dt[, calibration_bin := cut(predicted_probability, breaks = br, include.lowest = TRUE, labels = FALSE)]
  dt[!is.na(calibration_bin), .(
    n = .N,
    mean_predicted_probability = mean(predicted_probability),
    observed_pcr_rate = mean(endpoint_binary_pcr1_rd0),
    n_pcr = sum(endpoint_binary_pcr1_rd0 == 1L)
  ), by = calibration_bin]
}

discovery_label <- "GSE25066_GPL96"
disc_models <- list()
model_rows <- list()
for (cid in sort(unique(valid_scores$classifier_id))) {
  train <- valid_scores[matrix_label == discovery_label & classifier_id == cid]
  fit_obj <- fit_glm_safe(train)
  direction_sign <- ifelse(is.finite(fit_obj$beta) && fit_obj$beta < 0, -1, 1)
  disc_models[[cid]] <- fit_obj
  disc_models[[cid]]$direction_sign <- direction_sign
  model_rows[[length(model_rows) + 1L]] <- data.table(
    classifier_id = cid,
    discovery_matrix_label = discovery_label,
    n_discovery = nrow(train),
    n_discovery_pcr = sum(train$endpoint_binary_pcr1_rd0 == 1L),
    n_discovery_rd = sum(train$endpoint_binary_pcr1_rd0 == 0L),
    intercept = fit_obj$intercept,
    beta = fit_obj$beta,
    discovery_direction_sign = direction_sign,
    fit_status = fit_obj$status,
    score_status_values = paste(sort(unique(train$score_status)), collapse = ";")
  )
}
model_dt <- rbindlist(model_rows, fill = TRUE)
fwrite(model_dt, file.path(step5_dir, "discovery_logistic_models_v1.csv"))

perf_rows <- list()
cal_rows <- list()
pred_rows <- list()
cal_curve_rows <- list()

groups <- unique(valid_scores[, .(
  dataset_accession, matrix_label, platform, analysis_set, classifier_id,
  classifier_score_name, score_status
)])
setorder(groups, classifier_id, matrix_label, score_status)

for (i in seq_len(nrow(groups))) {
  g <- groups[i]
  dt <- valid_scores[
    dataset_accession == g$dataset_accession &
      matrix_label == g$matrix_label &
      platform == g$platform &
      analysis_set == g$analysis_set &
      classifier_id == g$classifier_id &
      score_status == g$score_status
  ]
  y <- dt$endpoint_binary_pcr1_rd0
  s <- dt$score_value
  model <- disc_models[[g$classifier_id]]
  sign <- if (is.null(model$direction_sign)) 1 else model$direction_sign
  auc_raw <- safe_auc(y, s)
  auc_aligned <- safe_auc(y, sign * s)
  ci <- auc_boot_ci(y, s, sign = sign, b = 500L)
  ap <- average_precision(y, sign * s)
  p_w <- wilcox_auc_p(y, sign * s)
  perf_rows[[length(perf_rows) + 1L]] <- cbind(g, data.table(
    n = nrow(dt),
    n_pcr = sum(y == 1L),
    n_rd = sum(y == 0L),
    pcr_prevalence = round(mean(y == 1L), 4),
    n_unique_patients = uniqueN(dt$patient_id),
    raw_auc_high_score_predicts_pcr = auc_raw,
    discovery_aligned_auc = auc_aligned,
    auc_ci_low_bootstrap = ci[1],
    auc_ci_high_bootstrap = ci[2],
    auc_bootstrap_valid_replicates = as.integer(ci[3]),
    average_precision_pr_auc = ap,
    pr_auc_baseline_pcr_prevalence = round(mean(y == 1L), 4),
    pr_auc_minus_pcr_prevalence = ap - mean(y == 1L),
    wilcoxon_p_value_aligned_score = p_w,
    discovery_direction_sign = sign,
    performance_role = ifelse(g$matrix_label == discovery_label, "discovery_calibration_base",
                              ifelse(g$analysis_set == "core_full_benchmark", "external_core_validation",
                                     ifelse(g$analysis_set == "large_stress_test", "large_stress_test",
                                            "limited_or_subtype_exploratory")))
  ))
  if (!is.null(model$fit) && is.finite(model$intercept) && is.finite(model$beta)) {
    lp <- model$intercept + model$beta * s
    pred <- 1 / (1 + exp(-lp))
    pred <- pmin(pmax(pred, 1e-6), 1 - 1e-6)
    pred_dt <- copy(dt)
    pred_dt[, predicted_probability := pred]
    pred_dt[, discovery_model_status := model$status]
    pred_dt[, discovery_beta := model$beta]
    pred_rows[[length(pred_rows) + 1L]] <- pred_dt
    intercept_only <- safe_cal_intercept(y, pred)
    cal_pair <- safe_cal_slope(y, pred)
    observed <- sum(y == 1L)
    expected <- sum(pred)
    cal_rows[[length(cal_rows) + 1L]] <- cbind(g, data.table(
      n = nrow(dt),
      n_pcr = observed,
      n_rd = sum(y == 0L),
      mean_predicted_probability = mean(pred),
      observed_pcr_rate = mean(y == 1L),
      brier_score = mean((pred - y)^2),
      calibration_intercept_offset = intercept_only,
      calibration_model_intercept = cal_pair[1],
      calibration_slope = cal_pair[2],
      observed_expected_ratio = ifelse(expected > 0, observed / expected, NA_real_),
      discovery_model_status = model$status,
      calibration_role = ifelse(g$matrix_label == discovery_label, "apparent_discovery",
                                ifelse(g$analysis_set == "core_full_benchmark", "external_core_validation",
                                       ifelse(g$analysis_set == "large_stress_test", "large_stress_test",
                                              "limited_or_subtype_exploratory")))
    ))
    bins <- bin_calibration(pred_dt, n_bins = 5L)
    if (nrow(bins)) {
      bins <- cbind(g, bins)
      cal_curve_rows[[length(cal_curve_rows) + 1L]] <- bins
    }
  } else {
    cal_rows[[length(cal_rows) + 1L]] <- cbind(g, data.table(
      n = nrow(dt), n_pcr = sum(y == 1L), n_rd = sum(y == 0L),
      mean_predicted_probability = NA_real_, observed_pcr_rate = mean(y == 1L),
      brier_score = NA_real_, calibration_intercept_offset = NA_real_,
      calibration_model_intercept = NA_real_, calibration_slope = NA_real_,
      observed_expected_ratio = NA_real_,
      discovery_model_status = ifelse(is.null(model$status), "no_discovery_model", model$status),
      calibration_role = "not_calibrated"
    ))
  }
}

performance <- rbindlist(perf_rows, fill = TRUE)
calibration <- rbindlist(cal_rows, fill = TRUE)
predictions <- if (length(pred_rows)) rbindlist(pred_rows, fill = TRUE) else data.table()
cal_curve <- if (length(cal_curve_rows)) rbindlist(cal_curve_rows, fill = TRUE) else data.table()

performance[, wilcoxon_q_value_aligned_score := p.adjust(wilcoxon_p_value_aligned_score, method = "BH")]

fwrite(performance, file.path(step5_dir, "classifier_performance_by_dataset_v1.csv"))
fwrite(calibration, file.path(step5_dir, "classifier_calibration_by_dataset_v1.csv"))
fwrite(predictions, file.path(step5_dir, "classifier_predicted_probabilities_v1.csv"))
fwrite(cal_curve, file.path(step5_dir, "calibration_curve_points_v1.csv"))

external_core <- performance[performance_role == "external_core_validation" & score_status == "scored_preferred_coverage"]
meta_summary <- external_core[, .(
  n_external_core_datasets = .N,
  median_external_auc = median(discovery_aligned_auc, na.rm = TRUE),
  min_external_auc = min(discovery_aligned_auc, na.rm = TRUE),
  max_external_auc = max(discovery_aligned_auc, na.rm = TRUE),
  n_auc_ge_0_60 = sum(discovery_aligned_auc >= 0.60, na.rm = TRUE),
  n_auc_ge_0_65 = sum(discovery_aligned_auc >= 0.65, na.rm = TRUE),
  median_pr_auc = median(average_precision_pr_auc, na.rm = TRUE),
  median_pcr_prevalence = median(pcr_prevalence, na.rm = TRUE),
  datasets = paste(matrix_label, collapse = ";")
), by = .(classifier_id, classifier_score_name)]
meta_summary[, external_transferability_tier := fifelse(
  n_external_core_datasets >= 3 & median_external_auc >= 0.65, "strong_candidate",
  fifelse(n_external_core_datasets >= 3 & median_external_auc >= 0.60, "moderate_candidate",
          fifelse(n_external_core_datasets >= 2 & median_external_auc >= 0.55, "weak_or_context_dependent", "poor_or_insufficient"))
)]
setorder(meta_summary, -median_external_auc, -n_auc_ge_0_60)
fwrite(meta_summary, file.path(step5_dir, "classifier_meta_summary_v1.csv"))

rank_dt <- copy(meta_summary)
rank_dt[, rank_external_core := frank(-median_external_auc, ties.method = "min")]
setorder(rank_dt, rank_external_core, classifier_id)
fwrite(rank_dt, file.path(step5_dir, "classifier_external_validation_rank_v1.csv"))

failure_mode <- performance[, .(
  n_comparisons = .N,
  median_auc = median(discovery_aligned_auc, na.rm = TRUE),
  min_auc = min(discovery_aligned_auc, na.rm = TRUE),
  max_auc = max(discovery_aligned_auc, na.rm = TRUE),
  n_auc_below_0_55 = sum(discovery_aligned_auc < 0.55, na.rm = TRUE),
  n_low_coverage = sum(score_status == "scored_low_coverage_exploratory"),
  n_not_primary = sum(performance_role %in% c("large_stress_test", "limited_or_subtype_exploratory")),
  n_total_samples = sum(n)
), by = .(classifier_id, platform, analysis_set, score_status)]
failure_mode[, failure_mode_flag := fifelse(
  n_auc_below_0_55 > 0 | n_low_coverage > 0 | analysis_set == "limited_subtype_panel_or_gene_space",
  "potential_failure_mode", "no_major_flag_in_step5"
)]
fwrite(failure_mode, file.path(step5_dir, "classifier_failure_mode_audit_v1.csv"))

subtype_dt <- valid_scores[
  score_status == "scored_preferred_coverage" &
    subtype_group %in% c(
      "TNBC",
      "TNBC_by_study",
      "TNBC_or_ER_PR_HER2_negative",
      "ER_negative_or_HR_negative_HER2_negative",
      "HR_negative_HER2_unknown",
      "HER2_positive",
      "HR_positive_HER2_negative"
    )
]
subtype_rows <- list()
if (nrow(subtype_dt)) {
  sg <- unique(subtype_dt[, .(dataset_accession, matrix_label, platform, analysis_set, subtype_group, classifier_id, classifier_score_name)])
  for (i in seq_len(nrow(sg))) {
    g <- sg[i]
    dt <- subtype_dt[
      dataset_accession == g$dataset_accession &
        matrix_label == g$matrix_label &
        subtype_group == g$subtype_group &
        classifier_id == g$classifier_id
    ]
    if (length(unique(dt$endpoint_binary_pcr1_rd0)) < 2L || nrow(dt) < 15L) next
    sign <- disc_models[[g$classifier_id]]$direction_sign
    subtype_rows[[length(subtype_rows) + 1L]] <- cbind(g, data.table(
      n = nrow(dt),
      n_pcr = sum(dt$endpoint_binary_pcr1_rd0 == 1L),
      n_rd = sum(dt$endpoint_binary_pcr1_rd0 == 0L),
      pcr_prevalence = mean(dt$endpoint_binary_pcr1_rd0 == 1L),
      discovery_aligned_auc = safe_auc(dt$endpoint_binary_pcr1_rd0, sign * dt$score_value),
      average_precision_pr_auc = average_precision(dt$endpoint_binary_pcr1_rd0, sign * dt$score_value)
    ))
  }
}
subtype_perf <- if (length(subtype_rows)) rbindlist(subtype_rows, fill = TRUE) else data.table()
fwrite(subtype_perf, file.path(step5_dir, "classifier_subtype_performance_v1.csv"))

write_plan <- function() {
  plan <- c(
    "# Step 5 external validation analysis plan v1",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## Objective",
    "",
    "Evaluate whether pre-specified transcriptomic classifier scores from Step 4 transfer across independent breast cancer neoadjuvant chemotherapy cohorts for pCR/RD prediction.",
    "",
    "## Primary endpoint",
    "",
    "`endpoint_binary_pcr1_rd0`, where pCR = 1 and residual disease = 0.",
    "",
    "## Primary score rule",
    "",
    "Discrimination is reported as discovery-aligned AUC. Direction is learned only in the discovery/calibration base GSE25066; external cohorts are not used to flip score direction.",
    "",
    "## Primary analysis set",
    "",
    "External core validation cohorts with `analysis_set == core_full_benchmark` and `score_status == scored_preferred_coverage`.",
    "",
    "## Secondary/exploratory analysis sets",
    "",
    "- GSE194040 platforms are reported as large stress tests.",
    "- Limited gene-space/panel datasets are reported separately.",
    "- Low-coverage scores are flagged as exploratory.",
    "",
    "## Metrics",
    "",
    "- ROC-AUC with bootstrap 95% CI.",
    "- Average precision / PR-AUC with pCR prevalence as baseline.",
    "- Discovery-trained logistic calibration: Brier score, calibration intercept, calibration slope, observed/expected ratio.",
    "- Failure-mode summaries by platform, analysis set, and score coverage status.",
    "",
    "## Leakage control",
    "",
    "No external validation outcome is used for direction flipping, score scaling, cutoff selection, calibration fitting, or model tuning."
  )
  writeLines(plan, file.path(step5_dir, "validation_analysis_plan_v1.md"), useBytes = TRUE)
}
write_plan()

make_figures <- function() {
  figure_rows <- list()
  fig_idx <- 1L
  register_fig <- function(filename, status, notes) {
    figure_rows[[length(figure_rows) + 1L]] <<- data.table(
      figure_file = filename,
      exists = file.exists(file.path(fig_dir, filename)),
      status = status,
      notes = notes
    )
  }
  save_figure_pair <- function(plot, stem, width, height, notes) {
    png_file <- paste0(stem, ".png")
    pdf_file <- paste0(stem, ".pdf")
    ggsave(file.path(fig_dir, png_file), plot, width = width, height = height, dpi = 300)
    register_fig(png_file, "generated", notes)
    ggsave(file.path(fig_dir, pdf_file), plot, width = width, height = height, device = "pdf")
    register_fig(pdf_file, "generated", paste0(notes, " PDF version"))
  }

  heat_dt <- performance[score_status %in% c("scored_preferred_coverage", "scored_low_coverage_exploratory")]
  if (nrow(heat_dt)) {
    heat_dt <- add_classifier_display(heat_dt)
    heat_dt[, dataset_display := paste0(matrix_label, "\n", analysis_set)]
    p <- ggplot(heat_dt, aes(x = dataset_display, y = classifier_axis, fill = discovery_aligned_auc)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = ifelse(is.finite(discovery_aligned_auc), sprintf("%.2f", discovery_aligned_auc), "")), size = 2.6) +
      scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0.5,
                           limits = c(0.25, 0.85), na.value = "grey90", name = "AUC") +
      labs(title = "External validation discrimination by cohort and classifier",
           subtitle = "AUC direction aligned using GSE25066 only; low-coverage cells are included but flagged in tables",
           x = "Dataset / analysis set", y = "Classifier") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
            panel.grid = element_blank())
    save_figure_pair(p, "fig_auc_heatmap_v1", 13, 6.8, "AUC heatmap generated from classifier_performance_by_dataset_v1.csv")
  } else {
    register_fig("fig_auc_heatmap_v1.png", "not_generated", "No performance rows available")
    register_fig("fig_auc_heatmap_v1.pdf", "not_generated", "No performance rows available")
  }

  pr_dt <- performance[score_status %in% c("scored_preferred_coverage", "scored_low_coverage_exploratory")]
  if (nrow(pr_dt)) {
    pr_dt <- add_classifier_display(pr_dt)
    pr_dt[, dataset_display := paste0(matrix_label, "\n", analysis_set)]
    pr_delta_limit <- max(abs(pr_dt$pr_auc_minus_pcr_prevalence[is.finite(pr_dt$pr_auc_minus_pcr_prevalence)]), na.rm = TRUE)
    if (!is.finite(pr_delta_limit) || pr_delta_limit == 0) pr_delta_limit <- 0.05
    p <- ggplot(pr_dt, aes(x = dataset_display, y = classifier_axis, fill = pr_auc_minus_pcr_prevalence)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = ifelse(is.finite(pr_auc_minus_pcr_prevalence), sprintf("%+.2f", pr_auc_minus_pcr_prevalence), "")), size = 2.6) +
      scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#238B45", midpoint = 0,
                           limits = c(-pr_delta_limit, pr_delta_limit), na.value = "grey90",
                           name = "PR-AUC -\npCR prevalence") +
      labs(title = "Precision-recall gain over cohort pCR prevalence",
           subtitle = "Average precision minus cohort-specific pCR prevalence; zero indicates the prevalence baseline",
           x = "Dataset / analysis set", y = "Classifier") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
            panel.grid = element_blank())
    save_figure_pair(p, "fig_pr_auc_heatmap_v1", 13, 6.8, "PR-AUC minus pCR prevalence heatmap generated from classifier_performance_by_dataset_v1.csv")
  } else {
    register_fig("fig_pr_auc_heatmap_v1.png", "not_generated", "No PR-AUC rows available")
    register_fig("fig_pr_auc_heatmap_v1.pdf", "not_generated", "No PR-AUC rows available")
  }

  forest_dt <- performance[
    performance_role %in% c("external_core_validation", "large_stress_test") &
      score_status == "scored_preferred_coverage" &
      is.finite(discovery_aligned_auc)
  ]
  if (nrow(forest_dt)) {
    forest_dt <- add_classifier_display(forest_dt)
    forest_dt[, dataset_display := paste0(matrix_label, " (", analysis_set, ")")]
    p <- ggplot(forest_dt, aes(x = discovery_aligned_auc, y = reorder(dataset_display, discovery_aligned_auc))) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey45") +
      geom_errorbarh(aes(xmin = auc_ci_low_bootstrap, xmax = auc_ci_high_bootstrap), height = 0.15, color = "grey35") +
      geom_point(aes(color = analysis_set, size = n), alpha = 0.9) +
      facet_wrap(~ classifier_facet, scales = "free_y", ncol = 2) +
      scale_x_continuous(limits = c(0.25, 0.9), breaks = seq(0.3, 0.9, by = 0.1)) +
      labs(title = "Transferability forest plot",
           subtitle = "Preferred-coverage external core and large stress-test cohorts only; 95% bootstrap CI",
           x = "Discovery-aligned ROC-AUC", y = "Dataset") +
      theme_bw(base_size = 9) +
      theme(legend.position = "bottom")
    save_figure_pair(p, "fig_transferability_forest_v1", 12, 10, "Forest plot generated from classifier_performance_by_dataset_v1.csv")
  } else {
    register_fig("fig_transferability_forest_v1.png", "not_generated", "No preferred-coverage external/stress AUC rows available")
    register_fig("fig_transferability_forest_v1.pdf", "not_generated", "No preferred-coverage external/stress AUC rows available")
  }

  cal_plot_dt <- cal_curve[
    analysis_set %in% c("core_full_benchmark", "large_stress_test") &
      score_status == "scored_preferred_coverage"
  ]
  if (nrow(cal_plot_dt)) {
    cal_plot_dt <- add_classifier_display(cal_plot_dt)
    p <- ggplot(cal_plot_dt, aes(x = mean_predicted_probability, y = observed_pcr_rate,
                                 group = matrix_label, color = analysis_set)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
      geom_line(alpha = 0.45) +
      geom_point(aes(size = n), alpha = 0.75) +
      facet_wrap(~ classifier_facet, ncol = 2) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(title = "Calibration curves from GSE25066-trained logistic models",
           subtitle = "Points are within-cohort quantile bins; external cohorts are not recalibrated",
           x = "Mean predicted pCR probability", y = "Observed pCR rate") +
      theme_bw(base_size = 9) +
      theme(legend.position = "bottom")
    save_figure_pair(p, "fig_calibration_grid_v1", 11, 10, "Calibration grid generated from calibration_curve_points_v1.csv")
  } else {
    register_fig("fig_calibration_grid_v1.png", "not_generated", "No calibration curve points available")
    register_fig("fig_calibration_grid_v1.pdf", "not_generated", "No calibration curve points available")
  }

  fail_dt <- failure_mode[is.finite(median_auc)]
  if (nrow(fail_dt)) {
    fail_dt <- add_classifier_display(fail_dt)
    p <- ggplot(fail_dt, aes(x = platform, y = classifier_axis, fill = median_auc)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = sprintf("%.2f", median_auc)), size = 2.7) +
      facet_grid(score_status ~ analysis_set, scales = "free_x", space = "free_x") +
      scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0.5,
                           limits = c(0.25, 0.85), na.value = "grey90", name = "Median AUC") +
      labs(title = "Failure-mode audit by platform, analysis set, and coverage status",
           x = "Platform", y = "Classifier") +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank())
    save_figure_pair(p, "fig_failure_mode_heatmap_v1", 13, 8, "Failure-mode heatmap generated from classifier_failure_mode_audit_v1.csv")
  } else {
    register_fig("fig_failure_mode_heatmap_v1.png", "not_generated", "No failure-mode rows available")
    register_fig("fig_failure_mode_heatmap_v1.pdf", "not_generated", "No failure-mode rows available")
  }

  fig_manifest <- rbindlist(figure_rows, fill = TRUE)
  fwrite(fig_manifest, file.path(step5_dir, "figure_manifest_v1.csv"))
  fig_manifest
}

fig_manifest <- make_figures()

summary_lines <- c(
  "# Step 5 external validation execution notes v1",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Data profile",
  "",
  paste0("- Score rows loaded: ", nrow(scores)),
  paste0("- Valid score rows for performance analysis: ", nrow(valid_scores)),
  paste0("- Datasets/matrices profiled: ", nrow(data_profile)),
  "",
  "## Outputs",
  "",
  "- validation_analysis_plan_v1.md",
  "- data_profile_step5_v1.csv",
  "- discovery_logistic_models_v1.csv",
  "- classifier_performance_by_dataset_v1.csv",
  "- classifier_calibration_by_dataset_v1.csv",
  "- classifier_predicted_probabilities_v1.csv",
  "- calibration_curve_points_v1.csv",
  "- classifier_meta_summary_v1.csv",
  "- classifier_external_validation_rank_v1.csv",
  "- classifier_failure_mode_audit_v1.csv",
  "- classifier_subtype_performance_v1.csv",
  "- figure_manifest_v1.csv",
  "- figures/*.png",
  "",
  "## Figure generation status",
  "",
  paste(capture.output(print(fig_manifest)), collapse = "\n"),
  "",
  "## External core ranking",
  "",
  if (nrow(rank_dt)) paste(capture.output(print(rank_dt[, .(rank_external_core, classifier_id, n_external_core_datasets, median_external_auc, min_external_auc, max_external_auc, external_transferability_tier)])), collapse = "\n") else "No external core rank rows.",
  "",
  "## Integrity note",
  "",
  "All metrics and figures in Step 5 were generated by this R script from Step 4 outputs. No result values were manually filled."
)
writeLines(summary_lines, file.path(step5_dir, "step5_external_validation_notes_v1.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(step5_dir, "R_session_info_step5_v1.txt"))

cat("\nOutput rows:\n")
cat("data_profile_step5_v1.csv:", nrow(data_profile), "\n")
cat("discovery_logistic_models_v1.csv:", nrow(model_dt), "\n")
cat("classifier_performance_by_dataset_v1.csv:", nrow(performance), "\n")
cat("classifier_calibration_by_dataset_v1.csv:", nrow(calibration), "\n")
cat("classifier_predicted_probabilities_v1.csv:", nrow(predictions), "\n")
cat("calibration_curve_points_v1.csv:", nrow(cal_curve), "\n")
cat("classifier_meta_summary_v1.csv:", nrow(meta_summary), "\n")
cat("classifier_external_validation_rank_v1.csv:", nrow(rank_dt), "\n")
cat("classifier_failure_mode_audit_v1.csv:", nrow(failure_mode), "\n")
cat("classifier_subtype_performance_v1.csv:", nrow(subtype_perf), "\n")
cat("figure_manifest_v1.csv:", nrow(fig_manifest), "\n")

missing_generated <- fig_manifest[status == "not_generated"]
if (nrow(missing_generated)) {
  cat("\nSome figures were not generated. See figure_manifest_v1.csv.\n")
} else {
  cat("\nAll planned Step 5 figures were generated by code.\n")
}
cat("\nDone.\n")
