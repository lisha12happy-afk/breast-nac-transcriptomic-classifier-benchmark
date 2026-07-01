args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
step5_dir <- file.path(project_dir, "manifests", "step 5")
step6_dir <- file.path(project_dir, "manifests", "step 6")
fig_dir <- file.path(step6_dir, "figures")
dir.create(step6_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("build_sensitivity_step6_v1.R\n")
cat("Project directory:", project_dir, "\n")
cat("Step 5 directory:", step5_dir, "\n")
cat("Step 6 directory:", step6_dir, "\n")
cat("Figure directory:", fig_dir, "\n")
cat("R version:", R.version.string, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for PDF figures")
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

input_files <- c(
  performance = file.path(step5_dir, "classifier_performance_by_dataset_v1.csv"),
  calibration = file.path(step5_dir, "classifier_calibration_by_dataset_v1.csv"),
  subtype = file.path(step5_dir, "classifier_subtype_performance_v1.csv"),
  step5_rank = file.path(step5_dir, "classifier_external_validation_rank_v1.csv"),
  failure_mode = file.path(step5_dir, "classifier_failure_mode_audit_v1.csv")
)

cat("Input file check:\n")
for (nm in names(input_files)) cat(nm, file.exists(input_files[[nm]]), input_files[[nm]], "\n")
if (!all(file.exists(input_files))) stop("Missing one or more Step 5 input files.")

perf <- fread(input_files[["performance"]], encoding = "UTF-8", showProgress = FALSE)
cal <- fread(input_files[["calibration"]], encoding = "UTF-8", showProgress = FALSE)
subtype <- fread(input_files[["subtype"]], encoding = "UTF-8", showProgress = FALSE)
step5_rank <- fread(input_files[["step5_rank"]], encoding = "UTF-8", showProgress = FALSE)
failure_mode <- fread(input_files[["failure_mode"]], encoding = "UTF-8", showProgress = FALSE)

numeric_cols_perf <- c("n", "n_pcr", "n_rd", "pcr_prevalence", "discovery_aligned_auc",
                       "auc_ci_low_bootstrap", "auc_ci_high_bootstrap", "average_precision_pr_auc",
                       "pr_auc_baseline_pcr_prevalence")
for (cc in intersect(numeric_cols_perf, names(perf))) perf[, (cc) := suppressWarnings(as.numeric(get(cc)))]
for (cc in c("n", "n_pcr", "n_rd", "brier_score", "calibration_intercept_offset",
             "calibration_model_intercept", "calibration_slope", "observed_expected_ratio",
             "observed_pcr_rate", "mean_predicted_probability")) {
  if (cc %in% names(cal)) cal[, (cc) := suppressWarnings(as.numeric(get(cc)))]
}
for (cc in c("n", "n_pcr", "n_rd", "pcr_prevalence", "discovery_aligned_auc", "average_precision_pr_auc")) {
  if (cc %in% names(subtype)) subtype[, (cc) := suppressWarnings(as.numeric(get(cc)))]
}

cat("\nLoaded rows:\n")
cat("performance:", nrow(perf), "\n")
cat("calibration:", nrow(cal), "\n")
cat("subtype:", nrow(subtype), "\n")
cat("step5_rank:", nrow(step5_rank), "\n")
cat("failure_mode:", nrow(failure_mode), "\n")

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}
safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

scenario_defs <- data.table(
  scenario_id = c(
    "primary_core_preferred",
    "core_include_low_coverage",
    "core_microarray_only_preferred",
    "core_exclude_gpl96_preferred",
    "stress_test_gse194040_preferred",
    "all_external_non_discovery_preferred",
    "all_external_non_discovery_with_low_coverage"
  ),
  scenario_label = c(
    "Primary external core, preferred coverage",
    "External core, include low-coverage exploratory",
    "External core microarray only, preferred coverage",
    "External core excluding GPL96, preferred coverage",
    "GSE194040 large stress-test only, preferred coverage",
    "All non-discovery matrices, preferred coverage",
    "All non-discovery matrices, include low coverage"
  ),
  scenario_definition = c(
    "performance_role == external_core_validation; score_status == scored_preferred_coverage",
    "performance_role == external_core_validation; score_status in preferred or low-coverage exploratory",
    "performance_role == external_core_validation; score_status == preferred; platform in GPL96/GPL570/GPL571",
    "performance_role == external_core_validation; score_status == preferred; platform != GPL96",
    "performance_role == large_stress_test; score_status == preferred",
    "performance_role != discovery_calibration_base; score_status == preferred",
    "performance_role != discovery_calibration_base; score_status in preferred or low-coverage exploratory"
  ),
  primary_for_claim = c("yes", rep("no", 6))
)
fwrite(scenario_defs, file.path(step6_dir, "sensitivity_analysis_plan_v1.csv"))

scenario_filter <- function(dt, scenario_id) {
  if (scenario_id == "primary_core_preferred") {
    return(dt[performance_role == "external_core_validation" & score_status == "scored_preferred_coverage"])
  }
  if (scenario_id == "core_include_low_coverage") {
    return(dt[performance_role == "external_core_validation" & score_status %in% c("scored_preferred_coverage", "scored_low_coverage_exploratory")])
  }
  if (scenario_id == "core_microarray_only_preferred") {
    return(dt[performance_role == "external_core_validation" & score_status == "scored_preferred_coverage" & platform %in% c("GPL96", "GPL570", "GPL571")])
  }
  if (scenario_id == "core_exclude_gpl96_preferred") {
    return(dt[performance_role == "external_core_validation" & score_status == "scored_preferred_coverage" & platform != "GPL96"])
  }
  if (scenario_id == "stress_test_gse194040_preferred") {
    return(dt[performance_role == "large_stress_test" & score_status == "scored_preferred_coverage"])
  }
  if (scenario_id == "all_external_non_discovery_preferred") {
    return(dt[performance_role != "discovery_calibration_base" & score_status == "scored_preferred_coverage"])
  }
  if (scenario_id == "all_external_non_discovery_with_low_coverage") {
    return(dt[performance_role != "discovery_calibration_base" & score_status %in% c("scored_preferred_coverage", "scored_low_coverage_exploratory")])
  }
  dt[0]
}

scenario_rows <- list()
for (sid in scenario_defs$scenario_id) {
  dt <- scenario_filter(perf, sid)
  if (!nrow(dt)) {
    scenario_rows[[length(scenario_rows) + 1L]] <- data.table(
      scenario_id = sid,
      classifier_id = character(),
      classifier_score_name = character()
    )
    next
  }
  sm <- dt[, .(
    n_datasets = .N,
    n_total_samples = sum(n, na.rm = TRUE),
    n_total_pcr = sum(n_pcr, na.rm = TRUE),
    median_auc = safe_median(discovery_aligned_auc),
    min_auc = safe_min(discovery_aligned_auc),
    max_auc = safe_max(discovery_aligned_auc),
    iqr_auc = ifelse(sum(is.finite(discovery_aligned_auc)) >= 2,
                     diff(as.numeric(quantile(discovery_aligned_auc, c(0.25, 0.75), na.rm = TRUE))), NA_real_),
    median_pr_auc = safe_median(average_precision_pr_auc),
    median_pcr_prevalence = safe_median(pcr_prevalence),
    n_auc_ge_0_60 = sum(discovery_aligned_auc >= 0.60, na.rm = TRUE),
    n_auc_ge_0_65 = sum(discovery_aligned_auc >= 0.65, na.rm = TRUE),
    datasets = paste(matrix_label, collapse = ";")
  ), by = .(classifier_id, classifier_score_name)]
  sm[, scenario_id := sid]
  sm <- merge(sm, scenario_defs[, .(scenario_id, scenario_label, scenario_definition, primary_for_claim)], by = "scenario_id", all.x = TRUE)
  sm[, scenario_rank := frank(-median_auc, ties.method = "min")]
  sm[, evaluable := ifelse(n_datasets >= 2 & is.finite(median_auc), "yes", "limited")]
  scenario_rows[[length(scenario_rows) + 1L]] <- sm
}
sensitivity <- rbindlist(scenario_rows, fill = TRUE)

primary_rank <- sensitivity[scenario_id == "primary_core_preferred", .(classifier_id, primary_rank = scenario_rank, primary_median_auc = median_auc)]
sensitivity <- merge(sensitivity, primary_rank, by = "classifier_id", all.x = TRUE)
sensitivity[, rank_shift_vs_primary := scenario_rank - primary_rank]
setorder(sensitivity, scenario_id, scenario_rank, classifier_id)
fwrite(sensitivity, file.path(step6_dir, "sensitivity_auc_summary_v1.csv"))

primary_result_lock <- sensitivity[scenario_id == "primary_core_preferred", .(
  classifier_id, classifier_score_name,
  locked_primary_analysis_set = "external_core_validation",
  locked_score_status = "scored_preferred_coverage",
  n_external_core_datasets = n_datasets,
  median_external_auc = median_auc,
  min_external_auc = min_auc,
  max_external_auc = max_auc,
  n_auc_ge_0_60,
  n_auc_ge_0_65,
  median_pr_auc,
  median_pcr_prevalence,
  primary_rank = scenario_rank,
  primary_interpretation = fifelse(
    n_datasets >= 5 & median_auc >= 0.65, "primary_strong_candidate",
    fifelse(n_datasets >= 5 & median_auc >= 0.60, "primary_moderate_candidate",
            fifelse(n_datasets >= 3 & median_auc >= 0.55, "primary_weak_or_context_dependent", "primary_poor_or_insufficient"))
  ),
  locked_conclusion_allowed = fifelse(n_datasets >= 5, "yes", "limited_evidence")
)]
setorder(primary_result_lock, primary_rank)
fwrite(primary_result_lock, file.path(step6_dir, "primary_result_lock_v1.csv"))

final_rank <- copy(primary_result_lock)
stability <- sensitivity[scenario_id != "primary_core_preferred" & is.finite(median_auc), .(
  n_sensitivity_scenarios_evaluable = .N,
  median_sensitivity_auc = safe_median(median_auc),
  min_sensitivity_auc = safe_min(median_auc),
  max_sensitivity_auc = safe_max(median_auc),
  max_abs_rank_shift = ifelse(all(is.na(rank_shift_vs_primary)), NA_real_, max(abs(rank_shift_vs_primary), na.rm = TRUE)),
  n_scenarios_auc_ge_0_60 = sum(median_auc >= 0.60, na.rm = TRUE),
  n_scenarios_auc_ge_0_65 = sum(median_auc >= 0.65, na.rm = TRUE)
), by = classifier_id]
final_rank <- merge(final_rank, stability, by = "classifier_id", all.x = TRUE)
final_rank[, robustness_tier := fifelse(
  primary_interpretation == "primary_strong_candidate" & n_scenarios_auc_ge_0_60 >= 4, "robust",
  fifelse(primary_interpretation %in% c("primary_strong_candidate", "primary_moderate_candidate") & n_scenarios_auc_ge_0_60 >= 3, "moderately_robust",
          fifelse(primary_interpretation %in% c("primary_weak_or_context_dependent") | median_external_auc < 0.60, "fragile_or_context_dependent", "uncertain"))
)]
setorder(final_rank, primary_rank)
fwrite(final_rank, file.path(step6_dir, "final_classifier_ranking_locked_v1.csv"))

cal_primary <- cal[calibration_role == "external_core_validation" & score_status == "scored_preferred_coverage"]
cal_summary <- cal_primary[, .(
  n_external_core_datasets = .N,
  median_brier = safe_median(brier_score),
  min_brier = safe_min(brier_score),
  max_brier = safe_max(brier_score),
  median_abs_calibration_intercept = safe_median(abs(calibration_intercept_offset)),
  median_calibration_slope = safe_median(calibration_slope),
  min_calibration_slope = safe_min(calibration_slope),
  max_calibration_slope = safe_max(calibration_slope),
  median_observed_expected_ratio = safe_median(observed_expected_ratio),
  n_abs_intercept_gt_0_5 = sum(abs(calibration_intercept_offset) > 0.5, na.rm = TRUE),
  n_slope_outside_0_5_1_5 = sum(calibration_slope < 0.5 | calibration_slope > 1.5, na.rm = TRUE),
  n_oe_outside_0_8_1_25 = sum(observed_expected_ratio < 0.8 | observed_expected_ratio > 1.25, na.rm = TRUE)
), by = .(classifier_id, classifier_score_name)]
cal_summary[, calibration_failure_class := fifelse(
  n_abs_intercept_gt_0_5 >= 3 | n_slope_outside_0_5_1_5 >= 3 | n_oe_outside_0_8_1_25 >= 3,
  "frequent_calibration_failure",
  fifelse(n_abs_intercept_gt_0_5 + n_slope_outside_0_5_1_5 + n_oe_outside_0_8_1_25 > 0,
          "some_calibration_instability", "no_major_calibration_flag")
)]
setorder(cal_summary, classifier_id)
fwrite(cal_summary, file.path(step6_dir, "calibration_failure_summary_v1.csv"))

subtype_summary <- data.table()
if (nrow(subtype)) {
  subtype[, subtype_story_group := fifelse(grepl("TNBC|ER_negative|HR_negative", subtype_group, ignore.case = TRUE),
                                           "TNBC_or_ER_negative",
                                           fifelse(grepl("HER2_positive", subtype_group, ignore.case = TRUE),
                                                   "HER2_positive",
                                                   fifelse(grepl("HR_positive", subtype_group, ignore.case = TRUE),
                                                           "HR_positive_HER2_negative", "other")))]
  subtype_summary <- subtype[subtype_story_group %in% c("TNBC_or_ER_negative", "HER2_positive"), .(
    n_dataset_subgroups = .N,
    n_total_samples = sum(n, na.rm = TRUE),
    n_total_pcr = sum(n_pcr, na.rm = TRUE),
    median_auc = safe_median(discovery_aligned_auc),
    min_auc = safe_min(discovery_aligned_auc),
    max_auc = safe_max(discovery_aligned_auc),
    median_pr_auc = safe_median(average_precision_pr_auc),
    n_auc_ge_0_60 = sum(discovery_aligned_auc >= 0.60, na.rm = TRUE),
    dataset_subgroups = paste(paste(matrix_label, subtype_group, sep = ":"), collapse = ";")
  ), by = .(subtype_story_group, classifier_id, classifier_score_name)]
  subtype_summary[, evaluable := ifelse(n_dataset_subgroups >= 2, "yes", "limited")]
  subtype_summary[, subtype_display := factor(
    fifelse(subtype_story_group == "TNBC_or_ER_negative", "TNBC/ER-negative", "HER2-positive"),
    levels = c("TNBC/ER-negative", "HER2-positive")
  )]
  setorder(subtype_summary, subtype_display, classifier_id)
}
fwrite(subtype_summary, file.path(step6_dir, "subtype_sensitivity_summary_v1.csv"))

claim_rows <- list(
  data.table(
    claim_id = "C1",
    manuscript_claim = "Published transcriptomic scores show modest but non-uniform external transferability for NAC pCR/RD prediction.",
    evidence_file = "primary_result_lock_v1.csv; sensitivity_auc_summary_v1.csv",
    evidence_summary = paste0("Primary median external AUC range: ",
                              sprintf("%.3f", safe_min(primary_result_lock$median_external_auc)), " to ",
                              sprintf("%.3f", safe_max(primary_result_lock$median_external_auc)), "."),
    claim_strength = "supported_with_moderate_wording",
    limitation = "Retrospective public-cohort benchmark; no prospective clinical validation."
  ),
  data.table(
    claim_id = "C2",
    manuscript_claim = "Immune-related and intrinsic-subtype/proliferation-associated scores rank highest in the locked external core analysis.",
    evidence_file = "final_classifier_ranking_locked_v1.csv",
    evidence_summary = paste(final_rank[primary_rank <= 3, paste0(primary_rank, ":", classifier_id, " AUC=", sprintf("%.3f", median_external_auc))], collapse = "; "),
    claim_strength = "supported",
    limitation = "AUC is still moderate; do not claim clinical utility."
  ),
  data.table(
    claim_id = "C3",
    manuscript_claim = "Calibration is less stable than discrimination across external cohorts.",
    evidence_file = "calibration_failure_summary_v1.csv",
    evidence_summary = paste(cal_summary[, paste0(classifier_id, ":", calibration_failure_class)], collapse = "; "),
    claim_strength = "supported_if_reported_as_failure_mode",
    limitation = "Calibration used GSE25066-trained single-score logistic models; no external recalibration was performed."
  ),
  data.table(
    claim_id = "C4",
    manuscript_claim = "Limited-panel or low-gene-coverage settings should be treated as exploratory rather than primary evidence.",
    evidence_file = "sensitivity_auc_summary_v1.csv; classifier_failure_mode_audit_v1.csv",
    evidence_summary = "Low-coverage and limited/panel scenarios are separated from the locked primary analysis.",
    claim_strength = "supported_by_design",
    limitation = "This is a predefined analytic boundary, not a biological result."
  ),
  data.table(
    claim_id = "C5",
    manuscript_claim = "The study should be framed as external validation and failure-mode analysis, not new signature development.",
    evidence_file = "validation_analysis_plan_v1.md; primary_result_lock_v1.csv",
    evidence_summary = "No new feature selection, cutoff tuning, or validation-cohort recalibration was used for primary conclusions.",
    claim_strength = "supported_by_protocol",
    limitation = "Approximate scores for commercial assays are not equivalent to proprietary clinical test outputs."
  )
)
claim_evidence <- rbindlist(claim_rows, fill = TRUE)
fwrite(claim_evidence, file.path(step6_dir, "claim_evidence_table_v1.csv"))

plan_lines <- c(
  "# Step 6 sensitivity analysis and primary result lock v1",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Lock the main external-validation conclusion and test whether it changes under prespecified sensitivity scenarios.",
  "",
  "## Locked primary analysis",
  "",
  "- `performance_role == external_core_validation`",
  "- `score_status == scored_preferred_coverage`",
  "- Discovery base GSE25066 is excluded from external validation ranking.",
  "- GSE194040 is reported as stress-test only.",
  "- Low-coverage exploratory scores and limited-panel datasets are not mixed into the primary result.",
  "",
  "## Sensitivity scenarios",
  "",
  paste0("- ", scenario_defs$scenario_label, ": ", scenario_defs$scenario_definition),
  "",
  "## Interpretation boundary",
  "",
  "Step 6 can support a paper claim about transferability, instability, and failure modes. It cannot support a claim of clinical utility, prospective predictive validity, or superiority over commercial assays."
)
writeLines(plan_lines, file.path(step6_dir, "sensitivity_analysis_plan_v1.md"), useBytes = TRUE)

figure_rows <- list()
register_fig <- function(filename, status, notes) {
  figure_rows[[length(figure_rows) + 1L]] <<- data.table(
    figure_file = filename,
    exists = file.exists(file.path(fig_dir, filename)),
    status = status,
    notes = notes
  )
}

if (nrow(sensitivity)) {
  plot_dt <- sensitivity[is.finite(median_auc)]
  plot_dt <- add_classifier_display(plot_dt)
  plot_dt[, scenario_label := factor(scenario_label, levels = scenario_defs$scenario_label)]
  p <- ggplot(plot_dt, aes(x = scenario_label, y = classifier_axis, fill = median_auc)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", median_auc)), size = 2.6) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0.60,
                         limits = c(0.45, 0.75), na.value = "grey90", name = "Median AUC") +
    labs(title = "Sensitivity analysis of classifier ranking",
         subtitle = "Median AUC across prespecified scenario-specific datasets",
         x = "Sensitivity scenario", y = "Classifier") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  fn <- "fig_sensitivity_rank_shift_v1.pdf"
  ggsave(file.path(fig_dir, fn), p, width = 12, height = 6, device = "pdf")
  register_fig(fn, "generated", "Generated from sensitivity_auc_summary_v1.csv")
} else {
  register_fig("fig_sensitivity_rank_shift_v1.pdf", "not_generated", "No sensitivity rows")
}

if (nrow(subtype_summary)) {
  plot_dt <- subtype_summary[is.finite(median_auc)]
  plot_dt <- add_classifier_display(plot_dt)
  plot_dt[, subtype_display := factor(as.character(subtype_display), levels = c("TNBC/ER-negative", "HER2-positive"))]
  p <- ggplot(plot_dt, aes(x = subtype_display, y = classifier_axis, fill = median_auc)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", median_auc)), size = 2.7) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0.60,
                         limits = c(0.45, 0.8), na.value = "grey90", name = "Median AUC") +
    labs(title = "Subtype sensitivity summary",
         subtitle = "Median AUC by broad subtype story group",
         x = "Subtype story group", y = "Classifier") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())
  fn <- "fig_subtype_auc_heatmap_v1.pdf"
  ggsave(file.path(fig_dir, fn), p, width = 9, height = 6, device = "pdf")
  register_fig(fn, "generated", "Generated from subtype_sensitivity_summary_v1.csv")
} else {
  register_fig("fig_subtype_auc_heatmap_v1.pdf", "not_generated", "No subtype sensitivity rows")
}

if (nrow(cal_summary)) {
  plot_dt <- rbindlist(list(
    cal_summary[, .(
      classifier_id,
      classifier_score_name,
      calibration_metric = "Median |intercept| <= 0.5",
      metric_pass = is.finite(median_abs_calibration_intercept) & median_abs_calibration_intercept <= 0.5
    )],
    cal_summary[, .(
      classifier_id,
      classifier_score_name,
      calibration_metric = "Median slope 0.5-1.5",
      metric_pass = is.finite(median_calibration_slope) &
        median_calibration_slope >= 0.5 & median_calibration_slope <= 1.5
    )],
    cal_summary[, .(
      classifier_id,
      classifier_score_name,
      calibration_metric = "Median O/E 0.8-1.25",
      metric_pass = is.finite(median_observed_expected_ratio) &
        median_observed_expected_ratio >= 0.8 & median_observed_expected_ratio <= 1.25
    )]
  ), fill = TRUE)
  plot_dt <- add_classifier_display(plot_dt)
  plot_dt[, calibration_metric := factor(
    calibration_metric,
    levels = c("Median |intercept| <= 0.5", "Median slope 0.5-1.5", "Median O/E 0.8-1.25")
  )]
  plot_dt[, calibration_flag := fifelse(
    metric_pass,
    "Within acceptable range",
    "Failed"
  )]
  plot_dt[, calibration_flag := factor(
    calibration_flag,
    levels = c("Within acceptable range", "Failed")
  )]
  p <- ggplot(plot_dt, aes(x = calibration_metric, y = classifier_axis, fill = calibration_flag)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = ifelse(metric_pass, "OK", "Fail")), size = 2.7) +
    scale_fill_manual(
      values = c("Within acceptable range" = "#DCECCB", "Failed" = "#D24B40"),
      drop = FALSE,
      name = "Flag"
    ) +
    labs(
      title = "Calibration failure flags",
      subtitle = "External core validation only; discrete flags avoid mixing intercept, slope, and O/E value scales",
      x = "Calibration acceptability rule",
      y = "Classifier"
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
  fn <- "fig_calibration_failure_summary_v1.pdf"
  ggsave(file.path(fig_dir, fn), p, width = 11, height = 6.5, device = "pdf")
  register_fig(fn, "generated", "Generated from calibration_failure_summary_v1.csv")
} else {
  register_fig("fig_calibration_failure_summary_v1.pdf", "not_generated", "No calibration summary rows")
}

fig_manifest <- rbindlist(figure_rows, fill = TRUE)
fwrite(fig_manifest, file.path(step6_dir, "figure_manifest_v1.csv"))

notes_lines <- c(
  "# Step 6 execution notes v1",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Output summary",
  "",
  paste0("- Sensitivity rows: ", nrow(sensitivity)),
  paste0("- Primary locked classifier rows: ", nrow(primary_result_lock)),
  paste0("- Final locked ranking rows: ", nrow(final_rank)),
  paste0("- Calibration summary rows: ", nrow(cal_summary)),
  paste0("- Subtype sensitivity rows: ", nrow(subtype_summary)),
  paste0("- Claim-evidence rows: ", nrow(claim_evidence)),
  "",
  "## Primary locked ranking",
  "",
  paste(capture.output(print(final_rank[, .(primary_rank, classifier_id, median_external_auc, robustness_tier, primary_interpretation)])), collapse = "\n"),
  "",
  "## Figure generation status",
  "",
  paste(capture.output(print(fig_manifest)), collapse = "\n"),
  "",
  "## Integrity statement",
  "",
  "All Step 6 tables and figures were generated by build_sensitivity_step6_v1.R from Step 5 result tables. No result value or figure was manually filled."
)
writeLines(notes_lines, file.path(step6_dir, "step6_sensitivity_notes_v1.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(step6_dir, "R_session_info_step6_v1.txt"))

cat("\nOutput rows:\n")
cat("sensitivity_auc_summary_v1.csv:", nrow(sensitivity), "\n")
cat("primary_result_lock_v1.csv:", nrow(primary_result_lock), "\n")
cat("final_classifier_ranking_locked_v1.csv:", nrow(final_rank), "\n")
cat("calibration_failure_summary_v1.csv:", nrow(cal_summary), "\n")
cat("subtype_sensitivity_summary_v1.csv:", nrow(subtype_summary), "\n")
cat("claim_evidence_table_v1.csv:", nrow(claim_evidence), "\n")
cat("figure_manifest_v1.csv:", nrow(fig_manifest), "\n")

if (any(fig_manifest$status != "generated")) {
  cat("\nSome Step 6 figures were not generated. Inspect figure_manifest_v1.csv.\n")
} else {
  cat("\nAll planned Step 6 PDF figures were generated by code.\n")
}

cat("\nDone.\n")
