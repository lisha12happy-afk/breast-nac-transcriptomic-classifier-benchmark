options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_root <- file.path(project_root, "manifests", "final_article_tables")
desktop_package <- file.path(out_root, "desktop_table_package")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(desktop_package, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_root, "R_raw_run_log_final_article_tables_v1.txt")
session_file <- file.path(out_root, "R_session_info_final_article_tables_v1.txt")
sink(log_file, split = TRUE)

cat("Build final manuscript-ready tables, excluding already-finished Table 1\n")
cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Project root:", project_root, "\n")
cat("Output root:", out_root, "\n")

step1 <- file.path(project_root, "manifests", "step 1")
step2 <- file.path(project_root, "manifests", "step 2")
step3 <- file.path(project_root, "manifests", "step 3")
step5 <- file.path(project_root, "manifests", "step 5")
step6 <- file.path(project_root, "manifests", "step 6")

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.csv(path, na.strings = c("", "NA"), check.names = FALSE)
}

write_one_table <- function(tbl, label, folder_name, file_stub) {
  local_dir <- file.path(out_root, folder_name)
  desktop_dir <- file.path(desktop_package, folder_name)
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(desktop_dir, recursive = TRUE, showWarnings = FALSE)

  csv_path <- file.path(local_dir, paste0(file_stub, ".csv"))
  csv_desktop_path <- file.path(desktop_dir, paste0(file_stub, ".csv"))
  xlsx_path <- file.path(local_dir, paste0(file_stub, ".xlsx"))
  xlsx_desktop_path <- file.path(desktop_dir, paste0(file_stub, ".xlsx"))

  write.csv(tbl, csv_path, row.names = FALSE, na = "")
  write.csv(tbl, csv_desktop_path, row.names = FALSE, na = "")

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, label)
    openxlsx::writeData(wb, label, label, startRow = 1, startCol = 1)
    openxlsx::writeData(wb, label, tbl, startRow = 3, startCol = 1)
    header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9EAF7", border = "TopBottomLeftRight", wrapText = TRUE, valign = "center")
    body_style <- openxlsx::createStyle(border = "TopBottomLeftRight", wrapText = TRUE, valign = "top")
    title_style <- openxlsx::createStyle(textDecoration = "bold", fontSize = 12)
    openxlsx::addStyle(wb, label, title_style, rows = 1, cols = 1)
    openxlsx::addStyle(wb, label, header_style, rows = 3, cols = seq_len(ncol(tbl)), gridExpand = TRUE)
    openxlsx::addStyle(wb, label, body_style, rows = 4:(nrow(tbl) + 3), cols = seq_len(ncol(tbl)), gridExpand = TRUE)
    openxlsx::setColWidths(wb, label, cols = seq_len(ncol(tbl)), widths = "auto")
    openxlsx::freezePane(wb, label, firstActiveRow = 4)
    openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
    openxlsx::saveWorkbook(wb, xlsx_desktop_path, overwrite = TRUE)
  } else if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(setNames(list(tbl), label), xlsx_path)
    writexl::write_xlsx(setNames(list(tbl), label), xlsx_desktop_path)
  } else {
    cat("No XLSX writer installed; CSV only for ", file_stub, "\n", sep = "")
  }

  cat("Wrote ", label, ": ", nrow(tbl), " rows x ", ncol(tbl), " cols\n", sep = "")
  data.frame(
    table_label = label,
    rows = nrow(tbl),
    columns = ncol(tbl),
    local_csv = csv_path,
    local_xlsx = if (file.exists(xlsx_path)) xlsx_path else NA_character_,
    desktop_package_folder = desktop_dir,
    stringsAsFactors = FALSE
  )
}

write_workbook <- function(sheet_list, label, folder_name, file_stub) {
  local_dir <- file.path(out_root, folder_name)
  desktop_dir <- file.path(desktop_package, folder_name)
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(desktop_dir, recursive = TRUE, showWarnings = FALSE)

  # Also write one CSV per sheet so each section can be copied independently.
  for (nm in names(sheet_list)) {
    write.csv(sheet_list[[nm]], file.path(local_dir, paste0(file_stub, "_", nm, ".csv")), row.names = FALSE, na = "")
    write.csv(sheet_list[[nm]], file.path(desktop_dir, paste0(file_stub, "_", nm, ".csv")), row.names = FALSE, na = "")
  }

  xlsx_path <- file.path(local_dir, paste0(file_stub, ".xlsx"))
  xlsx_desktop_path <- file.path(desktop_dir, paste0(file_stub, ".xlsx"))
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    for (nm in names(sheet_list)) {
      tbl <- sheet_list[[nm]]
      openxlsx::addWorksheet(wb, nm)
      openxlsx::writeData(wb, nm, paste(label, nm), startRow = 1, startCol = 1)
      openxlsx::writeData(wb, nm, tbl, startRow = 3, startCol = 1)
      header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9EAF7", border = "TopBottomLeftRight", wrapText = TRUE, valign = "center")
      body_style <- openxlsx::createStyle(border = "TopBottomLeftRight", wrapText = TRUE, valign = "top")
      openxlsx::addStyle(wb, nm, header_style, rows = 3, cols = seq_len(ncol(tbl)), gridExpand = TRUE)
      openxlsx::addStyle(wb, nm, body_style, rows = 4:(nrow(tbl) + 3), cols = seq_len(ncol(tbl)), gridExpand = TRUE)
      openxlsx::setColWidths(wb, nm, cols = seq_len(ncol(tbl)), widths = "auto")
      openxlsx::freezePane(wb, nm, firstActiveRow = 4)
    }
    openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
    openxlsx::saveWorkbook(wb, xlsx_desktop_path, overwrite = TRUE)
  } else if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheet_list, xlsx_path)
    writexl::write_xlsx(sheet_list, xlsx_desktop_path)
  } else {
    cat("No XLSX writer installed; CSV only for ", file_stub, "\n", sep = "")
  }
  cat("Wrote ", label, ": ", length(sheet_list), " sheets\n", sep = "")
  data.frame(
    table_label = label,
    rows = paste(vapply(sheet_list, nrow, integer(1)), collapse = "; "),
    columns = paste(vapply(sheet_list, ncol, integer(1)), collapse = "; "),
    local_csv = paste(file.path(local_dir, paste0(file_stub, "_", names(sheet_list), ".csv")), collapse = "; "),
    local_xlsx = if (file.exists(xlsx_path)) xlsx_path else NA_character_,
    desktop_package_folder = desktop_dir,
    stringsAsFactors = FALSE
  )
}

fmt3 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", formatC(x, digits = 3, format = "f"))
}
fmt1 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", formatC(x, digits = 1, format = "f"))
}
clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- gsub(";", "; ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}
yn <- function(x) {
  x <- as.character(x)
  ifelse(tolower(x) %in% c("yes", "true", "1"), "Yes",
         ifelse(tolower(x) %in% c("no", "false", "0"), "No", clean_text(x)))
}

manifest <- read_csv(file.path(step1, "analysis_ready_manifest_v1.csv"))
expr <- read_csv(file.path(step2, "expression_matrix_audit_v1.csv"))
dict <- read_csv(file.path(step3, "classifier_dictionary_v1.csv"))
risk <- read_csv(file.path(step3, "classifier_reproducibility_risk_v1.csv"))
coverage <- read_csv(file.path(step3, "classifier_gene_coverage_v1.csv"))
performance <- read_csv(file.path(step5, "classifier_performance_by_dataset_v1.csv"))
calibration <- read_csv(file.path(step5, "classifier_calibration_by_dataset_v1.csv"))
failure <- read_csv(file.path(step5, "classifier_failure_mode_audit_v1.csv"))
rank <- read_csv(file.path(step6, "final_classifier_ranking_locked_v1.csv"))
cal_fail <- read_csv(file.path(step6, "calibration_failure_summary_v1.csv"))
sensitivity <- read_csv(file.path(step6, "sensitivity_auc_summary_v1.csv"))
subtype <- read_csv(file.path(step6, "subtype_sensitivity_summary_v1.csv"))

cat("Loaded source tables.\n")

# ---- Main Table 2: classifier reconstruction and reproducibility ----
dict_small <- dict[, intersect(c("classifier_id", "original_endpoint_or_use"), names(dict)), drop = FALSE]
t2 <- merge(risk, dict_small, by = "classifier_id", all.x = TRUE, sort = FALSE)
rank_small <- rank[, intersect(c("classifier_id", "primary_rank", "robustness_tier", "locked_conclusion_allowed"), names(rank)), drop = FALSE]
t2 <- merge(t2, rank_small, by = "classifier_id", all.x = TRUE, sort = FALSE)
t2$analysis_decision <- ifelse(
  !is.na(t2$primary_rank),
  paste0("Scored in primary benchmark; locked rank ", t2$primary_rank, "."),
  paste0("Not scored in primary benchmark; ", clean_text(t2$recommended_step4_action), ".")
)
t2$primary_reference <- ifelse(
  !is.na(t2$primary_reference_doi) & nzchar(t2$primary_reference_doi),
  paste0("PMID ", t2$primary_reference_pmid, "; DOI ", t2$primary_reference_doi),
  paste0("PMID ", t2$primary_reference_pmid)
)
t2 <- t2[order(is.na(t2$primary_rank), t2$primary_rank, t2$classifier_id), ]
table2 <- data.frame(
  `Classifier` = t2$classifier_name,
  `Original endpoint or use` = t2$original_endpoint_or_use,
  `pCR-specific` = yn(t2$pcr_specific),
  `Gene-list status` = clean_text(t2$gene_list_status),
  `Formula/score reconstruction status` = clean_text(t2$formula_reconstruction_status),
  `Reproducibility risk` = clean_text(t2$reproducibility_risk),
  `Final analysis decision` = t2$analysis_decision,
  `Primary reference` = t2$primary_reference,
  check.names = FALSE
)

# ---- Main Table 3: locked external validation ranking ----
rank <- rank[order(rank$primary_rank), ]
interpretation_map <- c(
  primary_strong_candidate = "Strongest transferable signal in this benchmark",
  primary_moderate_candidate = "Moderate transferable signal",
  primary_weak_or_context_dependent = "Weak or context-dependent transferability"
)
table3 <- data.frame(
  `Rank` = rank$primary_rank,
  `Classifier` = rank$classifier_score_name,
  `External core datasets, n` = rank$n_external_core_datasets,
  `Median external AUC (range)` = paste0(fmt3(rank$median_external_auc), " (", fmt3(rank$min_external_auc), "-", fmt3(rank$max_external_auc), ")"),
  `Median PR-AUC` = fmt3(rank$median_pr_auc),
  `Datasets with AUC >=0.60` = paste0(rank$n_auc_ge_0_60, "/", rank$n_external_core_datasets),
  `Datasets with AUC >=0.65` = paste0(rank$n_auc_ge_0_65, "/", rank$n_external_core_datasets),
  `Median sensitivity AUC (range)` = paste0(fmt3(rank$median_sensitivity_auc), " (", fmt3(rank$min_sensitivity_auc), "-", fmt3(rank$max_sensitivity_auc), ")"),
  `Maximum absolute rank shift` = rank$max_abs_rank_shift,
  `Robustness tier` = clean_text(rank$robustness_tier),
  `Interpretation` = ifelse(rank$primary_interpretation %in% names(interpretation_map), interpretation_map[rank$primary_interpretation], clean_text(rank$primary_interpretation)),
  check.names = FALSE
)

# ---- Main Table 4: calibration and failure-mode summary ----
t4 <- merge(rank[, c("classifier_id", "classifier_score_name", "primary_rank", "median_external_auc", "robustness_tier")],
            cal_fail, by = "classifier_id", all.x = TRUE, sort = FALSE)
t4 <- t4[order(t4$primary_rank), ]
calibration_interpretation <- function(robustness, fail_class) {
  if (grepl("fragile", robustness)) {
    return("Discrimination and calibration were context-dependent; interpret as exploratory.")
  }
  if (grepl("frequent", fail_class)) {
    return("Discrimination was evaluated externally, but calibration was unstable; recalibration would be required before clinical use.")
  }
  "Calibration should still be checked before use outside the benchmark."
}
table4 <- data.frame(
  `Rank` = t4$primary_rank,
  `Classifier` = t4$classifier_score_name.x,
  `Median external AUC` = fmt3(t4$median_external_auc),
  `Median Brier score` = fmt3(t4$median_brier),
  `Median |calibration intercept|` = fmt3(t4$median_abs_calibration_intercept),
  `Median calibration slope` = fmt3(t4$median_calibration_slope),
  `Datasets with slope outside 0.5-1.5` = paste0(t4$n_slope_outside_0_5_1_5, "/", t4$n_external_core_datasets),
  `Datasets with O/E outside 0.8-1.25` = paste0(t4$n_oe_outside_0_8_1_25, "/", t4$n_external_core_datasets),
  `Calibration failure class` = clean_text(t4$calibration_failure_class),
  `Robustness tier` = clean_text(t4$robustness_tier),
  `Manuscript-use interpretation` = mapply(calibration_interpretation, t4$robustness_tier, t4$calibration_failure_class, USE.NAMES = FALSE),
  check.names = FALSE
)

# ---- Supplementary Table S1: cohort harmonization and expression availability ----
manifest$key <- paste(manifest$dataset_accession, manifest$platform, sep = "||")
expr$key <- paste(expr$dataset_accession, expr$platform, sep = "||")
manifest$endpoint_binary_clean <- suppressWarnings(as.integer(manifest$endpoint_binary_pcr1_rd0))
split_manifest <- split(manifest, manifest$key, drop = TRUE)
endpoint_by_platform <- do.call(rbind, lapply(names(split_manifest), function(k) {
  dat <- split_manifest[[k]]
  valid <- dat$endpoint_binary_clean %in% c(0L, 1L)
  data.frame(
    key = k,
    endpoint_mapped_n = sum(valid, na.rm = TRUE),
    pCR_n = sum(dat$endpoint_binary_clean == 1L, na.rm = TRUE),
    RD_n = sum(dat$endpoint_binary_clean == 0L, na.rm = TRUE),
    include_primary_values = paste(sort(unique(na.omit(dat$include_primary_analysis))), collapse = "; "),
    overlap_flags = paste(sort(unique(na.omit(dat$overlap_flag))), collapse = "; "),
    stringsAsFactors = FALSE
  )
}))
idx <- match(expr$key, endpoint_by_platform$key)
source_citation <- c(
  GSE25066 = "GEO accession GSE25066; GEO-linked citation(s): PMID 21558518, DOI 10.1001/jama.2011.593; PMID 24337596, DOI 10.1007/s10549-013-2763-z; PMID 36293478, DOI 10.3390/ijms232012625.",
  GSE20271 = "GEO accession GSE20271; GEO-linked citation(s): PMID 20829329, DOI 10.1158/1078-0432.CCR-10-1265; PMID 23185353, DOI 10.1371/journal.pone.0049529.",
  GSE32646 = "GEO accession GSE32646; GEO-linked citation(s): PMID 22320227, DOI 10.1111/j.1349-7006.2012.02231.x.",
  GSE41998 = "GEO accession GSE41998; GEO-linked citation(s): PMID 23340299, DOI 10.1158/1078-0432.CCR-12-1359.",
  GSE50948 = "GEO accession GSE50948; GEO-linked citation(s): PMID 24443618, DOI 10.1158/1078-0432.CCR-13-0239.",
  GSE66305 = "GEO accession GSE66305; GEO-linked citation(s): PMID 26245675, DOI 10.1634/theoncologist.2015-0138.",
  GSE163882 = "GEO accession GSE163882; GEO-linked citation(s): PMID 35935976, DOI 10.3389/fimmu.2022.948601.",
  GSE194040 = "GEO accession GSE194040; GEO-linked citation(s): PMID 35623341, DOI 10.1016/j.ccell.2022.05.005; PMID 40905691, DOI 10.1158/1078-0432.CCR-25-0553.",
  GSE106977 = "GEO accession GSE106977; GEO-linked citation(s): PMID 29899867, DOI 10.18632/oncotarget.25413.",
  GSE109710 = "GEO accession GSE109710; GEO-linked citation(s): PMID 29788230, DOI 10.1093/jnci/djy076.",
  GSE130786 = "GEO accession GSE130786; GEO-linked citation(s): PMID 33951110, DOI 10.1371/journal.pone.0251163; preprint DOI 10.1101/2021.01.26.21250250."
)
source_citation_vec <- unname(source_citation[as.character(expr$dataset_accession)])
source_citation_vec[is.na(source_citation_vec) | !nzchar(source_citation_vec)] <- paste0("GEO accession ", expr$dataset_accession, "; source publication not mapped in script.")
s1 <- data.frame(
  `Dataset accession` = expr$dataset_accession,
  `Platform` = expr$platform,
  `Source publication / dataset citation` = source_citation_vec,
  `Matrix label` = expr$matrix_label,
  `Expression source type` = clean_text(expr$expression_source_type),
  `Matrix samples, n` = expr$n_samples_in_matrix,
  `Matched manifest samples, n` = expr$n_samples_matched_to_manifest,
  `Endpoint-mapped samples, n` = endpoint_by_platform$endpoint_mapped_n[idx],
  `pCR, n` = endpoint_by_platform$pCR_n[idx],
  `RD, n` = endpoint_by_platform$RD_n[idx],
  `Unique genes after annotation, n` = expr$n_unique_genes_after_annotation,
  `Expression scale guess` = clean_text(expr$expression_scale_guess),
  `Usable for full benchmark` = yn(expr$usable_for_full_benchmark),
  `Usable for subtype-only validation` = yn(expr$usable_for_subtype_only_validation),
  `Planned role` = clean_text(expr$planned_role),
  `Story scope` = clean_text(expr$story_scope),
  `Raw data availability` = clean_text(expr$raw_data_availability),
  `Overlap or limitation note` = ifelse(!is.na(endpoint_by_platform$overlap_flags[idx]) & nzchar(endpoint_by_platform$overlap_flags[idx]), endpoint_by_platform$overlap_flags[idx], expr$notes),
  check.names = FALSE
)

# ---- Supplementary Table S2: classifier reconstruction and gene coverage ----
s2 <- coverage[, c(
  "classifier_id", "classifier_name", "dataset_accession", "platform", "analysis_set",
  "n_harmonized_genes", "n_genes_present", "coverage_pct_of_harmonized",
  "coverage_grade", "missing_genes", "usable_for_full_benchmark",
  "usable_for_subtype_only_validation", "expression_scale_guess"
)]
names(s2) <- c(
  "Classifier ID", "Classifier", "Dataset accession", "Platform", "Analysis set",
  "Expected harmonized genes, n", "Genes present, n", "Coverage, %",
  "Coverage grade", "Missing genes", "Usable for full benchmark",
  "Usable for subtype-only validation", "Expression scale guess"
)
s2$`Coverage, %` <- fmt1(s2$`Coverage, %`)
s2$`Analysis set` <- clean_text(s2$`Analysis set`)
s2$`Coverage grade` <- clean_text(s2$`Coverage grade`)
s2$`Usable for full benchmark` <- yn(s2$`Usable for full benchmark`)
s2$`Usable for subtype-only validation` <- yn(s2$`Usable for subtype-only validation`)
s2$`Expression scale guess` <- clean_text(s2$`Expression scale guess`)
s2 <- s2[order(s2$`Dataset accession`, s2$Platform, s2$`Classifier ID`), ]

# ---- Supplementary Table S3: full external validation performance ----
s3 <- performance[, c(
  "dataset_accession", "platform", "analysis_set", "classifier_id", "classifier_score_name",
  "score_status", "n", "n_pcr", "n_rd", "pcr_prevalence",
  "discovery_aligned_auc", "auc_ci_low_bootstrap", "auc_ci_high_bootstrap",
  "average_precision_pr_auc", "pr_auc_baseline_pcr_prevalence",
  "wilcoxon_p_value_aligned_score", "wilcoxon_q_value_aligned_score",
  "performance_role"
)]
names(s3) <- c(
  "Dataset accession", "Platform", "Analysis set", "Classifier ID", "Classifier score",
  "Score status", "Samples, n", "pCR, n", "RD, n", "pCR prevalence",
  "AUC", "AUC 95% CI low", "AUC 95% CI high",
  "PR-AUC", "PR-AUC baseline", "Wilcoxon P", "Wilcoxon q", "Performance role"
)
s3$`pCR prevalence` <- fmt3(s3$`pCR prevalence`)
s3$AUC <- fmt3(s3$AUC)
s3$`AUC 95% CI low` <- fmt3(s3$`AUC 95% CI low`)
s3$`AUC 95% CI high` <- fmt3(s3$`AUC 95% CI high`)
s3$`PR-AUC` <- fmt3(s3$`PR-AUC`)
s3$`PR-AUC baseline` <- fmt3(s3$`PR-AUC baseline`)
s3$`Wilcoxon P` <- signif(suppressWarnings(as.numeric(s3$`Wilcoxon P`)), 3)
s3$`Wilcoxon q` <- signif(suppressWarnings(as.numeric(s3$`Wilcoxon q`)), 3)
s3$`Analysis set` <- clean_text(s3$`Analysis set`)
s3$`Score status` <- clean_text(s3$`Score status`)
s3$`Performance role` <- clean_text(s3$`Performance role`)
s3 <- s3[order(s3$`Analysis set`, s3$`Dataset accession`, s3$`Classifier ID`), ]

# ---- Supplementary Table S4: robustness, subtype, calibration, failure modes ----
s4a <- sensitivity[, intersect(c(
  "sensitivity_scenario", "classifier_id", "classifier_score_name", "n_datasets",
  "median_auc", "min_auc", "max_auc", "rank", "rank_shift_vs_primary",
  "n_auc_ge_0_60", "n_auc_ge_0_65"
), names(sensitivity)), drop = FALSE]
names(s4a) <- clean_text(names(s4a))
num_cols_s4a <- intersect(c("median auc", "min auc", "max auc"), names(s4a))
for (cc in num_cols_s4a) s4a[[cc]] <- fmt3(s4a[[cc]])

s4b <- subtype[, intersect(c(
  "subtype_analysis_set", "subtype_group", "classifier_id", "classifier_score_name",
  "n_datasets", "median_auc", "min_auc", "max_auc", "median_pr_auc",
  "n_auc_ge_0_60", "n_auc_ge_0_65"
), names(subtype)), drop = FALSE]
names(s4b) <- clean_text(names(s4b))
num_cols_s4b <- intersect(c("median auc", "min auc", "max auc", "median pr auc"), names(s4b))
for (cc in num_cols_s4b) s4b[[cc]] <- fmt3(s4b[[cc]])

s4c <- cal_fail[, c(
  "classifier_id", "classifier_score_name", "n_external_core_datasets",
  "median_brier", "median_abs_calibration_intercept", "median_calibration_slope",
  "median_observed_expected_ratio", "n_abs_intercept_gt_0_5",
  "n_slope_outside_0_5_1_5", "n_oe_outside_0_8_1_25", "calibration_failure_class"
)]
names(s4c) <- c(
  "Classifier ID", "Classifier score", "External core datasets, n",
  "Median Brier score", "Median absolute calibration intercept", "Median calibration slope",
  "Median observed/expected ratio", "Datasets with |intercept| >0.5",
  "Datasets with slope outside 0.5-1.5", "Datasets with O/E outside 0.8-1.25", "Calibration failure class"
)
for (cc in c("Median Brier score", "Median absolute calibration intercept", "Median calibration slope", "Median observed/expected ratio")) {
  s4c[[cc]] <- fmt3(s4c[[cc]])
}
s4c$`Calibration failure class` <- clean_text(s4c$`Calibration failure class`)

s4d <- failure[, intersect(c(
  "classifier_id", "platform", "analysis_set", "score_status",
  "n_comparisons", "median_auc", "min_auc", "max_auc",
  "n_auc_below_0_55", "n_low_coverage", "n_not_primary",
  "n_total_samples", "failure_mode_flag"
), names(failure)), drop = FALSE]
names(s4d) <- clean_text(names(s4d))
for (cc in intersect(c("median auc", "min auc", "max auc"), names(s4d))) {
  s4d[[cc]] <- fmt3(s4d[[cc]])
}

manifest_table <- data.frame()
manifest_table <- rbind(manifest_table, write_one_table(
  table2,
  "Table 2. Published transcriptomic classifiers and reconstruction status",
  "Table 2",
  "Table_2_classifier_reconstruction_status_final_v1"
))
manifest_table <- rbind(manifest_table, write_one_table(
  table3,
  "Table 3. Locked external-validation transferability ranking",
  "Table 3",
  "Table_3_locked_external_validation_ranking_final_v1"
))
manifest_table <- rbind(manifest_table, write_one_table(
  table4,
  "Table 4. Calibration instability and failure-mode summary",
  "Table 4",
  "Table_4_calibration_failure_summary_final_v1"
))
manifest_table <- rbind(manifest_table, write_one_table(
  s1,
  "Supplementary Table S1. Cohort harmonization and expression availability",
  file.path("supplementary table", "Supplementary Table S1"),
  "Supplementary_Table_S1_cohort_harmonization_expression_availability_final_v1"
))
manifest_table <- rbind(manifest_table, write_one_table(
  s2,
  "Supplementary Table S2. Classifier gene coverage by dataset",
  file.path("supplementary table", "Supplementary Table S2"),
  "Supplementary_Table_S2_classifier_gene_coverage_final_v1"
))
manifest_table <- rbind(manifest_table, write_one_table(
  s3,
  "Supplementary Table S3. Full external validation performance",
  file.path("supplementary table", "Supplementary Table S3"),
  "Supplementary_Table_S3_full_external_validation_performance_final_v1"
))
manifest_table <- rbind(manifest_table, write_workbook(
  list(
    "S4a_sensitivity" = s4a,
    "S4b_subtype" = s4b,
    "S4c_calibration" = s4c,
    "S4d_failure_modes" = s4d
  ),
  "Supplementary Table S4. Robustness, subtype, calibration, and failure-mode results",
  file.path("supplementary table", "Supplementary Table S4"),
  "Supplementary_Table_S4_robustness_subtype_calibration_failure_modes_final_v1"
))

manifest_out <- file.path(out_root, "final_article_table_file_manifest_v1.csv")
write.csv(manifest_table, manifest_out, row.names = FALSE, na = "")
write.csv(manifest_table, file.path(desktop_package, "final_article_table_file_manifest_v1.csv"), row.names = FALSE, na = "")

writeLines(capture.output(sessionInfo()), session_file)
cat("Wrote output manifest:", manifest_out, "\n")
cat("Wrote session info:", session_file, "\n")
cat("Completed final table build.\n")

sink()
