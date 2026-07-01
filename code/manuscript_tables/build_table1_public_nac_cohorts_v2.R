options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "manifests", "table1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "R_raw_run_log_table1_v2.txt")
session_file <- file.path(out_dir, "R_session_info_table1_v2.txt")
sink(log_file, split = TRUE)

cat("Build Table 1: public NAC pretreatment transcriptome cohorts\n")
cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Project root:", project_root, "\n")

step1_dir <- file.path(project_root, "manifests", "step 1")
step2_dir <- file.path(project_root, "manifests", "step 2")

manifest_path <- file.path(step1_dir, "analysis_ready_manifest_v1.csv")
summary_path <- file.path(step1_dir, "manifest_summary_v1.csv")
expr_path <- file.path(step2_dir, "expression_matrix_audit_v1.csv")

manifest <- read.csv(manifest_path, na.strings = c("", "NA"), check.names = FALSE)
summary_tbl <- read.csv(summary_path, na.strings = c("", "NA"), check.names = FALSE)
expr <- read.csv(expr_path, na.strings = c("", "NA"), check.names = FALSE)

cat("Loaded analysis-ready manifest rows:", nrow(manifest), "\n")
cat("Loaded dataset-level summary rows:", nrow(summary_tbl), "\n")
cat("Loaded expression-audit rows:", nrow(expr), "\n")

required_manifest <- c(
  "dataset_accession", "platform", "endpoint_binary_pcr1_rd0",
  "planned_role", "story_scope", "include_primary_analysis", "overlap_flag",
  "curation_notes"
)
required_expr <- c(
  "dataset_accession", "platform", "n_samples_matched_to_manifest",
  "planned_role", "story_scope", "usable_for_full_benchmark",
  "usable_for_subtype_only_validation"
)

missing_manifest <- setdiff(required_manifest, names(manifest))
missing_expr <- setdiff(required_expr, names(expr))
if (length(missing_manifest) > 0) {
  stop("Missing columns in analysis_ready_manifest_v1.csv: ", paste(missing_manifest, collapse = ", "))
}
if (length(missing_expr) > 0) {
  stop("Missing columns in expression_matrix_audit_v1.csv: ", paste(missing_expr, collapse = ", "))
}

manifest$key <- paste(manifest$dataset_accession, manifest$platform, sep = "||")
expr$key <- paste(expr$dataset_accession, expr$platform, sep = "||")
manifest$endpoint_binary_clean <- suppressWarnings(as.integer(manifest$endpoint_binary_pcr1_rd0))

split_manifest <- split(manifest, manifest$key, drop = TRUE)
summary_by_platform <- do.call(
  rbind,
  lapply(names(split_manifest), function(k) {
    dat <- split_manifest[[k]]
    endpoint_valid <- dat$endpoint_binary_clean %in% c(0L, 1L)
    data.frame(
      key = k,
      dataset_accession = dat$dataset_accession[1],
      platform = dat$platform[1],
      manifest_samples_n = nrow(dat),
      endpoint_mapped_samples_n = sum(endpoint_valid, na.rm = TRUE),
      pCR_n = sum(dat$endpoint_binary_clean == 1L, na.rm = TRUE),
      RD_n = sum(dat$endpoint_binary_clean == 0L, na.rm = TRUE),
      include_values_observed = paste(sort(unique(na.omit(dat$include_primary_analysis))), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
)

role_label <- function(role, scope, dataset, platform) {
  if (dataset == "GSE25066") return("Discovery/model-reconstruction base")
  if (dataset == "GSE194040") return("Large stress test (I-SPY2)")
  if (dataset == "GSE163882") return("External validation 6 (RNA-seq FFPE)")
  if (dataset == "GSE106977") return("Subtype-only validation: TNBC/ER-negative")
  if (dataset %in% c("GSE109710", "GSE130786")) return("Subtype-only validation: HER2-positive")
  if (dataset == "GSE32646") return("External validation 2; TNBC/ER-negative subset")
  if (dataset == "GSE50948") return("External validation 4; HER2-positive subset")
  if (dataset == "GSE66305") return("External validation 5; HER2-positive subset")
  if (dataset == "GSE20271") return("External validation 1")
  if (dataset == "GSE41998") return("External validation 3")
  gsub("_", " ", role)
}

scope_label <- function(scope) {
  if (is.na(scope) || !nzchar(scope)) return("")
  x <- scope
  x <- gsub("full_benchmark_stress", "Full benchmark stress", x)
  x <- gsub("full_benchmark", "Full benchmark", x)
  x <- gsub("subtype_TNBC_ER_negative", "TNBC/ER-negative subset", x)
  x <- gsub("subtype_HER2_positive", "HER2-positive subset", x)
  x <- gsub(";", "; ", x)
  x <- gsub("_", " ", x)
  x
}

primary_label <- function(dataset, usable_full, endpoint_mapped, manifest_n) {
  if (dataset == "GSE25066") {
    return("Discovery/model reconstruction; endpoint-missing samples excluded")
  }
  if (dataset == "GSE194040") return("Stress test only")
  if (dataset %in% c("GSE106977", "GSE109710", "GSE130786")) return("No; subtype-only exploratory")
  if (!is.na(endpoint_mapped) && !is.na(manifest_n) && endpoint_mapped < manifest_n) {
    return("Yes; endpoint-missing samples excluded")
  }
  if (!is.na(usable_full) && usable_full == "yes") return("Yes")
  "No"
}

key_note <- function(dataset, endpoint_mapped, manifest_n) {
  if (dataset == "GSE25066") return("Selected over GSE25055/GSE25065; endpoint-missing samples excluded.")
  if (dataset == "GSE194040") return("I-SPY2 stress-test cohort; two expression platforms kept separate; one overlapping patient was flagged across platforms.")
  if (dataset == "GSE106977") return("Subtype-only exploratory TNBC/ER-negative validation; limited gene space.")
  if (dataset == "GSE109710") return("Subtype-only exploratory HER2-positive validation; NanoString/limited panel.")
  if (dataset == "GSE130786") return("Subtype-only exploratory HER2-positive validation; baseline samples only.")
  if (!is.na(endpoint_mapped) && !is.na(manifest_n) && endpoint_mapped < manifest_n) return("Endpoint-missing samples excluded from performance analyses.")
  ""
}

idx <- match(expr$key, summary_by_platform$key)
endpoint_mapped <- summary_by_platform$endpoint_mapped_samples_n[idx]
pcr_n <- summary_by_platform$pCR_n[idx]
rd_n <- summary_by_platform$RD_n[idx]
manifest_n <- summary_by_platform$manifest_samples_n[idx]
transcriptome_n <- suppressWarnings(as.integer(expr$n_samples_matched_to_manifest))
transcriptome_n[is.na(transcriptome_n)] <- manifest_n[is.na(transcriptome_n)]
pcr_rate <- ifelse(endpoint_mapped > 0, round(100 * pcr_n / endpoint_mapped, 1), NA_real_)

table1 <- data.frame(
  `Dataset accession` = expr$dataset_accession,
  `Platform` = expr$platform,
  `Analysis role` = mapply(role_label, expr$planned_role, expr$story_scope, expr$dataset_accession, expr$platform, USE.NAMES = FALSE),
  `Analysis scope` = vapply(expr$story_scope, scope_label, character(1)),
  `Transcriptome samples, n` = transcriptome_n,
  `Endpoint-mapped samples, n` = endpoint_mapped,
  `pCR, n` = pcr_n,
  `RD, n` = rd_n,
  `pCR rate, %` = pcr_rate,
  `Included in primary analysis` = mapply(primary_label, expr$dataset_accession, expr$usable_for_full_benchmark, endpoint_mapped, manifest_n, USE.NAMES = FALSE),
  `Key note` = mapply(key_note, expr$dataset_accession, endpoint_mapped, manifest_n, USE.NAMES = FALSE),
  check.names = FALSE
)

preferred_order <- c(
  "GSE25066||GPL96",
  "GSE20271||GPL96",
  "GSE32646||GPL570",
  "GSE41998||GPL571",
  "GSE50948||GPL570",
  "GSE66305||GPL570",
  "GSE163882||GPL18573",
  "GSE194040||GPL20078",
  "GSE194040||GPL30493",
  "GSE106977||GPL17586",
  "GSE109710||GPL24546",
  "GSE130786||GPL6480"
)
table1$key <- expr$key
table1 <- table1[order(match(table1$key, preferred_order)), ]
table1$key <- NULL

if (any(is.na(table1$`Endpoint-mapped samples, n`))) {
  stop("Endpoint summary missing for: ", paste(table1$`Dataset accession`[is.na(table1$`Endpoint-mapped samples, n`)], collapse = ", "))
}
if (any(table1$`Endpoint-mapped samples, n` != table1$`pCR, n` + table1$`RD, n`)) {
  bad <- table1$`Dataset accession`[table1$`Endpoint-mapped samples, n` != table1$`pCR, n` + table1$`RD, n`]
  stop("Endpoint count mismatch for: ", paste(bad, collapse = ", "))
}

audit <- data.frame(
  item = c(
    "source_manifest",
    "source_dataset_summary",
    "source_expression_audit",
    "rows_in_final_table",
    "gse194040_rows",
    "endpoint_count_check"
  ),
  value = c(
    manifest_path,
    summary_path,
    expr_path,
    nrow(table1),
    paste(table1$Platform[table1$`Dataset accession` == "GSE194040"], collapse = "; "),
    "PASS: Endpoint-mapped samples equal pCR + RD for every row"
  ),
  stringsAsFactors = FALSE
)

csv_out <- file.path(out_dir, "Table_1_public_NAC_transcriptome_cohorts_final_v1.csv")
xlsx_out <- file.path(out_dir, "Table_1_public_NAC_transcriptome_cohorts_final_v1.xlsx")
audit_out <- file.path(out_dir, "Table_1_source_audit_v1.csv")

write.csv(table1, csv_out, row.names = FALSE, na = "")
write.csv(audit, audit_out, row.names = FALSE, na = "")

cat("Wrote CSV:", csv_out, "\n")
cat("Wrote audit:", audit_out, "\n")

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Table 1")
  openxlsx::writeData(
    wb, "Table 1",
    "Table 1. Public neoadjuvant chemotherapy breast cancer pretreatment transcriptome cohorts included in the external-validation transferability benchmark.",
    startRow = 1, startCol = 1
  )
  openxlsx::writeData(wb, "Table 1", table1, startRow = 3, startCol = 1)
  title_style <- openxlsx::createStyle(textDecoration = "bold", fontSize = 12)
  header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9EAF7", border = "TopBottomLeftRight", wrapText = TRUE, valign = "center")
  body_style <- openxlsx::createStyle(border = "TopBottomLeftRight", wrapText = TRUE, valign = "top")
  num_style <- openxlsx::createStyle(border = "TopBottomLeftRight", numFmt = "0", valign = "top")
  pct_style <- openxlsx::createStyle(border = "TopBottomLeftRight", numFmt = "0.0", valign = "top")
  openxlsx::addStyle(wb, "Table 1", title_style, rows = 1, cols = 1, gridExpand = TRUE)
  openxlsx::addStyle(wb, "Table 1", header_style, rows = 3, cols = seq_len(ncol(table1)), gridExpand = TRUE)
  openxlsx::addStyle(wb, "Table 1", body_style, rows = 4:(nrow(table1) + 3), cols = seq_len(ncol(table1)), gridExpand = TRUE)
  numeric_cols <- which(names(table1) %in% c("Transcriptome samples, n", "Endpoint-mapped samples, n", "pCR, n", "RD, n"))
  openxlsx::addStyle(wb, "Table 1", num_style, rows = 4:(nrow(table1) + 3), cols = numeric_cols, gridExpand = TRUE, stack = TRUE)
  pct_col <- which(names(table1) == "pCR rate, %")
  openxlsx::addStyle(wb, "Table 1", pct_style, rows = 4:(nrow(table1) + 3), cols = pct_col, gridExpand = TRUE, stack = TRUE)
  openxlsx::setColWidths(wb, "Table 1", cols = 1:ncol(table1), widths = c(16, 12, 34, 34, 16, 18, 10, 10, 12, 34, 55))
  openxlsx::freezePane(wb, "Table 1", firstActiveRow = 4)
  openxlsx::addWorksheet(wb, "Source audit")
  openxlsx::writeData(wb, "Source audit", audit, startRow = 1, startCol = 1)
  openxlsx::setColWidths(wb, "Source audit", cols = 1:2, widths = c(26, 120))
  openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)
  cat("Wrote XLSX with openxlsx:", xlsx_out, "\n")
} else if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list("Table 1" = table1, "Source audit" = audit), xlsx_out)
  cat("Wrote XLSX with writexl:", xlsx_out, "\n")
} else {
  cat("XLSX not written by R: neither openxlsx nor writexl is installed.\n")
}

writeLines(capture.output(sessionInfo()), session_file)
cat("Wrote session info:", session_file, "\n")
cat("Completed Table 1 build.\n")

sink()
