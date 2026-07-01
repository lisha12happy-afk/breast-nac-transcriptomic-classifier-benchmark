args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
step2_dir <- file.path(project_dir, "manifests", "step 2")
step3_dir <- file.path(project_dir, "manifests", "step 3")

dir.create(step3_dir, recursive = TRUE, showWarnings = FALSE)

cat("build_classifier_coverage_audit_v1.R\n")
cat("Project directory:", project_dir, "\n")
cat("Step 2 directory:", step2_dir, "\n")
cat("Step 3 directory:", step3_dir, "\n")
cat("R version:", R.version.string, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

required_packages <- c("data.table")
for (pkg in required_packages) {
  cat("Package", pkg, "available:", requireNamespace(pkg, quietly = TRUE), "\n")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package data.table is required but not installed.")
}
library(data.table)

gene_availability_file <- file.path(step2_dir, "gene_availability_by_dataset_v1.csv")
expression_audit_file <- file.path(step2_dir, "expression_matrix_audit_v1.csv")
genefu_tarball <- file.path(step3_dir, "genefu_2.44.0.tar.gz")

cat("\nInput files:\n")
cat("gene_availability_file exists:", file.exists(gene_availability_file), "\n")
cat("expression_audit_file exists:", file.exists(expression_audit_file), "\n")
cat("genefu_tarball exists:", file.exists(genefu_tarball), "\n")

if (!file.exists(gene_availability_file)) {
  stop("Missing gene availability file: ", gene_availability_file)
}
if (!file.exists(expression_audit_file)) {
  stop("Missing expression audit file: ", expression_audit_file)
}
if (!file.exists(genefu_tarball)) {
  stop("Missing genefu source tarball: ", genefu_tarball)
}

gene_availability <- fread(gene_availability_file)
expression_audit <- fread(expression_audit_file)

matrix_cols <- setdiff(names(gene_availability), c("gene_symbol", "n_matrices_present"))
cat("\nLoaded gene availability rows:", nrow(gene_availability), "\n")
cat("Matrix columns:", paste(matrix_cols, collapse = ", "), "\n")
cat("Loaded expression audit rows:", nrow(expression_audit), "\n")

normalize_gene_symbol <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y %in% c("", "NA", "N/A", "NULL", "NONE", "<NA>")] <- NA_character_
  y
}

apply_alias_map <- function(x) {
  y <- normalize_gene_symbol(x)
  alias_map <- c(
    "CTSL2" = "CTSV",
    "ORC6L" = "ORC6",
    "GOLPH2" = "GOLM1",
    "SEPT9" = "SEPTIN9"
  )
  hit <- !is.na(y) & y %in% names(alias_map)
  y[hit] <- unname(alias_map[y[hit]])
  y
}

extract_genefu_object <- function(tarfile, rda_path) {
  tmp <- tempfile("genefu_src_")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(tarfile, files = rda_path, exdir = tmp)
  f <- file.path(tmp, rda_path)
  if (!file.exists(f)) {
    stop("Could not extract ", rda_path)
  }
  env <- new.env(parent = emptyenv())
  loaded <- load(f, envir = env)
  if (length(loaded) != 1) {
    stop("Expected one object in ", rda_path, " but loaded: ", paste(loaded, collapse = ", "))
  }
  get(loaded[[1]], envir = env)
}

make_rows <- function(classifier_id, classifier_name, genes, source_object,
                      gene_role = NA_character_, source_probe = NA_character_,
                      weight = NA_real_, source_note = NA_character_) {
  n <- length(genes)
  if (length(gene_role) == 1) gene_role <- rep(gene_role, n)
  if (length(source_probe) == 1) source_probe <- rep(source_probe, n)
  if (length(weight) == 1) weight <- rep(weight, n)
  if (length(source_note) == 1) source_note <- rep(source_note, n)
  data.table(
    classifier_id = classifier_id,
    classifier_name = classifier_name,
    gene_symbol_original = as.character(genes),
    gene_symbol_harmonized = apply_alias_map(genes),
    gene_role = as.character(gene_role),
    weight = as.numeric(weight),
    source_probe = as.character(source_probe),
    source_object = source_object,
    source_note = as.character(source_note)
  )
}

cat("\nExtracting genefu signature objects from source tarball...\n")
pam50 <- extract_genefu_object(genefu_tarball, "genefu/data/pam50.rda")
sig_gene70 <- extract_genefu_object(genefu_tarball, "genefu/data/sig.gene70.rda")
sig_oncotypedx <- extract_genefu_object(genefu_tarball, "genefu/data/sig.oncotypedx.rda")
sig_ggi <- extract_genefu_object(genefu_tarball, "genefu/data/sig.ggi.rda")
sig_gene76 <- extract_genefu_object(genefu_tarball, "genefu/data/sig.gene76.rda")
sig_endo <- extract_genefu_object(genefu_tarball, "genefu/data/sig.endoPredict.rda")

gene_lists <- rbindlist(list(
  make_rows(
    "PAM50",
    "PAM50 intrinsic subtype classifier",
    pam50$centroids.map$probe,
    "genefu::pam50$centroids.map",
    source_probe = pam50$centroids.map$probe,
    source_note = "50-gene centroid classifier from genefu source object"
  ),
  make_rows(
    "MAMMAPRINT_70",
    "MammaPrint / 70-gene signature",
    sig_gene70$HUGO.gene.symbol,
    "genefu::sig.gene70",
    gene_role = "70-gene prognosis signature",
    source_probe = sig_gene70$probe,
    source_note = ifelse(is.na(sig_gene70$HUGO.gene.symbol), "No HUGO symbol in genefu object; cannot audit by gene symbol", "HUGO symbol from genefu object")
  ),
  make_rows(
    "ONCOTYPE_DX_21",
    "Oncotype DX 21-gene recurrence score gene set",
    sig_oncotypedx$symbol,
    "genefu::sig.oncotypedx",
    gene_role = sig_oncotypedx$group,
    source_probe = sig_oncotypedx$probe.affy,
    weight = sig_oncotypedx$weight,
    source_note = "Gene set and group weights from genefu object; commercial assay algorithm is not reconstructed"
  ),
  make_rows(
    "GGI_128_PROBE",
    "Genomic Grade Index implementation in genefu",
    sig_ggi$HUGO.gene.symbol,
    "genefu::sig.ggi",
    gene_role = paste0("grade_", sig_ggi$grade),
    source_probe = sig_ggi$probe,
    source_note = ifelse(is.na(sig_ggi$HUGO.gene.symbol), "No HUGO symbol in genefu object; cannot audit by gene symbol", "HUGO symbol from genefu object")
  ),
  make_rows(
    "GENE76",
    "Wang 76-gene distant-metastasis signature",
    sig_gene76$HUGO.gene.symbol,
    "genefu::sig.gene76",
    gene_role = ifelse(sig_gene76$er == 1, "ER_positive_component", "ER_negative_component"),
    source_probe = sig_gene76$probe,
    weight = sig_gene76$std.cox.coefficient,
    source_note = ifelse(is.na(sig_gene76$HUGO.gene.symbol), "No HUGO symbol in genefu object; cannot audit by gene symbol", "HUGO symbol from genefu object")
  ),
  make_rows(
    "ENDOPREDICT_11",
    "EndoPredict 11-gene set",
    sig_endo$symbol,
    "genefu::sig.endoPredict",
    gene_role = sig_endo$group,
    source_probe = sig_endo$probe.affy,
    weight = sig_endo$weight,
    source_note = "Gene set and weights from genefu object; clinical assay calibration is not reconstructed"
  ),
  make_rows(
    "IFNG_18",
    "IFN-gamma-related 18-gene immune profile",
    c("STAT1", "IDO1", "CXCL10", "CXCL9", "HLA-DRA", "GZMB", "IFNG", "LAG3", "CD8A",
      "CCL5", "CXCL13", "HLA-E", "NKG7", "TIGIT", "PSMB10", "CD27", "CD274", "PDCD1LG2"),
    "manual_from_primary_publication",
    gene_role = "immune_response",
    source_note = "Manually encoded from the published IFN-gamma-related mRNA profile"
  ),
  make_rows(
    "CYTOLYTIC_ACTIVITY_2",
    "Cytolytic activity score",
    c("GZMA", "PRF1"),
    "manual_from_primary_publication",
    gene_role = "cytolytic_effector",
    source_note = "Two-gene cytolytic activity score"
  )
), fill = TRUE)

gene_lists[, row_in_classifier_source := seq_len(.N), by = classifier_id]

dictionary <- data.table(
  classifier_id = c(
    "PAM50", "MAMMAPRINT_70", "ONCOTYPE_DX_21", "GGI_128_PROBE", "GENE76", "ENDOPREDICT_11",
    "IFNG_18", "CYTOLYTIC_ACTIVITY_2", "HESS_PCR_30", "HATZIS_PCR_PREDICTOR",
    "HER2DX", "TNBCTYPE"
  ),
  classifier_name = c(
    "PAM50 intrinsic subtype classifier",
    "MammaPrint / 70-gene signature",
    "Oncotype DX 21-gene recurrence score gene set",
    "Genomic Grade Index implementation in genefu",
    "Wang 76-gene distant-metastasis signature",
    "EndoPredict 11-gene set",
    "IFN-gamma-related 18-gene immune profile",
    "Cytolytic activity score",
    "Hess 30-gene pharmacogenomic pCR predictor",
    "Hatzis taxane-anthracycline genomic response/survival predictor",
    "HER2DX genomic test",
    "TNBCtype expression subtype classifier"
  ),
  original_endpoint_or_use = c(
    "Intrinsic subtype / ROR-related biology; not originally a pCR classifier",
    "Distant metastasis / prognosis; not originally a pCR classifier",
    "Recurrence score in ER-positive breast cancer; not originally a pCR classifier",
    "Histologic grade / prognosis; not originally a pCR classifier",
    "Distant metastasis in lymph-node-negative breast cancer; not originally a pCR classifier",
    "Distant recurrence in ER-positive HER2-negative breast cancer; not originally a pCR classifier",
    "Immune response profile linked to anti-PD-1 response; mechanism comparator for pCR biology",
    "Tumor immune cytolytic activity; mechanism comparator for pCR biology",
    "Pathologic complete response to paclitaxel followed by FAC in breast cancer",
    "pCR and survival after taxane-anthracycline chemotherapy for invasive breast cancer",
    "HER2-positive early breast cancer genomic test; pCR/biology associations reported",
    "Triple-negative breast cancer expression subtype assignment"
  ),
  pcr_specific = c("no", "no", "no", "no", "no", "no", "no", "no", "yes", "yes", "partly", "no"),
  primary_scope_for_this_project = c(
    "full_benchmark_and_subtype_context",
    "full_benchmark",
    "full_benchmark_with_HR_positive_caveat",
    "full_benchmark",
    "full_benchmark",
    "HR_positive_context_only",
    "TNBC_HER2_positive_immune_sensitivity",
    "TNBC_HER2_positive_immune_sensitivity",
    "literature_anchor_not_implemented_v1",
    "literature_anchor_not_implemented_v1",
    "HER2_positive_context_not_implemented_v1",
    "TNBC_context_not_implemented_v1"
  ),
  gene_list_status = c(
    "available_from_genefu_source",
    "available_from_genefu_source_with_unmapped_original_contigs",
    "available_from_genefu_source",
    "available_from_genefu_source_with_unmapped_probes",
    "available_from_genefu_source_with_unmapped_probes",
    "available_from_genefu_source",
    "manual_public_gene_list",
    "manual_public_gene_list",
    "not_encoded_v1_requires_supplement_formula_review",
    "not_encoded_v1_requires_supplement_formula_review",
    "not_encoded_v1_proprietary_or_not_fully_public",
    "not_encoded_v1_requires_classifier_centroid_software_review"
  ),
  formula_reconstruction_status = c(
    "centroid_classifier_possible_but_requires_preprocessing_decision",
    "gene_set_level_possible_clinical_cutoff_not_reconstructed",
    "gene_set_level_possible_commercial_score_not_reconstructed",
    "score_possible_from_public/genefu_implementation_after_preprocessing_decision",
    "score_possible_from_public/genefu_implementation_after_preprocessing_decision",
    "gene_set_level_possible_clinical_calibration_not_reconstructed",
    "simple_gene_set_score_possible",
    "simple_two_gene_score_possible",
    "not_reconstructed_in_v1",
    "not_reconstructed_in_v1",
    "not_reconstructed_in_v1",
    "not_reconstructed_in_v1"
  ),
  primary_reference_pmid = c(
    "19204204", "11823860", "15591335", "16478745", "15721472", "21807638",
    "28650338", "25594174", "16896004", "21558518", "41324567", "21633166"
  ),
  primary_reference_doi = c(
    "10.1200/JCO.2008.18.1370",
    "10.1038/415530a",
    "10.1056/NEJMoa041588",
    "10.1093/jnci/djj052",
    "10.1016/S0140-6736(05)17947-1",
    "10.1158/1078-0432.CCR-11-0926",
    "10.1172/JCI91190",
    "10.1016/j.cell.2014.12.033",
    "10.1200/JCO.2006.05.6861",
    "10.1001/jama.2011.593",
    "10.1158/1078-0432.CCR-25-3123",
    "10.1172/JCI45014"
  ),
  primary_reference_url = paste0("https://pubmed.ncbi.nlm.nih.gov/", c(
    "19204204", "11823860", "15591335", "16478745", "15721472", "21807638",
    "28650338", "25594174", "16896004", "21558518", "41324567", "21633166"
  ), "/"),
  gene_list_source_url = c(
    rep("https://bioconductor.org/packages/release/bioc/html/genefu.html", 6),
    "https://pubmed.ncbi.nlm.nih.gov/28650338/",
    "https://pubmed.ncbi.nlm.nih.gov/25594174/",
    NA_character_, NA_character_, NA_character_, NA_character_
  ),
  include_in_step3_coverage = c(rep("yes", 8), rep("no", 4)),
  notes = c(
    "Useful as subtype/biology transferability benchmark; do not claim pCR predictor unless performance supports it.",
    "High clinical familiarity; use as transferability stress test, not as de novo pCR signature.",
    "Commercial score not reconstructed; gene-set availability still useful for platform coverage audit.",
    "Captures proliferation/grade signal; likely strong pCR correlate but not pCR-specific.",
    "Legacy prognostic signature; include mainly as robustness/proliferation comparator.",
    "ER-positive/HER2-negative recurrence context; not central for all-NAC pCR benchmark.",
    "Immune comparator relevant to TNBC/HER2+ pCR biology; not breast-NAC-specific classifier.",
    "Minimal immune cytotoxicity comparator; highly reproducible but too small as standalone classifier.",
    "Important prior pCR predictor; do not implement until exact gene/probe list, preprocessing and classifier rule are recovered.",
    "Important prior pCR predictor; do not implement until exact model object/rule and preprocessing are recovered.",
    "Recent HER2+ assay context; likely proprietary/commercial details limit independent reconstruction.",
    "Subtype classifier may be useful for TNBC subset, but needs separate centroid/software reproducibility review."
  )
)

coverage_gene_lists <- gene_lists[!is.na(gene_symbol_harmonized)]
coverage_gene_lists_unique <- unique(
  coverage_gene_lists[, .(
    classifier_id,
    classifier_name,
    gene_symbol_harmonized
  )]
)

gene_availability_long <- melt(
  gene_availability,
  id.vars = "gene_symbol",
  measure.vars = matrix_cols,
  variable.name = "matrix_label",
  value.name = "present"
)
gene_availability_long[, gene_symbol_harmonized := apply_alias_map(gene_symbol)]
gene_availability_long[, present := as.integer(present)]

coverage_expanded <- copy(coverage_gene_lists_unique)
coverage_expanded[, cartesian_join_key := 1L]
matrix_grid <- data.table(matrix_label = matrix_cols, cartesian_join_key = 1L)
coverage_expanded <- merge(
  coverage_expanded,
  matrix_grid,
  by = "cartesian_join_key",
  allow.cartesian = TRUE
)
coverage_expanded[, cartesian_join_key := NULL]
coverage_expanded <- merge(
  coverage_expanded,
  gene_availability_long[, .(matrix_label, gene_symbol_harmonized, present)],
  by = c("matrix_label", "gene_symbol_harmonized"),
  all.x = TRUE
)
coverage_expanded[is.na(present), present := 0L]

coverage_by_matrix <- coverage_expanded[, .(
  n_harmonized_genes = uniqueN(gene_symbol_harmonized),
  n_genes_present = sum(present == 1L),
  missing_genes = paste(sort(unique(gene_symbol_harmonized[present != 1L])), collapse = ";")
), by = .(classifier_id, classifier_name, matrix_label)]

total_rows <- gene_lists[, .(
  n_source_rows_total = .N,
  n_source_rows_with_harmonized_gene_symbol = sum(!is.na(gene_symbol_harmonized)),
  n_unique_harmonized_genes = uniqueN(gene_symbol_harmonized[!is.na(gene_symbol_harmonized)])
), by = .(classifier_id, classifier_name)]

coverage_by_matrix <- merge(
  coverage_by_matrix,
  total_rows,
  by = c("classifier_id", "classifier_name"),
  all.x = TRUE
)
coverage_by_matrix[, coverage_pct_of_harmonized := round(100 * n_genes_present / n_harmonized_genes, 1)]
coverage_by_matrix[, coverage_grade := fifelse(
  coverage_pct_of_harmonized >= 95, "excellent",
  fifelse(coverage_pct_of_harmonized >= 80, "usable_with_minor_missingness",
          fifelse(coverage_pct_of_harmonized >= 50, "limited", "poor"))
)]

coverage_by_matrix <- merge(
  coverage_by_matrix,
  expression_audit[, .(
    dataset_accession, matrix_label, platform, planned_role, story_scope,
    usable_for_full_benchmark, usable_for_subtype_only_validation, expression_scale_guess
  )],
  by = "matrix_label",
  all.x = TRUE
)

coverage_by_matrix <- coverage_by_matrix[
  order(classifier_id, dataset_accession, matrix_label)
]

coverage_by_matrix[, analysis_set := fifelse(
  usable_for_full_benchmark == "yes", "core_full_benchmark",
  fifelse(usable_for_full_benchmark == "stress_test_only", "large_stress_test",
          fifelse(grepl("limited", usable_for_subtype_only_validation, ignore.case = TRUE), "limited_subtype_panel_or_gene_space",
                  "subtype_validation_only"))
)]

not_computable_coverage <- dictionary[include_in_step3_coverage == "no", .(
  matrix_label = NA_character_,
  classifier_id,
  classifier_name,
  n_harmonized_genes = NA_integer_,
  n_genes_present = NA_integer_,
  missing_genes = NA_character_,
  n_source_rows_total = NA_integer_,
  n_source_rows_with_harmonized_gene_symbol = NA_integer_,
  n_unique_harmonized_genes = NA_integer_,
  coverage_pct_of_harmonized = NA_real_,
  coverage_grade = "not_computable_v1",
  dataset_accession = NA_character_,
  platform = NA_character_,
  planned_role = NA_character_,
  story_scope = NA_character_,
  usable_for_full_benchmark = NA_character_,
  usable_for_subtype_only_validation = NA_character_,
  expression_scale_guess = NA_character_
)]

coverage_final <- rbindlist(list(coverage_by_matrix, not_computable_coverage), fill = TRUE)

risk <- dictionary[, .(
  classifier_id,
  classifier_name,
  pcr_specific,
  include_in_step3_coverage,
  gene_list_status,
  formula_reconstruction_status,
  endpoint_match_for_this_project = fifelse(
    pcr_specific == "yes", "direct_pcr_rd_endpoint_match",
    fifelse(pcr_specific == "partly", "subtype_or_context_specific_pcr_relevance", "indirect_biology_or_prognostic_comparator")
  ),
  reproducibility_risk = fifelse(
    include_in_step3_coverage == "no", "high",
    fifelse(grepl("commercial|cutoff|calibration|preprocessing", formula_reconstruction_status, ignore.case = TRUE), "medium", "low_to_medium")
  ),
  recommended_step4_action = fifelse(
    include_in_step3_coverage == "no",
    "do_not_score_in_main_analysis_until_exact_formula_and_inputs_are_recovered",
    fifelse(classifier_id %in% c("IFNG_18", "CYTOLYTIC_ACTIVITY_2"),
            "use_as_secondary_mechanistic_gene_set_score_after locking normalization",
            "eligible_for_platform_coverage_pass; scoring requires outcome-blind preprocessing rule")
  ),
  primary_reference_pmid,
  primary_reference_doi,
  primary_reference_url,
  notes
)]

summary_by_classifier <- coverage_by_matrix[, .(
  n_matrices_audited = .N,
  min_coverage_pct = min(coverage_pct_of_harmonized, na.rm = TRUE),
  median_coverage_pct = median(coverage_pct_of_harmonized, na.rm = TRUE),
  n_matrices_excellent = sum(coverage_grade == "excellent"),
  n_matrices_usable_or_better = sum(coverage_pct_of_harmonized >= 80),
  matrices_below_80pct = paste(matrix_label[coverage_pct_of_harmonized < 80], collapse = ";")
), by = .(classifier_id, classifier_name)]
summary_by_classifier <- merge(
  summary_by_classifier,
  total_rows,
  by = c("classifier_id", "classifier_name"),
  all.x = TRUE
)
summary_by_classifier <- summary_by_classifier[order(classifier_id)]

summary_by_analysis_set <- coverage_by_matrix[, .(
  n_matrices_audited = .N,
  min_coverage_pct = min(coverage_pct_of_harmonized, na.rm = TRUE),
  median_coverage_pct = median(coverage_pct_of_harmonized, na.rm = TRUE),
  n_matrices_usable_or_better = sum(coverage_pct_of_harmonized >= 80),
  matrices_below_80pct = paste(matrix_label[coverage_pct_of_harmonized < 80], collapse = ";")
), by = .(classifier_id, classifier_name, analysis_set)]
summary_by_analysis_set <- summary_by_analysis_set[order(classifier_id, analysis_set)]

write_output <- function(x, filename) {
  out <- file.path(step3_dir, filename)
  fwrite(x, out)
  cat("Wrote:", out, "rows:", nrow(x), "\n")
}

write_output(dictionary, "classifier_dictionary_v1.csv")
write_output(gene_lists, "classifier_gene_lists_v1.csv")
write_output(coverage_final, "classifier_gene_coverage_v1.csv")
write_output(risk, "classifier_reproducibility_risk_v1.csv")
write_output(summary_by_classifier, "classifier_coverage_summary_v1.csv")
write_output(summary_by_analysis_set, "classifier_coverage_by_analysis_set_v1.csv")

decision_lines <- c(
  "# Step 3 classifier selection and gene-coverage decision v1",
  "",
  paste0("Run timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Objective",
  "",
  "This step converts the project from ordinary signature construction to an external-validation and transferability benchmark. The goal is to identify published breast-cancer transcriptomic classifiers or mechanistic scores that are sufficiently transparent to audit across the Step 2 expression matrices.",
  "",
  "## Main decision",
  "",
  "Use a two-tier classifier set:",
  "",
  "1. Main computable benchmark set: PAM50, MammaPrint/70-gene, OncotypeDX 21-gene set, Genomic Grade Index, Wang 76-gene signature, EndoPredict, IFN-gamma 18-gene immune profile, and cytolytic activity score.",
  "2. Literature-anchor but not-yet-computable set: Hess pCR 30-gene predictor, Hatzis pCR/survival predictor, HER2DX, and TNBCtype. These are important for the Introduction/Discussion and for showing that pCR prediction has been attempted, but they should not be scored in the main analysis until their exact public formula, feature preprocessing, and cutoffs are recovered.",
  "",
  "## Why this is not a new-signature project",
  "",
  "No new pCR signature is selected or optimized here. Step 3 only asks whether published classifiers can be transferred across the already-audited public NAC cohorts. This protects the paper frame: external validation, calibration, transferability, and failure modes.",
  "",
  "## Immediate Step 4 consequence",
  "",
  "Before any AUC/calibration analysis, lock an outcome-blind preprocessing rule per platform: probe-to-gene collapsing, log2 handling, z-score/quantile standardization choice, and missing-gene rule. Do not tune these choices against pCR/RD in validation cohorts.",
  "",
  "## Coverage summary",
  "",
  paste(capture.output(print(summary_by_classifier)), collapse = "\n"),
  "",
  "## Key reproducibility caveats",
  "",
  "- Commercial or clinical-use assays may have public gene lists but non-public weights, calibration, or cutoff rules. Treat those as gene-set or approximate-score benchmarks unless the exact algorithm is public.",
  "- Legacy signatures contain unmapped probes/contigs or old gene symbols. Coverage percentages are calculated against harmonized gene symbols available in Step 2, not against original probe-level assay content.",
  "- Hess and Hatzis are directly relevant pCR papers, but Step 3 v1 does not reconstruct them because the exact scoring rule was not recovered into an auditable local gene list during this pass.",
  "- GSE106977 and GSE109710 are expected stress points because Step 2 already found limited gene space / panel constraints.",
  "",
  "## Output files",
  "",
  "- classifier_dictionary_v1.csv",
  "- classifier_gene_lists_v1.csv",
  "- classifier_gene_coverage_v1.csv",
  "- classifier_coverage_summary_v1.csv",
  "- classifier_coverage_by_analysis_set_v1.csv",
  "- classifier_reproducibility_risk_v1.csv",
  "- build_classifier_coverage_audit_v1.R",
  "- R_raw_run_log_step3_v1.txt",
  "- R_session_info_step3_v1.txt"
)
writeLines(decision_lines, file.path(step3_dir, "classifier_selection_decision_v1.md"), useBytes = TRUE)
cat("Wrote:", file.path(step3_dir, "classifier_selection_decision_v1.md"), "\n")

session_out <- file.path(step3_dir, "R_session_info_step3_v1.txt")
capture.output(sessionInfo(), file = session_out)
cat("Wrote:", session_out, "\n")

cat("\nClassifier coverage summary:\n")
print(summary_by_classifier)

cat("\nDone.\n")
