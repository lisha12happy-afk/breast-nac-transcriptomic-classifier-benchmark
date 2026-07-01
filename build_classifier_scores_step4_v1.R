args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
step1_dir <- file.path(project_dir, "manifests", "step 1")
step2_dir <- file.path(project_dir, "manifests", "step 2")
step3_dir <- file.path(project_dir, "manifests", "step 3")
step4_dir <- file.path(project_dir, "manifests", "step 4")

dir.create(step4_dir, recursive = TRUE, showWarnings = FALSE)

cat("build_classifier_scores_step4_v1.R\n")
cat("Project directory:", project_dir, "\n")
cat("Step 1 directory:", step1_dir, "\n")
cat("Step 2 directory:", step2_dir, "\n")
cat("Step 3 directory:", step3_dir, "\n")
cat("Step 4 directory:", step4_dir, "\n")
cat("R version:", R.version.string, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package data.table is required but not installed.")
}
library(data.table)

normalize_token <- function(x) {
  if (length(x) == 0 || is.na(x)) return("")
  x <- trimws(as.character(x))
  x <- gsub("\\s+", " ", x, perl = TRUE)
  x <- gsub('^"|"$', "", x)
  trimws(x)
}

split_tab <- function(line) {
  vapply(strsplit(line, "\t", fixed = TRUE)[[1]], normalize_token, character(1))
}

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

sample_aliases_for_row <- function(row) {
  vals <- unique(c(row$gsm_sample_id, row$sample_title, row$patient_id))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  extra <- character()
  for (v in vals) {
    if (grepl(":", v, fixed = TRUE)) extra <- c(extra, sub(":.*$", "", v))
    if (grepl("^ISPY2_", v)) extra <- c(extra, sub("^ISPY2_", "", v))
    if (grepl("-GPL", v)) extra <- c(extra, sub("-GPL.*$", "", v))
  }
  unique(c(vals, extra))
}

normalize_matrix_sample_id <- function(x) {
  x <- normalize_token(x)
  x <- sub("^ISPY2_", "", x)
  x <- sub("-GPL.*$", "", x)
  x
}

safe_download_gz <- function(url, label) {
  cache_dir <- file.path(step4_dir, "download_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, paste0(gsub("[^A-Za-z0-9]+", "_", label), ".gz"))
  if (file.exists(cache_file) && file.info(cache_file)$size > 0) {
    cat("[download_cache_hit]", label, file.info(cache_file)$size, "bytes", cache_file, "\n")
    return(cache_file)
  }
  tmp <- tempfile(pattern = paste0(gsub("[^A-Za-z0-9]+", "_", label), "_"), fileext = ".gz")
  cat("[download_start]", label, url, "\n")
  utils::download.file(url, tmp, mode = "wb", quiet = TRUE, method = "auto")
  cat("[download_done]", label, file.info(tmp)$size, "bytes", tmp, "\n")
  file.copy(tmp, cache_file, overwrite = TRUE)
  unlink(tmp)
  cache_file
}

read_series_matrix_values <- function(url, label) {
  tmp <- safe_download_gz(url, paste0(label, "_series_matrix"))
  con <- gzfile(tmp, open = "rt")
  on.exit(close(con), add = TRUE)
  begin_line <- NA_integer_
  line_no <- 0L
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    line_no <- line_no + 1L
    if (identical(line, "!series_matrix_table_begin")) {
      begin_line <- line_no
      break
    }
  }
  if (is.na(begin_line)) stop("series_matrix_table_begin not found in ", url)
  cat("[matrix_fread_start]", label, "skip", begin_line, "\n")
  dt <- fread(tmp, skip = begin_line, header = TRUE, fill = TRUE, data.table = TRUE, showProgress = FALSE)
  first_col <- names(dt)[1]
  dt <- dt[!grepl("^!", as.character(get(first_col)))]
  setnames(dt, first_col, "feature_id")
  dt[, feature_id := vapply(feature_id, normalize_token, character(1))]
  sample_cols <- setdiff(names(dt), "feature_id")
  new_names <- c("feature_id", vapply(sample_cols, normalize_token, character(1)))
  setnames(dt, names(dt), new_names)
  cat("[matrix_fread_done]", label, "features", nrow(dt), "samples", ncol(dt) - 1L, "\n")
  dt
}

read_plain_expression_values <- function(url, label, sep = "\t", header_has_feature_col = FALSE, annotation_col = NULL) {
  tmp <- safe_download_gz(url, paste0(label, "_supp_matrix"))
  cat("[supp_fread_start]", label, "\n")
  if (header_has_feature_col) {
    dt <- fread(tmp, sep = sep, header = TRUE, fill = TRUE, data.table = TRUE, showProgress = FALSE)
    first_col <- names(dt)[1]
    setnames(dt, first_col, "feature_id")
    if (!is.null(annotation_col) && annotation_col %in% names(dt)) {
      dt[, (annotation_col) := NULL]
    }
  } else {
    con <- gzfile(tmp, open = "rt")
    first_line <- readLines(con, n = 1, warn = FALSE)
    close(con)
    header_values <- split_tab(first_line)
    dt <- fread(tmp, sep = sep, skip = 1, header = FALSE, fill = TRUE, data.table = TRUE, showProgress = FALSE)
    setnames(dt, names(dt)[1], "feature_id")
    sample_cols <- setdiff(names(dt), "feature_id")
    if (length(header_values) == ncol(dt)) {
      candidate_names <- header_values[-1]
    } else if (length(header_values) == length(sample_cols)) {
      candidate_names <- header_values
    } else if (length(header_values) == length(sample_cols) + 1L && !nzchar(header_values[1])) {
      candidate_names <- header_values[-1]
    } else {
      candidate_names <- sample_cols
    }
    candidate_names <- vapply(candidate_names, normalize_token, character(1))
    candidate_names[!nzchar(candidate_names)] <- sample_cols[!nzchar(candidate_names)]
    setnames(dt, sample_cols, make.unique(candidate_names, sep = "_dup"))
  }
  dt[, feature_id := vapply(feature_id, normalize_token, character(1))]
  dt <- dt[nzchar(feature_id)]
  cat("[supp_fread_done]", label, "features", nrow(dt), "samples", ncol(dt) - 1L, "\n")
  dt
}

get_reader_config <- function(dataset, platform, source_url, source_type) {
  if (dataset == "GSE163882") {
    return(list(reader = "plain", sep = ",", header_has_feature_col = TRUE, annotation_col = "annotation"))
  }
  if (dataset == "GSE109710") {
    return(list(reader = "plain", sep = "\t", header_has_feature_col = FALSE, annotation_col = NULL))
  }
  if (dataset == "GSE194040") {
    return(list(reader = "plain", sep = "\t", header_has_feature_col = FALSE, annotation_col = NULL))
  }
  list(reader = "series", sep = "\t", header_has_feature_col = TRUE, annotation_col = NULL)
}

extract_genefu_object <- function(tarfile, rda_path) {
  tmp <- tempfile("genefu_src_")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(tarfile, files = rda_path, exdir = tmp)
  f <- file.path(tmp, rda_path)
  if (!file.exists(f)) stop("Could not extract ", rda_path)
  env <- new.env(parent = emptyenv())
  loaded <- load(f, envir = env)
  get(loaded[[1]], envir = env)
}

aggregate_weights_by_gene <- function(genes, weights, roles = NULL) {
  dt <- data.table(
    gene = apply_alias_map(genes),
    weight = suppressWarnings(as.numeric(weights)),
    role = if (is.null(roles)) NA_character_ else as.character(roles)
  )
  dt <- dt[!is.na(gene)]
  dt[, .(
    weight = if (all(is.na(weight))) NA_real_ else mean(weight, na.rm = TRUE),
    role = paste(sort(unique(role[!is.na(role) & nzchar(role)])), collapse = ";")
  ), by = gene]
}

zscore_rows <- function(mat) {
  if (!nrow(mat) || !ncol(mat)) return(mat)
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, stats::sd, na.rm = TRUE)
  keep <- is.finite(mu) & is.finite(sdv) & sdv > 0
  z <- mat[keep, , drop = FALSE]
  z <- sweep(z, 1, mu[keep], "-")
  z <- sweep(z, 1, sdv[keep], "/")
  z
}

weighted_score <- function(zmat, gene_weights) {
  genes <- intersect(rownames(zmat), gene_weights$gene)
  if (!length(genes)) return(rep(NA_real_, ncol(zmat)))
  w <- gene_weights[match(genes, gene)]$weight
  valid <- is.finite(w)
  genes <- genes[valid]
  w <- w[valid]
  if (!length(genes)) return(rep(NA_real_, ncol(zmat)))
  score <- as.numeric(crossprod(w, zmat[genes, , drop = FALSE])) / sum(abs(w))
  names(score) <- colnames(zmat)
  score
}

mean_score <- function(zmat, genes) {
  genes <- intersect(rownames(zmat), apply_alias_map(genes))
  if (!length(genes)) return(rep(NA_real_, ncol(zmat)))
  score <- colMeans(zmat[genes, , drop = FALSE], na.rm = TRUE)
  names(score) <- colnames(zmat)
  score
}

mean_diff_score <- function(zmat, up_genes, down_genes) {
  up <- intersect(rownames(zmat), apply_alias_map(up_genes))
  down <- intersect(rownames(zmat), apply_alias_map(down_genes))
  up_score <- if (length(up)) colMeans(zmat[up, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(zmat))
  down_score <- if (length(down)) colMeans(zmat[down, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(zmat))
  score <- up_score - down_score
  names(score) <- colnames(zmat)
  score
}

row_max_finite <- function(mat) {
  apply(mat, 1, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    max(x)
  })
}

safe_summary_stat <- function(x, fun) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  fun(x)
}

make_score_rows <- function(meta, classifier_id, classifier_score_name, score, n_required,
                            n_available, preferred_min_pct, absolute_min_pct,
                            score_direction, note, auxiliary_class = NA_character_) {
  coverage_pct <- if (n_required > 0) round(100 * n_available / n_required, 1) else NA_real_
  score_status <- ifelse(
    is.na(coverage_pct) || coverage_pct < absolute_min_pct,
    "not_scored_insufficient_coverage",
    ifelse(coverage_pct < preferred_min_pct, "scored_low_coverage_exploratory", "scored_preferred_coverage")
  )
  if (score_status[1] == "not_scored_insufficient_coverage") {
    score[] <- NA_real_
  }
  score[!is.finite(score)] <- NA_real_
  data.table(
    dataset_accession = meta$dataset_accession,
    matrix_label = meta$matrix_label,
    matrix_file = meta$matrix_file,
    platform = meta$platform,
    gsm_sample_id = meta$gsm_sample_id,
    patient_id = meta$patient_id,
    endpoint_standard = meta$endpoint_standard,
    endpoint_binary_pcr1_rd0 = meta$endpoint_binary_pcr1_rd0,
    er_status = meta$er_status,
    pr_status = meta$pr_status,
    her2_status = meta$her2_status,
    hr_status = meta$hr_status,
    subtype_group = meta$subtype_group,
    treatment_arm = meta$treatment_arm,
    planned_role = meta$planned_role,
    story_scope = meta$story_scope,
    include_primary_analysis = meta$include_primary_analysis,
    analysis_set = meta$analysis_set,
    classifier_id = classifier_id,
    classifier_score_name = classifier_score_name,
    score_value = as.numeric(score),
    score_direction_higher_means = score_direction,
    score_auxiliary_class = auxiliary_class,
    n_genes_required = n_required,
    n_genes_available = n_available,
    coverage_pct = coverage_pct,
    score_status = score_status,
    score_note = note
  )
}

input_files <- c(
  sample_manifest = file.path(step1_dir, "analysis_ready_manifest_v1.csv"),
  expression_audit = file.path(step2_dir, "expression_matrix_audit_v1.csv"),
  classifier_dictionary = file.path(step3_dir, "classifier_dictionary_v1.csv"),
  genefu_tarball = file.path(step3_dir, "genefu_2.44.0.tar.gz")
)
cat("\nInput file check:\n")
for (nm in names(input_files)) cat(nm, file.exists(input_files[[nm]]), input_files[[nm]], "\n")
if (!all(file.exists(input_files))) {
  stop("One or more required input files are missing.")
}

sample_manifest <- fread(input_files[["sample_manifest"]], encoding = "UTF-8", showProgress = FALSE)
expression_audit <- fread(input_files[["expression_audit"]], encoding = "UTF-8", showProgress = FALSE)
classifier_dictionary <- fread(input_files[["classifier_dictionary"]], encoding = "UTF-8", showProgress = FALSE)

cat("\nLoaded rows:\n")
cat("sample_manifest:", nrow(sample_manifest), "\n")
cat("expression_audit:", nrow(expression_audit), "\n")
cat("classifier_dictionary:", nrow(classifier_dictionary), "\n")

genefu_tarball <- input_files[["genefu_tarball"]]
pam50 <- extract_genefu_object(genefu_tarball, "genefu/data/pam50.rda")
sig_gene70 <- extract_genefu_object(genefu_tarball, "genefu/data/sig.gene70.rda")
sig_oncotypedx <- extract_genefu_object(genefu_tarball, "genefu/data/sig.oncotypedx.rda")
sig_ggi <- extract_genefu_object(genefu_tarball, "genefu/data/sig.ggi.rda")
sig_gene76 <- extract_genefu_object(genefu_tarball, "genefu/data/sig.gene76.rda")
sig_endo <- extract_genefu_object(genefu_tarball, "genefu/data/sig.endoPredict.rda")

pam50_genes <- apply_alias_map(pam50$centroids.map$probe)
pam50_centroids <- pam50$centroids
rownames(pam50_centroids) <- apply_alias_map(rownames(pam50_centroids))

gene70_weights <- aggregate_weights_by_gene(sig_gene70$HUGO.gene.symbol, -sig_gene70$correlation)
oncotype_weights <- aggregate_weights_by_gene(sig_oncotypedx$symbol, sig_oncotypedx$weight, sig_oncotypedx$group)
ggi_dt <- data.table(gene = apply_alias_map(sig_ggi$HUGO.gene.symbol), grade = sig_ggi$grade)
ggi_dt <- unique(ggi_dt[!is.na(gene) & grade %in% c(1, 3)])
gene76_weights <- aggregate_weights_by_gene(sig_gene76$HUGO.gene.symbol, sig_gene76$std.cox.coefficient, ifelse(sig_gene76$er == 1, "ER_positive_component", "ER_negative_component"))
endo_weights <- aggregate_weights_by_gene(sig_endo$symbol, sig_endo$weight, sig_endo$group)
ifng_genes <- apply_alias_map(c("STAT1", "IDO1", "CXCL10", "CXCL9", "HLA-DRA", "GZMB", "IFNG", "LAG3", "CD8A",
                                "CCL5", "CXCL13", "HLA-E", "NKG7", "TIGIT", "PSMB10", "CD27", "CD274", "PDCD1LG2"))
cyto_genes <- apply_alias_map(c("GZMA", "PRF1"))

scoring_spec <- data.table(
  classifier_id = c("PAM50", "MAMMAPRINT_70", "ONCOTYPE_DX_21", "GGI_128_PROBE", "GENE76", "ENDOPREDICT_11", "IFNG_18", "CYTOLYTIC_ACTIVITY_2"),
  classifier_score_name = c(
    "PAM50 non-luminal affinity score",
    "MammaPrint approximate poor-prognosis score",
    "OncotypeDX approximate weighted expression score",
    "GGI grade-3 minus grade-1 expression score",
    "Wang 76-gene weighted risk score",
    "EndoPredict approximate weighted expression score",
    "IFN-gamma 18-gene mean z score",
    "Cytolytic activity mean z score"
  ),
  n_required_genes = c(
    length(unique(pam50_genes[!is.na(pam50_genes)])),
    uniqueN(gene70_weights$gene),
    uniqueN(oncotype_weights[is.finite(weight)]$gene),
    uniqueN(ggi_dt$gene),
    uniqueN(gene76_weights[is.finite(weight)]$gene),
    uniqueN(endo_weights[is.finite(weight)]$gene),
    length(unique(ifng_genes[!is.na(ifng_genes)])),
    length(unique(cyto_genes[!is.na(cyto_genes)]))
  ),
  preferred_min_coverage_pct = c(80, 80, 80, 80, 80, 80, 80, 100),
  absolute_min_coverage_pct = c(50, 50, 50, 50, 50, 50, 50, 100),
  score_direction_higher_means = c(
    "greater basal/HER2-enriched than luminal/normal centroid affinity; approximate pCR-favorable biology",
    "higher approximate poor-prognosis expression pattern; not a pCR-specific direction",
    "higher approximate recurrence/proliferation weighted expression; not the commercial recurrence score",
    "higher histologic-grade/proliferation-like expression",
    "higher distant-metastasis risk-like expression; not pCR-specific",
    "higher recurrence-risk-like expression in ER-positive/HER2-negative context",
    "higher interferon-gamma immune activation",
    "higher local cytolytic activity"
  ),
  scoring_rule = c(
    "Spearman correlation to PAM50 centroids using available genes; score=max(Basal,Her2)-max(LumA,LumB,Normal); subtype=max centroid.",
    "Weighted z-score using negative genefu gene70 correlation weights; clinical MammaPrint cutoff is not reconstructed.",
    "Weighted z-score using non-reference genefu OncotypeDX weights; commercial algorithm is not reconstructed.",
    "Mean z-score of grade 3 genes minus mean z-score of grade 1 genes from genefu GGI object.",
    "Weighted z-score using genefu Wang 76-gene Cox coefficients; ER-specific clinical rule is not reconstructed.",
    "Weighted z-score using genefu EndoPredict GOI weights; clinical calibration is not reconstructed.",
    "Mean z-score of the 18 IFN-gamma profile genes.",
    "Mean z-score of GZMA and PRF1; requires both genes."
  ),
  main_analysis_use = c(
    "eligible_if_preferred_coverage_met",
    "eligible_if_preferred_coverage_met_else_exploratory",
    "eligible_if_preferred_coverage_met",
    "eligible_if_preferred_coverage_met",
    "eligible_if_preferred_coverage_met",
    "secondary_HR_positive_context",
    "secondary_mechanistic_score",
    "secondary_mechanistic_score"
  )
)

fwrite(scoring_spec, file.path(step4_dir, "classifier_scoring_spec_v1.csv"))
cat("Wrote classifier_scoring_spec_v1.csv rows:", nrow(scoring_spec), "\n")

protocol <- c(
  "# Step 4 outcome-blind preprocessing and classifier scoring protocol v1",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Fixed analysis boundary",
  "",
  "- This step does not build a new pCR signature and does not select genes using pCR/RD.",
  "- The purpose is to convert each public expression matrix into deterministic published-classifier or mechanistic scores.",
  "- GSE194040 platforms are kept separate.",
  "- GSE106977 and GSE109710 are retained as limited gene-space/panel stress checks but should not drive main conclusions.",
  "",
  "## Expression preprocessing",
  "",
  "1. Expression matrices are read from the public URLs already audited in Step 2.",
  "2. Manifest sample matching uses GSM/sample title/patient aliases from Step 1 and the same normalization rule as Step 2.",
  "3. Features are mapped to Step 2 feature-index gene symbols.",
  "4. Probes mapping to multiple genes are excluded from scoring to avoid ambiguous feature attribution.",
  "5. Duplicate single-gene probes are collapsed by the arithmetic mean across probes for each sample. This is deterministic and outcome-blind.",
  "6. Matrices labelled as likely unlogged intensity/count/TPM, or with sampled maxima >50, are transformed as log2(x+1) after truncating negative values to zero. Matrices already on a log-like continuous scale are kept as-is.",
  "7. Gene expression is standardized within each matrix across matched manifest samples: z = (expression - gene mean) / gene SD.",
  "8. Genes with zero or non-finite SD are dropped before scoring.",
  "",
  "## Leakage control",
  "",
  "- pCR/RD, receptor status, treatment arm, and planned-role labels are never used to choose preprocessing parameters.",
  "- Coverage thresholds are pre-specified in classifier_scoring_spec_v1.csv.",
  "- Low-coverage scores are flagged as exploratory or not scored; they are not silently mixed into the main benchmark.",
  "",
  "## Scoring outputs",
  "",
  "- `classifier_score_matrix_v1.csv`: long sample-by-classifier score table.",
  "- `pam50_subtype_assignments_v1.csv`: PAM50 centroid correlations and subtype calls.",
  "- `classifier_score_qc_v1.csv`: per-dataset per-classifier coverage and scoring status.",
  "- `dataset_preprocessing_decisions_v1.csv`: scale transform, matched samples, gene-level matrix dimensions.",
  "- `sample_score_match_qc_v1.csv`: manifest-to-expression sample matching status."
)
writeLines(protocol, file.path(step4_dir, "preprocessing_protocol_v1.md"), useBytes = TRUE)

dataset_decisions <- list()
sample_match_qc <- list()
matrix_qc <- list()
score_rows <- list()
pam50_rows <- list()
errors <- list()

append_list <- function(object_name, value) {
  current <- get(object_name, envir = .GlobalEnv)
  current[[length(current) + 1L]] <- value
  assign(object_name, current, envir = .GlobalEnv)
  invisible(TRUE)
}

for (i in seq_len(nrow(expression_audit))) {
  ea <- expression_audit[i]
  dataset <- ea$dataset_accession
  platform <- ea$platform
  matrix_label <- ea$matrix_label
  matrix_file <- ea$matrix_file
  matrix_file_value <- ea$matrix_file
  source_url <- ea$expression_matrix_source_url
  source_type <- ea$expression_source_type
  cat("\n[dataset_start]", matrix_label, dataset, platform, "\n")
  tryCatch({
    cfg <- get_reader_config(dataset, platform, source_url, source_type)
    expr_dt <- if (identical(cfg$reader, "series")) {
      read_series_matrix_values(source_url, matrix_label)
    } else {
      read_plain_expression_values(
        source_url, matrix_label,
        sep = cfg$sep,
        header_has_feature_col = cfg$header_has_feature_col,
        annotation_col = cfg$annotation_col
      )
    }
    feature_index_file <- ea$feature_index_file
    if (!file.exists(feature_index_file)) stop("Missing feature index: ", feature_index_file)
    feature_index <- fread(feature_index_file, encoding = "UTF-8", showProgress = FALSE)
    feature_index[, gene_symbol_single := ifelse(grepl(";", gene_symbol, fixed = TRUE), NA_character_, apply_alias_map(gene_symbol))]

    sample_cols <- setdiff(names(expr_dt), "feature_id")
    manifest_rows <- sample_manifest[dataset_accession == dataset & matrix_file == matrix_file_value]
    norm_sample_cols <- vapply(sample_cols, normalize_matrix_sample_id, character(1))
    sample_lookup <- data.table(matrix_sample_id = sample_cols, matrix_sample_norm = norm_sample_cols)
    sample_lookup <- sample_lookup[nzchar(matrix_sample_norm)]
    sample_lookup <- sample_lookup[!duplicated(matrix_sample_norm)]

    match_rows <- list()
    for (j in seq_len(nrow(manifest_rows))) {
      aliases <- unique(vapply(sample_aliases_for_row(manifest_rows[j]), normalize_matrix_sample_id, character(1)))
      aliases <- aliases[nzchar(aliases)]
      hit <- sample_lookup[matrix_sample_norm %in% aliases]
      status <- if (nrow(hit)) "matched" else "not_matched"
      matrix_sample_id <- if (nrow(hit)) hit$matrix_sample_id[1] else NA_character_
      match_rows[[j]] <- data.table(
        dataset_accession = dataset,
        matrix_label = matrix_label,
        matrix_file = matrix_file,
        platform = platform,
        gsm_sample_id = manifest_rows$gsm_sample_id[j],
        patient_id = manifest_rows$patient_id[j],
        matrix_sample_id = matrix_sample_id,
        match_status = status
      )
    }
    match_dt <- rbindlist(match_rows, fill = TRUE)
    append_list("sample_match_qc", match_dt)
    matched_dt <- match_dt[match_status == "matched"]
    if (!nrow(matched_dt)) stop("No matched samples for ", matrix_label)

    keep_cols <- matched_dt$matrix_sample_id
    value_dt <- expr_dt[, c("feature_id", keep_cols), with = FALSE]
    for (cc in keep_cols) {
      value_dt[, (cc) := suppressWarnings(as.numeric(get(cc)))]
    }
    mat <- as.matrix(value_dt[, ..keep_cols])
    rownames(mat) <- value_dt$feature_id
    storage.mode(mat) <- "double"
    raw_min <- suppressWarnings(min(mat, na.rm = TRUE))
    raw_max <- suppressWarnings(max(mat, na.rm = TRUE))
    if (!is.finite(raw_min)) raw_min <- NA_real_
    if (!is.finite(raw_max)) raw_max <- NA_real_

    scale_guess <- ea$expression_scale_guess
    transform_rule <- "keep_as_provided"
    if (grepl("linear|unlogged|count", scale_guess, ignore.case = TRUE) || (!is.na(raw_max) && raw_max > 50)) {
      mat <- log2(pmax(mat, 0) + 1)
      transform_rule <- "log2_pmax0_plus1"
    }
    transformed_min <- suppressWarnings(min(mat, na.rm = TRUE))
    transformed_max <- suppressWarnings(max(mat, na.rm = TRUE))

    map <- feature_index[match(rownames(mat), feature_id)]
    gene <- map$gene_symbol_single
    keep_feature <- !is.na(gene) & nzchar(gene)
    mat_single <- mat[keep_feature, , drop = FALSE]
    gene_single <- gene[keep_feature]
    if (!nrow(mat_single)) stop("No single-gene features available for ", matrix_label)

    na_mask <- is.na(mat_single)
    mat_zero <- mat_single
    mat_zero[na_mask] <- 0
    sums <- rowsum(mat_zero, group = gene_single, reorder = FALSE)
    counts <- rowsum((!na_mask) + 0, group = gene_single, reorder = FALSE)
    gene_mat <- sums / counts
    gene_mat[counts == 0] <- NA_real_
    rownames(gene_mat) <- rownames(sums)
    colnames(gene_mat) <- keep_cols
    zmat <- zscore_rows(gene_mat)

    analysis_set <- if (ea$usable_for_full_benchmark == "yes") {
      "core_full_benchmark"
    } else if (ea$usable_for_full_benchmark == "stress_test_only") {
      "large_stress_test"
    } else if (grepl("limited", ea$usable_for_subtype_only_validation, ignore.case = TRUE)) {
      "limited_subtype_panel_or_gene_space"
    } else {
      "subtype_validation_only"
    }

    meta <- merge(
      matched_dt,
      manifest_rows,
      by = c("dataset_accession", "matrix_file", "platform", "gsm_sample_id", "patient_id"),
      all.x = TRUE,
      sort = FALSE
    )
    meta[, matrix_label := matrix_label]
    meta[, analysis_set := analysis_set]
    meta[, sample_order := match(matrix_sample_id, keep_cols)]
    setorder(meta, sample_order)
    meta[, sample_order := NULL]

    append_list("dataset_decisions", data.table(
      dataset_accession = dataset,
      matrix_label = matrix_label,
      matrix_file = matrix_file,
      platform = platform,
      expression_source_type = source_type,
      expression_matrix_source_url = source_url,
      n_manifest_samples = nrow(manifest_rows),
      n_matched_samples = nrow(matched_dt),
      n_features_raw = nrow(expr_dt),
      n_single_gene_features_used = nrow(mat_single),
      n_unique_genes_after_collapse = nrow(gene_mat),
      n_genes_after_zscore_nonzero_sd = nrow(zmat),
      raw_min = raw_min,
      raw_max = raw_max,
      expression_scale_guess_step2 = scale_guess,
      transform_rule = transform_rule,
      transformed_min = transformed_min,
      transformed_max = transformed_max,
      duplicate_probe_collapse_rule = "mean_across_single_gene_probes_outcome_blind",
      multi_gene_probe_rule = "excluded_from_scoring",
      analysis_set = analysis_set
    ))

    add_qc <- function(classifier_id, n_required, n_available, status_detail = "") {
      append_list("matrix_qc", data.table(
        dataset_accession = dataset,
        matrix_label = matrix_label,
        matrix_file = matrix_file,
        platform = platform,
        analysis_set = analysis_set,
        classifier_id = classifier_id,
        n_genes_required = n_required,
        n_genes_available = n_available,
        coverage_pct = if (n_required > 0) round(100 * n_available / n_required, 1) else NA_real_,
        status_detail = status_detail
      ))
    }

    # PAM50 centroid classifier and non-luminal affinity score
    pam_genes_avail <- intersect(rownames(zmat), rownames(pam50_centroids))
    pam_required <- length(unique(rownames(pam50_centroids)))
    pam_available <- length(pam_genes_avail)
    add_qc("PAM50", pam_required, pam_available)
    pam_scores <- rep(NA_real_, ncol(zmat))
    pam_class <- rep(NA_character_, ncol(zmat))
    names(pam_scores) <- colnames(zmat)
    pam_corr_dt <- NULL
    if (pam_available >= ceiling(0.5 * pam_required)) {
      cent <- pam50_centroids[pam_genes_avail, , drop = FALSE]
      expr <- zmat[pam_genes_avail, , drop = FALSE]
      corr_mat <- matrix(NA_real_, nrow = ncol(expr), ncol = ncol(cent))
      rownames(corr_mat) <- colnames(expr)
      colnames(corr_mat) <- colnames(cent)
      for (s in colnames(expr)) {
        for (ct in colnames(cent)) {
          corr_mat[s, ct] <- suppressWarnings(cor(expr[, s], cent[, ct], method = "spearman", use = "pairwise.complete.obs"))
        }
      }
      pam_class <- colnames(corr_mat)[max.col(corr_mat, ties.method = "first")]
      basal_cols <- grep("basal", colnames(corr_mat), ignore.case = TRUE, value = TRUE)
      her2_cols <- grep("her2", colnames(corr_mat), ignore.case = TRUE, value = TRUE)
      luma_cols <- grep("luma", colnames(corr_mat), ignore.case = TRUE, value = TRUE)
      lumb_cols <- grep("lumb", colnames(corr_mat), ignore.case = TRUE, value = TRUE)
      normal_cols <- grep("normal", colnames(corr_mat), ignore.case = TRUE, value = TRUE)
      nonlum <- row_max_finite(corr_mat[, c(basal_cols, her2_cols), drop = FALSE])
      lum <- row_max_finite(corr_mat[, c(luma_cols, lumb_cols, normal_cols), drop = FALSE])
      pam_scores <- nonlum - lum
      pam_scores[!is.finite(pam_scores)] <- NA_real_
      pam_corr_dt <- as.data.table(corr_mat, keep.rownames = "matrix_sample_id")
      pam_corr_dt <- cbind(meta[, .(dataset_accession, matrix_label, matrix_file, platform, gsm_sample_id, patient_id, analysis_set)], pam_corr_dt)
      pam_corr_dt[, pam50_subtype := pam_class]
      pam_corr_dt[, n_genes_required := pam_required]
      pam_corr_dt[, n_genes_available := pam_available]
      pam_corr_dt[, coverage_pct := round(100 * pam_available / pam_required, 1)]
      append_list("pam50_rows", pam_corr_dt)
    }
    append_list("score_rows", make_score_rows(
      meta, "PAM50", "PAM50 non-luminal affinity score", pam_scores,
      pam_required, pam_available, 80, 50,
      "higher basal/HER2-enriched centroid affinity relative to luminal/normal",
      "Approximate transferability score derived from PAM50 centroids; not an original pCR model.",
      pam_class
    ))

    # MammaPrint/gene70 approximate poor-prognosis score
    req <- uniqueN(gene70_weights$gene)
    avail <- length(intersect(rownames(zmat), gene70_weights$gene))
    add_qc("MAMMAPRINT_70", req, avail)
    append_list("score_rows", make_score_rows(
      meta, "MAMMAPRINT_70", "MammaPrint approximate poor-prognosis score",
      weighted_score(zmat, gene70_weights), req, avail, 80, 50,
      "higher approximate poor-prognosis expression pattern",
      "Clinical MammaPrint classification/cutoff is not reconstructed; use as approximate expression score."
    ))

    # OncotypeDX approximate weighted score
    onc_w <- oncotype_weights[is.finite(weight)]
    req <- uniqueN(onc_w$gene)
    avail <- length(intersect(rownames(zmat), onc_w$gene))
    add_qc("ONCOTYPE_DX_21", req, avail)
    append_list("score_rows", make_score_rows(
      meta, "ONCOTYPE_DX_21", "OncotypeDX approximate weighted expression score",
      weighted_score(zmat, onc_w), req, avail, 80, 50,
      "higher approximate recurrence/proliferation weighted expression",
      "Commercial 21-gene recurrence score algorithm is not reconstructed."
    ))

    # GGI grade score
    ggi_up <- ggi_dt[grade == 3, unique(gene)]
    ggi_down <- ggi_dt[grade == 1, unique(gene)]
    req <- length(unique(c(ggi_up, ggi_down)))
    avail <- length(intersect(rownames(zmat), unique(c(ggi_up, ggi_down))))
    add_qc("GGI_128_PROBE", req, avail)
    append_list("score_rows", make_score_rows(
      meta, "GGI_128_PROBE", "GGI grade-3 minus grade-1 expression score",
      mean_diff_score(zmat, ggi_up, ggi_down), req, avail, 80, 50,
      "higher grade/proliferation-like expression",
      "Computed from genefu GGI grade 3 versus grade 1 gene groups."
    ))

    # Gene76 weighted score
    g76_w <- gene76_weights[is.finite(weight)]
    req <- uniqueN(g76_w$gene)
    avail <- length(intersect(rownames(zmat), g76_w$gene))
    add_qc("GENE76", req, avail)
    append_list("score_rows", make_score_rows(
      meta, "GENE76", "Wang 76-gene weighted risk score",
      weighted_score(zmat, g76_w), req, avail, 80, 50,
      "higher distant-metastasis risk-like expression",
      "ER-specific clinical rule is not reconstructed; score is an approximate weighted expression benchmark."
    ))

    # EndoPredict approximate score
    endo_w <- endo_weights[is.finite(weight)]
    req <- uniqueN(endo_w$gene)
    avail <- length(intersect(rownames(zmat), endo_w$gene))
    add_qc("ENDOPREDICT_11", req, avail)
    append_list("score_rows", make_score_rows(
      meta, "ENDOPREDICT_11", "EndoPredict approximate weighted expression score",
      weighted_score(zmat, endo_w), req, avail, 80, 50,
      "higher recurrence-risk-like expression in ER-positive/HER2-negative context",
      "Clinical EndoPredict calibration is not reconstructed."
    ))

    # IFNG 18
    req <- length(unique(ifng_genes[!is.na(ifng_genes)]))
    avail <- length(intersect(rownames(zmat), ifng_genes))
    add_qc("IFNG_18", req, avail)
    append_list("score_rows", make_score_rows(
      meta, "IFNG_18", "IFN-gamma 18-gene mean z score",
      mean_score(zmat, ifng_genes), req, avail, 80, 50,
      "higher interferon-gamma immune activation",
      "Secondary mechanistic immune score."
    ))

    # Cytolytic activity, require both genes
    req <- length(unique(cyto_genes[!is.na(cyto_genes)]))
    avail <- length(intersect(rownames(zmat), cyto_genes))
    add_qc("CYTOLYTIC_ACTIVITY_2", req, avail)
    append_list("score_rows", make_score_rows(
      meta, "CYTOLYTIC_ACTIVITY_2", "Cytolytic activity mean z score",
      mean_score(zmat, cyto_genes), req, avail, 100, 100,
      "higher local cytolytic activity",
      "Secondary two-gene immune score; both genes required."
    ))

    rm(expr_dt, value_dt, mat, mat_single, gene_mat, zmat)
    gc(verbose = FALSE)
    cat("[dataset_done]", matrix_label, "\n")
  }, error = function(e) {
    cat("[dataset_error]", matrix_label, conditionMessage(e), "\n")
    append_list("errors", data.table(
      dataset_accession = dataset,
      matrix_label = matrix_label,
      matrix_file = matrix_file,
      platform = platform,
      error = conditionMessage(e)
    ))
  })
}

dataset_decisions_dt <- if (length(dataset_decisions)) rbindlist(dataset_decisions, fill = TRUE) else data.table()
sample_match_qc_dt <- if (length(sample_match_qc)) rbindlist(sample_match_qc, fill = TRUE) else data.table()
matrix_qc_dt <- if (length(matrix_qc)) rbindlist(matrix_qc, fill = TRUE) else data.table()
score_dt <- if (length(score_rows)) rbindlist(score_rows, fill = TRUE) else data.table()
pam50_dt <- if (length(pam50_rows)) rbindlist(pam50_rows, fill = TRUE) else data.table()
errors_dt <- if (length(errors)) rbindlist(errors, fill = TRUE) else data.table(dataset_accession = character(), matrix_label = character(), matrix_file = character(), platform = character(), error = character())

if (nrow(score_dt)) {
  score_dt[, endpoint_binary_pcr1_rd0 := suppressWarnings(as.integer(endpoint_binary_pcr1_rd0))]
}

score_qc <- if (nrow(score_dt)) {
  score_dt[, .(
    n_samples = .N,
    n_scores_nonmissing = sum(is.finite(score_value)),
    n_scores_missing = sum(!is.finite(score_value)),
    n_unique_patients = uniqueN(patient_id),
    score_mean = safe_summary_stat(score_value, mean),
    score_sd = safe_summary_stat(score_value, stats::sd),
    score_min = safe_summary_stat(score_value, min),
    score_median = safe_summary_stat(score_value, median),
    score_max = safe_summary_stat(score_value, max),
    score_status_values = paste(sort(unique(score_status)), collapse = ";")
  ), by = .(dataset_accession, matrix_label, platform, analysis_set, classifier_id, classifier_score_name, n_genes_required, n_genes_available, coverage_pct)]
} else {
  data.table()
}
if (nrow(score_qc)) {
  for (cc in c("score_mean", "score_sd", "score_min", "score_median", "score_max")) {
    score_qc[!is.finite(get(cc)), (cc) := NA_real_]
  }
}

fwrite(dataset_decisions_dt, file.path(step4_dir, "dataset_preprocessing_decisions_v1.csv"))
fwrite(sample_match_qc_dt, file.path(step4_dir, "sample_score_match_qc_v1.csv"))
fwrite(matrix_qc_dt, file.path(step4_dir, "gene_level_matrix_qc_v1.csv"))
fwrite(score_dt, file.path(step4_dir, "classifier_score_matrix_v1.csv"))
fwrite(score_qc, file.path(step4_dir, "classifier_score_qc_v1.csv"))
fwrite(pam50_dt, file.path(step4_dir, "pam50_subtype_assignments_v1.csv"))
fwrite(errors_dt, file.path(step4_dir, "classifier_score_errors_v1.csv"))

cat("\nOutput rows:\n")
cat("dataset_preprocessing_decisions_v1.csv:", nrow(dataset_decisions_dt), "\n")
cat("sample_score_match_qc_v1.csv:", nrow(sample_match_qc_dt), "\n")
cat("gene_level_matrix_qc_v1.csv:", nrow(matrix_qc_dt), "\n")
cat("classifier_score_matrix_v1.csv:", nrow(score_dt), "\n")
cat("classifier_score_qc_v1.csv:", nrow(score_qc), "\n")
cat("pam50_subtype_assignments_v1.csv:", nrow(pam50_dt), "\n")
cat("classifier_score_errors_v1.csv:", nrow(errors_dt), "\n")

notes <- c(
  "# Step 4 classifier scoring execution notes v1",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Completed outputs",
  "",
  "- preprocessing_protocol_v1.md",
  "- classifier_scoring_spec_v1.csv",
  "- dataset_preprocessing_decisions_v1.csv",
  "- sample_score_match_qc_v1.csv",
  "- gene_level_matrix_qc_v1.csv",
  "- classifier_score_matrix_v1.csv",
  "- classifier_score_qc_v1.csv",
  "- pam50_subtype_assignments_v1.csv",
  "- classifier_score_errors_v1.csv",
  "- build_classifier_scores_step4_v1.R",
  "- R_raw_run_log_step4_v1.txt",
  "- R_session_info_step4_v1.txt",
  "",
  "## Important interpretation",
  "",
  "These scores are pre-specified transferability benchmarks, not newly trained pCR predictors. Commercial or clinical classifiers are represented only by transparent gene-level approximate scores when the exact clinical algorithm is not public.",
  "",
  "## Error summary",
  "",
  if (nrow(errors_dt)) paste(capture.output(print(errors_dt)), collapse = "\n") else "No dataset-level scoring errors recorded.",
  "",
  "## Score QC summary",
  "",
  if (nrow(score_qc)) paste(capture.output(print(score_qc[, .(dataset_accession, matrix_label, classifier_id, n_scores_nonmissing, coverage_pct, score_status_values)])), collapse = "\n") else "No score QC rows produced."
)
writeLines(notes, file.path(step4_dir, "step4_classifier_scoring_notes_v1.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(step4_dir, "R_session_info_step4_v1.txt"))

cat("\nSession info written.\n")
if (nrow(errors_dt)) {
  cat("\nErrors were recorded; inspect classifier_score_errors_v1.csv.\n")
} else {
  cat("\nDone without dataset-level scoring errors.\n")
}
