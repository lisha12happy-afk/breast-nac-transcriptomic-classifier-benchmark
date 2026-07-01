options(stringsAsFactors = FALSE)

out_dir <- file.path(getwd(), "manifests", "final_article_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "qc_user_arranged_figures_tables_v1.txt")
sink(log_file, split = TRUE)

cat("QC user-arranged manuscript figure/table files\n")
cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("readxl is not installed")
}

xlsx <- "manifests/final_article_tables/supplementary table/Supplementary Table S4/Supplementary_Table_S4_robustness_subtype_calibration_failure_modes_final_v1.xlsx"
cat("Supplementary Tables path:", xlsx, "\n")
cat("Exists:", file.exists(xlsx), "\n")

sheet_names <- readxl::excel_sheets(xlsx)
cat("Sheets:", paste(sheet_names, collapse = " | "), "\n")

rows <- lapply(sheet_names, function(s) {
  x <- readxl::read_excel(xlsx, sheet = s)
  empty_cols <- names(x)[vapply(x, function(z) all(is.na(z) | z == ""), logical(1))]
  data.frame(
    sheet = s,
    rows = nrow(x),
    cols = ncol(x),
    empty_cols_n = length(empty_cols),
    empty_cols = paste(empty_cols, collapse = "; "),
    headers = paste(names(x), collapse = " || "),
    stringsAsFactors = FALSE
  )
})
summary_tbl <- do.call(rbind, rows)
print(summary_tbl[, c("sheet", "rows", "cols", "empty_cols_n")], row.names = FALSE)

write.csv(summary_tbl, file.path(out_dir, "qc_user_supplementary_tables_summary_v1.csv"), row.names = FALSE)
cat("Wrote QC summary CSV.\n")

sink()
