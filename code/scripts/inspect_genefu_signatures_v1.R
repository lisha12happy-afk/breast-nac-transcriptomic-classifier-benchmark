args <- commandArgs(trailingOnly = TRUE)
step3_dir <- if (length(args) >= 1) args[[1]] else file.path(getwd(), "manifests", "step 3")
tarfile <- file.path(step3_dir, "genefu_2.44.0.tar.gz")

cat("inspect_genefu_signatures_v1.R\n")
cat("Working directory:", getwd(), "\n")
cat("Step 3 directory:", step3_dir, "\n")
cat("Tarfile:", tarfile, "\n")
cat("Tarfile exists:", file.exists(tarfile), "\n")
cat("R version:", R.version.string, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

if (!file.exists(tarfile)) {
  stop("Missing genefu source tarball: ", tarfile)
}

tmp <- tempfile("genefu_src_")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

targets <- c(
  "genefu/data/pam50.rda",
  "genefu/data/pam50.robust.rda",
  "genefu/data/sig.gene70.rda",
  "genefu/data/sig.oncotypedx.rda",
  "genefu/data/sig.ggi.rda",
  "genefu/data/sig.gene76.rda",
  "genefu/data/sig.endoPredict.rda",
  "genefu/data/sig.tamr13.rda"
)

utils::untar(tarfile, files = targets, exdir = tmp)

for (target in targets) {
  f <- file.path(tmp, target)
  cat("\n===== ", basename(target), " =====\n", sep = "")
  if (!file.exists(f)) {
    cat("MISSING\n")
    next
  }
  env <- new.env(parent = emptyenv())
  loaded <- load(f, envir = env)
  cat("Loaded objects:", paste(loaded, collapse = ", "), "\n")
  for (obj_name in loaded) {
    obj <- get(obj_name, envir = env)
    cat("Object:", obj_name, "\n")
    cat("Class:", paste(class(obj), collapse = ", "), "\n")
    cat("Dimensions:", paste(dim(obj), collapse = " x "), "\n")
    if (is.data.frame(obj) || is.matrix(obj)) {
      cat("Column names:", paste(colnames(obj), collapse = ", "), "\n")
      print(utils::head(obj, 10))
    } else {
      utils::str(obj, max.level = 2)
    }
  }
}

cat("\nSession info:\n")
print(sessionInfo())
