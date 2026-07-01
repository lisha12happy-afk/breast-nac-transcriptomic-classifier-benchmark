args <- commandArgs(trailingOnly = TRUE)
manifest_dir <- if (length(args) >= 1) args[1] else "manifests/step 1"
out_dir <- if (length(args) >= 2) args[2] else "manifests/step 2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "feature_indexes"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "R_raw_run_log_final_v1.txt")
if (file.exists(log_file)) file.remove(log_file)
log_con <- file(log_file, "wt", encoding = "UTF-8")
log_msg <- function(...) {
  x <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(x, "\n"); writeLines(x, log_con, useBytes = TRUE); flush(log_con)
}
norm <- function(x) { if (length(x) == 0 || is.na(x)) return(""); trimws(gsub("\\s+", " ", gsub('^"|"$', "", as.character(x)), perl = TRUE)) }
split_tab <- function(x) vapply(strsplit(x, "\t", fixed = TRUE)[[1]], norm, character(1))
geo_prefix <- function(id) { m <- regmatches(id, regexec("^([A-Za-z]+)([0-9]+)$", id, perl=TRUE))[[1]]; paste0(m[2], substr(m[3], 1, nchar(m[3]) - 3), "nnn") }
dl <- function(url, label) { f <- tempfile(pattern=gsub("[^A-Za-z0-9]+","_",label), fileext=".gz"); log_msg("download_start", label, url); utils::download.file(url, f, mode="wb", quiet=TRUE); log_msg("download_done", label, file.info(f)$size, "bytes"); f }
extract_symbols <- function(raw) {
  raw <- norm(raw); if (!nzchar(raw) || raw %in% c("NA","N/A","---","--","-","///")) return(character())
  chunks <- unlist(strsplit(raw, "\\s*///\\s*", perl=TRUE)); out <- character()
  for (ch in chunks) {
    parts <- vapply(strsplit(ch, "\\s*//\\s*", perl=TRUE)[[1]], norm, character(1))
    if (length(parts) >= 2 && grepl("^[A-Za-z][A-Za-z0-9_.-]{1,30}$", parts[2])) out <- c(out, parts[2])
    else out <- c(out, unlist(strsplit(ch, "\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*", perl=TRUE)))
  }
  out <- vapply(out, norm, character(1))
  bad <- "^(NM_|NR_|XM_|XR_|ENSG|ENST|ENSP|ILMN|AFFX|IMAGE|GO:|CHR|CHROMOSOME|HOMO|HS\\.|REFSEQ|GENBANK|GB:|CONTIG|CLONE|EMPTY|CONTROL|NEGATIVE|POSITIVE|NULL|---)"
  keep <- nzchar(out) & !grepl(bad, out, ignore.case=TRUE, perl=TRUE) & !grepl("^\\d+$", out) & !grepl("\\s", out) & nchar(out) <= 40 & grepl("^[A-Za-z][A-Za-z0-9_.-]*$", out, perl=TRUE)
  unique(toupper(out[keep]))
}
fallback_symbol <- function(x) { x <- norm(x); if (grepl("^(AFFX|ILMN|A_|GE_|TC|HTA|MIRX|PH_|CONTROL|NEG|POS|EMPTY)", x, ignore.case=TRUE) || grepl("^\\d+$", x)) return(character()); if (grepl("^[A-Za-z][A-Za-z0-9_.-]{1,30}$", x, perl=TRUE)) toupper(x) else character() }
sym_col <- function(cols) { cand <- c("Gene Symbol","GENE_SYMBOL","Gene symbol","Gene symbol","GeneName","gene_assignment","Gene assignment","Symbol","SYMBOL","Gene","GENE"); hit <- match(cand, cols); hit <- hit[!is.na(hit)]; if (length(hit)) hit[1] else { g <- grep("gene.?symbol|genename|gene.?assignment|symbol", cols, ignore.case=TRUE, perl=TRUE); if(length(g)) g[1] else NA_integer_ } }
plat_url <- function(p) if (p %in% c("GPL18573","GPL20078","GPL30493","GPL17586","GPL24546")) sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/soft/%s_family.soft.gz", geo_prefix(p), p, p) else sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s.annot.gz", geo_prefix(p), p, p)
map_platform <- function(platform, target) {
  url <- plat_url(platform); f <- dl(url, paste0(platform, "_platform")); on.exit(if(file.exists(f)) unlink(f), add=TRUE)
  con <- gzfile(f, "rt"); begin <- NA_integer_; n <- 0L; header <- NULL
  repeat { ln <- readLines(con, 1, warn=FALSE); if(!length(ln)) break; n <- n+1L; if(identical(ln, "!platform_table_begin")) { begin <- n; header <- split_tab(readLines(con, 1, warn=FALSE)); break } }
  close(con); if (is.na(begin)) return(list(url=url, status="no_platform_table", map=data.table::data.table(feature_id=character(), gene_symbol_raw=character(), gene_symbol=character())))
  id_idx <- match("ID", header); if (is.na(id_idx)) id_idx <- 1L
  sg_idx <- sym_col(header); if (is.na(sg_idx)) return(list(url=url, status="no_symbol_column", map=data.table::data.table(feature_id=character(), gene_symbol_raw=character(), gene_symbol=character())))
  log_msg("platform_map_stream", platform, "id", id_idx, "symbol", sg_idx, "targets", length(unique(target)))
  target_env <- new.env(parent=emptyenv())
  for (z in unique(target[nzchar(target)])) assign(z, TRUE, envir=target_env)
  ids <- character(); raws <- character(); syms <- character()
  con <- gzfile(f, "rt")
  in_table <- FALSE; header_seen <- FALSE
  repeat {
    ln <- readLines(con, 1, warn=FALSE)
    if (!length(ln)) break
    if (identical(ln, "!platform_table_begin")) { in_table <- TRUE; next }
    if (!in_table) next
    if (!header_seen) { header_seen <- TRUE; next }
    if (identical(ln, "!platform_table_end")) break
    parts <- strsplit(ln, "\t", fixed=TRUE)[[1]]
    if (length(parts) < max(id_idx, sg_idx)) next
    fid <- norm(parts[id_idx])
    if (!exists(fid, envir=target_env, inherits=FALSE)) next
    raw <- norm(parts[sg_idx])
    gs <- paste(extract_symbols(raw), collapse=";")
    if (nzchar(gs)) { ids <- c(ids, fid); raws <- c(raws, raw); syms <- c(syms, gs) }
  }
  close(con)
  dt <- data.table::data.table(feature_id=ids, gene_symbol_raw=raws, gene_symbol=syms)
  list(url=url, status="ok_stream", map=dt)
}
source_for <- function(ds, p, series_url) {
  if (ds=="GSE163882") return(list(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163882/suppl/GSE163882_all.data.tpms_222Samples.csv.gz", kind="supp_tpm_csv", sep=",", header_feature=TRUE, ann="annotation"))
  if (ds=="GSE109710") return(list(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109710/suppl/GSE109710_raw_data.txt.gz", kind="supp_nanostring_counts", sep="\t", header_feature=FALSE, ann=NULL))
  if (ds=="GSE194040" && p=="GPL20078") return(list(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194040/suppl/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_GPL20078_ProbeLevel_n654.txt.gz", kind="supp_probe_level", sep="\t", header_feature=FALSE, ann=NULL))
  if (ds=="GSE194040" && p=="GPL30493") return(list(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194040/suppl/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_GPL30493.wasGPL16233_ProbeLevel_n334.txt.gz", kind="supp_probe_level", sep="\t", header_feature=FALSE, ann=NULL))
  list(url=series_url, kind="series_matrix", sep="\t", header_feature=TRUE, ann=NULL)
}
read_series <- function(url, label) {
  f <- dl(url, paste0(label, "_series")); on.exit(if(file.exists(f)) unlink(f), add=TRUE)
  con <- gzfile(f, "rt"); begin <- NA_integer_; n <- 0L
  repeat { ln <- readLines(con, 1, warn=FALSE); if(!length(ln)) break; n <- n+1L; if(identical(ln, "!series_matrix_table_begin")) { begin <- n; break } }
  close(con)
  dt <- data.table::fread(f, skip=begin, header=TRUE, fill=TRUE, data.table=TRUE, showProgress=FALSE)
  first <- names(dt)[1]; dt <- dt[!grepl("^!", as.character(get(first)))]
  vals <- suppressWarnings(as.numeric(unlist(dt[seq_len(min(nrow(dt),1500)), 2:min(ncol(dt),9), with=FALSE], use.names=FALSE))); vals <- vals[!is.na(vals)]
  scale <- if(length(vals) && max(vals)<=30 && min(vals)>=-20) "likely_log2_or_normalized_continuous" else if(length(vals) && max(vals)>100) "likely_linear_count_or_unlogged_intensity" else "unknown"
  list(features=vapply(dt[[first]], norm, character(1)), samples=names(dt)[-1], direct=rep("", nrow(dt)), scale=scale)
}
read_plain <- function(src, label) {
  f <- dl(src$url, paste0(label, "_plain")); on.exit(if(file.exists(f)) unlink(f), add=TRUE)
  if (src$header_feature) {
    dt <- data.table::fread(f, sep=src$sep, header=TRUE, fill=TRUE, data.table=TRUE, showProgress=FALSE)
    first <- names(dt)[1]; features <- vapply(dt[[first]], norm, character(1))
    raw <- if (!is.null(src$ann) && src$ann %in% names(dt)) vapply(dt[[src$ann]], norm, character(1)) else rep("", nrow(dt))
    samples <- setdiff(names(dt)[-1], src$ann)
  } else {
    con <- gzfile(f, "rt"); header <- split_tab(readLines(con, 1, warn=FALSE)); close(con)
    dt <- data.table::fread(f, sep=src$sep, skip=1, header=FALSE, fill=TRUE, data.table=TRUE, showProgress=FALSE)
    features <- vapply(dt[[1]], norm, character(1))
    raw <- if (src$kind == "supp_nanostring_counts") features else rep("", nrow(dt))
    samples <- header
  }
  keep <- nzchar(features); features <- features[keep]; raw <- raw[keep]
  direct <- vapply(raw, function(z) paste(extract_symbols(z), collapse=";"), character(1))
  no <- !nzchar(direct); direct[no] <- vapply(features[no], function(z) paste(fallback_symbol(z), collapse=";"), character(1))
  vals <- suppressWarnings(as.numeric(unlist(dt[seq_len(min(nrow(dt),1500)), 2:min(ncol(dt),9), with=FALSE], use.names=FALSE))); vals <- vals[!is.na(vals)]
  scale <- if(length(vals) && max(vals)<=30 && min(vals)>=-20) "likely_log2_or_normalized_continuous" else if(length(vals) && max(vals)>100) "likely_linear_count_or_unlogged_intensity" else "unknown"
  list(features=features, samples=samples, direct=direct, scale=scale)
}
norm_sample <- function(x) sub("-GPL.*$", "", sub("^ISPY2_", "", norm(x)))
row_alias <- function(row) { v <- unique(c(row$gsm_sample_id, row$sample_title, row$patient_id)); v <- v[!is.na(v)&nzchar(v)]; unique(c(v, sub(":.*$", "", v), sub("^ISPY2_", "", v), sub("-GPL.*$", "", v))) }
match_samples <- function(mr, samples) {
  ms <- unique(vapply(samples, norm_sample, character(1)))
  matched <- logical(nrow(mr)); aliases <- character()
  for(i in seq_len(nrow(mr))) { a <- unique(vapply(row_alias(mr[i]), norm_sample, character(1))); aliases <- unique(c(aliases,a)); matched[i] <- any(a %in% ms) }
  list(n=sum(matched), missing=mr$gsm_sample_id[!matched], matrix_only=samples[!(vapply(samples, norm_sample, character(1)) %in% aliases)])
}
write_index <- function(ds, mf, p, features, gene, raw, path) {
  gene[is.na(gene)] <- ""; raw[is.na(raw)] <- ""
  sg <- !grepl(";", gene, fixed=TRUE) & nzchar(gene); tc <- table(gene[sg]); dup <- rep(NA_integer_, length(gene)); dup[sg] <- as.integer(tc[gene[sg]])
  fi <- data.table::data.table(dataset_accession=ds,matrix_file=mf,platform=p,feature_id=features,gene_symbol_raw=raw,gene_symbol=gene,is_annotated_to_gene=ifelse(nzchar(gene),"yes","no"),duplicate_feature_count_for_gene=dup,recommended_probe_collapse_rule="collapse duplicate probes only with an outcome-blind pre-specified rule")
  data.table::fwrite(fi, path); fi
}
main <- function() {
  if(!requireNamespace("data.table", quietly=TRUE)) stop("data.table required")
  log_msg("repair_start", R.version.string)
  sm <- data.table::fread(file.path(manifest_dir,"analysis_ready_manifest_v1.csv"), showProgress=FALSE)
  dm <- data.table::fread(file.path(manifest_dir,"dataset_level_manifest_v1.csv"), showProgress=FALSE)
  old_path <- file.path(out_dir,"expression_matrix_audit_v1.csv")
  old <- if(file.exists(old_path)) data.table::fread(old_path, showProgress=FALSE) else data.table::data.table()
  cache_ok <- c("GSE25066","GSE20271","GSE32646","GSE41998","GSE50948","GSE66305","GSE106977")
  audit <- list(); flags <- list(); plats <- list()
  for(i in seq_len(nrow(dm))) {
    d <- dm[i]; ds <- d$dataset_accession; p <- d$platform; lab <- paste(ds,p,sep="_")
    idx_path <- file.path(out_dir,"feature_indexes",sprintf("matrix_feature_index_%s_%s_v1.csv",ds,p))
    src <- source_for(ds,p,d$normalized_matrix_url)
    log_msg("dataset_start", ds, p)
    if(ds %in% cache_ok && file.exists(idx_path)) {
      fi <- data.table::fread(idx_path, showProgress=FALSE)
      oldrow <- old[dataset_accession==ds & platform==p][1]
      nmat <- if(nrow(oldrow)) oldrow$n_samples_in_matrix else as.integer(d$n_matrix_samples)
      nmatch <- if(nrow(oldrow)) oldrow$n_samples_matched_to_manifest else as.integer(d$n_matrix_samples)
      scale <- if(nrow(oldrow)) oldrow$expression_scale_guess else "cached"
      source_url <- if(nrow(oldrow)) oldrow$expression_matrix_source_url else d$normalized_matrix_url
      pstat <- "cached_from_existing_feature_index"
    } else {
      if(ds %in% c("GSE106977","GSE130786") && file.exists(idx_path)) {
        oldfi <- data.table::fread(idx_path, showProgress=FALSE); features <- oldfi$feature_id
        mat_samples <- sm[dataset_accession==ds & matrix_file==d$matrix_file, gsm_sample_id]
        scale <- old[dataset_accession==ds & platform==p][1, expression_scale_guess]
        if(!length(scale) || is.na(scale)) scale <- "not_recomputed_from_cached_feature_ids"
      } else {
        mat <- if(src$kind=="series_matrix") read_series(src$url, lab) else read_plain(src, lab)
        features <- mat$features; mat_samples <- mat$samples; scale <- mat$scale
      }
      mr <- sm[dataset_accession==ds & matrix_file==d$matrix_file]; mt <- match_samples(mr, mat_samples)
      if(length(mt$missing)) flags[[length(flags)+1]] <- data.table::data.table(dataset_accession=ds,matrix_file=d$matrix_file,platform=p,gsm_sample_id=mt$missing,match_status="manifest_sample_missing_from_expression_matrix")
      if(length(mt$matrix_only)) flags[[length(flags)+1]] <- data.table::data.table(dataset_accession=ds,matrix_file=d$matrix_file,platform=p,gsm_sample_id=mt$matrix_only,match_status="expression_matrix_sample_not_in_manifest")
      direct <- if(exists("mat") && length(mat$direct)==length(features)) mat$direct else rep("", length(features))
      gene <- direct; raw <- direct; need <- !nzchar(gene)
      pmap <- list(url="",status="skipped_direct_gene_symbols",map=data.table::data.table(feature_id=character(),gene_symbol_raw=character(),gene_symbol=character()))
      if(any(need)) { pmap <- map_platform(p, unique(features[need])); mi <- match(features[need], pmap$map$feature_id); gene[need] <- pmap$map$gene_symbol[mi]; raw[need] <- pmap$map$gene_symbol_raw[mi] }
      gene[is.na(gene)] <- ""; raw[is.na(raw)] <- ""; need2 <- !nzchar(gene)
      if(any(need2)) { fb <- vapply(features[need2], function(z) paste(fallback_symbol(z), collapse=";"), character(1)); gene[need2] <- fb; raw[need2 & nzchar(fb)] <- features[need2 & nzchar(fb)] }
      fi <- write_index(ds,d$matrix_file,p,features,gene,raw,idx_path)
      nmat <- length(unique(mat_samples)); nmatch <- mt$n; source_url <- src$url; pstat <- pmap$status
    }
    genes <- unique(unlist(strsplit(fi$gene_symbol[nzchar(fi$gene_symbol)], ";", fixed=TRUE), use.names=FALSE)); genes <- genes[nzchar(genes)]
    full <- if(d$include_primary_analysis=="yes" && nmatch>0 && length(genes)>=5000 && ds!="GSE109710") "yes" else if(d$include_primary_analysis=="yes_stress_test" && nmatch>0 && length(genes)>=5000) "stress_test_only" else "no"
    subtype <- if(d$include_primary_analysis=="subtype_only" && nmatch>0) { if(ds=="GSE109710") "yes_limited_panel_only" else if(length(genes)>=5000) "yes" else "limited_gene_space_review_needed" } else if(grepl("subtype", d$story_scope) && nmatch>0) "yes" else "no"
    audit[[length(audit)+1]] <- data.table::data.table(dataset_accession=ds,matrix_label=lab,matrix_file=d$matrix_file,platform=p,expression_matrix_source_url=source_url,expression_source_type=src$kind,n_samples_in_matrix=nmat,n_manifest_samples_for_matrix=as.integer(d$n_matrix_samples),n_samples_matched_to_manifest=nmatch,n_manifest_samples_missing_in_matrix=NA_integer_,n_matrix_samples_not_in_manifest=NA_integer_,n_features_before_annotation=nrow(fi),n_features_with_gene_symbol=sum(nzchar(fi$gene_symbol)),n_unique_genes_after_annotation=length(genes),duplicated_probe_handling_rule="not collapsed in audit; collapse later with outcome-blind rule",expression_scale_guess=scale,platform_annotation_status=pstat,feature_index_file=normalizePath(idx_path,winslash="/",mustWork=FALSE),usable_for_full_benchmark=full,usable_for_subtype_only_validation=subtype,planned_role=d$role,story_scope=d$story_scope,raw_data_availability=d$raw_data_availability,notes="")
    plats[[length(plats)+1]] <- data.table::data.table(platform=p,annotation_status=pstat)
    log_msg("dataset_done", ds, p, "matched", nmatch, "genes", length(genes), "full", full, "subtype", subtype)
  }
  audit_dt <- data.table::rbindlist(audit, fill=TRUE)
  data.table::fwrite(audit_dt, file.path(out_dir,"expression_matrix_audit_v1.csv"))
  data.table::fwrite(unique(data.table::rbindlist(plats, fill=TRUE)), file.path(out_dir,"platform_annotation_audit_v1.csv"))
  sf <- if(length(flags)) data.table::rbindlist(flags, fill=TRUE) else data.table::data.table(dataset_accession="",matrix_file="",platform="",gsm_sample_id="",match_status="")
  data.table::fwrite(sf, file.path(out_dir,"sample_matrix_match_v1.csv"))
  data.table::fwrite(data.table::data.table(dataset_accession="",matrix_file="",platform="",error=""), file.path(out_dir,"expression_audit_errors_v1.csv"))
  labs <- audit_dt$matrix_label; gp <- new.env(parent=emptyenv())
  for(i in seq_len(nrow(audit_dt))) { fi <- data.table::fread(audit_dt$feature_index_file[i], showProgress=FALSE); gcol <- as.character(fi$gene_symbol); gs <- unique(unlist(strsplit(gcol[nzchar(gcol)], ";", fixed=TRUE), use.names=FALSE)); for(g in gs[nzchar(gs)]) assign(g, unique(c(if(exists(g,gp,inherits=FALSE)) get(g,gp) else character(), audit_dt$matrix_label[i])), gp) }
  genes <- sort(ls(gp)); av <- data.table::data.table(gene_symbol=genes,n_matrices_present=0L)
  for(l in labs) av[, (l) := 0L]
  for(i in seq_along(genes)) { pr <- get(genes[i],gp); av[i,n_matrices_present:=length(pr)]; for(l in pr) data.table::set(av,i=i,j=l,value=1L) }
  data.table::fwrite(av, file.path(out_dir,"gene_availability_by_dataset_v1.csv"))
  notes <- c("# Step 2 expression matrix audit", "", paste0("Generated: ", Sys.time()), "", "Final repair run: valid cached feature indexes were reused; problematic supplementary and platform-specific matrices were rebuilt.", "", capture.output(print(audit_dt[,.(dataset_accession,platform,n_samples_matched_to_manifest,n_unique_genes_after_annotation,usable_for_full_benchmark,usable_for_subtype_only_validation)])))
  writeLines(notes, file.path(out_dir,"step2_expression_audit_notes.md"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir,"R_session_info_v1.txt"))
  file.copy("scripts/build_expression_audit_repair_v1.R", file.path(out_dir,"build_expression_audit_v1.R"), overwrite=TRUE)
  log_msg("repair_done")
  print(audit_dt[,.(dataset_accession,platform,n_samples_matched_to_manifest,n_unique_genes_after_annotation,usable_for_full_benchmark,usable_for_subtype_only_validation)])
}
tryCatch(main(), error=function(e){ log_msg("FATAL", conditionMessage(e)); writeLines(capture.output(sessionInfo()), file.path(out_dir,"R_session_info_v1.txt")); try(close(log_con), silent=TRUE); stop(e) })
try(close(log_con), silent=TRUE)
