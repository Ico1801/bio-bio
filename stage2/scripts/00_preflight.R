#!/usr/bin/env Rscript
# Fast preflight: no packages required and does not load the full expression matrix.
expression_path <- "stage2/data/GSE256068_processed_data.csv.gz"
series_path <- "stage2/data/GSE256068_series_matrix.txt.gz"
out <- "stage2/results"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
if (!file.exists(expression_path) || !file.exists(series_path)) {
  stop("Missing input(s) in stage2/data/. Check file names and working directory.")
}
get_lines <- function(path, n=-1L) {
  con <- gzfile(path, open="rt")
  on.exit(close(con))
  readLines(con, n=n, warn=FALSE)
}
meta_lines <- get_lines(series_path)
get_meta <- function(key) {
  x <- grep(paste0("^!", key, "\t"), meta_lines, value=TRUE)
  if (!length(x)) return(NULL)
  vals <- strsplit(x[1], "\t", fixed=TRUE)[[1]][-1]
  sub('"$', "", sub('^"', "", vals))
}
titles <- get_meta("Sample_title")
gsm <- get_meta("Sample_geo_accession")
descs <- get_meta("Sample_description")
if (is.null(titles) || is.null(gsm)) stop("Missing sample titles/accessions in GEO matrix")
if (is.null(descs) || length(descs)!=length(titles)) descs <- rep("",length(titles))
heads <- get_lines(expression_path, 2L)
expr_names <- scan(text=heads[1], what=character(), sep="", quiet=TRUE, quote='"')
# Some matrices have a gene-ID header; others have an empty leading cell.
if (length(expr_names)==length(gsm)+1L) expr_names <- expr_names[-1L]
group <- ifelse(grepl("^TLE-HS_", titles, ignore.case=TRUE), "TLE-HS",
                ifelse(grepl("^ControlHippocampus_", titles, ignore.case=TRUE),
                       "ControlHippocampus", "other"))
idx <- which(group!="other")
match_ids <- match(gsm[idx],expr_names)
need <- which(is.na(match_ids))
if(length(need)) match_ids[need] <- match(titles[idx][need],expr_names)
need <- which(is.na(match_ids))
if(length(need)) match_ids[need] <- match(descs[idx][need],expr_names)
char_rows <- grep("^!Sample_characteristics_ch1\t", meta_lines, value=TRUE)
char_keys <- vapply(char_rows, function(line) {
  x <- strsplit(line, "\t", fixed=TRUE)[[1]]
  if (length(x)<2L) return("(empty)")
  tolower(trimws(sub(":.*$", "", x[2L])))
}, character(1))
report <- c(
  "STAGE 2 PREFLIGHT REPORT — NO DIFFERENTIAL EXPRESSION COMPUTED",
  paste("Working directory:",getwd()),
  paste("Expression file size (MB):",round(file.info(expression_path)$size/1e6,2)),
  paste("Series-matrix size (MB):",round(file.info(series_path)$size/1e6,2)),
  paste("GEO samples:",length(gsm)),
  paste("Expression header sample columns:",length(expr_names)),
  paste("TLE-HS selected:",sum(group=="TLE-HS")),
  paste("Hippocampal controls selected:",sum(group=="ControlHippocampus")),
  paste("Selected expression IDs matched:",sum(!is.na(match_ids)),"/",length(idx)),
  paste("Sample characteristics keys:",paste(char_keys,collapse=", ")),
  paste("Expression header first fields:",paste(head(expr_names,5),collapse=", ")),
  paste("GEO sample title first fields:",paste(head(titles,5),collapse=", ")),
  paste("GEO sample description first fields:",paste(head(descs,5),collapse=", ")),
  paste("Expression first data row (truncated):",substr(heads[2],1,240)),
  paste("Input expression MD5:",unname(tools::md5sum(expression_path))),
  paste("Input GEO matrix MD5:",unname(tools::md5sum(series_path)))
)
if (anyNA(match_ids)) report <- c(report,
  paste("UNMATCHED sample titles:",paste(head(titles[idx][is.na(match_ids)],10),collapse=" | ")))
writeLines(report, file.path(out,"preflight_report.txt"))
cat(paste(report,collapse="\n"),"\n")
if (sum(group=="TLE-HS")!=64L || sum(group=="ControlHippocampus")!=13L)
  warning("Expected cohort sizes 64 TLE-HS / 13 hippocampal controls do not match")
if (anyNA(match_ids)) warning("Expression sample matching failed: do not run DEG script yet")
