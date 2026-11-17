#!/usr/bin/env Rscript
# GSE256068 Stage 2: reproducible human TLE-HS differential expression.
# Input is submitter-processed log-scale expression, NOT integer raw counts.
suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly=TRUE)
arg <- function(flag, default=NULL) {
  pos <- match(flag, args)
  if (is.na(pos)) return(default)
  if (pos == length(args)) stop("Missing value for ", flag)
  args[pos+1L]
}
expression_path <- arg("--expression", "stage2/data/GSE256068_processed_data.csv.gz")
series_path <- arg("--series-matrix", "stage2/data/GSE256068_series_matrix.txt.gz")
out <- arg("--out", "stage2/results")
expect_case <- as.integer(arg("--expected-case", "64"))
expect_control <- as.integer(arg("--expected-control", "13"))
skip_pathways <- "--skip-pathways" %in% args
dir.create(out, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "figures"), recursive=TRUE, showWarnings=FALSE)
stopifnot(file.exists(expression_path), file.exists(series_path))

read_metadata <- function(path) {
  lines <- readLines(gzfile(path), warn=FALSE)
  field <- function(key) {
    matches <- grep(paste0("^!", key, "\t"), lines, value=TRUE)
    if (length(matches) != 1L) stop("Expected one GEO metadata row for ", key)
    values <- strsplit(matches, "\t", fixed=TRUE)[[1]][-1L]
    sub('"$', "", sub('^"', "", values))
  }
  gsm <- field("Sample_geo_accession")
  title <- field("Sample_title")
  if (length(gsm) != length(title)) stop("GEO sample metadata lengths differ")
  meta <- data.frame(gsm=gsm, title=title, age=NA_real_,
                     sex=NA_character_, tissue=NA_character_,
                     disease_state=NA_character_, stringsAsFactors=FALSE)
  char_lines <- grep("^!Sample_characteristics_ch1\t", lines, value=TRUE)
  for (line in char_lines) {
    values <- strsplit(line, "\t", fixed=TRUE)[[1]][-1L]
    values <- sub('"$', "", sub('^"', "", values))
    if (length(values) != nrow(meta)) stop("GEO characteristics row length mismatch")
    for (i in seq_along(values)) {
      entry <- values[i]
      if (!grepl(":", entry, fixed=TRUE)) next
      key <- tolower(trimws(sub(":.*$", "", entry)))
      value <- trimws(sub("^[^:]*:", "", entry))
      if (key == "age") meta$age[i] <- suppressWarnings(as.numeric(value))
      if (key == "sex") meta$sex[i] <- value
      if (key == "tissue") meta$tissue[i] <- value
      if (key == "disease state") meta$disease_state[i] <- value
    }
  }
  meta$group <- ifelse(grepl("^TLE-HS_", meta$title), "TLE-HS",
                       ifelse(grepl("^ControlHippocampus_", meta$title),
                              "ControlHippocampus", NA_character_))
  meta <- meta[!is.na(meta$group), , drop=FALSE]
  meta$case <- as.integer(meta$group == "TLE-HS")
  meta$sex <- ifelse(tolower(meta$sex) %in% c("female", "f"), "Female",
                     ifelse(tolower(meta$sex) %in% c("male", "m"), "Male", NA_character_))
  meta$sex <- factor(meta$sex, levels=c("Female", "Male"))
  if (anyNA(meta$age) || anyNA(meta$sex)) stop("Age/sex missing: inspect GEO metadata")
  if (anyDuplicated(meta$gsm) || anyDuplicated(meta$title)) stop("Duplicate GEO sample ID/title")
  meta
}

meta <- read_metadata(series_path)
nc <- sum(meta$case == 1L); nh <- sum(meta$case == 0L)
if (nc != expect_case || nh != expect_control)
  stop(sprintf("Cohort mismatch: TLE-HS=%d, controls=%d; expected %d/%d",
               nc, nh, expect_case, expect_control))
if (length(unique(meta$sex)) < 2L) stop("Sex factor has only one observed level")

message("Reading processed expression matrix...")
tab <- read.csv(gzfile(expression_path), check.names=FALSE,
                stringsAsFactors=FALSE)
if (ncol(tab) < nrow(meta)+1L) stop("Expression matrix has too few sample columns")
genes <- sub("\\.[0-9]+$", "", as.character(tab[[1L]]))
if (mean(grepl("^ENSG[0-9]+", genes)) < 0.8)
  stop("First expression column does not appear to contain Ensembl gene IDs")
if (anyDuplicated(genes)) stop("Duplicate Ensembl IDs after removing version suffix")
expr_names <- colnames(tab)[-1L]
index <- match(meta$gsm, expr_names)
if (anyNA(index)) index <- match(meta$title, expr_names)
if (anyNA(index)) {
  missing <- paste(head(meta$title[is.na(index)], 5L), collapse=", ")
  stop("Could not align expression columns to GEO sample IDs or titles: ", missing)
}
if (anyDuplicated(index)) stop("Duplicate sample-to-expression mapping")
selected <- tab[, index+1L, drop=FALSE]
y <- as.matrix(selected)
storage.mode(y) <- "double"
rownames(y) <- genes
colnames(y) <- meta$gsm
if (any(!is.finite(y))) stop("Nonfinite processed expression values")
if (all(abs(y-round(y)) < 1e-8) && all(y >= 0))
  stop("Input appears integer counts, not expected processed log expression")
meta$expression_column <- expr_names[index]
write.csv(meta, file.path(out, "selected_77_metadata.csv"), row.names=FALSE)
write.csv(data.frame(gsm=meta$gsm, group=meta$group, age=meta$age,
                     sex=meta$sex, median_log_expression=apply(y,2,median)),
          file.path(out, "sample_qc.csv"), row.names=FALSE)

# Age and sex adjustment is essential, but does not eliminate tissue,
# postmortem-versus-resection, batch, medication or cell-composition bias.
design <- model.matrix(~ case + age + sex, data=meta)
if (qr(design)$rank != ncol(design)) stop("Primary design matrix is rank deficient")
fit <- eBayes(lmFit(y, design))
tt <- topTable(fit, coef="case", number=Inf, sort.by="P", adjust.method="BH")
tt$ensembl_gene <- rownames(tt)
tt <- tt[, c("ensembl_gene", setdiff(names(tt), "ensembl_gene"))]
tt$gene_symbol <- NA_character_
if (requireNamespace("org.Hs.eg.db", quietly=TRUE) &&
    requireNamespace("AnnotationDbi", quietly=TRUE)) {
  tt$gene_symbol <- unname(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys=tt$ensembl_gene,
    column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
}
tt$FDR_significant <- !is.na(tt$adj.P.Val) & tt$adj.P.Val < 0.05
tt$FDR_and_effect <- tt$FDR_significant & abs(tt$logFC) >= log2(1.5)
write.csv(tt, file.path(out, "DEG_limma_age_sex_adjusted.csv"), row.names=FALSE)

fit0 <- eBayes(lmFit(y, model.matrix(~ case, data=meta)))
tt0 <- topTable(fit0, coef="case", number=Inf, sort.by="P", adjust.method="BH")
tt0$ensembl_gene <- rownames(tt0)
write.csv(tt0, file.path(out, "DEG_limma_unadjusted_sensitivity.csv"), row.names=FALSE)

# Restrict to the shared observed age interval as a sensitivity analysis.
cases <- meta$age[meta$case==1L]; controls <- meta$age[meta$case==0L]
lo <- max(min(cases), min(controls))
hi <- min(max(cases), max(controls))
overlap <- meta$age >= lo & meta$age <= hi
overlap_note <- sprintf("Age overlap: %g to %g years; TLE-HS=%d, controls=%d",
                        lo, hi, sum(overlap & meta$case==1L),
                        sum(overlap & meta$case==0L))
if (sum(overlap & meta$case==1L) >= 5L &&
    sum(overlap & meta$case==0L) >= 5L) {
  m2 <- droplevels(meta[overlap, , drop=FALSE])
  formula2 <- if (length(unique(m2$sex)) == 2L) ~ case + age + sex else ~ case + age
  d2 <- model.matrix(formula2, data=m2)
  if (qr(d2)$rank == ncol(d2)) {
    f2 <- eBayes(lmFit(y[,overlap,drop=FALSE], d2))
    t2 <- topTable(f2, coef="case", number=Inf, adjust.method="BH")
    t2$ensembl_gene <- rownames(t2)
    write.csv(t2, file.path(out, "DEG_age_overlap_sensitivity.csv"), row.names=FALSE)
  } else overlap_note <- paste(overlap_note, "(design rank deficient)")
} else overlap_note <- paste(overlap_note, "(insufficient group sizes)")

# Base R figures, no extra graphics dependencies.
group_col <- ifelse(meta$case == 1L, "#D55E00", "#0072B2")
png(file.path(out, "figures", "PCA_selected_77.png"), width=1300, height=1000, res=150)
v <- apply(y, 1, var)
ix <- order(v, decreasing=TRUE)[seq_len(min(1000L, nrow(y)))]
pc <- prcomp(t(y[ix,,drop=FALSE]), center=TRUE, scale.=FALSE)
pct <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
plot(pc$x[,1], pc$x[,2], col=group_col, pch=19,
     xlab=paste0("PC1 (",pct[1],"%)"), ylab=paste0("PC2 (",pct[2],"%)"),
     main="GSE256068: TLE-HS vs hippocampus controls")
legend("topright", c("TLE-HS","Control hippocampus"),
       col=c("#D55E00","#0072B2"), pch=19, bty="n")
dev.off()
write.csv(data.frame(gsm=meta$gsm, group=meta$group,
                     PC1=pc$x[,1], PC2=pc$x[,2]),
          file.path(out,"PCA_coordinates.csv"), row.names=FALSE)

png(file.path(out, "figures", "Volcano_age_sex_adjusted.png"),
    width=1300, height=1000, res=150)
volcano_col <- ifelse(tt$FDR_and_effect, "#D55E00", "#777777")
plot(tt$logFC, -log10(pmax(tt$adj.P.Val, 1e-300)),
     pch=16, cex=0.5, col=volcano_col,
     xlab="log2 fold change (TLE-HS vs control)",
     ylab="-log10(BH adjusted p-value)",
     main="Exploratory limma differential expression")
abline(v=c(-log2(1.5),log2(1.5)), lty=2, col="gray50")
abline(h=-log10(0.05), lty=2, col="gray50")
dev.off()

top <- head(tt$ensembl_gene, min(30L,nrow(tt)))
z <- y[top,,drop=FALSE]
z <- t(scale(t(z)))
z[!is.finite(z)] <- 0
png(file.path(out, "figures", "Top30_heatmap.png"),
    width=1800, height=1100, res=150)
heatmap(z, scale="none", ColSideColors=group_col,
        col=colorRampPalette(c("#2166AC","white","#B2182B"))(101),
        margins=c(7,8), cexCol=0.3, cexRow=0.65,
        main="Top 30 genes by adjusted model p-value")
dev.off()

# Competitive gene-set analysis; optional because msigdbr is externally sourced.
pathway_notes <- character()
if (!skip_pathways && requireNamespace("msigdbr", quietly=TRUE)) {
  collections <- list(
    Hallmark=list(collection="H"),
    Reactome=list(collection="C2", subcollection="CP:REACTOME"),
    GO_BP=list(collection="C5", subcollection="GO:BP"),
    KEGG=list(collection="C2", subcollection="CP:KEGG_LEGACY")
  )
  contrast <- as.numeric(colnames(design) == "case")
  for (nm in names(collections)) {
    tryCatch({
      ms <- do.call(msigdbr::msigdbr, c(list(species="Homo sapiens"), collections[[nm]]))
      if (!all(c("gs_name","ensembl_gene") %in% colnames(ms)))
        stop("MSigDB data missing gene-set/Ensembl columns")
      gene_sets <- split(ms$ensembl_gene[!is.na(ms$ensembl_gene)],
                         ms$gs_name[!is.na(ms$ensembl_gene)])
      idx <- limma::ids2indices(lapply(gene_sets, unique), rownames(y))
      idx <- idx[vapply(idx,length,integer(1)) >= 10L &
                 vapply(idx,length,integer(1)) <= 500L]
      if (length(idx) == 0L) stop("No gene sets with 10-500 measured genes")
      enrichment <- limma::camera(y, index=idx, design=design, contrast=contrast)
      enrichment$gene_set <- rownames(enrichment)
      write.csv(enrichment, file.path(out, paste0("Pathway_camera_",nm,".csv")),
                row.names=FALSE)
      pathway_notes <- c(pathway_notes,
                         sprintf("%s: %d sets tested",nm,nrow(enrichment)))
    }, error=function(e) {
      pathway_notes <<- c(pathway_notes, paste(nm, "FAILED:", conditionMessage(e)))
    })
  }
} else pathway_notes <- "Pathway analysis skipped or msigdbr unavailable"
writeLines(pathway_notes, file.path(out,"pathway_status.txt"))

write.csv(data.frame(
  metric=c("n_TLE_HS","n_control_hippocampus","n_genes",
           "n_FDR_lt_0.05","n_FDR_and_abs_logFC_ge_log2_1.5",
           "age_median_TLE_HS","age_median_control"),
  value=c(nc,nh,nrow(y),sum(tt$FDR_significant),
          sum(tt$FDR_and_effect),median(cases),median(controls))),
  file.path(out,"analysis_summary.csv"), row.names=FALSE)
writeLines(c(
  "Dataset: GEO GSE256068, submitter-processed expression matrix.",
  "Primary comparison: TLE-HS (case=1) vs ControlHippocampus (case=0).",
  "Primary model: limma lmFit + eBayes, covariates age and sex.",
  "Sensitivity: unadjusted group-only limma and age-overlap subset.",
  "Processed data are noninteger and can be negative: NEVER run DESeq2 on this file.",
  "Original voom observation-level weights are not supplied: these are not exact original-study DE results.",
  "This is observational cross-sectional human tissue, not EV treatment or causal target validation.",
  "TLE-HS GEO source-name 'Temporal' is generic; original publication describes hippocampal resections.",
  "Age, postmortem/resection handling, medication, batch, and cellular composition remain confounders.",
  "Differential expression and gene-set analysis are distinct from original paper's differential coexpression.",
  overlap_note,
  paste("Input expression MD5:", unname(tools::md5sum(expression_path))),
  paste("Input GEO metadata MD5:", unname(tools::md5sum(series_path)))
), file.path(out, "methods_and_limitations.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
message("Stage 2 limma pipeline completed: ", normalizePath(out))
