#!/usr/bin/env Rscript
# Stage 2: inspect statistical outputs without calling DEGs validated targets.
out <- "stage2/results"
must <- c("analysis_summary.csv","DEG_limma_age_sex_adjusted.csv",
          "DEG_limma_unadjusted_sensitivity.csv","selected_77_metadata.csv")
missing <- must[!file.exists(file.path(out,must))]
if(length(missing)) stop("Run Stage 2 pipeline first; missing: ",paste(missing,collapse=", "))
a <- read.csv(file.path(out,"DEG_limma_age_sex_adjusted.csv"),
              stringsAsFactors=FALSE,check.names=FALSE)
u <- read.csv(file.path(out,"DEG_limma_unadjusted_sensitivity.csv"),
              stringsAsFactors=FALSE,check.names=FALSE)
meta <- read.csv(file.path(out,"selected_77_metadata.csv"),
                 stringsAsFactors=FALSE)
if (!all(c("ensembl_gene","logFC","adj.P.Val","P.Value") %in% names(a)))
  stop("Primary DEG output lacks expected columns")
if (!all(c("ensembl_gene","logFC","adj.P.Val") %in% names(u)))
  stop("Unadjusted DEG output lacks expected columns")
if (anyDuplicated(a$ensembl_gene)) stop("Duplicated primary gene IDs")
if (any(!is.finite(a$logFC)) || any(!is.finite(a$P.Value))) stop("Nonfinite DEG statistics")
if (any(a$adj.P.Val<0 | a$adj.P.Val>1,na.rm=TRUE)) stop("FDR outside 0..1")
if (any(a$P.Value<0 | a$P.Value>1,na.rm=TRUE)) stop("P-value outside 0..1")
a$significant <- !is.na(a$adj.P.Val) & a$adj.P.Val < 0.05
a$significant_and_effect <- a$significant & abs(a$logFC)>=log2(1.5)
m <- merge(a[,c("ensembl_gene","logFC","adj.P.Val")],
           u[,c("ensembl_gene","logFC","adj.P.Val")],
           by="ensembl_gene",suffixes=c("_adjusted","_unadjusted"))
direction_agreement <- mean(sign(m$logFC_adjusted)==sign(m$logFC_unadjusted))
corr <- suppressWarnings(cor(m$logFC_adjusted,m$logFC_unadjusted,
                             method="spearman",use="complete.obs"))
top20 <- head(a[order(a$adj.P.Val,a$P.Value),,drop=FALSE],20L)
write.csv(top20,file.path(out,"top20_exploratory_DEGs_for_review.csv"),
          row.names=FALSE)
age_case <- meta$age[meta$group=="TLE-HS"]
age_ctrl <- meta$age[meta$group=="ControlHippocampus"]
lines <- c(
  "STAGE 2 REVIEW — EXPLORATORY DIFFERENTIAL EXPRESSION",
  paste("Genes tested:",nrow(a)),
  paste("FDR < 0.05:",sum(a$significant)),
  paste("FDR < 0.05 AND |log2FC| >= log2(1.5):",sum(a$significant_and_effect)),
  paste("Adjusted vs unadjusted logFC Spearman rho:",round(corr,3)),
  paste("Adjusted vs unadjusted direction agreement:",round(direction_agreement*100,1),"%"),
  paste("Median age TLE-HS:",median(age_case),"years; controls:",median(age_ctrl),"years"),
  paste("TLE-HS:",length(age_case),"; hippocampal controls:",length(age_ctrl)),
  paste("Sex by group:",paste(capture.output(print(table(meta$group,meta$sex))),collapse=" | ")),
  "Caution: age/sex adjustment does NOT correct postmortem-vs-resection, batch, drugs or cell composition."
)
if (file.exists(file.path(out,"DEG_age_overlap_sensitivity.csv"))) {
  x <- read.csv(file.path(out,"DEG_age_overlap_sensitivity.csv"),
                stringsAsFactors=FALSE)
  common <- merge(a[,c("ensembl_gene","logFC")],
                  x[,c("ensembl_gene","logFC")],by="ensembl_gene",
                  suffixes=c("_full","_overlap"))
  lines <- c(lines,paste("Age-overlap subset exists; logFC Spearman rho:",
                        round(cor(common$logFC_full,common$logFC_overlap,
                                  method="spearman"),3)))
} else lines <- c(lines,"Age-overlap subset unavailable or not estimable; check methods_and_limitations.txt")
paths <- list.files(out,pattern="^Pathway_camera_.*\\.csv$",full.names=TRUE)
for (p in paths) {
  z <- read.csv(p)
  if ("FDR" %in% names(z)) lines <- c(lines,
    paste(basename(p),":",sum(z$FDR<0.05,na.rm=TRUE),
          "gene sets with camera FDR < 0.05"))
}
lines <- c(lines,
  "Top 20 exploratory genes (NOT independently validated therapeutic targets):",
  paste(capture.output(print(top20[,intersect(c("ensembl_gene","gene_symbol","logFC","P.Value","adj.P.Val"),names(top20))],
                             row.names=FALSE)),collapse="\n"),
  "Do not interpret FDR or pathway enrichment as causal evidence or EV cargo efficacy."
)
writeLines(lines,file.path(out,"review_summary.txt"))
cat(paste(lines,collapse="\n"),"\n")
