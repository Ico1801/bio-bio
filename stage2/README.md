# Stage 2 — GSE256068 human TLE-HS transcriptomics

This is an **exploratory, reproducible** Stage 2 analysis for the NSTC CM03 goat-milk/SCAP extracellular-vesicle epilepsy proposal. It does **not** establish EV efficacy, target causality, or successful docking.

## Source and cohort

GEO accession: [GSE256068](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256068). Stage 1 selected 64 TLE-HS cases and 13 hippocampal controls (77 total) from 162 human samples. All source files must be independently verified and kept out of GitHub.

## Input

Put the **original** downloaded files into `stage2/data/` (excluded from version control):

- `GSE256068_processed_data.csv.gz`
- `GSE256068_series_matrix.txt.gz`

Do not rename or manually edit headers. The script requires the expression columns to match either GEO GSM accessions or GEO sample titles and stops if this is not true. Do not substitute the 'raw_data.csv.gz' file without independently validating its scale and meaning.

## Run

Requires R >= 4.3, packages `limma` and optional `org.Hs.eg.db`, `AnnotationDbi`, `msigdbr`. Install Bioconductor packages via BiocManager.

```r
install.packages("BiocManager")
BiocManager::install(c("limma", "org.Hs.eg.db", "AnnotationDbi"))
install.packages("msigdbr")
```

```bash
Rscript stage2/scripts/01_differential_expression.R \
  --expression stage2/data/GSE256068_processed_data.csv.gz \
  --series-matrix stage2/data/GSE256068_series_matrix.txt.gz \
  --out stage2/results
```

Optional: add `--skip-pathways` to run DEG/QC without online MSigDB resources.

The workflow intentionally **does not** automatically run on GitHub Actions, because the GEO matrix is large, remote downloads can change, and the repository is public. Run locally or in a configured R environment. Results should only be interpreted after checking `methods_and_limitations.txt`, the PCA, and cohort metadata.

## Planned outputs

- `selected_77_metadata.csv`, `sample_qc.csv`
- `DEG_limma_age_sex_adjusted.csv` (primary)
- `DEG_limma_unadjusted_sensitivity.csv`
- `DEG_age_overlap_sensitivity.csv` (when feasible)
- `PCA_coordinates.csv`, `analysis_summary.csv`
- `figures/PCA_selected_77.png`, `Volcano_age_sex_adjusted.png`, `Top30_heatmap.png`
- `Pathway_camera_{Hallmark,Reactome,GO_BP,KEGG}.csv` (where available)
- `pathway_status.txt`, `methods_and_limitations.txt`, `R_sessionInfo.txt`

## Analysis principles

1. Primary design: `~ case + age + sex`; contrast TLE-HS versus hippocampal control. No DESeq2 on transformed/noninteger values.
2. Report BH-FDR < 0.05; `|log2FC| >= log2(1.5)` is an additional descriptive effect-size filter.
3. Use unadjusted and age-overlap analyses to assess robustness; age overlap may be too limited for meaningful inference.
4. Pathway tests use limma `camera` and measured genes as background; gene sets are externally curated and version-dependent.
5. The disease comparison is observational. Differences in age, sample processing, cellular composition, medication, and batch may dominate observed signatures.
6. **Do not** interpret DEGs as direct miRNA targets, proven drug targets, or molecular docking results. Subsequent miRNA–mRNA interactions require validated sequence evidence and orthogonal experimental validation.
7. AI-assisted prioritization is a separate later step, with a prespecified evidence rubric and documented external sources.

## Status

**Code prepared; analysis not yet executed or validated against the downloaded source files.** No numerical DEG, pathway, or target discovery result is claimed in this repository.
