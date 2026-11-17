# Stage 2 validation checklist

- [ ] Confirm the GEO publication, data-processing steps, and whether the processed values are log2-normalized.
- [ ] Inspect first 10 rows and headers of both input files.
- [ ] Confirm GEO series-matrix metadata keys and all 77 group assignments.
- [ ] Confirm exact mapping of the 77 GSM IDs to expression columns.
- [ ] Confirm the age and sex covariates and report their distribution by group.
- [ ] Review PCA for outliers and major non-disease axes.
- [ ] Check design matrix rank and age-overlap support.
- [ ] Compare adjusted and unadjusted DEG direction and significance.
- [ ] Audit Ensembl annotation version and mapping losses.
- [ ] Record R/package/MSigDB versions and input checksums.
- [ ] Distinguish exploratory differential expression from the source paper's differential coexpression.
- [ ] Interpret gene-set results and nominate targets only after actual outputs are available.
- [ ] Report no preliminary docking or miRNA target binding results until independently computed.

## Recommended subsequent stages

Stage 2B: transparent candidate-gene prioritization using DEG effect size, sensitivity robustness, pathway consistency, and published TLE evidence.

Stage 2C: map candidate genes to experimentally supported miRNA–mRNA interactions (miRTarBase/TarBase) and predicted sites separately. Do not dock miRNA against a PCR primer sequence or infer direct binding from expression correlation.

Stage 3: experimentally test cargo loading, EV uptake, and target modulation.
