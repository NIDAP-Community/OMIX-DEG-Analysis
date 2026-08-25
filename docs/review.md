# Review: paired DEG analysis baseline

## Confirmed behavior

- Accepts a gene-by-sample count matrix and a sample-metadata table.
- Requires a donor/subject identifier and uses `limma::duplicateCorrelation`
  to account for repeated samples from the same individual.
- Supports one or two contrast variables and explicit limma contrasts.
- Uses `edgeR::DGEList`, TMM-family normalization, `limma::voom`, and limma
  empirical-Bayes moderation.
- Produces comparison-prefixed FC, logFC, SE, t-statistic, nominal p-value,
  and Benjamini-Hochberg adjusted-p-value columns.
- Can write a TSV output if `output_file_path` is supplied.

## Decisions required before implementation

1. **Input format:** define whether the capsule accepts separate count and
   metadata files, an OMIX/MOSuite object, or both through a documented
   conversion boundary.
2. **Sample selection:** decide how the App Panel exposes samples, donor,
   contrast variable(s), covariates, and optional batch column without
   hard-coded demo names.
3. **Replicates:** the template retains the first row for each
   donor-by-condition combination (`slice(1)`). Confirm whether technical
   replicates should instead be aggregated, rejected with an explicit error,
   or retained under another model.
4. **Filtering and counts:** define required count type, low-expression
   filtering policy, and whether the existing `input_in_log_counts` pathway is
   in scope for production use.
5. **Outputs:** establish stable output filenames and a schema for the wide
   DEG table, design/contrast diagnostics, and corrected-expression values.
6. **Runtime:** identify the shared image family. The baseline requires limma,
   edgeR, dplyr, tidyr, stringr, tibble, ggplot2, gridExtra, reshape2, gtable,
   plotly, and scales.
7. **Validation:** select a small paired fixture and expected results for
   direct tests before introducing a Code Ocean App Panel.

## Promotion path after review

1. Establish the portable CLI and tests in the canonical
   `OMIX/modules/OMIX-DEG-Analysis/` module.
2. Export the released scientific function to this adapter under
   `code/functions/`.
3. Add Code Ocean-specific `code/main.R`, `code/run`, App Panel metadata, and
   an environment overlay.
4. Validate in Code Ocean using a representative paired dataset before release.
# First implementation: raw counts plus Limma-style model selection

The first review implementation is `code/functions/OMIX_DEG_Analysis.R`.
It deliberately combines:

1. raw count validation and duplicate-gene summarisation;
2. `edgeR::filterByExpr()` using biological groups;
3. TMM normalisation and `limma::voom()`;
4. the Limma template's selectable linear or `duplicateCorrelation` model;
5. a single wide output table whose right-most columns are downstream
   expression values.

The statistical model includes all selected technical batch columns. The
right-most expression block removes selected technical batches while preserving
comparison groups. In a repeated-measures analysis it also removes the donor
effect by default, matching the useful downstream behaviour of the paired
template. Neither adjustment changes the fitted DE statistics.

Technical replicate policy: a donor-by-comparison-group combination with more
than one sample is rejected. A future aggregation policy must be explicit and
validated; the implementation must never silently select one replicate.

## Pending validation

- Run `OMIX_TEMPLATE_LIBRARY=../template-library Rscript
  tests/manual/validate-legacy-equivalence.R` from the repository root. The
  fixed synthetic fixture currently matches `DEG_Analysis_v193.R` (unpaired)
  and `DEG_Paired_v53.R` (paired) to less than `1e-8` for logFC and p-values
  when filtering and batch adjustment are disabled for an equivalent model.
- Add a real, de-identified raw-count validation fixture before release.
- Agree numerical tolerances for retained genes, coefficients, p-values, FDR,
  and the downstream expression block.
- Decide whether donor removal should remain the default downstream expression
  behaviour for every repeated-measures use case.
