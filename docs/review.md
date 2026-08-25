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

## Future input-adapter plan

The DEG scientific function remains platform- and object-model-neutral: it
accepts validated raw counts plus sample metadata. OMIX Core will define that
standard input contract and shared validation helpers.

After the count-table capsule is released and validated, add optional bridges:

1. **MOSuite bridge:** use the MOSuite-supported count extractor to convert a
   selected `multiOmicDataSet` count type and its sample metadata into the
   OMIX Core input. MOSuite is required only for this input mode.
2. **Seurat bridge:** use SeuratObject to retrieve a selected assay's raw
   count layer and cell metadata, then aggregate cells into sample-level
   pseudobulk counts (usually within a selected cell type). The DEG function
   must receive pseudobulk samples, not individual cells as replicates.

Neither object-specific dependency belongs in the core DEG implementation or
the default portable runtime.

## App Panel next milestone

The first Code Ocean App Panel will expose the standard count-table path only:

- input counts and sample metadata;
- feature-ID and sample-ID columns;
- one or two biological comparison columns and requested contrasts;
- optional fixed covariates, technical batch column(s), and donor column;
- design-aware filtering and the downstream-expression output policy.

MOO and Seurat choices will be added only with their corresponding optional
bridges, rather than mixing object-specific controls into the initial panel.

The initial adapter is now defined by `.codeocean/app-panel.json` and
`code/main.R`. It supports direct counts/metadata uploads without an
environment-specific data-asset ID. The Code Ocean environment must provide
R with `edgeR`, `limma`, and `optparse` before this adapter can be run there.

## Paired-analysis test data

The CCBR bulk training dataset cannot use `Batch` as a donor/patient proxy:
for the A-versus-B subset, each batch contains a duplicate condition and
therefore represents technical replication, not a paired design. This is an
expected rejection case for the new technical-replicate guard.

Its `Replicate` column has one A and one B sample for each of three values. It
may be used only as a **mock paired computational smoke test** of the
`duplicateCorrelation` code path, clearly labeled as synthetic pairing. A
real paired validation requires de-identified data with genuine repeated
samples from at least three donors.
