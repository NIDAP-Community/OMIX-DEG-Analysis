# OMIX DEG Analysis

## Code Ocean deployment repository

This public repository is the **Code Ocean deployment repository** for OMIX DEG
Analysis. It contains the Code Ocean-ready entry point, app-panel definition,
runtime adapter, and tests needed to run the analysis as a capsule.

For portable local R use and the canonical scientific implementation, start
with the [OMIX monorepo](https://github.com/NIDAP-Community/OMIX). The two
repositories intentionally serve different roles: OMIX is platform-neutral;
this repository packages that analysis for Code Ocean.

OMIX DEG Analysis identifies genes whose expression differs between biological
groups from bulk RNA-seq count data. It performs low-expression filtering,
library-size normalization, limma-voom modelling, and optional adjustment for
technical batch or genuine repeated measurements.

Use raw or pre-filtered **integer-like counts**. Do not use normalized,
log-transformed, TPM, FPKM, or CPM expression values as input.

## Quick start

1. Provide a count table and its matching sample-metadata table, or one
   compatible MOSuite MOO object.
2. Set **Biological Group Column(s)** to the metadata column defining the
   groups to compare, for example `Group`.
3. Set **Contrasts** as `test-reference`; for example, `B-A` compares group B
   with group A.
4. Enter a **Technical Batch Column** only when one exists. Leave it blank
   otherwise.
5. Run the analysis and review `DEG_Analysis.csv` together with the diagnostic
   plots and run summary.

Each group in a contrast needs at least two samples. For paired or repeated
measurements, enter the column that identifies the biological donor or patient;
do not use a technical batch as a donor identifier.

## Provide input data

Choose one of the following approaches.

### Table mode: one combined Data Asset

Attach one Data Asset that contains exactly one count table and one matching
sample-metadata table. Leave **Raw Counts File** and **Sample Metadata File**
blank. OMIX discovers the pair automatically.

The DEG-compatible MOSuite low-count-filter output uses this form:

- `Filtered_Counts.csv`
- `Sample_Metadata.csv`

### Table mode: two separately selected files

Use the **Raw Counts File** and **Sample Metadata File** selectors to provide
both files explicitly. Explicit selections take precedence over automatic
discovery and are useful when the attached data contains more than one possible
table.

### MOO mode

Provide one compatible MOSuite `.rds` multiOmicDataSet object. OMIX extracts
the raw-count layer and sample metadata automatically. This option requires the
MOSuite bridge to be available in the runtime; ordinary table input does not.

### Single-cell pseudobulk input

For single-cell data, first aggregate raw counts by biological donor and
condition into a pseudobulk count table and matching metadata table. The
[OmixSeurat bridge](https://github.com/NIDAP-Community/OMIX/tree/main/bridges/seurat)
can prepare these portable tables from a Seurat object before the capsule is
run. Do not use individual cells as independent replicates in this bulk-style
model.

## Required table structure

Count and metadata tables may be CSV, TSV, TXT, or tabular RDS files.

### Count table

- The first column contains feature identifiers; the default name is
  `GeneName`.
- Each remaining column contains raw or pre-filtered counts for one sample.
- Column names must match the sample IDs in the metadata table.

### Sample metadata table

- One row represents one sample.
- The default sample-ID column is `Sample`.
- Include a column that defines the biological groups to compare, commonly
  `Group`.
- Include `Batch` only when a technical batch should be modelled.

The app parameters let you replace `GeneName`, `Sample`, `Group`, and `Batch`
with the exact column names in your data.

## Configure the analysis

### Contrasts

Write each contrast as `test-reference`:

- `B-A` tests B relative to A.
- `Treatment-Control` tests Treatment relative to Control.
- Separate multiple contrasts with commas.

The sign of `logFC` follows this order: a positive value for `B-A` means higher
expression in B than in A.

### Covariates, batch, and repeated measures

- **Additional Fixed Covariates**: metadata columns to include in the model,
  such as sex or a prespecified clinical variable.
- **Technical Batch Column(s)**: technical source of variation to adjust for,
  such as sequencing run or library-preparation batch.
- **Donor / Patient Column**: biological subject identifier for paired or
  repeated-measures samples. This is not a substitute for a technical batch.

## Choose normalization

The default method, **TMM**, is suitable for most bulk RNA-seq studies. It
adjusts for library-size and compositional differences before limma-voom
modelling.

| Situation | Recommended method | Notes |
| --- | --- | --- |
| Typical bulk RNA-seq experiment | **TMM** | Default; preserves biological differences while accounting for composition bias. |
| Severe, documented technical distribution mismatch | **TMM + quantile** | Use cautiously: stronger distribution alignment can also attenuate subtle biology. |
| Sensitivity analysis | Other available methods | Compare diagnostics and interpret differences in the context of the experimental design. |

The diagnostics are part of the decision: inspect the normalization plots and
the voom mean-variance trend rather than choosing an alternative method solely
because it is available.

## Results

The analysis writes the following primary results.

- `DEG_Analysis.csv`: differential-expression statistics for every requested
  contrast, with sample-level expression values used for downstream plotting.
- `Sample_Metadata.csv`: metadata for exactly the samples whose expression is
  included in `DEG_Analysis.csv`, in matching order.
- Run summary: model settings, input checks, and analysis provenance.
- Diagnostics: `normalization_boxplots.png`, `normalization_densities.png`,
  and `voom_mean_variance.png`. They are written by default and can be disabled
  with the normalization-diagnostics setting.

When batch or donor adjustment is included, the expression columns appended to
the right side of the exported table are adjusted voom-scale log-CPM values for
downstream visualization. They are not replacements for the statistical
results, which come from the complete fitted design.

`DEG_Analysis.csv` and `Sample_Metadata.csv` are one downstream data bundle.
Keep both files together when attaching the result to another analysis or
visualization step.

## Common input issues

| Message or symptom | What to check |
| --- | --- |
| Multiple counts or metadata candidates found | Attach one combined Data Asset with exactly one pair, remove unrelated tables, or select both files explicitly. |
| `GeneName` or `Sample` is missing | Set the corresponding parameter to the exact column name in your files. |
| `Batch` is missing | Leave the Technical Batch setting blank, or choose the correct metadata column. |
| Samples do not match | Count-table column names and metadata sample IDs must match exactly. |
| A group has fewer than two samples | Add biological replicates or revise the requested contrast. |

## Methods

OMIX uses design-aware filtering and edgeR's TMM normalization before
limma-voom modelling with empirical-Bayes moderation. Repeated measures are
handled through the donor-aware model path.

- Robinson MD, Oshlack A. A scaling normalization method for differential
  expression analysis of RNA-seq data. *Genome Biology* 2010;11:R25.
  [doi:10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25)
- Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model
  analysis tools for RNA-seq read counts. *Genome Biology* 2014;15:R29.
  [doi:10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29)

## For developers

This repository is the deployment adapter for the portable OMIX module. The
canonical scientific implementation is maintained in the
[OMIX monorepo](https://github.com/NIDAP-Community/OMIX). See
[OMIX_MODULE_SOURCE.md](OMIX_MODULE_SOURCE.md) for the source map and
[docs/review.md](docs/review.md) for review and implementation history.
