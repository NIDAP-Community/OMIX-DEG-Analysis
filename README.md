# OMIX DEG Analysis

## Canonical OMIX module

| Item | Location |
| --- | --- |
| Canonical module | [OMIX DEG Analysis](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-DEG-Analysis) |
| Interface contract | [schemas/interface.yml](https://github.com/NIDAP-Community/OMIX/blob/main/modules/OMIX-DEG-Analysis/schemas/interface.yml) |
| Development contract | [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md) |
| Released source reference | [OMIX_MODULE_SOURCE.md](OMIX_MODULE_SOURCE.md) |

The canonical module owns the portable scientific implementation, CLI, tests,
and reusable input/output contract. This public repository is its deployment
adapter.

## What this deployment adds

- The Code Ocean entry point and App Panel.
- Attached-data and upload discovery for count tables, metadata, and MOO files.
- Capsule result layout and runtime configuration.

OMIX DEG Analysis identifies genes whose expression differs between biological
groups using one of three explicit matrix contracts:

| Analysis mode | Valid matrix | Model |
| --- | --- | --- |
| `raw_counts` | Raw or pre-filtered integer-like bulk/summed-pseudobulk counts | Design-aware filtering, edgeR normalization, and limma-voom |
| `harmony_mean_expression` | Declared log2-scale Harmony-corrected cell-type × donor means | Direct limma; no TMM, voom, or second batch correction |
| `sct_mean_expression` | Log2-scale SCT donor means written from Seurat `SCT/data` | Direct limma; no TMM or voom |

The mode is part of the scientific input contract. Do not analyze normalized
or log-transformed expression as `raw_counts`, and do not send raw counts to a
continuous-expression mode.

## Quick start

1. Keep **Input Type** at its default, `table`, for a pseudobulk bundle or a
   matrix-and-metadata pair. Select `moo` only for a compatible raw-count MOO
   object.
2. Select **Analysis Mode** to match the matrix, then attach the corresponding
   input.
3. Set **Biological Group Column(s)** to the metadata column defining the
   groups to compare, for example `Group`.
4. Set **Contrasts** as `test-reference`; for example, `B-A` compares group B
   with group A.
5. Keep processing controls on **Auto** unless the selected mode and study
   design justify an explicit override.
6. Run the analysis and review `DEG_Analysis.csv`, `Sample_Metadata.csv`, and
   `run_summary.txt`. Raw-count mode also writes normalization diagnostics by
   default.

Each group in a contrast needs at least two samples. For paired or repeated
measurements, enter the column that identifies the biological donor or patient;
do not use a technical batch as a donor identifier.

## Provide input data

Choose one of the following approaches.

### Pseudobulk bundle: automatic discovery

Attach one result or Data Asset containing:

- `Pseudobulk_Manifest.dcf`;
- `Pseudobulk_Sample_Metadata.csv`; and
- exactly one of `Pseudobulk_Counts.csv`, `Harmony_Mean_Expression.csv`, or
  `SCT_Mean_Log2_Expression.csv`.

Choose `table` as **Input Type**, select the matching **Analysis Mode**, and
leave all four file selectors blank.
The adapter finds exactly one manifest recursively, validates its `matrix_type`,
and resolves the matrix and metadata from that manifest's folder. If several
manifests are attached, select **Pseudobulk Manifest** under Advanced.

### Table mode: one combined Data Asset

Attach one Data Asset that contains exactly one mode-appropriate matrix and one
matching sample-metadata table. Keep **Input Type** at `table`; leave **Matrix
File** and **Sample Metadata File** blank. OMIX discovers the pair automatically
from mode-specific names.

The DEG-compatible MOSuite low-count-filter output uses this form:

- `Filtered_Counts.csv`
- `Sample_Metadata.csv`

### Table mode: two separately selected files

Use the **Matrix File** and **Sample Metadata File** selectors to provide
both files explicitly. Explicit selections take precedence over automatic
discovery and are useful when attached data contains more than one possible
table. Two explicit tables are authoritative; attach the manifest explicitly
when its matrix-type validation is also required.

### MOO mode

Provide one compatible `.rds` multiOmicDataSet object. OMIX extracts the
raw-count layer and sample metadata automatically. Select `moo` as **Input
Type** and `raw_counts` as **Analysis Mode**. This requires the lightweight
supported object bridge in the runtime; ordinary table input does not.

### Single-cell pseudobulk input

For single-cell data, first aggregate separately for each selected cell type,
biological donor, and condition with
[OMIX Seurat Pseudobulk](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-Seurat-Pseudobulk).
Its complete output bundle can be attached directly. Do not use individual
cells as independent replicates.

## Required table structure

Matrix and metadata tables may be CSV, TSV, TXT, or tabular RDS files.

### Matrix table

- The first column contains feature identifiers; the default name is
  `GeneName`.
- Each remaining column represents one donor/sample. Values must be
  integer-like counts for `raw_counts`, declared log2 Harmony-corrected means
  for `harmony_mean_expression`, or log2 SCT donor means for
  `sct_mean_expression`.
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
- **Technical Batch Column(s)**: for raw counts, `Auto` models `Batch` only
  when that metadata column exists; otherwise it uses no batch term. An
  explicitly named missing batch column remains an error. Harmony input
  forbids a second batch correction. For SCT means, put an unhandled technical
  batch in **Additional Fixed Covariates** instead.
- **Donor / Patient Column**: biological subject identifier for paired or
  repeated-measures samples. This is not a substitute for a technical batch.

## Mode-aware processing defaults

The App Panel defaults processing controls to `Auto`. This preserves the
historical raw-count defaults while preventing count-only operations from being
applied to continuous expression.

| Control | `raw_counts` | Harmony/SCT donor means |
| --- | --- | --- |
| Low-expression filtering | Enabled | Disabled |
| Normalization | TMM | None |
| Normalization diagnostics | Written | Not applicable |
| Duplicate feature IDs | Sum | Mean |
| Adjusted expression block | Batch/donor-adjusted when requested | Supplied expression returned unchanged |

### Raw-count normalization

The mode-aware default method, **TMM**, is suitable for most bulk RNA-seq studies. It
adjusts for library-size and compositional differences before limma-voom
modelling.

| Situation | Recommended method | Notes |
| --- | --- | --- |
| Typical bulk RNA-seq experiment | **TMM** | Default; preserves biological differences while accounting for composition bias. |
| High technical noise or a severe, documented distribution mismatch | **Quantile** | Second-choice profile. It force-aligns log-CPM distributions and can suppress severe technical variation, but may attenuate subtle or global biology. |
| Hybrid sensitivity analysis after Quantile | **TMM + Quantile** | Third-choice profile. Retains TMM library-size normalization before voom-scale quantile normalization. |
| Other documented special cases | Remaining available methods | `TMM + Scale`, `TMM + Cyclic Loess`, `TMMwsp`, `RLE`, and `Upper Quartile` should be compared with diagnostics and interpreted in the experimental-design context. |

The diagnostics are part of the decision: inspect the normalization plots and
the voom mean-variance trend rather than choosing an alternative method solely
because it is available.

Every raw-count voom pass uses the fixed smoothing policy `span = 0.5` and
`adaptive.span = FALSE`. These are reproducibility constants, not App Panel
parameters; they preserve the published-runtime behavior across supported
limma releases.

## Results

The analysis writes the following primary results.

- `DEG_Analysis.csv`: differential-expression statistics for every requested
  contrast, with sample-level expression values used for downstream plotting.
- `Sample_Metadata.csv`: metadata for exactly the samples whose expression is
  included in `DEG_Analysis.csv`, in matching order.
- `run_summary.txt`: resolved input paths, manifest use, analysis mode, matrix
  type, model, normalization, fixed raw-count voom span/adaptive-span status,
  and output provenance.
- Raw-count diagnostics: `normalization_boxplots.png`,
  `normalization_densities.png`, and `voom_mean_variance.png`. They are written
  by default and can be disabled. They are not produced for Harmony/SCT modes.

For raw counts, when batch or donor adjustment is included, the expression
columns appended to the right side of the exported table are adjusted
voom-scale log-CPM values for downstream visualization. Continuous-expression
modes append their supplied donor means unchanged. Neither expression block
replaces the statistical results from the fitted design.

`DEG_Analysis.csv` and `Sample_Metadata.csv` are one downstream data bundle.
Keep both files together when attaching the result to another analysis or
visualization step.

## Common input issues

| Message or symptom | What to check |
| --- | --- |
| Manifest `matrix_type` does not match Analysis Mode | Choose the mode declared by the bundle; do not override the manifest semantics. |
| Multiple manifests found | Attach one pseudobulk result or select **Pseudobulk Manifest** explicitly. |
| Multiple matrix or metadata candidates found | Attach one combined Data Asset with exactly one pair, remove unrelated tables, or select both files explicitly. |
| `GeneName` or `Sample` is missing | Set the corresponding parameter to the exact column name in your files. |
| An explicitly selected batch column is missing | Add the named metadata column or clear the explicit value. In raw-count mode, Auto uses `Batch` only when that column exists and otherwise runs without a batch term. |
| A count-only option is rejected | Restore that option to Auto; Harmony/SCT modes do not use TMM, voom, edgeR filtering, or downstream batch removal. |
| Samples do not match | Count-table column names and metadata sample IDs must match exactly. |
| A group has fewer than two samples | Add biological replicates or revise the requested contrast. |

## Methods

In `raw_counts`, OMIX uses design-aware filtering and edgeR's TMM normalization
before limma-voom modelling with empirical-Bayes moderation. Voom trend
smoothing is fixed at span `0.5` with adaptive span disabled. In the two
continuous-expression modes it fits limma directly to declared donor means and
does not repeat count normalization. Repeated measures use the donor-aware
`duplicateCorrelation` path.

- Robinson MD, Oshlack A. A scaling normalization method for differential
  expression analysis of RNA-seq data. *Genome Biology* 2010;11:R25.
  [doi:10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25)
- Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model
  analysis tools for RNA-seq read counts. *Genome Biology* 2014;15:R29.
  [doi:10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29)

## For developers

Read [AGENTS.md](AGENTS.md) and [OMIX_MODULE_SOURCE.md](OMIX_MODULE_SOURCE.md)
before editing. Deployment fixes that affect scientific behavior, defaults, or
the reusable input/output interface must be backported to the canonical OMIX
module before release.

This repository is the deployment adapter for the portable OMIX module. The
canonical scientific implementation is maintained in the
[OMIX monorepo](https://github.com/NIDAP-Community/OMIX). See
[OMIX_MODULE_SOURCE.md](OMIX_MODULE_SOURCE.md) for the source map and
[docs/review.md](docs/review.md) for review and implementation history.
