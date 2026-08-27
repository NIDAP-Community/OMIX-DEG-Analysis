# OMIX DEG Analysis

## Canonical module

The reusable scientific implementation is maintained in the OMIX monorepo at
[modules/OMIX-DEG-Analysis](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-DEG-Analysis).
This repository is the deployment adapter: it owns the App Panel, input
discovery, and runtime configuration, while released scientific changes are
promoted from the canonical module.

## Current review implementation

`code/functions/OMIX_DEG_Analysis.R` is the new raw-count implementation under
review. It uses the selectable linear / repeated-measures model structure of
the Limma template, but it front-loads raw-count validation, design-aware
filtering, edgeR TMM normalisation, and limma-voom preparation.

The returned DEG table places statistical results first and one downstream
expression column per sample at the far right. Those expression values are
batch-adjusted voom values when a valid adjustment is requested; DE statistics
are always fit using the full voom model, not the adjusted matrix.

## Summary of normalization best practices

| Dataset type | Recommended workflow | Reason |
| --- | --- | --- |
| **Typical dataset** (standard knockouts, treatments, diverse human tissues) | **Normalization Method: `TMM`** (default) [1, 2] | Preserves the observed biological range while correcting RNA-composition bias. |
| **High technical noise** (for example, varying platforms or batch-heavy historical data) | **Normalization Method: `TMM + Quantile`** | Forces log-CPM distributions to align and can reduce severe, non-linear technical variation, at the risk of attenuating subtle or global biology. [2] |

`TMM + Scale` equalizes sample medians, while `TMM + Cyclic Loess` is a
slower, pairwise alternative that can be more robust when differential
expression is unbalanced. `Quantile` omits TMM and applies voom-scale quantile
normalization alone; use it only when that deliberate choice fits the study.

### Choosing a profile

- Start with **`TMM`** for a new bulk RNA-seq analysis.
- Try **`TMM + Scale`** when QC shows a modest, sample-wide median shift after
  TMM and there is good reason to regard it as technical.
- Use **`TMM + Quantile`** only when QC shows a severe technical mismatch of
  whole distributions. It is intentionally stronger and can obscure genuine
  global shifts.
- Use **`TMM + Cyclic Loess`** as a sensitivity analysis for nonlinear,
  sample-specific distribution differences, particularly with unbalanced
  differential expression.
- Use **`Quantile`**, **`TMMwsp`**, **`RLE`**, or **`Upper Quartile`** for a
  deliberate legacy/reproducibility comparison or a documented special case;
  they are not routine first choices.

Normalization does not replace modelling known technical batch variables.

### Normalization diagnostics

The runtime writes `normalization_boxplots.png`,
`normalization_densities.png`, and `voom_mean_variance.png` by default. The
first two compare filtered log-CPM values before normalization with the final
voom values used for modelling. The mean-variance plot shows the voom trend
used to estimate precision weights and is useful for checking low-count
filtering and model preparation. Disable these outputs with
`--write_normalization_diagnostics false`. When calling `omix_deg_analysis()`
directly, set `normalization_diagnostics = TRUE` and supply
`diagnostics_output_dir` to write the same files.

### References

1. Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*. 2010;11:R25. [doi:10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25)
2. [limma documentation: `voom` and `normalizeBetweenArrays`](https://bioconductor.org/packages/devel/bioc/manuals/limma/man/limma.pdf)

The older paired template remains unchanged in this repository as a validation
reference. The unpaired and paired templates in `template-library/templates/`
will be used for regression comparisons before this becomes a capsule or an
OMIX monorepo module.

Review-stage repository for paired and repeated-measures differential gene
expression analysis. The baseline implementation is the supplied
`DEG_Paired_v53.R` template, preserved unchanged in
`code/functions/DEG_Paired_v53.R`.

## Current status

This is **not yet a production Code Ocean capsule or an OMIX monorepo module**.
It has a review-stage App Panel and runtime entry point, but the runtime image,
MOO dependency overlay, and release behavior remain to be finalized.

## Input contract

The standard path is one raw-count table plus one sample-metadata table. The
adapter accepts CSV, TSV, TXT, or tabular RDS files and requires the selected
sample-ID column to match the count-table sample columns exactly.

The optional MOSuite path accepts one `.rds` containing a
`MOSuite::multiOmicDataSet` (MOO). In this mode `code/main.R` extracts the
object's **raw** count layer and its embedded sample metadata through the
optional `OmixMOSuite` bridge. It then passes the same portable counts-plus-
metadata contract to the scientific DEG function. No MOO filename convention
is required: automatic discovery accepts exactly one `.rds` and validates its
class; otherwise the App Panel's MOO file selector must be used.

MOO mode requires MOSuite, OMIX Core, and the `OmixMOSuite` bridge in the
runtime. Ordinary table input does not install or load any of those packages.
The `filt` MOO layer is useful for testing conversion but must not be used for
the raw-count edgeR/voom DEG model.

## Optional Seurat pseudobulk source

The platform-neutral
[OmixSeurat bridge](https://github.com/NIDAP-Community/OMIX/tree/main/bridges/seurat)
can convert a Seurat object into donor-by-condition pseudobulk count and
metadata tables. That conversion occurs before this adapter: select **Table**
input and supply the generated files. Direct Seurat RDS input is not yet an
adapter option, so the capsule runtime does not install Seurat or SeuratObject
for ordinary DEG use.

## Scientific scope

The template performs repeated-measures differential-expression analysis with
`edgeR`, `limma::voom`, and `limma::duplicateCorrelation`. It models the donor
or subject identifier as a blocking factor, supports one or two contrast
variables, accommodates fixed-effect covariates, and can return expression
values corrected for donor and optional technical-batch effects.

The output is a wide differential-expression table with comparison-specific
`_FC`, `_logFC`, `_SE`, `_tstat`, `_pval`, and `_adjpval` columns. This shape is
compatible with the downstream OMIX GSEA and L2P workflows after its interface
is finalized.

Each comparison group requires at least two biological samples. For single-cell
data, aggregate raw counts by donor and condition before modelling; individual
cells are not independent replicates.

## Baseline provenance

- Source: `/Users/maggiec/GitHub/Maggie/NIDAP/template-library/templates/DEG_Paired_v53.R`
- SHA-256: `417684f4330fc28171084a950bfda651b7302fd367bc94659f2609c679837537`
- Imported without scientific or formatting changes.

See [docs/review.md](docs/review.md) for the decisions required before this is
promoted to a Code Ocean capsule and an OMIX module.
