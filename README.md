# OMIX DEG Analysis

Review-stage repository for paired and repeated-measures differential gene
expression analysis. The baseline implementation is the supplied
`DEG_Paired_v53.R` template, preserved unchanged in
`code/functions/DEG_Paired_v53.R`.

## Current status

This is **not yet a production Code Ocean capsule or an OMIX monorepo module**.
It intentionally has no App Panel, runtime entry point, or environment
definition until the review establishes the supported input contract and
release behavior.

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

## Baseline provenance

- Source: `/Users/maggiec/GitHub/Maggie/NIDAP/template-library/templates/DEG_Paired_v53.R`
- SHA-256: `417684f4330fc28171084a950bfda651b7302fd367bc94659f2609c679837537`
- Imported without scientific or formatting changes.

See [docs/review.md](docs/review.md) for the decisions required before this is
promoted to a Code Ocean capsule and an OMIX module.
