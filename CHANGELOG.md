# Changelog

## Unreleased

- Synchronize the byte-identical managed scientific export with canonical OMIX
  DEG Analysis `0.4.2`, interface version `2`, at commit
  `11bda7fb2cbfd94c0641f864d1c3a49a678f5af1`.
- Preserve the canonical fixed raw-count voom smoothing policy (`span = 0.5`,
  `adaptive.span = FALSE`) and record both settings in the adapter-generated
  run summary. No App Panel control or scientific default was added.
- Synchronize the managed scientific export with canonical OMIX DEG Analysis
  `0.4.1`, interface version `2`, at commit
  `b25430f52171c9c1e5fea1b093b70452e71f792a`.
- Resolve automatic raw-count batch selection after metadata is loaded: use
  `Batch` only when that column exists, otherwise run without a batch term;
  keep explicitly requested missing columns as errors.
- Add Code Ocean translation and discovery for raw-count, Harmony-mean, and
  SCTransform-mean pseudobulk bundles.
- Align the adapter CLI and App Panel with canonical input names, types,
  choices, order, and defaults while keeping unambiguous mounted-input
  discovery in the adapter layer.
- Isolate the byte-identical canonical R export under `code/functions/` and
  move adapter discovery and legacy validation code outside that managed tree.
- Add interface, input-discovery, three-mode execution, source-provenance, and
  legacy-equivalence tests.
- Record the canonical public `r-statistics` tag and immutable GHCR digest
  separately from the still-pending Code Ocean environment identity and run.

Code Ocean validation, the next adapter version and tag, the platform release
identifier, and the selected Code Ocean environment digest remain pending.

## 0.1.3

- Record the previously validated raw-count DEG adapter and its canonical
  OMIX `0.1.3` scientific source.
