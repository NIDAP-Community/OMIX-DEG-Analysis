# Changelog

## Unreleased

- Synchronize the managed scientific export with canonical OMIX DEG Analysis
  `0.4.0`, interface version `2`, at commit
  `018ef0e10bf29c69166b38ee6d81d4d71d18aeca`.
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
