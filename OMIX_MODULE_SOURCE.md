# Canonical OMIX Module Source

## Canonical module

- **Module:** [OMIX DEG Analysis](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-DEG-Analysis)
- **Canonical path:** `modules/OMIX-DEG-Analysis/`
- **Canonical module version:** `0.4.0`
- **Canonical interface version:** `2`
- **Canonical release tag:** **Pending**; canonical `0.4.0` has not yet been tagged.
- **Canonical source reference:** [`018ef0e10bf29c69166b38ee6d81d4d71d18aeca`](https://github.com/NIDAP-Community/OMIX/commit/018ef0e10bf29c69166b38ee6d81d4d71d18aeca)
- **Interface schema:** [schemas/interface.yml](https://github.com/NIDAP-Community/OMIX/blob/main/modules/OMIX-DEG-Analysis/schemas/interface.yml)
- **Module contract:** [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md)

## Exported scientific files

| Canonical file | Adapter copy | Canonical SHA-256 | Adapter SHA-256 | Purpose |
| --- | --- | --- | --- | --- |
| `R/OMIX_DEG_Analysis.R` | `code/functions/OMIX_DEG_Analysis.R` | `e0e64e4de04b91419c77d1d44ac607f62fa172bc1b99fd0936ebb3e6ae8eb13e` | `e0e64e4de04b91419c77d1d44ac607f62fa172bc1b99fd0936ebb3e6ae8eb13e` | Raw-count edgeR/limma-voom plus Harmony-mean and SCT-mean direct-limma implementation. |

This is the complete canonical `R/` tree at the source reference above. The
listed export was verified byte-for-byte; no adapter-side scientific edits are
present under `code/functions/`.

## Intentional adapter-only translation

- The adapter CLI and App Panel use the canonical input names `counts`,
  `metadata`, `pseudobulk_manifest`, and `moo`. Blank table selectors still
  trigger unambiguous mounted-input discovery; `input_type` retains the
  canonical default `table` and choices `table` or `moo`.
- Mounted-input discovery and pseudobulk-manifest translation live in
  `code/adapter/DEG_Adapter.R`. They resolve platform paths before calling the
  canonical scientific function and are not part of the managed scientific
  export.
- Mode-aware `auto` values are passed through the canonical interface: the
  effective raw-count defaults remain filtering + TMM + diagnostics + summed
  duplicate IDs, while continuous-expression modes disable count processing
  and average duplicate IDs.
- The canonical internal `output_dir` control is translated by `code/main.R`
  to `/results` when that Code Ocean mount exists; it is not exposed as an App
  Panel scientific control.

<!-- omix-adapter-contract: {"aliases": {}, "hidden_inputs": []} -->

## Adapter release record

| Field | Recorded value |
| --- | --- |
| Adapter version | Current released version: `0.1.3`; next version **Pending** an explicit release decision and platform validation. |
| Adapter release tag | Current release: [`v0.1.3`](https://github.com/NIDAP-Community/OMIX-DEG-Analysis/tree/v0.1.3); interface-2 adapter tag **Pending**. |
| Platform release | Interface-2 Code Ocean validation **Pending**. Historical v0.1.3 validation succeeded on 2026-09-04; a published release ID was not recorded. |
| Runtime identity | Immutable runtime digest **Pending** under [OMIX issue #43](https://github.com/NIDAP-Community/OMIX/issues/43); no digest is inferred here. |
| Syncweaver mapping | `.syncweaver-lock.json` **Pending** under [OMIX issue #44](https://github.com/NIDAP-Community/OMIX/issues/44). This pull request records the interim Harbor mapping and hashes. |

See the [OMIX versioning and release policy](https://github.com/NIDAP-Community/OMIX/blob/main/docs/versioning-and-releases.md). The source commit, adapter tag, platform release, and runtime identity are separate records.

## Adapter-only historical reference

`tests/manual/legacy/DEG_Paired_v53.R` is retained as a legacy validation
reference outside the managed export. It is not the canonical OMIX
implementation and must not receive permanent scientific changes in this
deployment repository.

## Ownership and synchronization

The canonical module owns scientific functions, portable CLI behavior, schemas,
tests, and scientific documentation. This adapter owns deployment UI, input
discovery, result paths, runtime setup, and the platform entry point.

Make reusable changes in the canonical module, update its tests and interface,
record the next canonical version and immutable source reference here, and
then re-export the listed file without unreviewed behavioral changes. The
adapter-only `code/adapter/DEG_Adapter.R` discovers mounted tables and
pseudobulk bundles; `code/main.R` translates platform parameters into the
canonical interface. Validate all three analysis modes in Code Ocean before
release.
