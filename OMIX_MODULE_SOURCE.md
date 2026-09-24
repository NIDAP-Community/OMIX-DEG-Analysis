# Canonical OMIX Module Source

## Canonical module

- **Module:** [OMIX DEG Analysis](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-DEG-Analysis)
- **Canonical path:** `modules/OMIX-DEG-Analysis/`
- **Canonical module version:** `0.4.0`
- **Canonical interface version:** `2`
- **Canonical release tag:** **Pending**; canonical `0.4.0` has not yet been tagged.
- **Canonical source reference:** [`163e37499c3daaa2a3e924072e7b6b86288ebef5`](https://github.com/NIDAP-Community/OMIX/commit/163e37499c3daaa2a3e924072e7b6b86288ebef5)
- **Interface schema:** [schemas/interface.yml](https://github.com/NIDAP-Community/OMIX/blob/main/modules/OMIX-DEG-Analysis/schemas/interface.yml)
- **Module contract:** [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md)

## Exported scientific files

| Canonical file | Adapter copy | Purpose |
| --- | --- | --- |
| `R/OMIX_DEG_Analysis.R` | `code/functions/OMIX_DEG_Analysis.R` | Raw-count edgeR/limma-voom plus Harmony-mean and SCT-mean direct-limma implementation. |

The listed export was verified byte-for-byte against the canonical source
reference above (SHA-256 `e0e64e4de04b91419c77d1d44ac607f62fa172bc1b99fd0936ebb3e6ae8eb13e`).

## Intentional adapter-only translation

- The canonical CLI defaults `input_type` to `table`; this adapter retains its
  historical `auto` choice so one unambiguous table bundle is preferred and a
  lone `.rds` can fall back to MOO validation.
- Mounted-input discovery and the optional explicit `pseudobulk_manifest`
  selector live in `code/functions/DEG_Adapter.R`. They resolve platform paths
  before calling the canonical scientific function.
- Mode-aware `auto` values are passed through the canonical interface: the
  effective raw-count defaults remain filtering + TMM + diagnostics + summed
  duplicate IDs, while continuous-expression modes disable count processing
  and average duplicate IDs.

## Adapter release record

| Field | Recorded value |
| --- | --- |
| Adapter version | Current released version: `0.1.3`; next version **Pending** an explicit release decision and platform validation. |
| Adapter release tag | Current release: [`v0.1.3`](https://github.com/NIDAP-Community/OMIX-DEG-Analysis/tree/v0.1.3); interface-2 adapter tag **Pending**. |
| Platform release | Interface-2 Code Ocean validation **Pending**. Historical v0.1.3 validation succeeded on 2026-09-04; a published release ID was not recorded. |
| Runtime identity | Not yet recorded as an immutable image digest or lockfile reference. |

See the [OMIX versioning and release policy](https://github.com/NIDAP-Community/OMIX/blob/main/docs/versioning-and-releases.md). The source commit, adapter tag, platform release, and runtime identity are separate records.

## Adapter-only historical reference

`code/functions/DEG_Paired_v53.R` is retained as a legacy validation reference.
It is not the canonical OMIX implementation and must not receive permanent
scientific changes in this deployment repository.

## Ownership and synchronization

The canonical module owns scientific functions, portable CLI behavior, schemas,
tests, and scientific documentation. This adapter owns deployment UI, input
discovery, result paths, runtime setup, and the platform entry point.

Make reusable changes in the canonical module, update its tests and interface,
record the next canonical version and immutable source reference here, and
then re-export the listed file without unreviewed behavioral changes. The
adapter-only `code/functions/DEG_Adapter.R` discovers mounted tables and
pseudobulk bundles; `code/main.R` translates platform parameters into the
canonical interface. Validate all three analysis modes in Code Ocean before
release.
