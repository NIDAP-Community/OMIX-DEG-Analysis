# Canonical OMIX Module Source

## Canonical module

- **Module:** [OMIX DEG Analysis](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-DEG-Analysis)
- **Canonical path:** `modules/OMIX-DEG-Analysis/`
- **Canonical module version:** `0.1.3`
- **Canonical interface version:** `1`
- **Canonical release tag:** [`module/omix-deg-analysis/v0.1.3`](https://github.com/NIDAP-Community/OMIX/tree/module/omix-deg-analysis/v0.1.3)
- **Canonical source reference:** [`2987296a45a2a11e6c4b3454c813a63ae6db6ac0`](https://github.com/NIDAP-Community/OMIX/commit/2987296a45a2a11e6c4b3454c813a63ae6db6ac0)
- **Interface schema:** [schemas/interface.yml](https://github.com/NIDAP-Community/OMIX/blob/main/modules/OMIX-DEG-Analysis/schemas/interface.yml)
- **Module contract:** [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md)

## Exported scientific files

| Canonical file | Adapter copy | Purpose |
| --- | --- | --- |
| `R/OMIX_DEG_Analysis.R` | `code/functions/OMIX_DEG_Analysis.R` | Raw-count differential-expression implementation. |

The listed export was verified byte-for-byte against the canonical source
reference above.

## Adapter release record

| Field | Recorded value |
| --- | --- |
| Adapter version | `0.1.3` |
| Adapter release tag | [`v0.1.3`](https://github.com/NIDAP-Community/OMIX-DEG-Analysis/tree/v0.1.3) |
| Platform release | Successful Code Ocean validation run on 2026-09-04; a published release ID was not recorded. |
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
record the next canonical version and immutable source reference here, and then re-export the listed
file without unreviewed behavioral changes. Validate the adapter with
representative deployment inputs before release.
