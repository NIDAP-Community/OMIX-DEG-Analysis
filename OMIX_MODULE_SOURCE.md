# Canonical OMIX Module Source

## Canonical module

- **Module:** [OMIX DEG Analysis](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-DEG-Analysis)
- **Canonical path:** `modules/OMIX-DEG-Analysis/`
- **Released source reference:** [`2987296a45a2a11e6c4b3454c813a63ae6db6ac0`](https://github.com/NIDAP-Community/OMIX/commit/2987296a45a2a11e6c4b3454c813a63ae6db6ac0)
- **Interface schema:** [schemas/interface.yml](https://github.com/NIDAP-Community/OMIX/blob/main/modules/OMIX-DEG-Analysis/schemas/interface.yml)
- **Module contract:** [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md)

## Exported scientific files

| Canonical file | Adapter copy | Purpose |
| --- | --- | --- |
| `R/OMIX_DEG_Analysis.R` | `code/functions/OMIX_DEG_Analysis.R` | Raw-count differential-expression implementation. |

The listed export was verified byte-for-byte against the released source
reference above.

## Adapter-only historical reference

`code/functions/DEG_Paired_v53.R` is retained as a legacy validation reference.
It is not the canonical OMIX implementation and must not receive permanent
scientific changes in this deployment repository.

## Ownership and synchronization

The canonical module owns scientific functions, portable CLI behavior, schemas,
tests, and scientific documentation. This adapter owns deployment UI, input
discovery, result paths, runtime setup, and the platform entry point.

Make reusable changes in the canonical module, update its tests and interface,
record the next released source reference here, and then re-export the listed
file without unreviewed behavioral changes. Validate the adapter with
representative deployment inputs before release.
