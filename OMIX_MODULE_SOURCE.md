# Canonical OMIX Module Source

## Canonical module

- **Module:** [OMIX DEG Analysis](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-DEG-Analysis)
- **Canonical path:** `modules/OMIX-DEG-Analysis/`
- **Canonical module version:** `0.4.2`
- **Canonical interface version:** `2`
- **Canonical release tag:** **Pending**; canonical `0.4.2` has not yet been tagged.
- **Canonical source reference:** [`11bda7fb2cbfd94c0641f864d1c3a49a678f5af1`](https://github.com/NIDAP-Community/OMIX/commit/11bda7fb2cbfd94c0641f864d1c3a49a678f5af1)
- **Interface schema:** [schemas/interface.yml](https://github.com/NIDAP-Community/OMIX/blob/11bda7fb2cbfd94c0641f864d1c3a49a678f5af1/modules/OMIX-DEG-Analysis/schemas/interface.yml)
- **Module contract:** [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md)

## Exported scientific files

| Canonical file | Adapter copy | Canonical SHA-256 | Adapter SHA-256 | Purpose |
| --- | --- | --- | --- | --- |
| `R/OMIX_DEG_Analysis.R` | `code/functions/OMIX_DEG_Analysis.R` | `44cfe51ea649ed72aa5c7c55ddfa6aed907517cad408c52a34a96ca2c692a7f2` | `44cfe51ea649ed72aa5c7c55ddfa6aed907517cad408c52a34a96ca2c692a7f2` | Raw-count edgeR/limma-voom with fixed trend smoothing plus Harmony-mean and SCT-mean direct-limma implementation. |

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
- Mode-aware `auto` values are resolved with the canonical implementation
  after inputs are loaded. In raw-count mode, automatic batch selection uses
  `Batch` only when that metadata column exists; other raw-count defaults
  remain filtering + TMM + diagnostics + summed duplicate IDs. Continuous-
  expression modes disable count processing and average duplicate IDs.
- The canonical internal `output_dir` control is translated by `code/main.R`
  to `/results` when that Code Ocean mount exists; it is not exposed as an App
  Panel scientific control.
- The adapter writes the canonical raw-count voom trend provenance to
  `run_summary.txt`: `span = 0.5` and `adaptive.span = FALSE`. Both fields are
  reported as not applicable for continuous-expression modes.

<!-- omix-adapter-contract: {"aliases": {}, "hidden_inputs": []} -->

<!-- omix-runtime-provenance: {"runtime_profile":"r-statistics","published_image":"ghcr.io/nidap-community/omix-r-statistics","published_tag":"r4.4.3-bioconductor3.20-v1","published_digest":"sha256:1325722877fec5167d171aa766ddf7bbfd056e4999bf40fd8c1eabee495da667","published_immutable_reference":"ghcr.io/nidap-community/omix-r-statistics@sha256:1325722877fec5167d171aa766ddf7bbfd056e4999bf40fd8c1eabee495da667","adapter_base_image":"codeocean/omix-r-statistics:r4.4.3-bioconductor3.20-v1","adapter_base_image_digest":null,"adapter_digest_status":"pending_code_ocean_resolution","runtime_profile_alignment_status":"aligned_by_profile_and_readable_tag_only","code_ocean_environment_id":null,"capsule_run":null} -->

## Adapter release record

| Field | Recorded value |
| --- | --- |
| Adapter version | Current released version: `0.1.3`; next version **Pending** an explicit release decision and platform validation. |
| Adapter release tag | Current release: [`v0.1.3`](https://github.com/NIDAP-Community/OMIX-DEG-Analysis/tree/v0.1.3); interface-2 adapter tag **Pending**. |
| Platform release | Canonical-0.4.2/interface-2 Code Ocean validation **Pending**. Historical v0.1.3 validation succeeded on 2026-09-04; a published release ID was not recorded. |
| Canonical runtime profile | `r-statistics`, as declared by the canonical module. |
| Published OMIX runtime | Tag `ghcr.io/nidap-community/omix-r-statistics:r4.4.3-bioconductor3.20-v1`; immutable public identity `ghcr.io/nidap-community/omix-r-statistics@sha256:1325722877fec5167d171aa766ddf7bbfd056e4999bf40fd8c1eabee495da667`. |
| Adapter-selected base image | `codeocean/omix-r-statistics:r4.4.3-bioconductor3.20-v1` from `.codeocean/environment.json`; the readable profile and tag align with the published OMIX runtime, but the Code Ocean environment's own digest is **Pending Code Ocean resolution** and is not inferred from GHCR. |
| Runtime build evidence | [Starter Environments run 32928780557](https://github.com/NIDAP-Community/OMIX/actions/runs/32928780557) built, verified, inventoried, and published the public `r-statistics` image; [GHCR version 1172700304](https://github.com/orgs/NIDAP-Community/packages/container/omix-r-statistics/1172700304) independently maps the version tag to the recorded digest. |
| Code Ocean environment identity | **Pending** — no immutable environment-release identifier is recorded. |
| Capsule run | Canonical-0.4.2/interface-2 named-parameter Code Ocean run **Pending**; GitHub or local checks are not platform validation. |
| Syncweaver mapping | `.syncweaver-lock.json` **Pending** under [OMIX issue #44](https://github.com/NIDAP-Community/OMIX/issues/44). This pull request records the interim Harbor mapping and hashes. |

The public image identity above is also recorded in the canonical
[`release-manifest.json`](https://github.com/NIDAP-Community/OMIX/blob/eee4433cd74d0ce9bd1daba1571ca6312cba7ea2/starter-environments/release-manifest.json).
It establishes the intended portable OMIX runtime; it does not prove that the
Code Ocean starter-environment import has the same digest. The source commit,
adapter tag, public image tag and digest, Code Ocean environment identity,
capsule run, and platform release remain separate records.

See the [OMIX versioning and release policy](https://github.com/NIDAP-Community/OMIX/blob/main/docs/versioning-and-releases.md).

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
