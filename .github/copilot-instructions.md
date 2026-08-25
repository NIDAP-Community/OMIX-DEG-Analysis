# OMIX DEG Analysis Review Instructions

This repository is in review before a canonical OMIX module or Code Ocean
capsule interface is created.

- Preserve `code/functions/DEG_Paired_v53.R` as the identified template
  baseline unless a reviewed scientific change is explicitly requested.
- Do not create App Panel metadata, a runtime Dockerfile, or a `code/main.R`
  wrapper until the input, output, replicate, and runtime decisions in
  `docs/review.md` are resolved.
- Keep data, results, credentials, and package inventories out of Git.
- Once approved, create the platform-neutral OMIX module first and make this
  repository its Code Ocean deployment adapter.
