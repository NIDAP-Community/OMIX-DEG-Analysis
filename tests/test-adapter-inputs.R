source(file.path("code", "adapter", "DEG_Adapter.R"))

root <- tempfile("omix-deg-adapter-inputs-")
dir.create(root, recursive = TRUE)
old_roots <- getOption("omix.deg.data_roots")
options(omix.deg.data_roots = root)

write_bundle <- function(name, analysis_mode) {
  bundle <- file.path(root, name)
  dir.create(bundle)
  matrix_type <- omix_deg_expected_manifest_type(analysis_mode)
  matrix_filename <- omix_deg_expected_matrix_filename(analysis_mode)
  write.dcf(data.frame(matrix_type = matrix_type), file.path(bundle, "Pseudobulk_Manifest.dcf"))
  utils::write.csv(
    data.frame(GeneName = "Gene1", A1 = 1, B1 = 2, check.names = FALSE),
    file.path(bundle, matrix_filename),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(Sample = c("A1", "B1"), Group = c("A", "B")),
    file.path(bundle, "Pseudobulk_Sample_Metadata.csv"),
    row.names = FALSE
  )
  bundle
}

modes <- c("raw_counts", "harmony_mean_expression", "sct_mean_expression")
for (mode in modes) {
  unlink(list.files(root, full.names = TRUE), recursive = TRUE)
  bundle <- write_bundle(mode, mode)
  resolved <- omix_deg_resolve_table_inputs("", "", "", mode, file.path(root, "code"))
  stopifnot(
    identical(basename(resolved$matrix), omix_deg_expected_matrix_filename(mode)),
    identical(basename(resolved$metadata), "Pseudobulk_Sample_Metadata.csv"),
    identical(basename(resolved$manifest$path), "Pseudobulk_Manifest.dcf"),
    identical(dirname(resolved$matrix), normalizePath(bundle))
  )
}

unlink(list.files(root, full.names = TRUE), recursive = TRUE)
raw_bundle <- write_bundle("raw", "raw_counts")
mismatch <- tryCatch(
  omix_deg_resolve_table_inputs(
    "",
    "",
    file.path(raw_bundle, "Pseudobulk_Manifest.dcf"),
    "sct_mean_expression",
    file.path(root, "code")
  ),
  error = conditionMessage
)
stopifnot(grepl("matrix_type", mismatch, fixed = TRUE))

harmony_bundle <- write_bundle("harmony", "harmony_mean_expression")
ambiguous <- tryCatch(
  omix_deg_resolve_table_inputs("", "", "", "raw_counts", file.path(root, "code")),
  error = conditionMessage
)
stopifnot(grepl("Multiple pseudobulk manifest candidates", ambiguous, fixed = TRUE))

explicit <- omix_deg_resolve_table_inputs(
  file.path(raw_bundle, "Pseudobulk_Counts.csv"),
  file.path(raw_bundle, "Pseudobulk_Sample_Metadata.csv"),
  "",
  "raw_counts",
  file.path(root, "code")
)
stopifnot(is.null(explicit$manifest))

override <- omix_deg_resolve_table_inputs(
  "",
  "",
  file.path(harmony_bundle, "Pseudobulk_Manifest.dcf"),
  "harmony_mean_expression",
  file.path(root, "code")
)
stopifnot(identical(dirname(override$matrix), normalizePath(harmony_bundle)))

options(omix.deg.data_roots = old_roots)
unlink(root, recursive = TRUE)
message("OMIX DEG adapter discovery and manifest tests passed")
