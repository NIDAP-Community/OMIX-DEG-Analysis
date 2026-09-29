# Code Ocean input discovery and manifest translation for OMIX DEG Analysis.
# Scientific modeling remains in OMIX_DEG_Analysis.R.

omix_deg_canonical_data_path <- function(path) {
  path <- normalizePath(path, mustWork = TRUE, winslash = "/")
  sub("^/{2,}", "/", path)
}

omix_deg_data_roots <- function(code_dir) {
  configured_roots <- getOption("omix.deg.data_roots", NULL)
  if (!is.null(configured_roots)) {
    configured_roots <- as.character(configured_roots)
    return(unique(configured_roots[dir.exists(configured_roots)]))
  }
  roots <- c("/data", file.path(dirname(code_dir), "data"))
  unique(roots[dir.exists(roots)])
}

omix_deg_find_files <- function(code_dir, extensions) {
  roots <- omix_deg_data_roots(code_dir)
  if (length(roots) == 0L) return(character())
  extension_pattern <- paste0("[.](", paste(extensions, collapse = "|"), ")$")
  files <- unique(unlist(lapply(roots, function(root) {
    list.files(
      root,
      pattern = extension_pattern,
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
  }), use.names = FALSE))
  unique(vapply(files, omix_deg_canonical_data_path, character(1)))
}

omix_deg_validate_upload <- function(path, label, extensions) {
  if (!file.exists(path)) stop(label, " upload does not exist: ", path, call. = FALSE)
  extension <- tolower(tools::file_ext(path))
  if (!extension %in% extensions) {
    stop(label, " must use one of these extensions: ", paste(extensions, collapse = ", "), ".", call. = FALSE)
  }
  omix_deg_canonical_data_path(path)
}

omix_deg_unique_candidate <- function(candidates, label, guidance) {
  candidates <- sort(unique(candidates))
  if (length(candidates) == 1L) return(candidates)
  if (length(candidates) == 0L) {
    stop("No ", label, " was uploaded or uniquely discovered. ", guidance, call. = FALSE)
  }
  stop(
    "Multiple ", label, " candidates found. Select the intended file explicitly: ",
    paste(candidates, collapse = ", "),
    call. = FALSE
  )
}

omix_deg_expected_manifest_type <- function(analysis_mode) {
  switch(
    analysis_mode,
    raw_counts = "raw_integer_counts",
    harmony_mean_expression = "harmony_corrected_mean_expression",
    sct_mean_expression = "sctransform_mean_log2_expression",
    stop("Unsupported analysis mode: ", analysis_mode, call. = FALSE)
  )
}

omix_deg_expected_matrix_filename <- function(analysis_mode) {
  switch(
    analysis_mode,
    raw_counts = "Pseudobulk_Counts.csv",
    harmony_mean_expression = "Harmony_Mean_Expression.csv",
    sct_mean_expression = "SCT_Mean_Log2_Expression.csv",
    stop("Unsupported analysis mode: ", analysis_mode, call. = FALSE)
  )
}

omix_deg_read_manifest <- function(path, analysis_mode) {
  manifest <- tryCatch(
    base::read.dcf(path),
    error = function(error) stop("Could not read pseudobulk manifest: ", conditionMessage(error), call. = FALSE)
  )
  if (nrow(manifest) != 1L || !"matrix_type" %in% colnames(manifest)) {
    stop("Pseudobulk manifest must contain exactly one record with matrix_type.", call. = FALSE)
  }
  expected_type <- omix_deg_expected_manifest_type(analysis_mode)
  observed_type <- unname(manifest[1L, "matrix_type"])
  if (!identical(observed_type, expected_type)) {
    stop(
      "Pseudobulk manifest matrix_type is '", observed_type,
      "' but analysis mode '", analysis_mode,
      "' requires '", expected_type, "'.",
      call. = FALSE
    )
  }
  list(path = path, values = as.list(manifest[1L, , drop = TRUE]))
}

omix_deg_resolve_manifest <- function(
  manifest_upload,
  analysis_mode,
  code_dir,
  explicit_table_pair = FALSE
) {
  if (nzchar(manifest_upload)) {
    path <- omix_deg_validate_upload(manifest_upload, "Pseudobulk manifest", "dcf")
    return(omix_deg_read_manifest(path, analysis_mode))
  }
  # Two explicit tables are authoritative. This avoids unrelated attached
  # manifests changing an otherwise complete upload.
  if (isTRUE(explicit_table_pair)) return(NULL)
  candidates <- omix_deg_find_files(code_dir, "dcf")
  candidates <- candidates[tolower(basename(candidates)) == "pseudobulk_manifest.dcf"]
  if (length(candidates) == 0L) return(NULL)
  path <- omix_deg_unique_candidate(
    candidates,
    "pseudobulk manifest",
    "Attach one pseudobulk result bundle or select Pseudobulk Manifest explicitly."
  )
  omix_deg_read_manifest(path, analysis_mode)
}

omix_deg_manifest_file <- function(manifest, filename, label) {
  directory <- dirname(manifest$path)
  candidates <- list.files(directory, full.names = TRUE, recursive = FALSE)
  candidates <- candidates[tolower(basename(candidates)) == tolower(filename)]
  omix_deg_unique_candidate(
    candidates,
    label,
    paste0("The manifest bundle must contain ", filename, ".")
  )
}

omix_deg_matrix_candidates <- function(code_dir, analysis_mode) {
  files <- omix_deg_find_files(code_dir, c("csv", "tsv", "txt", "rds"))
  names <- basename(files)
  expected <- omix_deg_expected_matrix_filename(analysis_mode)
  exact <- files[tolower(names) == tolower(expected)]
  if (length(exact) > 0L) return(exact)
  pattern <- switch(
    analysis_mode,
    raw_counts = "count",
    harmony_mean_expression = "harmony.*mean.*expression",
    sct_mean_expression = "sct.*mean.*(log2.*)?expression"
  )
  files[grepl(pattern, names, ignore.case = TRUE)]
}

omix_deg_metadata_candidates <- function(code_dir) {
  files <- omix_deg_find_files(code_dir, c("csv", "tsv", "txt", "rds"))
  exact <- files[tolower(basename(files)) == "pseudobulk_sample_metadata.csv"]
  if (length(exact) > 0L) return(exact)
  files[grepl("metadata|sample", basename(files), ignore.case = TRUE)]
}

omix_deg_resolve_table_inputs <- function(
  matrix_upload,
  metadata_upload,
  manifest_upload,
  analysis_mode,
  code_dir
) {
  table_extensions <- c("csv", "tsv", "txt", "rds")
  explicit_pair <- nzchar(matrix_upload) && nzchar(metadata_upload)
  manifest <- omix_deg_resolve_manifest(
    manifest_upload,
    analysis_mode,
    code_dir,
    explicit_table_pair = explicit_pair
  )

  if (nzchar(matrix_upload)) {
    matrix_path <- omix_deg_validate_upload(matrix_upload, "Expression matrix", table_extensions)
  } else if (!is.null(manifest)) {
    matrix_path <- omix_deg_manifest_file(
      manifest,
      omix_deg_expected_matrix_filename(analysis_mode),
      "pseudobulk expression matrix"
    )
  } else {
    matrix_path <- omix_deg_unique_candidate(
      omix_deg_matrix_candidates(code_dir, analysis_mode),
      "expression matrix",
      paste0(
        "Attach one unambiguous ", analysis_mode,
        " matrix or select Matrix File explicitly."
      )
    )
  }

  if (nzchar(metadata_upload)) {
    metadata_path <- omix_deg_validate_upload(metadata_upload, "Sample metadata", table_extensions)
  } else if (!is.null(manifest)) {
    metadata_path <- omix_deg_manifest_file(
      manifest,
      "Pseudobulk_Sample_Metadata.csv",
      "pseudobulk sample metadata"
    )
  } else {
    metadata_path <- omix_deg_unique_candidate(
      omix_deg_metadata_candidates(code_dir),
      "sample metadata",
      "Attach one unambiguous metadata table or select Sample Metadata File explicitly."
    )
  }

  list(matrix = matrix_path, metadata = metadata_path, manifest = manifest)
}

omix_deg_find_moo <- function(upload_path, code_dir) {
  if (nzchar(upload_path)) {
    return(omix_deg_validate_upload(upload_path, "MOO", "rds"))
  }
  omix_deg_unique_candidate(
    omix_deg_find_files(code_dir, "rds"),
    "MOO .rds file",
    "Attach one MOO object or select MOO File explicitly."
  )
}
