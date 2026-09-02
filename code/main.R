#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

script_path <- function() {
  candidates <- commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))]
  if (length(candidates) == 0L) return(getwd())
  dirname(normalizePath(sub("^--file=", "", candidates[[1L]]), mustWork = TRUE))
}

split_csv <- function(value) {
  if (is.null(value) || !nzchar(trimws(value))) return(character())
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

as_logical <- function(value, argument) {
  normalized <- tolower(trimws(value))
  if (normalized %in% c("true", "t", "1", "yes")) return(TRUE)
  if (normalized %in% c("false", "f", "0", "no")) return(FALSE)
  stop("--", argument, " must be true or false.")
}

option_list <- list(
  make_option("--input_type", type = "character", default = "auto"),
  make_option("--counts_file", type = "character", default = ""),
  make_option("--metadata_file", type = "character", default = ""),
  make_option("--moo_file", type = "character", default = ""),
  make_option("--gene_names_column", type = "character", default = "GeneName"),
  make_option("--sample_names_column", type = "character", default = "Sample"),
  make_option("--samples_to_include", type = "character", default = ""),
  make_option("--contrast_variable_columns", type = "character", default = "Group"),
  make_option("--contrasts", type = "character", default = "B-A"),
  make_option("--covariate_columns", type = "character", default = ""),
  make_option("--batch_effect_columns", type = "character", default = "Batch"),
  make_option("--donor_variable_column", type = "character", default = ""),
  make_option("--filter_low_expression", type = "character", default = "true"),
  make_option("--return_batch_corrected_values", type = "character", default = "true"),
  make_option("--remove_donor_effect_for_downstream", type = "character", default = "true"),
  make_option("--normalization_method", type = "character", default = "TMM"),
  make_option("--write_normalization_diagnostics", type = "character", default = "true"),
  make_option("--summarization_method", type = "character", default = "sum"),
  make_option("--output_dir", type = "character", default = "")
)

opt <- parse_args(OptionParser(option_list = option_list))
code_dir <- script_path()
source(file.path(code_dir, "functions", "OMIX_DEG_Analysis.R"))

read_tabular_input <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (extension == "rds") {
    object <- readRDS(path)
    if (is.matrix(object)) return(as.data.frame(object, check.names = FALSE))
    if (is.data.frame(object)) return(object)
    stop("RDS input must contain a data frame or matrix: ", path)
  }
  if (extension == "csv") {
    return(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (extension %in% c("tsv", "txt")) {
    return(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  stop("Unsupported input extension for ", path, ". Use CSV, TSV, TXT, or RDS.")
}

canonical_data_path <- function(path) {
  # Code Ocean can expose the same mounted data directory as both /data and
  # //data. Collapse repeated leading slashes so one physical file is not
  # treated as two discovered inputs.
  path <- normalizePath(path, mustWork = TRUE, winslash = "/")
  sub("^/{2,}", "/", path)
}

find_tabular_candidates <- function(filename_pattern) {
  search_roots <- c("/data", file.path(dirname(code_dir), "data"))
  files <- unique(unlist(lapply(search_roots[dir.exists(search_roots)], function(root) {
    list.files(root, pattern = "[.](csv|tsv|txt|rds)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  })))
  files <- unique(vapply(files, canonical_data_path, character(1)))
  files[grepl(filename_pattern, basename(files), ignore.case = TRUE)]
}

resolve_input <- function(upload_path, filename_pattern, label) {
  if (nzchar(upload_path)) {
    if (!file.exists(upload_path)) stop(label, " upload does not exist: ", upload_path)
    return(canonical_data_path(upload_path))
  }
  candidates <- find_tabular_candidates(filename_pattern)
  if (length(candidates) == 1L) return(candidates)
  if (length(candidates) == 0L) {
    stop("No ", label, " file was uploaded or uniquely discovered. Upload it explicitly or attach one file named with '", filename_pattern, "'.")
  }
  stop("Multiple ", label, " candidates found. Upload the intended file explicitly: ", paste(candidates, collapse = ", "))
}

find_moo_input <- function(upload_path) {
  if (nzchar(upload_path)) {
    if (!file.exists(upload_path)) stop("MOO upload does not exist: ", upload_path)
    if (tolower(tools::file_ext(upload_path)) != "rds") {
      stop("MOO input must be an .rds file: ", upload_path)
    }
    return(normalizePath(upload_path, mustWork = TRUE))
  }
  search_roots <- c("/data", file.path(dirname(code_dir), "data"))
  candidates <- unique(unlist(lapply(search_roots[dir.exists(search_roots)], function(root) {
    list.files(root, pattern = "[.]rds$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  })))
  if (length(candidates) == 1L) return(candidates)
  if (length(candidates) == 0L) {
    stop("No MOO .rds file was uploaded or uniquely discovered. Upload it explicitly.")
  }
  stop("Multiple .rds files were found. Upload the intended MOO explicitly: ", paste(candidates, collapse = ", "))
}

input_type <- tolower(trimws(opt$input_type))
if (!input_type %in% c("auto", "table", "moo")) {
  stop("--input_type must be one of: auto, table, moo.")
}
if (nzchar(opt$moo_file) && (nzchar(opt$counts_file) || nzchar(opt$metadata_file))) {
  stop("Use either --moo_file or --counts_file/--metadata_file, not both.")
}

has_explicit_table_input <- nzchar(opt$counts_file) || nzchar(opt$metadata_file)
use_moo <- identical(input_type, "moo") || nzchar(opt$moo_file)
if (identical(input_type, "auto") && !use_moo && !has_explicit_table_input) {
  count_candidates <- find_tabular_candidates("count")
  metadata_candidates <- find_tabular_candidates("metadata|sample")
  use_moo <- length(count_candidates) == 0L && length(metadata_candidates) == 0L
}

if (use_moo) {
  if (!requireNamespace("OmixMOSuite", quietly = TRUE)) {
    stop(
      "MOO input requires the optional OmixMOSuite bridge and MOSuite runtime. ",
      "Use table input or install the supported bridge before running this mode."
    )
  }
  moo_path <- find_moo_input(opt$moo_file)
  omix_input <- OmixMOSuite::omix_read_mosuite_rds(moo_path, count_type = "raw")
  counts <- omix_input$counts
  metadata <- omix_input$metadata
  counts_path <- paste0("MOO raw count layer: ", moo_path)
  metadata_path <- paste0("MOO embedded sample metadata: ", moo_path)
  opt$gene_names_column <- omix_input$feature_id_column
  opt$sample_names_column <- omix_input$sample_id_column
  message("MOO input: extracted raw counts and aligned sample metadata from ", moo_path)
} else {
  if (identical(input_type, "moo")) stop("MOO input mode requires a MOO .rds file.")
  counts_path <- resolve_input(opt$counts_file, "count", "counts")
  metadata_path <- resolve_input(opt$metadata_file, "metadata|sample", "metadata")
  counts <- read_tabular_input(counts_path)
  metadata <- read_tabular_input(metadata_path)
}

selected_samples <- split_csv(opt$samples_to_include)
if (length(selected_samples) == 0L) {
  if (!opt$sample_names_column %in% colnames(metadata)) {
    stop("Sample ID column '", opt$sample_names_column, "' was not found in metadata.")
  }
  selected_samples <- intersect(as.character(metadata[[opt$sample_names_column]]), colnames(counts))
}
if (length(selected_samples) == 0L) stop("No selected metadata sample IDs occur in the count table.")

output_dir <- opt$output_dir
if (!nzchar(output_dir)) output_dir <- if (dir.exists("/results")) "/results" else file.path(dirname(code_dir), "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- omix_deg_analysis(
  Dataset = counts,
  Metadata_Table = metadata,
  sample_names_column = opt$sample_names_column,
  samples_to_include = selected_samples,
  gene_names_column = opt$gene_names_column,
  contrast_variable_columns = split_csv(opt$contrast_variable_columns),
  contrasts = split_csv(opt$contrasts),
  covariate_columns = split_csv(opt$covariate_columns),
  donor_variable_column = split_csv(opt$donor_variable_column),
  batch_effect_columns = split_csv(opt$batch_effect_columns),
  filter_low_expression = as_logical(opt$filter_low_expression, "filter_low_expression"),
  normalization_method = opt$normalization_method,
  normalization_diagnostics = as_logical(opt$write_normalization_diagnostics, "write_normalization_diagnostics"),
  diagnostics_output_dir = output_dir,
  return_batch_corrected_values = as_logical(opt$return_batch_corrected_values, "return_batch_corrected_values"),
  remove_donor_effect_for_downstream = as_logical(opt$remove_donor_effect_for_downstream, "remove_donor_effect_for_downstream"),
  summarization_method = opt$summarization_method
)

utils::write.csv(results, file.path(output_dir, "DEG_Analysis.csv"), row.names = FALSE, na = "")
output_samples <- selected_samples[selected_samples %in% names(results)]
if (length(output_samples) == 0L) {
  stop("The DEG result does not contain any modeled sample-expression columns for metadata export.")
}
output_metadata <- metadata[
  match(output_samples, as.character(metadata[[opt$sample_names_column]])),
  ,
  drop = FALSE
]
utils::write.csv(
  output_metadata,
  file.path(output_dir, "Sample_Metadata.csv"),
  row.names = FALSE,
  na = ""
)
run_summary <- attr(results, "omix_deg_run")
writeLines(c(
  "OMIX DEG Analysis run summary",
  paste("counts input:", counts_path),
  paste("metadata input:", metadata_path),
  paste("model type:", run_summary$model_type),
  paste("normalization profile:", opt$normalization_method),
  paste("library-size normalization:", run_summary$library_size_normalization),
  paste("voom-scale normalization:", run_summary$voom_scale_normalization),
  paste("normalization diagnostics:", if (length(run_summary$normalization_diagnostic_files) == 0L) "not written" else paste(basename(run_summary$normalization_diagnostic_files), collapse = ", ")),
  paste("genes before filtering:", run_summary$genes_before_filtering),
  paste("genes modelled:", run_summary$genes_modelled),
  paste("expression output:", run_summary$expression_output),
  paste("adjusted columns:", paste(run_summary$adjusted_columns, collapse = ", "))
), file.path(output_dir, "DEG_Analysis_run_summary.txt"))
message("Wrote results to ", output_dir)
