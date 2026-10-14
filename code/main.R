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

as_auto_logical <- function(value, argument) {
  if (is.null(value) || !nzchar(trimws(value)) || identical(tolower(trimws(value)), "auto")) {
    return(NULL)
  }
  as_logical(value, argument)
}

read_tabular_input <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (extension == "rds") {
    object <- readRDS(path)
    if (is.matrix(object)) return(as.data.frame(object, check.names = FALSE))
    if (is.data.frame(object)) return(object)
    stop("RDS table input must contain a data frame or matrix: ", path)
  }
  if (extension == "csv") {
    return(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (extension %in% c("tsv", "txt")) {
    return(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  stop("Unsupported table input extension for ", path, ". Use CSV, TSV, TXT, or RDS.")
}

option_list <- list(
  make_option("--counts", type = "character", default = "", help = "Optional explicit matrix table; otherwise discover from attached data"),
  make_option("--metadata", type = "character", default = "", help = "Optional explicit sample-metadata table"),
  make_option("--pseudobulk_manifest", type = "character", default = "", help = "Optional explicit Pseudobulk_Manifest.dcf; otherwise discover exactly one"),
  make_option("--moo", type = "character", default = "", help = "Optional compatible MOO .rds; raw_counts only"),
  make_option("--input_type", type = "character", default = "table", help = "table or moo [default: %default]"),
  make_option("--analysis_mode", type = "character", default = "raw_counts", help = "raw_counts, harmony_mean_expression, or sct_mean_expression [default: %default]"),
  make_option("--gene_names_column", type = "character", default = "GeneName", help = "Matrix feature-ID column [default: %default]"),
  make_option("--sample_names_column", type = "character", default = "Sample", help = "Metadata sample-ID column [default: %default]"),
  make_option("--samples_to_include", type = "character", default = "", help = "Optional comma-separated sample IDs"),
  make_option("--contrast_variable_columns", type = "character", default = "Group", help = "One or two comma-separated biological group columns [default: %default]"),
  make_option("--contrasts", type = "character", default = "B-A", help = "Comma-separated numerator-denominator contrasts [default: %default]"),
  make_option("--covariate_columns", type = "character", default = "", help = "Optional comma-separated fixed covariates"),
  make_option("--batch_effect_columns", type = "character", default = "auto", help = "auto uses Batch for raw_counts only when metadata contains it; otherwise supply blank or comma-separated technical batch columns [default: %default]"),
  make_option("--donor_variable_column", type = "character", default = "", help = "Optional donor/patient column for repeated measures"),
  make_option("--filter_low_expression", type = "character", default = "auto", help = "auto, true, or false [default: %default]"),
  make_option("--return_batch_corrected_values", type = "character", default = "auto", help = "auto, true, or false [default: %default]"),
  make_option("--remove_donor_effect_for_downstream", type = "character", default = "auto", help = "auto, true, or false [default: %default]"),
  make_option("--normalization_method", type = "character", default = "auto", help = "Mode-aware normalization profile [default: %default]"),
  make_option("--write_normalization_diagnostics", type = "character", default = "auto", help = "auto, true, or false [default: %default]"),
  make_option("--summarization_method", type = "character", default = "auto", help = "auto, sum, mean, or max [default: %default]"),
  make_option("--output_dir", type = "character", default = "", help = "Output directory; defaults to /results in Code Ocean")
)

opt <- parse_args(OptionParser(option_list = option_list))
code_dir <- script_path()
source(file.path(code_dir, "functions", "OMIX_DEG_Analysis.R"))
source(file.path(code_dir, "adapter", "DEG_Adapter.R"))

input_type <- tolower(trimws(opt$input_type))
if (!input_type %in% c("table", "moo")) {
  stop("--input_type must be one of: table, moo.")
}
analysis_mode <- tolower(trimws(opt$analysis_mode))
if (!analysis_mode %in% c("raw_counts", "harmony_mean_expression", "sct_mean_expression")) {
  stop("--analysis_mode must be raw_counts, harmony_mean_expression, or sct_mean_expression.")
}
if (nzchar(opt$moo) && (nzchar(opt$counts) || nzchar(opt$metadata) || nzchar(opt$pseudobulk_manifest))) {
  stop("Use either --moo or table/manifest inputs, not both.")
}
if (identical(input_type, "moo") && (nzchar(opt$counts) || nzchar(opt$metadata) || nzchar(opt$pseudobulk_manifest))) {
  stop("--counts, --metadata, and --pseudobulk_manifest are only valid with --input_type table.")
}

use_moo <- identical(input_type, "moo")

if (use_moo) {
  if (!identical(analysis_mode, "raw_counts")) {
    stop("MOO input only supports --analysis_mode raw_counts.")
  }
  if (!requireNamespace("OmixMOSuite", quietly = TRUE)) {
    stop(
      "MOO input requires the optional OmixMOSuite bridge and its supported object runtime. ",
      "Use table input or install the bridge before running this mode."
    )
  }
  moo_path <- omix_deg_find_moo(opt$moo, code_dir)
  omix_input <- OmixMOSuite::omix_read_mosuite_rds(moo_path, count_type = "raw")
  counts <- omix_input$counts
  metadata <- omix_input$metadata
  opt$gene_names_column <- omix_input$feature_id_column
  opt$sample_names_column <- omix_input$sample_id_column
  input_summary <- c(
    "input type: MOO",
    paste("MOO file:", moo_path),
    paste("MOO count layer:", omix_input$provenance$count_type),
    "pseudobulk manifest: not supplied"
  )
  message("MOO input: extracted raw counts and aligned sample metadata from ", moo_path)
} else {
  if (nzchar(opt$moo)) stop("--moo is only valid with --input_type moo.")
  resolved <- omix_deg_resolve_table_inputs(
    matrix_upload = opt$counts,
    metadata_upload = opt$metadata,
    manifest_upload = opt$pseudobulk_manifest,
    analysis_mode = analysis_mode,
    code_dir = code_dir
  )
  counts <- read_tabular_input(resolved$matrix)
  metadata <- read_tabular_input(resolved$metadata)
  input_summary <- c(
    "input type: table",
    paste("matrix input:", resolved$matrix),
    paste("metadata input:", resolved$metadata),
    if (is.null(resolved$manifest)) {
      "pseudobulk manifest: not supplied"
    } else {
      paste("pseudobulk manifest:", resolved$manifest$path)
    }
  )
}

selected_samples <- split_csv(opt$samples_to_include)
if (length(selected_samples) == 0L) {
  if (!opt$sample_names_column %in% colnames(metadata)) {
    stop("Sample ID column '", opt$sample_names_column, "' was not found in metadata.")
  }
  selected_samples <- intersect(as.character(metadata[[opt$sample_names_column]]), colnames(counts))
}
if (length(selected_samples) == 0L) stop("No selected metadata sample IDs occur in the expression matrix.")

output_dir <- opt$output_dir
if (!nzchar(output_dir)) {
  output_dir <- if (dir.exists("/results")) "/results" else file.path(dirname(code_dir), "results")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

batch_effect_columns <- if (identical(tolower(trimws(opt$batch_effect_columns)), "auto")) {
  .omix_deg_auto_batch_effect_columns(metadata, analysis_mode)
} else {
  split_csv(opt$batch_effect_columns)
}
normalization_diagnostics <- if (
  identical(tolower(trimws(opt$write_normalization_diagnostics)), "auto") &&
    identical(analysis_mode, "raw_counts")
) {
  TRUE
} else {
  as_auto_logical(opt$write_normalization_diagnostics, "write_normalization_diagnostics")
}

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
  batch_effect_columns = batch_effect_columns,
  analysis_mode = analysis_mode,
  filter_low_expression = as_auto_logical(opt$filter_low_expression, "filter_low_expression"),
  normalization_method = opt$normalization_method,
  normalization_diagnostics = normalization_diagnostics,
  diagnostics_output_dir = output_dir,
  return_batch_corrected_values = as_auto_logical(
    opt$return_batch_corrected_values,
    "return_batch_corrected_values"
  ),
  remove_donor_effect_for_downstream = as_auto_logical(
    opt$remove_donor_effect_for_downstream,
    "remove_donor_effect_for_downstream"
  ),
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
  input_summary,
  paste("analysis mode:", run_summary$analysis_mode),
  paste("input matrix type:", run_summary$input_matrix_type),
  paste("model type:", run_summary$model_type),
  paste("normalization profile:", run_summary$normalization_method),
  paste("library-size normalization:", run_summary$library_size_normalization),
  paste("voom-scale normalization:", run_summary$voom_scale_normalization),
  paste(
    "voom trend span:",
    if (identical(run_summary$analysis_mode, "raw_counts")) run_summary$voom_trend_span else "not applicable"
  ),
  paste(
    "voom adaptive span:",
    if (identical(run_summary$analysis_mode, "raw_counts")) {
      tolower(as.character(run_summary$voom_adaptive_span))
    } else {
      "not applicable"
    }
  ),
  paste(
    "normalization diagnostics:",
    if (length(run_summary$normalization_diagnostic_files) == 0L) {
      "not written"
    } else {
      paste(basename(run_summary$normalization_diagnostic_files), collapse = ", ")
    }
  ),
  paste("genes before filtering:", run_summary$genes_before_filtering),
  paste("genes modelled:", run_summary$genes_modelled),
  paste("expression output:", run_summary$expression_output),
  paste("adjusted columns:", paste(run_summary$adjusted_columns, collapse = ", "))
), file.path(output_dir, "run_summary.txt"))
message("Wrote results to ", normalizePath(output_dir, mustWork = TRUE))
