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
  make_option("--counts_file", type = "character", default = ""),
  make_option("--metadata_file", type = "character", default = ""),
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

resolve_input <- function(upload_path, filename_pattern, label) {
  if (nzchar(upload_path)) {
    if (!file.exists(upload_path)) stop(label, " upload does not exist: ", upload_path)
    return(normalizePath(upload_path, mustWork = TRUE))
  }
  search_roots <- c("/data", file.path(dirname(code_dir), "data"))
  files <- unique(unlist(lapply(search_roots[dir.exists(search_roots)], function(root) {
    list.files(root, pattern = "[.](csv|tsv|txt|rds)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  })))
  candidates <- files[grepl(filename_pattern, basename(files), ignore.case = TRUE)]
  if (length(candidates) == 1L) return(candidates)
  if (length(candidates) == 0L) {
    stop("No ", label, " file was uploaded or uniquely discovered. Upload it explicitly or attach one file named with '", filename_pattern, "'.")
  }
  stop("Multiple ", label, " candidates found. Upload the intended file explicitly: ", paste(candidates, collapse = ", "))
}

counts_path <- resolve_input(opt$counts_file, "count", "counts")
metadata_path <- resolve_input(opt$metadata_file, "metadata|sample", "metadata")
counts <- read_tabular_input(counts_path)
metadata <- read_tabular_input(metadata_path)

selected_samples <- split_csv(opt$samples_to_include)
if (length(selected_samples) == 0L) {
  if (!opt$sample_names_column %in% colnames(metadata)) {
    stop("Sample ID column '", opt$sample_names_column, "' was not found in metadata.")
  }
  selected_samples <- intersect(as.character(metadata[[opt$sample_names_column]]), colnames(counts))
}
if (length(selected_samples) == 0L) stop("No selected metadata sample IDs occur in the count table.")

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
  return_batch_corrected_values = as_logical(opt$return_batch_corrected_values, "return_batch_corrected_values"),
  remove_donor_effect_for_downstream = as_logical(opt$remove_donor_effect_for_downstream, "remove_donor_effect_for_downstream"),
  summarization_method = opt$summarization_method
)

output_dir <- opt$output_dir
if (!nzchar(output_dir)) output_dir <- if (dir.exists("/results")) "/results" else file.path(dirname(code_dir), "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, file.path(output_dir, "DEG_Analysis.csv"), row.names = FALSE, na = "")
run_summary <- attr(results, "omix_deg_run")
writeLines(c(
  "OMIX DEG Analysis run summary",
  paste("counts input:", counts_path),
  paste("metadata input:", metadata_path),
  paste("model type:", run_summary$model_type),
  paste("genes before filtering:", run_summary$genes_before_filtering),
  paste("genes modelled:", run_summary$genes_modelled),
  paste("expression output:", run_summary$expression_output),
  paste("adjusted columns:", paste(run_summary$adjusted_columns, collapse = ", "))
), file.path(output_dir, "DEG_Analysis_run_summary.txt"))
message("Wrote results to ", output_dir)
