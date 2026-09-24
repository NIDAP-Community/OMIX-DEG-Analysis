# Manual regression check against the legacy unpaired and paired templates.
#
# Run from the OMIX-DEG-Analysis repository root:
#   OMIX_TEMPLATE_LIBRARY=../template-library Rscript tests/manual/validate-legacy-equivalence.R
#
# This check deliberately disables filtering and batch adjustment so it tests
# the common TMM + voom + limma model core. Filtering and adjusted-expression
# behaviour are covered separately by the implementation tests.

template_library <- Sys.getenv("OMIX_TEMPLATE_LIBRARY", unset = "../template-library")
template_library <- normalizePath(template_library, mustWork = TRUE)

source(file.path("code", "functions", "OMIX_DEG_Analysis.R"))
source(file.path(template_library, "templates", "DEG_Analysis_v193.R"))
source(file.path(template_library, "templates", "DEG_Paired_v53.R"))

make_fixture <- function() {
  set.seed(42)
  sample_ids <- c("A1", "A2", "A3", "B1", "B2", "B3")
  counts <- matrix(stats::rpois(400L * length(sample_ids), lambda = 80), ncol = length(sample_ids))
  counts[seq_len(60L), sample_ids %in% c("B1", "B2", "B3")] <-
    counts[seq_len(60L), sample_ids %in% c("B1", "B2", "B3")] + 80
  colnames(counts) <- sample_ids
  list(
    dataset = data.frame(Gene = paste0("Gene", seq_len(nrow(counts))), counts, check.names = FALSE),
    metadata = data.frame(
      Sample = sample_ids,
      Condition = rep(c("A", "B"), each = 3L),
      stringsAsFactors = FALSE
    )
  )
}

compare_statistics <- function(new_results, reference_results, label, tolerance = 1e-8) {
  new_subset <- new_results[, c("GeneName", "B-A_logFC", "B-A_pval")]
  colnames(new_subset)[1L] <- "Gene"
  shared <- merge(
    new_subset,
    reference_results[, c("Gene", "B-A_logFC", "B-A_pval")],
    by = "Gene",
    suffixes = c("_new", "_reference")
  )
  logfc_delta <- max(abs(shared$`B-A_logFC_new` - shared$`B-A_logFC_reference`))
  pvalue_delta <- max(abs(shared$`B-A_pval_new` - shared$`B-A_pval_reference`))
  message(
    label,
    ": ", nrow(shared), " shared genes; max |logFC delta| = ", format(logfc_delta, digits = 4),
    ", max |p-value delta| = ", format(pvalue_delta, digits = 4)
  )
  if (logfc_delta > tolerance || pvalue_delta > tolerance) {
    stop(label, " regression comparison exceeded tolerance ", tolerance, ".")
  }
}

fixture <- make_fixture()
scratch <- tempfile("omix-deg-legacy-validation-")
dir.create(scratch)
original_working_directory <- getwd()
on.exit({
  setwd(original_working_directory)
  unlink(scratch, recursive = TRUE)
}, add = TRUE)

unpaired_new <- omix_deg_analysis(
  Dataset = fixture$dataset,
  Metadata_Table = fixture$metadata,
  sample_names_column = "Sample",
  samples_to_include = fixture$metadata$Sample,
  gene_names_column = "Gene",
  contrast_variable_columns = "Condition",
  contrasts = "B-A",
  filter_low_expression = FALSE,
  return_batch_corrected_values = FALSE
)
setwd(scratch)
unpaired_reference <- DEG_Analysis(
  counts_matrix = fixture$dataset,
  sample_metadata = fixture$metadata,
  feature_gene_names_column = "Gene",
  sample_names_column = "Sample",
  columns_to_include = c("Gene", fixture$metadata$Sample),
  covariates_column_s_ = "Condition",
  contrast_variable_column = "Condition",
  contrasts = "B-A",
  return_normalized_counts = TRUE,
  normalization_method = "TMM"
)
setwd(original_working_directory)
compare_statistics(unpaired_new, unpaired_reference, "Unpaired legacy equivalence")

paired_metadata <- data.frame(
  Sample = c("D1_A", "D1_B", "D2_A", "D2_B", "D3_A", "D3_B"),
  Condition = rep(c("A", "B"), 3L),
  Donor = rep(paste0("D", 1:3), each = 2L),
  stringsAsFactors = FALSE
)
paired_counts <- fixture$dataset
colnames(paired_counts)[-1L] <- paired_metadata$Sample

paired_new <- omix_deg_analysis(
  Dataset = paired_counts,
  Metadata_Table = paired_metadata,
  sample_names_column = "Sample",
  samples_to_include = paired_metadata$Sample,
  gene_names_column = "Gene",
  contrast_variable_columns = "Condition",
  contrasts = "B-A",
  donor_variable_column = "Donor",
  filter_low_expression = FALSE,
  return_batch_corrected_values = FALSE
)
setwd(scratch)
paired_reference <- deg_paired(
  counts_matrix = paired_counts,
  sample_metadata = paired_metadata,
  gene_names_column = "Gene",
  sample_names_column = "Sample",
  samples_to_include = paired_metadata$Sample,
  contrast_variable_columns = "Condition",
  contrasts = "B-A",
  donor_variable_column = "Donor",
  return_batch_corrected_values = FALSE
)
setwd(original_working_directory)
compare_statistics(paired_new, paired_reference, "Paired legacy equivalence")

message("Legacy-equivalence validation passed")
