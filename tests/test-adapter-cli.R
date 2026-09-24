required <- c("optparse", "edgeR", "limma")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Adapter CLI test requires: ", paste(missing, collapse = ", "))
}

root <- tempfile("omix-deg-cli-")
dir.create(root, recursive = TRUE)
set.seed(42)
samples <- c("A1", "A2", "A3", "B1", "B2", "B3")
metadata <- data.frame(
  Sample = samples,
  Group = rep(c("A", "B"), each = 3L),
  Batch = rep(c("Run1", "Run2"), 3L),
  stringsAsFactors = FALSE
)
counts <- matrix(stats::rpois(300L * length(samples), lambda = 60), ncol = length(samples))
counts[seq_len(40L), 4:6] <- counts[seq_len(40L), 4:6] + 70
colnames(counts) <- samples

mode_specs <- list(
  raw_counts = list(filename = "Pseudobulk_Counts.csv", matrix_type = "raw_integer_counts", values = counts),
  harmony_mean_expression = list(
    filename = "Harmony_Mean_Expression.csv",
    matrix_type = "harmony_corrected_mean_expression",
    values = log2(counts + 1)
  ),
  sct_mean_expression = list(
    filename = "SCT_Mean_Log2_Expression.csv",
    matrix_type = "sctransform_mean_log2_expression",
    values = log2(counts + 1)
  )
)

for (mode in names(mode_specs)) {
  spec <- mode_specs[[mode]]
  bundle <- file.path(root, mode)
  output <- file.path(root, paste0(mode, "-results"))
  dir.create(bundle)
  matrix_table <- data.frame(GeneName = paste0("Gene", seq_len(nrow(counts))), spec$values, check.names = FALSE)
  utils::write.csv(matrix_table, file.path(bundle, spec$filename), row.names = FALSE)
  utils::write.csv(metadata, file.path(bundle, "Pseudobulk_Sample_Metadata.csv"), row.names = FALSE)
  write.dcf(data.frame(matrix_type = spec$matrix_type), file.path(bundle, "Pseudobulk_Manifest.dcf"))

  log <- suppressWarnings(system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = c(
      file.path("code", "main.R"),
      "--input_type", "table",
      "--analysis_mode", mode,
      "--pseudobulk_manifest", file.path(bundle, "Pseudobulk_Manifest.dcf"),
      "--output_dir", output
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(log, "status")
  if (!is.null(status) && status != 0L) {
    stop("Adapter CLI failed for ", mode, ":\n", paste(log, collapse = "\n"))
  }
  stopifnot(all(file.exists(file.path(output, c(
    "DEG_Analysis.csv",
    "Sample_Metadata.csv",
    "run_summary.txt"
  )))))
  summary <- readLines(file.path(output, "run_summary.txt"))
  stopifnot(
    any(summary == paste("analysis mode:", mode)),
    any(summary == paste("input matrix type:", spec$matrix_type)),
    any(grepl("Pseudobulk_Manifest.dcf", summary, fixed = TRUE))
  )
  diagnostics <- file.path(output, c(
    "normalization_boxplots.png",
    "normalization_densities.png",
    "voom_mean_variance.png"
  ))
  stopifnot(if (identical(mode, "raw_counts")) all(file.exists(diagnostics)) else !any(file.exists(diagnostics)))
}

unlink(root, recursive = TRUE)
message("OMIX DEG adapter CLI passed raw-count, Harmony-mean, and SCT-mean bundle tests")
