#' DEG Analysis for Paired Samples (Linear Mixed Model)
#'
#' @description
#' Performs Differential Expression of Genes (DEG) analysis for repeated
#' measures using a linear mixed model. This approach accounts for the
#' correlation between multiple samples from the same individual, recognizing
#' that they are not independent. It leverages the
#' \code{\link[limma]{duplicateCorrelation}}
#' function from the \pkg{limma} package to model these correlations
#' effectively.
#'
#' By default, donor- and batch-corrected values are returned to ensure
#' accurate results. A Donor variable (e.g., Patient ID) is **required**.
#'
#' @details
#' **Instructions for Data Preparation**
#'
#' **1. Handling the Donor Variable**
#' - Do **not** include the Donor variable in \code{covariate_columns}.
#' - It is handled separately via \code{donor_variable_column} as a random
#' effect.
#'
#' **2. Adding Contrast Variables**
#' - Up to **two** contrast variables are supported in
#' \code{contrast_variable_columns}.
#' - Two contrasts: combine as \code{factor1.factor2}
#' (e.g., \code{Treated.Day3 - Treated.Day1}).
#' - One contrast: use the format \code{Day3 - Day1}.
#' - Do **not** include contrast variable(s) in \code{covariate_columns}.
#'
#' **3. Defining Covariate Columns**
#' - Include all variables **except** the Donor and Contrast variables.
#' - May be left as \code{NULL} if no additional covariates are needed.
#'
#' **4. Batch Effect Removal**
#' - Enabled by default (\code{return_batch_corrected_values = TRUE}), but
#' optional and performed separately from the main statistical analysis.
#' - Do **not** include the Donor variable in \code{batch_column}; it is
#' handled automatically.
#'
#' This template is a Beta version undergoing active development. Contact
#' CCBR at \email{NCICCBRNIDAP@@mail.nih.gov} if you encounter problems.
#'
#' @param counts_matrix A data frame containing raw (or filtered) counts. Must
#' have one column for gene names and the remaining columns containing numeric
#' counts per sample.
#' @param sample_metadata A data frame containing sample group information.
#' Must include a sample names column and all grouping/covariate variables
#' used in the analysis.
#' @param gene_names_column
#' Character. Name of the column in \code{Counts_Matrix}
#' containing gene identifiers.
#' @param sample_names_column Character. Name of the column in
#' \code{Sample_Metadata} containing sample identifiers.
#' @param samples_to_include Character vector. Names of the sample columns from
#' \code{Counts_Matrix} to include in the analysis. Only numeric columns
#' should be selected.
#' @param contrast_variable_columns Character vector (length 1 or 2). Column
#' name(s) in \code{Sample_Metadata} containing the grouping variable(s) for
#' differential expression contrasts.
#' @param contrasts Character vector. Contrast specifications in the format
#' \code{"group1-group2"} (e.g., \code{"treated-control"}). For two-factor
#' analyses, use dot notation: \code{"Day10.40mg-Day3.20mg"}.
#' @param donor_variable_column Character. Name of the column in
#' \code{Sample_Metadata} specifying the donor/subject ID. Used as a random
#' effect via \code{duplicateCorrelation}. Required.
#' @param covariate_columns Character vector or \code{NULL}. Column names in
#' \code{Sample_Metadata} representing additional fixed-effect covariates.
#' Do **not** include the Donor or Contrast variable(s) here. Can be
#' \code{NULL}.
#' @param batch_column Character or \code{NULL}. Column name in
#' \code{Sample_Metadata} representing a technical batch variable to remove
#' (in addition to donor) during batch correction. Can be \code{NULL}.
#' @param input_in_log_counts
#' Logical. Set to \code{TRUE} if \code{Counts_Matrix}
#' is already in log counts; counts will be back-transformed before creating
#' the \code{DGEList}. Default: \code{FALSE}.
#' @param normalization_method Character. Normalization method applied to logCPM
#' values by \code{\link[limma]{voom}}. One of \code{"none"}, \code{"scale"},
#' \code{"quantile"}, \code{"cyclicloess"}, \code{"TMM"},
#' \code{"TMMwzp"}, \code{"RLE"}, or \code{"upperquartile"}.
#' Default: \code{"TMM"}.
#' @param return_batch_corrected_values Logical. If \code{TRUE}, returns
#' donor- and batch-corrected expression values alongside statistics. If
#' \code{FALSE}, returns normalized (voom) expression values. Default:
#' \code{TRUE}.
#' @param fold_change_threshold Numeric. Fold-change threshold used when
#' tabulating the number of differentially expressed genes. Default:
#' \code{1.2}.
#' @param first_pvalue_threshold Numeric. First p-value threshold for the
#' DEG summary table. Default: \code{0.05}.
#' @param second_pvalue_threshold Numeric. Second (more stringent) p-value
#' threshold for the DEG summary table. Default: \code{0.01}.
#' @param row_number Integer. Maximum number of rows per printed table chunk.
#' Default: \code{30}.
#' @param column_number Integer. Maximum number of columns per printed table
#' chunk. Default: \code{10}.
#' @param display_warnings Integer. Passed to \code{options(warn = ...)}.
#' Use \code{0} to display warnings; \code{-1} (default) suppresses them.
#' Default: \code{-1}.
#'
#' @return A data frame of differential expression results for the
#' paired/repeated-measures design, including fold change, p-value,
#' and adjusted p-value for each contrast.
#'
#' @importFrom limma .
#' @importFrom edgeR .
#' @importFrom dplyr .
#' @importFrom stringr .
#' @importFrom tidyr .
#' @param output_file_path Character or \code{NULL}. Optional output path for
#' writing the DEG result table as a TSV file.
#'
#' @export
deg_paired <- function(
  counts_matrix,
  sample_metadata,
  gene_names_column,
  sample_names_column,
  samples_to_include,
  contrast_variable_columns,
  contrasts,
  donor_variable_column,
  covariate_columns = NULL,
  batch_column = NULL,
  input_in_log_counts = FALSE,
  normalization_method = "TMM",
  return_batch_corrected_values = TRUE,
  fold_change_threshold = 1.2,
  first_pvalue_threshold = 0.05,
  second_pvalue_threshold = 0.01,
  row_number = 30,
  column_number = 10,
  display_warnings = -1,
  output_file_path = NULL
) {
  options(warn = display_warnings)

  ## --------- ##
  ## Libraries ##
  ## --------- ##

  library(limma)
  library(edgeR)
  library(ggplot2)
  library(stringr)
  library(gridExtra)
  library(reshape2)
  library(grid)
  library(gtable)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(plotly)

  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##

  if (is.character(counts_matrix) && length(counts_matrix) == 1) {
    if (!file.exists(counts_matrix)) {
      stop("ERROR: `counts_matrix` file does not exist: ", counts_matrix)
    }
    counts_matrix <- read.delim(
      counts_matrix,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (is.character(sample_metadata) && length(sample_metadata) == 1) {
    if (!file.exists(sample_metadata)) {
      stop("ERROR: `sample_metadata` file does not exist: ", sample_metadata)
    }
    sample_metadata <- read.delim(
      sample_metadata,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (!is.data.frame(counts_matrix)) {
    stop("ERROR: `counts_matrix` must be a data frame or a valid file path.")
  }

  if (!is.data.frame(sample_metadata)) {
    stop("ERROR: `sample_metadata` must be a data frame or a valid file path.")
  }

  infer_sample_names_column <- function(counts_df, metadata_df) {
    sample_columns <- colnames(counts_df)[vapply(
      counts_df,
      is.numeric,
      logical(1)
    )]
    candidates <- colnames(metadata_df)[vapply(
      metadata_df,
      function(x) is.character(x) || is.factor(x),
      logical(1)
    )]

    if (length(candidates) == 0) {
      stop("ERROR: Could not infer sample names column from metadata.")
    }

    scores <- vapply(
      candidates,
      function(col_name) {
        sum(as.character(metadata_df[[col_name]]) %in% sample_columns)
      },
      numeric(1)
    )

    best_index <- which.max(scores)
    if (scores[best_index] == 0) {
      stop(
        "ERROR: Could not match metadata sample names to count matrix columns."
      )
    }

    candidates[best_index]
  }

  if (
    is.null(sample_names_column) ||
      !nzchar(sample_names_column) ||
      !sample_names_column %in% colnames(sample_metadata)
  ) {
    sample_names_column <- infer_sample_names_column(
      counts_matrix,
      sample_metadata
    )
    message(paste0("Using inferred sample names column: ", sample_names_column))
  }

  if (sample_names_column %in% colnames(sample_metadata)) {
    numeric_sample_columns <- colnames(counts_matrix)[vapply(
      counts_matrix,
      is.numeric,
      logical(1)
    )]
    provided_overlap <- sum(
      as.character(sample_metadata[[sample_names_column]]) %in%
        numeric_sample_columns
    )
    inferred_column <- infer_sample_names_column(counts_matrix, sample_metadata)
    inferred_overlap <- sum(
      as.character(sample_metadata[[inferred_column]]) %in%
        numeric_sample_columns
    )

    if (
      provided_overlap == 0 &&
        inferred_overlap > 0 &&
        inferred_column != sample_names_column
    ) {
      message(paste0(
        "Provided sample names column '",
        sample_names_column,
        "' has no overlap with count sample names; using '",
        inferred_column,
        "' instead."
      ))
      sample_names_column <- inferred_column
    }
  }

  if (!(gene_names_column %in% colnames(counts_matrix))) {
    fallback_feature_candidates <- c(
      "Gene",
      "GeneID",
      "gene",
      "gene_id",
      "feature_id"
    )
    fallback_feature_column <- fallback_feature_candidates[
      fallback_feature_candidates %in% colnames(counts_matrix)
    ][1]

    if (is.na(fallback_feature_column)) {
      non_numeric_columns <- colnames(counts_matrix)[vapply(
        counts_matrix,
        function(x) !(is.numeric(x) || is.integer(x)),
        logical(1)
      )]
      if (length(non_numeric_columns) > 0) {
        fallback_feature_column <- non_numeric_columns[1]
      }
    }

    if (!is.na(fallback_feature_column) && nzchar(fallback_feature_column)) {
      message(paste0(
        "Provided gene_names_column '",
        gene_names_column,
        "' not found; using '",
        fallback_feature_column,
        "' instead."
      ))
      gene_names_column <- fallback_feature_column
    } else {
      stop("ERROR: Could not determine a valid gene names column.")
    }
  }

  if (is.null(samples_to_include) || length(samples_to_include) == 0) {
    samples_to_include <- colnames(counts_matrix)[vapply(
      counts_matrix,
      is.numeric,
      logical(1)
    )]
  }

  if (
    is.null(contrast_variable_columns) || length(contrast_variable_columns) == 0
  ) {
    contrast_variable_columns <- "condition"
  }

  contrast_variable_columns <- contrast_variable_columns[
    contrast_variable_columns %in% colnames(sample_metadata)
  ]
  if (length(contrast_variable_columns) == 0) {
    stop(
      "ERROR: No valid `contrast_variable_columns` found in sample metadata."
    )
  }

  if (
    is.null(donor_variable_column) ||
      !nzchar(donor_variable_column) ||
      !donor_variable_column %in% colnames(sample_metadata)
  ) {
    donor_candidates <- grep(
      "donor|patient|subject",
      colnames(sample_metadata),
      ignore.case = TRUE,
      value = TRUE
    )
    if (length(donor_candidates) == 0) {
      stop("ERROR: Could not infer `donor_variable_column`.")
    }
    donor_variable_column <- donor_candidates[1]
    message(paste0("Using donor variable column: ", donor_variable_column))
  }

  if (is.null(covariate_columns)) {
    covariate_columns <- character(0)
  }
  covariate_columns <- covariate_columns[
    covariate_columns %in% colnames(sample_metadata)
  ]

  if (is.null(contrasts) || length(contrasts) == 0 || all(!nzchar(contrasts))) {
    if (length(contrast_variable_columns) == 1) {
      contrast_levels <- unique(
        as.character(sample_metadata[[contrast_variable_columns[1]]])
      )
      contrast_levels <- contrast_levels[
        !is.na(contrast_levels) & contrast_levels != "NA"
      ]
      if (length(contrast_levels) >= 2) {
        contrasts <- paste0(contrast_levels[2], "-", contrast_levels[1])
      } else {
        stop("ERROR: Could not infer `contrasts` from metadata.")
      }
    } else {
      stop(
        paste0(
          "ERROR: Please provide `contrasts` when using multiple ",
          "contrast variables."
        )
      )
    }
  }

  if (is.null(batch_column) || length(batch_column) == 0) {
    batch_column <- NULL
  } else {
    batch_column <- batch_column[1]
    if (!batch_column %in% colnames(sample_metadata)) {
      message(paste0(
        "Provided batch_column '",
        batch_column,
        "' not found; skipping additional batch removal."
      ))
      batch_column <- NULL
    }

    if (!is.null(batch_column)) {
      donor_batch_summary <- sample_metadata %>%
        dplyr::filter(
          !is.na(.data[[donor_variable_column]]),
          !is.na(.data[[batch_column]])
        ) %>%
        dplyr::summarise(
          n_batches = dplyr::n_distinct(.data[[batch_column]]),
          .by = .data[[donor_variable_column]]
        )

      if (
        nrow(donor_batch_summary) > 0 &&
          all(donor_batch_summary$n_batches == 1)
      ) {
        message(paste0(
          "Batch column '",
          batch_column,
          "' is confounded with donor column '",
          donor_variable_column,
          paste0(
            "' (each donor appears in only one batch); skipping ",
            "additional batch removal."
          )
        ))
        batch_column <- NULL
      }
    }
  }

  ## --------------- ##
  ## Error Messages  ##
  ## --------------- ##

  if (sum(make.names(colnames(counts_matrix)) != colnames(counts_matrix)) > 0) {
    print("Error: The following counts matrix column names are not valid:\n")
    print(colnames(counts_matrix)[
      make.names(colnames(counts_matrix)) != colnames(counts_matrix)
    ])
    print(
      paste0(
        "Likely causes are columns starting with numbers or other ",
        "special characters eg spaces."
      )
    )
    stop("Bad column names.")
  }

  ## --------- ##
  ## Functions ##
  ## --------- ##

  plot_sankey <- function(
    sample_metadata,
    contrast_variable_columns,
    covariate_columns = NULL
  ) {
    # Prepare the data for the Sankey plot
    all_columns <- c(contrast_variable_columns, covariate_columns)

    if (length(all_columns) == 1) {
      message(
        paste0(
          "Cannot draw a Sankey plot with only one contrast variable ",
          "column and no covariate columns."
        )
      )
      return(NULL)
    } else {
      # Transform the data into a long format
      long_data <- sample_metadata %>%
        pivot_longer(
          cols = all_of(all_columns),
          names_to = "Category",
          values_to = "Value"
        )

      # Create nodes
      nodes <- data.frame(name = unique(long_data$Value))

      # Create links dynamically
      links <- NULL
      for (i in 1:(length(all_columns) - 1)) {
        temp_links <- sample_metadata %>%
          dplyr::group_by(across(all_of(all_columns[c(i, i + 1)]))) %>%
          dplyr::summarise(value = dplyr::n(), .groups = "drop") %>%
          dplyr::mutate(
            source = match(.data[[all_columns[i]]], nodes$name) - 1,
            target = match(.data[[all_columns[i + 1]]], nodes$name) - 1
          ) %>%
          dplyr::select(source, target, value)
        links <- bind_rows(links, temp_links)
      }

      # Create the plotly Sankey plot
      sankey <- plot_ly(
        type = "sankey",
        orientation = "h",
        node = list(
          label = nodes$name,
          pad = 15,
          thickness = 20,
          line = list(color = "black", width = 0.5)
        ),
        link = list(
          source = links$source,
          target = links$target,
          value = links$value
        )
      ) %>%
        layout(
          title = "Sankey Diagram",
          font = list(size = 10)
        )
      return(sankey)
    }
  }

  plot_sankey_with_confounders <- function(
    data,
    contrast_columns,
    covariate_column,
    confounder_names
  ) {
    # Transform data into long format
    long_data <- data %>%
      pivot_longer(
        cols = c(contrast_columns, covariate_column),
        names_to = "Category",
        values_to = "Value"
      )

    # Create nodes
    nodes <- data.frame(name = unique(long_data$Value))

    # Create links
    long_data <- long_data %>%
      dplyr::group_by(.data[[sample_names_column]]) %>%
      dplyr::arrange(Category) %>%
      dplyr::mutate(Source = lag(Value), Target = Value) %>%
      dplyr::filter(!is.na(Source)) %>%
      dplyr::ungroup()

    links <- long_data %>%
      dplyr::group_by(Source, Target) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(
        source = match(Source, nodes$name) - 1,
        target = match(Target, nodes$name) - 1
      ) %>%
      dplyr::select(source, target, count)

    # Highlight confounded nodes
    if (length(confounder_names) > 0) {
      nodes$group <- ifelse(
        nodes$name %in% confounder_names,
        "confounded",
        "non_confounded"
      )
      nodes$color <- ifelse(
        nodes$name %in% confounder_names,
        "red",
        "lightblue"
      )
    } else {
      # give different colors to each name
      nodes$group <- nodes$name
      nodes$color <- scales::hue_pal()(length(unique(nodes$name)))
    }

    # Highlight confounded links
    links$group <- ifelse(
      nodes$name[links$source + 1] %in%
        confounder_names |
        nodes$name[links$target + 1] %in% confounder_names,
      "confounded",
      "non_confounded"
    )
    links$color <- ifelse(
      nodes$name[links$source + 1] %in%
        confounder_names |
        nodes$name[links$target + 1] %in% confounder_names,
      "red",
      "lightblue"
    )

    # Check if there are any confounders
    if (length(confounder_names) == 0) {
      # If no confounders, use grey color for all links
      links$color <- "#A9A9A9"
    }

    # Create the plotly Sankey plot
    sankey <- plot_ly(
      type = "sankey",
      orientation = "h",
      node = list(
        label = nodes$name,
        color = nodes$color,
        pad = 15,
        thickness = 20,
        line = list(color = "black", width = 0.5)
      ),
      link = list(
        source = links$source,
        target = links$target,
        value = links$count,
        color = links$color
      )
    ) %>%
      layout(
        title = "Sankey Diagram",
        font = list(size = 10)
      )
    return(sankey)
  }

  # Function to check for confounders without using design matrix
  check_confounders <- function(data, contrast_vars, covariate_column) {
    confounders <- list()
    confounder_names <- list()

    for (conf_var in contrast_vars) {
      for (level in levels(data[[covariate_column]])) {
        contingency_table <- table(
          data[[conf_var]],
          data[[covariate_column]] == level
        )
        for (row_name in rownames(contingency_table)) {
          if (
            contingency_table[row_name, 1] == 0 &
              contingency_table[row_name, 2] > 0
          ) {
            detailed_message <- paste(
              row_name,
              "is completely confounded with",
              covariate_column,
              level
            )
            confounders[[conf_var]] <- c(
              confounders[[conf_var]],
              detailed_message
            )
            confounder_names[[conf_var]] <- c(
              confounder_names[[conf_var]],
              row_name
            )
          }
        }
      }
    }

    return(list(confounders = confounders, confounder_names = confounder_names))
  }

  # Function to extract unique elements inside parentheses
  extract_unique_groups <- function(contrasts) {
    # Extract elements within parentheses
    elements_in_parentheses <- str_extract_all(contrasts, "(?<=\\().+?(?=\\))")

    # Extract elements without parentheses
    elements_no_parentheses <- str_extract_all(contrasts, "^[^-()]+-[^-()]+$")

    # Combine all extracted elements
    all_elements <- c(
      unlist(elements_in_parentheses),
      unlist(elements_no_parentheses)
    )

    # Split the elements by "-" and flatten the list
    split_elements <- unlist(lapply(all_elements, function(x) {
      unlist(strsplit(x, "-"))
    }))

    # Get unique elements
    unique_elements <- unique(split_elements)

    # Return the unique elements
    return(unique_elements)
  }

  # Function to draw an error message with text wrapping and respecting newlines
  draw_error_message <- function(message, color, width = 120) {
    # Split the message by newlines
    segments <- str_split(message, "\n")[[1]]
    # Wrap each segment
    wrapped_segments <- lapply(segments, str_wrap, width = width)
    # Combine the wrapped segments with newlines
    wrapped_message <- paste(unlist(wrapped_segments), collapse = "\n")

    grid.newpage()
    grid.rect(gp = gpar(fill = color, col = NA))
    grid.text(
      wrapped_message,
      x = 0.5,
      y = 0.5,
      gp = gpar(fontface = "italic", cex = 3, col = "black"),
      just = "center"
    )
  }

  # Function to get gene lists based on fold change and p-value thresholds
  getgenelists <- function(FClimit, pvallimit, pval) {
    upreggenes <- list()
    downreggenes <- list()
    for (i in 1:length(contrasts)) {
      if (pval == "pval") {
        results %>%
          dplyr::filter(
            .data[[colnames(FC)[i]]] > FClimit &
              .data[[colnames(pvals)[i]]] < pvallimit
          ) %>%
          pull(Gene) %>%
          length() -> upreggenes[[i]]
        results %>%
          dplyr::filter(
            .data[[colnames(FC)[i]]] < -FClimit &
              .data[[colnames(pvals)[i]]] < pvallimit
          ) %>%
          pull(Gene) %>%
          length() -> downreggenes[[i]]
      } else {
        results %>%
          dplyr::filter(
            .data[[colnames(FC)[i]]] > FClimit &
              .data[[colnames(pvaladj)[i]]] < pvallimit
          ) %>%
          pull(Gene) %>%
          length() -> upreggenes[[i]]
        results %>%
          dplyr::filter(
            .data[[colnames(FC)[i]]] < -FClimit &
              .data[[colnames(pvaladj)[i]]] < pvallimit
          ) %>%
          pull(Gene) %>%
          length() -> downreggenes[[i]]
      }
    }
    names(upreggenes) <- contrasts
    names(downreggenes) <- contrasts
    allreggenes <- rbind(unlist(upreggenes), unlist(downreggenes))
    rownames(allreggenes) <- c(
      paste0("upreg>", FClimit, ", ", pval, "<", pvallimit),
      paste0("downreg<-", FClimit, ", ", pval, "<", pvallimit)
    )
    return(allreggenes)
  }

  # Custom function to calculate standard error
  se <- function(x) {
    n <- length(x[!is.na(x)]) # Only non-NA values should be counted
    sd(x, na.rm = TRUE) / sqrt(n)
  }

  # Function to split data frame rows
  split_data_frame_rows <- function(data, chunk_size) {
    data <- as.data.frame(data)
    split(data, ceiling(seq_along(1:nrow(data)) / chunk_size))
  }

  # Function to split data frame columns
  split_data_frame_cols <- function(data, chunk_size) {
    col_indices <- split(
      seq_len(ncol(data)),
      ceiling(seq_along(seq_len(ncol(data))) / chunk_size)
    )
    lapply(col_indices, function(cols) data[, cols, drop = FALSE])
  }

  # Function to print chunks with titles
  print_chunk <- function(chunk, index, title) {
    for (i in seq_along(colnames(chunk))) {
      if (class(chunk[[i]]) == "list") {
        chunk[[i]] <- unlist(lapply(chunk[[i]], function(x) {
          paste(x, collapse = ",")
        }))
      }
    }

    # Create a grob with title
    title_grob <- textGrob(title, gp = gpar(fontsize = 15, fontface = "bold"))
    table_grob <- tableGrob(chunk, theme = ttheme_default(base_size = 10))

    # Measure the size of the table grob
    table_width <- convertWidth(grobWidth(table_grob), "npc", valueOnly = TRUE)
    table_height <- convertHeight(
      grobHeight(table_grob),
      "npc",
      valueOnly = TRUE
    )

    # Adjust font size based on table size relative to the graphical window
    scaling_factor <- min(1 / table_width, 1 / table_height)
    adjusted_font_size <- max(1 * scaling_factor, 6)
    # Ensure minimum font size is 6

    # Recreate the table grob with the adjusted font size
    table_grob <- tableGrob(
      chunk,
      theme = ttheme_default(base_size = adjusted_font_size)
    )

    # Arrange the title and table in a vertical layout
    combined_grob <- arrangeGrob(
      title_grob,
      table_grob,
      ncol = 1,
      heights = c(1, 10)
    )

    grid.newpage()
    pushViewport(viewport(width = 0.75, height = 0.75))
    # Adjust to fit within 75% of the window
    grid.draw(combined_grob)
    popViewport()
    return(NULL)
  }

  # Function to split and print a table
  split_and_print_table <- function(
    data,
    row_chunk_size,
    col_chunk_size,
    tabtitle
  ) {
    # Split rows into chunks
    row_chunks <- suppressWarnings(split_data_frame_rows(data, row_chunk_size))
    row_index <- 1
    for (row_chunk in row_chunks) {
      row_start <- row_index * row_chunk_size - row_chunk_size + 1
      row_stop <- row_start + nrow(row_chunk) - 1

      # Split columns into chunks for each row chunk
      col_chunks <- suppressWarnings(split_data_frame_cols(
        row_chunk,
        col_chunk_size
      ))
      col_index <- 1
      for (col_chunk in col_chunks) {
        col_start <- col_index * col_chunk_size - col_chunk_size + 1
        col_stop <- col_start + ncol(col_chunk) - 1
        title <- paste(
          tabtitle,
          "Subset Rows",
          row_start,
          "-",
          row_stop,
          "Cols",
          col_start,
          "-",
          col_stop
        )
        print_chunk(
          col_chunk,
          paste0(
            "Row ",
            row_start,
            "-",
            row_stop,
            ", Col ",
            col_start,
            "-",
            col_stop
          ),
          title
        )
        col_index <- col_index + 1
      }
      row_index <- row_index + 1
    }
  }

  # Function to wrap lines
  wraplines <- function(y) {
    j <- unlist(strsplit(y, "-"))
    k <- strwrap(j, width = 10)
    l <- paste(k, collapse = "\n-")
    return(l)
  }

  # Function to calculate group means
  calculate_group_means <- function(exp, design) {
    contgroups <- colnames(design)
    group_means <- list()
    group_se <- list()
    for (group_name in contgroups) {
      group_indices <- which(design[, group_name] == 1)
      group_means[[group_name]] <- rowMeans(exp[, group_indices, drop = FALSE])
      group_se[[group_name]] <- apply(exp[, group_indices, drop = FALSE], 1, se)
    }
    return(list(group_means, group_se))
  }

  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##

  # Limit sample metadata to only the samples we want to include
  sample_metadata.filt <- sample_metadata %>%
    filter(.data[[sample_names_column]] %in% samples_to_include)
  rownames(sample_metadata.filt) <- sample_metadata.filt[[sample_names_column]]

  # Make sure the counts matrix matches sample names
  counts_matrix <- counts_matrix %>%
    select(any_of(c(gene_names_column, samples_to_include)))
  counts_matrix <- counts_matrix %>%
    group_by(.data[[gene_names_column]]) %>%
    summarise_all(sum) %>%
    as.data.frame()
  counts_matrix.filt <- counts_matrix[, sample_metadata.filt[[
    sample_names_column
  ]]]
  rownames(counts_matrix.filt) <- counts_matrix[[gene_names_column]]

  # Ensure that there is only 1 sample per condition per subject,
  # and print out the groups
  if (length(donor_variable_column) > 0) {
    sample_metadata.filt.dedup <- sample_metadata.filt %>%
      group_by(across(all_of(c(
        donor_variable_column,
        contrast_variable_columns
      )))) %>%
      slice(1) %>%
      ungroup()
    counts_matrix.filt.dedup <- counts_matrix.filt[,
      sample_metadata.filt.dedup[[sample_names_column]]
    ]
    rownames(counts_matrix.filt.dedup) <- rownames(counts_matrix.filt)

    removed_samples <- sample_metadata.filt[[sample_names_column]][
      !sample_metadata.filt[[sample_names_column]] %in%
        sample_metadata.filt.dedup[[sample_names_column]]
    ]
    removedlength <- length(removed_samples)
    remainlength <- dim(sample_metadata.filt.dedup)[1]
    if (removedlength > 0) {
      removed_samples <- paste(removed_samples, collapse = ", ")
      cat(sprintf(
        paste0(
          "The following %s samples were removed due to multiple ",
          "samples per condition per subject: %s leaving %s samples."
        ),
        removedlength,
        removed_samples,
        remainlength
      ))
    } else {
      cat(sprintf(
        paste0(
          "No samples were removed due to multiple samples per ",
          "condition per subject."
        )
      ))
    }
  } else {
    stop("No donor variable column provided")
  }

  # Create another contrast column called contmerge for 1 and 2 contrast columns
  if (length(contrast_variable_columns) > 1) {
    sample_metadata.batch <- sample_metadata.filt.dedup %>%
      filter(if_all(
        all_of(contrast_variable_columns),
        ~ !is.na(.) & . != "NA"
      )) %>%
      dplyr::mutate(
        contmerge = paste0(
          .data[[contrast_variable_columns[1]]],
          ".",
          .data[[contrast_variable_columns[2]]]
        )
      ) %>%
      as.data.frame()
    sample_metadata.final <- sample_metadata.filt.dedup %>%
      dplyr::mutate(
        contmerge = paste0(
          .data[[contrast_variable_columns[1]]],
          ".",
          .data[[contrast_variable_columns[2]]]
        )
      ) %>%
      as.data.frame()
  } else {
    sample_metadata.batch <- sample_metadata.filt.dedup %>%
      filter(
        !is.na(!!sym(contrast_variable_columns)) &
          !!sym(contrast_variable_columns) != "NA"
      ) %>%
      dplyr::mutate(contmerge = .data[[contrast_variable_columns]]) %>%
      as.data.frame()
    sample_metadata.final <- sample_metadata.filt.dedup %>%
      dplyr::mutate(contmerge = .data[[contrast_variable_columns]]) %>%
      as.data.frame()
  }

  cat(sprintf(
    "\n\n Contmerge factor created from %s \n\n",
    paste(contrast_variable_columns, collapse = ", ")
  ))
  sample_metadata.final <- sample_metadata.final %>%
    mutate(across(all_of(covariate_columns), as.factor))
  rownames(sample_metadata.final) <- sample_metadata.final[[
    sample_names_column
  ]]

  # Print out individual groups
  sample_metadata.groups <- sample_metadata.final %>% filter(!is.na(contmerge))
  grouped_summary <- sample_metadata.groups %>%
    group_by(across(all_of(c(
      donor_variable_column,
      contrast_variable_columns
    )))) %>%
    dplyr::summarise(
      Sample_Names = paste(.data[[sample_names_column]], collapse = ", "),
      .groups = "drop"
    )
  wide_table <- grouped_summary %>%
    pivot_wider(
      names_from = contrast_variable_columns,
      values_from = Sample_Names
    )
  wide_table <- replace(wide_table, is.na(wide_table), "")
  wide_table <- as.data.frame(wide_table)

  counts_matrix.final <- counts_matrix.filt.dedup[, rownames(
    sample_metadata.groups
  )]

  # Make design matrix
  if (length(covariate_columns) > 0) {
    dm.formula <- as.formula(paste(
      "~0 + contmerge + ",
      paste(covariate_columns, sep = "+", collapse = "+")
    ))
  } else {
    dm.formula <- as.formula("~0 + contmerge")
  }
  cat("\n\nstatistical formula:\n")
  print(dm.formula)
  design <- model.matrix(dm.formula, sample_metadata.groups)

  for (i in 1:length(contrast_variable_columns)) {
    colnames(design) <- gsub(contrast_variable_columns[i], "", colnames(design))
  }
  colnames(design) <- gsub(":", ".", colnames(design))
  colnames(design) <- gsub("contmerge", "", colnames(design))

  # Identify non-estimable factors
  nonest0 <- colnames(design)[colSums(design) == 0]
  nonest1 <- colnames(design)[colSums(design) == 1]

  # Combine non-estimable factors
  nonest <- c(nonest0, nonest1)

  # Print out non-estimable factors
  if (length(nonest) > 0) {
    # Print initial dimensions
    print("Variables before removing non-estimable factors:")
    print(colnames(design))

    # Remove non-estimable factors from the design matrix
    design <- design[, !colnames(design) %in% nonest]

    # Remove uncategorized samples from sample metadata and counts matrix
    counts_matrix.final <- (counts_matrix.final[, rownames(design)]) %>%
      as.data.frame()
    sample_metadata.groups <- sample_metadata.groups %>%
      filter(.data[[sample_names_column]] %in% rownames(design))

    # Print dimensions after removing non-estimable factors
    print("Variables after removing non-estimable factors:")
    print(colnames(design))

    message(
      paste0(
        "The following non-estimable factors have been removed from ",
        "the design matrix as they do not provide sufficient ",
        "variability:\n"
      ),
      "Factors with all zero counts: ",
      paste(nonest0, collapse = ", "),
      "\n",
      "Factors with only one non-zero count: ",
      paste(nonest1, collapse = ", "),
      "\n"
    )
  } else {
    message(
      paste0(
        "All factors in the design matrix provide sufficient ",
        "variability. No non-estimable factors were found."
      )
    )
  }

  # Perform QR decomposition
  qr_check <- qr(design)

  # Extract the rank from the QR object
  qr_rank <- qr_check$rank

  # Get the number of columns in the design matrix
  num_columns <- ncol(design)

  # Check for collinearity
  if (qr_rank < num_columns) {
    cat("Collinearity detected: The design matrix does not have full rank.\n")
    confounded_covariates <- colnames(design)[qr_check$pivot[
      (qr_rank + 1):num_columns
    ]]
    print(confounded_covariates)
  } else {
    cat("No collinearity detected: The design matrix has full rank.\n")
    confounded_covariates <- NULL
  }

  if (!is.null(confounded_covariates)) {
    # Identify confounded variables - column where confounding
    # variables are found
    confounding_variable_column <- covariate_columns[grepl(
      covariate_columns,
      confounded_covariates
    )][1]
    confounder_info <- check_confounders(
      sample_metadata.groups,
      contrast_variable_columns,
      confounding_variable_column
    )
    confounders <- confounder_info$confounders
    confounder_names <- unique(unlist(confounder_info$confounder_names))

    # Draw Sankey Plot to show the experimental design
    sankeyplot <- plot_sankey_with_confounders(
      sample_metadata.groups,
      c(contrast_variable_columns),
      covariate_columns,
      confounder_names
    )
    message <- paste0("Error: The following covariates are confounded:\n")

    for (covariate in names(confounders)) {
      if (length(confounders[[covariate]]) > 0) {
        message <- paste0(
          message,
          "\n",
          covariate,
          "\n",
          paste(confounders[[covariate]], collapse = "\n "),
          "\n"
        )
      } else {
        message <- paste0(
          message,
          covariate,
          " (no specific confounding variable identified)\n"
        )
      }
    }

    message <- paste0(
      message,
      paste0(
        "\n\nThis indicates that your experimental design contains ",
        "confounding variables, which can bias your results.\n"
      ),
      paste0(
        "Consider revising your experimental design to ensure that ",
        "each covariate provides unique information.\n"
      ),
      paste0(
        "You may need to collect more data or adjust the grouping ",
        "of your samples."
      )
    )
    draw_error_message(message, "tomato3")
    print(sankeyplot)
    return(NULL)
  } else {
    sankeyplot <- plot_sankey(
      sample_metadata.groups,
      c(contrast_variable_columns),
      covariate_columns
    )
    message <- "No collinearity was detected in the design matrix."
    draw_error_message(message, "lightgreen")
    print(sankeyplot)

    # Splitting and printing the pivoted table for donor samples
    if (!is.null(wide_table)) {
      split_and_print_table(
        wide_table,
        row_number,
        column_number,
        "Pivot Table of Different Groups by Donor"
      )
    }

    # Make DGEList with all samples and perform normalization
    gene_names <- NULL
    gene_names$GeneID <- rownames(counts_matrix.final)

    if (input_in_log_counts == TRUE) {
      x <- DGEList(counts = 2^counts_matrix.final, genes = gene_names)
    } else {
      x <- DGEList(counts = counts_matrix.final, genes = gene_names)
    }

    if (
      normalization_method %in%
        c("TMM", "TMMwzp", "RLE", "upperquartile")
    ) {
      x <- edgeR::calcNormFactors(x, method = normalization_method)
      voom_normalize_method <- "none"
    } else {
      x <- edgeR::calcNormFactors(x, method = "TMM")
      voom_normalize_method <- normalization_method
    }

    # Run Voom (first round)
    vtemp <- voom(x, design, normalize.method = voom_normalize_method)
    corfit <- duplicateCorrelation(
      vtemp,
      design,
      block = sample_metadata.groups[[donor_variable_column]]
    )

    # Second round of voom incorporating duplicate correlation (as recommended)
    v <- voom(
      x,
      design,
      normalize.method = voom_normalize_method,
      block = sample_metadata.groups[[donor_variable_column]],
      correlation = corfit$consensus.correlation
    )
    corfit <- duplicateCorrelation(
      v,
      design,
      block = sample_metadata.groups[[donor_variable_column]]
    )
    fit <- lmFit(
      v,
      design,
      block = sample_metadata.groups[[donor_variable_column]],
      correlation = corfit$consensus.correlation
    )
    cm <- makeContrasts(contrasts = contrasts, levels = design)
    fit <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit)

    # Put together the results dataframe
    tstat <- fit2$t
    pvals <- fit2$p.value
    logFC <- fit2$coefficients
    FC <- 2^fit2$coefficients
    FC <- apply(FC, MARGIN = c(1, 2), function(x) ifelse(x < 1, -1 / x, x))
    SE <- sqrt(fit2$s2.post) * fit2$stdev.unscaled

    colnames(logFC) <- paste(colnames(logFC), "logFC", sep = "_")
    colnames(FC) <- paste(colnames(FC), "FC", sep = "_")
    colnames(SE) <- paste(colnames(SE), "SE", sep = "_")
    colnames(tstat) <- paste(colnames(tstat), "tstat", sep = "_")
    colnames(pvals) <- paste(colnames(pvals), "pval", sep = "_")
    pvaladj <- apply(pvals, 2, function(x) p.adjust(x, "BH"))
    colnames(pvaladj) <- paste(colnames(fit$coefficients), "adjpval", sep = "_")

    # Format group sample numbers
    sampsize <- colSums(design)
    titleval <- "Please note Sample size:"
    titletext <- paste(names(sampsize), sampsize, sep = "=", collapse = " \n ")
    sampletitle <- textGrob(
      paste(titleval, "\n", titletext, "\n\n\n"),
      gp = gpar(fontsize = 10)
    )
    contrast <- colnames(cm)
    connames <- strsplit(contrast, "-")
    connames <- lapply(connames, function(x) {
      gsub("\\(", "", gsub("\\)", "", x))
    })
    contrastsize <- lapply(connames, function(x) sampsize[unlist(x)])
    contrasttext <- paste(contrast, contrastsize, sep = " : ", collapse = "\n")
    contrasttext <- textGrob(
      paste("\n\n\nContrasts:\n", contrasttext),
      gp = gpar(fontsize = 10)
    )

    grid.newpage()
    grid.draw(sampletitle)
    grid.newpage()
    grid.draw(contrasttext)

    # Extract unique groups
    unique_groups <- extract_unique_groups(contrasts)
    rownames(sample_metadata.batch) <- sample_metadata.batch[[
      sample_names_column
    ]]

    if (return_batch_corrected_values == TRUE) {
      design.batch <- model.matrix(dm.formula, sample_metadata.batch)
      colnames(design.batch) <- gsub(":", ".", colnames(design.batch))
      colnames(design.batch) <- gsub("contmerge", "", colnames(design.batch))
      # Filter out batch columns with only one sample - no variance to adjust
      if (any(colSums(design.batch) == 1)) {
        # Identify columns in the design matrix with colSum == 1
        columns_to_remove <- which(colSums(design.batch) == 1)
        # Check if there are any columns to remove
        if (length(columns_to_remove) > 0) {
          # Get the names of the affected columns
          column_names <- colnames(design.batch)[columns_to_remove]
          # Find the corresponding sample(s)
          affected_samples <- apply(
            design.batch[, columns_to_remove, drop = FALSE],
            2,
            function(col) {
              rownames(design.batch)[col == 1]
            }
          )
          # Print out a message for each affected column
          for (i in seq_along(column_names)) {
            cat(
              "Removing batch variable:",
              column_names[i],
              "associated with sample:",
              affected_samples[i],
              "\n"
            )
          }
          # Remove these columns from the design matrix
          design.batch <- design.batch[, -columns_to_remove]
        }
      } else {
        cat(
          paste0(
            "No batch variables with only one sample detected in the ",
            "design matrix. All batch groups have sufficient samples ",
            "for proper batch effect removal.\n"
          )
        )
      }

      # Identify sparse batches (batches with only one sample)
      sparse_batches <- names(which(
        table(sample_metadata.batch[[donor_variable_column]]) == 1
      ))
      if (length(sparse_batches) > 0) {
        # Get the specific samples associated with these sparse batches
        sparse_samples <- rownames(sample_metadata.batch)[
          sample_metadata.batch[[donor_variable_column]] %in% sparse_batches
        ]
        # Construct the warning message
        warning_message <- paste0(
          paste0(
            "Warning: The following batches have only one sample and ",
            "are not fully estimable:\n"
          ),
          paste(sparse_batches, collapse = ", "),
          ".\n",
          "These batches correspond to the following samples:\n",
          paste(sparse_samples, collapse = ", "),
          ".\n",
          paste0(
            "Batch effects for these samples may not be accurately ",
            "removed and results are only estimates based on other ",
            "samples."
          )
        )

        # Display the warning message
        grid.newpage()
        grid.text(
          warning_message,
          x = 0.5,
          y = 0.5,
          gp = gpar(fontsize = 12, col = "red"),
          just = "center"
        )
      }
      if (!is.null(batch_column)) {
        design_sub <- design.batch[,
          !grepl(batch_column, colnames(design.batch))
        ] # Removes batch column from the design
        batch_corrected_values <- removeBatchEffect(
          v$E[, sample_metadata.batch[[sample_names_column]]],
          batch = sample_metadata.batch[[batch_column]],
          batch2 = sample_metadata.batch[[donor_variable_column]],
          design = design_sub
        )
      } else {
        batch_corrected_values <- removeBatchEffect(
          v$E[, sample_metadata.batch[[sample_names_column]]],
          batch = sample_metadata.batch[[donor_variable_column]],
          design = design.batch
        )
      }
      group_stats <- calculate_group_means(batch_corrected_values, design.batch)
      sample_exp <- batch_corrected_values
    } else {
      group_stats <- calculate_group_means(v$E, design)
      sample_exp <- v$E
    }

    group_means <- group_stats[[1]]
    group_means <- do.call(data.frame, group_means)
    colnames(group_means) <- paste(colnames(group_means), "Mean", sep = "_")

    group_se <- group_stats[[2]]
    group_se <- do.call(data.frame, group_se)
    colnames(group_se) <- paste(colnames(group_se), "SE", sep = "_")

    results <- cbind(
      group_means,
      group_se,
      FC,
      logFC,
      SE,
      tstat,
      pvals,
      pvaladj,
      sample_exp
    )
    results <- as.data.frame(results) %>% rownames_to_column("Gene")
    colnames(results) <- gsub(" - ", "-", colnames(results))

    # Format genelist sizes for typical cutoffs and
    # print number of DEGs per cutoff
    FCpval1 <- getgenelists(
      FClimit = fold_change_threshold,
      pvallimit = first_pvalue_threshold,
      "pval"
    )
    FCpval2 <- getgenelists(
      FClimit = fold_change_threshold,
      pvallimit = second_pvalue_threshold,
      "pval"
    )
    FCadjpval1 <- getgenelists(
      FClimit = fold_change_threshold,
      pvallimit = first_pvalue_threshold,
      "adjpval"
    )
    FCadjpval2 <- getgenelists(
      FClimit = fold_change_threshold,
      pvallimit = second_pvalue_threshold,
      "adjpval"
    )
    pvaltab <- rbind(FCpval1, FCpval2, FCadjpval1, FCadjpval2)
    colnames(pvaltab) <- sapply(colnames(pvaltab), function(x) wraplines(x))

    split_and_print_table(
      pvaltab,
      row_number,
      column_number,
      "Number of contrasts"
    )

    colnames(results) <- gsub("\\(", "", gsub("\\)", "", colnames(results)))

    if (!is.null(output_file_path) && nzchar(trimws(output_file_path))) {
      output_dir <- dirname(output_file_path)
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      write.table(
        results,
        file = output_file_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      message(paste0("Output written to: ", output_file_path))
    }

    return(results)
  }
}
