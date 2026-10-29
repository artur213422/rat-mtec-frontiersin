
read_counts_file <- function(counts_file){
  if (!file.exists(counts_file)) {
    stop(sprintf("Counts file does not exist: %s", counts_file))
  }

  df <- read.csv(
    counts_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_columns <- c(
    "gene_symbol",
    "sample_id",
    "count"
  )

  missing_columns <- setdiff(required_columns, names(df))

  if (length(missing_columns) > 0L) {
    stop(sprintf("Counts file is missing required columns: %s", paste(missing_columns, collapse = ", ")))
  }

  if (nrow(df) == 0L) {
    stop("Counts file is empty.")
  }

  if (any(is.na(df$gene_symbol)) || any(df$gene_symbol == "")) {
    stop("Counts file contains missing or empty gene_symbol values.")
  }

  if (any(is.na(df$sample_id)) || any(df$sample_id == "")) {
    stop("Counts file contains missing or empty sample_id values.")
  }

  if (!is.numeric(df$count)) {
    stop("Counts column 'count' must be numeric.")
  }

  if (any(is.na(df$count))) {
    stop("Counts column contains NA values.")
  }

  if (any(!is.finite(df$count))) {
    stop("Counts column contains non-finite values.")
  }

  if (any(df$count < 0)) {
    stop("Counts cannot be negative.")
  }

  duplicated_rows <- duplicated(
    df[c("gene_symbol", "sample_id")]
  )

  if (any(duplicated_rows)) {
    stop("Counts file contains duplicate gene_symbol + sample_id combinations.")
  }

  return(df)
}

read_metadata_file <- function(metadata_file) {
  if (!file.exists(metadata_file)) {
    stop(sprintf("Metadata file does not exist: %s", metadata_file))
  }

  df <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_columns <- c(
    "sample_id",
    "experiment",
    "condition"
  )

  missing_columns <- setdiff(required_columns, names(df))

  if (length(missing_columns) > 0L) {
    stop(sprintf("Metadata file is missing required columns: %s",  paste(missing_columns, collapse = ", ")))
  }

  if (nrow(df) == 0L) {
    stop("Metadata file is empty.")
  }

  for (column in required_columns) {
    if (any(is.na(df[[column]]))) {
      stop(sprintf("Metadata column '%s' contains NA values.", column))
    }

    if (any(trimws(as.character(df[[column]])) == "")) {
      stop(sprintf("Metadata column '%s' contains empty values.", column))
    }
  }

  if (anyDuplicated(df$sample_id)) {
    stop("Metadata contains duplicate sample_id values.")
  }

  return(df)
}

validate_counts_and_metadata_df <- function(counts_df, metadata_df){
  # Only 2 conditions allowed
  conditions <- unique(metadata_df$condition)

  if (length(conditions) < 2L) {
    stop(sprintf("Metadata must contain at least 2 conditions. Found: %s", paste(conditions, collapse = ", ")))
  }

  # sample_id values must match exactly between counts and metadata
  counts_ids <- unique(counts_df$sample_id)
  metadata_ids <- unique(metadata_df$sample_id)

  if (!setequal(counts_ids, metadata_ids)) {
    only_in_counts <- setdiff(counts_ids, metadata_ids)
    only_in_metadata <- setdiff(metadata_ids, counts_ids)

    stop(
      sprintf(
        paste0("sample_id mismatch between counts and metadata. Only in counts: [%s]. Only in metadata: [%s]."),
        paste(only_in_counts, collapse = ", "),
        paste(only_in_metadata, collapse = ", ")
      )
    )
  }

  # This is not allowed:
  #gene_symbol  sample_id  count
  #GENE1        ctrl1      20
  #GENE1        ko1        50
  #GENE2        ctrl1      30
  # gene_symbols must match accross samples
  # Every gene must have one measurement for every sample.
  samples_per_gene <- table(counts_df$gene_symbol)
  expected_samples <- nrow(metadata_df)

  incomplete_genes <- names(samples_per_gene)[
    samples_per_gene != expected_samples
  ]

  if (length(incomplete_genes) > 0L) {
    stop(sprintf(
      paste0("Counts are missing gene/sample combinations for %d genes. Each gene must have one row per sample, including explicit zeros. Examples: %s"),
      length(incomplete_genes),
      paste(head(incomplete_genes, 10L), collapse = ", ")
    ))
  }
}


validate_reference_and_target_condition <- function(reference_condition, target_condition, metadata_df){
 if (reference_condition == target_condition) {
    stop("reference_condition and target_condition must be different.")
  }

  conditions <- unique(metadata_df$condition)

  if (!(reference_condition %in% conditions)) {
    stop(sprintf("reference_condition '%s' is not present in metadata.", reference_condition))
  }

  if (!(target_condition %in% conditions)) {
    stop(sprintf("target_condition '%s' is not present in metadata.", target_condition))
  }

  if (!setequal(conditions, c(reference_condition, target_condition))) {
    stop(sprintf("Metadata conditions [%s] do not match reference/target conditions [%s, %s].", paste(conditions, collapse = ", "), reference_condition, target_condition))
  }

}