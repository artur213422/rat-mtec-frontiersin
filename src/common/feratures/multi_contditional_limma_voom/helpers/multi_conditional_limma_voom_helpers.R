
read_contrasts_file <- function(contrasts_file) {
  if (!file.exists(contrasts_file)) {
    stop(sprintf("Contrasts file does not exist: %s", contrasts_file))
  }

  contrasts_df <- read.csv(
    contrasts_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"contrast" %in% colnames(contrasts_df)) {
    stop("Contrasts file must contain a 'contrast' column.")
  }

  if (nrow(contrasts_df) == 0L) {
    stop("Contrasts file is empty.")
  }

  if (any(is.na(contrasts_df$contrast)) || any(trimws(contrasts_df$contrast) == "")) {
    stop("Contrasts file contains missing or empty contrast names.")
  }

  if (anyDuplicated(contrasts_df$contrast)) {
    stop("Contrast names must be unique.")
  }

  coefficient_columns <- setdiff(colnames(contrasts_df), "contrast")

  if (length(coefficient_columns) < 2L) {
    stop("Contrasts file must contain at least two condition columns.")
  }

  for (column in coefficient_columns) {
    if (!is.numeric(contrasts_df[[column]])) {
      stop(sprintf("Contrast coefficient column '%s' must be numeric.", column))
    }

    if (any(is.na(contrasts_df[[column]])) || any(!is.finite(contrasts_df[[column]]))) {
      stop(sprintf("Contrast coefficient column '%s' contains invalid values.", column))
    }
  }

  coefficient_matrix <- as.matrix(
    contrasts_df[
      ,
      coefficient_columns,
      drop = FALSE
    ]
  )

  if (any(rowSums(abs(coefficient_matrix)) == 0)) {
    stop("Every contrast must contain at least one non-zero coefficient.")
  }

  # Comparisons between condition means should sum to zero.
  if (any(abs(rowSums(coefficient_matrix)) > 1e-8)) {
    stop(paste0("Every contrast must sum to zero. This is required for comparisons between condition means."))
  }

  return(contrasts_df)
}


prepare_matrix_and_design <- function(counts_df, metadata_df) {
  counts_matrix <- xtabs(
    count ~ gene_symbol + sample_id,
    data = counts_df
  )

  counts_matrix <- matrix(
    as.numeric(counts_matrix),
    nrow = nrow(counts_matrix),
    ncol = ncol(counts_matrix),
    dimnames = dimnames(counts_matrix)
  )

  metadata_df <- metadata_df[
    match(
      colnames(counts_matrix),
      metadata_df$sample_id
    ),
    ,
    drop = FALSE
  ]

  stopifnot(
    identical(
      colnames(counts_matrix),
      metadata_df$sample_id
    )
  )

  metadata_df$condition <- factor(metadata_df$condition)

  if (nlevels(metadata_df$condition) < 2L) {
    stop("At least two conditions are required.")
  }

  # One model coefficient per condition.
  design <- model.matrix(
    ~ 0 + condition,
    data = metadata_df
  )

  colnames(design) <- levels(metadata_df$condition)

  condition_columns <- colnames(design)
  metadata_df$experiment <- factor(metadata_df$experiment)
  if (nlevels(metadata_df$experiment) > 1L) {
    # Reference-coded experiment effects; the condition columns already span
    # the intercept. Do not add another intercept or all experiment indicators.
    experiment_design <- model.matrix(~ experiment, data = metadata_df)[, -1, drop = FALSE]
    if (any(colnames(experiment_design) %in% condition_columns)) {
      stop("Experiment coefficient names collide with condition names.")
    }
    design <- cbind(design, experiment_design)
  }

  rownames(design) <- metadata_df$sample_id

  if (qr(design)$rank < ncol(design)) {
    stop("The design matrix is not full rank.")
  }

  residual_df <- nrow(design) - qr(design)$rank

  if (residual_df <= 0L) {
    stop("The design has no residual degrees of freedom.")
  }

  dge <- edgeR::DGEList(counts = counts_matrix)

  keep <- edgeR::filterByExpr(dge, design = design)

  if (!any(keep)) {
    stop("No genes passed low-expression filtering.")
  }

  message(sprintf("Expression filtering: keeping %d of %d genes.", sum(keep), length(keep)))

  dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- edgeR::calcNormFactors(dge)

  list(
    design = design,
    condition_columns = condition_columns,
    dge = dge,
    metadata = metadata_df
  )
}


prepare_contrast_matrix <- function(
  contrasts_df,
  design,
  condition_columns = colnames(design)
) {
  contrast_names <- contrasts_df$contrast

  coefficient_df <- contrasts_df[
    ,
    setdiff(colnames(contrasts_df), "contrast"), drop = FALSE
  ]

  missing_design_columns <- setdiff(condition_columns, colnames(coefficient_df))

  extra_contrast_columns <- setdiff(colnames(coefficient_df), condition_columns)

  if (length(missing_design_columns) > 0L) {
    stop(sprintf(
      "Contrasts are missing design conditions: %s",
      paste(
        missing_design_columns,
        collapse = ", "
      )
    ))
  }

  if (length(extra_contrast_columns) > 0L) {
    stop(sprintf("Contrasts contain unknown conditions: %s", paste(extra_contrast_columns, collapse = ", ")))
  }

  # Experiment effects adjust the fit but cancel in condition comparisons.
  for (column in setdiff(colnames(design), condition_columns)) {
    coefficient_df[[column]] <- 0
  }

  # Match contrast coefficients to design-column order.
  coefficient_df <- coefficient_df[
    ,
    colnames(design),
    drop = FALSE
  ]

  contrast_matrix <- t(as.matrix(coefficient_df))

  storage.mode(contrast_matrix) <- "double"

  rownames(contrast_matrix) <- colnames(design)
  colnames(contrast_matrix) <- contrast_names

  contrast_matrix
}
