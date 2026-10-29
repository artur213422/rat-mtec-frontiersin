

prepare_matrix_and_design <- function(counts_df, metadata_df){
  counts_matrix <- xtabs(
    count ~ gene_symbol + sample_id,
    data = counts_df
  )

  # Convert the contingency table to a plain numeric matrix for voom.
  counts_matrix <- matrix(
    as.numeric(counts_matrix),
    nrow = nrow(counts_matrix),
    ncol = ncol(counts_matrix),
    dimnames = dimnames(counts_matrix)
  )

  # Make metadata order match count-matrix columns
  metadata_df <- metadata_df[
    match(colnames(counts_matrix), metadata_df$sample_id),
  ]

  stopifnot(
    identical(
      colnames(counts_matrix),
      metadata_df$sample_id
    )
  )

  # Explicitly define comparison direction
  metadata_df$condition <- factor(
    metadata_df$condition,
    levels = c(
      reference_condition,
      target_condition
    )
  )

  metadata_df$condition <- factor(metadata_df$condition)
  metadata_df$experiment <- factor(metadata_df$experiment)


  # Build statistical design
  if (nlevels(metadata_df$experiment) > 1L) {
    design <- model.matrix(
      ~ experiment + condition,
      data = metadata_df
    )
  } else {
    design <- model.matrix(
      ~ condition,
      data = metadata_df
    )
  }

  design_rank <- qr(design)$rank

  # Can the model distinguish all requested effects?
  if (design_rank < ncol(design)) {
    stop(paste0("Invalid design: model effects cannot be separated. Check whether experiment and condition are confounded."))
  }

  # Is there information left to estimate sample variability?
  residual_df <- nrow(design) - design_rank

  if (residual_df <= 0L) {
    stop(paste0("Invalid design: no residual degrees of freedom. More independent samples or a simpler justified model are required."))
  }

  dge <- edgeR::DGEList(
    counts = counts_matrix
  )

  # Keep genes with sufficient counts in enough samples.
  keep <- edgeR::filterByExpr(dge, design = design)

  if (!any(keep)) {
    stop("No genes passed low-expression filtering.")
  }

  message(sprintf(
    "Expression filtering: keeping %d of %d genes.",
    sum(keep),
    length(keep)
  ))

  # Recalculate library sizes after removing genes.
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  dge <- edgeR::calcNormFactors(dge)
    return(
      list(
        design = design,
        dge = dge
      )
    )
}