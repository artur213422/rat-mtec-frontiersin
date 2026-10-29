

source(here::here("src", "common", "feratures", "R", "json_config_handler.R"))
source(here::here("src", "common", "feratures", "multi_contditional_limma_voom", "helpers", "read_inputs.R"))
source(here::here("src", "common", "feratures", "multi_contditional_limma_voom", "helpers", "multi_conditional_limma_voom_helpers.R"))

###############################################################################
# PREPARE DATA
################################################################################

config <- read_json_config()

counts_file <- get_required_string(config, "counts_file")
metadata_file <- get_required_string(config, "metadata_file")
contrasts_file <- get_required_string(config, "contrasts_file")
output_file <- get_required_string(config, "output_file")
induction_log2_fc <- config$induction_log2_fc
specific_log2_fc <- config$specific_log2_fc

counts_df <- read_counts_file(counts_file)
metadata_df <- read_metadata_file(metadata_file)
contrasts_df  <- read_contrasts_file(contrasts_file)


validate_counts_and_metadata_df(counts_df, metadata_df)

###############################################################################
# PREPARE COUNTS AND DESIGN
################################################################################
matrix_and_design <- prepare_matrix_and_design(counts_df, metadata_df)

design <- matrix_and_design$design
dge <- matrix_and_design$dge

contrast_matrix <- prepare_contrast_matrix(
  contrasts_df = contrasts_df,
  design = design,
  condition_columns = matrix_and_design$condition_columns
)

###############################################################################
# RUN VOOM
################################################################################

voom_result <- limma::voom(dge, design, plot = FALSE)


###############################################################################
# RUN ANALYSIS & EXTRACT CONTRASTS
###############################################################################
result_list <- vector(mode = "list", length = ncol(contrast_matrix))

contrast_names <- colnames(contrast_matrix)

for (contrast_index in seq_along(contrast_names)) {
  contrast_name <- contrast_names[[contrast_index]]

  comparison_design <- limma::contrastAsCoef(
    design,
    contrast_matrix[, contrast_index, drop = FALSE]
  )

  fit <- limma::lmFit(
    voom_result,
    comparison_design$design
  )

  # Use a larger threshold for IFN-versus-control induction.
  # All other contrasts are target-IFN versus another IFN.
  if (grepl("_vs_ctrl$", contrast_name)) {
    treat_lfc <- induction_log2_fc
  } else {
    treat_lfc <- specific_log2_fc

  }

  contrast_fit <- limma::treat(
    fit,
    lfc = treat_lfc
  )

  contrast_results <- limma::topTable(
    contrast_fit,
    coef = contrast_name,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )

  stat_columns <- c(
    "logFC",
    "t",
    "P.Value",
    "adj.P.Val"
  )

  if (
    nrow(contrast_results) == 0L ||
    any(!is.finite(
      as.matrix(
        contrast_results[
          ,
          stat_columns,
          drop = FALSE
        ]
      )
    ))
  ) {
    stop(sprintf("Limma returned empty or non-finite results for contrast '%s'.", contrast_name))
  }

  result_list[[contrast_index]] <- data.frame(
    gene_symbol = rownames(contrast_results),
    contrast = contrast_name,
    log2_fold_change = contrast_results$logFC,
    p_val = contrast_results$P.Value,
    t_statistic = contrast_results$t,
    average_expression = contrast_results$AveExpr,
    p_adj = contrast_results$adj.P.Val,
    test = "limma_voom_contrast",
    row.names = NULL,
    check.names = FALSE
  )
}

results <- do.call(rbind, result_list)

###############################################################################
# SAVE
###############################################################################


write.csv(results, output_file, row.names = FALSE)

message("Conditional limma voom is done! Saved result at: ", output_file)
