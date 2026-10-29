source(here::here("src", "common", "feratures", "R", "json_config_handler.R"))
source(here::here("src", "common", "feratures", "limma_voom", "helpers", "read_inputs.R"))
source(here::here("src", "common", "feratures", "limma_voom", "helpers", "voom_helpers.R"))

###############################################################################
# PREPARE DATA
################################################################################

config <- read_json_config()

counts_file <- get_required_string(config, "counts_file")
metadata_file <- get_required_string(config, "metadata_file")

counts_df <- read_counts_file(counts_file)
metadata_df <- read_metadata_file(metadata_file)

output_file <- get_required_string(config, "output_file")

reference_condition <- get_required_string(config, "reference_condition")
target_condition <- get_required_string(config, "target_condition")

validate_counts_and_metadata_df(counts_df, metadata_df)
validate_reference_and_target_condition(reference_condition, target_condition, metadata_df)

###############################################################################
# RUN VOOM
################################################################################
prepared <- prepare_matrix_and_design(counts_df, metadata_df)

voom_result <- limma::voom(
  prepared$dge,
  prepared$design,
  plot = FALSE
)

###############################################################################
# RUN LIMMA
################################################################################

fit <- limma::lmFit(voom_result, prepared$design)

fit <- limma::eBayes(fit)

condition_coef <- grep(
  "^condition",
  colnames(prepared$design)
)

results <- limma::topTable(
  fit,
  coef = condition_coef,
  number = Inf,
  sort.by = "none"
)

stat_columns <- c("logFC", "t", "P.Value", "adj.P.Val")

if (nrow(results) == 0L || any(!is.finite(as.matrix(results[, stat_columns, drop = FALSE])))) {
  stop("Limma returned empty results or non-finite statistics.")
}

###############################################################################
# SAVE
###############################################################################

results <- data.frame(
  gene_symbol = rownames(results),
  log2_fold_change = results$logFC,
  target_condition = target_condition,
  reference_condition = reference_condition,
  p_val = results$P.Value,
  t_statistic = results$t,
  average_expression = results$AveExpr,
  p_adj = results$adj.P.Val,
  test = "limma_voom",
  row.names = NULL
)

write.csv(results, output_file, row.names = FALSE)

message("Limma voom is done! Saved result at: ", output_file)