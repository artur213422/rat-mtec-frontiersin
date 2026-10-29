args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript gsea.R '<config_json>'")
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Required R package is not installed: jsonlite")
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop("Required R package is not installed: fgsea")
}


get_required <- function(config, name) {
  if (!(name %in% names(config))) {
    stop(sprintf("Missing required config field: %s", name))
  }

  config[[name]]
}


get_integer <- function(config, name, minimum) {
  value <- get_required(config, name)

  if (length(value) != 1 || is.na(value) || value < minimum || value != as.integer(value)) {
    stop(sprintf("Config field '%s' must be an integer >= %d", name,minimum))
  }

  as.integer(value)
}


config <- jsonlite::fromJSON(args[[1]])

min_size <- get_integer(
  config,
  "min_size",
  minimum = 1
)

max_size <- get_integer(
  config,
  "max_size",
  minimum = min_size
)

permutation_num <- get_integer(
  config,
  "permutation_num",
  minimum = 1
)

rnk_file <- get_required(config, "rnk_file")
gene_sets_file <- get_required(config, "gene_sets_file")
result_file <- get_required(config, "result_file")
mode <- get_required(config, "mode")


# ----------------------------------------------------------------------
# Read ranked statistics
# ----------------------------------------------------------------------

rank_df <- read.csv(
  rnk_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

if (!identical(names(rank_df), c("gene", "score"))) {
  stop("Rank file must have exactly the columns: gene, score")
}

if (anyNA(rank_df$gene) || any(rank_df$gene == "") || anyDuplicated(rank_df$gene)) {
  stop("Rank file contains missing, empty, or duplicate genes")
}

if (anyNA(rank_df$score) || any(!is.finite(rank_df$score))) {
  stop("Rank file contains missing or non-finite scores")
}

ranked_statistics <- rank_df$score
names(ranked_statistics) <- rank_df$gene

ranked_statistics <- sort(
  ranked_statistics,
  decreasing = TRUE
)


# ----------------------------------------------------------------------
# Read GMT gene sets
# ----------------------------------------------------------------------

pathways <- fgsea::gmtPathways(gene_sets_file)

if (length(pathways) == 0) {
  stop("No gene sets were read from the GMT file")
}


# ----------------------------------------------------------------------
# Run negative one-tailed GSEA
# ----------------------------------------------------------------------

set.seed(123)

fgsea_result <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = ranked_statistics,
  minSize = min_size,
  maxSize = max_size,
  gseaParam = 1,
  scoreType = mode,
  nPermSimple = permutation_num,
  eps = 0
)


# ----------------------------------------------------------------------
# Create a Python-compatible result
# ----------------------------------------------------------------------

output_df <- data.frame(
  term = fgsea_result$pathway,
  es = fgsea_result$ES,
  nes = fgsea_result$NES,
  nominal_p_value = fgsea_result$pval,
  fdr_q_value = fgsea_result$padj,
  leading_edge_genes = vapply(
    fgsea_result$leadingEdge,
    paste,
    collapse = ";",
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  output_df,
  file = result_file,
  row.names = FALSE,
  quote = TRUE,
  na = "",
  fileEncoding = "UTF-8"
)

message("Done: ", result_file)