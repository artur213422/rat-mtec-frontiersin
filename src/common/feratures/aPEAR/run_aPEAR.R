source(here::here("src", "common", "feratures", "R", "json_config_handler.R"))
source(here::here("src", "common", "feratures", "aPEAR", "helpers", "aPEAR_helpers.R"))

###############################################################################
# READ CONFIG
###############################################################################

config <- read_json_config()

gsea_file <- get_required_string(config, "gsea_file")
output_svg <- get_required_string(config, "output_svg")
output_clusters <- get_required_string(config, "output_clusters")

gsea_df <- read_gsea(gsea_file)

###############################################################################
# SELECT PATHWAYS
###############################################################################

selected_df <- gsea_df[
  !is.na(gsea_df$nes) &
    is.finite(gsea_df$nes) &
    !is.na(gsea_df$fdr_q_value) &
    is.finite(gsea_df$fdr_q_value) &
    gsea_df$fdr_q_value <= 0.05 &
    !is.na(gsea_df$leading_edge_genes) &
    nzchar(trimws(gsea_df$leading_edge_genes)),
  ,
  drop = FALSE
]

# aPEAR expects genes within a pathway to be separated by "/".
selected_df$pathwayGenes <- gsub(
  pattern = ";",
  replacement = "/",
  x = selected_df$leading_edge_genes,
  fixed = TRUE
)

selected_df$leadingEdgeSize <- vapply(
  selected_df$pathwayGenes,
  count_pathway_genes,
  FUN.VALUE = integer(1)
)

selected_df <- selected_df[
  selected_df$leadingEdgeSize > 0L,
  ,
  drop = FALSE
]

enrichment_df <- data.frame(
  Description = selected_df$term,
  pathwayGenes = selected_df$pathwayGenes,
  NES = selected_df$nes,
  leadingEdgeSize = selected_df$leadingEdgeSize,
  pvalue = selected_df$nominal_p_value,
  p.adjust = selected_df$fdr_q_value,
  stringsAsFactors = FALSE
)


###############################################################################
# HANDLE TOO FEW PATHWAYS
###############################################################################

if (nrow(enrichment_df) < 2L) {
  handle_too_few_pathways()
}

###############################################################################
# RUN aPEAR
###############################################################################

set.seed(123)

cluster_result <- tryCatch(
  aPEAR::findPathClusters(
    enrichment = enrichment_df,
    similarity = "jaccard",
    cluster = "markov",
    clusterName = "pagerank",
    minClusterSize = 2,
    verbose = TRUE
  ),
  error = function(error) {
    if (identical(conditionMessage(error), "No clusters found.")) {
      return(NULL)
    }

    stop(error)
  }
)

if (is.null(cluster_result)) {
  cluster_output <- data.frame(
    term = enrichment_df$Description,
    nes = enrichment_df$NES,
    nominal_p_value = enrichment_df$pvalue,
    fdr_q_value = enrichment_df$p.adjust,
    leading_edge_gene_count = enrichment_df$leadingEdgeSize,
    cluster = NA_character_,
    stringsAsFactors = FALSE
  )

  write.csv(
    cluster_output,
    output_clusters,
    row.names = FALSE,
    quote = TRUE,
    na = "",
    fileEncoding = "UTF-8"
  )

  message("Skipping aPEAR: no clusters met the minimum cluster size.")
  quit(save = "no", status = 0L)
}

clusters <- as.data.frame(
  cluster_result$clusters,
  stringsAsFactors = FALSE
)

###############################################################################
# CREATE NETWORK FIGURE
###############################################################################

if (nrow(clusters) == 0L) {
  message("No clusters met aPEAR's minimum cluster size.")
} else {
  set.seed(123)

    # ggrepel fails for some tiny or single-cluster networks.
  repel_labels <- (
    nrow(clusters) > 2L &&
      length(unique(clusters$Cluster)) > 1L
  )

  network_plot <- aPEAR::plotPathClusters(
    enrichment = enrichment_df,
    sim = cluster_result$similarity,
    clusters = cluster_result$clusters,
    colorBy = "NES",
    colorType = "nes",
    nodeSize = "leadingEdgeSize",
    innerCutoff = 0.1,
    outerCutoff = 0.5,
    drawEllipses = FALSE,
    repelLabels = repel_labels,
    fontSize = 3,
    verbose = TRUE
  )

  svglite::svglite(
    output_svg,
    width = 16,
    height = 10
  )

  print(network_plot)

  grDevices::dev.off()
}