

handle_too_few_pathways <- function (){
   cluster_output <- data.frame(
    term = enrichment_df$Description,
    cluster = rep(
      NA_character_,
      nrow(enrichment_df)
    ),
    nes = enrichment_df$NES,
    nominal_p_value = enrichment_df$pvalue,
    fdr_q_value = enrichment_df$p.adjust,
    leading_edge_gene_count = enrichment_df$leadingEdgeSize,
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

  message("Skipping aPEAR: only ", nrow(enrichment_df), " pathway(s) passed filtering.")

  quit(save = "no", status = 0L)
}

read_gsea <- function (gsea_file) {
  gsea_df <- read.csv(
    gsea_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8"
  )

  required_columns <- c(
    "term",
    "nes",
    "nominal_p_value",
    "fdr_q_value",
    "leading_edge_genes"
  )

  missing_columns <- setdiff(required_columns, names(gsea_df))

  if (length(missing_columns) > 0L) {
    stop("GSEA input is missing required column(s): ", paste(missing_columns, collapse = ", "))
  }

  if ("dataset_id" %in% names(gsea_df) && length(unique(gsea_df$dataset_id)) > 1L) {
    stop("aPEAR must be run separately for each dataset. The input contains multiple dataset_id values.")
  }

  if (anyNA(gsea_df$term) || any(!nzchar(trimws(gsea_df$term)))) {
    stop("The GSEA input contains missing or empty pathway names.")
  }

  if (anyDuplicated(gsea_df$term)) {
    duplicate_terms <- unique(gsea_df$term[duplicated(gsea_df$term)])

    stop("The GSEA input contains duplicate pathway names: ", paste(head(duplicate_terms, 20L), collapse = ", "))
  }

  return(gsea_df)
}


count_pathway_genes <- function(pathway_genes) {
  if (
    is.na(pathway_genes) ||
    !nzchar(trimws(pathway_genes))
  ) {
    return(0L)
  }

  genes <- strsplit(
    pathway_genes,
    split = "/",
    fixed = TRUE
  )[[1]]

  genes <- unique(
    trimws(genes[nzchar(trimws(genes))])
  )

  length(genes)
}
