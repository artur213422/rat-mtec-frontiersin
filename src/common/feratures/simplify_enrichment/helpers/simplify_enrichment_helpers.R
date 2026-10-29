source(here::here("src", "common", "feratures", "R", "json_config_handler.R"))


get_required_ontology <- function(config, name) {
  value <- get_required(config, name)

  valid_ontologies <- c("BP", "CC", "MF")

  if (!is.character(value) || length(value) != 1L || is.na(value) || !value %in% valid_ontologies) {
    stop(sprintf(
      "Config field '%s' must be one of: %s. Got: %s",
      name,
      paste(valid_ontologies, collapse = ", "),
      paste(value, collapse = ", ")
    ))
  }

  return(value)
}


get_database <- function(reference_species) {
    if (reference_species == "mouse") {
        return("org.Mm.eg.db")
    } else if (reference_species == "rat") {
        return("org.Rn.eg.db")
    } else {
        stop(
            "reference_species must be one of: 'mouse', 'rat'"
        )
    }
}


create_go_lists <- function(gsea_df) {
    required_columns <- c(
        "dataset_id",
        "term",
        "fdr_q_value"
    )

    missing_columns <- setdiff(required_columns, names(gsea_df))

    if (length(missing_columns) > 0L) {
        stop(sprintf("GSEA data is missing required columns: %s", paste(missing_columns, collapse = ", ")))
    }

    go_lists <- lapply(
        split(gsea_df, gsea_df$dataset_id),
        function(dataset_df) {
            go_ids <- trimws(
                sub(
                    "\\|.*$",
                    "",
                    dataset_df$term
                )
            )

            invalid_ids <- !grepl("^GO:[0-9]{7}$", go_ids)

            if (any(invalid_ids)) {
                stop(sprintf("Invalid GO terms found: %s", paste(unique(dataset_df$term[invalid_ids]), collapse = ", ")))
            }

            setNames(dataset_df$fdr_q_value, go_ids
            )
        }
    )

    return(go_lists)
}