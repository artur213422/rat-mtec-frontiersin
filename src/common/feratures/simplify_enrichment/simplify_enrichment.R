source(here::here("src", "common", "feratures", "R", "json_config_handler.R"))
source(here::here("src", "common", "feratures", "simplify_enrichment", "helpers", "simplify_enrichment_helpers.R"))

###############################################################################
# PREPARE DATA
################################################################################

config <- read_json_config()

padj_cutoff <- get_required_double(config, "padj_cutoff")
reference_species <- get_required_string(config, "reference_species")
ontology <- get_required_ontology(config, "ontology")
gsea_file <- get_required_string(config, "gsea_file")
output_svg <- get_required_string(config, "output_svg")
output_clusters <- get_required_string(config, "output_clusters")
analysis_mode <- get_required_string(config, "analysis_mode")

db <- get_database(reference_species)

gsea_df <- read.csv(gsea_file)

go_lists <- create_go_lists(gsea_df)

###############################################################################
# RUN ANALYSE
################################################################################

svglite::svglite(
    output_svg,
    width = 16,
    height = 10
)

if (analysis_mode == "single") {

    significant_df <- gsea_df[
        !is.na(gsea_df$fdr_q_value) &
        gsea_df$fdr_q_value < padj_cutoff,
    ]

    go_ids <- unique(
        sub(
            "^\\s*(GO:[0-9]{7}).*$",
            "\\1",
            significant_df$term
        )
    )

    go_ids <- go_ids[
        grepl("^GO:[0-9]{7}$", go_ids)
    ]

    if (length(go_ids) < 2) {

        message(sprintf(
            "Skipping clustering: only %d significant GO term(s).",
            length(go_ids)
        ))

        plot.new()
        text(
            0.5,
            0.5,
            sprintf(
                "Not enough significant GO terms for clustering (n = %d)",
                length(go_ids)
            )
        )

        clusters <- data.frame(
            id = as.character(go_ids),
            cluster = rep(NA_integer_, length(go_ids)),
            stringsAsFactors = FALSE
        )

    } else {

        similarity_matrix <- simplifyEnrichment::GO_similarity(
            go_id = go_ids,
            ont = ontology,
            db = db
        )

        remaining_go_ids <- rownames(similarity_matrix)

        if (nrow(similarity_matrix) < 2) {

            message(sprintf(
                "Skipping clustering: only %d valid GO term(s) remained.",
                length(remaining_go_ids)
            ))

            plot.new()
            text(
                0.5,
                0.5,
                sprintf(
                    "Not enough valid GO terms for clustering (n = %d)",
                    length(remaining_go_ids)
                )
            )

            clusters <- data.frame(
                id = as.character(remaining_go_ids),
                cluster = rep(
                    NA_integer_,
                    length(remaining_go_ids)
                ),
                stringsAsFactors = FALSE
            )

        } else {

            clusters <- simplifyEnrichment::simplifyGO(
                similarity_matrix,
                min_term = 2,
                use_raster = FALSE
            )
        }
    }

} else if (analysis_mode == "comparison") {

    significant_go_ids <- as.character(
        unique(unlist(
            lapply(
                go_lists,
                function(x) {
                    names(x)[
                        !is.na(x) &
                        x < padj_cutoff
                    ]
                }
            ),
            use.names = FALSE
        ))
    )

    if (length(significant_go_ids) < 2) {

        message(sprintf(
            "Skipping comparison clustering: only %d significant GO term(s).",
            length(significant_go_ids)
        ))

        plot.new()
        text(
            0.5,
            0.5,
            sprintf(
                "Not enough significant GO terms for comparison (n = %d)",
                length(significant_go_ids)
            )
        )

        clusters <- data.frame(
            id = significant_go_ids,
            cluster = rep(
                NA_integer_,
                length(significant_go_ids)
            ),
            stringsAsFactors = FALSE
        )

    } else {

        clusters <- simplifyEnrichment::simplifyGOFromMultipleLists(
            lt = go_lists,
            padj_cutoff = padj_cutoff,
            ont = ontology,
            db = db,
            min_term = 2,
            heatmap_param = list(
                transform = function(x) {
                    -log10(pmax(
                        x,
                        .Machine$double.xmin
                    ))
                },
                name = "-log10(FDR)"
            ),
            show_barplot = TRUE,
            use_raster = FALSE
        )
    }

} else {

    stop(sprintf(
        "Unsupported analysis_mode: %s",
        analysis_mode
    ))
}

dev.off()


###############################################################################
# SAVE
###############################################################################

write.csv(
    clusters,
    output_clusters,
    row.names = FALSE
)