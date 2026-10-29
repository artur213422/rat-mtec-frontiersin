source(here::here("src", "common", "feratures", "R", "json_config_handler.R"))
source(here::here("src", "common", "feratures", "rrvgo", "helpers", "rrvgo_helpers.R"))
source(here::here("src", "common", "feratures", "simplify_enrichment", "helpers", "simplify_enrichment_helpers.R"))

###############################################################################
# PREPARE DATA
################################################################################

config <- read_json_config()

gsea_file <- get_required_string(config, "gsea_file")
reference_species <- get_required_string(config, "reference_species")
ontology <- get_required_ontology(config, "ontology")
padj_cutoff <- get_required_double(config, "padj_cutoff")
output_folder <- get_required_string(config, "output_folder")

db <- get_database(reference_species)
gsea_df <- read.csv(gsea_file)

output_reduced_terms <- file.path(output_folder, "output_reduced_terms.csv")
output_heatmap_svg <- file.path(output_folder, "output_heatmap.svg")
output_wordcloud_svg <- file.path(output_folder, "output_wordcloud.svg")
output_treemap_svg <- file.path(output_folder, "output_treemap.svg")
output_scatter_svg <- file.path(output_folder, "output_scatter.svg")

###############################################################################
# PREPARE SIGNIFICANT GO TERMS
###############################################################################


significant_df <- gsea_df[
    !is.na(gsea_df$fdr_q_value) &
    gsea_df$fdr_q_value <= padj_cutoff,
]

go_ids <- sub(
    "^\\s*(GO:[0-9]{7}).*$",
    "\\1",
    significant_df$term
)

valid <- grepl(
    "^GO:[0-9]{7}$",
    go_ids
)

significant_df <- significant_df[valid, ]
go_ids <- go_ids[valid]

################################################################################
# CREATE RRVGO SCORES
###############################################################################

scores <- setNames(
    -log10(
        pmax(
            significant_df$fdr_q_value,
            .Machine$double.xmin
        )
    ),
    go_ids
)

# If the same GO term occurs multiple times,
# keep its strongest significance score.
scores <- tapply(
    scores,
    names(scores),
    max
)

go_ids <- names(scores)


###############################################################################
# CALCULATE SIMILARITY MATRIX
###############################################################################

if (length(go_ids) < 2) {
    stop(sprintf(
        "Not enough significant GO terms for rrvgo (n = %d).",
        length(go_ids)
    ))
}

simMatrix <- rrvgo::calculateSimMatrix(
    go_ids,
    orgdb = db,
    ont = ontology,
    method = "Rel"
)


###############################################################################
# REDUCE GO TERMS
###############################################################################

# Keep scores synchronized with GO terms in the similarity matrix.
scores <- scores[rownames(simMatrix)]

reducedTerms <- rrvgo::reduceSimMatrix(
    simMatrix,
    scores = scores,
    threshold = 0.7,
    orgdb = db
)


###############################################################################
# PLOT - HEATMAP
###############################################################################

svglite::svglite(output_heatmap_svg, width = 16, height = 10)

rrvgo::heatmapPlot(
    simMatrix,
    reducedTerms,
    annotateParent = TRUE,
    annotationLabel = "parentTerm",
    fontsize = 6
)

dev.off()

###############################################################################
# PLOT - WORD CLOUD
###############################################################################

svglite::svglite(output_wordcloud_svg, width = 16, height = 10)

rrvgo::wordcloudPlot(reducedTerms, min.freq = 1, colors = "black")

dev.off()

###############################################################################
# PLOT - TREEMAP
###############################################################################

svglite::svglite(output_treemap_svg, width = 16, height = 10)

rrvgo::treemapPlot(reducedTerms)

dev.off()

###############################################################################
# PLOT - SCATTER
###############################################################################

if (nrow(simMatrix) >= 3) {
    svglite::svglite(output_scatter_svg, width = 16, height = 10)
    print(rrvgo::scatterPlot(simMatrix, reducedTerms))
    dev.off()
} else {
    message(sprintf(
        "Skipping scatter plot: at least 3 GO terms are required (n = %d).",
        nrow(simMatrix)
    ))
}

###############################################################################
# SAVE
###############################################################################

write.csv(reducedTerms, output_reduced_terms, row.names = FALSE)