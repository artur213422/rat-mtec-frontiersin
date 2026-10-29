args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("Expected one JSON configuration argument.")
}


suppressPackageStartupMessages({
    library(jsonlite)
    library(limma)
})

config <- fromJSON(args[[1]])

input_file <- config$input_file
output_file <- config$output_file
data_type <- config$data_type
pseudocount <- as.numeric(config$pseudocount)

control_columns <- unlist(
    config$control_columns,
    use.names = FALSE
)

ko_columns <- unlist(
    config$ko_columns,
    use.names = FALSE
)

data <- read.csv(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

expression_columns <- c(
    control_columns,
    ko_columns
)

missing_columns <- setdiff(
    expression_columns,
    colnames(data)
)

if (length(missing_columns) > 0) {
    stop(
        paste(
            "Missing expression columns:",
            paste(missing_columns, collapse = ", ")
        )
    )
}

expression_matrix <- as.matrix(
    data[, expression_columns, drop = FALSE]
)

storage.mode(expression_matrix) <- "double"

if (any(!is.finite(expression_matrix))) {
    stop("Expression matrix contains non-finite values.")
}

if (data_type == "rma") {
    transformed_expression <- expression_matrix
    use_trend <- FALSE
} else if (data_type == "rpkm") {
    if (any(expression_matrix < 0)) {
        stop("RPKM values cannot be negative.")
    }

    transformed_expression <- log2(
        expression_matrix + pseudocount
    )

    use_trend <- TRUE
} else {
    stop(
        paste(
            "Unsupported data type:",
            data_type
        )
    )
}

rownames(transformed_expression) <- paste0(
    "row_",
    seq_len(nrow(transformed_expression))
)

condition <- factor(
    c(
        rep("control", length(control_columns)),
        rep("ko", length(ko_columns))
    ),
    levels = c("control", "ko")
)

design <- model.matrix(
    ~ 0 + condition
)

colnames(design) <- levels(condition)

fit <- lmFit(
    transformed_expression,
    design
)

contrast_matrix <- makeContrasts(
    ko - control,
    levels = design
)

fit <- contrasts.fit(
    fit,
    contrast_matrix
)

fit <- eBayes(
    fit,
    trend = use_trend,
    robust = TRUE
)

result <- data.frame(
    gene_symbol = data$gene_symbol,
    log2_fold_change = fit$coefficients[, 1],
    t_statistic = fit$t[, 1],
    p_val = fit$p.value[, 1],
    p_adj = p.adjust(
        fit$p.value[, 1],
        method = "BH"
    ),
    average_log2_expression = fit$Amean,
    stringsAsFactors = FALSE
)

write.csv(
    result,
    output_file,
    row.names = FALSE
)