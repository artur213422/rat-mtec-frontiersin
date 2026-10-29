


load_expression_matrix <- function(control_group, ko_group, expression_matrix_file){
  sample_names <- c(control_group, ko_group)

  if (anyDuplicated(sample_names)) {
    stop("Duplicate sample names found in control_group/ko_group.")
  }

  df <- read.csv(
    expression_matrix_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(df) < 2) {
    stop("Expression matrix must contain a gene ID column and at least one sample column.")
  }

   gene_ids <- as.character(df[[1]])

    if (anyNA(gene_ids) || any(gene_ids == "")) {
        stop("Expression matrix contains missing or empty gene IDs.")
    }

    if (anyDuplicated(gene_ids)) {
        stop("Expression matrix contains duplicated gene IDs.")
    }

    missing_samples <- setdiff(sample_names, colnames(df))

    if (length(missing_samples) > 0) {
        stop("Samples missing from expression matrix: ", paste(missing_samples, collapse = ", "))
    }

    expression_df <- df[, sample_names, drop = FALSE]

    if (!all(vapply(expression_df, is.numeric, logical(1)))) {
        stop("Some selected expression columns are not numeric.")
    }

    expression_matrix <- as.matrix(expression_df)

    rownames(expression_matrix) <- gene_ids

    if (anyNA(expression_matrix)) {
        stop("Expression matrix contains NA values.")
    }

    return(expression_matrix)
}

load_gene_sets <- function(gmt_file){
  lines <- readLines(
    gmt_file,
    warn = FALSE,
    encoding = "UTF-8"
  )

  lines <- lines[nzchar(trimws(lines))]

  if (length(lines) == 0) {
    stop("GMT file is empty.")
  }

  gene_sets <- lapply(
    lines,
    function(line) {
      fields <- strsplit(
        line,
        "\t",
        fixed = TRUE
      )[[1]]

      if (length(fields) < 3) {
        stop("Invalid GMT line. Expected: gene_set<TAB>description<TAB>gene1...")
      }

      genes <- fields[3:length(fields)]
      genes <- genes[nzchar(genes)]

      unique(genes)
    }
  )

  gene_set_names <- vapply(
    lines,
    function(line) {
      strsplit(
        line,
        "\t",
        fixed = TRUE
      )[[1]][1]
    },
    character(1)
  )

  if (anyDuplicated(gene_set_names)) {
    stop("GMT file contains duplicated gene-set names.")
  }

  names(gene_sets) <- gene_set_names

  return(gene_sets)
}

load_padog_group <- function(control_group, ko_group){
  if (length(control_group) == 0) {
    stop("control_group must contain at least one sample.")
  }

  if (length(ko_group) == 0) {
    stop("ko_group must contain at least one sample.")
  }

  if (any(grepl("ko", control_group, ignore.case = TRUE))) {
    stop(
      "Sanity test failed: control_group contains a sample name with 'ko': ",
      paste(
        control_group[
          grepl("ko", control_group, ignore.case = TRUE)
        ],
        collapse = ", "
      )
    )
  }

  if (any(grepl("control|ctrl", ko_group, ignore.case = TRUE))) {
    stop(
      "Sanity test failed: ko_group contains a sample name with 'control' or 'ctrl': ",
      paste(
        ko_group[
          grepl("ko", control_group, ignore.case = TRUE)
        ],
        collapse = ", "
      )
    )
  }

  if (length(intersect(control_group, ko_group)) > 0) {
    stop("Samples occur in both control_group and ko_group: ", paste(intersect(control_group, ko_group), collapse = ", "))
  }

  return(
    c(
      rep("c", length(control_group)),
      rep("d", length(ko_group))
    )
  )
}