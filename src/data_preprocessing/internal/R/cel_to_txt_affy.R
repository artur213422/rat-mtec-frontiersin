args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript cel_to_txt_affy.R '<config_json>'")
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop(
    "Required R package is not installed: jsonlite\n",
    "Install it with: install.packages('jsonlite')"
  )
}

config_json <- args[[1]]
config <- jsonlite::fromJSON(config_json)

required_base_packages <- c(
  "affy",
  "Biobase"
)

for (pkg in required_base_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste(
      "Required R package is not installed:",
      pkg,
      "\nInstall it before running this script."
    ))
  }
}

required_fields <- c(
  "cel_files",
  "sample_names",
  "output_file",
  "platform_package"
)

for (field in required_fields) {
  if (is.null(config[[field]])) {
    stop(paste("Missing required field in JSON:", field))
  }
}

if (!is.character(config$cel_files)) {
  stop("'cel_files' must be a character vector")
}

if (!is.character(config$sample_names)) {
  stop("'sample_names' must be a character vector")
}

if (length(config$cel_files) != length(config$sample_names)) {
  stop("'cel_files' and 'sample_names' must have the same length")
}

if (length(config$cel_files) < 2) {
  stop("RMA should be run on at least two CEL files. Provide all CEL files for the dataset.")
}

if (!is.character(config$output_file) || length(config$output_file) != 1) {
  stop("'output_file' must be one string")
}

if (!is.character(config$platform_package) || length(config$platform_package) != 1) {
  stop("'platform_package' must be one string")
}

missing_files <- config$cel_files[!file.exists(config$cel_files)]

if (length(missing_files) > 0) {
  stop(paste(
    "The following CEL files do not exist:",
    paste(missing_files, collapse = "\n"),
    sep = "\n"
  ))
}

if (!requireNamespace(config$platform_package, quietly = TRUE)) {
  stop(paste(
    "Required CDF/platform package is not installed:",
    config$platform_package,
    "\nInstall it before running this script."
  ))
}

suppressPackageStartupMessages({
  library(affy)
  library(config$platform_package, character.only = TRUE)
})

message("Reading CEL files...")
raw_data <- affy::ReadAffy(filenames = config$cel_files)

Biobase::sampleNames(raw_data) <- config$sample_names

message("Running RMA...")
eset <- affy::rma(raw_data)

message("Extracting expression matrix...")
expr_matrix <- Biobase::exprs(eset)

colnames(expr_matrix) <- config$sample_names

output_df <- data.frame(
  ID_REF = rownames(expr_matrix),
  expr_matrix,
  check.names = FALSE
)

output_dir <- dirname(config$output_file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

quote_output <- TRUE

if (!is.null(config$quote_output)) {
  quote_output <- isTRUE(config$quote_output)
}

message("Writing output file...")
write.table(
  output_df,
  file = config$output_file,
  sep = "\t",
  quote = quote_output,
  row.names = FALSE,
  col.names = TRUE,
  fileEncoding = "UTF-8"
)

message("Done: ", config$output_file)