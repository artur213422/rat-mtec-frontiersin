args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop("Usage: Rscript attach_gene_location.R '<config_json_file>'")
}

required_packages <- c("jsonlite")

for (package_name in required_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(
      "Required R package is not installed: ", package_name,
      "\nInstall it before running this script."
    )
  }
}

config <- jsonlite::fromJSON(args[[1]])
required_fields <- c("source_file", "output_file", "platform_package")

for (field in required_fields) {
  value <- config[[field]]

  if (is.null(value)) {
    stop("Missing required field in JSON: ", field)
  }

  if (!is.character(value) || length(value) != 1L || !nzchar(value)) {
    stop("'", field, "' must be one non-empty string")
  }
}

if (!file.exists(config$source_file)) {
  stop("Source file does not exist: ", config$source_file)
}

if (!requireNamespace(config$platform_package, quietly = TRUE)) {
  stop(
    "Required platform package is not installed: ",
    config$platform_package
  )
}

# Transcript-level platform-design packages contain a NetAffx table with
# genomic coordinates. This is the most direct source for ST-array probes.
coordinates_from_platform <- function(platform_package) {
  netaffx_file <- system.file(
    "extdata",
    "netaffxTranscript.rda",
    package = platform_package
  )

  if (!nzchar(netaffx_file)) {
    return(NULL)
  }

  data_environment <- new.env(parent = emptyenv())
  loaded_objects <- load(netaffx_file, envir = data_environment)

  if (!("netaffxTranscript" %in% loaded_objects)) {
    stop(
      "The platform package contains netaffxTranscript.rda but not ",
      "the expected 'netaffxTranscript' object: ", platform_package
    )
  }

  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop(
      "Required R package is not installed: Biobase\n",
      "It is needed to read coordinates from ", platform_package, "."
    )
  }

  netaffx <- Biobase::pData(data_environment$netaffxTranscript)
  required_columns <- c(
    "transcriptclusterid", "seqname", "strand", "start", "stop"
  )
  missing_columns <- setdiff(required_columns, colnames(netaffx))

  if (length(missing_columns) > 0L) {
    stop(
      "NetAffx metadata in ", platform_package,
      " is missing columns: ", paste(missing_columns, collapse = ", ")
    )
  }

  data.frame(
    ID_REF = as.character(netaffx$transcriptclusterid),
    chromosome = as.character(netaffx$seqname),
    strand = as.character(netaffx$strand),
    start = as.numeric(netaffx$start),
    stop = as.numeric(netaffx$stop),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# Classic Affymetrix CDF packages describe probe layout, not genomic
# locations. Their companion chip annotation packages expose CHR/CHRLOC.
annotation_package_for <- function(platform_package) {
  known_packages <- c(
    "mgu74av2cdf" = "mgu74av2.db",
    "pd.mogene.1.0.st.v1" = "mogene10sttranscriptcluster.db"
  )

  if (platform_package %in% names(known_packages)) {
    return(unname(known_packages[[platform_package]]))
  }

  if (grepl("cdf$", platform_package)) {
    return(sub("cdf$", ".db", platform_package))
  }

  stop(
    "Cannot determine an annotation package for platform package: ",
    platform_package
  )
}

coordinates_from_annotation_db <- function(platform_package) {
  annotation_package <- annotation_package_for(platform_package)

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop(
      "Required R package is not installed: AnnotationDbi\n",
      "It is needed to query ", annotation_package, "."
    )
  }

  if (!requireNamespace(annotation_package, quietly = TRUE)) {
    stop(
      "Required chip annotation package is not installed: ",
      annotation_package,
      "\nIt supplies genomic locations for ", platform_package, "."
    )
  }

  annotation_db <- get(
    annotation_package,
    envir = asNamespace(annotation_package)
  )
  available_columns <- AnnotationDbi::columns(annotation_db)
  required_columns <- c("CHR", "CHRLOC", "CHRLOCEND")

  if (all(required_columns %in% available_columns)) {
    locations <- suppressMessages(
      AnnotationDbi::select(
        annotation_db,
        keys = AnnotationDbi::keys(annotation_db, keytype = "PROBEID"),
        columns = required_columns,
        keytype = "PROBEID"
      )
    )

    locations <- locations[
      !is.na(locations$CHR) &
        !is.na(locations$CHRLOC) &
        !is.na(locations$CHRLOCEND),
      ,
      drop = FALSE
    ]

    # A probe can have multiple mappings. Keep the first complete mapping as
    # one coherent chromosome/start/end record.
    locations <- locations[!duplicated(locations$PROBEID), , drop = FALSE]
    probe_ids <- as.character(locations$PROBEID)
    chromosomes <- as.character(locations$CHR)
    raw_start <- as.numeric(locations$CHRLOC)
    raw_stop <- as.numeric(locations$CHRLOCEND)
  } else {
    # Current chip annotation packages may retain CHRLOC as legacy Bimap
    # objects while omitting it from columns(annotation_db).
    object_prefix <- sub("\\.db$", "", annotation_package)
    namespace <- asNamespace(annotation_package)
    start_object <- paste0(object_prefix, "CHRLOC")
    stop_object <- paste0(object_prefix, "CHRLOCEND")

    if (!exists(start_object, envir = namespace, inherits = FALSE) ||
        !exists(stop_object, envir = namespace, inherits = FALSE)) {
      stop(
        "Annotation package ", annotation_package,
        " does not expose CHRLOC and CHRLOCEND mappings."
      )
    }

    start_map <- as.list(get(start_object, envir = namespace))
    stop_map <- as.list(get(stop_object, envir = namespace))
    probe_ids <- intersect(names(start_map), names(stop_map))

    first_complete_mapping <- function(probe_id) {
      starts <- start_map[[probe_id]]
      stops <- stop_map[[probe_id]]

      if (length(starts) == 0L || length(stops) == 0L) {
        return(c(chromosome = NA_character_, start = NA, stop = NA))
      }

      start_chromosomes <- names(starts)
      stop_chromosomes <- names(stops)

      for (start_index in seq_along(starts)) {
        matching_stop <- which(
          stop_chromosomes == start_chromosomes[[start_index]] &
            !is.na(stops)
        )

        if (!is.na(starts[[start_index]]) && length(matching_stop) > 0L) {
          stop_index <- matching_stop[[1L]]
          return(c(
            chromosome = start_chromosomes[[start_index]],
            start = starts[[start_index]],
            stop = stops[[stop_index]]
          ))
        }
      }

      c(chromosome = NA_character_, start = NA, stop = NA)
    }

    mapping_rows <- lapply(probe_ids, first_complete_mapping)
    chromosomes <- vapply(mapping_rows, `[[`, character(1L), "chromosome")
    raw_start <- as.numeric(vapply(mapping_rows, `[[`, character(1L), "start"))
    raw_stop <- as.numeric(vapply(mapping_rows, `[[`, character(1L), "stop"))
  }

  data.frame(
    ID_REF = probe_ids,
    chromosome = ifelse(
      grepl("^chr", chromosomes, ignore.case = TRUE),
      chromosomes,
      paste0("chr", chromosomes)
    ),
    strand = ifelse(raw_start < 0 | raw_stop < 0, "-", "+"),
    start = pmin(abs(raw_start), abs(raw_stop)),
    stop = pmax(abs(raw_start), abs(raw_stop)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

message("Reading expression matrix: ", config$source_file)
expression_data <- read.csv(
  config$source_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!("ID_REF" %in% colnames(expression_data))) {
  stop("Expression matrix is missing required column: ID_REF")
}

if (anyDuplicated(expression_data$ID_REF)) {
  stop("Expression matrix contains duplicated ID_REF values")
}

expression_data$ID_REF <- as.character(expression_data$ID_REF)
location_columns <- c("chromosome", "strand", "start", "stop")
expression_data <- expression_data[
  setdiff(colnames(expression_data), location_columns)
]

message("Loading genomic coordinates for ", config$platform_package, "...")
coordinates <- coordinates_from_platform(config$platform_package)

if (is.null(coordinates)) {
  coordinates <- coordinates_from_annotation_db(config$platform_package)
}

if (anyDuplicated(coordinates$ID_REF)) {
  stop("Coordinate source contains duplicated ID_REF values after normalization")
}

coordinates$chromosome <- trimws(as.character(coordinates$chromosome))
coordinates$chromosome[coordinates$chromosome == ""] <- NA_character_
needs_chr_prefix <- !is.na(coordinates$chromosome) &
  !grepl("^chr", coordinates$chromosome, ignore.case = TRUE)
coordinates$chromosome[needs_chr_prefix] <- paste0(
  "chr", coordinates$chromosome[needs_chr_prefix]
)
coordinates$strand <- trimws(as.character(coordinates$strand))
coordinates$strand[!(coordinates$strand %in% c("+", "-"))] <- NA_character_
coordinates$start <- as.numeric(coordinates$start)
coordinates$stop <- as.numeric(coordinates$stop)

invalid_interval <-
  !is.na(coordinates$start) &
  !is.na(coordinates$stop) &
  coordinates$stop < coordinates$start

if (any(invalid_interval)) {
  interval_start <- coordinates$start[invalid_interval]
  coordinates$start[invalid_interval] <- coordinates$stop[invalid_interval]
  coordinates$stop[invalid_interval] <- interval_start
}

row_index <- match(expression_data$ID_REF, coordinates$ID_REF)
expression_data$chromosome <- coordinates$chromosome[row_index]
expression_data$strand <- coordinates$strand[row_index]
expression_data$start <- coordinates$start[row_index]
expression_data$stop <- coordinates$stop[row_index]

# Keep identifiers and gene metadata together, followed by coordinates and
# expression values. This also makes repeated runs idempotent.
leading_columns <- intersect(c("ID_REF", "gene_symbol"), colnames(expression_data))
remaining_columns <- setdiff(
  colnames(expression_data),
  c(leading_columns, location_columns)
)
expression_data <- expression_data[
  c(leading_columns, location_columns, remaining_columns)
]

output_directory <- dirname(config$output_file)

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

mapped_count <- sum(stats::complete.cases(
  expression_data[location_columns]
))
message(
  "Mapped complete genomic coordinates for ", mapped_count, " of ",
  nrow(expression_data), " probes."
)
message("Writing output file: ", config$output_file)
write.csv(
  expression_data,
  file = config$output_file,
  row.names = FALSE,
  quote = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

message("Done: ", config$output_file)
