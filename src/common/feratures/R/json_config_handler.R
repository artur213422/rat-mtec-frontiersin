

read_json_config <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1L) {
    stop("Expected JSON config as the first command-line argument.")
  }

  config <- jsonlite::fromJSON(args[[1]])

  return(config)
}

get_required <- function(config, name) {
  if (!(name %in% names(config))) {
    stop(sprintf("Missing required config field: %s", name))
  }

  return(config[[name]])
}

get_required_double <- function(config, name) {
  value <- get_required(config, name)

  if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop(sprintf("Config field '%s' must be a single character string, got: %s", name, paste(class(value), collapse = ", ")))
  }

  return(value)
}

get_required_string <- function(config, name) {
  value <- get_required(config, name)

  if (!is.character(value) || length(value) != 1L) {
    stop(sprintf("Config field '%s' must be a single character string, got: %s", name, paste(class(value), collapse = ", ")))
  }

  return(value)
}


get_required_character_vector <- function(config, name) {
  value <- get_required(config, name)

  if (!is.character(value) || length(value) < 1L) {
    stop(sprintf("Config field '%s' must be a non-empty character vector, got: %s", name, paste(class(value), collapse = ", ")))
  }

  return(value)
}