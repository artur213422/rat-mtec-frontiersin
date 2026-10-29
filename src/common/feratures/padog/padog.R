source(here::here("src", "common", "feratures", "R", "json_config_handler.R"))
source(here::here("src", "common", "feratures", "padog", "helpers", "padog_helpers.R"))

###############################################################################
# PREPARE DATA
###############################################################################

config <- read_json_config()

expression_matrix_file <- get_required_string(config, "expression_matrix_file")
gmt_file <- get_required_string(config, "gmt_file")

ko_group <- get_required_character_vector(config, "ko_group")
control_group <- get_required_character_vector(config, "control_group")

expression_matrix <- load_expression_matrix(control_group, ko_group, expression_matrix_file)
gene_sets <- load_gene_sets(gmt_file)
padog_group <- load_padog_group(control_group, ko_group)

result_file <- get_required_string(config, "result_file")

###############################################################################
# PADOG
###############################################################################

set.seed(123)

padog_result <- PADOG::padog(
  esetm = expression_matrix,
  group = padog_group,
  paired = FALSE,
  block = NULL,
  gslist = gene_sets,
  annotation = NULL,
  NI = 1000,
  Nmin = 3,
  plots = FALSE,
  verbose = TRUE,
  parallel = FALSE,
  dseed = 123
)

###############################################################################
# SAVE
###############################################################################

write.csv(
  padog_result,
  file = result_file,
  row.names = FALSE,
  quote = TRUE,
  na = "",
  fileEncoding = "UTF-8"
)

message("Done: ", result_file)