from typing import Literal





MINIMUM_EXPRESSION_LOG2 = 5.0 # Similar to my TSA_MEAN_ENRICHMENT_LOG2

DROP_GENES_WITH_NO_ORTHOLOGS: bool = True


PRIMARY_SELECTION_RULE: Literal[
    "largest_j",
    "smallest_j",
] = "largest_j"

MINIMUM_EXPRESSION_LINEAR: float = 2.0 ** MINIMUM_EXPRESSION_LOG2



##########################################################################################
# TSA
##########################################################################################

MAX_SELECTED_TSA_GROUPS = 5


##########################################################################################
# HOUSE KEEPERS
##########################################################################################

MIN_EXPRESSED_GROUP_FRACTION: float = 0.90

MAX_SELECTED_HOUSEKEEPING_GROUPS = 6
