from src.common.paths import PROJECT_ROOT

DATA_DIR = PROJECT_ROOT / "data/temp/genomics/nearest_neighbor_clustering/data"
TOP_GENES_FILE = DATA_DIR / "top_aire_genes.csv"
NEIGHBOR_COUNTS_FILE = DATA_DIR / "aire_neighbors_count.csv"
PERMUTATIONS_FILE = DATA_DIR / "random_nearest_neighbor_permutations.csv"

TOP_UP_GENES: int = 200 # 297 is max
TOP_DOWN_GENES: int = 200 # 297 is max


LOG2_FC_THRESHOLD: float = 0.5

PERMUTATIONS: int = 3000

NEIGHBOR_DISTANCE_KB: int = 200
