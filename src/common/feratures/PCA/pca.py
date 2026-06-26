import polars as pl

from src.PCA.models import ExpressionTable, PcaResult
from src.PCA.run_pca_analysis import prepare_pca_input, run_pca_analysis

file: str = r"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\shared_orthologs\mapped\rat\mtec_high_mapped_mouse_orthologs.csv"


data: pl.DataFrame = pl.read_csv(file)

expression_table: ExpressionTable = prepare_pca_input(data)


pca_result: PcaResult = run_pca_analysis(expression_table)

print("Done")
print()
print("PC names:")
print(pca_result.sample_coordinates.pc_names)

print()
print("Explained variance:")
for pc_name, ratio in zip(
    pca_result.explained_variance.pc_names,
    pca_result.explained_variance.explained_variance_ratio,
    strict=True,
):
    print(f"{pc_name}: {ratio * 100:.2f}%")

print()
print("Sample coordinates shape:")
print(pca_result.sample_coordinates.sample_by_pc_coordinates.shape)

print()
print("Feature loadings shape:")
print(pca_result.feature_loadings.gene_by_pc_loadings.shape)

print()
print("Sample coordinates:")
for sample_id, coordinates in zip(
    pca_result.sample_coordinates.sample_ids,
    pca_result.sample_coordinates.sample_by_pc_coordinates,
    strict=True,
):
    pc1: float = float(coordinates[0])
    pc2: float = float(coordinates[1])
    print(f"{sample_id}: PC1={pc1:.3f}, PC2={pc2:.3f}")