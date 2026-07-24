from dataclasses import dataclass

import polars as pl

from src.common.feratures.analysers.deg_models import RankedGene, TopGenesResult
from src.common.feratures.gsea.rnk import get_genotype_expression_values
from src.common.feratures.math.p_val import calculate_raw_p_val, calculate_log2_fold_change
from src.common.types.dataset import GeneReference
from src.common.types.enums import Genotypes, Colors


@dataclass
class DatasetDegs:
    dataset_id: str
    color: Colors
    upregulated: list[str]
    downregulated: list[str]


def is_deg(p_val: float, fold_change: float, threshold: float = 1.0, p_val_threshold: float = 0.05) -> bool:
    return p_val < p_val_threshold and abs(fold_change) >= threshold

def get_degs(dataset_id: str, tec_df: pl.DataFrame, color:Colors) -> DatasetDegs:
    upregulated: list[str] = []
    downregulated: list[str] = []

    for row in tec_df.iter_rows(named=True):
        ctrl_values: list[float] = get_genotype_expression_values(row=row, gt=Genotypes.CTRL)
        ko_values: list[float] = get_genotype_expression_values(row=row, gt=Genotypes.KO)

        p_val: float = calculate_raw_p_val(ctrl_values, ko_values)
        log2_fc: float = calculate_log2_fold_change(ctrl_values, ko_values)

        if not is_deg(p_val=p_val, fold_change=log2_fc):
            continue

        candidate: str = str(row["gene_symbol"])

        if log2_fc > 0.0:
            upregulated.append(candidate)

        elif log2_fc < 0.0:
            downregulated.append(candidate)

    return DatasetDegs(
        dataset_id=dataset_id,
        color= color,
        upregulated=upregulated,
        downregulated=downregulated,
    )


def extract_top_genes(top: int, data: pl.DataFrame, threshold: float = 1.0,  p_val_threshold: float = 0.05) -> TopGenesResult:
    upregulated: list[RankedGene] = []
    downregulated: list[RankedGene] = []

    for row in data.iter_rows(named=True):
        ctrl_values: list[float] = get_genotype_expression_values(row, Genotypes.CTRL)
        ko_values: list[float] = get_genotype_expression_values(row, Genotypes.KO)

        p_val: float = calculate_raw_p_val(ctrl_values, ko_values)
        log2_fold_change: float = calculate_log2_fold_change(ctrl_values, ko_values)

        if not is_deg(
            p_val=p_val,
            fold_change=log2_fold_change,
            threshold=threshold,
            p_val_threshold=p_val_threshold,
        ):
            continue

        ranked_gene = RankedGene(
            gene_reference=GeneReference.from_row(row),
            log2_fold_change=log2_fold_change,
            p_val=p_val,
        )

        if log2_fold_change > 0:
            upregulated.append(ranked_gene)
        elif log2_fold_change < 0:
            downregulated.append(ranked_gene)

    upregulated = sorted(
        upregulated,
        key=lambda gene: gene.log2_fold_change,
        reverse=True,
    )

    downregulated = sorted(
        downregulated,
        key=lambda gene: gene.log2_fold_change,
    )

    return TopGenesResult(
        top_upregulated=upregulated[:top],
        top_downregulated=downregulated[:top],
    )