import math
from dataclasses import dataclass

import polars as pl

from src.common.feratures.analysers.deg_models import RankedGene, TopGenesResult
from src.common.feratures.gsea.rnk import get_genotype_expression_values
from src.common.feratures.math.fold_change import calculate_log2_fold_change
from src.common.feratures.math.p_val import calculate_raw_p_val
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
    required_columns: set[str] = {
        "gene_symbol",
        "log2_fold_change",
        "p_val",
    }

    missing_columns: set[str] = required_columns.difference(tec_df.columns)

    if missing_columns:
        missing: str = ", ".join(sorted(missing_columns))
        raise ValueError(f"Missing required columns: {missing}")

    upregulated: list[str] = []
    downregulated: list[str] = []

    for row in tec_df.iter_rows(named=True):
        gene_symbol_raw: object = row["gene_symbol"]
        log2_fc_raw: object = row["log2_fold_change"]
        p_val: object = row["p_adj"]

        if gene_symbol_raw is None or log2_fc_raw is None or p_val is None:
            continue

        gene_symbol: str = str(gene_symbol_raw).strip()
        log2_fc: float = float(log2_fc_raw)
        p_val: float = float(p_val)

        if not gene_symbol or not math.isfinite(log2_fc) or not math.isfinite(p_val):
            continue

        if not is_deg(p_val=p_val, fold_change=log2_fc):
            continue

        if log2_fc > 0.0:
            upregulated.append(gene_symbol)
        elif log2_fc < 0.0:
            downregulated.append(gene_symbol)

    return DatasetDegs(
        dataset_id=dataset_id,
        color=color,
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


def print_deg_summary(datasets: list[DatasetDegs],show_genes: bool = False,
) -> None:
    print("\nDEG SUMMARY")
    print("-" * 60)

    for dataset in datasets:
        up_count: int = len(dataset.upregulated)
        down_count: int = len(dataset.downregulated)
        total_count: int = up_count + down_count

        print(f"\nDataset: {dataset.dataset_id}")
        print(f"  Upregulated:   {up_count}")
        print(f"  Downregulated: {down_count}")
        print(f"  Total:         {total_count}")

        if show_genes:
            print(
                "  Upregulated genes: "
                + (", ".join(dataset.upregulated) or "None")
            )
            print(
                "  Downregulated genes: "
                + (", ".join(dataset.downregulated) or "None")
            )