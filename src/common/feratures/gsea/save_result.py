from pathlib import Path
import polars as pl

from src.common.feratures.gsea.models import AnalysedGseaResult


def save_gsea_result(out: Path, result: list[AnalysedGseaResult]):
    out.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | float | int]] = []

    for analysed_result in result:
        for gsea_result in analysed_result.gsea_result:
            rows.append(
                {
                    "dataset_id": analysed_result.dataset_id,
                    "term": gsea_result.term,
                    "es": gsea_result.es,
                    "nes": gsea_result.nes,
                    "nominal_p_value": gsea_result.nominal_p_value,
                    "fdr_q_value": gsea_result.fdr_q_value,
                    "leading_edge_gene_count": len(gsea_result.leading_edge_genes),
                    "leading_edge_genes": ";".join(gsea_result.leading_edge_genes),
                }
            )

    schema = {
        "dataset_id": pl.String,
        "term": pl.String,
        "es": pl.Float64,
        "nes": pl.Float64,
        "nominal_p_value": pl.Float64,
        "fdr_q_value": pl.Float64,
        "leading_edge_gene_count": pl.Int64,
        "leading_edge_genes": pl.String,
    }

    result_df = pl.DataFrame(rows, schema=schema)

    result_df.write_csv(out)

    print(f"Saved GSEA results: {out}")
    print(f"Rows: {result_df.height}")