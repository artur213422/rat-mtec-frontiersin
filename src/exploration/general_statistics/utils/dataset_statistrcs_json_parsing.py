import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from src.common.feratures.analysers.deg_models import (
    RankedGene,
    TopGenesResult,
)
from src.common.types.dataset import GeneReference
from src.exploration.general_statistics.utils.types import (
    DatasetStatistics,
    DegStatistics,
    ExpressionStatistics,
    ReplicateStatistics,
)


STATISTICS_SCHEMA_VERSION: int = 1


##########################################################################################
#                                         PUBLIC
##########################################################################################


def save_statistics(statistics: list[DatasetStatistics], file: Path) -> None:
    """
    Save dataset statistics to a versioned JSON file.

    Raises ValueError if the data contain NaN or infinity.
    """
    file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    payload: dict[str, Any] = {
        "schema_version": STATISTICS_SCHEMA_VERSION,
        "statistics": [
            asdict(dataset_statistics)
            for dataset_statistics in statistics
        ],
    }

    json_text: str = json.dumps(
        payload,
        indent=2,
        ensure_ascii=False,
        allow_nan=False,
    )

    temporary_file: Path = file.with_suffix(
        file.suffix + ".tmp"
    )

    temporary_file.write_text(
        json_text,
        encoding="utf-8",
    )

    # Works on both Windows and Linux.
    temporary_file.replace(file)


def read_statistics(file: Path) -> list[DatasetStatistics]:
    """
    Read dataset statistics from a versioned JSON file.
    """
    if not file.is_file():
        raise FileNotFoundError(
            f"Statistics file does not exist: {file}"
        )

    try:
        raw_data: Any = json.loads(
            file.read_text(encoding="utf-8")
        )
    except json.JSONDecodeError as error:
        raise ValueError(
            f"Invalid JSON in statistics file: {file}"
        ) from error

    root: dict[str, Any] = _require_dict(
        raw_data,
        context="statistics file root",
    )

    schema_version: int = _require_int(
        root,
        key="schema_version",
        context="statistics file root",
    )

    if schema_version != STATISTICS_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported statistics schema version: "
            f"{schema_version}. "
            f"Expected {STATISTICS_SCHEMA_VERSION}."
        )

    raw_statistics: list[Any] = _require_list(
        root,
        key="statistics",
        context="statistics file root",
    )

    return [
        _dataset_statistics_from_dict(item)
        for item in raw_statistics
    ]


##########################################################################################
#                                  DATASET STATISTICS
##########################################################################################


def _dataset_statistics_from_dict(
    data: Any,
) -> DatasetStatistics:
    record: dict[str, Any] = _require_dict(
        data,
        context="dataset statistics",
    )

    return DatasetStatistics(
        dataset_id=_require_str(
            record,
            key="dataset_id",
            context="dataset statistics",
        ),
        degs=_deg_statistics_from_dict(
            _require_value(
                record,
                key="degs",
                context="dataset statistics",
            )
        ),
        total_unique_genes=_require_int(
            record,
            key="total_unique_genes",
            context="dataset statistics",
        ),
        ctrl=_expression_statistics_from_dict(
            _require_value(
                record,
                key="ctrl",
                context="dataset statistics",
            )
        ),
        ko=_expression_statistics_from_dict(
            _require_value(
                record,
                key="ko",
                context="dataset statistics",
            )
        ),
        combined=_expression_statistics_from_dict(
            _require_value(
                record,
                key="combined",
                context="dataset statistics",
            )
        ),
        replicate_statistics=_replicate_statistics_from_dict(
            _require_value(
                record,
                key="replicate_statistics",
                context="dataset statistics",
            )
        ),
        top_15=_top_genes_result_from_dict(
            _require_value(
                record,
                key="top_15",
                context="dataset statistics",
            )
        ),
    )


##########################################################################################
#                                 EXPRESSION STATISTICS
##########################################################################################


def _expression_statistics_from_dict(
    data: Any,
) -> ExpressionStatistics:
    record: dict[str, Any] = _require_dict(
        data,
        context="expression statistics",
    )

    return ExpressionStatistics(
        mean_expression=_require_float(
            record,
            key="mean_expression",
            context="expression statistics",
        ),
        median_expression=_require_float(
            record,
            key="median_expression",
            context="expression statistics",
        ),
        mean_gene_sd=_require_float(
            record,
            key="mean_gene_sd",
            context="expression statistics",
        ),
        median_gene_sd=_require_float(
            record,
            key="median_gene_sd",
            context="expression statistics",
        ),
        mean_gene_cv=_require_float(
            record,
            key="mean_gene_cv",
            context="expression statistics",
        ),
        median_gene_cv=_require_float(
            record,
            key="median_gene_cv",
            context="expression statistics",
        ),
        expression_q25=_require_float(
            record,
            key="expression_q25",
            context="expression statistics",
        ),
        expression_q75=_require_float(
            record,
            key="expression_q75",
            context="expression statistics",
        ),
    )


##########################################################################################
#                                    DEG STATISTICS
##########################################################################################


def _deg_statistics_from_dict(
    data: Any,
) -> DegStatistics:
    record: dict[str, Any] = _require_dict(
        data,
        context="DEG statistics",
    )

    return DegStatistics(
        upregulated_count=_require_int(
            record,
            key="upregulated_count",
            context="DEG statistics",
        ),
        downregulated_count=_require_int(
            record,
            key="downregulated_count",
            context="DEG statistics",
        ),
        total_count=_require_int(
            record,
            key="total_count",
            context="DEG statistics",
        ),
        upregulated_fraction=_require_float(
            record,
            key="upregulated_fraction",
            context="DEG statistics",
        ),
        downregulated_fraction=_require_float(
            record,
            key="downregulated_fraction",
            context="DEG statistics",
        ),
        total_fraction=_require_float(
            record,
            key="total_fraction",
            context="DEG statistics",
        ),
        median_abs_log2_fc=_require_float(
            record,
            key="median_abs_log2_fc",
            context="DEG statistics",
        ),
        mean_abs_log2_fc=_require_float(
            record,
            key="mean_abs_log2_fc",
            context="DEG statistics",
        ),
        max_abs_log2_fc=_require_float(
            record,
            key="max_abs_log2_fc",
            context="DEG statistics",
        ),
        median_up_log2_fc=_require_float(
            record,
            key="median_up_log2_fc",
            context="DEG statistics",
        ),
        median_down_log2_fc=_require_float(
            record,
            key="median_down_log2_fc",
            context="DEG statistics",
        ),
    )


##########################################################################################
#                                REPLICATE STATISTICS
##########################################################################################


def _replicate_statistics_from_dict(
    data: Any,
) -> ReplicateStatistics:
    record: dict[str, Any] = _require_dict(
        data,
        context="replicate statistics",
    )

    return ReplicateStatistics(
        mean_ctrl_correlation=_require_float(
            record,
            key="mean_ctrl_correlation",
            context="replicate statistics",
        ),
        minimum_ctrl_correlation=_require_float(
            record,
            key="minimum_ctrl_correlation",
            context="replicate statistics",
        ),
        mean_ko_correlation=_require_float(
            record,
            key="mean_ko_correlation",
            context="replicate statistics",
        ),
        minimum_ko_correlation=_require_float(
            record,
            key="minimum_ko_correlation",
            context="replicate statistics",
        ),
    )


##########################################################################################
#                                      TOP GENES
##########################################################################################


def _top_genes_result_from_dict(
    data: Any,
) -> TopGenesResult:
    record: dict[str, Any] = _require_dict(
        data,
        context="top genes result",
    )

    raw_upregulated: list[Any] = _require_list(
        record,
        key="top_upregulated",
        context="top genes result",
    )

    raw_downregulated: list[Any] = _require_list(
        record,
        key="top_downregulated",
        context="top genes result",
    )

    return TopGenesResult(
        top_upregulated=[
            _ranked_gene_from_dict(item)
            for item in raw_upregulated
        ],
        top_downregulated=[
            _ranked_gene_from_dict(item)
            for item in raw_downregulated
        ],
    )


def _ranked_gene_from_dict(
    data: Any,
) -> RankedGene:
    record: dict[str, Any] = _require_dict(
        data,
        context="ranked gene",
    )

    return RankedGene(
        gene_reference=_gene_reference_from_dict(
            _require_value(
                record,
                key="gene_reference",
                context="ranked gene",
            )
        ),
        log2_fold_change=_require_float(
            record,
            key="log2_fold_change",
            context="ranked gene",
        ),
        p_val=_require_float(
            record,
            key="p_val",
            context="ranked gene",
        ),
    )


def _gene_reference_from_dict(
    data: Any,
) -> GeneReference:
    record: dict[str, Any] = _require_dict(
        data,
        context="gene reference",
    )

    gene_symbol: str = _require_str(
        record,
        key="gene_symbol",
        context="gene reference",
    )

    raw_ortholog_id: Any = _require_value(
        record,
        key="ortholog_id",
        context="gene reference",
    )

    if raw_ortholog_id is None:
        ortholog_id: str | None = None
    elif isinstance(raw_ortholog_id, str):
        ortholog_id = raw_ortholog_id
    else:
        raise ValueError(
            "Field 'ortholog_id' in gene reference must be "
            f"a string or null, got {type(raw_ortholog_id).__name__}."
        )

    return GeneReference(
        gene_symbol=gene_symbol,
        ortholog_id=ortholog_id,
    )


##########################################################################################
#                                  JSON VALIDATION
##########################################################################################


def _require_value(
    record: dict[str, Any],
    *,
    key: str,
    context: str,
) -> Any:
    if key not in record:
        raise ValueError(
            f"Missing required field '{key}' in {context}."
        )

    return record[key]


def _require_dict(
    value: Any,
    *,
    context: str,
) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(
            f"Expected a JSON object for {context}, "
            f"got {type(value).__name__}."
        )

    return value


def _require_list(
    record: dict[str, Any],
    *,
    key: str,
    context: str,
) -> list[Any]:
    value: Any = _require_value(
        record,
        key=key,
        context=context,
    )

    if not isinstance(value, list):
        raise ValueError(
            f"Field '{key}' in {context} must be a JSON array, "
            f"got {type(value).__name__}."
        )

    return value


def _require_str(
    record: dict[str, Any],
    *,
    key: str,
    context: str,
) -> str:
    value: Any = _require_value(
        record,
        key=key,
        context=context,
    )

    if not isinstance(value, str):
        raise ValueError(
            f"Field '{key}' in {context} must be a string, "
            f"got {type(value).__name__}."
        )

    return value


def _require_int(
    record: dict[str, Any],
    *,
    key: str,
    context: str,
) -> int:
    value: Any = _require_value(
        record,
        key=key,
        context=context,
    )

    # bool is a subclass of int, so reject it explicitly.
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(
            f"Field '{key}' in {context} must be an integer, "
            f"got {type(value).__name__}."
        )

    return value


def _require_float(
    record: dict[str, Any],
    *,
    key: str,
    context: str,
) -> float:
    value: Any = _require_value(
        record,
        key=key,
        context=context,
    )

    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValueError(
            f"Field '{key}' in {context} must be numeric, "
            f"got {type(value).__name__}."
        )

    return float(value)