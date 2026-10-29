from enum import Enum

import polars as pl

class WindowSource:
    DATASET_ID: str = "dataset_id"
    ORTHOLOG_ID: str = "ortholog_id"
    GENE_SYMBOL: str = "gene_symbol"
    CHROMOSOME: str = "chromosome"
    START: str = "start"
    STOP: str = "stop"

    T_STATISTIC: str = "t_statistic"
    LOG2_FOLD_CHANGE: str = "log2_fold_change"
    P_VALUE: str = "p_val"

    SCHEMA: pl.Schema = pl.Schema({
        DATASET_ID: pl.String,
        ORTHOLOG_ID: pl.String,
        GENE_SYMBOL: pl.String,
        CHROMOSOME: pl.String,
        START: pl.Int64,
        STOP: pl.Int64,
        T_STATISTIC: pl.Float64,
        LOG2_FOLD_CHANGE: pl.Float64,
        P_VALUE: pl.Float64,
    })


class AireDirection(str, Enum):
    POSITIVE = "positive"
    NEGATIVE = "negative"

class Windows:
    DATASET_ID: str = "dataset_id"
    WINDOW_ID: str = "window_id"
    CHROMOSOME: str = "chromosome"

    START: str = "window_start"
    STOP: str = "window_stop"
    GENE_COUNT: str = "window_gene_count"
    SPAN_BP: str = "span_bp"

    ORTHOLOG_IDS: str = "ortholog_ids"
    GENES: str = "genes"

    T_STATISTICS: str = "t_statistics"
    LOG2_FOLD_CHANGES: str = "log2_fold_changes"
    P_VALUES: str = "p_val"

    SCHEMA: pl.Schema = pl.Schema({
        DATASET_ID: pl.String,
        WINDOW_ID: pl.String,
        CHROMOSOME: pl.String,

        START: pl.Int64,
        STOP: pl.Int64,
        GENE_COUNT: pl.Int64,
        SPAN_BP: pl.Int64,

        ORTHOLOG_IDS: pl.List(pl.String),
        GENES: pl.List(pl.String),

        T_STATISTICS: pl.List(pl.Float64),
        LOG2_FOLD_CHANGES: pl.List(pl.Float64),
        P_VALUES: pl.List(pl.Float64),
    })


class ScoredWindows:
    SCORE: str = "score"
    AIRE_DIRECTION: str = "aire_direction"

    PASSENGER_GENES: str = "passenger_genes"
    PASSENGER_COUNT: str = "passenger_count"
    CORE_GENES: str = "core_genes"


    SCHEMA: pl.Schema = pl.Schema({
        **Windows.SCHEMA,

        SCORE: pl.Float64,
        AIRE_DIRECTION: pl.Enum(AireDirection),
        PASSENGER_GENES: pl.List(pl.String),
        PASSENGER_COUNT: pl.Int64,
        CORE_GENES: pl.List(pl.String)
    })


class Clusters:

    SCHEMA: pl.Schema = pl.Schema({
        **ScoredWindows.SCHEMA
    })