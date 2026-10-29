
import polars as pl

from src.analysis.mtec_high.genomics.window_sliding.columns import Windows, WindowSource

def enumerate_genomic_windows(tec_df: pl.DataFrame, dataset_id: str, min_size: int, max_size: int) -> pl.DataFrame:
    if min_size < 1:
        raise ValueError("min_size must be >= 1.")

    if max_size < min_size:
        raise ValueError("max_size must be >= min_size.")

    windows: list[dict[str, object]] = []
    window_index: int = 0

    chromosome_dfs: list[pl.DataFrame] = tec_df.partition_by(WindowSource.CHROMOSOME, maintain_order=True)

    for chromosome_df in chromosome_dfs:
        chromosome_df = chromosome_df.sort(WindowSource.START, WindowSource.STOP,)

        chromosome: str = chromosome_df.get_column(WindowSource.CHROMOSOME).item(0)

        gene_count: int = chromosome_df.height

        for start_index in range(gene_count):
            for window_size in range(min_size, max_size + 1):
                stop_index: int = start_index + window_size

                if stop_index > gene_count:
                    break

                window_df: pl.DataFrame = chromosome_df.slice(start_index, window_size,)

                window_start: int = window_df.get_column(WindowSource.START).min()

                window_stop: int = window_df.get_column(WindowSource.STOP).max()

                windows.append({
                    Windows.DATASET_ID: dataset_id,
                    Windows.WINDOW_ID: f"{dataset_id}_{window_index}",
                    Windows.CHROMOSOME: chromosome,

                    Windows.START: window_start,
                    Windows.STOP: window_stop,
                    Windows.GENE_COUNT: window_size,
                    Windows.SPAN_BP: window_stop - window_start,

                    Windows.ORTHOLOG_IDS: window_df.get_column(WindowSource.ORTHOLOG_ID).to_list(),

                    Windows.GENES: window_df.get_column(WindowSource.GENE_SYMBOL).to_list(),

                    Windows.T_STATISTICS: window_df.get_column(WindowSource.T_STATISTIC).to_list(),

                    Windows.LOG2_FOLD_CHANGES: window_df.get_column(WindowSource.LOG2_FOLD_CHANGE).to_list(),

                    Windows.P_VALUES: window_df.get_column(WindowSource.P_VALUE).to_list(),
                })

                window_index += 1

    return pl.DataFrame(windows, schema=Windows.SCHEMA)