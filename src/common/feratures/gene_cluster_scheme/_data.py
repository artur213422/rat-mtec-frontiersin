from __future__ import annotations

import polars as pl


class _GeneClusterDataMixin:
    data: pl.DataFrame
    value_col: str
    gene_col: str
    chromosome_col: str
    start_col: str
    stop_col: str
    cluster_genes: list[str]

    def _validate_columns(self, data: pl.DataFrame) -> None:
        required_columns: set[str] = {
            self.value_col,
            self.gene_col,
            self.chromosome_col,
            self.start_col,
            self.stop_col,
        }

        missing_columns: set[str] = required_columns - set(data.columns)


        if missing_columns:
            raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    def _prepare_data(self, data: pl.DataFrame) -> pl.DataFrame:
        return (
            data
            .filter(pl.col(self.gene_col).is_in(self.cluster_genes))
            .with_columns(
                pl.min_horizontal(self.start_col, self.stop_col).alias("_plot_start"),
                pl.max_horizontal(self.start_col, self.stop_col).alias("_plot_stop"),
            )
            .with_columns(((pl.col("_plot_start") + pl.col("_plot_stop")) / 2).alias("_midpoint")
            )
            .sort(
                self.chromosome_col,
                "_plot_start",
            )
        )

    def _validate_cluster(self) -> None:
        if self.data.height == 0:
            raise ValueError("None of cluster_genes were found in data.")

        found_genes: set[str] = set(
            self.data
            .get_column(self.gene_col)
            .cast(pl.String)
            .to_list()
        )

        requested_genes: set[str] = set(self.cluster_genes)

        missing_genes: set[str] = requested_genes - found_genes

        if missing_genes:
            raise ValueError(f"Cluster genes missing from data: {sorted(missing_genes)}")

        if self.data.get_column(self.chromosome_col).null_count():
            raise ValueError("Chromosome values must not be null.")
        for column in ("_plot_start", "_plot_stop"):
            values = self.data.get_column(column).cast(pl.Float64)
            if values.null_count() or not values.is_finite().all():
                raise ValueError("Genomic coordinates must be finite and non-null.")

    def _chromosome_segments(self) -> list[tuple[str, list[pl.DataFrame]]]:
        """Split at uncovered gaps, tracking the end of overlapping genes."""
        result = []
        for chromosome in self.chromosomes:
            data = self.data.filter(pl.col(self.chromosome_col) == chromosome)
            segments = []
            first = 0
            end = float(data["_plot_stop"][0])
            for index, row in enumerate(data.iter_rows(named=True)):
                start = float(row["_plot_start"])
                if start - end > self.gap_break_threshold_kb * 1000:
                    segments.append(data.slice(first, index - first))
                    first = index
                end = max(end, float(row["_plot_stop"]))
            segments.append(data.slice(first))
            result.append((chromosome, segments))
        return result
