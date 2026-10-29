from __future__ import annotations

import polars as pl


_RESCUABLE_KALLIKREINS: frozenset[str] = frozenset(
    {
        "klk1b24",
        "klk1b26",
        "klk1b27",
    }
)


def _rescue_kallikrein_symbol(gene_symbol: str) -> str:
    """Rescue one specifically requested kallikrein from a mixed annotation."""
    if "///" not in gene_symbol:
        return gene_symbol

    matches = {
        symbol.strip()
        for symbol in gene_symbol.split("///")
        if symbol.strip().lower() in _RESCUABLE_KALLIKREINS
    }

    # Preserve genuinely ambiguous annotations, including probes that contain
    # multiple requested kallikreins.
    if len(matches) != 1:
        return gene_symbol

    return matches.pop()


def try_rescue_kallikreins(tec_df: pl.DataFrame) -> pl.DataFrame:
    """Rescue Klk1b24, Klk1b26, and Klk1b27 from multi-gene annotations.

    This intentionally does not rescue broad mixed probes such as
    ``Klk1 /// Klk1b1 /// Klk1b4 /// Klk1b7-ps /// Klk1b9``.
    """
    if "gene_symbol" not in tec_df.columns:
        raise ValueError("Expected a 'gene_symbol' column.")

    rescued_symbols = [
        _rescue_kallikrein_symbol(symbol)
        if symbol is not None
        else None
        for symbol in tec_df["gene_symbol"].to_list()
    ]

    return tec_df.with_columns(pl.Series("gene_symbol", rescued_symbols))
