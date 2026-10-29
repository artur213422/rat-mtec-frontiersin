from __future__ import annotations

import re

import polars as pl


_CHEMOKINE_SYMBOL = re.compile(
    r"^(?P<prefix>ccl|cxcl|cx3cl|xcl)(?P<number>\d+)[a-z]*$",
    flags=re.IGNORECASE,
)

_CANONICAL_PREFIXES: dict[str, str] = {
    "ccl": "Ccl",
    "cxcl": "Cxcl",
    "cx3cl": "Cx3cl",
    "xcl": "Xcl",
}


def _chemokine_family_name(gene_symbol: str) -> str:
    """Return a family label only for an ambiguous, single-family annotation."""
    if gene_symbol.strip().lower() in {"ccl27a", "ccl27b", "ccl27al"}:
        return "Ccl27_family"

    if "///" not in gene_symbol:
        return gene_symbol

    families: set[tuple[str, str]] = set()
    for symbol in gene_symbol.split("///"):
        match = _CHEMOKINE_SYMBOL.fullmatch(symbol.strip())
        if match is not None:
            families.add(
                (
                    match.group("prefix").lower(),
                    match.group("number"),
                )
            )

    # Do not guess when the annotation contains no chemokine, or more than one
    # chemokine family (for example, an annotation spanning Ccl19 and Ccl21).
    if len(families) != 1:
        return gene_symbol

    prefix, number = families.pop()
    return f"{_CANONICAL_PREFIXES[prefix]}{number}_family"


def try_rescue_chemokines(tec_df: pl.DataFrame) -> pl.DataFrame:
    """Replace ambiguous chemokine annotations with non-GO family labels.

    For example, ``Ccl21a /// Ccl21b /// Ccl21c /// Gm10591`` becomes
    ``Ccl21_family``. The underlying probe remains one family-level
    measurement and should therefore be excluded from ordinary GO GSEA.
    """
    if "gene_symbol" not in tec_df.columns:
        raise ValueError("Expected a 'gene_symbol' column.")

    rescued_symbols = [
        _chemokine_family_name(symbol)
        if symbol is not None
        else None
        for symbol in tec_df["gene_symbol"].to_list()
    ]

    return tec_df.with_columns(pl.Series("gene_symbol", rescued_symbols))
