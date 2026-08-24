from src.common.dataframe_ops import get_genotype_expression_values, get_custom_genotype_expression_values
from src.common.feratures.math.fold_change import calculate_log2_fold_change
from src.common.feratures.math.is_significant import is_significant
from src.common.feratures.math.logarithm import log2_transform
from src.common.feratures.math.p_val import calculate_raw_p_val
from src.common.types.enums import Genotypes

from dataclasses import dataclass

@dataclass(frozen=True, slots=True)
class CandidateResult:
    gene_symbol: str
    ifn: str
    p_value: float
    log2_fold_change: float


def is_candidate(row: dict[str, object], ifns: list[str]) -> CandidateResult | None:
    """
      Return the IFN name when exactly one IFN significantly upregulates the gene.

      Return None when:
      - no IFN is significant;
      - multiple IFNs are significant;
      - the only significant response is downregulated.
      """
    gene_value: object | None = row.get("gene_symbol", row.get("gene"))

    if gene_value is None:
        raise KeyError("The row must contain either 'gene_symbol' or 'gene'.")

    gene_symbol: str = str(gene_value)

    raw_ctrl_values: list[float] = get_genotype_expression_values(row=row, gt=Genotypes.CTRL)
    log2_ctrl_values: list[float] = log2_transform(raw_ctrl_values)

    candidate: CandidateResult | None = None

    for ifn in ifns:
        raw_ifn_values: list[float] = get_custom_genotype_expression_values(row=row, pattern=ifn)
        log2_ifn_values: list[float] = log2_transform(raw_ifn_values)

        log2_fold_change: float = calculate_log2_fold_change(ctrl_values=log2_ctrl_values, ko_values=log2_ifn_values)
        p_value: float = calculate_raw_p_val(ctrl_values=log2_ctrl_values,ko_values=log2_ifn_values)

        if not is_significant(log2_fold_change=log2_fold_change, p_value=p_value):
            continue

        if log2_fold_change <= 0.0:
            continue

        # A candidate was already found, so this gene responds
        # significantly to more than one IFN.
        if candidate is not None:
            return None

        candidate = CandidateResult(
            gene_symbol=gene_symbol,
            ifn=ifn,
            p_value=p_value,
            log2_fold_change=log2_fold_change,
        )

    return candidate