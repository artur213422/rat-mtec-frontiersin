import math

DEFAULT_MIN_ABS_LOG2_FOLD_CHANGE: float = 1.0
DEFAULT_MAX_P_VALUE: float = 0.05

def is_significant(log2_fold_change: float,
                   p_value: float,
                   min_abs_log2_fold_change: float = DEFAULT_MIN_ABS_LOG2_FOLD_CHANGE,
                   max_p_value: float = DEFAULT_MAX_P_VALUE) -> bool:

    if not (math.isfinite(log2_fold_change) and math.isfinite(p_value)):
        return False

    if not 0.0 <= p_value <= 1.0:
        raise ValueError(f"p_value must be between 0 and 1, received {p_value}.")

    if min_abs_log2_fold_change < 0.0:
        raise ValueError("min_abs_log2_fold_change cannot be negative.")

    return abs(log2_fold_change) >= min_abs_log2_fold_change and p_value <= max_p_value


