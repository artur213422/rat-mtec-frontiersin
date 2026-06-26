import math
from statistics import mean
from typing import Sequence

from scipy.stats import ttest_ind


def calculate_raw_p_val(ctrl_values: list[float], ko_values: list[float]) -> float:
    ctrl: list[float] = [
        float(value)
        for value in ctrl_values
        if math.isfinite(value)
    ]

    ko: list[float] = [
        float(value)
        for value in ko_values
        if math.isfinite(value)
    ]

    if len(ctrl) < 2:
        raise ValueError(f"Need at least 2 control values, got {len(ctrl)}.")

    if len(ko) < 2:
        raise ValueError(f"Need at least 2 KO values, got {len(ko)}.")

    result = ttest_ind(
        ctrl,
        ko,
        equal_var=False,   # Welch t-test
        alternative="two-sided",
    )

    p_value = float(result.pvalue)

    if math.isnan(p_value):
        ctrl_mean = sum(ctrl) / len(ctrl)
        ko_mean = sum(ko) / len(ko)

        if ctrl_mean == ko_mean:
            return 1.0

        return 0.0

    return p_value


def calculate_log2_fold_change(ctrl_values: list[float], ko_values: list[float]) -> float:
    if not ctrl_values:
        raise ValueError("Control values are empty.")

    if not ko_values:
        raise ValueError("KO values are empty.")

    return mean(ko_values) - mean(ctrl_values)