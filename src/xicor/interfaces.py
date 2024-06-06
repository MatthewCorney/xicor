from typing import TypedDict


class XiResult(TypedDict):
    """
    Interface for the xi result
    """
    fr: list
    xi: float
    CU: float


class XiStatistic(TypedDict):
    """
    Interface for the Xi statistic
    """
    xi: float
    sd: float
    pval: float