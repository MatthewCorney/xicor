import numpy as np
from typing import Literal, Union, TypedDict
from scipy.stats import norm
import logging

logging.basicConfig(format="[%(levelname)s] %(asctime)s %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


class XiResult(TypedDict):
    fr: list
    xi: float
    CU: float


class XiStatistic(TypedDict):
    xi: float
    sd: float
    pval: float


def calculate_xi_statistics(res: XiResult) -> XiStatistic:
    """

    :param res:
    :return:
    """
    fr = np.array(res['fr'])
    n = len(fr)

    # The following steps calculate the theoretical variance in the presence of ties:
    qfr = np.sort(fr)
    ind = np.arange(1, n + 1)
    ind2 = 2 * n - 2 * ind + 1
    ai = np.mean(ind2 * qfr * qfr) / n
    ci = np.mean(ind2 * qfr) / n
    cq = np.cumsum(qfr)
    m = (cq + (n - ind) * qfr) / n
    b = np.mean(m ** 2)
    v = (ai - 2 * b + ci ** 2) / (res['CU'] ** 2)

    # Return xi, standard deviation of xi, and P-value:
    xi = res['xi']
    sd_xi = np.sqrt(v / n)
    pval = 1 - norm.cdf(np.sqrt(n) * xi / np.sqrt(v))
    statistic = XiStatistic(xi=xi, sd=sd_xi, pval=pval)
    return statistic


def random_rank(vec: np.ndarray) -> np.ndarray:
    """
    From an array vec, sorts it in order breaking ties at random
    :param vec:
    :return:
    """
    np.random.seed(42)  # Seed here for consistency
    random_values = np.random.random(len(vec))
    combined = list(zip(vec, random_values))
    sorted_combined = sorted(combined)
    ranks = [sorted_combined.index(x) + 1 for x in combined]
    return np.array(ranks)


def calculate_xi(xvec: np.ndarray,
                 yvec: np.ndarray,
                 ) -> XiResult:
    """

    :param xvec:
    :param yvec:
    :return:
    """
    n = len(xvec)

    # PI is the rank vector for x, with ties broken at random
    PI = random_rank(xvec)
    # fr[i] is number of j s.t. y[j] <= y[i], divided by n.
    fr = np.array([np.sum(yvec <= yi) for yi in yvec]) / n

    # gr[i] is number of j s.t. y[j] >= y[i], divided by n.
    gr = np.array([np.sum(yvec >= yi) for yi in yvec]) / n

    # order of the x's, ties broken at random
    ord_ = np.argsort(PI)

    # Rearrange fr according to ord.
    fr = fr[ord_]

    # xi is calculated in the next three lines:
    A1 = np.sum(np.abs(fr[:-1] - fr[1:])) / (2 * n)
    CU = np.mean(gr * (1 - gr))[0]
    xi = 1 - A1 / CU
    xi_result = XiResult(xi=xi, fr=fr, CU=CU)
    return xi_result


def calculate_xi_permutation_test(yvec: np.ndarray,
                                  nperm: int,
                                  xi: float):
    """

    :param yvec:
    :param nperm:
    :param xi:
    :return:
    """
    n = len(yvec)

    # Initialize the permutation results array
    rp = np.zeros(nperm)

    # Perform permutations
    for i in range(nperm):
        np.random.seed(42 + i)
        x1 = np.random.uniform(0, 1, n)
        rp[i] = calculate_xi(x1, yvec)

    # Calculate the standard deviation and P-value based on permutation test
    sd_rp = np.sqrt(np.var(rp))
    pval = np.mean(rp > xi)[0]
    statistic = XiStatistic(xi=xi, sd=sd_rp, pval=pval)
    return statistic


def xicor(xvec: Union[np.ndarray, list],
          yvec: Union[np.ndarray, list],
          method: Literal['asymptotic', 'permutation'], nperm: int):
    """

    :param xvec:
    :param yvec:
    :param method:
    :param nperm:
    :return:
    """
    if isinstance(xvec, list):
        xvec = np.array(xvec)
    if isinstance(yvec, list):
        yvec = np.array(yvec)
    if method == "asymptotic":
        res = calculate_xi(xvec,
                           yvec,
                           )
        res = calculate_xi_statistics(res)
        return res
    elif method == 'permutation':
        res = calculate_xi(xvec,
                           yvec,
                           )
        res = calculate_xi_permutation_test(yvec=yvec, xi=res['xi'], nperm=nperm)
    else:
        logger.error('Method should be either "asymptotic" or "permutation"')
        raise Exception
    return res
