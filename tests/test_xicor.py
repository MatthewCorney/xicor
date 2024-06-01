import pytest
from src.xicor import xicor

anscombes_quartet = [([10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5],
                      [8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68]),
                     ([10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5],
                      [9.14, 8.14, 8.74, 8.77, 9.26, 8.1, 6.13, 3.1, 9.13, 7.26, 4.74]),
                     ([0, 0, 0, 0, 1, 1, 2, 3],
                      [0, 1, 2, 3, 4, 5, 6, 7]),
                     ([8, 8, 8, 8, 8, 8, 8, 19, 8, 8, 8],
                      [6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.5, 5.56, 7.91, 6.89]),
                     ]
expected_results_asym = [
    0.275,
    0.6,
    0.38,
    -0.25]


# Define pytest functions
@pytest.mark.parametrize("data, expected_results_asym", zip(anscombes_quartet, expected_results_asym))
def test_xicor_asymptotic(data, expected_results_asym):
    xvec, yvec = data
    result = xicor(xvec=xvec, yvec=yvec, method='asymptotic', nperm=10000)
    assert result['xi'] == pytest.approx(expected_results_asym, rel=1e-2), f"Failed for data {data}"


expected_results_perm = [
    0.275,
    0.6,
    0.38,
    -0.25]


@pytest.mark.parametrize("data, expected_results_perm", zip(anscombes_quartet, expected_results_perm))
def test_xicor_permutation(data, expected_results_perm):
    xvec, yvec = data
    result = xicor(xvec=xvec, yvec=yvec, method='permutation', nperm=100)
    assert result['xi'] == pytest.approx(expected_results_perm, rel=1e-2), f"Failed for data {data}"
