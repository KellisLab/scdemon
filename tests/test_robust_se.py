
import pytest
import numpy as np
import scdemon as scd


@pytest.fixture(params=[(1000, 30), (500, 200), (200, 10)])
def random_array_and_row(request):
    shape = request.param
    np.random.seed(42)
    array = np.random.random(shape).astype("f8")
    col_idx = np.random.randint(0, shape[1])
    return col_idx, array

def test_robust_se_X(random_array_and_row):
    from statsmodels.api import OLS
    from scdemon.robust_se import _robust_se_X
    i, X = random_array_and_row
    scd_tv = _robust_se_X(i, X)
    np.testing.assert_array_equal(scd_tv.shape, (X.shape[1],))
    ols_tv = np.zeros_like(scd_tv)
    for j in range(X.shape[1]):
        ols_tv[j] = OLS(X[:, j], X[:, i]).fit(cov_type="HC0").tvalues[0]
    assert np.argmax(scd_tv) == np.argmax(ols_tv)
    assert np.argmax(scd_tv) == i
    scd_tv[i] = 0 ## Numerical instability, set to 0 after check
    ols_tv[i] = 0
    np.testing.assert_allclose(scd_tv, ols_tv)
