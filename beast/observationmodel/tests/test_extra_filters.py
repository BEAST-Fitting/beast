import numpy as np
import pytest

from beast.observationmodel.extra_filters import make_integration_filter, make_top_hat_filter


@pytest.mark.parametrize(
    "lambda_start,lambda_finish,d_lambda",
    [(90., 913., 1.), (1000, 3000, 100)],
)
def test_extra_filters(lambda_start, lambda_finish, d_lambda):

    # create example integration filter
    f_int = make_integration_filter(
        lambda_start, lambda_finish, d_lambda, "QION", observatory="Pseudo", instrument="Fake"
    )

    # test bandwidth and name
    np.testing.assert_allclose(f_int.bandwidth, lambda_finish - lambda_start - d_lambda)
    if not f_int.name == "Pseudo_Fake_QION":
        raise AssertionError()

    # create example top hat filters
    f_top = make_top_hat_filter(
        lambda_start, lambda_finish, d_lambda, "TOP", observatory="Pseudo", instrument="Fake"
    )

    np.testing.assert_allclose(f_top.bandwidth, lambda_finish - lambda_start - d_lambda)
    if not f_top.name == "Pseudo_Fake_TOP":
        raise AssertionError()
