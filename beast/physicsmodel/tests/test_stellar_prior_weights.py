import numpy as np

from beast.physicsmodel.priormodel import PriorAgeModel


def test_age_prior_weights():
    """
    Test the age prior weights
    """

    log_age = np.array([6.0, 7.0, 8.0, 9.0, 10.0])

    # Test A(V) priors
    #  (with random input values)

    models = [
        {"name": "flat"},
        {"name": "flat", "amp": 0.1},
        {
            "name": "bins_histo",
            "x": [6.0, 7.0, 8.0, 9.0, 10.0],
            "values": [1.0, 2.0, 1.0, 5.0, 3.0],
        },
        {
            "name": "bins_interp",
            "x": [6.0, 7.0, 8.0, 9.0, 10.0],
            "values": [1.0, 2.0, 1.0, 5.0, 3.0],
        },
        {"name": "exponential", "a": 0.1},
    ]
    # fmt: off
    expected_vals = [
        [1.0, 1.0, 1.0, 1.0, 1.0],
        [0.1, 0.1, 0.1, 0.1, 0.1],
        [1.0, 2.0, 1.0, 5.0, 3.0],
        [1.0, 2.0, 1.0, 5.0, 3.0],
        [
            9.900498e-01,
            9.048374e-01,
            3.678794e-01,
            4.539993e-05,
            3.720076e-44,
        ]
    ]
    # fmt: on

    for cmod, cvals in zip(models, expected_vals):
        age_prior = PriorAgeModel(cmod)
        mname = cmod["name"]
        np.testing.assert_allclose(
            age_prior(log_age),
            cvals,
            atol=1e-6,
            err_msg=f"A problem occurred while setting the age priors to {mname}.",
        )
