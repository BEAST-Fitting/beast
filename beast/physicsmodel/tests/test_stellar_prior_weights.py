import numpy as np

from beast.physicsmodel.priormodel import (
    PriorAgeModel,
    PriorMassModel,
    PriorMetallicityModel,
    PriorDistanceModel,
)


def test_age_prior_weights():
    """
    Test the age prior weights
    """

    log_age = np.array([6.0, 7.0, 8.0, 9.0, 10.0])

    models = [
        {"name": "flat"},
        {"name": "flat", "amp": 0.1},
        {"name": "flat_log"},
        {
            "name": "bins_histo",
            "x": [6.0, 7.0, 8.0, 9.0, 10.0],
            "values": [1.0, 2.0, 1.0, 5.0, 3.0],
        },
        {  # case when values has one entry less than x -> assumed that last x has 0 SFR
            "name": "bins_histo",
            "x": [6.0, 7.0, 8.0, 9.0, 10.0],
            "values": [1.0, 2.0, 1.0, 5.0],
        },
        {
            "name": "bins_interp",
            "x": [6.0, 7.0, 8.0, 9.0, 10.0],
            "values": [1.0, 2.0, 1.0, 5.0, 3.0],
        },
        {"name": "exponential", "tau": 0.1},
    ]
    # fmt: off
    expected_vals = [
        [1.0, 1.0, 1.0, 1.0, 1.0],
        [0.1, 0.1, 0.1, 0.1, 0.1],
        [9.00009e-01, 9.00009e-02, 9.00009e-03, 9.00009e-04, 9.00009e-05],
        [1.0, 2.0, 1.0, 5.0, 3.0],
        [1.0, 2.0, 1.0, 5.0, 0.0],
        [1.0, 2.0, 1.0, 5.0, 3.0],
        [9.900498e-01, 9.048374e-01, 3.678794e-01, 4.539993e-05, 3.720076e-44]
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


def test_mass_prior_weights():
    """
    Test the mass prior weights
    """

    mass = np.array([1, 2, 3, 4, 5])

    models = [
        {"name": "flat"},
        {"name": "salpeter"},
        {"name": "salpeter", "slope": 2.35},
        {"name": "salpeter", "slope": 2.00},
        {"name": "kroupa"},
        {"name": "kroupa", "alpha0": 0.3, "alpha1": 1.3, "alpha2": 2.3, "alpha3": 2.3},
        {"name": "kroupa", "alpha0": 0.5, "alpha1": 1.0, "alpha2": 2.0, "alpha3": 2.5},
    ]
    # fmt: off
    expected_vals = [
        [1.0, 1.0, 1.0, 1.0, 1.0],
        [4.02338441, 0.58842044, 0.21633931, 0.10825509, 0.06360075],
        [4.02338441, 0.58842044, 0.21633931, 0.10825509, 0.06360075],
        [3.666667, 0.733333, 0.314286, 0.174603, 0.111111],
        [3.97740709, 0.60861986, 0.22874078, 0.11618704, 0.06904523],
        [3.97740709, 0.60861986, 0.22874078, 0.11618704, 0.06904523],
        [4.036514, 0.601346, 0.20694 , 0.098998, 0.056201],
    ]
    # fmt: on

    for cmod, cvals in zip(models, expected_vals):
        mass_prior = PriorMassModel(cmod)
        mname = cmod["name"]
        np.testing.assert_allclose(
            mass_prior(mass),
            cvals,
            atol=1e-6,
            err_msg=f"A problem occurred while setting the mass priors to {mname}.",
        )


def test_met_prior_weights():
    """
    Test the metallicity prior weights
    """

    met = np.array([0.01, 0.1, 1.0])

    models = [{"name": "flat"}]
    # fmt: off
    expected_vals = [
        [1.0, 1.0, 1.0]
    ]
    # fmt: on

    for cmod, cvals in zip(models, expected_vals):
        met_prior = PriorMetallicityModel(cmod)
        mname = cmod["name"]
        np.testing.assert_allclose(
            met_prior(met),
            cvals,
            atol=1e-6,
            err_msg=f"A problem occurred while setting the metallicity priors to {mname}.",
        )


def test_dist_prior_weights():
    """
    Test the distance prior weights
    """

    dist = np.array([10.0, 100.0, 1000.0])

    models = [{"name": "flat"}]
    # fmt: off
    expected_vals = [
        [1.0, 1.0, 1.0]
    ]
    # fmt: on

    for cmod, cvals in zip(models, expected_vals):
        dist_prior = PriorDistanceModel(cmod)
        mname = cmod["name"]
        np.testing.assert_allclose(
            dist_prior(dist),
            cvals,
            atol=1e-6,
            err_msg=f"A problem occurred while setting the distance priors to {mname}.",
        )
