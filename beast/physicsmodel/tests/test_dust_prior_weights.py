import numpy as np
import pytest
import astropy.units as u

from beast.physicsmodel.priormodel import PriorDustModel


def test_av_prior_weights():
    """
    Test the A(V) dust prior weights
    """

    av_vals = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])

    # Test A(V) priors

    models = [
        {"name": "flat"},
        {"name": "lognormal", "mean": 2.0, "sigma": 1.0},
        {
            "name": "two_lognormal",
            "mean1": 0.2,
            "mean2": 2.0,
            "sigma1": 1.0,
            "sigma2": 0.2,
            "N1_to_N2": 1.0 / 5.0,
        },
    ]
    # fmt: off
    expected_vals = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [0.0, 0.120985, 0.095149, 0.066168, 0.046282, 0.033133],
        [0.000000e+00, 4.969014e-01, 2.563834e-03, 3.719773e-04, 1.341047e-04, 5.743044e-05],
    ]
    # fmt: on

    for cmod, cvals in zip(models, expected_vals):
        av_prior = PriorDustModel(cmod)
        mname = cmod["name"]
        np.testing.assert_allclose(
            av_prior(av_vals),
            cvals,
            atol=1e-6,
            err_msg=f"A problem occurred while setting the A(V) priors to {mname}.",
        )

    # test error output when unsupported model is given
    with pytest.raises(NotImplementedError) as exc_av:
        PriorDustModel({"name": "exp"})
    assert str(exc_av.value) == "exp is not an allowed model"


def test_av_prior_weights_2d():
    """
    Test the A(V) dust prior weights for 2D priors
    """

    dist_vals = np.array([50.0, 59.0, 61.0, 70.0]) * 1e3
    av_vals = np.array([0.1, 0.9, 1.1, 2.0])

    # Test A(V) priors
    model = {
        "name": "step",
        "dist0": 60 * u.kpc,
        "amp1": 0.1,
        "damp2": 1.0,
        "lgsigma1": 0.05,
        "lgsigma2": 0.05,
    }

    # fmt: off
    expected_vals = [7.968878e+01, 0.000000e+00, 7.244435e+00, 6.544045e-31]
    # fmt: on

    av_prior = PriorDustModel(model)
    mname = model["name"]
    np.testing.assert_allclose(
        av_prior(av_vals, y=dist_vals),
        expected_vals,
        atol=1e-6,
        err_msg=f"A problem occurred while setting the 2D A(V) priors to {mname}.",
    )


def test_rv_prior_weights():
    """
    Test the R(V) dust prior weights
    """

    rv_vals = np.array([2.0, 3.0, 4.0, 5.0, 6.0])

    # Test R(V) priors

    models = [
        {"name": "flat"},
        {"name": "lognormal", "mean": 3.1, "sigma": 0.25},
        {
            "name": "two_lognormal",
            "mean1": 3.1,
            "mean2": 4.5,
            "sigma1": 0.1,
            "sigma2": 0.2,
            "N1_to_N2": 2.0 / 5.0,
        },
    ]
    # fmt: off
    expected_vals = [
        [1.0, 1.0, 1.0, 1.0, 1.0],
        [0.107331, 0.494654, 0.296693, 0.080183, 0.015239],
        [1.096692e-04, 3.915668e-01, 2.787845e-01, 2.737155e-01, 1.117566e-01],
    ]
    # fmt: on

    for cmod, cvals in zip(models, expected_vals):
        rv_prior = PriorDustModel(cmod)
        mname = cmod["name"]
        np.testing.assert_allclose(
            rv_prior(rv_vals),
            cvals,
            atol=1e-6,
            err_msg=f"A problem occurred while setting the R(V) priors to {mname}.",
        )


def test_rv_prior_weights_2d():
    """
    Test the R(V) dust prior weights for 2D priors
    """

    dist_vals = np.array([50.0, 59.0, 61.0, 70.0]) * 1e3
    rv_vals = np.array([3.1, 3.8, 4.6, 5.5])

    # Test A(V) priors
    model = {
        "name": "step",
        "dist0": 60 * u.kpc,
        "amp1": 3.1,
        "damp2": 1.5,
        "lgsigma1": 0.05,
        "lgsigma2": 0.05,
    }

    # fmt: off
    expected_vals = [2.570606e+00, 6.449322e-04, 1.732365e+00, 2.918564e-03]
    # fmt: on

    av_prior = PriorDustModel(model)
    mname = model["name"]
    np.testing.assert_allclose(
        av_prior(rv_vals, y=dist_vals),
        expected_vals,
        atol=1e-6,
        err_msg=f"A problem occurred while setting the 2D R(V) priors to {mname}.",
    )


def test_fA_prior_weights():
    """
    Test the f_A dust prior weights
    """

    fA_vals = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # Test f_A priors

    models = [
        {"name": "flat"},
        {"name": "lognormal", "mean": 0.8, "sigma": 0.1},
        {
            "name": "two_lognormal",
            "mean1": 0.1,
            "mean2": 0.8,
            "sigma1": 0.1,
            "sigma2": 0.2,
            "N1_to_N2": 2.0 / 5.0,
        },
    ]
    # fmt: off
    expected_vals = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [0.000000e+00, 9.205133e-42, 1.831235e-10, 7.916202e-02, 4.961907e+00, 4.115346e-01],
        [0.000000e+00, 8.506796e-10, 7.578321e-03, 1.092738e+00, 3.074674e+00, 1.650018e+00],
    ]
    # fmt: on

    for cmod, cvals in zip(models, expected_vals):
        fA_prior = PriorDustModel(cmod)
        mname = cmod["name"]
        np.testing.assert_allclose(
            fA_prior(fA_vals),
            cvals,
            atol=1e-6,
            err_msg=f"A problem occurred while setting the fA priors to {mname}.",
        )


def test_fA_prior_weights_2d():
    """
    Test the fA dust prior weights for 2D priors
    """

    dist_vals = np.array([50.0, 59.0, 61.0, 70.0]) * 1e3
    fA_vals = np.array([0.1, 0.5, 0.9, 0.3])

    # Test A(V) priors
    model = {
        "name": "step",
        "dist0": 60 * u.kpc,
        "amp1": 0.1,
        "damp2": 0.8,
        "lgsigma1": 0.05,
        "lgsigma2": 0.1,
    }

    # fmt: off
    expected_vals = [7.968878e+001, 8.158912e-224, 4.410584e+000, 2.728553e-026]
    # fmt: on

    av_prior = PriorDustModel(model)
    mname = model["name"]
    np.testing.assert_allclose(
        av_prior(fA_vals, y=dist_vals),
        expected_vals,
        atol=1e-6,
        err_msg=f"A problem occurred while setting the 2D f_A priors to {mname}.",
    )
