import numpy as np

from beast.physicsmodel.prior_weights_stars import (
    compute_distance_prior_weights,
    compute_age_prior_weights,
    compute_mass_prior_weights,
    compute_metallicity_prior_weights,
    imf_kroupa,
)


def test_flat_age_prior_weights():
    """
    Test for flat age prior
    """
    log_age = np.array([6.0, 7.0, 8.0, 9.0, 10.0])
    log_age_prior_model = {"name": "flat"}
    log_age_prior = compute_age_prior_weights(log_age, log_age_prior_model)
    expected_log_age_prior = [1, 1, 1, 1, 1]
    np.testing.assert_allclose(
        log_age_prior, expected_log_age_prior, err_msg=("Flat age prior error")
    )


def test_flat_log_age_prior_weights():
    """
    Test for flat log age prior
    """
    log_age = np.array([6.0, 7.0, 8.0, 9.0, 10.0])
    log_age_prior_model = {"name": "flat_log"}
    log_age_prior = compute_age_prior_weights(log_age, log_age_prior_model)
    expected_log_age_prior = [
        4.500045e00,
        4.500045e-01,
        4.500045e-02,
        4.500045e-03,
        4.500045e-04,
    ]
    np.testing.assert_allclose(
        log_age_prior, expected_log_age_prior, err_msg=("Flat log, log age prior error")
    )


def test_bins_histo_age_prior_weights():
    """
    Test for bin histogram age prior
    """
    log_age = np.array([7.0, 8.0, 9.0])
    log_age_prior_model = {
        "name": "bins_histo",
        "logages": [6.0, 7.0, 8.0, 9.0, 10.0],
        "values": [1.0, 2.0, 1.0, 5.0, 3.0],
    }
    log_age_prior = compute_age_prior_weights(log_age, log_age_prior_model)
    expected_log_age_prior = [0.75, 0.375, 1.875]
    np.testing.assert_allclose(
        log_age_prior,
        expected_log_age_prior,
        err_msg=("Bin histogram log age prior error"),
    )


def test_bins_interp_age_prior_weights():
    """
    Test for bin interpolation age prior
    """
    log_age = np.array([6.0, 7.0, 8.0, 9.0, 10.0])
    log_age_prior_model = {
        "name": "bins_interp",
        "logages": [6.0, 7.0, 8.0, 9.0, 10.0],
        "values": [1.0, 2.0, 1.0, 5.0, 3.0],
    }
    log_age_prior = compute_age_prior_weights(log_age, log_age_prior_model)
    expected_log_age_prior = [0.41666667, 0.83333333, 0.41666667, 2.08333333, 1.25]
    np.testing.assert_allclose(
        log_age_prior,
        expected_log_age_prior,
        err_msg=("Bin histogram log age prior error"),
    )


def test_exp_age_prior_weights():
    """
    Test for exponential age prior with a tau = 0.1
    """
    log_age = np.array([6.0, 7.0, 8.0, 9.0, 10.0])
    log_age_prior_model = {"name": "exp", "tau": 0.1}
    log_age_prior = compute_age_prior_weights(log_age, log_age_prior_model)
    expected_log_age_prior = [
        2.18765367e00,
        1.99936491e00,
        8.12881110e-01,
        1.00317499e-04,
        8.22002849e-44,
    ]
    np.testing.assert_allclose(
        log_age_prior,
        expected_log_age_prior,
        err_msg=("Exponential log age prior error"),
    )


def test_imf_kroupa():
    """
    Test for creating kroupa IMF
    """
    mass = np.array([0.1, 1, 2, 3, 4, 50])
    imf = imf_kroupa(mass)
    expected_imf = [
        3.99052463e01,
        1.00000000e00,
        2.03063099e-01,
        7.99136770e-02,
        4.12346222e-02,
        1.23699798e-04,
    ]
    np.testing.assert_allclose(
        imf, expected_imf, err_msg=("Kroupa IMF calculation error")
    )


def test_kroupa_mass_prior_weight():
    """
    Test the kroupa mass prior
    """
    mass = np.array([1, 2, 3, 4, 5])
    mass_prior_model = {"name": "kroupa"}
    weights = compute_mass_prior_weights(mass, mass_prior_model)
    expected_weights = [3.97740709, 0.60861986, 0.22874078, 0.11618704, 0.06904523]
    np.testing.assert_allclose(
        weights,
        expected_weights,
        err_msg=("Stellar mass prior weight error (kroupa IMF)"),
    )


def test_salpeter_mass_prior_weight():
    """
    Test the salpeter mass prior
    """
    mass = np.array([1, 2, 3, 4, 5])
    mass_prior_model = {"name": "salpeter"}
    weights = compute_mass_prior_weights(mass, mass_prior_model)
    expected_weights = [4.02338441, 0.58842044, 0.21633931, 0.10825509, 0.06360075]
    np.testing.assert_allclose(
        weights,
        expected_weights,
        err_msg=("Stellar mass prior weight error (salpeter IMF)"),
    )


def test_flat_mass_prior_weight():
    """
    Test the flat mass prior
    """
    mass = np.array([1, 2, 3, 4, 5])
    mass_prior_model = {"name": "flat"}
    weights = compute_mass_prior_weights(mass, mass_prior_model)
    np.testing.assert_allclose(
        weights,
        np.full((len(weights)), 1.0),
        err_msg=("Stellar mass prior weight error (flat IMF)"),
    )


def test_flat_metallicity_prior_weight():
    """
    Test the flat metallicity prior
    """
    z = [10.0, 100.0, 1000.0]
    z_prior_model = {"name": "flat"}
    weights = compute_metallicity_prior_weights(z, z_prior_model)
    np.testing.assert_allclose(
        weights,
        np.full((len(weights)), 1.0),
        err_msg=("Stellar flat metallicity prior weight error"),
    )


def test_flat_distance_prior_weight():
    """
    Test the flat distance prior
    """
    dists = [10.0, 100.0, 1000.0]
    dist_prior_model = {"name": "flat"}
    weights = compute_distance_prior_weights(dists, dist_prior_model)
    np.testing.assert_allclose(
        weights,
        np.full((len(weights)), 1.0),
        err_msg=("Stellar flat distance prior weight error"),
    )
