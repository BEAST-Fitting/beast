import numpy as np

from beast.physicsmodel.prior_weights_stars import *


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
        imf, expected_imf, err_msg=("\nKroupa IMF calculation error\n")
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
        err_msg=("\nStellar mass prior weight error (kroupa IMF)\n"),
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
        err_msg=("\nStellar mass prior weight error (salpeter IMF)\n"),
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
        err_msg=("\nStellar mass prior weight error (flat IMF)\n"),
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
        err_msg=("\nStellar flat metallicity prior weight error\n"),
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
        err_msg=("\nStellar flat distance prior weight error\n"),
    )
