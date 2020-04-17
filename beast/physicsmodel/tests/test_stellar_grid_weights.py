import numpy as np

from beast.physicsmodel.grid_weights_stars import (
    compute_distance_grid_weights,
    compute_age_grid_weights,
    compute_mass_grid_weights,
    compute_metallicity_grid_weights,
)


def test_flat_distance_grid_weight():
    """
    Test the flat distance grid weights
    """
    dists = [10.0, 100.0, 1000.0]
    expected_weights = [0.18181818, 1.0, 1.81818182]

    weight = compute_distance_grid_weights(dists)

    np.testing.assert_allclose(weight, expected_weights)


def test_mass_grid_weights():
    """
    Test mass grid weights
    """
    masses = np.array([1, 2, 5, 7, 10, 30])
    weights = compute_mass_grid_weights(masses)
    expected_weights = [
        0.15189873,
        0.30379747,
        0.37974684,
        0.37974684,
        1.74683544,
        3.03797468,
    ]
    np.testing.assert_allclose(
        weights, expected_weights, err_msg=("Stellar grid mass weights error")
    )


def test_metallicity_grid_weights():
    """
    Test metallicities grid weights
    """
    metallicities = np.array([0.03, 0.019, 0.008, 0.004])
    weights = compute_metallicity_grid_weights(metallicities)
    expected_weights = [1.31343284, 1.31343284, 0.89552239, 0.47761194]
    np.testing.assert_allclose(
        weights, expected_weights, err_msg=("Stellar grid metallicity weights error")
    )
