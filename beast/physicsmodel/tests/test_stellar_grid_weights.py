import numpy as np

from beast.physicsmodel.grid_weights_stars import compute_distance_grid_weights


def test_flat_distance_grid_weight():
    """
    Test the flat distance grid weights
    """
    dists = [10., 100., 1000.]
    expected_weights = [0.18181818, 1.0, 1.81818182]

    weight = compute_distance_grid_weights(dists)

    np.testing.assert_allclose(weight, expected_weights)
