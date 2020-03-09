import numpy as np

from beast.physicsmodel.prior_weights_stars import compute_distance_prior_weights


def test_flat_distance_prior_weight():
    """
    Test the flat distance prior
    """
    dists = [10., 100., 1000.]
    dist_prior_model = {"name": "flat"}
    weight = compute_distance_prior_weights(dists, dist_prior_model)

    np.testing.assert_allclose(weight, np.full((len(weight)), 1.0))
