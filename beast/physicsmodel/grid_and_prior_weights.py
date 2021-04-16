"""
Grid and Prior Weights
======================
The use of a non-uniformly spaced grid complicates the marginalization
step as the trick of summation instead of integration is used.  But this
trick only works when the grid is uniformaly spaced in all dimensions.

If the grid is not uniformally spaced, weights can be used to correct
for the non-uniform spacing.

Basically, we want the maginalization using these grid weights to provide
flat priors on all the fit parameters.  Non-flat priors will be implemented
with prior weights.
"""
import numpy as np

from beast.physicsmodel.grid_weights_stars import compute_distance_grid_weights
from beast.physicsmodel.grid_weights_stars import compute_age_grid_weights
from beast.physicsmodel.grid_weights_stars import compute_mass_grid_weights
from beast.physicsmodel.grid_weights_stars import compute_metallicity_grid_weights

from beast.physicsmodel.priormodel import (
    PriorAgeModel,
    PriorMassModel,
    PriorMetallicityModel,
    PriorDistanceModel,
)

__all__ = [
    "compute_age_mass_metallicity_weights",
    "compute_distance_age_mass_metallicity_weights",
]


def compute_distance_age_mass_metallicity_weights(
    _tgrid,
    distance_prior_model={"name": "flat"},
    age_prior_model={"name": "flat"},
    mass_prior_model={"name": "kroupa"},
    met_prior_model={"name": "flat"},
):
    """
    Computes the distance and age-mass-metallicity grid
    and prior weights on the BEAST model spectra grid

    Parameters
    ----------
    _tgrid : BEAST model spectra grid.

    distance_prior_model,
    age_prior_model,
    mass_prior_model,
    met_prior_model: dict
        dict including prior model name and parameters

    Returns
    -------
    Grid and prior weight columns updated by multiplying by the
    the distance and age-mass-metallicity weight.
    """

    # get the unique distances
    uniq_dists = np.unique(_tgrid["distance"])

    # setup the vector to hold the distance weight vectors
    n_dist = len(uniq_dists)
    total_dist_grid_weight = np.zeros((n_dist))
    total_dist_prior_weight = np.zeros((n_dist))
    total_dist_weight = np.zeros((n_dist))

    for dz, dist_val in enumerate(uniq_dists):
        print("computing the distance plus weights for dist = ", dist_val)
        (dindxs,) = np.where(_tgrid["distance"] == dist_val)
        compute_age_mass_metallicity_weights(
            _tgrid,
            dindxs,
            age_prior_model=age_prior_model,
            mass_prior_model=mass_prior_model,
            met_prior_model=met_prior_model,
        )
        total_dist_grid_weight[dz] = np.sum(_tgrid[dindxs]["grid_weight"])
        total_dist_prior_weight[dz] = np.sum(_tgrid[dindxs]["prior_weight"])
        total_dist_weight[dz] = np.sum(_tgrid[dindxs]["weight"])

    # ensure that the distance prior is uniform
    if n_dist > 1:
        # get the distance weights
        dist_grid_weights = compute_distance_grid_weights(uniq_dists)
        dist_grid_weights /= np.sum(dist_grid_weights)
        dist_prior = PriorDistanceModel(distance_prior_model)
        dist_prior_weights = dist_prior(uniq_dists)
        dist_prior_weights /= np.sum(dist_prior_weights)
        dist_weights = dist_grid_weights * dist_prior_weights

        # correct for any non-unformity in the number size of the
        # age-mass grids between metallicity points
        total_dist_grid_weight /= np.sum(total_dist_grid_weight)
        total_dist_prior_weight /= np.sum(total_dist_prior_weight)
        total_dist_weight /= np.sum(total_dist_weight)

        for i, dist_val in enumerate(uniq_dists):
            # get the grid for this distance
            (dindxs,) = np.where(_tgrid["distance"] == dist_val)
            _tgrid[dindxs]["grid_weight"] *= (
                dist_grid_weights[i] * total_dist_grid_weight[i]
            )
            _tgrid[dindxs]["prior_weight"] *= (
                dist_prior_weights[i] * total_dist_prior_weight[i]
            )
            _tgrid[dindxs]["weight"] *= dist_weights[i] * total_dist_weight[i]


def compute_age_mass_metallicity_weights(
    _tgrid,
    indxs,
    age_prior_model={"name": "flat"},
    mass_prior_model={"name": "kroupa"},
    met_prior_model={"name": "flat"},
    **kwargs
):
    """
    Computes the age-mass-metallicity grid and prior weights
    on the BEAST model spectra grid
    Grid and prior weight columns updated by multiplying by the
    age-mass-metallicity weight.

    Parameters
    ----------
    _tgrid : SpectralGrid
        BEAST models spectral grid
    age_prior_model : dict
        dict including prior model name and parameters
    mass_prior_model : dict
        dict including prior model name and parameters
    met_prior_model : dict
        dict including prior model name and parameters
    """

    # get the unique metallicities
    uniq_Zs = np.unique(_tgrid[indxs]["Z"])

    # setup the vector to hold the z weight vector
    total_z_grid_weight = np.zeros(len(uniq_Zs))
    total_z_prior_weight = np.zeros(len(uniq_Zs))
    total_z_weight = np.zeros(len(uniq_Zs))

    for az, z_val in enumerate(uniq_Zs):
        print("computing the age-mass-metallicity grid weight for Z = ", z_val)

        # get the grid for a single metallicity
        (zindxs,) = np.where(_tgrid[indxs]["Z"] == z_val)

        # get the unique ages for this metallicity
        zindxs = indxs[zindxs]
        uniq_ages = np.unique(_tgrid[zindxs]["logA"])

        # compute the age weights
        age_grid_weights = compute_age_grid_weights(uniq_ages)
        age_prior = PriorAgeModel(age_prior_model)
        age_prior_weights = age_prior(uniq_ages)

        for ak, age_val in enumerate(uniq_ages):
            # get the grid for a single age
            (aindxs,) = np.where(
                (_tgrid[indxs]["logA"] == age_val) & (_tgrid[indxs]["Z"] == z_val)
            )
            aindxs = indxs[aindxs]
            _tgrid_single_age = _tgrid[aindxs]

            # compute the mass weights
            if len(aindxs) > 1:
                cur_masses = _tgrid_single_age["M_ini"]
                mass_grid_weights = compute_mass_grid_weights(cur_masses)
                mass_prior = PriorMassModel(mass_prior_model)
                mass_prior_weights = mass_prior(cur_masses)
            else:
                # must be a single mass for this age,z combination
                # set mass weight to zero to remove this point from the grid
                mass_grid_weights = np.zeros(1)
                mass_prior_weights = np.zeros(1)

            # apply both the mass and age weights
            for i, k in enumerate(aindxs):
                comb_grid_weights = mass_grid_weights[i] * age_grid_weights[ak]
                comb_prior_weights = mass_prior_weights[i] * age_prior_weights[ak]
                _tgrid[k]["grid_weight"] *= comb_grid_weights
                _tgrid[k]["prior_weight"] *= comb_prior_weights
                _tgrid[k]["weight"] *= comb_grid_weights * comb_prior_weights

        # compute the current total weight at each metallicity
        total_z_grid_weight[az] = np.sum(_tgrid[zindxs]["grid_weight"])
        total_z_prior_weight[az] = np.sum(_tgrid[zindxs]["prior_weight"])
        total_z_weight[az] = np.sum(_tgrid[zindxs]["weight"])

    # ensure that the metallicity prior is uniform
    if len(uniq_Zs) > 1:
        # get the metallicity weights
        met_grid_weights = compute_metallicity_grid_weights(uniq_Zs)
        met_grid_weights /= np.sum(met_grid_weights)
        met_prior = PriorMetallicityModel(met_prior_model)
        met_prior_weights = met_prior(uniq_Zs)
        met_prior_weights /= np.sum(met_prior_weights)
        met_weights = met_grid_weights * met_prior_weights

        # correct for any non-unformity in the number size of the
        # age-mass grids between metallicity points
        total_z_grid_weight /= np.sum(total_z_grid_weight)
        total_z_prior_weight /= np.sum(total_z_prior_weight)
        total_z_weight /= np.sum(total_z_weight)

        for i, z_val in enumerate(uniq_Zs):
            # get the grid for this metallicity
            (zindxs,) = np.where(_tgrid[indxs]["Z"] == z_val)
            zindxs = indxs[zindxs]
            _tgrid[zindxs]["grid_weight"] *= (
                met_grid_weights[i] * total_z_grid_weight[i]
            )
            _tgrid[zindxs]["prior_weight"] *= (
                met_prior_weights[i] * total_z_prior_weight[i]
            )
            _tgrid[zindxs]["weight"] *= met_weights[i] * total_z_weight[i]
