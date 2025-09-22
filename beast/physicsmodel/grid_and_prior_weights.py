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

from beast.physicsmodel.grid_weights import compute_grid_weights

from beast.physicsmodel.priormodel import (
    PriorAgeModel,
    PriorMassModel,
    PriorMetallicityModel,
    PriorDistanceModel,
    PriorDustModel,
)

__all__ = [
    "compute_age_mass_metallicity_weights",
    "compute_distance_age_mass_metallicity_weights",
    "compute_av_rv_fA_prior_weights",
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

    distance_prior_model, age_prior_model, mass_prior_model, met_prior_model: dict
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

    if n_dist > 1:
        # get the distance weights
        dist_grid_weights = compute_grid_weights(uniq_dists)
        dist_grid_weights /= np.sum(dist_grid_weights)
        dist_prior = PriorDistanceModel(distance_prior_model)
        dist_prior_weights = dist_prior(uniq_dists)
        dist_prior_weights /= np.sum(dist_prior_weights)
        dist_weights = dist_grid_weights * dist_prior_weights

        # correct for any non-uniformity in the number size of the
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


def compute_age_mass_weights(
    _tgrid,
    z_val,
    indxs,
    age_prior_model={"name": "flat"},
    mass_prior_model={"name": "kroupa"},
):
    """
    Compute the age and meteallicity prior.  Handles two cases.
    The age on a grid and mass non-uniform and the mass on a grid and the age non-uniform.

    _tgrid : SpectralGrid
        BEAST models spectral grid
    z_val : float
        value of metallicity to compute priors
    age_prior_model : dict
        dict including prior model name and parameters
    mass_prior_model : dict
        dict including prior model name and parameters
    """

    # get the grid for a single metallicity
    (zindxs,) = np.where(_tgrid[indxs]["Z"] == z_val)

    # determine which of the two cases based on the variable on a uniform grid will
    # have many fewer unique elements
    ages = np.unique(_tgrid[indxs]["logA"])
    masses = np.unique(_tgrid[indxs]["M_ini"])
    if len(ages) > len(masses):
        pname1 = "M_ini"
        pname1_logval = False
        prior1 = mass_prior_model
        prior1mod = PriorMassModel
        pname2 = "logA"
        pname2_logval = True
        prior2 = age_prior_model
        prior2mod = PriorAgeModel
    else:
        pname1 = "logA"
        pname1_logval = True
        prior1 = age_prior_model
        prior1mod = PriorAgeModel
        pname2 = "M_ini"
        pname2_logval = False
        prior2 = mass_prior_model
        prior2mod = PriorMassModel

    # get the unique values of pname1 for this metallicity
    zindxs = indxs[zindxs]
    uniq_pname1 = np.unique(_tgrid[zindxs][pname1])

    pname1_grid_weights = compute_grid_weights(uniq_pname1, log=pname1_logval)
    if isinstance(prior1, dict):
        pname1_prior = prior1mod(prior1)
    else:
        pname1_prior = prior1
    pname1_prior_weights = pname1_prior(uniq_pname1)

    for ak, pname1_val in enumerate(uniq_pname1):
        # get the grid for a single age
        (aindxs,) = np.where(
            (_tgrid[indxs][pname1] == pname1_val) & (_tgrid[indxs]["Z"] == z_val)
        )
        aindxs = indxs[aindxs]
        _tgrid_single_pname1 = _tgrid[aindxs]

        # compute the mass weights
        if len(aindxs) > 1:
            if isinstance(prior2, dict):
                pname2_prior = PriorMassModel(prior2)
            else:
                pname2_prior = prior2

            # deal with repeat masses or ages - happens for MegaBEAST
            # and have discovered this can happen even for a standard BEAST run
            # as sometimes two masses in an isochrone are exactly the same
            #   new code for MegaBEAST is more correct as then the grid weight
            #   will be correctly set for any repeated masses
            cur_pname2 = np.unique(_tgrid_single_pname1[pname2])
            n_pname2 = len(_tgrid_single_pname1[pname2])
            if len(cur_pname2) < n_pname2:
                upname2_grid_weights = compute_grid_weights(
                    cur_pname2, log=pname2_logval
                )
                upname2_prior_weights = pname2_prior(cur_pname2)
                pname2_grid_weights = np.zeros(n_pname2, dtype=float)
                pname2_prior_weights = np.zeros(n_pname2, dtype=float)
                for k, cpname2 in enumerate(cur_pname2):
                    gvals = _tgrid_single_pname1[pname2] == cpname2
                    pname2_grid_weights[gvals] = upname2_grid_weights[k]
                    pname2_prior_weights[gvals] = upname2_prior_weights[k]
            else:
                cur_pname2 = _tgrid_single_pname1[pname2]
                pname2_grid_weights = compute_grid_weights(
                    cur_pname2, log=pname2_logval
                )
                pname2_prior_weights = pname2_prior(cur_pname2)

        else:
            # must be a single mass for this age,z combination
            # set mass weight to zero to remove this point from the grid
            pname2_grid_weights = np.zeros(1)
            pname2_prior_weights = np.zeros(1)

        # apply both the mass and age weights
        for i, k in enumerate(aindxs):
            comb_grid_weights = pname2_grid_weights[i] * pname1_grid_weights[ak]
            comb_prior_weights = pname2_prior_weights[i] * pname1_prior_weights[ak]
            _tgrid[k]["grid_weight"] *= comb_grid_weights
            _tgrid[k]["prior_weight"] *= comb_prior_weights
            _tgrid[k]["weight"] *= comb_grid_weights * comb_prior_weights


def compute_age_mass_metallicity_weights(
    _tgrid,
    indxs,
    age_prior_model={"name": "flat"},
    mass_prior_model={"name": "kroupa"},
    met_prior_model={"name": "flat"},
    **kwargs,
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

        # compute the age and mass weights
        compute_age_mass_weights(
            _tgrid, z_val, indxs, age_prior_model, mass_prior_model
        )

        # compute the current total weight at each metallicity
        (zindxs,) = np.where(_tgrid[indxs]["Z"] == z_val)
        total_z_grid_weight[az] = np.sum(_tgrid[zindxs]["grid_weight"])
        total_z_prior_weight[az] = np.sum(_tgrid[zindxs]["prior_weight"])
        total_z_weight[az] = np.sum(_tgrid[zindxs]["weight"])

    # ensure that the metallicity prior is uniform
    if len(uniq_Zs) > 1:
        # get the metallicity weights
        met_grid_weights = compute_grid_weights(uniq_Zs)
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


def compute_av_rv_fA_prior_weights(
    Av,
    Rv,
    f_A,
    dists,
    av_prior_model={"name": "flat"},
    rv_prior_model={"name": "flat"},
    fA_prior_model={"name": "flat"},
):
    """
    Computes the av, rv, f_A grid and prior weights
    on the BEAST model spectra grid
    Grid and prior weight columns updated by multiplying by the
    existing weights

    Parameters
    ----------
    Av : vector
        A(V) values
    Rv : vector
        R(V) values
    f_A : vector
        f_A values
    dists : vector
        distance values
    av_prior_model : dict
        dict including prior model name and parameters
    rv_prior_model : dict
        dict including prior model name and parameters
    fA_prior_model :dict
        dict including prior model name and parameters
    """
    av_prior = PriorDustModel(av_prior_model)
    rv_prior = PriorDustModel(rv_prior_model)
    fA_prior = PriorDustModel(fA_prior_model)
    if av_prior_model["name"] == "step":
        av_weights = av_prior(np.full((len(dists)), Av), y=dists)
    else:
        av_weights = av_prior(Av)
    if rv_prior_model["name"] == "step":
        rv_weights = rv_prior(np.full((len(dists)), Rv), y=dists)
    else:
        rv_weights = rv_prior(Rv)
    if fA_prior_model["name"] == "step":
        f_A_weights = fA_prior(np.full((len(dists)), f_A), y=dists)
    else:
        f_A_weights = fA_prior(f_A)

    dust_prior = av_weights * rv_weights * f_A_weights

    # normalize to control for numerical issues
    # dust_prior /= np.max(dust_prior)

    return dust_prior
