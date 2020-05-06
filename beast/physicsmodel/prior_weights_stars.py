"""
Prior Weights
=============
The priors on age, mass, and metallicty are computed as weights to use
in the posterior calculations.
"""
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

from beast.physicsmodel.grid_weights_stars import (
    compute_bin_boundaries,
    compute_age_grid_weights,
)

__all__ = [
    "compute_distance_prior_weights",
    "compute_age_prior_weights",
    "compute_mass_prior_weights",
    "compute_metallicity_prior_weights",
    "imf_kroupa",
]


def compute_age_prior_weights(logages, age_prior_model):
    """
    Computes the age prior for the specified model

    Parameters
    ----------
    logages : numpy vector
       log(ages)

    age_prior_model: dict
        dict including prior model name and parameters

    Returns
    -------
    age_weights : numpy vector
       weights needed according to the prior model
    """
    if age_prior_model["name"] == "flat" or age_prior_model["name"] == "flat_linear":
        age_weights = np.full(len(logages), 1.0)
    elif age_prior_model["name"] == "flat_log":
        # flat in log space means use the native log(age) grid spacing
        # thus the priors weights are the inverse of the grid weights
        # assumes the logace spacing is uniform
        age_weights = 1.0 / compute_age_grid_weights(logages)
    elif age_prior_model["name"] == "bins_histo":
        # check if all ages within interpolation range
        if np.all(
            [np.max(logages) <= x <= np.min(logages) for x in age_prior_model["values"]]
        ):
            raise ValueError(
                "Age prior weight error: Requested ages outside of model range"
            )

        # interpolate according to bins, assuming SFR constant from i to i+1
        # and allow for bin edges input
        if len(age_prior_model["values"]) == len(age_prior_model["logages"]) - 1:
            age_prior_model["values"].append(0.0)
        ageND = interp1d(
            age_prior_model["logages"], age_prior_model["values"], kind="zero"
        )
        age_weights = ageND(logages)
    elif age_prior_model["name"] == "bins_interp":
        # interpolate model to grid ages
        age_weights = np.interp(
            logages,
            np.array(age_prior_model["logages"]),
            np.array(age_prior_model["values"]),
        )
    elif age_prior_model["name"] == "exp":
        # assumes SFR(t) \propto e**(-t/tau) \propto e**(age/tau)
        # where age \propto -t (for age=t0-t) and tau in Gyr
        vals = (10 ** logages) / (age_prior_model["tau"] * 1e9)
        age_weights = np.exp(-1.0 * vals)
    else:
        raise NotImplementedError(
            "input age prior ''{}'' function not supported".format(
                age_prior_model["name"]
            )
        )

    # normalize to avoid numerical issues (too small or too large)
    age_weights /= np.average(age_weights)

    return age_weights


def imf_kroupa(in_x):
    """
    Compute a Kroupa IMF

    Parameters
    ----------
    in_x : numpy vector
      masses

    Returns
    -------
    imf : numpy vector
      unformalized IMF
    """
    # allows for single float or an array
    x = np.atleast_1d(in_x)

    m1 = 0.08
    m2 = 0.5
    alpha0 = -0.3
    alpha1 = -1.3
    alpha2 = -2.3
    imf = np.full((len(x)), 0.0)

    indxs, = np.where(x >= m2)
    if len(indxs) > 0:
        imf[indxs] = x[indxs] ** alpha2

    indxs, = np.where((x >= m1) & (x < m2))
    fac1 = (m2 ** alpha2) / (m2 ** alpha1)
    if len(indxs) > 0:
        imf[indxs] = (x[indxs] ** alpha1) * fac1

    indxs, = np.where(x < m1)
    fac2 = fac1 * ((m1 ** alpha1) / (m1 ** alpha0))
    if len(indxs) > 0:
        imf[indxs] = (x[indxs] ** alpha0) * fac2

    return imf


def imf_salpeter(x):
    """
    Compute a Salpeter IMF

    Parameters
    ----------
    x : numpy vector
      masses

    Returns
    -------
    imf : numpy vector
      unformalized IMF
    """
    return x ** (-2.35)


def imf_flat(x):
    """
    Compute a flat IMF (useful for simulations, not for normal BEAST runs)

    Parameters
    ----------
    x : numpy vector
      masses

    Returns
    -------
    imf : numpy vector
      unformalized IMF
    """
    return 1.0


def compute_mass_prior_weights(masses, mass_prior_model):
    """
    Compute the mass prior for the specificed model

    Parameters
    ----------
    masses : numpy vector
        masses

    mass_prior_model: dict
        dict including prior model name and parameters

    Returns
    -------
    mass_weights : numpy vector
      Unnormalized IMF integral for each input mass
      integration is done between each bin's boundaries
    """

    # sort the initial mass along this isochrone
    sindxs = np.argsort(masses)

    # Compute the mass bin boundaries
    mass_bounds = compute_bin_boundaries(masses[sindxs])

    # compute the weights = mass bin widths
    mass_weights = np.empty(len(masses))

    # integrate the IMF over each bin
    if mass_prior_model["name"] == "kroupa":
        imf_func = imf_kroupa
    elif mass_prior_model["name"] == "salpeter":
        imf_func = imf_salpeter
    elif mass_prior_model["name"] == "flat":
        imf_func = imf_flat
    else:
        raise NotImplementedError("input mass prior function not supported")

    # calculate the average prior in each mass bin
    for i in range(len(masses)):
        mass_weights[sindxs[i]] = (quad(imf_func, mass_bounds[i], mass_bounds[i + 1]))[
            0
        ] / (mass_bounds[i + 1] - mass_bounds[i])

    # normalize to avoid numerical issues (too small or too large)
    mass_weights /= np.average(mass_weights)

    return mass_weights


def compute_metallicity_prior_weights(mets, met_prior_model):
    """
    Computes the metallicity prior for the specified model

    Parameters
    ----------
    mets : numpy vector
        metallicities
    met_prior_model: dict
        dict including prior model name and parameters

    Returns
    -------
    metallicity_weights : numpy vector
       weights to provide the requested prior model
    """
    if met_prior_model["name"] == "flat":
        met_weights = np.full(len(mets), 1.0)
    else:
        raise NotImplementedError("input metallicity prior function not supported")

    # normalize to avoid numerical issues (too small or too large)
    met_weights /= np.average(met_weights)

    return met_weights


def compute_distance_prior_weights(dists, dist_prior_model):
    """
    Computes the distance prior for the specified model

    Parameters
    ----------
    dists : numpy vector
        distances
    dist_prior_model: dict
        dict including prior model name and parameters

    Returns
    -------
    dists_weights : numpy vector
       weights to provide the requested prior model
    """
    if dist_prior_model["name"] == "flat":
        dists_weights = np.full(len(dists), 1.0)
    else:
        raise NotImplementedError("input distance prior function not supported")

    # normalize to avoid numerical issues (too small or too large)
    dists_weights /= np.average(dists_weights)

    return dists_weights
