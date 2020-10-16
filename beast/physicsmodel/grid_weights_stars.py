"""
Grid Weights
============
The use of a non-uniformly spaced grid complicates the marginalization
step as the trick of summation instead of integration is used.  But this
trick only works when the grid is uniformly spaced in all dimensions.

If the grid is not uniformly spaced, weights can be used to correct
for the non-uniform spacing.
"""
import numpy as np

__all__ = [
    "compute_distance_grid_weights",
    "compute_age_grid_weights",
    "compute_mass_grid_weights",
    "compute_metallicity_grid_weights",
    "compute_bin_boundaries",
]


def compute_bin_boundaries(tab):
    """
    Computes the boundaries of bins

    The bin boundaries are defined as the midpoint between each value in tab.
    At the two edges, 1/2 of the bin width is subtractted/added to the
    min/max of tab.

    Parameters
    ----------
    tab : numpy array
       centers of each bin

    Returns
    -------
    tab2 : numpy array
       boundaries of the bins
    """
    temp = tab[1:] - np.diff(tab) / 2.0
    tab2 = np.zeros(len(tab) + 1)
    tab2[0] = tab[0] - np.diff(tab)[0] / 2.0
    tab2[-1] = tab[-1] + np.diff(tab)[-1] / 2.0
    tab2[1:-1] = temp
    return tab2


def compute_age_grid_weights(logages):
    """
    Computes the age weights to set a uniform prior on linear SFR

    Parameters
    ----------
    logages : numpy vector
       log(ages)

    Returns
    -------
    age_weights : numpy vector
       total masses at each age for a constant SFR in linear age
    """
    # ages need to be monotonically increasing
    aindxs = np.argsort(logages)

    # Computes the bin boundaries in log
    logages_bounds = compute_bin_boundaries(logages[aindxs])

    # initialize the age weights
    age_weights = np.full(len(aindxs), 0.0)

    # Returns the age weight as a numpy array
    age_weights[aindxs] = np.diff(10 ** (logages_bounds))

    # normalize to avoid numerical issues (too small or too large)
    age_weights /= np.average(age_weights)

    # return in the order that logages was passed
    return age_weights


def compute_mass_grid_weights(masses):
    """
    Computes the mass weights to set a uniform prior on linear mass

    Parameters
    ----------
    masses : numpy vector
        masses

    Returns
    -------
    mass_weights : numpy vector
       weights to provide a constant SFR in linear age
    """
    # sort the initial mass along this isochrone
    sindxs = np.argsort(masses)

    # Compute the mass bin boundaries
    masses_bounds = compute_bin_boundaries(masses[sindxs])

    # compute the weights = bin widths
    mass_weights = np.zeros(len(masses))
    mass_weights[sindxs] = np.diff(masses_bounds)

    # normalize to avoid numerical issues (too small or too large)
    mass_weights /= np.average(mass_weights)

    return mass_weights


def compute_metallicity_grid_weights(mets):
    """
    Computes the metallicity weights to set a uniform prior on linear metallicity

    Parameters
    ----------
    mets : numpy vector
        metallicities

    Returns
    -------
    metallicity_weights : numpy vector
       weights to provide a flat metallicity
    """
    # sort the initial mass along this isochrone
    sindxs = np.argsort(mets)

    # Compute the mass bin boundaries
    mets_bounds = compute_bin_boundaries(mets[sindxs])

    # compute the weights = bin widths
    mets_weights = np.zeros(len(mets))
    mets_weights[sindxs] = np.diff(mets_bounds)

    # normalize to avoid numerical issues (too small or too large)
    mets_weights /= np.average(mets_weights)

    return mets_weights


def compute_distance_grid_weights(dists):
    """
    Computes the distance weights to set a uniform prior on linear distance

    Parameters
    ----------
    dists : numpy vector
        distances

    Returns
    -------
    dist_weights : numpy vector
       weights to provide a flat distance
    """
    # sort
    tdists = np.array(dists)
    sindxs = np.argsort(tdists)

    # Compute the bin boundaries
    dists_bounds = compute_bin_boundaries(tdists[sindxs])

    # compute the weights = bin widths
    dists_weights = np.zeros(len(tdists))
    dists_weights[sindxs] = np.diff(dists_bounds)

    # normalize to avoid numerical issues (too small or too large)
    dists_weights /= np.average(dists_weights)

    return dists_weights
