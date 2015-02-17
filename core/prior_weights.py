"""
Priors as weights
================

The priors on age and mass are computed as weights to be used in the
likelihood computation.  This code was created by Heddy Arab and
integrated into the beast core by Karl Gordon.
"""
from __future__ import print_function

import numpy as np
from scipy.integrate import quad

from .grid import FileSEDGrid
from .grid import SpectralGrid
from ..external.eztables import Table

# compute the bin edges
# approximate the edge bins by adding 1/2 the adjacent bin width
def compute_bin_boundaries(tab):
    """ Computes the bin boundaries.

    Keywords
    ----------
    tab : table of float values

    Returns
    -------
    table of bounardies to the input tab

    Note
    ----
    The bin boundaries are defined as the midpoint between each value in tab.
    At the two edges, 1/2 of the bin width is subtractted/added to the min/max of tab.
    """
    temp = tab[1:]-np.diff(tab)/2.
    tab2 = np.empty(len(tab)+1)
    tab2[0] = tab[0]-np.diff(tab)[0]/2.
    tab2[-1] = tab[-1]+np.diff(tab)[-1]/2.
    tab2[1:-1] = temp
    return tab2

# Kroupa IMF
def imf_kroupa(x):
    """ Computes a Kroupa IMF

    Keywords
    ----------
    x : masses

    Returns
    -------
    Unnormalized IMF based on the input masses
    """
    m0 = 0.01
    m1 = 0.08
    m2 = 0.5
    alpha0 = -0.3
    alpha1 = -1.3
    alpha2 = -2.3
    if (x < m1):
        return x**alpha0
    elif (x >= m1) and (x < m2):
        return x**alpha1
    elif (x>=m2):
        return x**alpha2
    
# Salpeter IMF
def imf_salpeter(x):
    """ Computes a Salpeter IMF

    Keywords
    ----------
    x : masses

    Returns
    -------
    Unnormalized IMF based on the input masses
    """
    return x**(-2.35) # Salpeter IMF

# compute the age weights for a constant SFR in linear age
def compute_age_weights(logages):
    """ Computes the age weights assuming a constant star formation rate (linear age)

    Keywords
    --------
    logages : vector of log(ages)

    Returns
    -------
    Unnormalized total masses at each age assuming a constant SFR.
    """
    aindxs = np.argsort(logages)   # ages need to be monotonically increasing
    logages2 = compute_bin_boundaries(logages[aindxs])    # Computes the bin boundaries in log
    age_weights = np.full(len(aindxs),0.0)
    age_weights[aindxs] = np.diff(10**(logages2))           # Returns the age weight as a numpy array
    return age_weights    # return in the order that logages was passed

# compute the mass weights at a constant age
# uses an assumed IMF to generate the weights
#  IMF norm is the integral of the assumed IMF over the full possible mass range
#  must be precomputed as this information is not usual available for a specific isochrone age
def compute_mass_weights(masses, full_imf_integral):
    """ Computes the mass weights for a Kroupa IMF.

    Keywords
    --------
    masses : vector of masses

    Returns
    -------
    Unnormalized integral of the IMF for each input mass.  Integration is done between the
    bin min/max boundary.

    Notes
    -----
    Should make the IMF function used a passed variable.
    """
    d = np.zeros(len(masses))
        
    isoc = np.sort(masses)               # sort the initial mass along this isochrone
    index_isoc = np.argsort(masses)
    
    isoc2 = compute_bin_boundaries(isoc)      # Compute the initial mass bin boundaries

    I1 = np.empty(len(isoc))
    for ik, uk in enumerate(isoc2[:-1]):
        res = quad(imf_kroupa, isoc2[ik], isoc2[ik+1]) # integrate according to the prior on the mass bin
        I1[index_isoc[ik]] = res[0]/full_imf_integral      # compute the integrated weight
    return I1

# compute age-mass-metallicity prior weights
# age prior is a constant SFR in linear age
# mass prior is a Kroupa IMF (need to update code to allow user to pick the function)
# metallicity prior is flat
def compute_age_mass_metallicity_prior_weights(_tgrid):
    """ Computes the age-mass-metallicity weights on a BEAST model spectra grid.

    Keywords
    --------
    _tgrid : BEAST model spectra grid.

    Returns
    -------
    Nothing.
    Weight column in the grid is updated by multiplying by the age-mass-metallicity
    weight.

    """

    uniq_Zs = np.unique(_tgrid['Z'])  # get the unique metallicities
    total_z_weight = np.zeros(len(uniq_Zs))
    for az, z_val in enumerate(uniq_Zs):
        print('working computing the age-mass prior for Z = ', z_val)
        
        zindxs, = np.where(_tgrid['Z'] == z_val)   # get the grid for a single metallicity
        uniq_ages = np.unique(_tgrid[zindxs]['logA']) # get the unique ages for this metallicity
        age_weights = compute_age_weights(uniq_ages)  # compute the age weights for a constant SFR in linear age

        # get the integral over the IMF
        # assuming the min/max masses at any age for a single metallicity are the assumed min/max masses
        #min_mass = np.min(_tgrid[zindxs]['M_ini'])
        #max_mass = np.max(_tgrid[zindxs]['M_ini'])
        #imf_norm = quad(imf_kroupa, min_mass, max_mass)   # integrate according to the desired IMF along the isochrone
        imf_norm = [1.0]  # effectively ignore this term - assumes same mass range for all metallicities

        for ak, age_val in enumerate(uniq_ages):
            aindxs, = np.where((_tgrid['logA'] == age_val) & (_tgrid['Z'] == z_val))   # get the grid for a single age
            _tgrid_single_age = _tgrid[aindxs]
            if len(aindxs) > 1:
                mass_weights = compute_mass_weights(_tgrid_single_age['M_ini'],imf_norm[0])
            else:
                # must be a single mass for this age,z combination
                # set mass weight to zero to remove this point from the grid
                mass_weights = np.zeros(1)

            for i, k in enumerate(aindxs):
                _tgrid[k]['weight'] *= mass_weights[i]*age_weights[ak]

        total_z_weight[az] = np.sum(_tgrid[zindxs]['weight'])
        
    # ensure that the metallicity prior is uniform
    if len(uniq_Zs) > 1:
        z_boundaries = compute_bin_boundaries(uniq_Zs)
        z_widths = np.diff(z_boundaries)

        total_z_weight *= z_widths   # very simple integration
        z_weights = total_z_weight/np.sum(total_z_weight)
        for az, z_val in enumerate(uniq_Zs):
            zindxs, = np.where(_tgrid['Z'] == z_val)   # get the grid for a single metallicity
            _tgrid[zindxs]['weight'] *= z_weights[az]

    # Add index for use later with the SED grid
    # =================================================
    # useful for looking up the best fit spectrum from a SED fit
    _tgrid[:]['specgrid_indx'] = np.arange(len(_tgrid), dtype=np.int64)
