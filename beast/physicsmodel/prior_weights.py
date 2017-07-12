"""
Prior Weights
=============
The priors on age, mass, and metallicty are computed as weights to use
in the posterior calculations.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from scipy.integrate import quad

from .grid_weights import compute_bin_boundaries

__all__ = ['compute_age_prior_weights',
           'compute_mass_prior_weights',
           'compute_metallicity_prior_weights']

def compute_age_prior_weights(logages):
    """ Computes the age weights to provide constant star formation rate
    (in linear age)

    Keywords
    --------
    logages : numpy vector
       log(ages)

    Returns
    -------
    age_weights : numpy vector
       total masses at each age for a constant SFR in linear age
    """
    # initialize the age weights to one
    #   for a flat prior, nothing else is needed
    #   non-uniform grid spacing is handled in the grid_weights code
    age_weights = np.full(len(logages),1.0)

    # code will be needed here for non-flat priors
    
    # return in the order that logages was passed
    return age_weights    

def imf_kroupa(x):
    """ Computes a Kroupa IMF

    Keywords
    ----------
    x : numpy vector
      masses

    Returns
    -------
    imf : numpy vector
      unformalized IMF
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
    
def imf_salpeter(x):
    """ Computes a Salpeter IMF

    Keywords
    ----------
    x : numpy vector
      masses

    Returns
    -------
    imf : numpy vector
      unformalized IMF
    """
    return x**(-2.35)

def compute_mass_prior_weights(masses):
    """ Computes the mass weights for a kroupa IMF

    Keywords
    --------
    masses : numpy vector
        masses

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
    for i in range(len(masses)):
        mass_weights[sindxs[i]] = (quad(imf_kroupa,
                                        mass_bounds[i],
                                        mass_bounds[i+1]))[0]

    return mass_weights

def compute_metallicity_prior_weights(mets):
    """ Computes the metallicity weights to provide a default flat prior

    Keywords
    --------
    mets : numpy vector
        metallicities

    Returns
    -------
    metallicity_weights : numpy vector
       weights to provide a flat metallicity
    """
    # initialize the metalicity weights to one
    #   for a flat prior, nothing else is needed
    #   non-uniform grid spacing is handled in the grid_weights code
    met_weights = np.full(len(mets),1.0)

    return met_weights
