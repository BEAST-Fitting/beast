"""
Grid Weights
============
The use of a non-uniformly spaced grid complicates the marginalization 
step as the trick of summation instead of integration is used.  But this 
trick only works when the grid is uniformly spaced in all dimensions.

If the grid is not uniformly spaced, weights can be used to correct
for the non-uniform spacing.  

Basically, we want the maginalization using these grid weights to provide
flat priors on all the fit parameters.  Non-flat priors will be implemented 
with prior weights.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

__all__ = ['compute_age_grid_weights',
           'compute_mass_grid_weights',
           'compute_metallicity_grid_weights']

def compute_bin_boundaries(tab):
    """ Computes the bin boundaries.

    Parameters
    ----------
    tab : numpy array
       centers of each bin

    Returns
    -------
    tab2 : numpy array
       boundaries of the bins

    Note
    ----
    The bin boundaries are defined as the midpoint between each value in tab.
    At the two edges, 1/2 of the bin width is subtractted/added to the
    min/max of tab.
    """
    temp = tab[1:]-np.diff(tab)/2.
    tab2 = np.empty(len(tab)+1)
    tab2[0] = tab[0]-np.diff(tab)[0]/2.
    tab2[-1] = tab[-1]+np.diff(tab)[-1]/2.
    tab2[1:-1] = temp
    return tab2

def compute_age_grid_weights(logages, constantSFR=True):
    """ Computes the age weights to set prior on star formation history

    Keywords
    --------
    logages : numpy vector
       log(ages)

    constantSFR : boolean
       Sets assumption of star formation history: flat in log or linear age
       Default = True (constant in linear age)

    Returns
    -------
    age_weights : numpy vector
       total masses at each age for a constant SFR in linear age
    """
    if constantSFR:
        # ages need to be monotonically increasing
        aindxs = np.argsort(logages)   
        
        # Computes the bin boundaries in log
        logages_bounds = compute_bin_boundaries(logages[aindxs])    
        
        # initialize the age weights
        age_weights = np.full(len(aindxs),0.0)
        
        # Returns the age weight as a numpy array
        age_weights[aindxs] = np.diff(10**(logages_bounds))

    else:
        # Returns the age weight as a numpy array
        age_weights = np.ones(len(logages))

    # return in the order that logages was passed
    return age_weights    

def compute_mass_grid_weights(masses):
    """ Computes the mass weights to provide a uniform IMF

    Keywords
    --------
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

    # compute the weights = mass bin widths
    mass_weights = np.empty(len(masses))
    mass_weights[sindxs] = np.diff(masses_bounds)

    return mass_weights

def compute_metallicity_grid_weights(mets):
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
    # sort the initial mass along this isochrone
    sindxs = np.argsort(mets)
    
    # Compute the mass bin boundaries
    mets_bounds = compute_bin_boundaries(mets[sindxs])      

    # compute the weights = mass bin widths
    mets_weights = np.empty(len(mets))
    mets_weights[sindxs] = np.diff(mets_bounds)

    return mets_weights

