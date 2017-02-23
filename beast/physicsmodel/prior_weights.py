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

from grid import FileSEDGrid
from grid import SpectralGrid
from ..external.eztables import Table

from grid_weights import compute_bin_boundaries

__all__ = ['compute_age_mass_metallicity_prior_weights']

def compute_age_weights(logages):
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

def compute_mass_weights(masses):
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

# compute age-mass-metallicity prior weights
# age prior is a constant SFR in linear age
# mass prior is a Kroupa IMF (need to update code to allow user to pick
#   the function)
# metallicity prior is flat
def compute_age_mass_metallicity_prior_weights(_tgrid):
    """ Computes the age-mass-metallicity weights on a BEAST model spectra grid.

    Keywords
    --------
    _tgrid : BEAST model spectra grid.

    Returns
    -------
    Weight column in the grid is updated by multiplying by the 
    age-mass-metallicity weight.
    """

    # get the unique metallicities
    uniq_Zs = np.unique(_tgrid['Z'])  
    total_z_weight = np.zeros(len(uniq_Zs))
    for az, z_val in enumerate(uniq_Zs):
        print('working computing the age-mass prior for Z = ', z_val)
        
        # get the grid for a single metallicity
        zindxs, = np.where(_tgrid['Z'] == z_val)   
        # get the unique ages for this metallicity
        uniq_ages = np.unique(_tgrid[zindxs]['logA']) 
        # compute the age weights for a constant SFR in linear age
        age_weights = compute_age_weights(uniq_ages)  

        # assumes same mass range for all metallicities
        imf_norm = [1.0]  

        for ak, age_val in enumerate(uniq_ages):
            # get the grid for a single age
            aindxs, = np.where((_tgrid['logA'] == age_val) &
                               (_tgrid['Z'] == z_val))   
            _tgrid_single_age = _tgrid[aindxs]
            if len(aindxs) > 1:
                mass_weights = compute_mass_weights(_tgrid_single_age['M_ini'])
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
            # get the grid for a single metallicity
            zindxs, = np.where(_tgrid['Z'] == z_val)   
            _tgrid[zindxs]['weight'] *= z_weights[az]

