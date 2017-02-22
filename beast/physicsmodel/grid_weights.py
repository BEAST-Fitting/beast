"""
Grid Weights
============
The use of a non-uniformly spaced grid complicates the marginalization 
step as the trick of summation instead of integration is used.  But this 
trick only works when the grid is uniformaly spaced in all dimensions.

If the grid is not uniformally spaced, weights can be used to correct
for the non-uniform spacing.  

Basically, we want the maginalization using these grid weights to provide
flat priors on all the fit parameters.  Non-flat priors will be implemented 
with prior weights.
"""
from __future__ import print_function

import numpy as np
from scipy.integrate import quad

from grid import FileSEDGrid
from grid import SpectralGrid
from ..external.eztables import Table

__all__ = ['compute_age_weights','compute_mass_weights',
           'compute_age_mass_metallicity_grid_weights']

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
    # ages need to be monotonically increasing
    aindxs = np.argsort(logages)   

    # Computes the bin boundaries in log
    logages_bounds = compute_bin_boundaries(logages[aindxs])    

    # initialize the age weights
    age_weights = np.full(len(aindxs),0.0)

    # Returns the age weight as a numpy array
    age_weights[aindxs] = np.diff(10**(logages_bounds))

    # return in the order that logages was passed
    return age_weights    

def compute_mass_weights(masses):
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

def compute_metallicity_weights(mets):
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

def compute_age_mass_metallicity_grid_weights(_tgrid):
    """ 
    Computes the age-mass-metallicity grid weights the BEAST model spectra grid

    Keywords
    --------
    _tgrid : BEAST model spectra grid.

    Returns
    -------
    Grid weight column is updated by multiplying by the 
       age-mass-metallicity weight.
    """

    # get the unique metallicities
    uniq_Zs = np.unique(_tgrid['Z'])  

    # setup the vector to hold the z weight vector
    total_z_weight = np.zeros(len(uniq_Zs))

    for az, z_val in enumerate(uniq_Zs):
        print('computing the age-mass grid weight for Z = ', z_val)
        
        # get the grid for a single metallicity
        zindxs, = np.where(_tgrid['Z'] == z_val)   

        # get the unique ages for this metallicity
        uniq_ages = np.unique(_tgrid[zindxs]['logA']) 

        # compute the age weights
        age_weights = compute_age_weights(uniq_ages)  

        for ak, age_val in enumerate(uniq_ages):
            # get the grid for a single age
            aindxs, = np.where((_tgrid['logA'] == age_val) &
                               (_tgrid['Z'] == z_val))   
            _tgrid_single_age = _tgrid[aindxs]

            # compute the mass weights
            if len(aindxs) > 1:
                mass_weights = compute_mass_weights(_tgrid_single_age['M_ini'])
            else:
                # must be a single mass for this age,z combination
                # set mass weight to zero to remove this point from the grid
                mass_weights = np.zeros(1)

            # apply both the mass and age weights  
            for i, k in enumerate(aindxs):
                _tgrid[k]['grid_weight'] *= mass_weights[i]*age_weights[ak]

        # compute the current total weight at each metallicity
        total_z_weight[az] = np.sum(_tgrid[zindxs]['weight'])
        
    # ensure that the metallicity prior is uniform
    if len(uniq_Zs) > 1:
        # get the metallicity weights
        met_weights = compute_metallicity_weights(uniq_Zs)

        total_z_weight *= met_widths   # very simple integration
        z_weights = total_z_weight/np.sum(total_z_weight)

        for az, z_val in enumerate(uniq_Zs):
            # get the grid for this metallicity
            zindxs, = np.where(_tgrid['Z'] == z_val)   
            _tgrid[zindxs]['grid_weight'] *= z_weights[az]

    # Add index for use later with the SED grid
    # =================================================
    # useful for looking up the best fit spectrum from a SED fit
    _tgrid[:]['specgrid_indx'] = np.arange(len(_tgrid), dtype=np.int64)
