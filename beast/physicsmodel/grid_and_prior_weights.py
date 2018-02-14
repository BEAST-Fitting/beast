"""
Grid and Prior Weights
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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from .grid_weights import compute_age_grid_weights
from .grid_weights import compute_mass_grid_weights
from .grid_weights import compute_metallicity_grid_weights

from .prior_weights import compute_age_prior_weights
from .prior_weights import compute_mass_prior_weights
from .prior_weights import compute_metallicity_prior_weights

__all__ = ['compute_age_mass_metallicity_weights',
           'compute_bin_boundaries']

def compute_age_mass_metallicity_weights(_tgrid, **kwargs):
    """ 
    Computes the age-mass-metallicity grid and prior weights 
    on the BEAST model spectra grid

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
    total_z_grid_weight = np.zeros(len(uniq_Zs))
    total_z_prior_weight = np.zeros(len(uniq_Zs))
    total_z_weight = np.zeros(len(uniq_Zs))

    for az, z_val in enumerate(uniq_Zs):
        print('computing the age-mass-metallicity grid weight for Z = ',
              z_val)
        
        # get the grid for a single metallicity
        zindxs, = np.where(_tgrid['Z'] == z_val)   

        # get the unique ages for this metallicity
        uniq_ages = np.unique(_tgrid[zindxs]['logA']) 

        # compute the age weights
        age_grid_weights = compute_age_grid_weights(uniq_ages, **kwargs)  
        age_prior_weights = compute_age_prior_weights(uniq_ages)  

        for ak, age_val in enumerate(uniq_ages):
            # get the grid for a single age
            aindxs, = np.where((_tgrid['logA'] == age_val) &
                               (_tgrid['Z'] == z_val))   
            _tgrid_single_age = _tgrid[aindxs]

            # compute the mass weights
            if len(aindxs) > 1:
                cur_masses = _tgrid_single_age['M_ini']
                mass_grid_weights = compute_mass_grid_weights(cur_masses)
                mass_prior_weights = compute_mass_prior_weights(cur_masses)
            else:
                # must be a single mass for this age,z combination
                # set mass weight to zero to remove this point from the grid
                mass_grid_weights = np.zeros(1)
                mass_prior_weights = np.zeros(1)

            # apply both the mass and age weights  
            for i, k in enumerate(aindxs):
                _tgrid[k]['grid_weight'] *= mass_grid_weights[i] \
                                            *age_grid_weights[ak]
                _tgrid[k]['prior_weight'] *= mass_prior_weights[i] \
                                             *age_prior_weights[ak]
                _tgrid[k]['weight'] *= mass_prior_weights[i] \
                                       *age_grid_weights[ak]

        # compute the current total weight at each metallicity
        total_z_grid_weight[az] = np.sum(_tgrid[zindxs]['grid_weight'])
        total_z_prior_weight[az] = np.sum(_tgrid[zindxs]['prior_weight'])
        total_z_weight[az] = np.sum(_tgrid[zindxs]['weight'])
        
    # ensure that the metallicity prior is uniform
    if len(uniq_Zs) > 1:
        # get the metallicity weights
        met_grid_weights = compute_metallicity_grid_weights(uniq_Zs)
        met_grid_weights /= np.sum(met_grid_weights)
        met_prior_weights = compute_metallicity_prior_weights(uniq_Zs)
        met_prior_weights /= np.sum(met_prior_weights)
        met_weights = met_grid_weights*met_prior_weights

        # correct for any non-unformity in the number size of the
        # age-mass grids between metallicity points
        total_z_grid_weight /= np.sum(total_z_grid_weight)
        total_z_prior_weight /= np.sum(total_z_prior_weight)
        total_z_weight /= np.sum(total_z_weight)
        
        for i, z_val in enumerate(uniq_Zs):
            # get the grid for this metallicity
            zindxs, = np.where(_tgrid['Z'] == z_val)   
            _tgrid[zindxs]['grid_weight'] *= met_grid_weights[i] \
                                             *total_z_grid_weight[i]
            _tgrid[zindxs]['prior_weight'] *= met_prior_weights[i] \
                                              *total_z_prior_weight[i]
            _tgrid[zindxs]['weight'] *= met_weights[i] \
                                        *total_z_weight[i]

