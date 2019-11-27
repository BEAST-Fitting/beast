"""
Functions for interacting with the BEAST model
"""

# system imports
from __future__ import (absolute_import, division, print_function)

# other package imports
import numpy as np
import h5py
from tqdm import tqdm


def read_lnp_data(filename, nstars=None, shift_lnp=True):
    """
    Read in the sparse lnp for all the stars in the hdf5 file

    Parameters
    ----------
    filename : string
       name of the file with the sparse lnp values

    nstars : int (default=None)
       if you want to check that the number of lnp values is correct, set this
       to the number of stars expected in the file

    shift_lnp : boolean (default=True)
        if True, shift lnp values to have a max of 0.0

    Returns
    -------
    lnp_data: dictonary
       contains arrays of the lnp values and indexs to the BEAST model grid
    """


    with h5py.File(filename, 'r') as lnp_hdf:

        if nstars is not None:
            if len(lnp_hdf.keys()) != nstars:
                raise ValueError(
                    "Error: number of stars not equal between nstars image and lnp file"
                )

        # initialize arrays
        # - find the length of the sparse likelihoods
        lnp_sizes = [lnp_hdf[sname]['lnp'].value.shape[0] for sname in lnp_hdf.keys()]
        # - set arrays to the maximum size
        lnp_vals = np.zeros((np.max(lnp_sizes), nstars), dtype=float) - np.inf
        lnp_indxs = np.zeros((np.max(lnp_sizes), nstars), dtype=int)

        # loop over all the stars (groups)
        for k, sname in enumerate(lnp_hdf.keys()):
            lnp_vals[:lnp_sizes[k], k] = lnp_hdf[sname]['lnp'].value
            lnp_indxs[:lnp_sizes[k], k] = np.int64(np.array(lnp_hdf[sname]['idx'].value))

        if shift_lnp:
            # shift the log(likelihood) values to have a max of 0.0
            #  ok if the same shift is applied to all stars in a pixel
            #  avoids numerical issues later when we go to intergrate probs
            lnp_vals -= np.max(lnp_vals)

    return {'vals': lnp_vals, 'indxs': lnp_indxs}


def read_noise_data(
    filename,
    param_list=['bias', 'completeness', 'error'],
    filter_col=None,
):
    """
    Read some or all of the noise model parameters, for one or all of the filters

    Parameters
    ----------
    filename : string
       name of the file with the BEAST observationmodel grid

    param_list : list of strings
       the set of parameters to extract

    filter_col : int (default=None)
        if set, only return the data for this column number

    noise_data
    -------
    beast_data: dictonary
       contains arrays of the noise parameters
    """
    noise_data = {}

    # open files for reading
    with h5py.File(noise_filename, 'r') as noise_hdf:

        # get beast physicsmodel params
        for param in tqdm(param_list, desc='reading beast data'):
            if filter_col is None:
                noise_data[cparam] = np.array(noise_hdf[cparam])
            else:
                noise_data[cparam] = noise_hdf[cparam][:,filter_col]

    return noise_data


def read_sed_data(
    filename,
    param_list=['Av', 'Rv', 'f_A', 'M_ini', 'logA', 'Z', 'distance'],
    verbose=True
):
    """
    Read in the beast data needed by all the pixels

    Parameters
    ----------
    filename : string
       name of the file with the BEAST physicsmodel grid

    param_list : list of strings
       the set of parameters to extract
       default = [Av, Rv, f_A, M_ini, logA, Z, distance]
       If set to None, return the list of possible parameters

    Returns
    -------
    grid_param_list : list of strings
        if param_list is None, return the list of parameter options

    beast_data: dictonary
       contains arrays of the requested SED grid parameters
    """
    sed_data = {}

    # open files for reading
    with h5py.File(sed_filename, 'r') as sed_hdf:

        # get the possible list of parameters
        grid_param_list = list(h['grid'].value.dtype.names)
        # return that if the user is so inclined
        if param_list is None:
            return grid_param_list + ['seds', 'lamb']

        # get parameters
        for param in tqdm(param_list, desc='reading beast data'):
            # grid parameter
            if param in grid_param_list:
                sed_data[cparam] = sed_hdf['grid'][param]
            # wavelengths of the filters -or- SED photometry values
            elif (param == 'lamb') or (param == 'seds'):
                sed_data[param] = sed_hdf[param].value
            else:
                raise ValueError(
                    'parameter {0} not found in SED grid'.format(param)
                )



    return sed_data


def extract_beast_data(beast_data, lnp_data):
    """
    Read in the beast data for the locations where the lnp values
    were saved

    Parameters
    ----------
    beast_data: dictonary
       contains arrays of the beast parameters and priors

    lnp_data: dictonary
       contains arrays of the lnp values and indexs to the BEAST model grid

    Returns
    -------
    lnp_grid_vals: dictonary
       contains arrays of the beast parameters and priors for the sparse
       lnp saved model grid points
    """
    # get the keys in beast_data
    beast_params = beast_data.keys()

    # setup the output
    lnp_grid_vals = {}
    n_lnps, n_stars = lnp_data['indxs'].shape
    for cparam in beast_params:
        lnp_grid_vals[cparam] = np.empty((n_lnps, n_stars), dtype=float)

    # loop over the stars and extract the requested BEAST data
    # for k in tqdm(range(n_stars), desc='extracting beast data'):
    for k in range(n_stars):
        for cparam in beast_params:
            lnp_grid_vals[cparam][:, k] = \
                            beast_data[cparam][lnp_data['indxs'][:, k]]

    return lnp_grid_vals
