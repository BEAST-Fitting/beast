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
        lnp_hdf.close()

        if shift_lnp:
            # shift the log(likelihood) values to have a max of 0.0
            #  ok if the same shift is applied to all stars in a pixel
            #  avoids numerical issues later when we go to intergrate probs
            lnp_vals -= np.max(lnp_vals)

    return {'vals': lnp_vals, 'indxs': lnp_indxs}


def read_beast_data(
    sed_filename,
    noise_filename,
    beast_params=['Av', 'Rv', 'f_A', 'M_ini', 'logA', 'Z', 'distance', 'completeness'],
    verbose=True
):
    """
    Read in the beast data needed by all the pixels

    Parameters
    ----------
    sed_filename: string
       name of the file with the BEAST physicsmodel grid

    noise_filename: string
       name of the file with the BEAST observationmodel grid

    beast_params: list of strings
       contains the set of BEAST parameters to extract
       default = [completeness, Av, Rv, f_A, M_ini, logA, Z, distance]

    Returns
    -------
    beast_data: dictonary
       contains arrays of the beast parameters, priors, and completeness
    """
    beast_data = {}

    # open files for reading
    with h5py.File(noise_filename, 'r') as beast_noise_hdf, h5py.File(sed_filename, 'r') as beast_seds_hdf:

        # get beast physicsmodel params
        for cparam in tqdm(beast_params, desc='reading beast data'):
            if cparam == 'completeness':
                beast_data[cparam] = np.max(beast_noise_hdf[cparam], axis=1)
            else:
                beast_data[cparam] = beast_seds_hdf['grid'][cparam]


    return beast_data


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
