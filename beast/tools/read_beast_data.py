"""
Functions for interacting with the BEAST model
"""

import numpy as np
import h5py
from tqdm import tqdm


__all__ = ["read_lnp_data", "read_noise_data", "read_sed_data", "get_lnp_grid_vals"]


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
    lnp_data : dictonary
       contains arrays of the lnp values and indices to the BEAST model grid
    """

    with h5py.File(filename, "r") as lnp_hdf:

        # get keyword names for the stars (as opposed to filter info)
        star_key_list = [sname for sname in lnp_hdf.keys() if "star" in sname]
        tot_stars = len(star_key_list)

        if nstars is not None:
            if tot_stars != nstars:
                raise ValueError(
                    "Error: number of stars not equal between nstars image and lnp file"
                )

        # initialize arrays
        # - find the lengths of the sparse likelihoods
        lnp_sizes = [lnp_hdf[sname]["lnp"][()].shape[0] for sname in star_key_list]
        # - set arrays to the maximum size
        lnp_vals = np.full((np.max(lnp_sizes), tot_stars), -np.inf)
        lnp_indxs = np.full((np.max(lnp_sizes), tot_stars), -1)

        # loop over all the stars (groups)
        for k, sname in enumerate(star_key_list):
            lnp_vals[: lnp_sizes[k], k] = lnp_hdf[sname]["lnp"][()]
            lnp_indxs[: lnp_sizes[k], k] = np.array(lnp_hdf[sname]["idx"][()])

        if shift_lnp:
            # shift the log(likelihood) values to have a max of 0.0
            #  ok if the same shift is applied to all stars in a pixel
            #  avoids numerical issues later when we go to intergrate probs
            lnp_vals -= np.max(lnp_vals)

    return {"vals": lnp_vals, "indxs": lnp_indxs}


def read_noise_data(
    filename,
    param_list=["bias", "completeness", "error"],
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

    Returns
    -------
    noise_data : dictonary
       contains arrays of the noise parameters
    """
    noise_data = {}

    # open files for reading
    with h5py.File(filename, "r") as noise_hdf:

        # get beast physicsmodel params
        for param in tqdm(param_list, desc="reading noise data"):
            if filter_col is None:
                noise_data[param] = np.array(noise_hdf[param])
            else:
                noise_data[param] = noise_hdf[param][:, filter_col]

    return noise_data


def read_sed_data(
    filename,
    param_list=["Av", "Rv", "f_A", "M_ini", "logA", "Z", "distance"],
    return_params=False,
):
    """
    Read in the beast data needed by all the pixels

    Parameters
    ----------
    filename : string
        name of the file with the BEAST physicsmodel grid

    param_list : list of strings
        The set of parameters to extract (default: Av, Rv, f_A, M_ini, logA, Z,
        distance).  If set to 'all', extract all parameters and model fluxes in
        the grid.

    return_params : boolean (default=False)
        If True, return the list of keywords for all parameters and model fluxes
        in the grid.  Useful for checking what columns are present.

    Returns
    -------
    Two possible returns depending on param_list input

    sed_data : dictonary (param_list input as list of strings)
        contains arrays of the requested SED grid parameters
        if return_params is True, then also provides

    grid_param_list : list of strings (return_params is True)
        if param_list is None, return the list of parameter options
    """
    sed_data = {}

    # open files for reading
    with h5py.File(filename, "r") as sed_hdf:

        # get the possible list of parameters
        grid_param_list = list(sed_hdf["grid"][()].dtype.names)
        # return that if the user is so inclined
        if return_params:
            return grid_param_list + ["seds", "lamb", "filters"]
        if param_list == "all":
            param_list = grid_param_list

        # get parameters
        for param in tqdm(param_list, desc="reading sed data"):
            # grid parameter
            if param in grid_param_list:
                sed_data[param] = sed_hdf["grid"][param]
            # wavelengths of the filters -or- SED photometry values
            elif (param == "lamb") or (param == "seds"):
                sed_data[param] = sed_hdf[param][()]
            elif param == "filters":
                filtstr = sed_hdf["grid"].attrs["filters"]
                if isinstance(filtstr, bytes):
                    filtstr = filtstr.decode()
                sed_data[param] = filtstr.split(" ")
            else:
                raise ValueError("parameter {0} not found in SED grid".format(param))

    return sed_data


def get_lnp_grid_vals(sed_data, lnp_data, verbose=False):
    """
    Acquire the SED parameter values for the locations where the lnp values
    were saved

    Parameters
    ----------
    sed_data : dictonary or string
       if dictionary: contains arrays of the beast parameters (output from
       read_sed_data)
       if string: name of the file with the BEAST physicsmodel grid, which will
       be used in read_sed_data to get default parameters

    lnp_data : dictonary or string
       if dictionary: contains arrays of the lnp values and indices to the BEAST
       model grid (output from read_lnp_data)
       if string: name of the file with the sparse lnp values, which will be
       used in read_lnp_data with default parameters

    Returns
    -------
    lnp_grid_vals : dictonary
        arrays of the SED grid parameters for the points in the lnp lists
    """

    if isinstance(sed_data, str):
        sed_data = read_sed_data(sed_data)
    if isinstance(lnp_data, str):
        lnp_data = read_lnp_data(lnp_data)

    # get the keys in beast_data
    param_list = sed_data.keys()

    # setup the output
    lnp_grid_vals = {}
    n_lnps, n_stars = lnp_data["indxs"].shape
    for param in param_list:
        lnp_grid_vals[param] = np.zeros((n_lnps, n_stars), dtype=float)

    # loop over the stars and extract the requested BEAST data
    for k in tqdm(
        range(n_stars), desc="extracting params for each lnP", disable=not verbose
    ):
        lnp_inds = lnp_data["indxs"][:, k]
        good_inds = np.isfinite(lnp_inds)
        for param in param_list:
            lnp_grid_vals[param][good_inds, k] = sed_data[param][
                lnp_inds[good_inds].astype(int)
            ]

    return lnp_grid_vals
