import os
from multiprocessing import Pool

import h5py
import numpy as np

from ..observationmodel.noisemodel.generic_noisemodel import get_noisemodelcat
from ..physicsmodel import grid
from ..external.eztables import Table


def split_grid(grid_fname, num_subgrids):
    """
    Splits a spectral or sed grid (they are the same class actually)
    according to grid point index (so basically, arbitrarily).

    Parameters
    ----------
    grid_fname: string
        file name of the existing grid to be split up

    num_subgrids: integer
        the number of parts the grid should be split into

    Returns
    -------
    list of string
        the names of the newly created subgrid files

    """

    # With h5py we can choose which data we want to load to memory by
    # providing a slice
    h5grid = h5py.File(grid_fname)
    lamb = h5grid['lamb']
    seds = h5grid['seds']
    gr = h5grid['grid']

    fnames = []

    num_seds = seds.shape[0]
    q = num_seds // num_subgrids
    r = num_seds % num_subgrids
    for i in range(num_subgrids):

        subgrid_fname = grid_fname.replace('.hd5', 'sub{}.hd5'.format(i))
        fnames.append(subgrid_fname)
        if os.path.isfile(subgrid_fname):
            print('{} already exists. Skipping.'.format(subgrid_fname))
            continue
        else:
            print('constructing subgrid ' + str(i))

        # First, do strides of q+1
        if i < r:
            start = i * (q + 1)
            stop = start + q + 1
        # After the remainder has been taken care of, do strides of q
        else:
            start = r * (q + 1) + (i - r) * q
            stop = start + q

        # Load a slice as a SpectralGrid object
        slc = slice(start, stop)
        g = grid.SpectralGrid(lamb, seds=seds[slc], grid=Table(gr[slc]),
                              backend='memory')

        # Save it to a new file
        g.writeHDF(subgrid_fname, append=False)

    return fnames


def merge_grids(seds_fname, sub_names):
    """
    Merges a set of grids into one big grid. The grids need to have the
    same columns

    Parameters
    ----------
    seds_fname: string
        path for the output file

    sub_names: list of strings
        paths for the input grids
    """

    if not os.path.isfile(seds_fname):
        for n in sub_names:
            print('Appending {} to {}'.format(n, seds_fname))
            g = grid.FileSEDGrid(n)
            g.writeHDF(seds_fname, append=True)
    else:
        print('{} already exists'.format(seds_fname))


def gather_mins_maxes(grid_fname, noise_fname=None):
    """
    Generates a list of mins and maxes of all the quantities in the given grid

    Parameters
    ----------
    grid_fname: string
        path to a beast grid file (hd5 format)

    noise_fname: string
        Path to the noise model file for the given grid (hd5 format)
        (optional). If this is given, the mins/maxes for the full model
        fluxes are added too, under the name 'log'+filter+'_wd_bias'
        (needs to conform to the name used in fit.py).

    Returns
    -------
    min_max_dict: dictionary
        {name of quantity [string]: (minimum, maximum), ...}
    """

    # Use the HDFStore (pytables) backend
    sedgrid = grid.FileSEDGrid(grid_fname, backend='hdf')
    seds = sedgrid.seds

    min_max_dict = {}

    qnames = sedgrid.keys()
    for q in qnames:
        qvals = sedgrid[q]
        qmin = np.amin(qvals)
        qmax = np.amax(qvals)
        min_max_dict[q] = (qmin, qmax)

    if noise_fname is not None:
        # This code is more or less copied from fit.py
        noisemodel = get_noisemodelcat(noise_fname)
        full_model_flux = seds[:] + noisemodel.root.bias[:]
        indxs = np.where(full_model_flux > 0)
        full_model_flux[indxs] = np.log10(full_model_flux[indxs])
        full_model_flux[np.where(full_model_flux <= 0)] = -100.

        filters = sedgrid.filters
        for i, f in enumerate(filters):
            f_fluxes = full_model_flux[:, i]
            qmin = np.amin(f_fluxes)
            qmax = np.amax(f_fluxes)

            q = 'log'+f+'_wd_bias'
            qnames.append(q)

            min_max_dict[q] = (qmin, qmax)

    print('Got minima and maxima of {} for {}'.format(qnames, grid_fname))
    return min_max_dict

def reduce_mins_maxes(grid_fnames, noise_fnames=None, nprocs=1):
    """
    Computes the total minimum and maximum of the necessary quantities
    across all the subgrids. Can run in parallel.

    Parameters
    ----------
    grid_fnames: list of str
        subgrid file paths

    noise_fnames: list of str (optional)
        noise file for each subgrid

    nprocs: int
        Number of processes to use

    Returns
    -------
    min_max_dict: dictionary
        {name of quantity: (min, max), ...}
    """
    # Gather the mins and maxes for the subgrid
    if noise_fnames is None:
        arguments = [(g, None) for g in grid_fnames]
    else:
        arguments = list(zip(grid_fnames, noise_fnames))

    parallel = nprocs > 1
    if (parallel):
        p = Pool(nprocs)
        min_max_dicts = p.starmap(gather_mins_maxes, arguments)
    else:
        min_max_dicts = []
        for a in arguments:
            min_max_dicts.append(gather_mins_maxes(*a))

    # Assume that all dicts have the same keys. Copy the (min,max)
    # tuples from the first dict as a starting point.
    result_min_max_dict = min_max_dicts[0]

    # Then, reduce the values over the rest of the dicts
    for individual_min_max_dict in min_max_dicts:
        for q in result_min_max_dict:
            currentmin, currentmax = result_min_max_dict[q]
            othermin, othermax = individual_min_max_dict[q]
            result_min_max_dict[q] = (min(currentmin, othermin),
                                          max(currentmax, othermax))

    return result_min_max_dict
