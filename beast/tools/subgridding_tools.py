import os
from multiprocessing import Pool

import h5py
import numpy as np
from astropy.io import fits
from astropy.table import Table

from ..observationmodel.noisemodel.generic_noisemodel import get_noisemodelcat
from ..physicsmodel import grid
from ..external import eztables
from ..fitting.fit import save_pdf1d


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
        g = grid.SpectralGrid(lamb, seds=seds[slc], grid=eztables.Table(gr[slc]),
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


def subgrid_info(grid_fname, noise_fname=None):
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
    info_dict: dictionary
        {name of quantity [string]: {'min': min, 'max': max, 'unique': unique values}}
    """

    # Use the HDFStore (pytables) backend
    sedgrid = grid.FileSEDGrid(grid_fname, backend='hdf')
    seds = sedgrid.seds

    info_dict = {}

    qnames = sedgrid.keys()
    for q in qnames:
        qvals = sedgrid[q]
        qmin = np.amin(qvals)
        qmax = np.amax(qvals)
        qunique = np.unique(qvals)
        info_dict[q] = {}
        info_dict[q]['min'] = qmin
        info_dict[q]['max'] = qmax
        info_dict[q]['unique'] = qunique

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
            qunique = np.unique(qvals)

            q = 'log'+f+'_wd_bias'
            info_dict[q] = {}
            info_dict[q]['min'] = qmin
            info_dict[q]['max'] = qmax
            info_dict[q]['unique'] = qunique

    print('Gathered grid info for {}'.format(grid_fname))
    return info_dict

def unpack_and_subgrid_info(x):
    """
    Utility to call this function in parallel, with multiple arguments
    """
    return subgrid_info(*x)

def reduce_grid_info(grid_fnames, noise_fnames=None, nprocs=1):
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
    info_dict: dictionary
        {name of quantity: (min, max), ...}
    """
    # Gather the mins and maxes for the subgrid
    if noise_fnames is None:
        arguments = [(g, None) for g in grid_fnames]
    else:
        arguments = list(zip(grid_fnames, noise_fnames))

    # Use generators here for memory efficiency
    parallel = nprocs > 1
    if (parallel):
        p = Pool(nprocs)
        info_dicts_generator = p.imap(unpack_and_subgrid_info, arguments)
    else:
        info_dicts_generator = (subgrid_info(*a) for a in arguments)

    # Then, reduce the values over the rest of the dicts
    result_dict = {}

    # Assume that all info dicts have the same keys
    first_info_dict = next(info_dicts_generator)
    for q in first_info_dict:
        # Combine the values of the first subgrid
        union_min = first_info_dict[q]['min']
        union_max = first_info_dict[q]['max']
        union_unique = first_info_dict[q]['unique']

        # And all the other subgrids
        for individual_dict in info_dicts_generator:
            other_min = individual_dict[q]['min']
            union_min = min(union_min, other_min)

            other_max = individual_dict[q]['max']
            union_max = max(union_max, other_max)

            other_unique = individual_dict[q]['unique']
            union_unique = np.union1d(union_unique, other_unique)

        result_dict[q] = {}
        result_dict[q]['min'] = union_min
        result_dict[q]['max'] = union_max
        result_dict[q]['num_unique'] = len(union_unique)

    return result_dict

def merge_pdf1d_stats(subgrid_pdf1d_fnames, subgrid_stats_fnames):
    """
    Merge a set of 1d pdfs that were generated by fits on different
    grids. It is nessecary (and checked) that all the 1d pdfs have the
    same limits, bin values, and number of bins.

    The stats files are also combined; some values for the total grid
    can be calculated by simply comparing them across all the grids,
    others are recalculated after obtaining the new 1dpdfs.

    Parameters
    ----------
    subgrid_pdf1d_fnames: list of string
        file names of all the pdf1d fits files

    subgrid_stats_fnames: list of string
        file names of the stats files. Should be in the same order as
        subgrid_pdf1d_fnames. These files are needed to help with
        averaging the pdf1d files as they contain the total weight of
        each subgrid.

    Returns
    -------
    merged_pdf1d_fname: string
        file name of the resulting pdf1d fits file (newly created by
        this function)
    """

    nsubgrids = len(subgrid_pdf1d_fnames)
    assert(len(subgrid_stats_fnames) == nsubgrids)

    # Open the first one, and check the contents of the others against
    # it. Gather some useful information along the way.
    nbins = {}
    with fits.open(subgrid_pdf1d_fnames[0]) as hdul_0:
        qnames = [hdu.name for hdu in hdul_0[1:]]
        nbins = {q: hdul_0[q].data.shape[1] for q in qnames}
        bincenters = {q: hdul_0[q].data[-1, :] for q in qnames}
        nobs = hdul_0[qnames[0]].data.shape[0] - 1

        # Check if all the bins are equal
        for pdf1d_f in subgrid_pdf1d_fnames[1:]:
            with fits.open(pdf1d_f) as hdul:
                for q in qnames:
                    pdf1d_0 = hdul_0[q].data
                    pdf1d = hdul[q].data
                    # the number of bins
                    assert(pdf1d_0.shape[1] == pdf1d.shape[1])
                    # the number of stars
                    assert(pdf1d_0.shape[0] == pdf1d.shape[0])
                    # the bin centers (stored in the last row of the
                    # image) should be equal (or both nan, in case there
                    # is only one bin; this is a glitch in the
                    # implementation, but no reason for alarm).
                    bin_centers_ok = \
                    np.isnan(pdf1d_0[-1,0]) and np.isnan(pdf1d[-1,0]) \
                    or (pdf1d_0[-1, :] == pdf1d[-1, :]).all()
                    assert(bin_centers_ok)

    # First, let's read the arrays of weights (each subgrid has an array
    # of weights, containing one weight for each source).
    logweight = np.zeros((nobs, nsubgrids))
    for i, stats_f in enumerate(subgrid_stats_fnames):
        stats = Table.read(stats_f)
        logweight[:, i] = stats['total_log_norm']

    # Take the maximum over all the grids. This will be used to
    # normalize the weights, per star.
    max_weights = np.amax(logweight, axis=1, keepdims=True)

    # Get linear weights for each object/grid. By using keepdims above,
    # max_weights is a column, and the subtraction will be done for each
    # column.
    weights = np.exp(logweight - max_weights)

    # We will try to reuse the save function defined in fit.py
    save_pdf1d_vals = []
    for i, q in enumerate(qnames):
        # Prepare the ouput array
        save_pdf1d_vals.append(np.zeros((nobs+1, nbins[q])))
        # Copy the bin centers
        save_pdf1d_vals[i][-1, :] = bincenters[q]

    # Now, go over all the pdf1d files, and sum the weighted pdf1d values
    for g, pdf1d_f in enumerate(subgrid_pdf1d_fnames):
        with fits.open(pdf1d_f) as hdul:
            for i, q in enumerate(qnames):
                pdf1d_g = hdul[q].data[:-1, :]
                weight_column = weights[:, [g]] # use [g] to keep dimension
                save_pdf1d_vals[i][:-1, :] += pdf1d_g * weight_column

    # Normalize all the pdfs of the final result
    for i in range(len(save_pdf1d_vals)):
        # sum for each source in a column
        norms_col = np.sum(save_pdf1d_vals[i][:-1, :], axis=1, keepdims=True)
        # non zero mask as 1d array
        nonzero = norms_col[:,0] > 0
        save_pdf1d_vals[i][:-1][nonzero, :] /= norms_col[nonzero]

    # Save the combined 1dpdf file
    save_pdf1d("combined_pdf1d.fits", save_pdf1d_vals, qnames)
