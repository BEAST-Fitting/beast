from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

import h5py
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from astropy.table import Column

from ..vega import Vega


def mag_limits(seds, limits, Nfilter=1):
    """
    Selects models which have at least N filter above the limits

    Parameters
    ----------
    seds:   np.array
            Magnitude array from BEAST grid

    limits: list
            List of limit magnitudes

    Nfilter: integer
             In how many filters, you want a fake star to be brighter
             than the limit

    Returns
    -------
    idx:    np.array
            Array of integers contining the indices of allowed models

    """
    flag = seds.copy()

    # flag is True if the models are brigter (=smaller number in mag)
    # than the limits
    for i, limit in enumerate(limits):
        flag[:, i] = seds[:, i] < limit

    # Keep index where model is brighter than the limit in N filters
    s = np.sum(flag, axis=1)
    idx, = np.where(s >= Nfilter)

    return idx


def pick_models_toothpick_style(sedgrid_fname, filters, mag_cuts, Nfilter,
                                N_fluxes, min_N_per_flux,
                                outfile=None, bins_outfile=None, mag_pad=.25):
    """
    Creates a fake star catalog from a BEAST model grid. The chosen seds
    are optimized for the toothpick model, by working with a given
    number of flux bins, and making sure that every flux bin is covered
    by at least a given number of models (for each filter individually,
    which is how the toothpick model works).

    Parameters
    ----------
    sedgrid_fname: string
        BEAST model grid from which the models are picked (hdf5 file)

    filters: list of string
        Names of the filters, to be used as columns of the output table

    mag_cuts: list of float
        List of magnitude limits for each filter

    Nfilter: integer
        In how many filters a fake star needs to be brighter than the
        mag_cut value

    N_fluxes: integer
        The number of flux bins into which the dynamic range of the
        model grid in each filter is divided

    min_N_per_flux: integer
        Minimum number of model seds that need to fall into each bin

    outfile: string
        Output path for the models (optional). If this file already
        exists, the chosen seds are loaded from this file instead.

    bins_outfile: string
        Output path for a file containing the flux bin limits for each
        filter, and the number of samples for each (optional)

    mag_pad: float
        The range that the lowest and highest bins should extend above
        and below the minimum and maximum magnitude. Negative values
        will shove models that fall outside of the clipped range into
        the outermost bins.

    Returns
    -------
    sedsMags: astropy Table
        A table containing the selected model seds (columns are named
        after the filters)

    """
    if outfile is not None and os.path.isfile(outfile):
        print('{} already exists. Will attempt to load SEDs for ASTs from there.'.format(
            outfile))
        t = Table.read(outfile, format='ascii')
        return t

    with Vega() as v:
        vega_f, vega_flux, lambd = v.getFlux(filters)

    gridf = h5py.File(sedgrid_fname)

    sedsMags = -2.5 * np.log10(gridf['seds'][:] / vega_flux)
    Nf = sedsMags.shape[1]

    idxs = mag_limits(sedsMags, mag_cuts, Nfilter=Nfilter)
    sedsMags_cut = sedsMags[idxs]

    # Note that i speak of fluxes, but I've recently modified this to
    # work with mags instead

    # Set up a number of flux bins for each filter
    maxes = np.amax(sedsMags_cut, axis=0) + mag_pad
    mins = np.amin(sedsMags_cut, axis=0) - mag_pad

    bin_edges = np.zeros((N_fluxes + 1, Nf))  # indexed on [fluxbin, nfilters]
    for f in range(Nf):
        bin_edges[:, f] = np.linspace(mins[f], maxes[f], N_fluxes + 1)
    bin_mins = bin_edges[:-1, :]
    bin_maxs = bin_edges[1:, :]
    assert(len(bin_mins) == N_fluxes)
    assert(len(bin_maxs) == N_fluxes)

    bin_count = np.zeros((N_fluxes, Nf))
    chosen_idxs = []
    counter = 0
    successes = 0
    include_mask = np.full(idxs.shape, True, dtype=bool)
    chunksize = 100000
    while True:
        counter += 1
        # pick some random models
        rand_idx = np.random.choice(idxs[include_mask], size=chunksize)
        randomseds = sedsMags[rand_idx, :]

        # Find in which bin each model belongs, for each filter
        fluxbins = np.zeros(randomseds.shape, dtype=int)
        for fltr in range(Nf):
            fluxbins[:, fltr] = np.digitize(randomseds[:, fltr],
                                            bin_maxs[:, fltr])

        # Clip in place (models of which the flux is equal to the max
        # are assigned bin nr N_fluxes. Move these down to bin nr
        # N_fluxes - 1)
        np.clip(fluxbins, a_min=0, a_max=N_fluxes - 1, out=fluxbins)

        add_these = np.full((len(rand_idx)), False, dtype=bool)
        for r in range(len(rand_idx)):
            # If any of the flux bins that this model falls into does
            # not have enough samples yet, add it to the list of model
            # spectra to be output
            if (bin_count[fluxbins[r, :], range(Nf)] < min_N_per_flux).any():
                bin_count[fluxbins[r, :], range(Nf)] += 1
                successes += 1
                add_these[r] = True

            # If all these bins are full...
            else:
                # ... do not include this model again, since we will reject it
                # anyway.
                include_mask[idxs == rand_idx] = False

        # Add the approved models
        chosen_idxs.extend(rand_idx[add_these])

        # If some of the randomly picked models were not added
        if not add_these.any():
            # ... check if we have enough samples everywhere, or if all
            # the models have been exhausted (and hence the bins are
            # impossible to fill).
            enough_samples = (bin_count.flatten() >= min_N_per_flux).all()
            still_models_left = include_mask.any()
            if enough_samples or not still_models_left:
                break

        if not counter % 10:
            print('Sampled {} models. {} successfull seds. Ratio = {}'.format(
                counter * chunksize, successes, successes / counter / chunksize))
            print('Bin array:')
            print(bin_count)

    # Gather the selected model seds in a table
    sedsMags = Table(sedsMags[chosen_idxs, :], names=filters)

    if outfile is not None:
        ascii.write(sedsMags, outfile, overwrite=True,
                    formats={k: '%.5f' for k in sedsMags.colnames})

    if bins_outfile is not None:
        bin_info_table = Table()
        col_bigarrays = [bin_mins, bin_maxs, bin_count]
        col_basenames = ['bin_mins_', 'bin_maxs_', 'bin_count_']
        for fltr, filter_name in enumerate(filters):
            for bigarray, basename in zip(col_bigarrays, col_basenames):
                bin_info_table.add_column(
                    Column(bigarray[:, fltr], name=basename + filter_name))
        ascii.write(bin_info_table, bins_outfile, overwrite=True)

    return sedsMags


def pick_models(sedgrid_fname, filters, mag_cuts, Nfilter=3, N_stars=70, Nrealize=20,
                outfile=None):
    """Creates a fake star catalog from a BEAST model grid

    Parameters
    ----------
    sedgrid_fname: string
        BEAST model grid from which the models are picked (hdf5 file)

    filters: list of string
        Names of the filters

    mag_cuts: list
        List of magnitude limits for each filter

    Nfilter: Integer
             In how many filters, you want a fake star to be brighter
             than the limit (mag_cut) (default = 3)

    N_stars: Integer
               Number of stellar models picked per a single log(age)
               (default=70)

    Nrealize: Integer
              Number of realization of each models (default = 20)

    outfile: str
        If a file name is given, the selected models will be written to
        disk

    Returns
    -------
    astropy Table of selected models
    - and optionally -
    ascii file: A list of selected models, written to 'outfile'
    """

    with Vega() as v:               # Get the vega fluxes
        vega_f, vega_flux, lamb = v.getFlux(filters)

    gridf = h5py.File(sedgrid_fname)

    # Convert to Vega mags
    sedsMags = -2.5 * np.log10(gridf['seds'][:] / vega_flux)

    # Select the models above the magnitude limits in N filters
    idxs = mag_limits(sedsMags, mag_cuts, Nfilter=Nfilter)
    grid_cut = gridf['grid'][list(idxs)]

    # Sample the model grid uniformly
    prime_params = np.column_stack(
        (grid_cut['logA'], grid_cut['M_ini'], grid_cut['Av']))
    search_age = np.unique(prime_params[:, 0])

    N_sample = N_stars
    models = []
    for iage in search_age:
        tmp, = np.where(prime_params[:, 0] == iage)
        models.append(np.random.choice(tmp, N_sample))

    index = np.repeat(idxs[np.array(models).reshape((-1))], Nrealize)
    sedsMags = Table(sedsMags[index, :], names=filters)

    if outfile is not None:
        ascii.write(sedsMags, outfile, overwrite=True,
                    formats={k: '%.5f' for k in sedsMags.colnames})

    return sedsMags
