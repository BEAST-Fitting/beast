from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse
import sys

import numpy as np

from astropy.table import Table
from astropy.table import Column
from astropy.io import ascii, fits
from astropy.wcs import WCS

from ..vega import Vega
from ...physicsmodel.grid import FileSEDGrid
from ...tools.pbar import Pbar


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


def pick_models_toothpick_style(sedgrid, filters, mag_cuts, Nfilter,
                                N_fluxes, min_N_per_flux, Nrealize,
                                outfile=None, bins_outfile=None, mag_pad=.25):
    with Vega() as v:
        vega_f, vega_flux, lambd = v.getFlux(filters)

    sedsMags = -2.5 * np.log10(sedgrid.seds[:] / vega_flux)
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
    while True:
        counter += 1
        # pick a random model
        rand_idx = np.random.choice(idxs[include_mask])

        # find which flux bin it belongs to for each filter
        fluxbins = [np.digitize(flux, bin_maxs[:, fltr])
            for fltr, flux in enumerate(sedsMags[rand_idx, :])]

        # If any of the flux bins that this model falls into does not
        # have enough samples yet, add it to the list of model spectra
        # to be output
        if (bin_count[fluxbins, range(Nf)] < min_N_per_flux).any():
            bin_count[fluxbins, range(Nf)] += 1
            successes += 1
            chosen_idxs.append(rand_idx)

        # If all these bins are full...
        else:
            # ... do not include this model again, since we will reject it
            # anyway.
            include_mask[idxs == rand_idx] = False
            # ... check if we have enough samples everywhere, or if all
            # the models have been exhausted (and hence the bins are
            # impossible to fill).
            enough_samples = (bin_count.flatten() >= min_N_per_flux).all()
            still_models_left = include_mask.any()
            if enough_samples or not still_models_left:
                break

        if not counter % 10000:
            print('Sampled {} models. {} successfull. Ratio = {}'.format(
                counter, successes, successes / counter))

    # Gather the selected model seds in a table
    sedsMags = Table(sedsMags[chosen_idxs, :], names=filters)

    if outfile is not None:
        ascii.write(sedsMags, outfile, overwrite=True,
                    formats={k: '%.5f' for k in sedsMags.colnames})

    if bins_outfile is not None:
        bin_info_table = Table()
        col_bigarrays = [bin_mins, bin_maxs, bin_count]
        col_basenames = ['bin_mins_', 'bin_maxs_', 'bin_count_']
        all_cols = []
        for fltr, filter_name in enumerate(filters):
            for bigarray, basename in zip(col_bigarrays, col_basenames):
                bin_info_table.add_column(
                    Column(bigarray[:, fltr], name=basename + filter_name))
        ascii.write(bin_info_table, bins_outfile, overwrite=True)

    return sedsMags


def pick_models(sedgrid, filters, mag_cuts, Nfilter=3, N_stars=70, Nrealize=20,
                outfile=None):
    """Creates a fake star catalog from a BEAST model grid

    Parameters
    ----------
    sedgrid: beast.grid
               BEAST model grid from which the models are picked

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

    # Convert to Vega mags
    sedsMags = -2.5 * np.log10(sedgrid.seds[:] / vega_flux)

    # Select the models above the magnitude limits in N filters
    idxs = mag_limits(sedsMags, mag_cuts, Nfilter=Nfilter)
    grid_cut = sedgrid.grid[idxs]

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
