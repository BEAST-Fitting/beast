from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import argparse
import numpy as np
from ..vega import Vega
from ...physicsmodel.grid import FileSEDGrid
from astropy.table import Table
from astropy.table import Column
from astropy.io import ascii
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


def pick_models_per_background(sedgrid, bg_map, N_bg_bins, filters, mag_cuts,
                               Nfilter=3, N_per_age=70, Nrealize=20, outfile=None):
    """
    Create a fake start catalog from a BEAST model grid. Makes sure the
    stars are evenly distributed across regions of different background
    intensity.

    Parameters
    ----------
    sedgrid: beast.grid
        Beast model grid from which the models are picked

    bg_map: str
        Path to a fits file containing a background map. Each row in the
        fits table should represent a tile of the map. The table should
        have columns describing for each tile: the minimum and maximum
        RA, the minimum and maximum DEC,and a value which represents the
        background density.

    N_bg_bins: int
        The number of bins for the range of background density values.
        The bins will be picked on a linear grid, rangin from the
        minimum to the maximum background value of the map. Then, each
        tile will be put in a bin, so that a set of tiles of the map is
        obtained for each range of background density values.

    filters: list of str
        List of filter names used in the catalog/grid

    mag_cuts: list
        List of magnitude limits for each filter

    Nfilter: int
        Model SEDs will be accepted when they are brighter than
        mag_cut in at least Nfilter filters

    N_per_age: int
        Number of stellar models picked for each log(age)

    Nrealize: int
        Number of times all these models are repeated

    Returns
    -------
    astropy Table: List of fake stars, with magnitudes and positions
    - optionally -
    ascii file of this table, written to outfile
    """
    # Get a set of seds, without writing them to file
    fake_star_seds = pick_models(sedgrid, filters, mag_cuts, Nfilter=Nfilter,
                                 N_stars=N_per_age, Nrealize=Nrealize)
    Nseds = len(fake_star_seds)

    # Load the background map
    bg = Table.read(bg_map)
    tile_bg_vals = bg['median_bg']
    min_bg = np.amin(tile_bg_vals)
    max_bg = np.amax(tile_bg_vals)

    # Create the background bins
    # [min, ., ., ., max]
    # 0 [1, 2, 3, 4, 5] 6
    bg_bins = np.linspace(min_bg, max_bg, N_bg_bins + 1)

    # Find which bin each tile belongs to
    bgbin_foreach_tile = np.digitize(tile_bg_vals, bg_bins)
    # Invert this
    tiles_foreach_bgbin = [np.nonzero(bgbin_foreach_tile == b + 1)[0]
                           for b in range(N_bg_bins)]
    print(tiles_foreach_bgbin)

    # Remove empty bins
    tile_sets = [tile_set for tile_set in tiles_foreach_bgbin if len(tile_set)]

    # For each set of tiles, repeat the seds and spread them evenly over
    # the tiles
    repeated_seds = np.repeat(fake_star_seds[:], len(tile_sets))

    out_table = Table(repeated_seds, names=filters)
    ras = np.zeros(len(out_table))
    decs = np.zeros(len(out_table))
    bin_indices = np.zeros(len(out_table))

    tile_ra_min = bg['min_ra']
    tile_dec_min = bg['min_dec']
    tile_ra_delta = bg['max_ra'] - tile_ra_min
    tile_dec_delta = bg['max_dec'] - tile_dec_min

    pbar = Pbar(len(tile_sets), desc='{} models per background bin'.format(Nseds))
    for bin_index, tile_set in pbar.iterover(enumerate(tile_sets)):

        start = bin_index * Nseds
        stop = start + Nseds
        bin_indices[start:stop] = bin_index
        for i in range(Nseds):
            j = bin_index * Nseds + i
            # Pick a random tile
            tile = np.random.choice(tile_set)
            # Within this tile, pick a random ra and dec
            ras[j] = tile_ra_min[tile] + \
                np.random.random_sample() * tile_ra_delta[tile]
            decs[j] = tile_dec_min[tile] + \
                np.random.random_sample() * tile_dec_delta[tile]

    # Add the positions to the table
    ra_col = Column(ras, name='RA')
    out_table.add_column(ra_col)
    dec_col = Column(decs, name='DEC')
    out_table.add_column(dec_col)
    bin_col = Column(bin_indices, name='bg_bin')
    out_table.add_column(bin_col)

    # Write out the table in ascii
    if outfile:
        formats = {k: '%.5f' for k in out_table.colnames}
        ascii.write(out_table, outfile, overwrite=True, formats=formats)

    return out_table


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
    idx = mag_limits(sedsMags, mag_cuts, Nfilter=Nfilter)
    grid_cut = sedgrid.grid[idx]

    # Sample the model grid uniformly
    prime_params = np.column_stack(
        (grid_cut['logA'], grid_cut['M_ini'], grid_cut['Av']))
    search_age = np.unique(prime_params[:, 0])

    N_sample = N_stars
    models = []
    for iage in search_age:
        tmp, = np.where(prime_params[:, 0] == iage)
        models.append(np.random.choice(tmp, N_sample))

    index = np.repeat(idx[np.array(models).reshape((-1))], Nrealize)
    sedsMags = Table(sedsMags[index, :], names=filters)

    if outfile:
        ascii.write(sedsMags, outfile, overwrite=True, formats={
            k: '%.5f' for k in sedsMags.colnames})

    return sedsMags
