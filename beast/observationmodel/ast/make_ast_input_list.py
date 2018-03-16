from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import argparse
import numpy as np
from ..vega import Vega
import datamodel as datamodel
from ...physicsmodel.grid import FileSEDGrid
from astropy.table import Table
from astropy.io import ascii


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
             In how many filters, you want a fake star to be brighter than the limit

    Returns
    -------
    idx:    np.array
            Array of integers contining the indices of allowed models

    """
    flag = seds.copy()

    # flag is True if the models are brigter (=smaller number in mag) than the limits
    for i, limit in enumerate(limits):
        flag[:, i] = seds[:, i] < limit

    # Keep index where model is brighter than the limit in N filters
    s = np.sum(flag, axis=1)
    idx, = np.where(s >= Nfilter)

    return idx


def pick_models_per_background(sedgrid, bg_map, Nbg_bins, N_per_bg, filters, mag_cuts,
                               Nfilter=3, N_per_age=70, Nrealize=20):
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

    bg_bins: int
        The number of bins for the range of background density values.
        The bins will be picked on a linear grid, rangin from the
        minimum to the maximum background value of the map. Then, each
        tile will be put in a bin, so that a set of tiles of the map is
        obtained for each range of background density values.

    filters: list of str
        List of filter names (helps with getting the right magnitudes
        from the grid to apply mag_cuts)

    mag_cuts: list
        List of magnitude limits for each filter

    Nfitler: int
        Model SEDs will be accepted when they are brighter than the
        mag_cut in at least Nfilter filters

    N_per_age: int
        Number of stellar models picked for each log(age)

    Nrealize: int
        Number of times all these models are repeated

    Returns
    -------
    nothing, but writes an ascii file
    """
    # Get the vega flux for each filter
    with Vega() as v:
        vega_f, vega_fluxes, lamb = v.getFlux(filters)

    # Convert model star seds to vega fluxes
    sedsMags = -2.5 * np.log10(sedgrid.seds[:] / vega_fluxes)

    # Select the models that are above the magnitude limits in at least
    # N filters. We will then sample from this grid.
    indices_above_cut = mag_limits(sedsMags, mag_cuts, Nfilter)
    grid_cut = sedgrid.grid[indices_above_cut]

    # Read in the background density map
    bg = Table.read(bg_map)
    tile_bg_vals = bg['median_bg']
    min_bg = np.amin(tile_bg_vals)
    max_bg = np.amax(tile_bg_vals)

    # Create the background bins
    # [min, ., ., ., max]
    # 0 [1, 2, 3, 4, 5] 6
    bg_bins = np.linspace(min_bg, max_bg, Nbg_bins + 1)

    # Find which bin each tile belongs to
    bgbin_foreach_tile = np.digitize(tile_bg_vals, bg_bins)
    # Invert this
    tiles_foreach_bgbin = [np.nonzero(bgbin_foreach_tile == b)[0]
                           for b in range(Nbg_bins)]

    print(tiles_foreach_bgbin)

def pick_models(sedgrid, mag_cuts, Nfilter=3, N_stars=70, Nrealize=20):
    """
    Creates a fake star catalog from a BEAST model grid

    Parameters
    ----------
    sedgrid: beast.grid
               BEAST model grid from which the models are picked

    mag_cuts: list
               List of magnitude limits

    Nfilter: Integer
             In how many filters, you want a fake star to be brighter than the limit
             (default = 3)

    N_stars: Integer
               Number of stellar models picked per a single log(age) (default=70)

    Nrealize: Integer
              Number of realization of each models (default = 20)

    Returns
    -------
    ascii file: A list of selected models
    """

    filters = datamodel.filters

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

    outfile = './' + datamodel.project + '/' + datamodel.project + '_inputAST.txt'

    ascii.write(sedsMags, outfile, overwrite=True, formats={
                k: '%.5f' for k in sedsMags.colnames})
