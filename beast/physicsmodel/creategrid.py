"""
Create extinguished grid more segmented dealing with large grids with enough
memory

All functions are now transformed into generators. As a result, any function
allows computation of a grid in an arbitrary number of chunks. This offers the
possibility to generate grids that cannot fit in memory.


.. note::

    * dependencies have also been updated accordingly.

    * likelihood computations need to be updated to allow computations even if
      the full grid does not fit in memory
"""
from __future__ import (absolute_import, print_function, division)

__version__ = '2.0dev'

import numpy as np
import copy

from astropy import units

from .stars import stellib
from .grid import SpectralGrid
from .prior_weights_dust import PriorWeightsDust
from ..external.eztables import Table
from ..tools.pbar import Pbar
from ..tools.helpers import generator
from ..tools import helpers

from ..observationmodel.noisemodel import absflux_covmat

__all__ = ['gen_spectral_grid_from_stellib_given_points',
           'make_extinguished_grid',
           'add_spectral_properties',
           'calc_absflux_cov_matrices']


@generator
def gen_spectral_grid_from_stellib_given_points(osl, pts,
                                                bounds=dict(dlogT=0.1,
                                                            dlogg=0.3),
                                                chunksize=0):
    """
    Generator that reinterpolates a given stellar spectral library on to
       an Isochrone grid

    It will iterate over a list of `pts` points and generate
       `chunksize` models until all the list of points is processed

    Parameters
    ----------
    osl: stellib.stellib
        a stellar library

    pts: dict like structure of points
        dictionary like or named data structure of points to interpolate at
        must contain logg, logT, logL, and Z

    bounds:  dict, optional (default={dlogT:0.1, dlogg:0.3})
        sensitivity to extrapolation (see grid.get_stellib_boundaries)

    chunksize: int, optional (default=0)
        number of models to generate at each cycle.
        If default <= 0, all models will be returned at once.

    Returns
    -------
    g: SpectralGrid
        Spectral grid (in memory) containing the requested list of stars
        and associated spectra
    """

    helpers.type_checker('osl', osl, stellib.Stellib)

    if chunksize <= 0:
        yield osl.gen_spectral_grid_from_given_points(pts, bounds=bounds)
    else:
        try:
            # Yield successive n-sized chunks from l, assuming we can take
            # slices of the iterator
            for chunk_slice in helpers.chunks(list(range(len(pts))),
                                              chunksize):
                chunk_pts = pts[chunk_slice]
                yield osl.gen_spectral_grid_from_given_points(chunk_pts,
                                                              bounds=bounds)
        except Exception as e:
            # chunks may not work on this as pts is most likely a Table
            print(e)
            for chunk_pts in helpers.chunks(pts, chunksize):
                yield osl.gen_spectral_grid_from_given_points(chunk_pts,
                                                              bounds=bounds)


def _make_dust_fA_valid_points_generator(it, min_Rv, max_Rv):
    """
    compute the allowed points based on the R(V) versus f_A plane
    duplicates effort for all A(V) values, but it is quick compared to
    other steps

    .. note::

        on 2.74: SMC extinction implies f_A = 0. and Rv = 2.74

    Parameters
    ----------
    it: an iterable
        an initial sequence of points that will be trimmed to only valid ones

    min_Rv: float
        lower Rv limit

    max_Rv: float
        upper Rv limit

    Returns
    -------
    npts: int
        the actual number of valid points

    pts: generator
        a generator that only produce valid points
    """
    itn = copy.copy(it)
    npts = 0

    def is_valid(ak, rk, fk):
        return (fk / max_Rv + (1. - fk) / 2.74 <= 1. / rk <= fk * 1. /
                min_Rv + (1. - fk) / 2.74)

    # explore the full list once
    # not very time consuming
    for ak, rk, fk in itn:
        if is_valid(ak, rk, fk):
            npts += 1

    # make the iterator
    pts = ((float(ak), float(rk), float(fk)) for ak, rk, fk in it if
           is_valid(ak, rk, fk))

    return npts, pts


def apply_distance_grid(specgrid, distances,
                        redshift=0):
    """
    Distances are applied to the spectral grid by copying the grid and
    applying a scaling factor.

    Parameters
    ----------

    project: str
        project name

    specgrid: grid.SpectralGrid object
        spectral grid to transform

    distances: list of float
        Distances at which models should be shifted
        0 means absolute magnitude.
        Expecting pc units

    redshift: float
        Redshift to which wavelengths should be shifted
        Default is 0 (rest frame)
    """
    g0 = specgrid

    # Current length of the grid
    N0 = len(g0.grid)
    N = N0 * len(distances)

    # Make singleton list if a single distance is given
    if not hasattr(distances, '__iter__'):
        _distances = [distances]
    else:
        _distances = distances

    # Add distance column if multiple distances are specified
    cols = {}
    cols['distance'] = np.empty(N, dtype=float)

    # Existing columns
    keys0 = list(g0.keys())
    for key in keys0:
        cols[key] = np.empty(N, dtype=float)

    n_sed_points = g0.seds.shape[1]
    new_seds = np.empty((N, n_sed_points), dtype=float)

    for count, distance in \
            Pbar(len(_distances),
                 desc='grid with distances').iterover(enumerate(_distances)):

        # The range where the current distance points will live
        distance_slice = slice(N0 * count, N0 * (count + 1))

        # The seds default to 10 pc.
        # Therefore, scale them with (d / (10 pc))**(-2).
        distance_pc = distance.to(units.pc).value
        new_seds[distance_slice, :] = g0.seds / (0.1 * distance_pc) ** 2

        # Fill in the distance in the distance column
        cols['distance'][distance_slice] = distance_pc

        # Copy the old columns
        for key in keys0:
            cols[key][distance_slice] = g0.grid[key]

    # apply redshift
    g0.lamb = g0.lamb * (1. + redshift)

    # New object
    g = SpectralGrid(g0.lamb, seds=new_seds, grid=Table(cols),
                     backend='memory')
    return g


@generator
def make_extinguished_grid(spec_grid, filter_names, extLaw,
                           avs, rvs, fAs=None,
                           av_prior_model={'name': 'flat'},
                           rv_prior_model={'name': 'flat'},
                           fA_prior_model={'name': 'flat'},
                           chunksize=0,
                           add_spectral_properties_kwargs=None,
                           absflux_cov=False,
                           filterLib=None):
    """
    Extinguish spectra and extract an SEDGrid through given series of filters
    (all wavelengths in stellar SEDs and filter response functions are assumed
    to be in Angstroms)

    Parameters
    ----------
    spec_grid: string or grid.SpectralGrid
        if string:
        spec_grid is the filename to the grid file with stellar spectra
        the backend to load this grid will be the minimal invasive: 'HDF'
        if possible, 'cache' otherwise.

        if not a string, expecting the corresponding SpectralGrid instance
        (backend already setup)

    filter_names: list
        list of filter names according to the filter lib

    Avs: sequence
        Av values to iterate over

    av_prior_model: list
        list including prior model name and parameters

    Rvs: sequence
        Rv values to iterate over

    rv_prior_model: list
        list including prior model name and parameters

    fAs: sequence (optional)
        f_A values to iterate over
        f_A can be omitted if the extinction Law does not use it or allow
        fixed values

    fA_prior_model: list
        list including prior model name and parameters

    chunksize: int, optional (default=0)
        number of extinction model variations to generate at each cycle.
        Note that this means len(spec_grid * chunksize)
        If default <= 0, all models will be returned at once.

    filterLib:  str
        full filename to the filter library hd5 file

    add_spectral_properties_kwargs: dict
        keyword arguments to call :func:`add_spectral_properties` at each
        iteration to add model properties from the spectra into the grid
        property table

    asbflux_cov: boolean
        set to calculate the absflux covariance matrices for each model
        (can be very slow!!!  But it is the right thing to do)

    Returns
    -------
    g: grid.SpectralGrid
        final grid of reddened SEDs and models
    """
    # Check inputs
    # ============
    # get the stellar grid (no dust yet)
    # if string is provided try to load the most memory efficient backend
    # otherwise use a cache-type backend (load only when needed)
    if type(spec_grid) == str:
        ext = spec_grid.split('.')[-1]
        if ext in ['hdf', 'hd5', 'hdf5']:
            g0 = SpectralGrid(spec_grid, backend='hdf')
        else:
            g0 = SpectralGrid(spec_grid, backend='cache')
    else:
        helpers.type_checker('spec_grid', spec_grid, SpectralGrid)
        g0 = spec_grid

    # Tag fA usage
    if fAs is None:
        with_fA = False
    else:
        with_fA = True

    # get the min/max R(V) values necessary for the grid point definition
    min_Rv = min(rvs)
    max_Rv = max(rvs)

    # Create the sampling mesh
    # ========================
    # basically the dot product from all input 1d vectors
    # setup interation over the full dust parameter grid
    if with_fA:
        dustpriors = PriorWeightsDust(avs, av_prior_model,
                                      rvs, rv_prior_model,
                                      fAs, fA_prior_model)

        it = np.nditer(np.ix_(avs, rvs, fAs))
        niter = np.size(avs) * np.size(rvs) * np.size(fAs)
        npts, pts = _make_dust_fA_valid_points_generator(it, min_Rv, max_Rv)

        # Pet the user
        print("""number of initially requested points = {0:d}
              number of valid points = {1:d} (based on restrictions in R(V)
                 versus f_A plane)
              """.format(niter, npts))

        if npts == 0:
            raise AttributeError('No valid points')
    else:
        dustpriors = PriorWeightsDust(avs, av_prior_model,
                                      rvs, rv_prior_model,
                                      [1.0], fA_prior_model)

        it = np.nditer(np.ix_(avs, rvs))
        npts = np.size(avs) * np.size(rvs)
        pts = ((float(ak), float(rk)) for ak, rk in it)

    # Generate the Grid
    # =================
    N0 = len(g0.grid)
    N = N0 * npts

    if chunksize <= 0:
        print('Generating a final grid of {0:d} points'.format(N))
    else:
        print('Generating a final grid of {0:d} points in {1:d}' +
              ' pieces'.format(N, int(float(N0) / chunksize + 1.)))

    if chunksize <= 0:
        chunksize = npts

    if add_spectral_properties_kwargs is not None:
        nameformat = add_spectral_properties_kwargs.pop('nameformat',
                                                        '{0:s}') + '_wd'

    for chunk_pts in helpers.chunks(pts, chunksize):
        # iter over chunks of models

        # setup chunk outputs
        cols = {'Av': np.empty(N, dtype=float),
                'Rv': np.empty(N, dtype=float)
                }

        if with_fA:
            cols['Rv_A'] = np.empty(N, dtype=float)
            cols['f_A'] = np.empty(N, dtype=float)

        keys = list(g0.keys())
        for key in keys:
            cols[key] = np.empty(N, dtype=float)

        n_filters = len(filter_names)
        _seds = np.empty((N, n_filters), dtype=float)
        if absflux_cov:
            n_offdiag = (((n_filters**2)-n_filters)/2)
            _cov_diag = np.empty((N, n_filters), dtype=float)
            _cov_offdiag = np.empty((N, n_offdiag), dtype=float)

        for count, pt in \
                Pbar(npts, desc='SED grid').iterover(enumerate(chunk_pts)):

            if with_fA:
                Av, Rv, f_A = pt
                dust_prior_weight = dustpriors.get_weight(Av, Rv, f_A)
                Rv_MW = extLaw.get_Rv_A(Rv, f_A)
                r = g0.applyExtinctionLaw(extLaw, Av=Av, Rv=Rv, f_A=f_A,
                                          inplace=False)
                # add extra "spectral bands" if requested
                if add_spectral_properties_kwargs is not None:
                    r = add_spectral_properties(
                        r, nameformat=nameformat,
                        filterLib=filterLib,
                        **add_spectral_properties_kwargs)
                temp_results = r.getSEDs(filter_names,
                                         filterLib=filterLib)
                # adding the dust parameters to the models
                cols['Av'][N0 * count: N0 * (count + 1)] = Av
                cols['Rv'][N0 * count: N0 * (count + 1)] = Rv
                cols['f_A'][N0 * count:N0 * (count + 1)] = f_A
                cols['Rv_A'][N0 * count: N0 * (count + 1)] = Rv_MW

            else:
                Av, Rv = pt
                dust_prior_weight = dustpriors.get_weight(Av, Rv, 1.0)
                r = g0.applyExtinctionLaw(extLaw, Av=Av, Rv=Rv, inplace=False)

                if add_spectral_properties_kwargs is not None:
                    r = add_spectral_properties(
                        r, nameformat=nameformat,
                        filterLib=filterLib,
                        **add_spectral_properties_kwargs)
                temp_results = r.getSEDs(filter_names,
                                         filterLib=filterLib)
                # adding the dust parameters to the models
                cols['Av'][N0 * count: N0 * (count + 1)] = Av
                cols['Rv'][N0 * count: N0 * (count + 1)] = Rv

            # get new attributes if exist
            for key in list(temp_results.grid.keys()):
                if key not in keys:
                    k1 = N0 * count
                    k2 = N0 * (count + 1)
                    cols.setdefault(key, np.empty(N, dtype=float))[k1:k2] = \
                        temp_results.grid[key]

            # compute the fractional absflux covariance matrices
            if absflux_cov:
                absflux_covmats = calc_absflux_cov_matrices(r, temp_results,
                                                            filter_names)
                _cov_diag[N0 * count: N0 * (count + 1)] = absflux_covmats[0]
                _cov_offdiag[N0 * count: N0 * (count + 1)] = absflux_covmats[1]

            # assign the extinguished SEDs to the output object
            _seds[N0 * count: N0 * (count + 1)] = temp_results.seds[:]

            # copy the rest of the parameters
            for key in keys:
                cols[key][N0 * count: N0 * (count + 1)] = g0.grid[key]

            # multiply existing prior weights by the dust prior weight
            cols['weight'][N0 * count: N0 * (count + 1)] \
                *= dust_prior_weight
            cols['prior_weight'][N0 * count: N0 * (count + 1)] \
                *= dust_prior_weight

            if count == 0:
                cols['lamb'] = temp_results.lamb[:]

        _lamb = cols.pop('lamb')

        # free the memory of temp_results
        # del temp_results
        # del tempgrid

        # Ship
        if absflux_cov:
            g = SpectralGrid(_lamb, seds=_seds,
                             cov_diag=_cov_diag, cov_offdiag=_cov_offdiag,
                             grid=Table(cols), backend='memory')
        else:
            g = SpectralGrid(_lamb, seds=_seds,
                             grid=Table(cols), backend='memory')

        g.grid.header['filters'] = ' '.join(filter_names)

        yield g


def add_spectral_properties(specgrid, filternames=None, filters=None,
                            callables=None, nameformat=None, filterLib=None):
    """ Addon spectral calculations to spectral grids to extract in the fitting
    routines

    Parameters
    ----------
    specgrid: SpectralGrid instance
        instance of the spectral grid

    filternames: sequence(str)
        compute the integrated values through given filters in the library

    filters: sequence(Filters)
        sequence of filter instances from which extract integrated values

    callables: sequence(callable)
        sequence of functions to apply onto the spectral grid assuming storing
        results is internally processed by the individual functions

    nameformat: str
        naming format to adopt for filternames and filters
        default value is '{0:s}_0' where the value will be the filter name

    filterLib:  str
        full filename to the filter library hd5 file

    Returns
    -------
    specgrid: SpectralGrid instance
        instance of the input spectral grid which will include more properties
    """
    if nameformat is None:
        nameformat = '{0:s}_0'

    if filternames is not None:
        temp = specgrid.getSEDs(filternames, extLaw=None, filterLib=filterLib)

        logtempseds = np.array(temp.seds)
        indxs = np.where(temp.seds > 0)
        if len(indxs) > 0:
            logtempseds[indxs] = np.log10(temp.seds[indxs])
        indxs = np.where(temp.seds <= 0)
        if len(indxs) > 0:
            logtempseds[indxs] = -100.

        for i, fk in enumerate(filternames):
            specgrid.grid.addCol('log'+nameformat.format(fk),
                                 logtempseds[:, i])
        del temp

    if filters is not None:
        temp = specgrid.getSEDs(filters, extLaw=None)

        logtempseds = np.array(temp.seds)
        indxs = np.where(temp.seds > 0)
        if len(indxs) > 0:
            logtempseds[indxs] = np.log10(temp.seds[indxs])

        indxs = np.where(temp.seds <= 0)
        if len(indxs) > 0:
            logtempseds[indxs] = -100.

        for i, fk in enumerate(filters):
            specgrid.grid.addCol('log'+nameformat.format(fk.name),
                                 logtempseds[:, i])
        del temp

    if callables is not None:
        for fn in callables:
            fn(specgrid)

    return specgrid


def calc_absflux_cov_matrices(specgrid, sedgrid, filter_names):
    """ Calculate the absflux covariance matrices for each model
    Must be done on the full spectrum of each model to account for
    the changing combined spectral response due to the model SED and
    the filter response curve.

    Parameters
    ----------
    specgrid: SpectralGrid instance
        instance of the spectral grid containing the full spectrum fluxes

    sedgrid: SpectralGrid instance
        instance of the spectral grid containing the band SED fluxes

    Returns
    -------
    absflux_covmat :
    """

    # get the fractional absflux covariance matrix
    absflux_cov_mats = absflux_covmat.hst_frac_matrix(
        filter_names, spectrum=(specgrid.lamb[:], specgrid.seds))

    # setup the output quantities
    n_models = specgrid.seds.shape[0]
    n_filters = len(filter_names)
    n_offdiag = (((n_filters**2)-n_filters)/2)
    cov_diag = np.empty((n_models, n_filters), dtype=np.float64)
    cov_offdiag = np.empty((n_models, n_offdiag), dtype=np.float64)

    # pack the resulting covariance matrices into diganonal and
    # non-diagnonal terms
    #   much more efficient for use later in combining with AST results
    #     and fitting
    #   also convert from fractional to physical flux units
    m = 0
    cov_diag[:, n_filters-1] = (absflux_cov_mats[:, n_filters-1, n_filters-1]
                                * np.square(sedgrid.seds[:, n_filters-1]))
    for k in range(n_filters-1):
        cov_diag[:, k] = (absflux_cov_mats[:, k, k]
                          * np.square(sedgrid.seds[:, k]))
        for l in range(k+1, n_filters):
            cov_offdiag[:, m] = (absflux_cov_mats[:, k, l]
                                 * sedgrid.seds[:, k]*sedgrid.seds[:, l])
            m += 1

    return (cov_diag, cov_offdiag)
