"""
Create extinguished grid more segmented dealing with large grids with enough memory

All functions are now transformed into generators. As a result, any function
allows computation of a grid in an arbitrary number of chunks. This offers the
possibility to generate grids that cannot fit in memory.


.. note::

    * dependencies have also been updated accordingly.

    * likelihood computations need to be updated to allow computations even if
      the full grid does not fit in memory
"""

__version__ = '2.0dev'

import numpy as np
import copy

from . import stellib
from . import isochrone
from . import extinction
from .grid import SpectralGrid
from ..external import ezunits
from ..external.eztables import Table
from ..tools.pbar import Pbar
from ..tools.helpers import generator
from ..tools import helpers
#from .priors import KDTreeDensityEstimator


@generator
def gen_spectral_grid_from_stellib_given_points(osl, pts, bounds=dict(dlogT=0.1, dlogg=0.3), chunksize=0):
    """ generator that reinterpolates a given stellar spectral library on to an Isochrone grid
    it will iterate over a list of `pts` points and generate `chunksize` models
    until all the list of points is processed.

    Parameters
    ----------
    osl: stellib.stellib
        a stellar library

    pts: dict like structure of points
        dictionary like or named data structure of points to interpolate at.
        must contain:
            * logg  surface gravity in log-scale
            * logT  log of effective temperatures (in Kelvins)
            * logL  log of luminosity in Lsun units
            * Z     metallicity

    bounds:  dict, optional (default={dlogT:0.1, dlogg:0.3})
        sensitivity to extrapolation (see grid.get_stellib_boundaries)

    chunksize: int, optional (default=0)
        number of models to generate at each cycle.
        If default <= 0, all models will be returned at once.

    Returns
    -------
    g: SpectralGrid
        Spectral grid (in memory) containing the requested list of stars and associated spectra
    """

    helpers.type_checker('osl', osl, stellib.Stellib)

    if chunksize <= 0:
        yield osl.gen_spectral_grid_from_given_points(pts, bounds=bounds)
    else:
        try:
            # Yield successive n-sized chunks from l, assuming we can take
            # slices of the iterator
            for chunk_slice in helpers.chunks(range(len(pts)), chunksize):
                chunk_pts = pts[chunk_slice]
                yield osl.gen_spectral_grid_from_given_points(chunk_pts, bounds=bounds)
        except Exception as e:
            #chunks may not work on this as pts is most likely a Table
            print(e)
            for chunk_pts in helpers.chunks(pts, chunksize):
                yield osl.gen_spectral_grid_from_given_points(chunk_pts, bounds=bounds)


@helpers.deprecated
def gen_spectral_grid_from_stellib(osl, oiso, ages=(1e7,), masses=(3,), Z=(0.02,),
                                   bounds=dict(dlogT=0.1, dlogg=0.3)):
    """ Reinterpolate a given stellar spectral library on to an Isochrone grid

    DEPRECATED: use gen_spectral_grid_from_stellib_given_points instead

    Parameters
    ----------
    osl: stellib.stellib
        a stellar library

    oiso: isochrone.Isochrone
        an isochrone library

    ages: iterable
        list of age points to include in the grid   (in Yr)

    masses: iterable
        list of mass points to include in the grid  (M/Msun)

    ages: iterable
        list of metallicity points to include in the grid (Z/Zsun)

    bounds: dict
        sensitivity to extrapolation (see get_stellib_boundaries)

    Returns
    -------
    g: SpectralGrid
        Spectral grid (in memory) containing the requested list of stars and associated spectra
    """
    helpers.type_checker('osl', osl, stellib.Stellib)
    helpers.type_checker('oiso', oiso, isochrone.Isochrone)

    # Step 0: generate interpolation points and prepare outputs
    # =========================================================

    # Grid points
    # -----------
    # Points are provided by age, mass, Z
    # We need to iterate over isochrones defined in (age, Z) space
    # and extract relevant information for specific masses

    _masses = np.asarray(masses)

    # make an iterator containing every point of the model grid
    # only one point is loaded at a time and not the full list
    niter = len(ages) * len(Z)
    it = np.nditer(np.ix_(ages, Z))

    # total number of points
    ndata = niter * len(masses)

    # prepare outputs
    # ---------------
    # Grid properties will be stored into a dictionary format until saved on disk
    # SEDs are kept into a ndarray

    _grid  = {}
    for k in oiso.data.keys():
        _grid[k] = np.empty(ndata, dtype=float )

    _grid['radius'] = np.empty(ndata, dtype=float )
    _grid['keep'] = np.empty(ndata, dtype=bool )

    lamb = osl.wavelength[:]
    specs = np.empty( (ndata, len(lamb)), dtype=float )

    # Loop over (age, Z) space
    # ========================

    # some constants
    kdata = 0
    rsun = ezunits.unit['Rsun'].to('m').magnitude  # 6.955e8 m

    for k, (_ak, _Zk) in Pbar(niter, desc='spectral grid').iterover(enumerate(it)):

        # Step 1: get isochrone points
        # ============================
        # get the isochrone of (age, Z) sampled at given masses
        r = oiso._get_isochrone(_ak, metal=_Zk, masses=np.log10(_masses))

        # keep array pointer mark
        start_idx = k * len(_masses)
        end_idx   = start_idx + len(r)

        # Step 2: Avoid Extrapolation
        # ===========================
        # check boundary conditions, keep the data but do not compute the sed if not needed
        bound_cond = osl.points_inside(zip(r['logg'], r['logT']))
        _grid['keep'][start_idx: end_idx] = bound_cond[:]

        # Step 3: radii
        # =============
        # Stellar library models are given in cm^-2  ( 4 pi R)
        # Compute radii of each point using log(T) and log(L)
        # get the isochrone of (age, Z) sampled at given masses
        radii = osl.get_radius(r['logL'], r['logT'])
        weights = 4. * np.pi * (radii * 1e2) ** 2  # denorm models are in cm**-2 (4*pi*rad)
        _grid['radius'][start_idx: end_idx] = radii / rsun

        # Step 4: Interpolation
        # =====================
        # Do the actual interpolation, avoiding exptrapolations
        for mk in range(r.nrows):
            if bound_cond[mk]:
                s = np.array( osl.interp(r['logT'][mk], r['logg'][mk], _Zk, 0.) ).T
                specs[kdata, :] = osl.genSpectrum(s) * weights[mk]
            else:
                specs[kdata, :] = np.zeros(len(lamb), dtype=float )
            kdata += 1

        # Step 4: Store properties
        # ========================
        for key in r.keys():
            _grid[key][start_idx: end_idx] = r[key]

    # Step 5: filter points without spectrum
    # ======================================
    #filter unbound values
    idx = np.array(_grid.pop('keep'))

    specs = specs.compress(idx, axis=0)
    for k in _grid.keys():
            _grid[k] = _grid[k].compress(idx, axis=0)

    # Step 5: Ship
    # ============
    header = {'stellib': osl.source,
              'isoch': oiso.source,
              'comment': 'radius in Rsun',
              'name': 'Reinterpolated stellib grid'}

    g = SpectralGrid(lamb, seds=specs, grid=Table(_grid), header=header, backend='memory')

    return g


def _make_dust_fbump_valid_points_generator(it, min_Rv, max_Rv):
    """
    compute the allowed points based on the R(V) versus f_bump plane
    duplicates effort for all A(V) values, but it is quick compared to
    other steps

    .. note::

        on 2.74: Bumpless extinction implies f_bump = 0. and Rv = 2.74

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
        return (fk / max_Rv + (1. - fk) / 2.74 <= 1. / rk <= fk * 1. / min_Rv + (1. - fk) / 2.74)

    # explore the full list once
    # not very time consuming
    for ak, rk, fk in itn:
        if is_valid(ak, rk, fk):
            npts += 1

    #make the iterator
    pts = ( (float(ak), float(rk), float(fk)) for ak, rk, fk in it if is_valid(ak, rk, fk) )

    return npts, pts


@generator
def make_extinguished_grid(spec_grid, filter_names, extLaw, avs, rvs,
                           fbumps=None, chunksize=0,
                           add_spectral_properties_kwargs=None):
    """
    Extinguish spectra and extract an SEDGrid through given series of filters
    (all wavelengths in stellar SEDs and filter response functions are assumed
    to be in Angstroms)

    keywords
    --------

    spec_grid: string or grid.SpectralGrid
        if string:
            spec_grid is the filename to the grid file with stellar spectra
            the backend to load this grid will be the minimal invasive: 'HDF' if
            possible, 'cache' otherwise.

        if not a string, expecting the corresponding SpectralGrid instance (backend already setup)

    filter_names: list
        list of filter names according to the filter lib

    Avs: sequence
        Av values to iterate over

    Rvs: sequence
        Rv values to iterate over

    fbumps: sequence (optional)
        f_bump values to iterate over
        f_bump can be omitted if the extinction Law does not use it or allow fixed values

    chunksize: int, optional (default=0)
        number of extinction model variations to generate at each cycle.
        Note that this means len(spec_grid * chunksize)
        If default <= 0, all models will be returned at once.

    add_spectral_properties_kwargs: dict
        keyword arguments to call :func:`add_spectral_properties` at each
        iteration to add model properties from the spectra into the grid property table

    returns
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

    # Tag fbump usage
    if fbumps is None:
        with_fb = False
    else:
        with_fb = True

    # get the min/max R(V) values necessary for the grid point definition
    min_Rv = min(rvs)
    max_Rv = max(rvs)

    # Create the sampling mesh
    # ========================
    # basically the dot product from all input 1d vectors
    # setup interation over the full dust parameter grid
    if with_fb:
        it = np.nditer(np.ix_(avs, rvs, fbumps))
        niter = np.size(avs) * np.size(rvs) * np.size(fbumps)
        npts, pts = _make_dust_fbump_valid_points_generator(it, min_Rv, max_Rv)

        # Pet the user
        print("""number of initially requested points = {0:d}
              number of valid points = {1:d} (based on restrictions in R(V) versus f_bump plane)
              """.format(niter, npts))

        if npts == 0:
            raise AttributeError('No valid points')
    else:
        it = np.nditer(np.ix_(avs, rvs))
        npts = np.size(avs) * np.size(rvs)
        pts = ( (float(ak), float(rk)) for ak, rk in it)

    # Generate the Grid
    # =================
    N0 = len(g0.grid)
    N = N0 * npts

    if chunksize <= 0:
        print('Generating a final grid of {0:d} points'.format(N))
    else:
        print('Generating a final grid of {0:d} points in {1:d} pieces'.format(N, int(float(N0) / chunksize + 1.)))

    if chunksize <= 0:
        chunksize = npts

    if add_spectral_properties_kwargs is not None:
        nameformat = add_spectral_properties_kwargs.pop('nameformat', '{0:s}') + '_1'

    for chunk_pts in helpers.chunks(pts, chunksize):
        # iter over chunks of models

        # setup chunk outputs
        cols = {'Av': np.empty(N, dtype=float),
                'Rv': np.empty(N, dtype=float)
                }

        if with_fb:
            cols['Rv_MW'] = np.empty(N, dtype=float)
            cols['f_bump'] = np.empty(N, dtype=float)

        keys = g0.keys()
        for key in keys:
            cols[key] = np.empty(N, dtype=float)

        _seds = np.empty( (N, len(filter_names)), dtype=float)

        for count, pt in Pbar(npts, desc='SED grid').iterover(enumerate(chunk_pts)):

            if with_fb:
                Av, Rv, f_bump = pt
                Rv_MW = extLaw.get_Rv_A(Rv, f_bump)
                r = g0.applyExtinctionLaw(extLaw, Av=Av, Rv=Rv, f_bump=f_bump, inplace=False)
                if add_spectral_properties_kwargs is not None:
                    r = add_spectral_properties(r, nameformat=nameformat, **add_spectral_properties_kwargs)
                temp_results = r.getSEDs(filter_names)
                # adding the dust parameters to the models
                cols['Av'][N0 * count: N0 * (count + 1)] = Av
                cols['Rv'][N0 * count: N0 * (count + 1)] = Rv
                cols['f_bump'][N0 * count:N0 * (count + 1)] = f_bump
                cols['Rv_MW'][N0 * count: N0 * (count + 1)] = Rv_MW
            else:
                Av, Rv = pt
                r = g0.applyExtinctionLaw(extLaw, Av=Av, Rv=Rv, inplace=False)
                if add_spectral_properties_kwargs is not None:
                    r = add_spectral_properties(r, nameformat=nameformat, **add_spectral_properties_kwargs)
                temp_results = r.getSEDs(filter_names)
                # adding the dust parameters to the models
                cols['Av'][N0 * count: N0 * (count + 1)] = Av
                cols['Rv'][N0 * count: N0 * (count + 1)] = Rv

            # get new attributes if exist
            for key in temp_results.grid.keys():
                if key not in keys:
                    cols.setdefault(key, np.empty(N, dtype=float))[N0 * count: N0 * (count + 1)] = temp_results.grid[key]

            # assign the extinguished SEDs to the output object
            _seds[N0 * count: N0 * (count + 1)] = temp_results.seds[:]

            # copy the rest of the parameters
            for key in keys:
                cols[key][N0 * count: N0 * (count + 1)] = g0.grid[key]

            if count == 0:
                cols['lamb'] = temp_results.lamb[:]

        #Adding Density
        #tempgrid = np.array([ cols[k] for k in 'logA M_ini M_act Av Rv f_bump Z'.split() if k in cols ]).T
        #tr = KDTreeDensityEstimator(tempgrid)
        #cols['Density'] = tr(tempgrid)
        _lamb = cols.pop('lamb')

        # free the memory of temp_results
        #del temp_results
        #del tempgrid

        # Ship
        g = SpectralGrid(_lamb, seds=_seds, grid=Table(cols), backend='memory')
        g.grid.header['filters'] = ' '.join(filter_names)
        yield g


def add_spectral_properties(specgrid, filternames=None, filters=None, callables=None, nameformat=None):
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

    Returns
    -------
    specgrid: SpectralGrid instance
        instance of the input spectral grid which will include more properties
    """
    if nameformat is None:
        nameformat = '{0:s}_0'

    if filternames is not None:
        temp = specgrid.getSEDs(filternames, extLaw=None)
        for i, fk in enumerate(filternames):
            specgrid.grid.addCol(nameformat.format(fk), temp.seds[:, i])
        del temp

    if filters is not None:
        temp = specgrid.getSEDs(filters, extLaw=None)
        for i, fk in enumerate(filters):
            specgrid.grid.addCol(nameformat.format(fk.name), temp.seds[:, i])
        del temp

    if callables is not None:
        for fn in callables:
            fn(specgrid)

    return specgrid


#=================== TESTUNITS ============================

def test_gen_spectral_grid_from_stellib_given_points():
    """test_gen_spectral_grid_from_stellib_given_points
    Make sure it runs and returns a grid
    """
    osl = stellib.Kurucz()
    oiso = isochrone.padova2010()
    chunksize = 10000

    #as it is an interator, list does the actual loop
    list(gen_spectral_grid_from_stellib_given_points(osl, oiso.data, chunksize=chunksize))


@helpers.deprecated
def test_gen_spectral_grid_from_stellib():
    """test_gen_spectral_grid_from_stellib
    Make sure it runs and returns a grid
    """
    osl = stellib.Kurucz()
    oiso = isochrone.padova2010()
    gen_spectral_grid_from_stellib(osl, oiso, ages=(1e7, 1e8), masses = (3., 4.), Z=(0.02, 0.004))


def test_make_extinguished_grid():
    """Make a grid from isochrone points and run the extinction to extract SEDs
    A spectral grid will be generated (but not saved), using the stellar parameters
    by directly using points from the isochrones.
    This spectral grid will be then reddened according to the Mixture models and photometry will be extracted.
    The final grid is saved on disk.
    """

    #select the stellar library and the isochrones to use
    osl = stellib.Kurucz()
    oiso = isochrone.PadovaWeb()
    oiso.data = oiso.get_t_isochrones(6.5, 8.0, 0.1, 0.02)
    oext = extinction.RvFbumpLaw()
    chunksize = 1e3

    # a spectral grid will be generated but not saved
    # using the stellar parameters by interpolation of the isochrones and the generation of spectra into the physical units
    grid_fname = 'test_kurucz2004.spectral.grid.hd5'

    # define filters for the grid
    filter_names  = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    # variable to ensure that range is fully covered in using np.arange
    __tiny_delta__ = 0.001
    # grid spacing for dust
    avs            = np.arange(0.0, 5.0 + __tiny_delta__, 0.5)
    rvs            = np.arange(1.0, 6.0 + __tiny_delta__, 1.0)
    fbumps         = np.asarray([1.0])

    #make the spectral grid
    g = gen_spectral_grid_from_stellib_given_points(osl, oiso.data, chunksize=chunksize)

    part = 0
    for gk in g:
        # make the grid
        gen_extgrid = make_extinguished_grid(gk, filter_names, oext, avs, rvs, fbumps, chunksize=chunksize)

        print('Calling make_extinguished_grid')
        for extgrid in gen_extgrid:
            part += 1
            # save grid to file
            _fname = grid_fname.split('.')
            _fname = '.'.join(_fname[:-1]) + 'part{0:d}'.format(part) + _fname[-1]
            print('Saving part {0:d}'.format(part))
            extgrid.writeHDF(grid_fname, append=True)
