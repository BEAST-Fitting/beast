"""
Create extinguished grid

Dec 2012: Originally written by Kirill T.
Jan 2013: Modified by Karl G. to separate the extinguished grid creation into a separate program.
          Added restriction of parameter space in R(V) and f_bump.

Major rewritting progress:
    functions return a grid object (MemoryBackend) that can be saved or not
    cleaning functions and copy and better documentation
    memory efficiency as a priority

TODO: as clearly visible with this version, most of this goes into stellib class unless recursive imports:
    especially:  * gen_spectral_grid_from_stellib_given_points
"""

__version__ = '0.8dev'

import numpy as np

from . import stellib
from . import isochrone
from .grid import SpectralGrid
from .gridhelpers import isNestedInstance
from ..external import ezunits
from ..external.eztables import Table
from ..tools import progressbar


def gen_spectral_grid_from_stellib_given_points(osl, pts, bounds=dict(dlogT=0.1, dlogg=0.3)):
    """ Reinterpolate a given stellar spectral library on to an Isochrone grid

    keywords
    --------
    osl: stellib.stellib
        a stellar library

    pts: dict like structure of points
        dictionary like or named data structure of points to interpolate at.
        must contain:
            * logg  surface gravity in log-scale
            * logT  log of effective temperatures (in Kelvins)
            * logL  log of luminosity in Lsun units
            * Z     metallicity

    bounds:  dict
        sensitivity to extrapolation (see grid.get_stellib_boundaries)
        default: {dlogT:0.1, dlogg:0.3}

    Returns
    -------
    g: SpectralGrid
        Spectral grid (in memory) containing the requested list of stars and associated spectra
    """
    assert(isNestedInstance(osl, stellib.Stellib) )
    return osl.gen_spectral_grid_from_given_points(pts, bounds=bounds)


def gen_spectral_grid_from_stellib(osl, oiso, ages=(1e7,), masses=(3,), Z=(0.02,), bounds=dict(dlogT=0.1, dlogg=0.3)):
    """ Reinterpolate a given stellar spectral library on to an Isochrone grid
    keywords
    --------
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

    returns
    -------
    g: SpectralGrid
        Spectral grid (in memory) containing the requested list of stars and associated spectra
    """
    assert(isNestedInstance(osl, stellib.Stellib) )
    assert(isNestedInstance(oiso, isochrone.Isochrone) )

    # Step 0: generate interpolation points and prepare outputs
    # =========================================================

    # Grid points
    # -----------
    # Points are provided by age, mass, Z
    # We need to iterate over isochrones defined in (age, Z) space
    # and extract relevant information for specific masses

    _masses = np.asarray(masses)

    # make an iterator
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

    with progressbar.PBar(niter, txt='spectral grid') as Pbar:
        for k, (_ak, _Zk) in enumerate(it):
            Pbar.update(k)

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


def make_extinguished_grid(spec_grid, filter_names, extLaw, avs, rvs, fbumps=None):
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
        g0 = spec_grid

    N0 = len(g0.grid)

    # Tag fbump usage
    if fbumps is None:
        with_fb = False
    else:
        with_fb = True

    # get the min/max R(V) values
    min_Rv = min(rvs)
    max_Rv = max(rvs)

    # Create the sampling mesh
    # ========================
    # basically the dot product from all input 1d vectors
    # setup interation over the full dust parameter grid
    if with_fb:
        it = np.nditer(np.ix_(avs, rvs, fbumps))
        niter = np.size(avs) * np.size(rvs) * np.size(fbumps)

        # compute the allowed points based on the R(V) versus f_bump plane
        # duplicates effort for all A(V) values, but it is quick compared to
        # other steps
        # Note on 2.74: Bumpless extinction implies f_bump = 0. and Rv = 2.74
        pts = [ (float(ak), float(rk), float(fk))
                for ak, rk, fk in it
                if fk / max_Rv + (1. - fk) / 2.74 <= 1. / rk <= fk * 1. / min_Rv + (1. - fk) / 2.74
                ]

        npts = len(pts)

        # Pet the user
        print('number of initially requested points = {0:d}'.format(niter))
        print('number of valid points = {0:d} (based on restrictions in R(V) versus f_bump plane)'.format(npts) )

        # setup of output
        N = N0 * npts
        cols = {'Av': np.empty(N, dtype=float),
                'Rv': np.empty(N, dtype=float),
                'Rv_MW': np.empty(N, dtype=float),
                'f_bump': np.empty(N, dtype=float)
                }
    else:
        it = np.nditer(np.ix_(avs, rvs))
        niter = np.size(avs) * np.size(rvs)

        pts = [ (float(ak), float(rk)) for ak, rk in it]
        npts = niter

        # setup of output
        N = N0 * niter
        cols = {'Av': np.empty(N, dtype=float),
                'Rv': np.empty(N, dtype=float)
                }

    # Generate the Grid
    # =================
    print('Generating a final grid of {0:d} points'.format(N))

    keys = g0.keys()
    for key in keys:
        cols[key] = np.empty(N, dtype=float)

    _seds = np.empty( (N, len(filter_names)), dtype=float)

    with progressbar.PBar(npts, txt='SED grid') as Pbar:
        for count, pt in enumerate(pts):
            # info showing program is running
            Pbar.update(count, force=True)

            if with_fb:
                Av, Rv, f_bump = pt
                Rv_MW = extLaw.get_Rv_A(Rv, f_bump)
                temp_results = g0.getSEDs(filter_names, extLaw=extLaw, Av=Av, Rv=Rv, f_bump=f_bump)
                # adding the dust parameters to the models
                cols['f_bump'][N0 * count:N0 * (count + 1)] = f_bump
                cols['Rv_MW'][N0 * count: N0 * (count + 1)] = Rv_MW
            else:
                Av, Rv = pt
                temp_results = g0.getSEDs(filter_names, extLaw=extLaw, Av=Av, Rv=Rv)
                # adding the dust parameters to the models
                cols['Av'][N0 * count: N0 * (count + 1)] = Av
                cols['Rv'][N0 * count: N0 * (count + 1)] = Rv

            # assign the extinguished SEDs to the output object
            _seds[N0 * count: N0 * (count + 1)] = temp_results.seds[:]

            # copy the rest of the parameters
            for key in keys:
                cols[key][N0 * count: N0 * (count + 1)] = g0.grid[key]

    _lamb = temp_results.lamb[:]
    # free the memory of temp_results
    del temp_results

    # Ship
    g = SpectralGrid(_lamb, seds=_seds, grid=Table(cols), backend='memory')
    g.grid.header['filters'] = ' '.join(filter_names)
    return g


#=================== TESTUNITS ============================

def test_gen_spectral_grid_from_stellib_given_points():
    """test_gen_spectral_grid_from_stellib_given_points
    Make sure it runs and returns a grid
    """
    osl = stellib.Kurucz()
    oiso = isochrone.padova2010()
    gen_spectral_grid_from_stellib_given_points(osl, oiso.data)


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
    oiso = isochrone.padova2010()

    # a spectral grid will be generated but not saved
    # using the stellar parameters by interpolation of the isochrones and the generation of spectra into the physical units
    grid_fname = 'kurucz2004.spectral.grid.hd5'

    # define filters for the grid
    filter_names  = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    # variable to ensure that range is fully covered in using np.arange
    __tiny_delta__ = 0.001
    # grid spacing for dust
    avs            = np.arange(0.0, 5.0 + __tiny_delta__, 0.1)
    rvs            = np.arange(1.0, 6.0 + __tiny_delta__, 0.5)
    fbumps         = np.asarray([1.0])

    #make the spectral grid
    g = gen_spectral_grid_from_stellib_given_points(osl, oiso.data)

    # make the grid
    extgrid = make_extinguished_grid(g, filter_names, avs, rvs, fbumps)

    # save grid to file
    extgrid.writeHDF(grid_fname)
