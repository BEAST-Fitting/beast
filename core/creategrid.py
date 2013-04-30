"""
Create extinguished grid

Dec 2012: Originally written by Kirill T.
Jan 2013: Modified by Karl G. to separate the extinguished grid creation into a separate program.
          Added restriction of parameter space in R(V) and f_bump.

"""

__version__ = '0.3dev'

import numpy
import pyfits
import numpy as np
from matplotlib.nxutils import points_inside_poly
from ..tools import progressbar
from . import grid
from . import stellib
from . import isochrone
from ..external import mytables
from ..external import eztables
from ..external import ezunits


def get_radius(logl, logt):
    """ Returns the radius of a star given its luminosity and temperature
            assuming a black body
            R ** 2 = L / ( 4 pi sig T ** 4 )
    INPUTS:
        logl   ndarray[float, ndim=1]   log luminosities from the isochrones *in Lsun*
        logt   ndarray[float, ndim=1]   log temperatures from the isochrones *in Kelvin*
    """
    lsun = 1. * ezunits.unit['lsun'].to('W').magnitude  # 3.839e26 W
    sig  = 5.67037321 * 1e-8 * ezunits.unit[' W * m**-2 * K**-4'].magnitude
    return np.sqrt( (10 ** logl) * lsun / (4.0 * np.pi * sig * ((10 ** logt) ** 4)) )


def gen_spectral_grid_from_stellib_given_points(outfile, osl, oiso, pts, bounds=dict(dlogT=0.1, dlogg=0.3)):
    """ Reinterpolate a given stellar spectral library on to an Isochrone grid
    INPUTS:
        outfile str                 fits file to export to
        osl     stellib.stellib     a stellar library
        oiso    isochrone.Isochrone an isochrone library
        pts     
    KEYWORDS:
        bounds  dict                sensitivity to extrapolation (see get_stellib_boundaries)

    OUTPUTS:
        None

        only write into outfile
    """
    assert(grid.isNestedInstance(osl, stellib.Stellib) )
    assert(grid.isNestedInstance(oiso, isochrone.Isochrone) )

    ndata = len(pts)
    _grid  = {}
    for k in oiso.data.keys():
        _grid[k] = np.empty(ndata, dtype=float )
    _grid['radius'] = np.empty(ndata, dtype=float )
    _grid['keep'] = np.empty(ndata, dtype=bool )

    specs = np.empty( (ndata + 1, len(osl.wavelength)), dtype=float )
    specs[-1] = osl.wavelength[:]

    _bounds = osl.get_boundaries(dlogT=bounds['dlogT'], dlogg=bounds['dlogg'], closed=True)

    # compute logg, Teff of the required points,
    # this means interp iso
    # compute the spectrum of the point as well
    # store all
    # then check the boundary conditions to filter them out
    radii = get_radius(pts['logL'], pts['logT'])
    weights = 4. * np.pi * (radii * 1e2) ** 2  # denorm models are in cm**-2 (4*pi*rad)
    #check boundaries, keep the data but do not compute the sed if not needed
    data = np.array([pts['logg'], pts['logT']]).T
    bound_cond = points_inside_poly(data, _bounds)
    del data
    _grid['keep'] = bound_cond[:]
    _grid['radius'] = radii[:]
    for key in pts.keys():
        _grid[key] = pts[key]
    with progressbar.PBar(ndata, txt='spectral grid') as Pbar:
        for mk, r in enumerate(pts):
            Pbar.update(mk)
            if bound_cond[mk]:
                try:
                    s = np.array( osl.interp(r['logT'], r['logg'], r['Z'], 0.) ).T
                    specs[mk, :] = osl.genSpectrum(s) * weights[mk]
                except:
                    #print "error"
                    specs[mk, :] = np.zeros(len(osl.wavelength), dtype=float )                
                    #assert(False)
                    
    #filter unbound values
    idx = np.array(_grid.pop('keep'))

    #specs = np.vstack( [specs, osl.wavelength] )
    specs = np.vstack( [specs.compress(idx, axis=0), osl.wavelength] )
    for k in _grid.keys():
            _grid[k] = _grid[k][idx]

    rsun = ezunits.unit['Rsun'].to('m').magnitude  # 6.955e8 m
    _grid['radius'] /= rsun

    pars  = mytables.Table(_grid, name='Reinterpolated stellib grid')
    pars.header['stellib'] = osl.source
    pars.header['isoch'] = oiso.source
    pars.setUnit('radius', 'Rsun')
    #assert(False)
    pyfits.writeto(outfile, specs, clobber=True)
    pars.write(outfile, append=True)


def gen_spectral_grid_from_stellib(outfile, osl, oiso, ages=(1e7,), masses=(3,), Z=(0.02,), bounds=dict(dlogT=0.1, dlogg=0.3)):
    """ Reinterpolate a given stellar spectral library on to an Isochrone grid
    INPUTS:
        outfile str                 fits file to export to
        osl     stellib.stellib     a stellar library
        oiso    isochrone.Isochrone an isochrone library
    KEYWORDS:
        ages    iterable            list of age points to include in the grid   (in Yr)
        masses  iterable            list of mass points to include in the grid  (M/Msun)
        ages    iterable            list of metallicity points to include in the grid (Z/Zsun)
        bounds  dict                sensitivity to extrapolation (see get_stellib_boundaries)

    OUTPUTS:
        None

        only write into outfile
    """
    assert(grid.isNestedInstance(osl, stellib.Stellib) )
    assert(grid.isNestedInstance(oiso, isochrone.Isochrone) )

    _ages, _Z = np.ix_(ages, Z)
    _masses = np.asarray(masses)
    it = np.nditer([_ages, _Z])

    niter = len(ages) * len(Z)
    ndata = niter * len(masses)

    _grid  = {}
    for k in oiso.data.keys():
        _grid[k] = np.empty(ndata, dtype=float )
    _grid['radius'] = np.empty(ndata, dtype=float )
    _grid['keep'] = np.empty(ndata, dtype=bool )

    specs = np.empty( (ndata + 1, len(osl.wavelength)), dtype=float )
    specs[-1] = osl.wavelength[:]

    _bounds = osl.get_boundaries(dlogT=bounds['dlogT'], dlogg=bounds['dlogg'], closed=True)

    # compute logg, Teff of the required points,
    # this means interp iso
    # compute the spectrum of the point as well
    # store all
    # then check the boundary conditions to filter them out
    kdata = 0
    with progressbar.PBar(niter, txt='spectral grid') as Pbar:
        for k, (_ak, _Zk) in enumerate(it):
            Pbar.update(k)
            r = oiso._get_isochrone(_ak, metal=_Zk, masses=np.log10(_masses))
            radii = get_radius(r['logL'], r['logT'])
            weights = 4. * np.pi * (radii * 1e2) ** 2  # denorm models are in cm**-2 (4*pi*rad)
            #check boundaries, keep the data but do not compute the sed if not needed
            data = np.array([r['logg'], r['logT']]).T
            bound_cond = points_inside_poly(data, _bounds)
            del data
            start_idx = k * len(masses)
            end_idx   = start_idx + r.nrows
            _grid['keep'][start_idx: end_idx] = bound_cond[:]
            _grid['radius'][start_idx: end_idx] = radii[:]
            for key in r.keys():
                _grid[key][start_idx: end_idx] = r[key]
            for mk in range(r.nrows):
                if bound_cond[mk]:
                    try:
                        s = np.array( osl.interp(r['logT'][mk], r['logg'][mk], _Zk, 0.) ).T
                        specs[kdata, :] = osl.genSpectrum(s) * weights[mk]
                    except:
                        print "error"
                        specs[kdata, :] = np.zeros(len(osl.wavelength), dtype=float )
                        #assert(False)
                else:
                    specs[kdata, :] = np.zeros(len(osl.wavelength), dtype=float )
                kdata += 1

    #filter unbound values
    idx = np.array(_grid.pop('keep'))

    #specs = np.vstack( [specs, osl.wavelength] )
    specs = np.vstack( [specs.compress(idx, axis=0), osl.wavelength] )
    for k in _grid.keys():
            _grid[k] = _grid[k][idx]

    rsun = ezunits.unit['Rsun'].to('m').magnitude  # 6.955e8 m
    _grid['radius'] /= rsun

    pars  = mytables.Table(_grid, name='Reinterpolated stellib grid')
    pars.header['stellib'] = osl.source
    pars.header['isoch'] = oiso.source
    pars.setUnit('radius', 'Rsun')

    pyfits.writeto(outfile, specs, clobber=True)
    pars.write(outfile, append=True)

def merge_spectral_grids(preferred_fname, alt_fname, outname):
    """
    Combine two spectral grids using points from pref_grid wherever available and
    points from alt_grid elsewhere.
    INPUTS:
        preferred_fname   string    Name of preferred spectral grid
        alt_fname         string    Name of spectral grid to use where preferred grid has no coverage
        outname           string    What to save the output as
    
    """
    alt_grid = grid.FileSpectralGrid(alt_fname)
    pref_grid = grid.FileSpectralGrid(preferred_fname)

    alt_inds = np.where((np.in1d(alt_grid.logL, pref_grid.logL)==False)*(np.in1d(alt_grid.logg, pref_grid.logg)==False))[0]

    tab = eztables.Table()
    for key in pref_grid.grid.keys():
        tab.addCol(key, np.hstack([alt_grid.grid[key][alt_inds], pref_grid.grid[key]]))

    merged = grid.MemoryGrid(pref_grid.lamb, seds=np.vstack([alt_grid.seds[alt_inds],pref_grid.seds]), grid=tab)
    merged.write(outname, clobber=True)


def make_extinguished_grid(stellar_filename, filter_names, extLaw, avs, rvs, fbumps):

    """
    Extinguish and extract fluxes through filters
       (all wavelengths in stellar SEDs and filter response functions assumed to be in Angstroms)
    INPUTS:
        stellar_filename string    FITS file with stellar SEDs (luminosities)
        filter_names   list        list of filter names according to the filter lib
    KEYWORDS:
        Av_vals     numpy array    Av values to iterate over
        Rv_vals     numpy array    Rv values to iterate over
        f_bump_vals numpy array    f_bump values to iterate over
    """
    # get the stellar grid (no dust)
    g0 = grid.FileSpectralGrid(stellar_filename)

    # reduce size of grid for debugging (TBR)
    #indx = numpy.arange(0, g0.grid.nrows, 100)
    #g0.grid = g0.grid.extract(indx)
    #g0.seds = g0.seds[indx]
    #g0.write('small_grid.fits', clobber=True)

    # define the extinction law be be used (fixed for now)
    #extLaw = extinction.RvFbumpLaw()

    # get the min/max R(V) values
    min_Rv = min(rvs)
    max_Rv = max(rvs)


    
    # create mesh from input 1d vectors
    Av_vals, Rv_vals, f_bump_vals = numpy.ix_(avs, rvs, fbumps)

    # setup interation over the full dust parameter grid
    it = numpy.nditer([Av_vals, Rv_vals, f_bump_vals])
    niter = Av_vals.size * Rv_vals.size * f_bump_vals.size

    # compute the allowed points based on the R(V) versus f_bump plane
    #  duplicates effort for all A(V) values, but it is quick compared to other steps
    # Note on 2.74: Bumpless extinction implies f_bump = 0. and Rv = 2.74
    pts = [ (float(ak), float(rk), float(fk)) for ak, rk, fk in it \
            if fk / max_Rv + (1. - fk) / 2.74 <= 1./rk <= fk * 1./min_Rv + (1. - fk)/2.74]
    npts = len(pts)

    print 'n possible = ', niter
    print 'n actual   = ', npts, ' (based on restrictions in R(V) versus f_bump plane'

    # setup of output
    cols = {'Av': numpy.empty(g0.grid.nrows * npts), 'Rv': numpy.empty(g0.grid.nrows * npts), 'f_bump': numpy.empty(g0.grid.nrows * npts)}

    for key in g0.keys():
        cols[key] = numpy.empty(g0.grid.nrows * npts)

    count = 0
    print 'Generating a final grid of %d points' % (g0.grid.nrows * npts)

    with progressbar.PBar(npts, txt='SED grid') as Pbar:
        for Av, Rv, f_bump in pts:
            # info showing program is running
            Pbar.update(count, force=True)

            # compute R(V)^MW to return the Rv requested
            if f_bump > 0.:
                Rv_MW = 1./(1./(Rv*f_bump) - (1. - f_bump)/(f_bump*2.74))
            else:
                Rv_MW = 2.74  # doesn't matter

            # apply extinction and integrate over band response functions
            temp_results = g0.getSEDs(filter_names, extLaw=extLaw, Av=Av, Rv=Rv_MW, f_bump=f_bump)
            # setup the output object
            #  must be done after the 1st extraction to get the wavelength vector
            if count is 0:
                results = grid.MemoryGrid(temp_results.lamb,
                                          seds=numpy.empty( (g0.grid.nrows * npts, len(filter_names)) ),
                                          grid=eztables.Table(cols) )

            # assign the extinguished SEDs to the output object
            results.seds[g0.grid.nrows * count: g0.grid.nrows * (count + 1)] = temp_results.seds

            # free the memory of temp_results
            del temp_results

            # adding the dust parameters to the models
            results.grid.data['Av'][g0.grid.nrows * count:g0.grid.nrows * (count + 1)] = Av
            results.grid.data['Rv'][g0.grid.nrows * count:g0.grid.nrows * (count + 1)] = Rv
            results.grid.data['f_bump'][g0.grid.nrows * count:g0.grid.nrows * (count + 1)] = f_bump

            # the rest of the parameters
            for key in g0.keys():
                results.grid.data[key][g0.grid.nrows * count:g0.grid.nrows * (count + 1)] = g0.grid[key]
            count += 1

    results.grid.header['filters'] = ' '.join(filter_names)
    return results


#=================== TESTUNITS ============================

def test_gen_spectral_grid_from_stellib():
    osl = stellib.Kurucz()
    oiso = isochrone.padova2010()
    gen_spectral_grid_from_stellib('tmp.fits', osl, oiso,
                                   ages=(1e7, 1e8),
                                   masses = (3., 4.),
                                   Z=(0.02, 0.004),
                                   bounds=dict(dlogT=0.1, dlogg=0.3))


def main():

    #select the stellar library and the isochrones to use
    osl = stellib.Kurucz()
    oiso = isochrone.padova2010()

    # a spectral grid will be generated using the stellar parameters by interpolation of the isochrones and the generation of spectra into the physical units
    spectral_grid_fname = 'kurucz2004.spectral.grid.fits'

    # a photometric grid precomputing attenuation values as well
    #    sed_grid_fname = spectral_grid_fname.replace('spectral', 'seds')
    sed_grid_fname = 'fbump_only.fits'

    # define filters for the grid
    filter_names  = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    # grid spacing for stars
    __tiny_delta__ = 0.001
    # ages           = (1e7, 1e8, 1e9)
    # masses         = (1., 2., 3., 4., 50.)
    # Z              = (0.02, 0.004)
    ages           = 10 ** numpy.arange(6., 9. + __tiny_delta__, 0.1)
    masses         = 10 ** numpy.arange(0.5, 20 + __tiny_delta__, 0.1)
    Z              = (0.02)
    # grid spacing for dust
    # variable to ensure that range is fully covered in using numpy.arange
    avs            = numpy.arange(0.0, 5.0 + __tiny_delta__, 0.1)
    rvs            = numpy.arange(1.0, 6.0 + __tiny_delta__, 0.5)
    fbumps         = numpy.asarray([1.0])
    #fbumps         = numpy.arange(0.0, 1. + __tiny_delta__, 0.1)
    #avs            = numpy.arange(0.0, 5.0 + __tiny_delta__, 1.)
    #rvs            = numpy.arange(1.0, 6.0 + __tiny_delta__, 1.)
    #fbumps         = numpy.arange(0.0, 1. + __tiny_delta__, 0.5)

    #filter extrapolations of the grid with given sensitivities in logg and logT
    bounds = dict(dlogT=0.1, dlogg=0.3)

    #make the spectral grid
    #gen_spectral_grid_from_stellib(spectral_grid_fname, osl, oiso, ages=ages, masses=masses, Z=Z, bounds=bounds)

    # make the grid
    extgrid = make_extinguished_grid(spectral_grid_fname, filter_names, avs, rvs, fbumps)

    # save grid to file
    extgrid.write(sed_grid_fname, clobber=True)


def generate_dense_SED_grid(Av_min=0.0, Av_max=5.0, Av_step=1.0,
                            Rv_min=2.0, Rv_max=6.0, Rv_step=2.0,
                            fb_min=0.0, fb_max=1.0, fb_step=0.7,
                            sed_grid_outname='sedfitter/libs/tlusty_kurucz.padova.sed.grid.fits',
                            filters='hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()):
    """
    Keywords:
       *_min, *_max, *_step   floats    Minimum, maximum, and step values for dust parameters
       filters                list      Filter names (need to match those in filter definition file)
       sed_grid_outname       string    Path and name from sedfitter main directory to save
                                        dense SED grid in

    """
    from . import extinction
    from .. import config

    extLaw = extinction.RvFbumpLaw()

    to_sedfitter_dir = config.os.path.realpath(config.__file__).replace('config.pyc','')

    spec_grid_kurucz_fname= to_sedfitter_dir + 'libs/kurucz.padova.spectral.grid.fits'
    spec_grid_tlusty_fname= to_sedfitter_dir + 'libs/tlusty.padova.spectral.grid.fits'
    spec_grid_combined_fname= to_sedfitter_dir + 'libs/tlusty_kurucz.padova.spectral.grid.fits'

    oiso = isochrone.padova2010()
    iso_select = oiso.data.selectWhere('*', '(Z==0.02)')
    osl = stellib.Kurucz()
    gen_spectral_grid_from_stellib_given_points(spec_grid_kurucz_fname, osl, oiso, iso_select)
    osl = stellib.Tlusty()
    gen_spectral_grid_from_stellib_given_points(spec_grid_tlusty_fname, osl, oiso, iso_select)
    del osl, oiso

    merge_spectral_grids(spec_grid_tlusty_fname, spec_grid_kurucz_fname, spec_grid_combined_fname)
    
    tiny_delta = np.min([Av_step, Rv_step, fb_step])/2
    Avs = np.arange(Av_min, Av_max+tiny_delta, Av_step)
    Rvs = np.arange(Rv_min, Rv_max+tiny_delta, Rv_step)
    fbs = np.arange(fb_min, fb_max+tiny_delta, fb_step)

    extgrid = make_extinguished_grid(spec_grid_combined_fname, filters, extLaw, Avs, Rvs, fbs)

    extgrid.write(sed_grid_outname, clobber=True)
