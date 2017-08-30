""" This module gives tools to generate Kurucz grid from original downloads """
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

#from anased import grid
from ..core import grid
from ..core import stellib, isochrone
from ..decorators import timeit
from ..external.eztables import Table

import pyfits
import numpy as np
#import glob
from matplotlib.nxutils import points_inside_poly
import sys


def __treatSingleFile__(fname='ckp00_10000.fits'):
    """ Grab the useful data from a single file

    Used in Kurucz_to_Stellib

    INPUTS:
        fname   str file to process
    OUTPUTS:
        logg, teff, logz, lamb, data
    """

    with pyfits.open(fname) as f:
        logz  = f[0].header['log_Z']
        teff  = f[0].header['TEFF']
        d     = f[1].data
        logg  = [ float(k[1:]) * 0.1 for k in d.dtype.names if k[0] == 'g']
        d     = np.array(d.tolist())
        lamb  = d[:, 0]
        data  = d[:, 1:].T
        nspec = data.shape[0]
        logz  = np.asarray([logz] * nspec)
        teff  = np.asarray([teff] * nspec)
    return logg, teff, logz, lamb, data


def Kurucz_to_Stellib(lst):
    """ Extract SED parameters and spectra from a single file and collapse
    them into a MemoryGrid object

    INPUTS:
        lst list    list of files to process
                e.g. lst = glob.glob('ck*/*fits')

    OUTPUTS:
        g   MemoryGrid  stellib grid in the original units
    """
    d = dict(logg=[], teff=[], logz=[])
    specs = []

    for kfile in lst:
        logg, teff, logz, lamb, data = __treatSingleFile__(kfile)
        d['logg'].append(logg)
        d['teff'].append(teff)
        d['logz'].append(logz)
        specs.append(data)

    #post process
    for k in d:
        d[k] = np.ravel(d[k])
    specs = np.vstack(specs)
    t     = Table(d)
    g     = grid.MemoryGrid(lamb, specs, t)
    return g


def gen_spectral_grid_from_kurucz(outfile, osl, oiso, Z=0.02):
    """ Reinterpolate a given stellar spectral library on to an Isochrone grid
    INPUTS:
        outfile     str         fits file to export to
        osl     stellib.stellib     a stellar library
        oiso        isochrone.Isochrone an isochrone library
        Z       float           metallicity to use

    OUTPUTS:
        None

        only write into outfile
    """
    assert(grid.isNestedInstance(osl, stellib.Stellib) )
    assert(grid.isNestedInstance(oiso, isochrone.Isochrone) )
    specs = np.empty( (oiso.data.nrows + 1, len(osl.wavelength)), dtype=float )
    specs[-1] = osl.wavelength[:]

    def get_radius(logl, logt):
        #get the radius of a star given its luminosity and temperature
        #(assuming a black body)
        lsun = 3.839e26  # W
        sig  = 5.67037321 * 1e-8   # W m**-2 K**-4
        return np.sqrt( (10 ** logl) * lsun / (4.0 * np.pi * sig * ((10 ** logt) ** 4)) )

    bounds = get_stellib_boundaries(osl, dlogT=0.1, dlogg=0.3, closed=True)

    progress = 0
    data = np.array([oiso.data['logg'], oiso.data['logT']]).T
    bound_cond = points_inside_poly(data, bounds)
    del data
    radii = get_radius(oiso.data['logL'], oiso.data['logT'])
    weights = 4. * np.pi * (radii * 1e2) ** 2  # denorm models are in cm**-2 (4 * pi * rad)
    with timeit('interpolation'):
        for k in range(oiso.data.nrows):
            p = int(100 * (k + 1) / oiso.data.nrows)
            if (progress < p):
                progress = p
                sys.stdout.write("progress... %d / 100\r" % progress )
            if bound_cond[k] is True:
                r = np.array( osl.interp(oiso.data['logT'][k], oiso.data['logg'][k], Z, 0.) ).T
                specs[k, :] = osl.genSpectrum(r) * weights[k]
            else:
                specs[k, :] = np.zeros(len(osl.wavelength), dtype=float )
        sys.stdout.write("progress... %d / 100" % progress )

    specs = specs[bound_cond, :]

    pyfits.writeto(outfile, specs)

    #copy pars
    data = {}
    for k in list(oiso.data.keys()):
        data[k] = oiso.data[k][bound_cond]
    data['radius'] = radii[bound_cond] / 6.955e8  # Rsun
    pars  = Table(data, name='Reinterpolated stellib grid')
    pars.header['stellib'] = osl.source
    pars.header['isoch'] = oiso.source

    pars.write(outfile, append=True)


def get_stellib_boundaries(s, dlogT=0.1, dlogg=0.3, closed=True):
    """ Returns the closed boundary polygon around the stellar library with
    given margins

    INPUTS:
        s   Stellib     Stellar library object

    KEYWORDS:
        dlogT   float       margin in logT
        dlogg   float       margin in logg
        closed  bool        if set, close the polygon

    OUTPUTS:
        b   ndarray[float, ndim=2]  (closed) boundary points: [logg, Teff]

    Note:
        use "points_inside_poly" to test wether a point is inside the limits
        >>> data = np.array([iso.data['logg'], iso.data['logT']]).T
        >>> aa = points_inside_poly(data, leftb)
    """
    leftb   = [(k, np.max(s.logT[s.logg == k]) + dlogT ) for k in np.unique(s.logg)]
    leftb  += [ (leftb[-1][0] + dlogg, leftb[-1][1]) ]
    leftb   = [ (leftb[0][0] - dlogg, leftb[0][1]) ] + leftb
    rightb  = [(k, np.min(s.logT[s.logg == k]) - dlogT ) for k in np.unique(s.logg)[::-1]]
    rightb += [ (rightb[-1][0] - dlogg, rightb[-1][1]) ]
    rightb  = [ (rightb[0][0] + dlogg, rightb[0][1]) ] + rightb
    b = leftb + rightb
    if closed:
        b += [b[0]]
    return np.array(b)


#=========================== TEST UNITS ==================================
def test_stellib_boundaries():
    """ Test get_stellib_boundaries function """
    import pylab as plt

    s   = stellib.Kurucz()
    iso = isochrone.padova2010()
    leftb = get_stellib_boundaries(s, .1, 0.3, True)

    plt.plot(s.Teff, s.logg, '.', color='k')
    plt.plot(leftb[:, 1], leftb[:, 0], 'o-')

    #leftb = [logg, teff]
    plt.plot(iso.data['logT'], iso.data['logg'], ',', color='r')
    data = np.array([iso.data['logg'], iso.data['logT']]).T
    aa   = points_inside_poly(data, leftb)
    plt.plot(iso.data['logT'][aa], iso.data['logg'][aa], ',', color='g')
    plt.xlabel('log(Teff)')
    plt.ylabel('log(g)')
    plt.xlim(plt.xlim()[::-1])
    plt.ylim(plt.ylim()[::-1])


def main_last_grid():
    s   = stellib.Kurucz()
    iso = isochrone.padova2010()
    gen_spectral_grid_from_kurucz('tmp.fits', s, iso)
