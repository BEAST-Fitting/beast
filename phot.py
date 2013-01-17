""" Trying to speed up the photometry part """
import numpy
import tables
import inspect
import os
import grid
import time
import numexpr
from decorators import timeit
from scipy.integrate import simps
localpath = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])

__default__      = localpath + '/libs/filters.hd5'
__default_vega__ = localpath + '/libs/vega.hd5'

# this is used to convert from bolometric luminosities to abs fluxes
# object to 10parsecs -- abs mag.
distc = 4. * numpy.pi * (3.1e19) ** 2


class filter(object):
    """Class filter
    Define a filter by its name, wavelength and transmission
    """
    #----------------------------------------------------------------------
    def info(self):
        """ display information about the current filter"""
        print "Filter object information:"
        print "   name: %s" % self.name
        print "   central wavelength: %f" % self.cl
        print "   norm: %f" % self.norm
        print "   pivot wavelength: %f" % self.lpivot
        print "   definition contains %d points" % self.transmit.size

    def __repr__(self):
        return "Filter: %s, %s" % (self.name, object.__repr__(self))

    def getFlux(self, slamb, sflux):
        """getFlux
        Integrate the flux within the filter and return the integrated energy
        If you consider applying the filter to many spectra, you might want to consider extractSEDs
        INPUTS:
           slamb: spectrum wavelength definition domain
           sflux: associated flux
        OUTPUTS:
           <float>: Energy of the spectrum within the filter
        """
        if True in numpy.isinf(sflux):
            indinf = numpy.where(numpy.isinf(sflux))
            indfin = numpy.where(numpy.isfinite(sflux))
            sflux[indinf] = numpy.interp(slamb[indinf], slamb[indfin], sflux[indfin])
        ifT = numpy.interp(slamb, self.wavelength, self.transmit, left=0., right=0.)
        if True in (ifT > 0.):
            ind = numpy.where(ifT > 0.)
            a = numpy.trapz( slamb[ind] * ifT[ind] * sflux[ind], slamb[ind] )
            b = numpy.trapz( slamb[ind] * ifT[ind], slamb[ind] )
            if (numpy.isinf(a) | numpy.isinf(b)):
                print self.name, "Warn for inf value"
            return a / b
        else:
            return 0.

    def __call__(self, slamb, sflux):
        return self.applyTo(slamb, sflux)

    def applyTo(self, slamb, sflux):
        """applyTo
        Apply filter to a spectrum
        INPUTS:
           slamb: spectrum wavelength definition domain
           sflux: associated flux
        OUTPUTS:
           [<float>]: new spectrum values accounting for the filter
        """
        ifT = numpy.interp(slamb, self.wavelength, self.transmit)
        return ifT * sflux

    def __init__(self, wavelength, transmit, name=''):
        """Constructor"""
        self.name       = name
        self.wavelength = wavelength
        self.transmit   = transmit
        self.norm       = simps(transmit, wavelength)
        self.lT         = simps(wavelength * transmit, wavelength)
        self.lpivot     = numpy.sqrt( self.lT / simps(transmit / wavelength, wavelength) )
        self.cl         = self.lT / self.norm


def __load__(fname, ftab, interp=True, lamb=None):
    """ Load a given filter from the library
        INPUTS:
            fname       str                     normalized names according to filtersLib
            ftab        hd5root                 root from the filter library hd5 file
        KEYWORDS:
            interp      bool                    reinterpolate the filters over given lambda points
            lamb        ndarray[float, ndim=1]  desired wavelength definition of the filter
        OUTPUTS:
            filter      filter                  filter object
    """
    fnode    = ftab.getNode('/filters/' + fname)
    flamb    = fnode[:]['WAVELENGTH']
    transmit = fnode[:]['THROUGHPUT']
    if interp & (lamb is not None):
        ifT = numpy.interp(lamb, flamb, transmit, left=0., right=0.)
        return filter( lamb, ifT, name=fnode.name )
    else:
        return filter( flamb, transmit, name=fnode.name )


def load_all_filters(interp=True, lamb=None, filterLib=__default__):
    """ load all filters from the library
        KEYWORDS:
            interp      bool                    reinterpolate the filters over given lambda points
            lamb        ndarray[float, ndim=1]  desired wavelength definition of the filter
            filterLib   path                    path to the filter library hd5 file
        OUTPUTS:
            filters     list[filter]            list of filter objects
    """
    with tables.openFile(filterLib, 'r') as ftab:
        filters = [ __load__(fname, ftab, interp=interp, lamb=lamb) for fname in ftab.root.content.cols.TABLENAME ]
    return(filters)


def load_filters(names, interp=True, lamb=None, filterLib=__default__):
    """ load a limited set of filters
        INPUTS:
            names       list[str]               normalized names according to filtersLib
        KEYWORDS:
            interp      bool                    reinterpolate the filters over given lambda points
            lamb        ndarray[float, ndim=1]  desired wavelength definition of the filter
            filterLib   path                    path to the filter library hd5 file
        OUTPUTS:
            filters     list[filter]            list of filter objects
    """
    with tables.openFile(filterLib, 'r') as ftab:
        filters = [ __load__(fname, ftab, interp=interp, lamb=lamb) for fname in names ]
    return(filters)


def extractPhotometry(lamb, spec, flist, absFlux=True):
    """ Extract seds from a one single spectrum

        INPUTS:
            lamb    ndarray[float,ndim=1]   wavelength of spec
            spec    ndarray[float, ndim=1]  spectrum
            flist   list[filter]            list of filter objects
        KEYWORDS:
            absflux bool                    return SEDs in absolute fluxes if set
        OUTPUT:
            cls     ndarray[float, ndim=1]  filters central wavelength
            seds    ndarray[float, ndim=1]  integrated sed
    """
    cls  = numpy.empty( len(flist), dtype=float)
    seds = numpy.empty( len(flist), dtype=float)
    for e, k in enumerate(flist):
        xl  = k.transmit > 0.
        tmp = lamb[xl] * k.transmit[xl]
        s0  = spec[:, xl]
        # apply absolute flux conversion if requested
        if absFlux:
            s0 /= distc
        a = simps( tmp[None, :] * s0, lamb[xl], axis=1 )
        seds[e] = a / k.lT
        cls[e]  = k.cl

    return cls, seds


def extractSEDs(g0, flist, absFlux=True):
    """ Extract seds from a grid

        INPUTS:
            g0      grid            Initial spectral grid
            flist   list[filter]    list of filter objects
        KEYWORDS:
            absflux bool            return SEDs in absolute fluxes if set
        OUTPUT:
            g       grid            SED grid object
    """
    lamb = g0.lamb
    seds = numpy.empty(( g0.grid.nrows, len(flist) ), dtype=float)
    cls  = numpy.empty( len(flist), dtype=float)
    for e, k in enumerate(flist):
        xl  = k.transmit > 0.
        tmp = lamb[xl] * k.transmit[xl]
        s0  = g0.seds[:, xl]
        # apply absolute flux conversion if requested
        if absFlux:
            s0 /= distc
        a = simps( tmp[None, :] * s0, lamb[xl], axis=1 )
        seds[:, e] = a / k.lT
        cls[e] = k.cl

    return grid.MemoryGrid(cls, seds, g0.grid)

def extractExtinguishedSEDs(g0, flist, extCurve, absFlux=True):
    """ Extract seds from a grid

        INPUTS:
            g0      grid            Initial spectral grid
            flist   list[filter]    list of filter objects
            extCurve np.exp(-tau)   extinction curve
        KEYWORDS:
            absflux bool            return SEDs in absolute fluxes if set
        OUTPUT:
            g       grid            SED grid object
    """
    lamb = g0.lamb
    seds = numpy.empty(( g0.grid.nrows, len(flist) ), dtype=float)
    for e, k in enumerate(flist):
        xl  = k.transmit > 0.
        tmp = lamb[xl] * k.transmit[xl]
        s0  = numexpr.evaluate('seds*extCurve',local_dict={'seds':g0.seds[:,xl],'extCurve':extCurve[xl]})
        a = simps( numexpr.evaluate('tmp*s0',local_dict={'tmp':tmp,'s0':s0}), lamb[xl], even='first', axis=1)
        seds[:, e] = a / k.lT
    #apply absolute flux conversion if needed
    if absFlux:
        seds /= distc
    return seds

def test(absFlux=True):
    """ Test units """
    gridfile = 'libs/stellib_kurucz2004_padovaiso.spectralgrid.fits'
    filter_names = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

	#Load the initial model grid
    with timeit('Loading Grid, %s' % gridfile):
        g0   = grid.FileSpectralGrid(gridfile)
        lamb = g0.lamb

    # Load the filters
    with timeit('Loading Filters (with interpolation)'):
        flist = load_filters(filter_names, interp=True, lamb=lamb)

    with timeit('SED integrations over %d ' % g0.grid.nrows):
        return extractSEDs(g0, flist, absFlux=absFlux)
