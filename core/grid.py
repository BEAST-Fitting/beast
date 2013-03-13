""" Manage Various SED/spectral grids is a generic way """
import numpy
import pyfits

#from . import stellib
from . import phot
#from . import isochrone
from . import extinction
from ..external.eztables import Table

#from tools.decorators import timeit

from copy import deepcopy


def isNestedInstance(obj, cl):
    """ Test for sub-classes types
        I could not find a universal test
    """
    tree = []
    for k in cl.__subclasses__():
        tree += k.__subclasses__()
    tree += cl.__subclasses__() + [ cl ]
    return  issubclass(obj.__class__, tuple(tree))


class ModelGrid(object):
    """ Generic class for a minimum update of future codes """
    def __init__(self, *args, **kwargs):
        self.lamb = None
        self.seds = None
        self.grid = None

    def keys(self):
        """ returns the grid dimension names """
        return []

    def getGridPoints(self, *args, **kwargs):
        """ Returns age, mass, logg, logT, logL... """
        pass

    def getPDF(self, Qname, lnp, *args, **kwargs):
        assert (Qname in self.keys() ), "Cannot find %s in the grid description" % Qname

    def write(self, fname, *args, **kwargs):
        if ( (self.lamb is not None) & (self.seds is not None) & (self.grid is not None) ):
            assert (isinstance(self.grid, Table)), 'Only eztables.Table are supported so far'
            r = numpy.vstack( [ numpy.copy(self.seds), self.lamb ])
            pyfits.writeto(fname, r, **kwargs)
            if getattr(self, 'filters', None) is not None:
                if ('FILTERS' not in self.grid.header.keys()):
                    self.grid.header['FILTERS'] = self.filters
            self.grid.write(fname, append=True)

    def copy(self):
        """ returns a copy of the object """
        return deepcopy(self)


class MemoryGrid(ModelGrid):
    """ Instanciate an grid object that has no physical storage
        Helps to create new grids on the fly. Because it deriveds from
        ModelGrid, this can be exported on disk too.
    """
    def __init__(self, lamb, seds=None, grid=None):
        """ MemoryGrid constructor
        INPUTS:
            lamb    ModelGrid or subclass   New ref to the given grid (debug purpose)
                lamb            wavelengths
                sed         seds
                grid            grid associated to seds
        """
        if isNestedInstance(lamb, ModelGrid):
            self.lamb = lamb.lamb
            self.seds = lamb.seds
            self.grid = lamb.grid
        else:
            assert ((seds is not None) & (grid is not None)), 'Wrong number of arguments'
            self.lamb = lamb
            self.seds = seds
            self.grid = grid


class SpectralGrid(ModelGrid):
    """ Generate a grid that contains spectra.
    It provides an access to integrated photometry function getSEDs """

    def getSEDs(self, filter_names, absFlux=True, extLaw=None, inplace=False, **kwargs):
        """
        Extract integrated fluxes through filters
        INPUTS:
            filter_names    list    list of filter names according to the filter lib
        KEYWORDS:
            absFlux         bool    returns absolute fluxes if set
            extLaw          extinction.ExtinctionLaw    apply extinction law if provided
            inplace         bool                        if set, do not copy the grid and apply on it

            **kwargs        extra keywords will be forwrded to extLaw
        """
        if type(filter_names[0]) == str:
            flist = phot.load_filters(filter_names, interp=True, lamb=self.lamb)
            _fnames = filter_names
        else:
            flist = filter_names
            _fnames = [ fk.name for fk in filter_names ]
        if extLaw is not None:
            if not inplace:
                r = self.applyExtinctionLaw(extLaw, inplace=inplace, **kwargs)
                return phot.extractSEDs(r, flist, absFlux=absFlux)
            else:
                self.applyExtinctionLaw(extLaw, inplace=inplace, **kwargs)
                r = self
        lamb, seds, grid = phot.extractSEDs(self, flist, absFlux=absFlux)
        memgrid = MemoryGrid(lamb, seds, grid)
        setattr(memgrid, 'filters', _fnames)
        return memgrid

    def applyExtinctionLaw(self, extLaw, inplace=False, **kwargs):
        """
        Apply an extinction law to the model grid
        INPUTS:
            extLaw          extinction.ExtinctionLaw    apply extinction law if provided
        KEYWORDS:
            inplace         bool                        if set, do not copy the grid and apply on it

            **kwargs        extra keywords will be forwrded to extLaw
        """
        assert( isinstance(extLaw, extinction.ExtinctionLaw)), 'Expecting ExtinctionLaw object got %s' % type(extLaw)
        extCurve = numpy.exp(-1. * extLaw.function(self.lamb[:], **kwargs))
        if not inplace:
            g = self.copy()
            g.seds *= extCurve[None, :]
            g.grid.header['ExtLaw'] = extLaw.name
            for k, v in kwargs.iteritems():
                g.grid.header[k] = v
            return g
        else:
            self.grid.header['ExtLaw'] = extLaw.name
            for k, v in kwargs.iteritems():
                self.grid.header[k] = v
            self.seds *= extCurve[None, :]


class FileSEDGrid(SpectralGrid):
    """ Generate a grid from a spectral library """

    def __init__(self, fname, *args, **kwargs):
        with pyfits.open(fname) as f:
            self.seds = f[0].data[:-1]
            self.lamb = f[0].data[-1]
        self.grid = Table(fname)
        self.filters = self.grid.header.get('FILTERS', None)
        if self.filters is not None:
            self.filters = self.filters.split()
        #lamb, seds = self.getSEDs(filters, self.osl.wavelength, self.osl.spectra)
        for k in self.grid.keys():
            self.__dict__[k] = self.grid[k]

        if self.filters is None:
            print 'Warning: Filter names where not found!'
            print '(This may happen if you loaded an older version of grid file)'
            print 'Please correct by setting the filters: self.filters = [...]'

    def keys(self):
        """ returns the grid dimension names """
        return self.grid.keys()


class FileSpectralGrid(SpectralGrid):
    """ Generate a grid from a spectral library """

    def __init__(self, fname, *args, **kwargs):
        with pyfits.open(fname) as f:
            self.seds = f[0].data[:-1]
            self.lamb = f[0].data[-1]
        self.grid = Table(fname)
        #lamb, seds = self.getSEDs(filters, self.osl.wavelength, self.osl.spectra)
        for k in self.grid.keys():
            self.__dict__[k] = self.grid[k]

    def keys(self):
        """ returns the grid dimension names """
        return self.grid.keys()


class StellibGrid(SpectralGrid):
    """ Generate a grid from a spectral library """

    def __init__(self, osl, filters, *args, **kwargs):
        self.osl = osl  # stellib.BaSeL()
        lamb, seds = self.getSEDs(filters, self.osl.wavelength, self.osl.spectra)
        self.lamb = lamb
        self.seds = seds
        self.filters = filters
        self.grid = self.osl.grid
        for k in self.grid.keys():
            self.__dict__[k] = self.grid[k]

    def keys(self):
        """ returns the grid dimension names """
        return self.grid.keys()


## Not used anywhere
#def generate_spectral_grid_from_isochrones(outfile, osl, oiso, Z=0.02):
#    """ Reinterpolate a given stellar spectral library on to an Isochrone grid
#    INPUTS:
#        outfile     str         fits file to export to
#        osl     stellib.stellib     a stellar library
#        oiso        isochrone.Isochrone an isochrone library
#        Z       float           metallicity to use
#
#    OUTPUTS:
#        None
#
#        only write into outfile
#    """
#
#    assert(isNestedInstance(osl, stellib.Stellib) )
#    assert(isNestedInstance(oiso, isochrone.Isochrone) )
#    specs = numpy.empty( (oiso.data.nrows + 1, len(osl.wavelength)), dtype=float )
#    specs[-1] = osl.wavelength[:]
#
#    progress = 0
#    with timeit('interpolation'):
#        for k in range(oiso.data.nrows):
#            if progress < int(100 * (k + 1) / oiso.data.nrows):
#                progress = int(100 * (k + 1) / oiso.data.nrows)
#                print "progress... %d / 100" % progress
#            r = numpy.array( osl.interp(oiso.data['logT'][k], oiso.data['logg'][k], Z, oiso.data['logL'][k]) ).T
#            specs[k, :] = osl.genSpectrum(r)
#    pyfits.writeto(outfile, specs)
#
#    #copy pars
#    data = {}
#    for k in oiso.data.keys():
#        data[k] = oiso.data[k]
#    pars  = Table(data, name='Reinterpolated stellib grid')
#    pars.header['stellib'] = osl.source
#    pars.header['isoch'] = oiso.source
#
#    pars.write(outfile, append=True)
