""" Manage Various SED/spectral grids is a generic way """

import numpy
import pyfits

import stellib
import phot
import isochrone
import mytables
import extinction

from decorators import timeit
from time import time

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
            assert (isinstance(self.grid, mytables.Table)), 'Only mytables.Table are supported so far'
            r = numpy.vstack( [ numpy.copy(self.seds), self.lamb ])
            pyfits.writeto(fname, r, **kwargs)
            self.grid.write(fname, append=True)


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

    def getSEDs(self, filter_names, absFlux=True):
        """
        Extract integrated fluxes through filters
        INPUTS:
            filter_names    list    list of filter names according to the filter lib
        """
        if type(filter_names[0]) == str:
            flist = phot.load_filters(filter_names, interp=True, lamb=self.lamb)
        else:
            flist = filter_names
        return phot.extractSEDs(self, flist, absFlux=absFlux)

    def getExtinguishedSEDs(self, filter_names, toPhotLamb=1.0e4,toFbumpLawLamb=1.0,absFlux=True, Av_vals=numpy.array([3.0]), Rv_vals = numpy.array([2.7]), f_bump_vals = numpy.array([0.0])):
        """
        Extinguish and extract fluxes through filters
        INPUTS:
              filter_names   list    list of filter names according to the filter lib
        KEYWORDS:
              toPhotLamb     float   conversion factor from grid wavelength units to angstroms
              toFbumpLawLamb float   conversion factor from grid wavelength units to microns
              Av_vals   numpy array  Av values to iterate over
              Rv_vals   numpy array  Rv values to iterate over
              f_bump_vals numpy array f_bump values to iterate over
        """
        extLaw = extinction.RvFbumpLaw()
        
        if type(filter_names[0]) == str:
            flist = phot.load_filters(filter_names, interp=True, lamb=self.lamb*toPhotLamb)
        else:
            flist = filter_names

        lamb  = numpy.empty( len(flist), dtype=float)
        for e,k in enumerate(flist):
            lamb[e] = k.cl
        
        Av_vals, Rv_vals, f_bump_vals = numpy.ix_(Av_vals, Rv_vals, f_bump_vals)
        it = numpy.nditer([Av_vals,Rv_vals,f_bump_vals])
        ext_lamb = self.lamb*toFbumpLawLamb

        niter = Av_vals.size*Rv_vals.size*f_bump_vals.size
        cols = {'Av':numpy.empty(self.grid.nrows*niter),'Rv':numpy.empty(self.grid.nrows*niter),'f_bump':numpy.empty(self.grid.nrows*niter)}
            
        for key in self.keys():
            cols[key] = numpy.empty(self.grid.nrows*niter)

        results = MemoryGrid(lamb,seds= numpy.empty((self.grid.nrows*niter,filter_names.__len__())),grid=mytables.Table(iterable=cols,header=self.grid.header))

        itime = time()
        count = 0
        print "Memory allocated, starting on " + str(niter) + " iterations"
        for Av,Rv,f_bump in it:
            if ((count+1) % 100) == 0:
                print 100.*count/(Av_vals.size*Rv_vals.size*f_bump_vals.size),time()-itime
                itime = time()
            ext = numpy.exp(-1*extLaw.function(ext_lamb,Av=Av,Rv=Rv,f_bump=f_bump))
            results.grid['Av'][self.grid.nrows*count:self.grid.nrows*(count+1)] = Av
            results.grid['Rv'][self.grid.nrows*count:self.grid.nrows*(count+1)] = Rv
            results.grid['f_bump'][self.grid.nrows*count:self.grid.nrows*(count+1)] = f_bump
            results.seds[self.grid.nrows*count:self.grid.nrows*(count+1)] = phot.extractExtinguishedSEDs(self,flist,ext,absFlux=absFlux)
            for key in self.keys():
                results.grid[key][self.grid.nrows*count:self.grid.nrows*(count+1)] = self.grid[key]
            count += 1

        return results

class FileSEDGrid(SpectralGrid):
    """ Generate a grid from a spectral library """

    def __init__(self, fname, *args, **kwargs):
        with pyfits.open(fname) as f:
            self.seds = f[0].data[:-1]
            self.lamb = f[0].data[-1]
        self.grid = mytables.load(fname)
        #lamb, seds = self.getSEDs(filters, self.osl.wavelength, self.osl.spectra)
        for k in self.grid.keys():
            self.__dict__[k] = self.grid[k]

    def keys(self):
        """ returns the grid dimension names """
        return self.grid.keys()


class FileSpectralGrid(SpectralGrid):
    """ Generate a grid from a spectral library """

    def __init__(self, fname, *args, **kwargs):
        with pyfits.open(fname) as f:
            self.seds = f[0].data[:-1]
            self.lamb = f[0].data[-1]
        self.grid = mytables.load(fname)
        #lamb, seds = self.getSEDs(filters, self.osl.wavelength, self.osl.spectra)
        for k in self.grid.keys():
            self.__dict__[k] = self.grid[k]

    def keys(self):
        """ returns the grid dimension names """
        return self.grid.keys()


class StellibGrid(SpectralGrid):
    """ Generate a grid from a spectral library """

    def __init__(self, filters, *args, **kwargs):
        self.osl = stellib.BaSeL()
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


def generate_spectral_grid_from_isochrones(outfile, osl, oiso, Z=0.02):
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

    assert(isNestedInstance(osl, stellib.Stellib) )
    assert(isNestedInstance(oiso, isochrone.Isochrone) )
    specs = numpy.empty( (oiso.data.nrows + 1, len(osl.wavelength)), dtype=float )
    specs[-1] = osl.wavelength[:]

    progress = 0
    with timeit('interpolation'):
        for k in range(oiso.data.nrows):
            if progress < int(100 * (k + 1) / oiso.data.nrows):
                progress = int(100 * (k + 1) / oiso.data.nrows)
                print "progress... %d / 100" % progress
            r = numpy.array( osl.interp(oiso.data['logT'][k], oiso.data['logg'][k], Z, oiso.data['logL'][k]) ).T
            specs[k, :] = osl.genSpectrum(r)
    pyfits.writeto(outfile, specs)

    #copy pars
    data = {}
    for k in oiso.data.keys():
        data[k] = oiso.data[k]
    pars  = mytables.Table(data, name='Reinterpolated stellib grid')
    pars.header['stellib'] = osl.source
    pars.header['isoch'] = oiso.source

    pars.write(outfile, append=True)
