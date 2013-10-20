""" Manage Various SED/spectral grids is a generic way

Major changes from the previous version of core.grid:
    Removed general write method
    added backend functions and direct access to properties
    MemoryGrid migrates to a function that creates a ModelGrid with a MemoryBackend
    FileSEDGrid, FileSpectralGrid migrated to functions as well

Currently no majors variation is expected as long as memory or cache backend types are used

More optimization can be done, especially in SpectralGrid.getSEDs

TODO: Check where any beast code uses eztable.Table's specific methods and
      implement equivalent in the backends for transparency in the case of HDFBackend
         * aliases
         * eval expression
         * selectWhere
         * readCoordinates (although should work already)
"""
import sys
from copy import deepcopy
import numpy as np

from . import phot
from . import extinction
from .gridbackends import MemoryBackend, CacheBackend, HDFBackend, GridBackend
from .gridhelpers import pretty_size_print


def find_backend(txt):
    """find_backend

    keywords
    --------

    txt: str
        name to find in the list

    returns
    -------

    b: GridBackend class or subclass
        corresponding backend class
    """

    maps = {'memory': MemoryBackend,
            'cache': CacheBackend,
            'hdf': HDFBackend,
            'generic': GridBackend
            }
    return maps.get(txt.lower(), None)


class ModelGrid(object):
    """ Generic class for a minimum update of future codes """
    def __init__(self, *args, **kwargs):
        """
        keywords
        --------
        *args and **kwargs are directly forwarded to the backend constructor

        lamb: ndarray or str or GridBackend
            if ndarray: wavelength of the SEDs (requires seds and grid arguments)
            if str: filename to the grid
            if backend: ref to the given grid

        seds: ndarray[dtype=float, ndim=2]
            array of seds

        grid: eztable.Table
            table of properties associated to each sed

        header: dict
            if provided, update the grid table header

        aliases: dict
            if provided, update the grid table aliases

        backend: str or GridBackend class or subclass
            corresponding backend class

            'memory': MemoryBackend,
            'cache': CacheBackend,
            'hdf': HDFBackend,
            'generic': GridBackend
        """
        backend = kwargs.pop('backend', None)
        if backend is None:
            self._backend = GridBackend(*args, **kwargs)
        elif type(backend) == str:
            self._backend = find_backend(backend)(*args, **kwargs)
        else:
            self._backend = backend(*args, **kwargs)

    @property
    def lamb(self):
        return self._backend.lamb

    @property
    def seds(self):
        return self._backend.seds

    @property
    def grid(self):
        return self._backend.grid

    def __repr__(self):
        txt = '{} ({})'
        return txt.format(object.__repr__(self), pretty_size_print(self.nbytes))

    @property
    def nbytes(self):
        """ return the number of bytes of the object """
        n = sum(k.nbytes if hasattr(k, 'nbytes') else sys.getsizeof(k) for k in self.__dict__.values())
        return n

    def keys(self):
        """ returns the grid dimension names """
        if hasattr(self.grid, 'keys'):
            return self.grid.keys()
        else:
            return []

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        elif hasattr(self._backend, name):
            return getattr(self._backend, name)
        elif name in self.keys():
            return self.grid[name]
        else:
            msg = "'{0}' object has no attribute '{1}'"
            raise AttributeError(msg.format(type(self).__name__, name))

    def copy(self):
        """ returns a copy of the object """
        return deepcopy(self)


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
                lamb, seds, grid = phot.extractSEDs(r, flist, absFlux=absFlux)
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
        extCurve = np.exp(-1. * extLaw.function(self.lamb[:], **kwargs))
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


class StellibGrid(SpectralGrid):
    """ Generate a grid from a spectral library """

    def __init__(self, osl, filters, header={}, aliases={}, *args, **kwargs):
        self.osl = osl
        lamb, seds = self.getSEDs(filters, self.osl.wavelength, self.osl.spectra)
        super(StellibGrid, self).__init__(lamb, seds=seds, grid=self.osl.grid, header=header, aliases=aliases)
        self.filters = filters


def MemoryGrid(lamb, seds=None, grid=None, header={}, aliases={}):
    """ Replace the MemoryGrid class for backwards compatibility

        Instanciate an grid object that has no physical storage
        Helps to create new grids on the fly. Because it deriveds from
        ModelGrid, this can be exported on disk too.

    keywords
    --------

    lamb: ndarray or GridBackend subclass
        if ndarray: wavelength of the SEDs (requires seds and grid arguments)
        if backend: ref to the given grid

    seds: ndarray[dtype=float, ndim=2]
        array of seds

    grid: eztable.Table
        table of properties associated to each sed

    header: dict
        if provided, update the grid table header

    aliases:
        if provided, update the grid table aliases

    returns
    -------
    g: ModelGrid
        grid of models with no physical storage (MemoryBackend)
    """
    return ModelGrid(lamb, seds=seds, grid=grid, header=header,
                     aliases=aliases, backend=MemoryBackend)


def FileSEDGrid(fname, header={}, aliases={}, backend='memory'):
    """ Replace the FileSEDGrid class for backwards compatibility
        Generates a grid from a spectral library on disk

    keywords
    --------

    lamb: ndarray or GridBackend subclass
        if ndarray: wavelength of the SEDs (requires seds and grid arguments)
        if backend: ref to the given grid

    seds: ndarray[dtype=float, ndim=2]
        array of seds

    grid: eztable.Table
        table of properties associated to each sed

    header: dict
        if provided, update the grid table header

    aliases: dict
        if provided, update the grid table aliases

    backend: str or GridBackend class or subclass
        corresponding backend class

        'memory': MemoryBackend,
        'cache': CacheBackend,
        'hdf': HDFBackend,
        'generic': GridBackend

    returns
    -------
    g: ModelGrid
        grid of models with no physical storage (MemoryBackend)
    """
    return SpectralGrid(fname, header=header, aliases=aliases, backend=backend)


# Backward compatibility
FileSpectralGrid = FileSEDGrid
