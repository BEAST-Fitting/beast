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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import numpy as np
from copy import deepcopy

from ..observationmodel import phot
from .dust import extinction
from .helpers.gridbackends import MemoryBackend, CacheBackend, HDFBackend, GridBackend
from .helpers.gridhelpers import pretty_size_print, isNestedInstance

try:
    unicode = unicode
except NameError:
    # 'unicode' is undefined, must be Python 3
    str = str
    unicode = str
    bytes = bytes
    basestring = (str,bytes)
else:
    # 'unicode' exists, must be Python 2
    str = str
    unicode = unicode
    bytes = str
    basestring = (str, unicode)

__all__ = ['ModelGrid', 'SpectralGrid','StellibGrid',
           'MemoryGrid','FileSEDGrid']

def find_backend(txt):
    """find_backend

    Parameters
    ----------

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
        Parameters
        ----------

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
        elif type(backend) in basestring:
            self._backend = find_backend(backend)(*args, **kwargs)
        elif isNestedInstance(backend, GridBackend):
            self._backend = backend
        else:
            self._backend = backend(*args, **kwargs)

    @property
    def lamb(self):
        return self._backend.lamb

    @lamb.setter
    def lamb(self, value):
        """ Allow temporary overriding properties """
        self._backend.lamb = value

    @property
    def header(self):
        return self._backend.header

    @header.setter
    def header(self, value):
        self._backend.header = value

    @property
    def seds(self):
        return self._backend.seds

    @seds.setter
    def seds(self, value):
        """ Allow temporary overriding properties """
        self._backend.seds = value

    @property
    def grid(self):
        return self._backend.grid

    @grid.setter
    def grid(self, value):
        """ Allow temporary overriding properties """
        self._backend.grid = value

    def __repr__(self):
        txt = '{} ({})'
        return txt.format(object.__repr__(self),
                          pretty_size_print(self.nbytes))

    @property
    def nbytes(self):
        """ return the number of bytes of the object """
        n = sum(k.nbytes if hasattr(k, 'nbytes') else \
                sys.getsizeof(k) for k in list(self.__dict__.values()))
        return n

    def keys(self):
        """ returns the grid dimension names """
        if hasattr(self.grid, 'keys'):
            return list(self.grid.keys())
        elif hasattr(self.grid, 'colnames'):
            return self.grid.colnames
        else:
            return []

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        elif hasattr(self._backend, name):
            return getattr(self._backend, name)
        #elif name in self.keys():
        #    return self.grid[name]
        elif hasattr(self.grid, 'keys'):
            return self.grid[name]
        else:
            msg = "'{0}' object has no attribute '{1}'"
            raise AttributeError(msg.format(type(self).__name__, name))

    def __getitem__(self, name):
        if hasattr(self.grid, 'read'):
            try:
                return self.grid.read(field=name)
            except TypeError:
                return self.grid[name]
        else:
            return self.grid[name]

    def copy(self):
        """ returns a copy of the object """
        return self.__class__(backend=self._backend.copy())


class SpectralGrid(ModelGrid):
    """ Generate a grid that contains spectra.
    It provides an access to integrated photometry function getSEDs """

    def getSEDs(self, filter_names, absFlux=True, extLaw=None, inplace=False,
                filterLib=None, **kwargs):
        """
        Extract integrated fluxes through filters

        Parameters
        ----------
        filter_names: list
            list of filter names according to the filter lib or filter
            instances
            (no mixing between name and instances)

        absFlux:bool
            returns absolute fluxes if set

        extLaw: extinction.ExtinctionLaw
            apply extinction law if provided

        inplace:bool
            if set, do not copy the grid and apply extinction on it

        filterLib:  str
            full filename to the filter library hd5 file

        **kwargs extra keywords will be forwrded to extLaw

        Returns
        -------
        memgrid: MemoryGrid instance
            grid of SEDs
        """
        if type(filter_names[0]) == str:
            flist = phot.load_filters(filter_names, interp=True,
                                      lamb=self.lamb, filterLib=filterLib)
            _fnames = filter_names
        else:
            flist = phot.load_Integrationfilters(filter_names, interp=True,
                                                 lamb=self.lamb)
            _fnames = [ fk.name for fk in filter_names ]
        if extLaw is not None:
            if not inplace:
                r = self.applyExtinctionLaw(extLaw, inplace=inplace, **kwargs)
                lamb, seds, grid = phot.extractSEDs(r, flist, absFlux=absFlux)
            else:
                self.applyExtinctionLaw(extLaw, inplace=inplace, **kwargs)
                r = self
                lamb, seds, grid = phot.extractSEDs(self, flist,
                                                    absFlux=absFlux)
        else:
            r = self
            lamb, seds, grid = phot.extractSEDs(self, flist, absFlux=absFlux)
        memgrid = MemoryGrid(lamb, seds, grid)
        setattr(memgrid, 'filters', _fnames)
        return memgrid

    def applyExtinctionLaw(self, extLaw, inplace=False, **kwargs):
        """
        Apply an extinction law to the model grid

        Parameters
        ----------
        extLaw: extinction.ExtinctionLaw
            apply extinction law if provided

        inplace: bool
            if set, do not copy the grid and apply on it

        **kwargs        extra keywords will be forwrded to extLaw

        Return
        ------
        g: ModelGrid instance or None
            if not inplace, returns a new ModelGrid instance. Otherwise returns
            nothing
        """
        assert( isinstance(extLaw, extinction.ExtinctionLaw)), \
                'Expecting ExtinctionLaw object got %s' % type(extLaw)
        extCurve = np.exp(-1. * extLaw.function(self.lamb[:], **kwargs))
        if not inplace:
            g = self.copy()
            g.seds = g.seds[:] * extCurve[None, :]
            g.header['ExtLaw'] = extLaw.name
            for k, v in kwargs.items():
                g.header[k] = v
            return g
        else:
            self.header['ExtLaw'] = extLaw.name
            for k, v in kwargs.items():
                self.header[k] = v
            self.seds = self.seds[:] * extCurve[None, :]


class StellibGrid(SpectralGrid):
    """ Generate a grid from a spectral library """

    def __init__(self, osl, filters, header={}, aliases={}, *args, **kwargs):
        self.osl = osl
        lamb, seds = self.getSEDs(filters,
                                  self.osl.wavelength,
                                  self.osl.spectra)
        super(StellibGrid, self).__init__(lamb,
                                          seds=seds,
                                          grid=self.osl.grid,
                                          header=header,
                                          aliases=aliases,
                                          backend=MemoryBackend)
        self.filters = filters

    def copy(self):
        g = super(StellibGrid, self).copy()
        g.osl = deepcopy(self.osl)


def MemoryGrid(lamb, seds=None, grid=None, header={}, aliases={}):
    """ Replace the MemoryGrid class for backwards compatibility

        Instanciate an grid object that has no physical storage
        Helps to create new grids on the fly. Because it deriveds from
        ModelGrid, this can be exported on disk too.

    Parameters
    ----------

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

    Parameters
    ----------

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
