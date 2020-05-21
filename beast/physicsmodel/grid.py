"""
SED/spectral grids
"""
import sys
import numpy as np

from beast.observationmodel import phot
from beast.physicsmodel.dust import extinction
from beast.physicsmodel.helpers.gridbackends import (
    MemoryBackend,
    CacheBackend,
    DiskBackend,
    GridBackend,
)
from beast.physicsmodel.helpers.gridhelpers import pretty_size_print, isNestedInstance

__all__ = ["ModelGrid", "SEDGrid", "SpectralGrid"]


def find_backend(txt):
    """
    Determine the needed background based on a text string

    Parameters
    ----------
    txt : str
        name to find in the list

    Returns
    -------
    b : :class:`~beast.physicsmodel.helpers.gridbackends.GridBackend` subclass
        corresponding backend class
    """

    maps = {
        "memory": MemoryBackend,
        "cache": CacheBackend,
        "disk": DiskBackend,
    }
    btype = maps.get(txt.lower(), None)
    if btype is None:
        raise ValueError(f"{txt} backend not supported")
    return btype


class ModelGrid(object):
    """
    Generic class
    """

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        lamb : ndarray or str or :class:`~beast.physicsmodel.helpers.gridbackends.GridBackend` subclass
            - if ndarray: wavelength of the SEDs (requires seds and
              grid arguments)
            - if str: filename to the grid
            - if backend: ref to the given grid

        seds : ndarray
            2D `float` array of the seds

        grid : :class:`~astropy.table.Table`
            table of properties associated to each sed

        header : dict
            if provided, update the grid table header

        aliases : dict
            if provided, update the grid table aliases

        backend : str or :class:`~beast.physicsmodel.helpers.gridbackends.GridBackend` subclass, optional
            if str corresponding backend class
            'memory' = MemoryBackend,
            'cache' = CacheBackend,
            'disk' = DiskBackend
        """
        backend = kwargs.pop("backend", None)
        if backend is None:
            self._backend = MemoryBackend(*args, **kwargs)
        elif isinstance(backend, (str, bytes)):
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
        self._backend.seds = value

    @property
    def grid(self):
        return self._backend.grid

    @grid.setter
    def grid(self, value):
        self._backend.grid = value

    def __repr__(self):
        txt = "{} ({})"
        return txt.format(object.__repr__(self), pretty_size_print(self.nbytes))

    def __len__(self):
        return self._backend.__len__()

    @property
    def nbytes(self):
        """ The number of bytes of the object """
        n = sum(
            k.nbytes if hasattr(k, "nbytes") else sys.getsizeof(k)
            for k in list(self.__dict__.values())
        )
        return n

    def keys(self):
        """ The grid column names """
        if hasattr(self._backend, "keys"):
            return self._backend.keys()
        else:
            return []

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        elif hasattr(self._backend, name):
            return getattr(self._backend, name)
        elif hasattr(self.grid, "keys"):
            return self.grid[name]
        else:
            msg = "'{0}' object has no attribute '{1}'"
            raise AttributeError(msg.format(type(self).__name__, name))

    def __getitem__(self, name):
        return self.grid[name]

    def copy(self):
        """ returns a copy of the object """
        return self.__class__(backend=self._backend.copy())


class SEDGrid(ModelGrid):
    """
    Generate a grid that the full observational model (SEDs).
    Currently a direct interface to ModelGrid.  Setup for later expansion.

    Attributes
    ----------
    seds : ndarray
        2D `float` array (# models, # bands) giving the seds

    lamb : ndarray
        1D `float` array of the wavelengths of the sed bands

    filters : list
        list of the filter names of the sed bands

    grid : :class:`~astropy.table.Table`
        table with columns providing the model parameters and other
        characteristics of the grid

    header : dict
        header information

    cov_diag, cov_offdiag : ndarray
        2D 'float' arrays with the covariance matrices of the absolute calibration
        uncertainties for each model
    """


class SpectralGrid(ModelGrid):
    """
    Generate a grid that contains spectra.
    It provides an access to integrated photometry function getSEDs.

    Attributes
    ----------
    seds : ndarray
        2D `float` array (# models, # bands) giving the seds

    lamb : ndarray
        1D `float` array of the wavelengths of the sed bands

    grid : :class:`~astropy.table.Table`
        table with columns providing the model parameters and other
        characteristics of the grid

    header : dict
        header information
    """

    def getSEDs(
        self,
        filter_names,
        absFlux=True,
        extLaw=None,
        inplace=False,
        filterLib=None,
        **kwargs
    ):
        """
        Extract integrated fluxes through filters

        Parameters
        ----------
        filter_names : list
            list of filter names according to the filter lib or filter
            instances (no mixing between name and instances)

        absFlux : bool, optional
            returns absolute fluxes if set
            [capability should be removed]

        extLaw : extinction.ExtinctionLaw, optional
            apply extinction law if provided

        inplace : bool, optional
            if set, do not copy the grid and apply extinction on it

        filterLib : str, optional
            full filename to the filter library hd5 file

        **kwargs extra keywords will be forworded to extLaw

        Returns
        -------
        memgrid : :class:`~beast.physicsmodel.helpers.grid.SEDGrid` instance
            grid info with memory backend
        """
        if isinstance(filter_names[0], str):
            flist = phot.load_filters(
                filter_names, interp=True, lamb=self.lamb, filterLib=filterLib
            )
            _fnames = filter_names
        else:
            flist = phot.load_Integrationfilters(
                filter_names, interp=True, lamb=self.lamb
            )
            _fnames = [fk.name for fk in filter_names]
        if extLaw is not None:
            if not inplace:
                r = self.applyExtinctionLaw(extLaw, inplace=inplace, **kwargs)
                lamb, seds, grid = phot.extractSEDs(r, flist, absFlux=absFlux)
            else:
                self.applyExtinctionLaw(extLaw, inplace=inplace, **kwargs)
                lamb, seds, grid = phot.extractSEDs(self, flist, absFlux=absFlux)
        else:
            lamb, seds, grid = phot.extractSEDs(self, flist, absFlux=absFlux)
        memgrid = SEDGrid(lamb, seds, grid, backend=MemoryBackend)

        setattr(memgrid, "filters", _fnames)
        return memgrid

    def applyExtinctionLaw(self, extLaw, inplace=False, **kwargs):
        """
        Apply an extinction law to the model grid

        Parameters
        ----------
        extLaw : extinction.ExtinctionLaw
            apply extinction law if provided

        inplace : bool
            if set, do not copy the grid and apply on it

        **kwargs
            extra keywords will be forwarded to extLaw

        Returns
        -------
        g : :class:`~beast.physicsmodel.helpers.grid.SEDGrid` instance or None
            if not inplace, returns a new :class:`~beast.physicsmodel.helpers.grid.SEDGrid`
            instance. Otherwise returns `None`
        """
        if not isinstance(extLaw, extinction.ExtinctionLaw):
            raise TypeError("Expecting ExtinctionLaw object got %s" % type(extLaw))
        extCurve = np.exp(-1.0 * extLaw.function(self.lamb[:], **kwargs))
        if not inplace:
            g = self.copy()
            g.seds = g.seds[:] * extCurve[None, :]
            g.header["ExtLaw"] = extLaw.name
            for k, v in kwargs.items():
                g.header[k] = v
            return g
        else:
            self.header["ExtLaw"] = extLaw.name
            for k, v in kwargs.items():
                self.header[k] = v
            self.seds = self.seds[:] * extCurve[None, :]
