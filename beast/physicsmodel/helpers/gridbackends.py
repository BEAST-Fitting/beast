"""
Backends to handle the model grids different ways
=================================================

Multiple backends are available to reduce the memory footprint for a
performance cost as small as possible. They allow grids to be stored into FITS
files but also into HDF5 format.

In principle, HDFBackend is the most optimized as it allow you to attack the
grid file directly without memory overhead thanks to the full use of pytables
package.

Implemented Backends
--------------------

MemoryBackend:
    Load everything into memory. Can initiate from variables, a filename,
    CacheBackend, and HDFBackend.

CacheBackend:
    Load data only at the first request. You can work using only seds or only
    model properties without the overhead of loading both. (works with FITS and
    HDF files). Offers also dropping part of the data.

HDFBackend:
    Works directly with an HDFStore support, ie., on disk. Cache and reading
    are allowed through any way offered by pytables, which becomes very handy
    for very low-memory tasks such as doing single star figures.

All backends are able to write on disk into FITS and HDF format.

TODO: add evalexpr into the HDFBackend grid
        needs a evalexpr into the backends directly to offer transparent access

TODO: add readCoordinates into all backends

TODO: check read(field=) exists into all backends.grid, give direct access
"""
from __future__ import (absolute_import, division, print_function)

import sys
import numpy
import astropy.io.fits as pyfits
import copy

from ...external.eztables import Table
from .hdfstore import HDFStore
from .gridhelpers import isNestedInstance, pretty_size_print

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
    basestring = (str, str)

__all__ = ['GridBackend', 'MemoryBackend', 'CacheBackend', 'HDFBackend']


class GridBackend(object):
    """GridBackend
    How the content of a grid is handled. The idea is to provide enough
    flexibility that low-memory footprint can be achieved if needed

    This class is a generic implementation that will be derived into
    more classes
    """
    def __init__(self, *args, **kwargs):
        self._filters = None
        self._header = None
        self.fname = None
        self._aliases = {}

    @property
    def nbytes(self):
        """ return the number of bytes of the object """
        n = sum(k.nbytes if hasattr(k, 'nbytes') else sys.getsizeof(k) \
                for k in list(self.__dict__.values()))
        return n

    @property
    def header(self):
        """header"""
        return self._header

    @property
    def filters(self):
        """filters"""
        return self._filters

    def __len__(self):
        return len(self.grid)

    def keys(self):
        """ returns the grid dimension names """
        if hasattr(self.grid, 'keys'):
            return list(self.grid.keys())
        else:
            return []

    def _get_type(self, fname):
        """_get_type -- guess the type of the file fname
        """
        types = {'fits': 'fits',
                 'hdf': 'hdf',
                 'hd5': 'hdf',
                 'hdf5': 'hdf5'
                 }
        #if fname.split('.')[-1] not in types:
        #    print(fname)
        #    try:
        #        hdulist = pyfits.open(fname)
        #        rtype = 'fits'
        #    except:
        #        print('An error occured trying to read the file.')            
        #else:
        rtype = types[fname.split('.')[-1]]
        return rtype

    def __repr__(self):
        """__repr__"""
        txt = '{}\n source: {}, \n current memory footprint: {}'
        return txt.format(object.__repr__(self),
                          self.fname, pretty_size_print(self.nbytes))

    def _from_HDFBackend(self, b):
        """_from_HDFBackend -- convert from HDFBackend

        Parameters
        ----------

        b: GridBackend or sub class
            backend to convert from
        """
        self.lamb = b.lamb.read()
        self.seds = b.seds.read()
        self.grid = Table(b.grid.read())
        self._filters = b._filters[:]
        self._header = b.header
        self._aliases = b._aliases

    def _from_GridBackend(self, b):
        """_from_GridBackend -- convert from generic backend

        Parameters
        ----------

        b: GridBackend or sub class
            backend to convert from
        """
        self.lamb = b.lamb
        self.seds = b.seds
        self.grid = b.grid
        self._filters = b._filters
        self._header = b.header
        self._aliases = b._aliases

    def copy(self):
        """ implement a copy method """
        g = GridBackend()
        g.lamb = copy.deepcopy(self.lamb)
        g.seds = copy.deepcopy(self.seds)
        g.grid = copy.deepcopy(self.grid)
        g._filters = copy.deepcopy(self._filters)
        g._header = copy.deepcopy(self._header)
        g._aliases = copy.deepcopy(self._aliases)
        return g


class MemoryBackend(GridBackend):
    """ Instanciate an grid object that has no physical storage
        Helps to create new grids on the fly. Because it deriveds from
        ModelGrid, this can be exported on disk too.
    """
    def __init__(self, lamb, seds=None, grid=None,
                 cov_diag=None, cov_offdiag=None,
                 header={}, aliases={}):
        """__init__

        Parameters
        ----------

        lamb: ndarray or GridBackend subclass
            if ndarray: wavelength of the SEDs (requires seds and grid
                                                arguments)
            if backend: ref to the given grid

        seds: ndarray[dtype=float, ndim=2]
            array of seds

        grid: eztable.Table
            table of properties associated to each sed

        header: dict
            if provided, update the grid table header

        aliases:
            if provided, update the grid table aliases

        """
        super(MemoryBackend, self).__init__()

        #read from various formats
        if isinstance(lamb, HDFBackend):
            self._fromHDFBackend(lamb)
        elif isNestedInstance(lamb, GridBackend):
            self._from_GridBackend(lamb)
        elif type(lamb) in basestring:
            self._from_File(lamb)
        else:
            if ((seds is None) | (grid is None)):
                raise ValueError('Wrong number of arguments')
            self.lamb = lamb
            self.seds = seds
            self.grid = grid
            if (cov_diag is not None) & (cov_offdiag is not None):
                print('including cov diag and offdiag')
                self.cov_diag = cov_diag
                self.cov_offdiag = cov_offdiag
            else:
                self.cov_diag = None
                self.cov_offdiag = None

        #update header
        if self._header is None:
            self._header = header
        else:
            for k, v in list(header.items()):
                self.grid.header[k] = v

        #update aliases
        self._aliases.update(aliases)
        self.fname = ':memory:'

    @property
    def filters(self):
        """filters"""
        r = self._header.get('filters', None) or self._header.get('FILTERS',
                                                                  None)
        if r is not None:
            r = r.split()
        return r

    @property
    def header(self):
        return self._header

    def _from_File(self, fname):
        """_from_File -- load the content of a FITS or HDF file

        Parameters
        ----------

        fname: str
            filename (incl. path) to read from
        """

        # load_seds - load wavelength and seds
        if self._get_type(fname) == 'fits':
            with pyfits.open(fname) as f:
                self.seds = f[0].data[:-1]
                self.lamb = f[0].data[-1]
            self.grid = Table(fname)

        elif self._get_type(fname) == 'hdf':
            with HDFStore(fname, mode='r') as s:
                self.seds = s['/seds'].read()
                self.lamb = s['/lamb'].read()
                try:
                    self.cov_diag = s['/covdiag'].read()
                except:
                    self.cov_diag = None
                try:
                    self.cov_offdiag = s['/covoffdiag'].read()
                except:
                    self.cov_offdiag = None
            self.grid = Table(fname, tablename='/grid')

        self._header = self.grid.header

    def writeFITS(self, fname, *args, **kwargs):
        """write -- export to fits file

        Parameters
        ----------

        fname: str
            filename (incl. path) to export to
        """
        if ( (self.lamb is not None) & (self.seds is not None) &
             (self.grid is not None) ):
            assert(isinstance(self.grid, Table)), \
                                    'Only eztables.Table are supported so far'
            r = numpy.vstack( [ self.seds, self.lamb ])
            pyfits.writeto(fname, r, **kwargs)
            if getattr(self, 'filters', None) is not None:
                if ('FILTERS' not in list(self.grid.header.keys())):
                    self.grid.header['FILTERS'] = ' '.join(self.filters)
            self.grid.write(fname, append=True)

    def writeHDF(self, fname, append=False, *args, **kwargs):
        """write -- export to HDF file

        Parameters
        ----------

        fname: str
            filename (incl. path) to export to

        append: bool, optional (default False)
            if set, it will append data to each Array or Table
        """
        if ( (self.lamb is not None) & (self.seds is not None) &
             (self.grid is not None) ):
            assert(isinstance(self.grid, Table)), \
                                    'Only eztables.Table are supported so far'
            with HDFStore(fname, mode='a') as hd:
                if not append:
                    hd['/seds'] = self.seds[:]
                    hd['/lamb'] = self.lamb[:]
                    if self.cov_diag is not None:
                        hd['/covdiag'] = self.cov_diag[:]
                    if self.cov_offdiag is not None:
                        hd['/covoffdiag'] = self.cov_offdiag[:]
                else:
                    try:
                        node = hd.get_node('/seds')
                        node.append(self.seds[:])
                    except:
                        hd['/seds'] = self.seds[:]
                        hd['/lamb'] = self.lamb[:]
                        if self.cov_diag is not None:
                            hd['/covdiag'] = self.cov_diag[:]
                        if self.cov_offdiag is not None:
                            hd['/covoffdiag'] = self.cov_offdiag[:]
            if getattr(self, 'filters', None) is not None:
                if ('FILTERS' not in list(self.grid.header.keys())):
                    self.grid.header['FILTERS'] = ' '.join(self.filters)
            self.grid.write(fname, tablename='grid', append=True)

    def copy(self):
        """ implement a copy method """
        g = MemoryBackend(copy.deepcopy(self.lamb),
                          seds=copy.deepcopy(self.seds),
                          grid=copy.deepcopy(self.grid),
                          cov_diag=copy.deepcopy(self.cov_diag),
                          cov_offdiag=copy.deepcopy(self.cov_offdiag))
        g._filters = copy.deepcopy(self._filters)
        g._header = copy.deepcopy(self._header)
        g._aliases = copy.deepcopy(self._aliases)
        return g


class CacheBackend(GridBackend):
    """CacheBackend -- Load content from a file only when needed

    The key idea is to be able to load the content only at the first query

    Currently the grid attribute is an eztable.Table object as it was before.
    """

    def __init__(self, fname, *args, **kwargs):
        """__init__

        Parameters
        ----------

        fname: str
            FITS or HD5 file containing the grid
        """
        super(CacheBackend, self).__init__()

        self.fname = fname
        self._type = self._get_type(fname)
        self.clear()

    def clear(self, attrname=None):
        """clear current cache

        Parameters
        ----------

        attrname: str in [lamb, filters, grid, header, lamb, seds]
            if provided clear only one attribute
            else all cache will be erased
        """
        if attrname is None:
            self._seds = None
            self._lamb = None
            self._filters = None
            self._grid = None
            self._header = None
        else:
            setattr(self, '_{0}'.format(attrname), None)

    def _load_seds(self, fname):
        """load_seds - load seds"""
        if (self._seds is None):
            if self._get_type(fname) == 'fits':
                with pyfits.open(self.fname) as f:
                    self._seds = f[0].data[:-1]

            elif self._get_type(fname) == 'hdf':
                with HDFStore(self.fname, mode='r') as s:
                    self._seds = s['/seds'].read()

    def _load_lamb(self, fname):
        """load_seds - load wavelength"""
        if self._lamb is None:
            if self._get_type(fname) == 'fits':
                with pyfits.open(self.fname) as f:
                    self._lamb = f[0].data[-1]

            elif self._get_type(fname) == 'hdf':
                with HDFStore(self.fname, mode='r') as s:
                    self._lamb = s['/lamb'].read()

    def _load_grid(self, fname):
        """load_grid - load grid table"""
        # load_seds - load wavelength and seds
        if (self._grid is None):
            if self._get_type(fname) == 'fits':
                self._grid = Table(self.fname)

            elif self._get_type(fname) == 'hdf':
                self._grid = Table(self.fname, tablename='/grid')

    def _load_filters(self, fname):
        """load_filters -- load only filters"""
        if self._filters is None:
            if self._type == 'fits':
                with pyfits.open(self.fname) as f:
                    self._filters = f[1].header.get('FILTERS', None) or f[1].header.get('filters', None)
                    if self._filters is not None:
                        self._filters = self._filters.split()
            elif self._type == 'hdf':
                self._filters = self.header.get('FILTERS', None) or self.header.get('filters', None)
                if self._filters is not None:
                    self._filters = self._filters.split()

    @property
    def seds(self):
        """seds - load in cache if needed """
        self._load_seds(self.fname)
        return self._seds

    @seds.setter
    def seds(self, value):
        """ replace seds value """
        self._seds = value

    @property
    def lamb(self):
        """lamb - load in cache if needed """
        self._load_lamb(self.fname)
        return self._lamb

    @lamb.setter
    def lamb(self, value):
        """ replace seds value """
        self._lamb = value

    @property
    def grid(self):
        """grid - load in cache if needed """
        self._load_grid(self.fname)
        return self._grid

    @grid.setter
    def grid(self, value):
        """ replace seds value """
        self._grid = value

    @property
    def header(self):
        """header - load in cache if needed """
        self._load_grid(self.fname)
        return self._grid.header

    @header.setter
    def header(self, value):
        """ replace seds value """
        self._header = value

    @property
    def filters(self):
        """filters - load in cache if needed """
        self._load_filters(self.fname)
        return self._filters

    @filters.setter
    def filters(self, value):
        """ replace seds value """
        self._filters = value

    def keys(self):
        """ return column names when possible, avoid loading when possible """
        if hasattr(self._grid, 'coldescrs'):
            return list(self._grid.coldescrs.keys())
        elif hasattr(self._grid, 'keys'):
            return list(self._grid.keys())
        elif hasattr(self.grid, 'keys'):
            return list(self.grid.keys())
        else:
            return []

    def writeFITS(self, fname, *args, **kwargs):
        """write -- export to fits file

        Parameters
        ----------

        fname: str
            filename (incl. path) to export to
        """
        if ( (self.lamb is not None) & (self.seds is not None) & (self.grid is not None) ):
            assert(isinstance(self.grid, Table)), 'Only eztables.Table are supported so far'
            r = numpy.vstack( [ self.seds, self.lamb ])
            pyfits.writeto(fname, r, **kwargs)
            del r
            if getattr(self, 'filters', None) is not None:
                if ('FILTERS' not in list(self.grid.header.keys())):
                    self.grid.header['FILTERS'] = ' '.join(self.filters)
            self.grid.write(fname, append=True)

    def writeHDF(self, fname, append=False, *args, **kwargs):
        """write -- export to HDF file

        Parameters
        ----------

        fname: str
            filename (incl. path) to export to

        append: bool, optional (default False)
            if set, it will append data to each Array or Table
        """
        if ( (self.lamb is not None) & (self.seds is not None) & (self.grid is not None) ):
            assert(isinstance(self.grid, Table)), 'Only eztables.Table are supported so far'
            with HDFStore(fname, mode='a') as hd:
                if not append:
                    hd['/seds'] = self.seds[:]
                    hd['/lamb'] = self.lamb[:]
                else:
                    try:
                        node = hd.get_node('/seds')
                        node.append(self.seds[:])
                    except:
                        hd['/seds'] = self.seds[:]
                        hd['/lamb'] = self.lamb[:]
            if getattr(self, 'filters', None) is not None:
                if ('FILTERS' not in list(self.grid.header.keys())):
                    self.grid.header['FILTERS'] = ' '.join(self.filters)
            self.grid.write(fname, tablename='grid', append=True)

    def copy(self):
        """ implement a copy method """
        g = CacheBackend(self.fname)
        g._aliases = copy.deepcopy(self._aliases)
        if (self._grid is not None):
            g._grid = copy.deepcopy(self._grid)
        if (self._seds is not None):
            g._seds = copy.deepcopy(self._seds)
        if (self._lamb is not None):
            g._lamb = copy.deepcopy(self._lamb)
        if (self._header is not None):
            g._header = copy.deepcopy(self._header)
        if (self._filters is not None):
            g._filters = copy.deepcopy(self._filters)

        return g


class HDFBackend(GridBackend):
    """HDFBackend -- Laziest grid
    Operations are optimized on disk through pytables
    """
    def __init__(self, fname, *args, **kwargs):
        super(HDFBackend, self).__init__()
        ftype = self._get_type(fname)
        if ftype != 'hdf':
            raise ValueError('Expecting HDF file got {0}'.format(ftype))

        self.fname = fname
        self.store = HDFStore(self.fname, mode='r')
        self.seds = self.store['/seds']
        self.lamb = self.store['/lamb']
        self.grid = self.store['/grid']
        self._filters = None
        self._header = None
        self._aliases = {}

    @property
    def header(self):
        if self._header is None:
            #update header & aliases
            exclude = ['NROWS', 'VERSION', 'CLASS', 'EXTNAME']
            header = {}
            for k in self.grid.attrs._v_attrnames:
                if (not k in exclude) & (k[:5] != 'FIELD') & (k[:5] != 'ALIAS'):
                    header[k] = self.grid.attrs[k]
                if (k[:5] == 'ALIAS'):
                    c0, c1 = self.grid.attrs[k].split('=')
                    self._aliases[c0] = c1

            empty_name = ['', 'None', 'Noname', None]
            if (header['NAME'] in empty_name) & (header.get('TITLE', None) not in empty_name):
                header['NAME'] = header['TITLE']
            self._header = header
        return self._header

    @property
    def filters(self):
        """filters - load in cache if needed """
        if self._filters is None:
                self._filters = self.header.get('FILTERS', None) or self.header.get('filters', None)
                if self._filters is not None:
                    self._filters = self._filters.split()
        return self._filters

    def keys(self):
        """ return column names when possible """
        if hasattr(self._grid, 'coldescrs'):
            return list(self._grid.coldescrs.keys())
        else:
            return []

    def writeHDF(self, fname, append=False, *args, **kwargs):
        """write -- export to HDF file

        Parameters
        ---------

        fname: str
            filename (incl. path) to export to

        append: bool, optional (default False)
            if set, it will append data to each Array or Table
        """
        if ( (self.lamb is not None) & (self.seds is not None) & (self.grid is not None) ):
            with HDFStore(fname, mode='a') as hd:
                if not append:
                    hd['/seds'] = self.seds[:]
                    hd['/lamb'] = self.lamb[:]
                else:
                    try:
                        node = hd.get_node('/seds')
                        node.append(self.seds[:])
                    except:
                        hd['/seds'] = self.seds[:]
                        hd['/lamb'] = self.lamb[:]
                hd.write(self.grid[:], group='/', tablename='grid', header=self.header, append=append)

    def copy(self):
        g = HDFBackend(self.fname)
        g._aliases = copy.deepcopy(self._aliases)
        return g
