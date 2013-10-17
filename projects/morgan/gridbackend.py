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
"""
import sys
import numpy
import pyfits

from beast.external.eztables import Table
from hdfstore import HDFStore


__all__ = ['GridBackend', 'MemoryBackend', 'CacheBackend', 'HDFBackend']


def isNestedInstance(obj, cl):
    """ Test for sub-classes types
        I could not find a universal test

        keywords
        --------
        obj: object instance
            object to test

        cl: Class
            top level class to test

        returns
        -------
        r: bool
            True if obj is indeed an instance or subclass instance of cl
    """
    tree = []
    for k in cl.__subclasses__():
        tree += k.__subclasses__()
    tree += cl.__subclasses__() + [ cl ]
    return issubclass(obj.__class__, tuple(tree))


def pretty_size_print(num_bytes):
    """
    Output number of bytes in a human readable format

    keywords
    --------
    num_bytes: int
        number of bytes to convert

    returns
    -------
    output: str
        string representation of the size with appropriate unit scale
    """
    if num_bytes is None:
        return

    KiB = 1024
    MiB = KiB * KiB
    GiB = KiB * MiB
    TiB = KiB * GiB
    PiB = KiB * TiB
    EiB = KiB * PiB
    ZiB = KiB * EiB
    YiB = KiB * ZiB

    if num_bytes > YiB:
        output = '%.3g YB' % (num_bytes / YiB)
    elif num_bytes > ZiB:
        output = '%.3g ZB' % (num_bytes / ZiB)
    elif num_bytes > EiB:
        output = '%.3g EB' % (num_bytes / EiB)
    elif num_bytes > PiB:
        output = '%.3g PB' % (num_bytes / PiB)
    elif num_bytes > TiB:
        output = '%.3g TB' % (num_bytes / TiB)
    elif num_bytes > GiB:
        output = '%.3g GB' % (num_bytes / GiB)
    elif num_bytes > MiB:
        output = '%.3g MB' % (num_bytes / MiB)
    elif num_bytes > KiB:
        output = '%.3g KB' % (num_bytes / KiB)
    else:
        output = '%.3g Bytes' % (num_bytes)

    return output


class GridBackend(object):
    """GridBackend
    How the content of a grid is handled. The idea is to provide enough
    flexibility that low-memory footprint can be achieved if needed

    This class is a generic implementation that will be derived into more classes
    """
    def __init__(self, *args, **kwargs):
        self._filters = None
        self._header = None
        self._aliases = {}

    @property
    def nbytes(self):
        """ return the number of bytes of the object """
        n = sum(k.nbytes if hasattr(k, 'nbytes') else sys.getsizeof(k) for k in self.__dict__.values())
        return n

    @property
    def header(self):
        """header"""
        return self._header

    @property
    def filters(self):
        """filters"""
        return self._filters

    def _get_type(self, fname):
        """_get_type -- guess the type of the file fname
        """
        types = {'fits': 'fits',
                 'hdf': 'hdf',
                 'hd5': 'hdf',
                 'hdf5': 'hdf5'
                 }
        return types[fname.split('.')[-1]]

    def __repr__(self):
        """__repr__"""
        txt = '{}\n source: {}, \n current memory footprint: {})'
        return txt.format(object.__repr__(self), self.fname, pretty_size_print(self.nbytes))

    def _from_HDFBackend(self, b):
        """_from_HDFBackend -- convert from HDFBackend

        keywords
        --------

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

        keywords
        --------

        b: GridBackend or sub class
            backend to convert from
        """
        self.lamb = b.lamb
        self.seds = b.seds
        self.grid = b.grid
        self._filters = b._filters
        self._header = b.header
        self._aliases = b._aliases


class MemoryBackend(GridBackend):
    """ Instanciate an grid object that has no physical storage
        Helps to create new grids on the fly. Because it deriveds from
        ModelGrid, this can be exported on disk too.
    """
    def __init__(self, lamb, seds=None, grid=None, header={}, aliases={}):
        """__init__

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

        """
        super(MemoryBackend, self).__init__()

        #read from various formats
        if isinstance(lamb, HDFBackend):
            self._fromHDFBackend(lamb)
        elif isNestedInstance(lamb, GridBackend):
            self._from_GridBackend(lamb)
        elif type(lamb) == str:
            self._from_File(lamb)
        else:
            if ((seds is None) | (grid is None)):
                raise ValueError('Wrong number of arguments')
            self.lamb = lamb
            self.seds = seds
            self.grid = grid

        #update header
        if self._header is None:
            self._header = header
        else:
            for k, v in header.items():
                self.grid.header[k] = v

        #update aliases
        self._aliases.update(aliases)
        self.fname = ':memory:'

    def _from_File(self, fname):
        """_from_File -- load the content of a FITS or HDF file

        keywords
        --------

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
            self.grid = Table(fname, tablename='/grid')

        self._header = self.grid.header

    def writeFITS(self, fname, *args, **kwargs):
        """write -- export to fits file

        keywords
        --------

        fname: str
            filename (incl. path) to export to
        """
        if ( (self.lamb is not None) & (self.seds is not None) & (self.grid is not None) ):
            assert(isinstance(self.grid, Table)), 'Only eztables.Table are supported so far'
            r = numpy.vstack( [ self.seds, self.lamb ])
            pyfits.writeto(fname, r, **kwargs)
            if getattr(self, 'filters', None) is not None:
                if ('FILTERS' not in self.grid.header.keys()):
                    self.grid.header['FILTERS'] = ' '.join(self.filters)
            self.grid.write(fname, append=True)

    def writeHDF(self, fname, *args, **kwargs):
        """write -- export to HDF file

        keywords
        --------

        fname: str
            filename (incl. path) to export to
        """
        if ( (self.lamb is not None) & (self.seds is not None) & (self.grid is not None) ):
            assert(isinstance(self.grid, Table)), 'Only eztables.Table are supported so far'
            with HDFStore(fname, mode='a') as hd:
                hd['/seds'] = self.seds
                hd['/lamb'] = self.lamb
            if getattr(self, 'filters', None) is not None:
                if ('FILTERS' not in self.grid.header.keys()):
                    self.grid.header['FILTERS'] = ' '.join(self.filters)
            self.grid.write(fname, tablename='grid', append=True)


class CacheBackend(GridBackend):
    """CacheBackend -- Load content from a file only when needed

    The key idea is to be able to load the content only at the first query

    Currently the grid attribute is an eztable.Table object as it was before.
    """

    def __init__(self, fname, *args, **kwargs):
        super(CacheBackend, self).__init__()

        self.fname = fname
        self._type = self._get_type(fname)
        self.clear()

    def clear(self, attrname=None):
        """clear current cache

        keywords
        --------

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
            setattr(self, '_{}'.format(attrname), None)

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

    @property
    def lamb(self):
        """lamb - load in cache if needed """
        self._load_lamb(self.fname)
        return self._lamb

    @property
    def grid(self):
        """grid - load in cache if needed """
        self._load_grid(self.fname)
        return self._grid

    @property
    def header(self):
        """header - load in cache if needed """
        self._load_grid(self.fname)
        return self._grid.header

    @property
    def filters(self):
        """filters - load in cache if needed """
        self._load_filters(self.fname)
        return self._filters

    def writeFITS(self, fname, *args, **kwargs):
        """write -- export to fits file

        keywords
        --------

        fname: str
            filename (incl. path) to export to
        """
        if ( (self.lamb is not None) & (self.seds is not None) & (self.grid is not None) ):
            assert(isinstance(self.grid, Table)), 'Only eztables.Table are supported so far'
            r = numpy.vstack( [ self.seds, self.lamb ])
            pyfits.writeto(fname, r, **kwargs)
            del r
            if getattr(self, 'filters', None) is not None:
                if ('FILTERS' not in self.grid.header.keys()):
                    self.grid.header['FILTERS'] = ' '.join(self.filters)
            self.grid.write(fname, append=True)

    def writeHDF(self, fname, *args, **kwargs):
        """write -- export to HDF file

        keywords
        --------

        fname: str
            filename (incl. path) to export to
        """
        if ( (self.lamb is not None) & (self.seds is not None) & (self.grid is not None) ):
            assert(isinstance(self.grid, Table)), 'Only eztables.Table are supported so far'
            with HDFStore(fname, mode='a') as hd:
                hd['/seds'] = self.seds
                hd['/lamb'] = self.lamb
            if getattr(self, 'filters', None) is not None:
                if ('FILTERS' not in self.grid.header.keys()):
                    self.grid.header['FILTERS'] = ' '.join(self.filters)
            self.grid.write(fname, tablename='grid', append=True)


class HDFBackend(GridBackend):
    """HDFBackend -- Laziest grid
    Operations are optimized on disk through pytables
    """
    def __init__(self, fname, *args, **kwargs):
        super(HDFBackend, self).__init__()
        ftype = self._get_type(fname)
        if ftype != 'hdf':
            raise ValueError('Expecting HDF file got {}'.format(ftype))

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

    def writeHDF(self, fname, *args, **kwargs):
        """write -- export to HDF file

        keywords
        --------

        fname: str
            filename (incl. path) to export to
        """
        if ( (self.lamb is not None) & (self.seds is not None) & (self.grid is not None) ):
            assert(isinstance(self.grid, Table)), 'Only eztables.Table are supported so far'
            with HDFStore(fname, mode='a') as hd:
                hd['/seds'] = self.seds[:]
                hd['/lamb'] = self.lamb[:]
                hd.write(self.grid[:], group='/', tablename='grid', header=self.header)
