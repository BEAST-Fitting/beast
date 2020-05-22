"""
Backends to handle the model grids different ways
=================================================

Multiple backends are available to reduce the memory footprint for a
performance cost as small as possible.

Implemented Backends
--------------------

MemoryBackend:
    Load everything into memory. Can initiate from variables, a filename,
    CacheBackend, and DiskBackend.

CacheBackend:
    Load data only at the first request. You can work using only seds or only
    model properties without the overhead of loading both. (works with FITS and
    HDF files). Offers also dropping part of the data.

DiskBackend:
    Works directly with an h5py support, ie., on disk. Cache and reading
    are allowed through any way offered by h5py, which becomes very handy
    for very low-memory tasks such as doing single star figures.
"""
import sys
from astropy.io import fits
import h5py
import copy
from astropy.table import Table
from astropy.io.misc.hdf5 import read_table_hdf5

from beast.physicsmodel.helpers.gridhelpers import pretty_size_print, isNestedInstance

__all__ = ["GridBackend", "MemoryBackend", "CacheBackend", "DiskBackend"]


def _decodebytestring(a):
    """
    Convert to string if input is a bytestring.

    Parameters
    ----------
    a : byte or str
        string or bytestring

    Returns
    -------
    str
        string version of input
    """
    if isinstance(a, bytes):
        return a.decode()
    else:
        return a


def _gethdfdatasetmeta(hdfds):
    """
    Extract the meta(header) information from the grid dataset in a hdf file.
    Done without reading the entire grid into memory.

    Parameters
    ----------
    hdfds : h5py dataset
        the hdf dataset

    Returns
    -------
    header : dict
        dictionary of header information
    """
    exclude = ["NROWS", "VERSION", "CLASS", "EXTNAME"]
    header = {}
    attrs = hdfds.attrs
    for k in hdfds.attrs.keys():
        if (k not in exclude) & (k[:5] != "FIELD") & (k[:5] != "ALIAS"):
            header[k] = _decodebytestring(attrs[k])
    return header


class GridBackend(object):
    """
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
        n = sum(
            k.nbytes if hasattr(k, "nbytes") else sys.getsizeof(k)
            for k in list(self.__dict__.values())
        )
        return n

    @property
    def header(self):
        return self._header

    def __len__(self):
        """ number of models in grid """
        return len(self.grid)

    def keys(self):
        """ returns the grid keys"""
        if hasattr(self.grid, "keys"):
            return list(self.grid.keys())
        else:
            return []

    def _get_type(self, fname):
        """ determine the type of the file fname
        """
        types = {"fits": "fits", "hdf": "hdf", "hd5": "hdf", "hdf5": "hdf"}
        # else:
        fext = fname.split(".")[-1]
        if fext not in types:
            raise ValueError(f"{fext} file type not supported")
        else:
            return types[fext]

    def __repr__(self):
        """ print the object and memory usage"""
        txt = "{}\n source: {}, \n current memory footprint: {}"
        return txt.format(
            object.__repr__(self), self.fname, pretty_size_print(self.nbytes)
        )

    def write(self, fname, append=False):
        """
        Save the file in a format based on the filename extension

        fname: str
            filename (incl. path)
        """
        # non supported types raise an error in self._get_type
        if self._get_type(fname) == "fits":
            self.writeFITS(fname)
        elif self._get_type(fname) == "hdf":
            self.writeHDF(fname, append=append)

    def writeFITS(self, fname, overwrite=False):
        """
        Save to fits file

        Parameters
        ----------
        fname: str
            filename (incl. path) to export to

        overwrite : bool, optional
            Set to overwrite the fits file
        """
        if (self.lamb is not None) & (self.seds is not None) & (self.grid is not None):
            if not isinstance(self.grid, Table):
                raise ValueError("Only astropy.Table are supported")

            hdulist = fits.HDUList()
            hdulist.append(fits.PrimaryHDU(self.lamb))
            hdulist.append(fits.ImageHDU(self.seds, name="seds"))
            if self.cov_diag is not None:
                hdulist.append(fits.ImageHDU(self.cov_diag, name="covdiag"))
            if self.cov_offdiag is not None:
                hdulist.append(fits.ImageHDU(self.cov_offdiag, name="covoffdiag"))
            hdulist.append(fits.BinTableHDU(self.grid, name="grid"))
            hdulist.writeto(fname, overwrite=overwrite)
        else:
            raise ValueError("Full data set not specified (lamb, seds, grid)")

    def writeHDF(self, fname, append=False):
        """
        Save to HDF file

        Parameters
        ----------
        fname : str
            filename (incl. path)

        append : bool, optional (default False)
            if set, it will append data to each Array or Table
        """
        if (self.lamb is not None) & (self.seds is not None) & (self.grid is not None):
            if not isinstance(self.grid, Table):
                raise ValueError("Only astropy.Table are supported")
            with h5py.File(fname, "w") as hd:
                if (not append) or ("seds" not in hd.keys()):
                    hd["seds"] = self.seds[:]
                    hd["lamb"] = self.lamb[:]
                    if self.cov_diag is not None:
                        hd["covdiag"] = self.cov_diag[:]
                    if self.cov_offdiag is not None:
                        hd["covoffdiag"] = self.cov_offdiag[:]
                else:
                    raise Exception("Appending to HDF5 file not supported")

            if getattr(self, "filters", None) is not None:
                if "filters" not in list(self.header.keys()):
                    self.header["filters"] = " ".join(self.filters)
            # append the table of the grid parameters to the file
            self.grid.meta = self.header
            self.grid.write(fname, path="grid", format="hdf5", append=True)
        else:
            raise ValueError("Full data set not specified (lamb, seds, grid)")

    def _from_GridBackend(self, b):
        """
        _from_GridBackend -- convert from generic backend
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

    def _from_DiskBackend(self, b):
        """
        convert from DiskBackend

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


class MemoryBackend(GridBackend):
    """ Instanciate a grid object that has no physical storage

        Helps to create new grids on the fly. Because it derives from
        ModelGrid, this can be exported on disk too.
    """

    def __init__(
        self,
        lamb,
        seds=None,
        grid=None,
        cov_diag=None,
        cov_offdiag=None,
        header={},
        aliases={},
    ):
        """
        Parameters
        ----------
        lamb : ndarray or str or GridBackend subclass
            if ndarray - 1D `float` wavelength of the SEDs (requires seds and grid arguments)
            if str - filename of grid on disk
            if backend - ref to the given grid

        seds : ndarray, optional
            2D `float` array of seds

        grid : astropy.Table, optional
            table of properties associated to each sed

        cov_diag : ndarray, optional
            2D `float` array (# models, # filters) of the diagonal elements
            of the absolute flux covariance matrix

        cov_offdiag : ndarray, optional
            2D `float` array (# models, # elements) of the off diagonal elements
            of the absolute flux covariance matrix

        header : dict, optional
            if provided, update the grid table header

        aliases : dict, optional
            if provided, update the grid table aliases
        """
        super().__init__()

        # read from various formats
        if isinstance(lamb, DiskBackend):
            self._from_DiskBackend(lamb)
        elif isNestedInstance(lamb, GridBackend):
            self._from_GridBackend(lamb)
        elif isinstance(lamb, (str, bytes)):
            self._from_File(lamb)
        else:
            if (seds is None) | (grid is None):
                raise ValueError("seds or grid not passed")
            self.lamb = lamb
            self.seds = seds
            self.grid = grid
            if (cov_diag is not None) & (cov_offdiag is not None):
                self.cov_diag = cov_diag
                self.cov_offdiag = cov_offdiag
            else:
                self.cov_diag = None
                self.cov_offdiag = None

        # update header
        if self._header is None:
            self._header = header
        else:
            for k, v in list(header.items()):
                self.header[k] = v

        # update aliases
        self._aliases.update(aliases)
        self.fname = ":memory:"

    @property
    def filters(self):
        """filters"""
        r = self._header.get("filters", None) or self._header.get("FILTERS", None)
        if r is None:
            return r
        r = r.split()
        return [_decodebytestring(tr) for tr in r]

    def _from_File(self, fname):
        """
        Load the content of a file

        Parameters
        ----------
        fname: str
            filename (incl. path) to read from
        """

        # load_seds - load wavelength and seds
        if self._get_type(fname) == "fits":
            with fits.open(fname) as f:
                extnames = [f[k].header["EXTNAME"].lower() for k in range(1, len(f))]
                if "seds" in extnames:  # new format
                    self.lamb = f[0].data
                    self.seds = f["seds"].data
                    if "covdiag" in extnames:
                        self.cov_diag = f["covdiag"].data
                    else:
                        self.cov_diag = None
                    if "covoffdiag" in extnames:
                        self.cov_offdiag = f["covoffdiag"].data
                    else:
                        self.cov_offdiag = None
                    self._header = f["grid"].header
                    self.grid = Table(f["grid"].data)
                else:  # old format (used for stellar atmosphere grids, remove when those updated)
                    with fits.open(fname) as f:
                        self.seds = f[0].data[:-1]
                        self.lamb = f[0].data[-1]
                        self._header = f[1].header
                    self.grid = Table.read(fname)

        elif self._get_type(fname) == "hdf":
            with h5py.File(fname, "r") as s:
                self.seds = s["seds"][()]
                self.lamb = s["lamb"][()]
                if "covdiag" in s.keys():
                    self.cov_diag = s["covdiag"][()]
                else:
                    self.cov_diag = None
                if "covdiag" in s.keys():
                    self.cov_offdiag = s["covoffdiag"][()]
                else:
                    self.cov_offdiag = None
                self.grid = read_table_hdf5(s["grid"])
                self._header = self.grid.meta

        if "filters" in self._header.keys():
            self._header["filters"] = _decodebytestring(self._header["filters"])

    def copy(self):
        """ implement a copy method """
        g = MemoryBackend(
            copy.deepcopy(self.lamb),
            seds=copy.deepcopy(self.seds),
            grid=copy.deepcopy(self.grid),
            cov_diag=copy.deepcopy(self.cov_diag),
            cov_offdiag=copy.deepcopy(self.cov_offdiag),
        )
        g._filters = copy.deepcopy(self._filters)
        g._header = copy.deepcopy(self._header)
        g._aliases = copy.deepcopy(self._aliases)
        return g


class CacheBackend(GridBackend):
    """
    Load content from a file only when needed
    """

    def __init__(self, fname, *args, **kwargs):
        """
        Parameters
        ----------
        fname : str
            name of file containing the grid
        """
        super().__init__(*args, **kwargs)

        self.fname = fname
        self._type = self._get_type(fname)
        self.clear()

    def clear(self, attrname=None):
        """clear current cache

        Parameters
        ----------
        attrname : str in [lamb, filters, grid, header, lamb, seds]
            if provided clear only one attribute
            else all cache will be erased
        """
        if attrname is None:
            self._seds = None
            self._lamb = None
            self._cov_diag = None
            self._cov_offdiag = None
            self._filters = None
            self._grid = None
            self._header = None
            self._filters = None
        else:
            setattr(self, "_{0}".format(attrname), None)

    def _load_seds(self):
        """
        Load in the SEDs from file if not present
        """
        if self._seds is None:
            if self._type == "fits":
                with fits.open(self.fname) as f:
                    self._seds = f["seds"].data

            elif self._type == "hdf":
                with h5py.File(self.fname, "r") as s:
                    self._seds = s["seds"][()]

    def _load_cov_diag(self):
        """
        Load in the cov_diag from file if not present
        """
        if self._cov_diag is None:
            if self._type == "fits":
                with fits.open(self.fname) as f:
                    self._cov_diag = f["covdiag"].data

            elif self._type == "hdf":
                with h5py.File(self.fname, "r") as s:
                    self._cov_diag = s["covdiag"][()]

    def _load_cov_offdiag(self):
        """
        Load in the cov_offdiag from file if not present
        """
        if self._cov_offdiag is None:
            if self._type == "fits":
                with fits.open(self.fname) as f:
                    self._cov_offdiag = f["covoffdiag"].data

            elif self._type == "hdf":
                with h5py.File(self.fname, "r") as s:
                    self._cov_offdiag = s["covoffdiag"][()]

    def _load_lamb(self):
        """
        Load in the wavelengths from file if not present
        """
        if self._lamb is None:
            if self._type == "fits":
                with fits.open(self.fname) as f:
                    self._lamb = f[0].data

            elif self._type == "hdf":
                with h5py.File(self.fname, mode="r") as s:
                    self._lamb = s["lamb"][()]

    def _load_grid(self):
        """
        Load in the grid from file if not present
        """
        if self._grid is None:
            if self._type == "fits":
                self._grid = Table.read(self.fname)

            elif self._type == "hdf":
                self._grid = Table.read(self.fname, path="grid", format="hdf5")

    def _load_header(self):
        """
        Load in the header of the grid if not present
        """
        if self._header is None:
            if self._type == "fits":
                with fits.open(self.fname) as f:
                    self._header = f["grid"].header
            elif self._type == "hdf":
                with h5py.File(self.fname, mode="r") as s:
                    self._header = _gethdfdatasetmeta(s["grid"])

    def _load_filters(self):
        """
        Load in the filters from file if not present
        """
        if self._filters is None:
            self._load_header()
            try1 = self._header.get("FILTERS", None)
            try2 = self._header.get("filters", None)
            self._filters = try1 or try2
            if self._filters is not None:
                self._filters = self._filters.split()

    @property
    def seds(self):
        self._load_seds()
        return self._seds

    @seds.setter
    def seds(self, value):
        self._seds = value

    @property
    def lamb(self):
        self._load_lamb()
        return self._lamb

    @lamb.setter
    def lamb(self, value):
        self._lamb = value

    @property
    def cov_diag(self):
        self._load_cov_diag()
        return self._cov_diag

    @cov_diag.setter
    def cov_diag(self, value):
        self._cov_diag = value

    @property
    def cov_offdiag(self):
        self._load_cov_offdiag()
        return self._cov_offdiag

    @cov_offdiag.setter
    def cov_offdiag(self, value):
        self._cov_offdiag = value

    @property
    def grid(self):
        self._load_grid()
        return self._grid

    @grid.setter
    def grid(self, value):
        self._grid = value

    @property
    def header(self):
        self._load_header()
        return self._header

    @header.setter
    def header(self, value):
        self._header = value

    @property
    def filters(self):
        self._load_filters()
        return self._filters

    @filters.setter
    def filters(self, value):
        self._filters = value

    def keys(self):
        """ return column names when possible, avoid loading when possible """
        if hasattr(self._grid, "keys"):
            return list(self._grid.keys())
        else:
            super().keys()

    def copy(self):
        """ implement a copy method """
        g = CacheBackend(self.fname)
        g._aliases = copy.deepcopy(self._aliases)
        if self._grid is not None:
            g._grid = copy.deepcopy(self._grid)
        if self._seds is not None:
            g._seds = copy.deepcopy(self._seds)
        if self._lamb is not None:
            g._lamb = copy.deepcopy(self._lamb)
        if self._cov_diag is not None:
            g._cov_diag = copy.deepcopy(self._cov_diag)
        if self._cov_offdiag is not None:
            g._cov_offdiag = copy.deepcopy(self._cov_offdiag)
        if self._header is not None:
            g._header = copy.deepcopy(self._header)
        if self._filters is not None:
            g._filters = copy.deepcopy(self._filters)

        return g


class DiskBackend(GridBackend):
    """
    Reads the data from disk when it is accessed.  This supports reading
    only a portion (e.g., slices/subsets) of the requsted data.
    This allows spectral and SED grids larger than can fit into memory.

    Only hdf files supported.
    """

    def __init__(self, fname, *args, **kwargs):
        super().__init__(*args, **kwargs)
        ftype = self._get_type(fname)
        if ftype != "hdf":
            raise ValueError("Expecting HDF file got {0}".format(ftype))

        self.fname = fname
        self.store = h5py.File(self.fname, mode="r")
        self.seds = self.store["seds"]
        self.lamb = self.store["lamb"]
        self.grid = self.store["grid"]
        self.cov_diag = None
        if "covdiag" in self.store.keys():
            self.cov_diag = self.store["covdiag"]
        if "covoffdiag" in self.store.keys():
            self.cov_offdiag = self.store["covoffdiag"]
        self._filters = None
        self._header = None
        self._aliases = {}

    @property
    def header(self):
        if self._header is None:
            self._header = _gethdfdatasetmeta(self.grid)
        return self._header

    @property
    def filters(self):
        """load in cache if needed """
        if self._filters is None:
            if self._header is None:
                self._header = _gethdfdatasetmeta(self.grid)
            self._filters = self._header.get("FILTERS", None) or self._header.get(
                "filters", None
            )
            if self._filters is not None:
                self._filters = self._filters.split()
        return self._filters

    def keys(self):
        """ return column names when possible """
        return list(self.grid.dtype.fields.keys())

    def copy(self):
        g = DiskBackend(self.fname)
        g._aliases = copy.deepcopy(self._aliases)
        return g
