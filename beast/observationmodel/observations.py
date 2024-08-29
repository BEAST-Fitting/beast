import numpy as np

from astropy.table import Table

from beast.observationmodel.vega import Vega

__all__ = ["Observations"]


class Observations(object):
    """
    A generic class that interfaces observation catalog in a standardized way

    Attributes
    ----------
    inputFile : str
        catalog source file
    filters : list
        list of filter names (internal standards)
    filter_aliases : dict
        alias of filter names between internal and external names
    desc : str
        description of the observations
    badvalue : float
        value that tags a bad measurement that should not be used in the
        fitting
    nObs : int
        number of observations in the catalog
    """

    def __init__(
        self, inputFile, filters, obs_colnames=None, vega_fname=None, desc=None
    ):
        """
        Parameters
        ----------
        inputFile : str
            observation file
        filters : list
            interal filter names of the data
        obs_colnames : list, optional
            filter names in the observed catalog
        vega_fname : str, optional
            name of the file with the vega model spectrum
        desc : str, optional
            description of the observations
        """
        if desc is None:
            self.desc = "GENERIC: %s" % inputFile
        else:
            self.desc = desc
        self.inputFile = inputFile
        self.setFilters(filters)
        self.filter_aliases = {}
        for ik, k in enumerate(filters):
            self.filter_aliases[k] = obs_colnames[ik]
        self.readData()
        self.setVegaFluxes(filters, vega_fname=vega_fname)
        # some bad values smaller than expected
        # in physical flux units
        self.setBadValue(6e-40)

    @property
    def nObs(self):
        return len(self.data)

    def __len__(self):
        return self.nObs

    def __call__(self):
        """Calling the object will show info"""
        self.info()

    def info(self):
        """Prints some information about the catalog"""
        txt = "Data read from {s.inputFile:s}\n"
        if self.desc is not None:
            txt += "Description: {s.desc:s}\n"
        txt += "Number of records: {s.nObs:d}\n\n"
        txt += "Dataset contains:"

        print("Data read from %s " % self.inputFile)
        if self.desc is not None:
            print("Description: %s" % self.desc)
            print("Number of records: %d" % self.nObs)
            print("")
            print("Dataset contains:")

        for k in list(self.data.keys()):
            txt += "\t {0:s}\n".format(k)

        if self.filters is None:
            txt += "\n No filters given yet!"
        else:
            txt += "\n Using Filters: {s.filters}\n"

        print(txt.format(s=self))

    def __getitem__(self, *args, **kwargs):
        """get item will generate a subsample"""
        return self.data.__getitem__(*args, **kwargs)

    def keys(self):
        """Returns dataset content names"""
        return self.data.keys()

    def setDescription(self, txt):
        self.desc = txt

    def setBadValue(self, val):
        self.badvalue = val

    def getFilters(self):
        return self.filters

    def setFilters(self, filters):
        self.filters = filters

    def setVegaFluxes(self, filters, vega_fname=None):
        """
        Set vega reference fluxes for conversions

        Parameters
        ----------
        filters : list
            list of filters using the internally normalized namings
        vega_fname : str, optional
            name of the file with the vega model spectrum
        """
        # for optimization purpose: pre-compute
        with Vega(source=vega_fname) as v:
            _, vega_flux, _ = v.getFlux(filters)
        self.vega_flux = vega_flux

    def getFlux(self, num, units=False):
        """
        Flux of an observation computed from normalized vega fluxes

        Parameters
        ----------
        num : int
            index of the star in the catalog to get measurement from
        units : bool
            if set returns the fluxes with units

        Returns
        -------
        flux : ndarray[dtype=float, ndim=1]
            Measured integrated flux values throughout the filters
            in erg/s/cm^2/A
        """
        if self.vega_flux is None:
            raise ValueError("vega_flux not set, can't return fluxes")

        # case for using '_flux' result
        d = self.data[num]

        flux = (
            np.array([d[self.filter_aliases[ok]] for ok in self.filters])
            * self.vega_flux
        )

        if units is True:
            return flux * units.erg / (units.s * units.cm * units.cm * units.angstrom)
        else:
            return flux

    def getFluxerr(self, num):
        """returns the error on the flux of an observation from the number of
        counts (not used in the analysis)"""

        fluxerr = np.zeros(len(self.filters), dtype=float)

        for ek, ok in enumerate(self.filters):
            fluxerr[ek] = self.data[ok + "_err"][num]

        return fluxerr

    def getObs(self, num=0):
        """returns the flux"""
        if self.filters is None:
            raise AttributeError("No filter set provided.")

        flux = self.getFlux(num)

        return flux

    def readData(self):
        """read the dataset from the original source file"""

        if isinstance(self.inputFile, str):
            self.data = Table.read(self.inputFile)
        else:
            self.data = self.inputFile

    def iterobs(self):
        """yield getObs"""
        for k in range(self.nObs):
            yield self.getObs(k)

    def enumobs(self):
        for k in range(self.nObs):
            yield k, self.getObs(k)
