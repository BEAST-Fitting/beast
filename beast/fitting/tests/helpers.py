import numpy as np

from beast.observationmodel.observations import Observations
from beast.observationmodel.vega import Vega


class GenFluxCatalog(Observations):
    """
    This class implements a direct access to the Generic HST measured fluxes.
    """

    def __init__(self, inputFile, filters, obs_colnames, vega_fname=None):
        """
        Parameters
        ----------
        inputFile : str
            observation file

        filters : list
            interal filter names of the data

        obs_colnames : list
            filter names in the observed catalog

        vega_fname : str, optional
            name of the file with the vega model spectrum
        """

        desc = "GENERIC star: %s" % inputFile
        Observations.__init__(self, inputFile, desc=desc)
        self.setFilters(filters, vega_fname=vega_fname)
        # some bad values smaller than expected
        # in physical flux units
        self.setBadValue(6e-40)

        # rate column needed as this is the *flux* column
        for ik, k in enumerate(filters):
            self.filter_aliases[k] = obs_colnames[ik]

    def getFlux(self, num, units=False):
        """returns the absolute flux of an observation

        Parameters
        ----------
        num : int
            index of the star in the catalog to get measurement from

        units : bool
            if set returns the fluxes with a unit capsule

        Returns
        -------
        flux : ndarray[dtype=float, ndim=1]
            Measured integrated flux values throughout the filters
            in erg/s/cm^2/A
        """

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

    def setFilters(self, filters, vega_fname=None):
        """
        Set the filters and update the vega reference for the conversions

        Parameters
        ----------
        filters : list
            list of filters using the internally normalized namings

        vega_fname : str, optional
            name of the file with the vega model spectrum
        """
        self.filters = filters

        # Data "rates" are normalized to Vega already, fits are not using vega

        # for optimization purpose: pre-compute
        #  getting vega mags, require to open and read the content of one file.
        #  since getObs, calls getFlux, for each star you need to do this
        #  expensive operation
        with Vega(source=vega_fname) as v:
            _, vega_flux, _ = v.getFlux(filters)

        self.vega_flux = vega_flux


def get_obscat(obsfile, filters, obs_colnames, vega_fname=None, *args, **kwargs):
    """
    Generates a data catalog object

    Parameters
    ----------
    obsfile : str
        observation file

    filters : list
        interal filter names of the data

    obs_colnames : list
        filter names in the observed catalog

    vega_fname : str, optional
        name of the file with the vega model spectrum

    Returns
    -------
    obs : GenFluxCatalog
        observation catalog
    """
    obs = GenFluxCatalog(obsfile, filters, obs_colnames, vega_fname=vega_fname)
    return obs
