from __future__ import (absolute_import, division, print_function)

import numpy as np

# BEAST imports
from beast.observationmodel.observations import Observations
from beast.observationmodel.vega import Vega



class GenFluxCatalog(Observations):
    """Generic n band filter photometry
    This class implements a direct access to the Generic HST measured fluxes.

    ..note::
        it does not implement uncertainties as in this model, the noise is
        given through artificial star tests
    """
    def __init__(self, inputFile, filters=filters):
        """ Construct the interface """
        desc = 'GENERIC star: %s' % inputFile
        Observations.__init__(self, inputFile, desc=desc)
        self.setFilters( filters )
        #some bad values smaller than expected
        # in physical flux units
        self.setBadValue(6e-40)

        # rate column needed as this is the *flux* column
        for ik,k in enumerate(filters):
            self.data.set_alias(k, obs_colnames[ik])

    def getFlux(self, num, units=False):
        """returns the absolute flux of an observation 

        Parameters
        ----------
        num: int
            index of the star in the catalog to get measurement from

        units: bool
            if set returns the fluxes with a unit capsule

        Returns
        -------
        flux: ndarray[dtype=float, ndim=1]
            Measured integrated flux values throughout the filters 
            in erg/s/cm^2/A
        """

        # case for using '_flux' result
        d = self.data[num]
        
        flux = np.array([ d[self.data.resolve_alias(ok)] 
                          for ok in self.filters ]) * self.vega_flux
        
        if units is True:
            return flux * units.erg / (units.s*units.cm*units.cm*units.angstrom)
        else:
            return flux

    def setFilters(self, filters):
        """ set the filters and update the vega reference for the conversions

        Parameters
        ----------
        filters: sequence
            list of filters using the internally normalized namings
        """
        self.filters = filters

        #Data "rates" are normalized to Vega already, fits are not using vega

        # for optimization purpose: pre-compute
        #   getting vega mags, require to open and read the content of one file.
        #   since getObs, calls getFlux, for each star you need to do this
        #   expensive operation
        with Vega() as v:
            _, vega_flux, _ = v.getFlux(filters)

        self.vega_flux = vega_flux


def get_obscat(obsfile=obsfile, filters=filters, *args, **kwargs):
    """ Function that generates a data catalog object with the correct
    arguments

    Parameters
    ----------
    obsfile: str, optional (default datamodel.obsfile)
        observation file

    filters: sequence(str), optional, datamodel.filters
        seaquence of filters of the data

    returns
    -------
    obs: GenFluxCatalog
        observation catalog
    """
    obs = GenFluxCatalog(obsfile, filters=filters)
    return obs
