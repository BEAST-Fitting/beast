""" Data Model interface v2.0
BEAST datamodel for the SMC SMIDGE data.
Karl G. - 3 Dec 2015
  - adapted from the BEAST datamodel for PHAT data
"""

import numpy as np

# BEAST imports
from beast.core import stellib
from beast.core import extinction
from beast.core.observations import Observations
from beast.core.vega import Vega
from beast.core.noisemodel import absflux_covmat
from beast.external.ezunits import unit

#---------------------------------------------------------
# User inputs                                   [sec:conf]
#---------------------------------------------------------
# Parameters that are required to make models
# and to fit the data
#---------------------------------------------------------
project = 'smidge_dec15_small'

filters = ['HST_WFC3_F225W', 'HST_WFC3_F275W', 'HST_WFC3_F336W', 
           'HST_ACS_WFC_F475W','HST_ACS_WFC_F550M', 'HST_ACS_WFC_F814W',
           'HST_WFC3_F110W', 'HST_WFC3_F160W']

# observations
obsfile = 'data/13659_SMC-F14.gst.fits'

# AST files (single camera ASTs)
astfile = ''

# name for noise model
noisefile = project + '/' + project + '_noisemodel.hd5'

# absflux calibration covariance matrix for HST specific filters

absflux_a_matrix = absflux_covmat.hst_frac_matrix(filters)

# distance to the SMC
distanceModulus = 18.9 * unit['mag']

### Stellar grid definition

# log10(Age) -- [min,max,step] to generate the isochrones
logt = [6.0, 10.13, 0.25]

#note: Mass is not sampled, use the isochrone def instead.

#Metallicity
#z = [0.03, 0.019, 0.008, 0.004]
z = 0.004

# Isochrone CMD version (2.3 for Girardi et al. (2010) or 2.7 for PARSECv1.2S)
trackVersion = 2.7

# Stellar Atmospheres library definition
osl = stellib.Tlusty() + stellib.Kurucz()

################

### Dust extinction grid definition
extLaw = extinction.RvFbumpLaw()

# A(V): dust column
avs = [0.0, 10.055, 0.5]

# R(V): dust average grain size
rvs = [2.0,6.0,0.5]

# fbump (should be f_A): mixture factor between "MW" and "SMCBar" extinction curves
fbumps = [0.0,1.0, 0.10]

################

## ..note::
##      in the above grid definitions the upper limit makes sure the last point of interest is
##      included


class SMIDGEFluxCatalog(Observations):
    """SMIDGE 8 filter photometry
    This class implements a direct access to the SMIDGE measured fluxes.

    ..note::
        it does not implement uncertainties as in this model, the noise is
        given through artificial star tests
    """
    def __init__(self, inputFile, distanceModulus=distanceModulus, filters=filters):
        """ Construct the interface """
        desc = 'HTTP star: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        #some bad values smaller than expected
        # in physical flux units
        self.setBadValue(6e-40)

        #hard code mapping directly with the interface to HTTP
        for k in filters:
            self.data.set_alias(k, k.split('_')[-1].upper() + '_RATE')

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
            Measured integrated flux values throughout the filters in erg/s/cm^2
        """

        # case for using '_flux' result
        d = self.data[num]
        flux = np.array([ d[self.data.resolve_alias(ok)] for ok in self.filters ]) * self.vega_flux
        
        if units is True:
            return flux * unit['erg/s/cm**2']
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
        #   since getObs, calls getFlux, for each star you need to do this expensive
        #   op.
        with Vega() as v:
            _, vega_flux, _ = v.getFlux(filters)

        self.vega_flux = vega_flux


def get_obscat(obsfile=obsfile, distanceModulus=distanceModulus,
               filters=filters, *args, **kwargs):
    """ Function that generates a data catalog object with the correct
    arguments

    Parameters
    ----------
    obsfile: str, optional (default datamodel.obsfile)
        observation file

    distanceModulus: float, optional (default datamodel.distanceModulus)
        distance modulus to correct the data from (in magitude)

    filters: sequence(str), optional, datamodel.filters
        seaquence of filters of the data

    returns
    -------
    obs: PHATFluxCatalog
        observation catalog
    """
    obs = SMIDGEFluxCatalog(obsfile, distanceModulus=distanceModulus, filters=filters)
    return obs
