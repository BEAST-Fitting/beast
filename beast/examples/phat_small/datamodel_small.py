""" Data Model interface v2.0
BEAST datamodel for the example based on PHAT data
KDG - 21 Dec 2016
"""

import numpy as np

# BEAST imports
from beast.physicsmodel.stars import stellib
from beast.physicsmodel.dust import extinction
from beast.observationmodel.observations import Observations
from beast.observationmodel.vega import Vega
from beast.observationmodel.noisemodel import absflux_covmat
from beast.external.ezunits import unit
#from extra_filters import make_integration_filter, make_top_hat_filter

#-----------------------------------------------------------------
# User inputs                                   [sec:conf]
#-----------------------------------------------------------------
# Parameters that are required to make models
# and to fit the data
#-----------------------------------------------------------------
# AC == authomatically created
# indicates where user's input change is NOT necessary/recommeded
#-----------------------------------------------------------------

# project : string
# the name of the output results directory
project = 'beast_example_phat'

# filters : list of strings
# full filter names in BEAST filter database
filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
           'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
# basefilters : list of strings
# short names for filters
basefilters = ['F275W','F336W','F475W',
               'F814W','F110W','F160W']
# obs_colnames : list of strings
# names of columns for filters in the observed catalog
# need to match column names in the observed catalog,
# input data MUST be in fluxes, NOT in magnitudes 
# fluxes MUST be in normalized Vega units
obs_colnames = [ f + '_rate' for f in basefilters ]
# ast_colnames : list of strings 
# names of columns for filters in the AST catalog (AC)
ast_colnames = np.array(basefilters)

# bright_limits_mag, sens_limits_mag : lists of floats
# sensitivity limits (used for AST input generation)
# units are Vega magnitudes
bright_limits_mag = [14., 14.5, 16., 15., 16., 14., 14.5, 14., 14.]
sens_limits_mag = [26., 26., 27., 29., 27.5, 28., 28.5, 27., 26.]

# obsfile : string 
# pathname of the observed catalog
obsfile = 'data/b15_4band_det_27_A.fits'

# astfile : string
# pathname of the AST files (single camera ASTs)
astfile = 'data/fake_stars_b15_27_all.hd5'

# noisefile : string
# create a name for the noise model
noisefile = project + '/' + project + '_noisemodel.hd5'

# absflux calibration covariance matrix for HST specific filters (AC)
absflux_a_matrix = absflux_covmat.hst_frac_matrix(filters)

# distance modulus to the galaxy
distanceModulus = 24.47 * unit['mag']

### Stellar grid definition

# log10(Age) -- [min,max,step] to generate the isochrones in years
# recommended [6.0, 10.13, 1.0]
logt = [6.0, 10.13, 1.0]

# note: Mass is not sampled, use the isochrone def instead.

# Metallicity : list of floats
# acceptable z > 0.0
# reasonable [min, max] = [0.003, 1.3] 
# can they be set as [min, max, step]?
z = [0.03, 0.019, 0.008, 0.004]

# Isochrone CMD version (2.3 for Girardi et al. (2010) or 2.7 for PARSECv1.2S)
trackVersion = 2.7

# Stellar Atmospheres library definition
osl = stellib.Tlusty() + stellib.Kurucz()

################

### Dust extinction grid definition
extLaw = extinction.Gordon16_RvFALaw()

# A(V): dust column in magnitudes
# acceptable avs > 0.0
# reasonable [min, max, step] = [0.0, 10.055, 1.0]
avs = [0.0, 10.055, 1.0]

# R(V): dust average grain size
# rasonable [min, max, step] = [2.0,6.0,1.0]
rvs = [2.0,6.0,1.0]

# fbump (should be f_A): mixture factor between
#     "MW" and "SMCBar" extinction curves
# reasonable [min, max, step] = [0.0,1.0, 0.25]
fbumps = [0.0,1.0, 0.25]

################

# add in the standard filters to enable output of stats and pdf1d values
# for the observed fitlers (AC)
add_spectral_properties_kwargs = dict(filternames=filters)

################
# The following code does not require user's attention (AC)
################

class PHATFluxCatalog(Observations):
    """PHAT 6 filter photometry
    This class implements a direct access to the PHAT measured fluxes.

    ..note::
        it does not implement uncertainties as in this model, the noise is
        given through artificial star tests
    """
    def __init__(self, inputFile, distanceModulus=distanceModulus, filters=filters):
        """ Construct the interface """
        desc = 'PHAT star: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        #some bad values smaller than expected
        # in physical flux units
        self.setBadValue(6e-40)

        #hard code mapping directly with the interface to HTTP
        for k in filters:
            self.data.set_alias(k, k.split('_')[-1].lower() + '_rate')

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
        flux = np.array([ d[self.data.resolve_alias(ok)] 
                          for ok in self.filters ]) * self.vega_flux
        
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
    obs = PHATFluxCatalog(obsfile, distanceModulus=distanceModulus, filters=filters)
    return obs
