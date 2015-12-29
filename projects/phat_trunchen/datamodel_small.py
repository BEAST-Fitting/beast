""" Data Model interface v2.0

Includes artificial star tests (ASTs) and related function to generate the
noise model

The noise model is then applied to the models while only data measurements are
used in constrasts with using measurement uncertainties.

with limited quantity units handling to help both reading and robustness.

..note::
    By default units on returned values from function calls are turned off to
    avoid breaking possible other scripts.
"""
import numpy as np

# BEAST imports
from beast.core import stellib
from beast.core import extinction
from beast.core.observations import Observations
from beast.core.vega import Vega
from beast.external.ezunits import unit
from beast.core.noisemodel import absflux_covmat
from extra_filters import make_integration_filter, make_top_hat_filter

#---------------------------------------------------------
# User inputs                                   [sec:conf]
#---------------------------------------------------------
# Parameters that are required to make models
# and to fit the data
#---------------------------------------------------------
project = 'b15_dec15_small'

filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
           'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
basefilters = ['F275W','F336W','F475W','F814W','F110W','F160W']
obs_colnames = [ f + '_rate' for f in basefilters ]
ast_colnames = np.array(basefilters)

# observations
obsfile = 'data_small/b15_4band_det_27_A.fits'

# AST file
astfile = 'ASTs_6band/12057_M31-B09-F02.gst.fake.fits'

# name for noise model
noisefile = project + '/' + project + '_noisemodel.hd5'

# absflux calibration covariance matrix for PHAT specific filters
absflux_cov = True
if not absflux_cov:
    absflux_a_matrix = absflux_covmat.hst_frac_matrix(filters)
    print(np.sqrt(absflux_a_matrix))

# distance to M31
distanceModulus = 24.47 * unit['mag']

### Stellar grid definition

# log10(Age) -- [min,max,step] to generate the isochrones
logt = [6.0, 10.13, 2.0]

#note: Mass is not sampled, use the isochrone def instead.

#Metallicity
#z = [0.03, 0.019, 0.008, 0.004]
z = 0.019

# Isochrone CMD version (2.3 for Girardi et al. (2010) or 2.7 for PARSECv1.2S)
trackVersion = 2.7

# Stellar Atmospheres library definition
osl = stellib.Tlusty() + stellib.Kurucz()

################

### Dust extinction grid definition
extLaw = extinction.RvFbumpLaw()

# A(V): dust column
avs = [0.0, 5.055, 10.0]

# R(V): dust average grain size
rvs = [2.0,6.0,5.0]

# fbump (should be f_A): mixture factor between "MW" and "SMCBar" extinction curves
fbumps = None
#fbumps = [0.0,1.0, 1.0]

################

### extra filters for specific projects

qion_filter = make_integration_filter(90., 912., 1, 'F_QION')
qion_filter.name = 'F_QION'  # getting rid of instrument etc
maria_filter = make_top_hat_filter(912., 2066., 1, 'F_UV_6_13e')
maria_filter.name = 'F_UV_6_13e'

additional_filters = ['GALEX_FUV', 'GALEX_NUV']
# note: remember to multiply by bandwidth to get the actual energy

add_spectral_properties_kwargs = dict(filternames=filters + additional_filters,
                                      filters=[qion_filter, maria_filter])

################

## ..note::
##      in the above grid definitions the upper limit makes sure the last point of interest is
##      included


class PHATFluxCatalog(Observations):
    """PHAT 6 filter photometry
    This class implements a direct access to the PHAT measured fluxes.

    ..note::
        it does not implement uncertainties as in this model, the noise is
        given through artificial star tests
    """
    def __init__(self, inputFile, distanceModulus=distanceModulus,
                 filters=filters):
        """ Construct the interface """
        desc = 'PHAT star: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        #some bad values smaller than expected
        # in RATE = flux units
        self.setBadValue(6e-30)

        #hard code mapping directly with the interface to PHAT
        for k in filters:
            self.data.set_alias(k, k.split('_')[-1].lower() + '_rate')

    def getFlux(self, num, units=False):
        """returns the absolute flux of an observation from the number of
        counts

        Parameters
        ----------
        num: int
            index of the star in the catalog to get measurement from

        units: bool
            if set returns the fluxes with a unit capsule

        Returns
        -------
        flux: ndarray[dtype=float, ndim=1]
            Measured integrated flux values throughout the filters in
            erg/s/cm^2
        """

        # case for using '_RATE' result instead of '_VEGA'
        d = self.data[num]
        flux = np.array([ d[self.data.resolve_alias(ok)] for
                          ok in self.filters ]) * self.vega_flux
        #flux = Observations.getFlux(self, num) * self.vega_flux
        
        # case for using '_VEGA' which are vega magnitudes.
        #d = self.data[num]
        #flux = 10 ** (-0.4 * np.array([ d[self.data.resolve_alias(ok)]
        #                 for ok in self.filters ])) * self.vega_flux

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
        #   getting vega mags, require to open and read the content of one
        #   file.  since getObs, calls getFlux, for each star you need to
        #   do this expensive op.
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
    obs = PHATFluxCatalog(obsfile, distanceModulus=distanceModulus,
                          filters=filters)
    return obs
