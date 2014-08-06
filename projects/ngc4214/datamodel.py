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
from special import make_top_hat_filter

#---------------------------------------------------------
# User inputs                                   [sec:conf]
#---------------------------------------------------------
# Parameters that are required to make models
# and to fit the data
#---------------------------------------------------------
project = 'mf_ngc4214_Aug2014'

obsfile = 'data/N4214_3band_detects.fits'
astfile = 'data/N4214_gst_fake.fits'

filters = ['HST_WFC3_F225W', 'HST_WFC3_F336W', 'HST_WFC3_F438W',
           'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']


# absflux calibration covariance matrix for NGC4214 specific filters
calibration_covariance = np.array(
    [[1.80e-4,  1.37e-4, 6.02e-5, 2.44e-5, 1.23e-6, -4.21e-6],
    [1.37e-4,  1.09e-4, 5.72e-5, 3.23e-5, 1.65e-5, 1.32e-5],
    [6.02e-5,  5.72e-5, 5.07e-5, 4.66e-5, 4.40e-5, 4.34e-5],
    [2.44e-5,  3.23e-5, 4.66e-5, 5.42e-5, 5.87e-5, 5.99e-5],
    [1.23e-6,  1.65e-5, 4.40e-5, 5.87e-5, 6.98e-5, 7.33e-5],
    [-4.21e-6, 1.32e-5, 4.34e-5, 5.99e-5, 7.33e-5, 7.81e-5] ])

noisefile = 'ngc4214/mf_phat_july14_noisemodel.hd5'

distanceModulus = 27.41 * unit['mag']

#Spectral grid definition
#log10(Age) -- [min,max,step] to generate the isochrones
logt = [6.0, 10.13, 0.15]

#note: Mass is not sampled, use the isochrone def instead.

#Metallicity
z = [0.03, 0.019, 0.008, 0.004]
z = 0.004  # Strong prior SMC metalliticy

#Stellar library definition
osl = stellib.Tlusty() + stellib.Kurucz()

#Extinction's law definition
extLaw = extinction.RvFbumpLaw()
avs = [0., 10.055,0.15]
rvs = [2.0, 6.1, 1.0]
fbumps = [0., 1.01, 0.25]
## ..note::
##      in the above the upper limit makes sure the last point of interest is
##      included

# extracting non-default spectral properties of the models
# -------------------------------------------------------

qion_filter = make_top_hat_filter(912., 2066., 1, 'F_QION')
qion_filter.name = 'F_QION'  # getting rid of instrument etc
# note: remember to multiply by bandwidth to get the actual energy

add_spectral_properties_kwargs = dict(filternames=filters, filters=[qion_filter])


class NGC4214_FluxCatalog(Observations):
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
        # in RATE = flux units
        self.setBadValue(6e-11)

        #hard code mapping directly with the interface
        for k in filters:
            self.data.set_alias(k.upper(), k.split('_')[-1].upper() + '_VEGA')

    def getFlux(self, num, units=False):
        """returns the absolute flux of an observation from the catalog

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
        # catalog uses _VEGA which are vega magnitudes.
        d = self.data[num]
        flux = 10 ** (-0.4 * np.array([ d[self.data.resolve_alias(ok)] for ok in self.filters ])) * self.vega_flux

        if units is True:
            return flux * unit['erg/s/cm**2']
        else:
            return flux

    def getFluxerr(self, num, units=False):
        """returns the error on the absolute flux of an observation from the
        number of counts (not used in the analysis)

        Parameters
        ----------
        num: int
            index of the star in the catalog to get measurement from

        units: bool
            if set returns the fluxes with a unit capsule

        Returns
        -------
        fluxerr: ndarray[dtype=float, ndim=1]
            Measured integrated flux uncertainties in erg/s/cm^2
        """

        fluxerr = 10 ** (-0.4 * Observations.getFluxerr(self, num)) * self.vega_flux
        if units is True:
            return fluxerr * unit['erg/s/cm**2']
        else:
            return fluxerr

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
    obs = NGC4214_FluxCatalog(obsfile, distanceModulus=distanceModulus, filters=filters)
    return obs
