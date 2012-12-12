"""
Handling PHAT SED Cluster analysis

   [Globals]	phat filters, distance, Rv
   [Data]	PHATData(Observation) class, vegamag to fluxes
   [Extinction] Extinction to use, incl. Av sampling
   [Models]	Which models to use

"""
import numpy
from observations import Observations
from anased import extinction
from vega import Vega, from_Vegamag_to_Flux


""" PHAT Globals """

PHAT_filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
		'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']

""" Assertions """
distanceModulus = 24.3

Av_law          = extinction.Cardelli()
##define Av with a number of points
Av = numpy.linspace(0, 1, 3)
Rv = numpy.array([3.1], dtype=float)
fb = numpy.array([0.5], dtype=float)     # only one value still needs to be array/list/tuple


""" Data handling """

#Data are in Vega magnitudes
#  Need to use Vega
with Vega() as v:
	vega_f, vega_mag, lamb = v.getMag(PHAT_filters)


class PhatData(Observations):
    """ PHAT catalog for clusters in M31 """
    def __init__(self, inputFile, distanceModulus=distanceModulus):
        desc = 'PHAT cluster catalog yr1: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( PHAT_filters )
        self.setBadValue(99.0)

    @from_Vegamag_to_Flux(lamb, vega_mag)
    def getObs(self, num):
        """ Using the decorator @from_Vegamag_to_Flux
            Hence, results are in flux (not in flux/flux_vega)

            Returns the fluxes, errors and mask of an observation.
        """
        return Observations.getObs(self, num)

    def getObsinMag(self, num):
        """ Returns the original catalog magnitudes """
        return Observations.getObs(self, num)

""" DATA """
#obs = PhatData(obsfile, distanceModulus)
#obs.setFilter(PHAT_filters)
