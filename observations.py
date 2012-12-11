""" Defines a generic interface to observation catalog
    This enables to handle non detections, (upper limits one day?), flux and
    magnitude conversions to avoid painful preparation of the dataset

UNDER DEV

TODO:
    convert magnitudes to fluxes
        * need to store the magnitudes type {vega, ab, st}
        * generate a converter especially for vega...

"""
import numpy


"""
STMAGS --- Convert an ST magnitude to erg/s/cm2/AA (Flambda)
      mag = -2.5*log10(F) - 21.10
      M0 = 21.10
      F0 = 3.6307805477010028e-09 erg/s/cm2/AA
"""


def STmag_to_flux( v ):
    v0 = 21.1
    return numpy.power(10., -0.4 * (v - v0) )


def STmag_from_flux( v ):
    v0 = 21.1
    return -2.5 * numpy.log10( v ) - v0


class Observations(object):

    def __init__(self, inputFile, distanceModulus=0., desc=None):
        """ Generate a data interface object """
        self.inputFile = inputFile
        self.filters   = None
        self.desc      = desc
        self.setDistanceModulus(distanceModulus)
        self.readData()
        self.badvalue  = None

    @property
    def nObs(self):
        return self.data.nrows

    def __len__(self):
        return self.nObs

    def __call__(self):
        """ Calling the object will show info """
        print "Data read from %s " % self.inputFile
        if self.desc is not None:
            print "Description: %s" % self.desc
            print "Number of records: %d" % self.nObs
            print ""
            print "Dataset contains:"

        for k in self.data.keys():
            print "\t %s" % k

        if self.filters is None:
            print "No filters set yet!"
        else:
            print "Using filters:", self.filters

    def __getitem__(self, *args, **kwargs):
        """ get item will generate a subsample """
        return self.data.__getitem__(*args, **kwargs)

    def keys(self):
        """ Returns dataset content names """
        return self.data.keys()

    def setDescription(self, txt):
        self.desc = txt

    def setDistanceModulus(self, val):
        """ Set the distance modulus to consider the dataset """
        self.distanceModulus = val
        self.distance = 10 ** ( (val - 25.) / 5. )

    def setDistance(self, val):
        """ Set observed object distance to X Megaparsecs
            this will update also the distance Modulus
        """
        self.distance = val
        self.distanceModulus = 5. * numpy.log10( val * 1e5 )

    def setBadValue(self, val):
        self.badvalue = val

    def getFilters(self):
        return self.filters

    def setFilters(self, filters):
        self.filters = filters

    def getMags(self, num, filters):
        return numpy.array([ self.data[tt][num] - self.distanceModulus for tt in filters])

    def getErrors(self, num, filters):
        return numpy.array([ self.data[tt + 'err'][num] for tt in filters])

    def getObs(self, num=0):
        """ returns the dictionnary used during the analysis """
        assert ( not self.filters is None), "No filter set."
        mags = self.getMags(num, self.filters)
        errs = self.getErrors(num, self.filters)
        if not self.badvalue is None:
            mask = (mags >= self.badvalue)
        else:
            mask = numpy.zeros(len(mags), dtype=bool)

        return mags, errs, mask

    def readData(self):
        """ read the dataset from the original source file """
        import mytables
        self.data = mytables.load(self.inputFile)

    def iterobs(self):
        for k in range(self.nObs):
            yield self.getObs(k)

    def enumobs(self):
        for k in range(self.nObs):
            yield k, self.getObs(k)


class FakeObs(Observations):

    def getObs(self, num=0, err=0.05):
        assert ( self.filters is not None), "No filter set."
        mags = self.getMags(num, self.filters)
        errs = numpy.ones(len(mags), dtype=float) * err
        if self.badvalue is not None:
            mask = (mags >= self.badvalue)
        else:
            mask = numpy.zeros(len(mags), dtype=bool)

        return mags, errs, mask

    def readData(self):
        """ read the dataset from the original source file """
        import mytables
        self.data = mytables.load(self.inputFile)
