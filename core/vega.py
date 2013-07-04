""" Handle vega spec/mags/fluxes manipulations """
from functools import wraps
import numpy
import tables
from ..config import __ROOT__


class Vega(object):
    """ Handle vega spec manipulations """

    def __init__(self, source='{0}/libs/vega.hd5'.format(__ROOT__)):
        self.source = source
        self.hdf = None

    def __enter__(self):
        if self.hdf is None:
            self.hdf = tables.openFile(self.source)
        return self

    def __exit__(self,  *exc_info):
        if not self.hdf is None:
            self.hdf.close()
            self.hdf = None
        return False

    def getFlux(self, filters):
        """ Return vega abs. fluxes in filters """
        with self as s:
            FNAME  = s.hdf.root.sed.cols.FNAME[:]
            LUM    = s.hdf.root.sed.cols.LUM[:]
            CWAVE  = s.hdf.root.sed.cols.CWAVE[:]
        idx = numpy.asarray([ numpy.where( FNAME == k) for k in filters ])
        return numpy.ravel(FNAME[idx]), numpy.ravel(LUM[idx]), numpy.ravel(CWAVE[idx])

    def getMag(self, filters):
        """ Return vega abs. magnitudes in filters """
        with self as s:
            FNAME  = s.hdf.root.sed.cols.FNAME[:]
            MAG    = s.hdf.root.sed.cols.MAG[:]
            CWAVE  = s.hdf.root.sed.cols.CWAVE[:]
        idx = numpy.asarray([ numpy.where( FNAME == k) for k in filters ])
        return numpy.ravel(FNAME[idx]), numpy.ravel(MAG[idx]), numpy.ravel(CWAVE[idx])


def testUnit():
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_WFC3_F475W',
           'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    with Vega() as v:
        vega_f, vega_mag, flamb = v.getSed(filters)
    print vega_f, vega_mag, flamb


def from_Vegamag_to_Flux(lamb, vega_mag):
    def deco(f):
        def vegamagtoFlux(mag, err, mask):
            f = numpy.power(10, -0.4 * (mag + vega_mag))
            e = f * ( 1. - numpy.power(10, -0.4 * err) )
            return f, e, mask

        @wraps(f)
        def wrapper(*args, **kwargs):
            mag, err, mask = f(args[0], args[1], **kwargs)
            return vegamagtoFlux( mag, err, mask )

        return wrapper
    return deco
