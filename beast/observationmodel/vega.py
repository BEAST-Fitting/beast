""" Handle vega spec/mags/fluxes manipulations """
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from functools import wraps
import numpy

import tables
from ..config import __ROOT__


__all__ = ['Vega', 'from_Vegamag_to_Flux']


class Vega(object):
    """ Class that handles vega spectrum and references.
    This class know where to find the Vega synthetic spectrum (Bohlin 2007) in
    order to compute fluxes and magnitudes in given filters

    An instance can be used as a context manager as::

        filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_WFC3_F475W', \
                   'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
        with Vega() as v:
            vega_f, vega_mag, flamb = v.getSed(filters)
        print vega_f, vega_mag, flamb

    """

    def __init__(self, source=None):
        """ Constructor """
        if source is None:
            source = '{0}/vega.hd5'.format(__ROOT__)
        self.source = source
        self.hdf = None

    def __enter__(self):
        """ Enter context """
        if self.hdf is None:
            self.hdf = tables.open_file(self.source)
        return self

    def __exit__(self,  *exc_info):
        """ end context """
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
        idx = numpy.asarray([ numpy.where( FNAME == k.encode('utf-8') )
                              for k in filters ])
        return numpy.ravel(FNAME[idx]), numpy.ravel(LUM[idx]), numpy.ravel(CWAVE[idx])

    def getMag(self, filters):
        """ Return vega abs. magnitudes in filters """
        with self as s:
            FNAME  = s.hdf.root.sed.cols.FNAME[:]
            MAG    = s.hdf.root.sed.cols.MAG[:]
            CWAVE  = s.hdf.root.sed.cols.CWAVE[:]
        idx = numpy.asarray([ numpy.where( FNAME == k) for k in filters ])
        return numpy.ravel(FNAME[idx]), numpy.ravel(MAG[idx]), numpy.ravel(CWAVE[idx])


def xxtestUnit():
    """ Unit test and example usage """
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_WFC3_F475W',
               'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    with Vega() as v:
        vega_f, vega_mag, flamb = v.getSed(filters)
    print(vega_f, vega_mag, flamb)


def from_Vegamag_to_Flux(lamb, vega_mag):
    """ function decorator that transforms vega magnitudes to fluxes (without vega reference) """
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


def from_Vegamag_to_Flux_SN_errors(lamb, vega_mag):
    """ function decorator that transforms vega magnitudes to fluxes (without vega reference) """
    def deco(f):
        def vegamagtoFlux(mag, errp, errm, mask):
            f = 10 ** (-0.4 * (mag + vega_mag))
            fp = 10 ** (-0.4 * (mag - errp + vega_mag))
            fm = 10 ** (-0.4 * (mag + errm + vega_mag))
            return f, fp - f, f - fm, mask

        @wraps(f)
        def wrapper(*args, **kwargs):
            mag, errp, errm, mask = f(args[0], args[1], **kwargs)
            return vegamagtoFlux( mag, errp, errm, mask )

        return wrapper
    return deco
