""" Defines a generic interface to observation catalog
    This enables to handle non detections, (upper limits one day?), flux and
    magnitude conversions to avoid painful preparation of the dataset

    Data model v2 with limited quantity units handling
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from scipy.interpolate import interp1d

__all__ = ['Observations', 'FakeObs', 'PhotCharact']

class Observations(object):
    """ A generic class that interfaces observation catalog in a standardized way

    Attributes
    ----------
    inputFile: str
        catalog source file

    filters: sequence(str)
        list of filter names (internal standards)

    desc: str, optional
        description of the observations

    badvalue: float, optional
        value that tags a bad measurement that should not be used in the
        fitting.

    nObs: int
        number of observations in the catalog
    """
    def __init__(self, inputFile, desc=None):
        """ Generate a data interface object """
        self.inputFile = inputFile
        self.filters   = None
        self.desc      = desc
        self.readData()
        self.badvalue  = None

    @property
    def nObs(self):
        return self.data.nrows

    def __len__(self):
        return self.nObs

    def __call__(self):
        """ Calling the object will show info """
        self.info()

    def info(self):
        """ Prints some information about the catalog """
        txt = "Data read from {s.inputFile:s}\n"
        if self.desc is not None:
            txt += "Description: {s.desc:s}\n"
        txt += "Number of records: {s.nObs:d}\n\n"
        txt += "Dataset contains:"

        print("Data read from %s " % self.inputFile)
        if self.desc is not None:
            print("Description: %s" % self.desc)
            print("Number of records: %d" % self.nObs)
            print("")
            print("Dataset contains:")

        for k in list(self.data.keys()):
            txt += "\t {0:s}\n".format(k)

        if self.filters is None:
            txt += '\n No filters given yet!'
        else:
            txt += '\n Using Filters: {s.filters}\n'

        print(txt.format(s=self))

    def __getitem__(self, *args, **kwargs):
        """ get item will generate a subsample """
        return self.data.__getitem__(*args, **kwargs)

    def keys(self):
        """ Returns dataset content names """
        return self.data.keys()

    def setDescription(self, txt):
        self.desc = txt

    def setBadValue(self, val):
        self.badvalue = val

    def getFilters(self):
        return self.filters

    def setFilters(self, filters):
        self.filters = filters

    def getMags(self, num, filters):
        raise Exception('Do not use as magnitudes')
        return np.array([ self.data[tt][num] for tt in filters])

    def getErrors(self, num, filters):
        raise Exception('Do not use as magnitudes')
        return np.array([ self.data[tt + 'err'][num] for tt in filters])

    def getFlux(self, num):
        """returns the flux of an observation from the number of counts"""

        flux = np.empty(len(self.filters), dtype=float)
        for ek, ok in enumerate(self.filters):
            flux[ek] = self.data[ok][num]

        return flux

    def getFluxerr(self, num):
        """returns the error on the flux of an observation from the number of counts (not used in the analysis)"""

        fluxerr = np.empty(len(self.filters), dtype=float)

        for ek, ok in enumerate(self.filters):
            fluxerr[ek] = self.data[ok + '_err'][num]

        return fluxerr

    def getObs(self, num=0):
        """ returns the flux"""
        if self.filters is None:
            raise AttributeError('No filter set provided.')

        flux = self.getFlux(num, self.filters)

        return flux

    def getObsWithUncertainties(self, num=0):
        """ returns the flux and uncertainties and the mask of bad values"""
        assert ( not self.filters is None), "No filter set."
        mags = self.getMags(num, self.filters)
        errs = self.getErrors(num, self.filters)

        if not self.badvalue is None:
            mask = (mags >= self.badvalue)
        else:
            mask = np.zeros(len(mags), dtype=bool)

        return mags, errs, mask

    def readData(self):
        """ read the dataset from the original source file """
        from ..external.eztables import AstroTable
        if type(self.inputFile) == str:
            self.data = AstroTable(self.inputFile)
        else:
            self.data = self.inputFile

    def iterobs(self):
        """ yield getObs """
        for k in range(self.nObs):
            yield self.getObs(k)

    def enumobs(self):
        for k in range(self.nObs):
            yield k, self.getObs(k)

#******************
# Code below is not tested/used for sometime (KDG - Jul 2017)
# Not clear if any of this code is needed any longer.
#******************
            
class FakeObs(Observations):

    def getObs(self, num=0, err=0.05):
        assert ( self.filters is not None), "No filter set."
        mags = self.getMags(num, self.filters)
        #errs = np.ones(len(mags), dtype=float) * err
        errs = self.getErrors(num, self.filters)
        if self.badvalue is not None:
            mask = (mags >= self.badvalue)
        else:
            mask = np.zeros(len(mags), dtype=bool)

        return mags, errs, mask

    def readData(self):
        """ read the dataset from the original source file """
        from ..external.eztables import Table
        self.data = Table(self.inputFile)


def gen_FakeObs_from_sedgrid(sedgrid, nrows, err=0.05, 
                             filters=None, save=False):
    from ..external.eztables import Table
    from . import grid
    if type(sedgrid) == str:
        sedgrid = grid.FileSEDGrid(sedgrid)

    inds = np.random.randint(0, high=sedgrid.grid.nrows, size=nrows)
    obsTab = Table()
    if filters is None:
        filters = sedgrid.grid.header.FILTERS.split()
    for e, filt in enumerate(filters):
        errs = np.random.normal(loc=0., scale=err, size=nrows)
        obsTab.addCol(filt, (1. + errs) * sedgrid.seds[inds, e])
        obsTab.addCol(filt + 'err', err * sedgrid.seds[inds, e])
    for key in list(sedgrid.grid.keys()):
        obsTab.addCol(key, sedgrid.grid[key][inds])

    if save is True:
        return obsTab
    else:
        obsTab.write(save, clobber=True, append=False)


class PhotCharact(object):
    def __init__(self, fname, filters):
        self.inputFile = fname
        self.filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W',
                        'HST_ACS_WFC_F475W', 'HST_ACS_WFC_F814W',
                        'HST_WFC3_F110W', 'HST_WFC3_F160W']

        self.pars = { 'magin': 'magin',
                      'comp': 'comp',
                      'bias': 'bias',
                      'berr': 'bias_err',
                      'join': '_'}

        self.interp_bias = { 'kind': 'linear',
                             'axis': -1,
                             'copy': False,
                             'bounds_error': False,
                             'fill_value': 0.0 }

        self.interp_bias_error = { 'kind': 'linear',
                                   'axis': -1,
                                   'copy': False,
                                   'bounds_error': False,
                                   'fill_value': 0.0 }

        self.interp_comp = { 'kind': 'linear',
                             'axis': -1,
                             'copy': False,
                             'bounds_error': False,
                             'fill_value': 0.0 }

        self.readData()
        self._funcs = [ self.getCharactFilterFunctions(k, output='flux') for k in self.filters ]

    def readData(self):
        """ read the dataset from the original source file """
        from ..external.eztables import AstroTable
        from .vega import Vega

        self.data = AstroTable(self.inputFile)

        #data_filters = [ k.split(self.pars['join'])[0] for k in self.data.keys() if k[-5:] == self.pars['magin'] ]

        #Data are in Vega magnitudes
        #  Need to use Vega
        with Vega() as v:
            #self.vega = vega_f, vega_mag, lamb = v.getMag(self.filters)
            self.vega = v.getMag(self.filters)

    def get_filter_index(self, fname):
        for e, k in enumerate(self.filters):
            if k == fname:
                return e

    def getCharactFilterFunctions(self, fname, output='flux'):
        if type(fname) == int:
            _fname = self.filters[fname]
        else:
            _fname = fname

        if output not in ['mag', 'flux']:
            raise ValueError("interfrom must be either mag or flux")

        join = self.pars['join']
        bkey = self.pars['bias']
        ekey = self.pars['berr']
        mkey = self.pars['magin']

        m_val = self.data[join.join([_fname, mkey])]
        b_val = self.data[join.join([_fname, bkey])]
        e_val = self.data[join.join([_fname, ekey])]

        if output == 'mag':
            #interp from mags
            bias_mag_fn = interp1d(m_val, b_val, **self.interp_bias)
            bias_err_mag_fn  = interp1d(m_val, e_val, **self.interp_bias_error)

            return (bias_mag_fn, bias_err_mag_fn)
        else:
            #vegamag to fluxes
            #b_val is a delta_mag, does not need to check the vega ref.
            vega_mag = self.vega[1][self.get_filter_index(_fname)]
            flux_in = np.power(10., -0.4 * (m_val + vega_mag))
            flux_bias = np.power(10., -0.4 * (b_val))
            flux_err = b_val * ( 1. - np.power(10., -0.4 * e_val) )
            bias_flux_fn = interp1d(flux_in, flux_bias, **self.interp_bias)
            bias_err_flux_fn  = interp1d(flux_in, flux_err, **self.interp_bias_error)
            self.flux_in = flux_in
            self.flux_bias = flux_bias
            self.flux_err = flux_err

            return (bias_flux_fn, bias_err_flux_fn)

    def get_bias_of_sed(self, sed, **kwargs):
        if np.ndim(sed) > 1:
            nlamb = np.shape(sed)[1]
            if nlamb != len(self.filters):
                raise ValueError('expecting {0} values per sed, got {1}'.format(len(self.filters), nlamb))
            biases = np.empty(sed.shape, dtype=float)
            errors = np.empty(sed.shape, dtype=float)

            for e, fk in enumerate(self._funcs):
                biases[e, :] = fk[0](sed[e, :])
                errors[e, :] = fk[1](sed[e, :])
        else:
            biases = np.empty(sed.shape, dtype=float)
            errors = np.empty(sed.shape, dtype=float)
            for e, fk in enumerate(self._funcs):
                biases[e] = fk[0](sed[e])
                errors[e] = fk[1](sed[e])

        return biases, errors
