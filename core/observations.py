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
import numpy as np
from scipy.interpolate import interp1d
from .vega import Vega


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


""" Some helpers """
def fluxToMag(flux):
    """ Return the magnitudes from flux values
    INPUTS:
        flux    np.ndarray[float, ndim=N]   array of fluxes
    OUTPUTS:
        mag np.ndarray[float, ndim=N]   array of magnitudes
    """
    return -2.5 * np.log10(flux)


def fluxErrTomag(flux, fluxerr):
    """ Return the magnitudes and associated errors from fluxes and flux error values
    INPUTS:
        flux    np.ndarray[float, ndim=1]   array of fluxes
        fluxerr np.ndarray[float, ndim=1]   array of flux errors
    OUTPUTS:
        mag np.ndarray[float, ndim=1]   array of magnitudes
        err np.ndarray[float, ndim=1]   array of magnitude errors
    """
    mag = fluxToMag(flux)
    return mag, -2.5 * np.log10( 1. - fluxerr / flux )


def magToFlux(mag):
    """ Return the flux from magnitude values
    INPUTS:
        mag np.ndarray[float, ndim=N]   array of magnitudes
    OUTPUTS:
        flux    np.ndarray[float, ndim=N]   array of fluxes
    """
    return 10 ** (-0.4 * mag)


def magErrToFlux(mag, err):
    """ Return the flux and associated errors from magnitude and mag error values
    INPUTS:
        mag np.ndarray[float, ndim=1]   array of magnitudes
        err np.ndarray[float, ndim=1]   array of magnitude errors
    OUTPUTS:
        flux    np.ndarray[float, ndim=1]   array of fluxes
        fluxerr np.ndarray[float, ndim=1]   array of flux errors
    """
    flux = magToFlux(mag)
    return flux, flux * ( 1. - magToFlux(err) )


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
        from ..external.eztables import AstroTable
        if type(self.inputFile) == str:
            self.data = AstroTable(self.inputFile)
        else:
            self.data = self.inputFile

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
        #errs = numpy.ones(len(mags), dtype=float) * err
        errs = self.getErrors(num, self.filters)
        if self.badvalue is not None:
            mask = (mags >= self.badvalue)
        else:
            mask = numpy.zeros(len(mags), dtype=bool)

        return mags, errs, mask

    def readData(self):
        """ read the dataset from the original source file """
        from ..external.eztables import Table
        self.data = Table(self.inputFile)


def gen_FakeObs_from_sedgrid(sedgrid, nrows, err=0.05, distanceModulus=0., filters=None, save=False):
    from ..external.eztables import Table
    from . import grid
    if type(sedgrid) == str:
        sedgrid = grid.FileSEDGrid(sedgrid)

    inds = numpy.random.randint(0, high=sedgrid.grid.nrows, size=nrows)
    obsTab = Table()
    if filters is None:
        filters = sedgrid.grid.header.FILTERS.split()
    for e, filt in enumerate(filters):
        errs = numpy.random.normal(loc=0., scale=err, size=nrows)
        obsTab.addCol(filt, (1. + errs) * sedgrid.seds[inds, e])
        obsTab.addCol(filt + 'err', err * sedgrid.seds[inds, e])
    for key in sedgrid.grid.keys():
        obsTab.addCol(key, sedgrid.grid[key][inds])

    if save is True:
        return obsTab
    else:
        obsTab.write(save, clobber=True, append=False)


class PhotCharact(object):
    def __init__(self, fname, filters):
        self.inputFile = fname
        self.filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W', 'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']

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
            flux_in = numpy.power(10., -0.4 * (m_val + vega_mag))
            flux_bias = numpy.power(10., -0.4 * (b_val))
            flux_err = b_val * ( 1. - numpy.power(10., -0.4 * e_val) )
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
                raise ValueError('expecting {} values per sed, got {}'.format(len(self.filters), nlamb))
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
