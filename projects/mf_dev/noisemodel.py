""" Noise model interface v0.1

The idea is to use artificial star tests (ASTs) to characterize the noise
introduced by crowding and selection function.

This version assumes that all bands are independent (which is wrong)
"""
import numpy as np
import tables
from scipy.interpolate import interp1d
from beast.external.eztables import AstroTable
from beast.core.vega import Vega
#from beast.external.ezunits import unit
#from beast.tools.helpers import val_in_unit


class OneD_ASTs_ModelGenerator(object):
    """
    Generate an interface to compute independent band biases, dispersions and
    completenesses given a set of ASTs

    .. note::
        This assumes that ASTs are given in mag!

    Attributes
    ----------
    inputFile: str
        AST input file

    pars: dict
        internal to external name mapping

    interp_bias: dict
        keyword arguments used to define the interpolator on bias values

    interp_bias_error: dict
        keyword arguments used to define the interpolator on bias dispersion values

    interp_comp: dict
        keyword arguments used to define the interpolator on completeness values

    .. note::
        interpolations are 1D interpolations (see: `:func:scipy.interp1d`)
    """
    def __init__(self, fname, filters):
        self.inputFile = fname
        self.setFilters(filters)

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

    def setBadValue(self, val):
        self.badvalue = val

    def getFilters(self):
        return self.filters

    def setFilters(self, filters):
        self.filters = filters

        #Data are in Vega magnitudes
        with Vega() as v:
            _, self.vega, _ = v.getMag(self.filters)

    @property
    def nObs(self):
        return self.data.nrows

    def __len__(self):
        return self.nObs

    def readData(self):
        """ read the dataset from the original source file """
        self.data = AstroTable(self.inputFile)

    def get_filter_index(self, fname):
        for e, k in enumerate(self.filters):
            if k == fname:
                return e

    def getCharactFilterFunctions(self, fname, output='flux'):
        """ Return interpolation functions

        Parameters
        ----------
        fname: str or int
            if str, filter name, filter index otherwise

        output: str in ['mag', 'flux']
            Interpolation function predicts either `mag` or `flux`

        Returns
        -------
        bias_fn: callable
            bias interpolation function

        bias_err_fn: callable
            bias dispersion interpolation function
        """
        if type(fname) == int:
            _fname = self.filters[fname]
        else:
            _fname = fname

        if output not in ['mag', 'flux']:
            raise ValueError("interpolation `output` must be either mag or flux")

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
            vega_mag = self.vega[self.get_filter_index(_fname)]

            flux_in = 10. ** (-0.4 * (m_val + vega_mag))
            flux_bias = 10 ** (-0.4 * (b_val))
            flux_err = b_val * ( 1. - 10 ** (-0.4 * e_val) )
            bias_flux_fn = interp1d(flux_in, flux_bias, **self.interp_bias)
            bias_err_flux_fn  = interp1d(flux_in, flux_err, **self.interp_bias_error)

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


def get_noisemodelcat(filename):
    """
    returns the noise model

    Parameters
    ----------
    filename: str
        file containing the outputs from OneD_ASTs_ModelGenerator

    Returns
    -------
    table: pytables.Table
        table containing the elements of the noise model
    """
    return tables.openFile(filename)
