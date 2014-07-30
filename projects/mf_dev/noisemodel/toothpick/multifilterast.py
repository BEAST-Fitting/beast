import numpy as np

from ..noisemodel import NoiseModel
from ...vega import Vega

from .helpers import _prepare_x, nearest_neighbors, toFlux, convert_dict_to_structured_ndarray


class MultiFilterASTs(NoiseModel):

    def __init__(self, astfile, filters, *args, **kwargs):
        NoiseModel.__init__(self, astfile, *args, **kwargs)
        self.setFilters(filters)
        self.set_data_mappings()

        self._fluxes = None
        self._biases = None
        self._sigmas = None
        self._compls = None

    def set_data_mappings(self):
        """ hard code mapping directly with the interface to PHAT-like ASTs

        .. note::

            it makes it trivial to update this function for other input formats
        """
        try:
            for k in self.filters:
                #self.data.set_alias(k, k.split('_')[-1].lower() + '_rate')
                self.data.set_alias(k + '_out', k.split('_')[-1].lower() + '_vega')
                self.data.set_alias(k + '_in', k.split('_')[-1].lower() + '_in')
        except Exception as e:
            print(e)
            print('Warning: Mapping failed. This could lead to wrong results')

    def _compute_stddev(self, mag_in, mag_out, k=10, eps=0,
                        completeness_mag_cut=80, name_prefix=None):
        """
        Computes standard deviation and store the result in a dictionary

        Parameters
        ----------
        mag_in: ndarray
            input magnitudes

        mag_out: ndarray
            output magnitudes

        completeness_mag_cut: float
            magnitude at which consider a star not recovered

        k: Integer
            Number of nearest neighbors taken in the standard deviation computation

        eps: non-negative float
            precision on the NN search

        name_prefix: str
            if set, all output names in the final structure will start with this
            prefix.

        Returns
        -------
        d: np.recarray
            named array containing the statistics


        Method
        ------
        Statistics are computed for each input artificial star.
        For each input star, we find the k-NN in input magnitude space and compute
        their mean and variance also in magnitude and completeness fractions.

        Finally, we return the flux conversion of the statistics: bias and stddev
        dispersion per input star averaged over k-NN.

        """
        n_asts = len(mag_in)

        if name_prefix is None:
            name_prefix = ''
        else:
            if name_prefix[-1] != '_':
                name_prefix += '_'

        mag_in_prep = _prepare_x(mag_in)
        NN_mag = nearest_neighbors(mag_in_prep, k=k, eps=eps)

        std_mag = np.zeros(n_asts, dtype=float)
        completeness = np.zeros(n_asts, dtype=float)

        for i in range(n_asts):
            nn_i = mag_out[NN_mag[i]]
            std_mag[i] = np.sqrt( (1. / k) * ( (nn_i - mag_in[i]) ** 2).sum(axis=0))
            completeness[i] = np.sum(nn_i < completeness_mag_cut, axis=0)

        completeness /= float(k)

        mag_inup = mag_in + std_mag
        mag_indown = mag_in - std_mag

        f_in = toFlux(mag_in)
        f_out = toFlux(mag_out)
        f_up = toFlux(mag_indown)
        f_down = toFlux(mag_inup)
        std_flux = 0.5 * (f_up - f_down)

        d = {name_prefix + 'MAG_STD': std_mag,
             name_prefix + 'MAG_BIAS': mag_out - mag_in,
             name_prefix + 'MAG_IN': mag_in,
             name_prefix + 'MAG_OUT': mag_out,
             name_prefix + 'FLUX_STD': std_flux,
             name_prefix + 'FLUX_BIAS': f_out - f_in,
             name_prefix + 'FLUX_IN': f_in,
             name_prefix + 'FLUX_OUT': f_out,
             name_prefix + 'COMPLETENESS': completeness}

        return convert_dict_to_structured_ndarray(d)

    def fit(self, k=10, eps=0, completeness_mag_cut=80):
        """
        Compute the necessary statistics before evaluating the noise model

        Parameters
        ----------
        k: Integer
            Number of nearest neighbors taken in the standard deviation computation

        eps: non-negative float
            precision on the NN search

        completeness_mag_cut: float
            magnitude at which consider a star not recovered

        .. see also: :func:`_compute_stddev`
        """

        shape = len(self.data), len(self.filters)

        self._biases = np.empty( shape, dtype=float)
        self._sigmas = np.empty( shape, dtype=float)
        self._compls = np.empty( shape, dtype=float)

        for e, k in enumerate(self.filters):

            mag_in = self.data[k + '_in']
            mag_out = self.data[k + '_out']

            d = self._compute_stddev(mag_in, mag_out, k=k, eps=eps,
                                     completeness_mag_cut=completeness_mag_cut)

            self._fluxes[:, e] = d['FLUX_IN'] / self.vega_flux[e]
            self._sigmas[:, e] = d['FLUX_STD'] / self.vega_flux[e]
            self._biases[:, e] = d['FLUX_BIAS'] / self.vega_flux[e]
            self._compls[:, e] = d['COMPLETENESS']
            del d

    def setFilters(self, filters):
        """ set the filters and update the vega reference for the conversions

        Parameters
        ----------
        filters: sequence
            list of filters using the internally normalized namings
        """
        self.filters = filters

        # ASTs inputs are in vega mag whereas models are not
        # for optimization purpose: pre-compute
        with Vega() as v:
            # name, vega_flux, lamb = v.getFlux(filters)
            _, vega_flux, _ = v.getFlux(filters)

        self.vega_flux = vega_flux

    def interpolate(self, sedgrid):
        """
        Interpolate the results of the ASTs on a model grid

        Parameters
        ----------
        sedgrid: beast.core.grid type
            model grid to interpolate AST results on

        Returns
        -------
        bias: ndarray
            bias table of the models

        sigma: ndarray
            dispersion table of the models

        comp: ndarray
            completeness table per model
        """
        flux = sedgrid.seds
        N, M = flux.shape

        if M != len(self.filter):
            raise AttributeError('the grid of models does not seem to be defined with the same number of filters')

        bias = np.empty((N, M), dtype=float)
        sigma = np.empty((N, M), dtype=float)
        compl = np.empty((N, M), dtype=float)

        for i in range(M):

            _fluxes = self._fluxes[:, i]
            arg_sort = np.argsort(self._fluxes[:, i])
            _fluxes = _fluxes[arg_sort]

            bias[:, i] = np.interp(flux[:, i], _fluxes, self._biases[arg_sort, i] )
            sigma[:, i] = np.interp(flux[:, i], _fluxes, self._sigmas[arg_sort, i])
            compl[:, i] = np.interp(flux[:, i], _fluxes, self._compls[arg_sort, i])

        return (bias, sigma, compl)
