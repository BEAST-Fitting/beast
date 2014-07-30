import numpy as np

from .. import NoiseModel
from .multifilterast import MultiFilterASTs


class perCamereASTs(NoiseModel):

    def __init__(self, astfiles, filters, *args, **kwargs):

        NoiseModel.__init__(self, astfiles, *args, **kwargs)
        self.models = [ MultiFilterASTs(astfile) for astfile in astfiles ]
        self.setFilters(filters)

    def set_data_mappings(self):
        """ hard code mapping directly with the interface to PHAT-like ASTs

        .. note::

            it makes it trivial to update this function for other input formats
        """
        #TODO: update the mapping to stick to the initial PHAT version
        pass

    def fit(self, k=10, eps=0, completeness_mag_cut=80):
        """
        Fit one model per camera
        """
        for model in self.models:
            model.compute_knn_statistics(k=k, eps=eps,
                                         completeness_mag_cut=completeness_mag_cut)

    def setFilters(self, filters):
        """ set the filters and update the vega reference for the conversions

        Parameters
        ----------
        filters: sequence
            list of filters using the internally normalized namings
        """
        self.filters = filters
        for model in self.models:
            model.setFilters(filters)

    @property
    def vega_flux(self):
        return self.models[0].vega_flux

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

    def __call__(self, sedgrid):
        return self.interpolate(sedgrid)

