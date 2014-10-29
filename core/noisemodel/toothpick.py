"""
Toothpick noise model assumes that every photometric band is independent from the others.

The following package implements two classes that corresponds to two variants of input AST data.

the first :class:`MultiFilterASTs` assumes that all AST information is compiled
into on single table, in which one entry corresponds to one artificial star and recovered values

the second :class:`perCameraASTs` assumes that the information is split into
multiple tables, and implements the equivalent of multiple instances of
:class:`MultiFilterASTs` in parallel to calculate the model.

Method
------
Statistics are computed for each input artificial star.  For each input star,
we find the k-NN in input magnitude space and compute their mean and variance
also in magnitude and completeness fractions.

The method was updated in late Oct 2014 by KDG.  The mean bias and standard deviation about
this bias.  In addition, as the calculation is done in magnitude space, sources that are
not recovered (magnitudes above the cutoff, e.g. all 99.99 mags) are used only to calculate
the completeness, nothing else.

Finally, we return the flux conversion of the statistics: bias and stddev
dispersion per input star averaged over k-NN.

TODO: +++ not recovered if delta_mag > 0.8 (factor of 2 in flux)
  +++ change to do the calcuation in flux space as this is more accurate (linear instead of log averaging)
"""
import numpy as np

from .noisemodel import NoiseModel
from ..vega import Vega

from .helpers import _prepare_x, nearest_neighbors, toFlux, convert_dict_to_structured_ndarray
from ...tools.pbar import Pbar


class MultiFilterASTs(NoiseModel):
    """ Implement a noise model for which input information of ASTs are
    provided as one single table

    Attributes
    ----------
    astfile: str
        file containing the ASTs

    filters: sequence(str)
        sequence of filter names
    """
    def __init__(self, astfile, filters, *args, **kwargs):
        NoiseModel.__init__(self, astfile, *args, **kwargs)
        self.setFilters(filters)
        if not 'pass_mapping' in kwargs:
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
                               completeness_mag_cut=80, name_prefix=None,
                               asarray=False):
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

        asarray: bool
            if set returns a structured ndarray instead of a dictionary

        Returns
        -------
        d: dict or np.recarray
            dictionary or named array containing the statistics


        Method
        ------
        Statistics are computed for each input artificial star.
        For each input star, we find the k-NN in input magnitude space and compute
        their mean and variance also in magnitude and completeness fractions.

        Finally, we return the flux conversion of the statistics: bias and stddev
        dispersion per input star averaged over k-NN.
        """
        if name_prefix is None:
            name_prefix = ''
        else:
            if name_prefix[-1] != '_':
                name_prefix += '_'

        # number of good ASTs (recovered)
        n_asts = len(mag_in)

        ### compute the completeness using the full set of ASTs
        # setup the nearest neighbors (only recovered ASTs)
        mag_in_prep = _prepare_x(mag_in)
        NN_mag = nearest_neighbors(mag_in_prep, k=k, eps=eps)

        # storage for the completeness
        completeness = np.zeros(n_asts, dtype=float)

        norm_val = float( 1. / np.sqrt( float(k) - 1.))
        for i in range(n_asts):
            nn_i = mag_out[NN_mag[i]]
            # only recovered stars are considered
            # recovered =      stars < completeness_mag_cut
            #             and  stars with |flux_out / flux_in| < 2
            #                         <=> |mag_out - mag_in| < 0.76 (0.756...)
            ind = (nn_i < completeness_mag_cut) & (np.abs(nn_i - mag_in[NN_mag[i]] < 0.76))
            completeness[i] = np.sum(ind)

        completeness /= float(k)

        # get the ASTs that were recovered
        #  only bias and stddev values for these sources will be calculated
        good_indxs,= np.where(mag_out < completeness_mag_cut)

        # number of good ASTs (recovered)
        n_asts = len(good_indxs)

        ### compute the bias and bias_stddev from the good ASTs (recovered)
        mag_in_good = mag_in[good_indxs]
        mag_out_good = mag_out[good_indxs]

        # setup the nearest neighbors (only recovered ASTs)
        mag_in_prep = _prepare_x(mag_in_good)
        NN_mag = nearest_neighbors(mag_in_prep, k=k, eps=eps)

        # setup variables to store the average bias and std dev around that bias
        ave_bias = np.zeros(n_asts, dtype=float)
        std_bias = np.zeros(n_asts, dtype=float)

        # compute the bias for each AST
        bias_mag = mag_out_good - mag_in_good

        norm_val = float( 1. / np.sqrt( float(k) - 1.))
        for i in range(n_asts):
            nn_i = mag_out_good[NN_mag[i]]
            bias_i = bias_mag[NN_mag[i]]
            ave_bias[i] = bias_i.sum(axis=0)/float(k)
            std_bias[i] = np.sqrt( ((bias_i - ave_bias[i]) ** 2).sum(axis=0))

        std_bias *= norm_val

        # compute the final values for output
        mag_inup = mag_in_good + std_bias
        mag_indown = mag_in_good - std_bias

        f_in = toFlux(mag_in_good)
        f_out = toFlux(mag_in_good + ave_bias)
        f_up = toFlux(mag_indown)
        f_down = toFlux(mag_inup)
        std_flux = 0.5 * (f_up - f_down)

        d = {name_prefix + 'MAG_STD': std_bias,
             name_prefix + 'MAG_BIAS': ave_bias,
             name_prefix + 'MAG_IN': mag_in_good,
             name_prefix + 'MAG_OUT': mag_out_good,
             name_prefix + 'FLUX_STD': std_flux,
             name_prefix + 'FLUX_BIAS': f_out - f_in,
             name_prefix + 'FLUX_IN': f_in,
             name_prefix + 'FLUX_OUT': f_out,
             name_prefix + 'COMPLETENESS': completeness[good_indxs]}

        if asarray:
            return convert_dict_to_structured_ndarray(d)
        else:
            return d

    def _compute_stddev_wrong(self, mag_in, mag_out, k=10, eps=0,
                              completeness_mag_cut=80, name_prefix=None,
                              asarray=False):
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

        asarray: bool
            if set returns a structured ndarray instead of a dictionary

        Returns
        -------
        d: dict or np.recarray
            dictionary or named array containing the statistics


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

        norm_std = float( 1. / np.sqrt( float(k) - 1.))
        for i in range(n_asts):
            nn_i = mag_out[NN_mag[i]]
            # only recovered stars are considered
            # recovered =      stars < completeness_mag_cut
            #             and  stars with |flux_out / flux_in| < 2
            #                         <=> |mag_out - mag_in| < 0.76 (0.756...)
            ind = (nn_i < completeness_mag_cut) & (np.abs(nn_i - mag_in[NN_mag[i]] < 0.76))
            std_mag[i] = np.sqrt(( (nn_i[ind] - mag_in[i]) ** 2).sum(axis=0))
            completeness[i] = np.sum(ind)

        std_mag *= norm_std
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

        if asarray:
            return convert_dict_to_structured_ndarray(d)
        else:
            return d

    def fit(self, k=10, eps=0, completeness_mag_cut=80, progress=True):
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

        progress: bool, optional
            if set, display a progress bar

        .. see also: :func:`_compute_stddev`
        """

        shape = len(self.data), len(self.filters)

        self._fluxes = np.empty( shape, dtype=float)
        self._biases = np.empty( shape, dtype=float)
        self._sigmas = np.empty( shape, dtype=float)
        self._compls = np.empty( shape, dtype=float)

        if progress is True:
            it = Pbar(desc='fitting model').iterover(self.filters)
        else:
            it = self.filters

        for e, filterk in enumerate(it):

            mag_in = self.data[filterk + '_in']
            mag_out = self.data[filterk + '_out']

            d = self._compute_stddev(mag_in, mag_out, k=k, eps=eps,
                                     completeness_mag_cut=completeness_mag_cut)

            self._fluxes[:, e] = d['FLUX_IN'] * self.vega_flux[e]
            self._sigmas[:, e] = d['FLUX_STD'] * self.vega_flux[e]
            self._biases[:, e] = d['FLUX_BIAS'] * self.vega_flux[e]
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

    def interpolate(self, sedgrid, progress=True):
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

        progress: bool, optional
            if set, display a progress bar
        """
        flux = sedgrid.seds
        N, M = flux.shape

        if M != len(self.filters):
            raise AttributeError('the grid of models does not seem to be defined with the same number of filters')

        bias = np.empty((N, M), dtype=float)
        sigma = np.empty((N, M), dtype=float)
        compl = np.empty((N, M), dtype=float)

        if progress is True:
            it = Pbar(desc='Evaluating model').iterover(range(M))
        else:
            it = range(M)

        for i in it:

            _fluxes = self._fluxes[:, i]
            arg_sort = np.argsort(self._fluxes[:, i])
            _fluxes = _fluxes[arg_sort]

            bias[:, i] = np.interp(flux[:, i], _fluxes, self._biases[arg_sort, i] )
            sigma[:, i] = np.interp(flux[:, i], _fluxes, self._sigmas[arg_sort, i])
            compl[:, i] = np.interp(flux[:, i], _fluxes, self._compls[arg_sort, i])

        return (bias, sigma, compl)

    def __call__(self, sedgrid, **kwargs):
        return self.interpolate(sedgrid, **kwargs)


class perCameraASTs(NoiseModel):
    """ Implement a noise model for which input information of ASTs are
    provided as multiple tables

    Attributes
    ----------
    astfiles: sequence(str)
        files containing the ASTs

    filters: sequence(sequence(str))
        sequence of sequence of filter names (one per astfile)
    """

    def __init__(self, astfiles, filters, *args, **kwargs):
        self.models = [ MultiFilterASTs(astfile, filts, pass_mapping=True) for astfile, filts in zip(astfiles, filters) ]
        self.set_data_mappings()

    def set_data_mappings(self):
        """ hard code mapping directly with the interface to PHAT-like ASTs

        .. note::

            it makes it trivial to update this function for other input formats
        """
        #TODO: update the mapping to stick to the initial PHAT version
        for filts, model in zip(self.filters, self.models):
            for k in filts:
                try:
                    #self.data.set_alias(k, k.split('_')[-1].lower() + '_rate')
                    model.data.set_alias(k + '_out', k.split('_')[-1].upper() + '_VEGA')
                    model.data.set_alias(k + '_in', k.split('_')[-1].upper() + '_IN')
                except Exception as e:
                    print('Warning: Mapping failed. This could lead to wrong results')
                    raise e

    def fit(self, k=10, eps=0, completeness_mag_cut=80, progress=True):
        """
        Fit one model per camera

        Parameters
        ----------
        k: Integer
            Number of nearest neighbors taken in the standard deviation computation

        eps: non-negative float
            precision on the NN search

        completeness_mag_cut: float
            magnitude at which consider a star not recovered

        progress: bool, optional
            if set, display a progress bar
        """
        if progress is True:
            it = Pbar(desc='fitting camera').iterover(self.models)
        else:
            it = self.models

        for model in it:
            model.fit(k=k, eps=eps, completeness_mag_cut=completeness_mag_cut,
                      progress=progress)

    def setFilters(self, filters):
        """ set the filters and update the vega reference for the conversions

        Parameters
        ----------
        filters: sequence or sequences
            list of filters per camera using the internally normalized namings
        """
        self.filters = filters
        for model, filts in zip(self.models, self.filters):
            model.setFilters(filts)

    @property
    def filters(self):
        return [model.filters for model in self.models]

    @property
    def astfile(self):
        return [model.astfile for model in self.models]

    @property
    def vega_flux(self):
        return self.models[0].vega_flux

    def interpolate(self, sedgrid, progress=True):
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

        nfilters = sum([len(model.filters) for model in self.models])

        if M != nfilters:
            raise AttributeError('the grid of models does not seem to be defined with the same number of filters')

        bias = np.empty((N, M), dtype=float)
        sigma = np.empty((N, M), dtype=float)
        compl = np.empty((N, M), dtype=float)

        if progress is True:
            pbar = Pbar(M, desc='Evaluating model')

        for j, model in enumerate(self.models):
            for i in range(len(model.filters)):

                if progress:
                    pbar.update(i + j)

                _fluxes = model._fluxes[:, i + j]
                arg_sort = np.argsort(model._fluxes[:, i + j])
                _fluxes = _fluxes[arg_sort]

                bias[:, i + j] = np.interp(flux[:, i], _fluxes, model._biases[arg_sort, i] )
                sigma[:, i + j] = np.interp(flux[:, i], _fluxes, model._sigmas[arg_sort, i])
                compl[:, i + j] = np.interp(flux[:, i], _fluxes, model._compls[arg_sort, i])

            return (bias, sigma, compl)

    def __call__(self, sedgrid, **kwargs):
        return self.interpolate(sedgrid, **kwargs)
