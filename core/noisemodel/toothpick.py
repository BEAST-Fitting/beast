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

Additonal method added that computes the noise model in equally spaced bins in log flux space to
avoid injecting noise when the ASTs grossly oversample the model space (e.g. in the case of
single band ASTs).

Finally, we return the flux conversion of the statistics: bias and stddev
dispersion per input star averaged over k-NN.

TODO:
  +++ perCameraASTs has not been updated - delete?  Does not work with PHAT single camera ASTs - column names duplicated
"""
import math

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

    def _compute_stddev(self, magflux_in, magflux_out, k=10, eps=0,
                               completeness_mag_cut=80, name_prefix=None,
                               asarray=False):
        """
        Computes standard deviation and store the result in a dictionary

        Parameters
        ----------
        magflux_in: ndarray
             input mag or flux

        kagflux_out: ndarray
             output mag or flux

        completeness_mag_cut: float
            magnitude at which consider a star not recovered
            set to -1 if the magflux_out is in fluxes (not magnitudes)

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
        For each input star, we find the k-NN in input flux space and compute
        their mean and variance in flux only.

        """
        if name_prefix is None:
            name_prefix = ''
        else:
            if name_prefix[-1] != '_':
                name_prefix += '_'

        # convert the AST input from magnitudes to fluxes
        # always convert the magflux_in to fluxes (the way the ASTs are reported)
        flux_in = 10 ** (-0.4*magflux_in)

        # convert the AST output from magnitudes to fluxes if needed
        #  this is designated by setting the completeness_mag_cut to a negative number
        #    good_indxs gives the list of recovered sources
        if completeness_mag_cut > 0:
            flux_out = 10 ** (-0.4*magflux_out)
            bad_indxs,= np.where(magflux_out >= completeness_mag_cut)
            flux_out[bad_indxs] = 0.0
        else:
            flux_out = magflux_out

        # number of ASTs (all)
        n_asts = len(flux_in)

        ### compute the completeness using the full set of ASTs
        # setup the nearest neighbors
        flux_in_prep = _prepare_x(flux_in)
        NN_flux = nearest_neighbors(flux_in_prep, k=k, eps=eps)

        # storage for the completeness
        completeness = np.zeros(n_asts, dtype=float)

        for i in range(n_asts):
            nn_i = flux_out[NN_flux[i]]
            # only recovered stars are considered with recovered stars w/ flux != 0
            ind = (nn_i != 0.)
            completeness[i] = np.sum(ind)

        completeness /= float(k)

        # get the ASTs that were recovered
        #  only bias and stddev values for these sources will be calculated 
        good_indxs,= np.where(flux_out != 0.0)
        flux_in = flux_in[good_indxs]
        flux_out = flux_out[good_indxs]

        # number of ASTs
        n_asts = len(flux_in)

        # setup the nearest neighbors
        flux_in_prep = _prepare_x(flux_in)
        NN_flux = nearest_neighbors(flux_in_prep, k=k, eps=eps)

        # setup variables to store the average bias and std dev around that bias
        ave_bias = np.zeros(n_asts, dtype=float)
        std_bias = np.zeros(n_asts, dtype=float)

        # compute the bias for each AST
        bias_flux = flux_out - flux_in

        ave_bias = np.mean(bias_flux[NN_flux], axis=1)
        std_bias = np.std(bias_flux[NN_flux], axis=1)

        d = {name_prefix + 'FLUX_STD': std_bias,
             name_prefix + 'FLUX_BIAS': ave_bias,
             name_prefix + 'FLUX_IN': flux_in,
             name_prefix + 'FLUX_OUT': flux_in + ave_bias,
             name_prefix + 'COMPLETENESS': completeness[good_indxs]}

        if asarray:
            return convert_dict_to_structured_ndarray(d)
        else:
            return d

    def _compute_stddev_bins(self, magflux_in, magflux_out, nbins=30,
                             completeness_mag_cut=80, name_prefix=None,
                             asarray=False):
        """
        Computes standard deviation and store the result in a dictionary

        Parameters
        ----------
        magflux_in: ndarray
             input mag or flux

        kagflux_out: ndarray
             output mag or flux

        completeness_mag_cut: float
            magnitude at which consider a star not recovered
            set to -1 if the magflux_out is in fluxes (not magnitudes)

        nbins: Integer
            Number of logrithmically spaced bins between the min/max values

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
        For each input star, we find the k-NN in input flux space and compute
        their mean and variance in flux only.

        """
        if name_prefix is None:
            name_prefix = ''
        else:
            if name_prefix[-1] != '_':
                name_prefix += '_'

        # convert the AST output from magnitudes to fluxes if needed
        #  this is designated by setting the completeness_mag_cut to a negative number
        #    good_indxs gives the list of recovered sources
        if completeness_mag_cut > 0:
            # first remove cases that have input magnitudes below the cut
            #   not sure why this is possible, but they exist and contain
            #   *no information* as mag_in = mag_out = 99.99
            good_in_indxs, = np.where(magflux_in < completeness_mag_cut)
            if len(good_in_indxs) < len(magflux_in):
                magflux_in = magflux_in[good_in_indxs]
                magflux_out = magflux_out[good_in_indxs]

            # now convert from input mags to normalized vega fluxes
            flux_out = 10 ** (-0.4*magflux_out)
            bad_indxs,= np.where(magflux_out >= completeness_mag_cut)
            flux_out[bad_indxs] = 0.0
        else:
            flux_out = magflux_out

        # convert the AST input from magnitudes to fluxes
        # always convert the magflux_in to fluxes (the way the ASTs are reported)
        flux_in = 10 ** (-0.4*magflux_in)

        # storage the storage of the results
        ave_flux_in = np.zeros(nbins, dtype=float)
        ave_bias = np.zeros(nbins, dtype=float)
        std_bias = np.zeros(nbins, dtype=float)
        completeness = np.zeros(nbins, dtype=float)
        good_bins = np.zeros(nbins, dtype=int)

        # get the indexs to the recovered fluxes
        good_indxs,= np.where(flux_out != 0.0)

        # setup the bins (done in log units due to dynamic range)
        #  add a very small value to the max to make sure all the data is included
        min_flux = math.log10(min(flux_in))
        max_flux = math.log10(max(flux_in)*1.000001)
        delta_flux = (max_flux - min_flux)/float(nbins)
        bin_min_vals = min_flux + np.arange(nbins)*delta_flux
        bin_max_vals = bin_min_vals + delta_flux
        bin_ave_vals = 0.5*(bin_min_vals + bin_max_vals)

        # convert the bin min/max value to linear space for computational ease
        bin_min_vals = 10 ** bin_min_vals
        bin_max_vals = 10 ** bin_max_vals
        bin_ave_vals = 10 ** bin_ave_vals

        for i in range(nbins):
            bindxs, = np.where((flux_in >= bin_min_vals[i]) & (flux_in < bin_max_vals[i]))
            n_bindxs = len(bindxs)
            if n_bindxs > 0:
                bin_flux_in = flux_in[bindxs]
                bin_flux_out = flux_out[bindxs]
                # compute completenss
                g_bindxs, = np.where(bin_flux_out != 0.0)
                n_g_bindxs = len(g_bindxs)
                completeness[i] = n_g_bindxs/float(n_bindxs)
                if n_g_bindxs > 5:
                    ave_flux_in[i] = np.mean(bin_flux_in)
                    bin_bias_flux = bin_flux_out[g_bindxs] - bin_flux_in[g_bindxs]
                    ave_bias[i] = np.mean(bin_bias_flux)
                    std_bias[i] = np.std(bin_bias_flux)
                    good_bins[i] = 1
                    
            #print(i,n_bindxs,ave_flux_in[i],ave_bias[i],std_bias[i],completeness[i])

        # only pass back the bins with non-zero results
        gindxs, = np.where(good_bins == 1)

        d = {name_prefix + 'FLUX_STD': std_bias[gindxs],
             name_prefix + 'FLUX_BIAS': ave_bias[gindxs],
             name_prefix + 'FLUX_IN': bin_ave_vals[gindxs],
             name_prefix + 'FLUX_OUT': bin_ave_vals[gindxs] + ave_bias[gindxs],
             #name_prefix + 'FLUX_IN': ave_flux_in,
             #name_prefix + 'FLUX_OUT': ave_flux_in + ave_bias,
             name_prefix + 'COMPLETENESS': completeness[gindxs]}

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
        self._nasts = np.empty(shape[1], dtype=long)

        if progress is True:
            it = Pbar(desc='fitting model').iterover(self.filters)
        else:
            it = self.filters

        for e, filterk in enumerate(it):

            mag_in = self.data[filterk + '_in']
            magflux_out = self.data[filterk + '_out']
            
            d = self._compute_stddev(mag_in, magflux_out, k=k, eps=eps,
                                     completeness_mag_cut=completeness_mag_cut)

            ncurasts = len(d['FLUX_IN'])
            self._fluxes[0:ncurasts, e] = d['FLUX_IN'] * self.vega_flux[e]
            self._sigmas[0:ncurasts, e] = d['FLUX_STD'] * self.vega_flux[e]
            self._biases[0:ncurasts, e] = d['FLUX_BIAS'] * self.vega_flux[e]
            self._compls[0:ncurasts, e] = d['COMPLETENESS']
            self._nasts[e] = ncurasts

            del d

    def fit_bins(self, nbins=30, completeness_mag_cut=80, progress=True):
        """
        Compute the necessary statistics before evaluating the noise model

        Parameters
        ----------
        completeness_mag_cut: float
            magnitude at which consider a star not recovered

        progress: bool, optional
            if set, display a progress bar

        .. see also: :func:`_compute_stddev`
        """

        shape = nbins, len(self.filters)

        self._fluxes = np.empty( shape, dtype=float)
        self._biases = np.empty( shape, dtype=float)
        self._sigmas = np.empty( shape, dtype=float)
        self._compls = np.empty( shape, dtype=float)
        self._nasts = np.empty(shape[1], dtype=long)

        if progress is True:
            it = Pbar(desc='fitting model').iterover(self.filters)
        else:
            it = self.filters

        for e, filterk in enumerate(it):

            mag_in = self.data[filterk + '_in']
            magflux_out = self.data[filterk + '_out']

            d = self._compute_stddev_bins(mag_in, magflux_out, nbins=nbins,
                                          completeness_mag_cut=completeness_mag_cut)

            ncurasts = len(d['FLUX_IN'])
            self._fluxes[0:ncurasts, e] = d['FLUX_IN'] * self.vega_flux[e]
            self._sigmas[0:ncurasts, e] = d['FLUX_STD'] * self.vega_flux[e]
            self._biases[0:ncurasts, e] = d['FLUX_BIAS'] * self.vega_flux[e]
            self._compls[0:ncurasts, e] = d['COMPLETENESS']
            self._nasts[e] = ncurasts

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

            ncurasts = self._nasts[i]
            _fluxes = self._fluxes[0:ncurasts, i]
            _biases = self._biases[0:ncurasts, i]
            _sigmas = self._sigmas[0:ncurasts, i]
            _compls = self._compls[0:ncurasts, i]

            arg_sort = np.argsort(_fluxes)
            _fluxes = _fluxes[arg_sort]

            bias[:, i] = np.interp(flux[:, i], _fluxes, _biases[arg_sort] )
            sigma[:, i] = np.interp(flux[:, i], _fluxes, _sigmas[arg_sort])
            compl[:, i] = np.interp(flux[:, i], _fluxes, _compls[arg_sort])

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
