"""
Trunchen version of noisemodel
Goal is to compute the full 6-band covariance matrix for each model
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math
import numpy as np

from scipy.spatial import cKDTree

#from .noisemodel_trunchen import NoiseModel
from .noisemodel import NoiseModel
from ..vega import Vega

from ...tools.pbar import Pbar

__all__ = ['MultiFilterASTs']

class MultiFilterASTs(NoiseModel):
    """ Implement a noise model where the ASTs are provided as a single table

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

        # needs updating
        self._input_fluxes = None
        self._biases = None
        self._completenesses = None
        self._cov_matrices = None
        self._corr_matrices = None

    def setFilters(self, filters):
        """ set the filters and update the vega reference for the conversions

        Parameters
        ----------
        filters: sequence
            list of filters using the internally normalized namings
        """
        self.filters = filters

        # ASTs inputs are in vega mag whereas models are in flux units
        #     for optimization purpose: pre-compute
        with Vega() as v:
            _, vega_flux, _ = v.getFlux(filters)

        self.vega_flux = vega_flux

    def _calc_ast_cov(self, indxs, filters, return_all=False):
        """
        The NxN-dimensional covariance matrix and N-dimensional bias vector are
        calculated from M independent ASTs computed for N bands

        Parameters
        ----------
        indxs : index array giving the ASTs assocaited with a single
                model SED
        filters : base filter names in the AST file

        Keywords
        --------
        return_all : True/False
        
        Returns
        -------
        if return_all = False
           (cov_mat, bias, compls)
        else
           (cov_mat, bias, stddevs, corr_mat, diffs, ifluxes, compls)

        cov_mat : NxN dim numpy array
                  covariance matrix in flux units
        bias : N dim numpy vector
               vector of the biases in each filter
        stddevs : N dim numpy vector
                  vector of standard deviations in each filter
        corr_mat : NxN dim numpy array
                   correlation matrix
        diffs : KxN dim numpy vector
                raw flux differences for N filters and K AST instances
        ifluxes : N dim numpy vector
                  input fluxes of the AST in each filter
        compl : float
                AST completeness for this model
        """

        # set the asts for this star using the input index array
        asts = self.data[indxs]

        # now check that the source was recovered in at least 1 band
        #   this replicates how the observed catalog is created
        n_asts = len(asts)
        gtindxs = np.full((n_asts),1)
        for k in range(n_asts):
            cgood = 0
            for cfilter in filters:
                if asts[cfilter+'_VEGA'][k] < 90:
                    cgood = cgood + 1
            gtindxs[k] = cgood

        indxs, = np.where(gtindxs > 0)
        n_indxs = len(indxs)
        if n_indxs <= 5:
            return False

        # completeness
        compl = float(n_indxs)/float(n_asts)

        # setup the variables for output
        n_filters = len(filters)
        ifluxes = np.empty((n_filters),dtype=np.float32)
        diffs = np.empty((n_filters,n_indxs),dtype=np.float32)
        biases = np.empty((n_filters),dtype=np.float32)
        cov_matrix = np.full((n_filters, n_filters),0.0,dtype=np.float32)
    
        for ck, cfilter in enumerate(filters):
            ifluxes[ck] = np.power(10.0,-0.4*asts[cfilter+'_IN'][indxs[0]])* \
                          self.vega_flux[ck]
            # compute the difference vector between the input and output fluxes
            #    note that the input fluxes are in magnitudes and the
            #    output fluxes in normalized vega fluxes
            diffs[ck,:] = asts[cfilter+'_RATE'][indxs]*self.vega_flux[ck] - \
                          ifluxes[ck]
            # compute the bias and standard deviations around said bias
            biases[ck] = np.mean(diffs[ck,:])

        # compute the covariance matrix
        for ck in range(n_filters):
            for dk in range(ck,n_filters):
                for ci in range(n_indxs):
                    cov_matrix[ck,dk] += (diffs[ck,ci] - biases[ck])* \
                                         (diffs[dk,ci] - biases[dk])
                # fill in the symmetric terms
                cov_matrix[dk,ck] = cov_matrix[ck,dk]

        cov_matrix /= (n_indxs - 1)
        stddevs = np.sqrt(np.diagonal(cov_matrix))

        # compute the corrleation matrix
        corr_matrix = np.array(cov_matrix)
        for ck in range(n_filters):
            for dk in range(ck,n_filters):
                if stddevs[ck]*stddevs[dk] > 0:
                    corr_matrix[ck,dk] /= stddevs[ck]*stddevs[dk]
                else:
                    corr_matrix[ck,dk] = 0.0
                # fill in the symmetric terms
                corr_matrix[dk,ck] = corr_matrix[ck,dk]

        if return_all:
            return (cov_matrix, biases, stddevs, corr_matrix, diffs, ifluxes,
                    compl)
        else:
            return (cov_matrix, biases, compl)

    def _calc_all_ast_cov(self, filters, progress=True):
        """
        The covariance matrices and biases are calculated for all the
        independent models in the AST file

        Parameters
        ----------
        filters : filter names for the AST data
        
        Keywords
        --------
        progress: bool, optional
            if set, display a progress bar

        Returns
        -------
        (cov_mats, biases, completenesses, corr_mats, ifluxes)

        cov_mats : KxNxN dim numpy array
                   K AST covariance matrices in flux units
        bias : KxN dim numpy vector
               K vectors of the biases in each filter
        completenesses : K dim numpy vector
                         completeness versus model
        corr_mats : KxNxN dim numpy array
                    K AST correlation matrices
        ifluxes : KxN dim numpy vector
                  K vectors of the input fluxes in each filter
        """

        # find the stars by using unique values of the magnitude values
        #   in filtername
        filtername = filters[-1] + '_IN'
        uvals, ucounts = np.unique(self.data[filtername], return_counts=True)
        n_models = len(uvals)

        # setup the output
        n_filters = len(filters)
        all_covs = np.empty((n_models,n_filters,n_filters),dtype=np.float64)
        all_corrs = np.empty((n_models,n_filters,n_filters),dtype=np.float32)
        all_biases = np.empty((n_models,n_filters),dtype=np.float64)
        all_ifluxes = np.empty((n_models,n_filters),dtype=np.float32)
        all_compls = np.empty((n_models),dtype=np.float32)

        ast_minmax = np.empty((2,n_filters),dtype=np.float64)
        ast_minmax[0,:] = 1e99
        ast_minmax[1,:] = 1e-99

        # loop over the unique set of models and
        # calculate the covariance matrix using the ASTs for this model
        good_asts = np.full((n_models),True)
        if progress is True:
            it = Pbar(desc='Calculating AST Covariance ' + \
                      'Matrices').iterover(list(range(n_models)))
        else:
            it = list(range(n_models))
        for i in it:
            # find all the ASTs for this model
            indxs, = np.where(self.data[filtername] == uvals[i])
            n_asts = len(indxs)

            if n_asts > 5:
                results = self._calc_ast_cov(indxs, filters,
                                             return_all=True)
                if results:
                    all_covs[i,:,:] = results[0]
                    all_biases[i,:] = results[1]
                    all_corrs[i,:,:] = results[3]
                    all_ifluxes[i,:] = results[5]
                    all_compls[i] = results[6]

                    for k in range(n_filters):
                        ast_minmax[0,k] = min(ast_minmax[0,k],all_ifluxes[i,k])
                        ast_minmax[1,k] = max(ast_minmax[1,k],all_ifluxes[i,k])
                else:
                    good_asts[i] = False

        indxs, = np.where(good_asts)

        return (all_covs[indxs,:,:], all_biases[indxs,:], all_compls[indxs], 
                all_corrs[indxs,:,:], all_ifluxes[indxs,:], ast_minmax)

    def process_asts(self, filters):
        """
        Process all the AST results creating average biases and
        covariance matrices for each model SED.
        Also, prep for the interpolation by setting up the kd-tree
        
        Parameters
        ----------
        filters : filter names for the AST data
        
        Returns
        -------
        N/A.
        """
        results = self._calc_all_ast_cov(filters)

        self._cov_matrices = results[0]
        self._biases = results[1]
        self._completenesses = results[2]
        self._corr_matrices = results[3]
        self._input_fluxes = results[4]
        self._minmax_asts = results[5]

        print('building kd-tree...')
        self._kdtree = cKDTree(np.log10(self._input_fluxes))
        print('...done')

    def __call__(self, sedgrid,
                 generic_absflux_a_matrix=None,
                 progress=True):
        """
        Interpolate the results of the ASTs on the model grid

        Parameters
        ----------
        sedgrid: beast.core.grid type
            model grid to interpolate AST results on

        Returns
        -------

        progress: bool, optional
            if set, display a progress bar
        """
        flux = sedgrid.seds
        if generic_absflux_a_matrix is not None:
            model_absflux_cov = False
            if generic_absflux_a_matrix is not None:
                print('using model indepdent absflux cov matrix')
            else:
                print('not using any absflux cov matrix')
        elif ((sedgrid.cov_diag is not None) &
              (sedgrid.cov_offdiag is not None)):
            model_absflux_cov = True
            absflux_cov_diag = sedgrid.cov_diag
            absflux_cov_offdiag = sedgrid.cov_offdiag
            print('using model dependent absflux cov matrix')
        else:
            model_absflux_cov = False
            
        n_models, n_filters = flux.shape
        n_offdiag = (((n_filters**2)-n_filters)/2)

        if n_filters != len(self.filters):
            raise AttributeError('the grid of models does not seem to' + 
                                 'be defined with the same number of filters') 

        biases = np.empty((n_models, n_filters), dtype=np.float64)
        sigmas = np.empty((n_models, n_filters), dtype=np.float64)
        cov_diag = np.empty((n_models, n_filters), dtype=np.float64)
        cov_offdiag = np.empty((n_models, n_offdiag), dtype=np.float64)
        icov_diag = np.empty((n_models, n_filters), dtype=np.float64)
        icov_offdiag = np.empty((n_models, n_offdiag), dtype=np.float64)
        q_norm = np.empty((n_models), dtype=np.float64)
        compls = np.empty((n_models), dtype=float)

        if progress is True:
            it = Pbar(desc='Evaluating model').iterover(list(range(n_models)))
        else:
            it = list(range(n_models))

        for i in it:
            # AST results are in vega fluxes
            cur_flux = flux[i,:]

            # find the 10 nearest neighbors to the model SED
            result = self._kdtree.query(np.log10(cur_flux),10)

            dist = result[0]
            indxs = result[1]

            # check if the distance is very small, set to a reasonable value
            tindxs, = np.where(dist < 0.01)
            if len(tindxs) > 0:
                dist[tindxs] = 0.01

            # compute the interpolated covariance matrix
            #    use the distances to generate weights for the sum
            dist_weights = 1.0/dist
            dist_weights /= np.sum(dist_weights)

            cur_cov_matrix = np.average(self._cov_matrices[indxs,:,:],
                                        axis=0,
                                        weights=dist_weights)

            # add in the absflux covariance matrix
            #   unpack off diagonal terms the same way they were packed
            if model_absflux_cov:
                m = 0
                cur_cov_matrix[n_filters-1,n_filters-1] += \
                                                absflux_cov_diag[i,n_filters-1]
                for k in range(n_filters-1):
                    cur_cov_matrix[k,k] += absflux_cov_diag[i,k]
                    for l in range(k+1,n_filters):
                        cur_cov_matrix[k,l] += absflux_cov_offdiag[i,m]
                        cur_cov_matrix[l,k] += absflux_cov_offdiag[i,m]
                        m += 1
            elif generic_absflux_a_matrix is not None:
                for k in range(n_filters):
                    for l in range(n_filters):
                        cur_cov_matrix[k,l] += (generic_absflux_a_matrix[k,l]*
                                                cur_flux[k]*cur_flux[l])

            # compute the interpolated biases
            biases[i,:] = np.average(self._biases[indxs,:],
                                     axis=0,
                                     weights=dist_weights)

            # compute the interpolated completeness
            compls[i] = np.average(self._completenesses[indxs],
                                  weights=dist_weights)

            # save the straight uncertainties
            sigmas[i,:] = np.sqrt(np.diagonal(cur_cov_matrix))

            # invert covariance matrix
            inv_cur_cov_matrix = np.linalg.inv(cur_cov_matrix)

            # save the diagnonal and packed version of non-diagonal terms
            m = 0
            icov_diag[i,n_filters-1] = inv_cur_cov_matrix[n_filters-1,
                                                          n_filters-1]
            cov_diag[i,n_filters-1] = cur_cov_matrix[n_filters-1,
                                                     n_filters-1]
            for k in range(n_filters-1):
                icov_diag[i,k] = inv_cur_cov_matrix[k,k]
                cov_diag[i,k] = cur_cov_matrix[k,k]
                for l in range(k+1,n_filters):
                    icov_offdiag[i,m] = inv_cur_cov_matrix[k,l]
                    cov_offdiag[i,m] = cur_cov_matrix[k,l]
                    m += 1

            # save the log of the determinat for normalization
            #   the ln(det) is calculated and saved as this is what will
            #   be used in the actual calculation
            #       norm = 1.0/sqrt(Q)
            det = np.linalg.slogdet(cur_cov_matrix)
            #print(det)
            if det[0] <= 0:
                print('something bad happened')
                print('determinant of covarinace matrix is zero or negative')
                print(det)
            q_norm[i] = -0.5*det[1]

        return (biases, sigmas, compls, q_norm, icov_diag, icov_offdiag,
                cov_diag, cov_offdiag)
        
