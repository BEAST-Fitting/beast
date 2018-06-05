# class to generate 1D PDFs for many objects all with
#  spare or full nD likelihoods on the same grid of models
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math
import numpy as np

class pdf1d():
    def __init__(self, gridvals, nbins, logspacing=False, minval=None, maxval=None):
        """
        Create an object which can be used to efficiently generate a 1D pdf for an observed object

        Parameters
        ----------

        gridvals: array-like
            values of the quantity for all the grid points

        nbins: int
            number of bins to use for the 1D pdf

        logspacing: bool
            whether to use logarithmic spacing for the bins

        minval, maxval: float (optional)
            override the range for the bins. this can be useful to make
            sure that the pdfs for different runs have the same bins
        """
        self.nbins = nbins
        self.n_gridvals = len(gridvals)
        self.logspacing = logspacing

        indxs = np.arange(self.n_gridvals)
        self.n_indxs = len(indxs)
            
        # storage of the grid values to consider
        tgridvals = np.array(gridvals[indxs])

        if len(tgridvals) <= 0:
            # this is a hack to just get the code to work when
            # all the possible values are negative and the requested
            # pdf is for log x values
            self.bad = True
            self.bin_vals = np.linspace(0.0,1.0, num=self.nbins)
        else:
            self.bad = False
            self.min_val = tgridvals.min() if minval is None else minval
            self.max_val = tgridvals.max() if maxval is None else maxval

            # if log spacing requested, make the transformation
            if logspacing:
                self.min_val = math.log10(self.min_val)
                self.max_val = math.log10(self.max_val)
                tgridvals = np.log10(tgridvals)
        
            self.bin_delta = (self.max_val - self.min_val)/(self.nbins-1)
            self.bin_vals = self.min_val + np.arange(self.nbins)*self.bin_delta
            self.bin_edges = self.min_val + \
                             (np.arange(self.nbins+1) - 0.5)*self.bin_delta

            # get in indices of the grid for each bin in the PDF
            _tpdf_indxs = np.digitize(tgridvals, self.bin_edges)
        
            # generate the reverse indices 
            # (like the IDL version returned by the histogram function)
            _tgrid_indxs = np.arange(self.n_indxs)
        
            self.bin_edges_indxs = np.zeros((self.nbins+1), dtype=np.uint64)
            for i in range(nbins):
                # find the indicies for the current bin
                _cur_indxs, = np.where(_tpdf_indxs == (i+1)) 
                _cur_indxs = indxs[_cur_indxs]

                self.bin_edges_indxs[i+1] = self.bin_edges_indxs[i] + \
                                            len(_cur_indxs)
                if len(_cur_indxs) > 0:
                    _tgrid_indxs[self.bin_edges_indxs[i]
                                 :self.bin_edges_indxs[i+1]] = _cur_indxs
                
            # transform the bin edges back to linear spacing if log spacing 
            #  was asked for
            if logspacing:
                self.bin_vals = np.power(10.0,self.bin_vals)
                self.bin_edges = np.power(10.0,self.bin_edges)

            self.tpdf_indxs = _tgrid_indxs

    def gen1d(self, gindxs, weights):
        
        if self.bad:
            return (self.bin_vals, np.zeros((self.nbins)))
        else:
            _tgrid = np.zeros(self.n_gridvals)
            _tgrid[gindxs] = weights
            _vals_1d = np.zeros(self.nbins)
            for i in range(self.nbins):
                if self.bin_edges_indxs[i] < self.bin_edges_indxs[i+1]:
                    _vals_1d[i] = np.sum(_tgrid[self.tpdf_indxs[
                        self.bin_edges_indxs[i]:self.bin_edges_indxs[i+1]]])

            return (self.bin_vals, _vals_1d)

    def gen1d_full(self, weights):

        if self.bad:
            return (self.bin_vals, np.zeros((self.nbins)))
        else:
            _vals_1d = np.zeros(self.nbins)
            for i in range(self.nbins):
                if self.bin_edges_indxs[i] < self.bin_edges_indxs[i+1]:
                    _vals_1d[i] = np.sum(weights[self.tpdf_indxs[
                        self.bin_edges_indxs[i]:self.bin_edges_indxs[i+1]]])

            return (self.bin_vals, _vals_1d)

