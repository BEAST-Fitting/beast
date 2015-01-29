# class to generate 1D PDFs for many objects all with
#  spare nD likelihoods on the same grid of models

import math
import numpy as np

class pdf1d():
    def __init__(self, gridvals, nbins):
        self.nbins = nbins
        self.n_gridvals = len(gridvals)

        self.min_val = gridvals.min()
        self.max_val = gridvals.max()
        self.bin_delta = (self.max_val - self.min_val)/(self.nbins-1)
        self.bin_vals = self.min_val + np.arange(self.nbins)*self.bin_delta
        self.bin_edges = self.min_val + (np.arange(self.nbins+1) - 0.5)*self.bin_delta

        # get in indices of the grid for each bin in the PDF
        _tpdf_indxs = np.digitize(gridvals, self.bin_edges)
        
        # generate the reverse indices (like the IDL version returned by the histogram function)
        _tgrid_indxs = np.arange(self.n_gridvals)

        self.bin_edges_indxs = np.zeros(self.nbins+1)
        for i in range(nbins):
            _cur_indxs, = np.where(_tpdf_indxs == (i+1)) # find the indicies for the current bin

            self.bin_edges_indxs[i+1] = self.bin_edges_indxs[i] + len(_cur_indxs)
            if len(_cur_indxs) > 0:
                _tgrid_indxs[self.bin_edges_indxs[i]:self.bin_edges_indxs[i+1]] = _cur_indxs

        self.tpdf_indxs = _tgrid_indxs

    def gen1d(self, gindxs, weights):

        _tgrid = np.zeros(self.n_gridvals)
        _tgrid[gindxs] = weights
        _vals_1d = np.zeros(self.nbins)
        for i in range(self.nbins):
            if self.bin_edges_indxs[i] < self.bin_edges_indxs[i+1]:
                _vals_1d[i] = np.sum(_tgrid[self.tpdf_indxs[self.bin_edges_indxs[i]:self.bin_edges_indxs[i+1]]])

        return (self.bin_vals, _vals_1d)

