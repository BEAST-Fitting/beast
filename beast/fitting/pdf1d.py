# class to generate 1D PDFs for many objects all with
#  spare or full nD likelihoods on the same grid of models
import math
import numpy as np


class pdf1d:
    def __init__(
        self, gridvals, nbins, logspacing=False, minval=None, maxval=None, uniqvals=None
    ):
        """
        Create an object which can be used to efficiently generate a 1D pdf
        for an observed object

        Parameters
        ----------
        gridvals : ndarray
            1D `float` array with the values of the quantity for all the grid points
        nbins : int
            number of bins to use for the 1D pdf
        logspacing : bool, optional
            whether to use logarithmic spacing for the bins
        minval, maxval : float, optional
            override the range for the bins. this can be useful to make
            sure that the pdfs for different runs have the same bins
        uniqvals : ndarray, optional
            unique values for the full physics grid
        """
        self.nbins = nbins
        self.n_gridvals = len(gridvals)
        self.logspacing = logspacing

        # grab copy of gridvals that can be edited without messing with original
        tgridvals = np.array(gridvals)
        self.n_indxs = len(gridvals)

        if uniqvals is None:
            uniqvals = np.unique(gridvals)

        if len(uniqvals) < 10:
            print(uniqvals)

        if len(tgridvals) <= 0:
            # this is a hack to just get the code to work when
            # all the possible values are negative and the requested
            # pdf is for log x values
            self.bad = True
            self.bin_vals = np.linspace(0.0, 1.0, num=self.nbins)
        else:
            self.bad = False
            # set bin ranges
            self.min_val = tgridvals.min() if minval is None else minval
            self.max_val = tgridvals.max() if maxval is None else maxval

            # if log spacing requested, make the transformation
            if logspacing:
                self.min_val = math.log10(self.min_val)
                self.max_val = math.log10(self.max_val)
                tgridvals = np.log10(tgridvals)
                uniqvals = np.log10(uniqvals)

            # set bin widths
            if self.nbins > 1:
                self.bin_delta = (self.max_val - self.min_val) / (self.nbins - 1)
            else:
                self.bin_delta = 1

            # set values for the bin middles/edges
            if len(uniqvals) > self.nbins:
                self.bin_vals = self.min_val + np.arange(self.nbins) * self.bin_delta
                self.bin_edges = (
                    self.min_val + (np.arange(self.nbins + 1) - 0.5) * self.bin_delta
                )
            else:
                self.bin_vals = np.array(uniqvals)
                self.bin_edges = np.zeros(self.nbins + 1)
                if self.nbins > 1:
                    self.bin_edges[1:-1] = 0.5 * (
                        self.bin_vals[0:-1] + self.bin_vals[1:]
                    )
                    self.bin_edges[0] = self.bin_vals[0] - (
                        self.bin_edges[1] - self.bin_vals[0]
                    )
                    self.bin_edges[-1] = self.bin_vals[-1] + (
                        self.bin_vals[-1] - self.bin_edges[-2]
                    )
                else:
                    self.bin_edges[0] = 0.95 * self.bin_vals[0]
                    self.bin_edges[1] = 1.05 * self.bin_vals[0]

            # get PDF bin associated with each grid val
            pdf_bin_num = np.digitize(tgridvals, self.bin_edges)

            # array to hold indices for each bin
            # (like the IDL version returned by the histogram function)
            pdf_bin_indxs = []

            used_nindxs = 0
            for i in range(nbins):
                # find the indicies for the current bin
                (cur_bin_indxs,) = np.where(pdf_bin_num == (i + 1))
                used_nindxs += len(cur_bin_indxs)

                # save them
                pdf_bin_indxs.append(cur_bin_indxs)

            # transform the bin edges back to linear spacing if log spacing
            #  was asked for
            if logspacing:
                self.bin_vals = np.power(10.0, self.bin_vals)
                self.bin_edges = np.power(10.0, self.bin_edges)

            self.pdf_bin_indxs = pdf_bin_indxs

            if used_nindxs != self.n_indxs:
                print(used_nindxs, self.n_indxs)
                raise ValueError(
                    "Not all the physics grid model points mapped to 1d ppdf bins - should not happen"
                )

    def gen1d(self, gindxs, weights):
        """
        Compute the 1D posterior PDFs based on the nD probabilities

        Parameters
        ----------
        gindxs : ndarray
            1D `int` array with the indxs of the weights in the full model grid
        weights : ndarray
            1D `float` array with the fit probabilities (likelihood*prior)
            at each grid point

        Returns
        -------
        bin_vals : ndarray
            1D `float` array giving the values at the bin centers
        vals_1d : ndarray
            1D `float` array giving the bin pPDF values
        """

        if self.bad:
            return (self.bin_vals, np.zeros((self.nbins)))
        else:
            _tgrid = np.zeros(self.n_gridvals)
            _tgrid[gindxs] = weights
            _vals_1d = np.zeros(self.nbins)
            for i in range(self.nbins):
                if len(self.pdf_bin_indxs[i]) > 0:
                    _vals_1d[i] = np.sum(_tgrid[self.pdf_bin_indxs[i]])

            return (self.bin_vals, _vals_1d)
