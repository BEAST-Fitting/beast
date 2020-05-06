# class to generate 1D PDFs for many objects all with
#  spare or full nD likelihoods on the same grid of models
import math
import numpy as np


class pdf2d:
    def __init__(
        self,
        gridvals_p1,
        gridvals_p2,
        nbins_p1,
        nbins_p2,
        logspacing_p1=False,
        logspacing_p2=False,
        minval_p1=None,
        maxval_p1=None,
        minval_p2=None,
        maxval_p2=None,
    ):
        """
        Create an object which can be used to efficiently generate a 2D pdf
        for an observed object

        Parameters
        ----------
        gridvals_p1, gridvals_p2 : ndarray
            1D `float` array with the values of the quantity for all the grid points
        nbins_p1, nbins_p2 : int
            number of bins to use for the 1D pdf
        logspacing_p1, logspacing_p2 : bool, optional
            whether to use logarithmic spacing for the bins
        minval_p1, maxval_p1, minval_p2, maxval_p2 : float, optional
            override the range for the bins. this can be useful to make
            sure that the pdfs for different runs have the same bins
        """

        # copy values over
        self.nbins_p1 = nbins_p1
        self.nbins_p2 = nbins_p2
        self.n_gridvals = len(gridvals_p1)  # same as len(gridvals_p2)
        self.logspacing_p1 = logspacing_p1
        self.logspacing_p2 = logspacing_p2

        # grab copies of gridvals that can be edited without messing with originals
        tgridvals_p1 = np.array(gridvals_p1)
        tgridvals_p2 = np.array(gridvals_p2)

        # set bin ranges
        self.min_val_p1 = tgridvals_p1.min() if minval_p1 is None else minval_p1
        self.max_val_p1 = tgridvals_p1.max() if maxval_p1 is None else maxval_p1
        self.min_val_p2 = tgridvals_p2.min() if minval_p2 is None else minval_p2
        self.max_val_p2 = tgridvals_p2.max() if maxval_p2 is None else maxval_p2

        # if log spacing requested, make the transformation
        if logspacing_p1:
            self.min_val_p1 = math.log10(self.min_val_p1)
            self.max_val_p1 = math.log10(self.max_val_p1)
            tgridvals_p1 = np.log10(tgridvals_p1)
        if logspacing_p2:
            self.min_val_p2 = math.log10(self.min_val_p2)
            self.max_val_p2 = math.log10(self.max_val_p2)
            tgridvals_p2 = np.log10(tgridvals_p2)

        # set bin widths
        if self.nbins_p1 > 1:
            self.bin_delta_p1 = (self.max_val_p1 - self.min_val_p1) / (
                self.nbins_p1 - 1
            )
        else:
            self.bin_delta_p1 = 1
        if self.nbins_p2 > 1:
            self.bin_delta_p2 = (self.max_val_p2 - self.min_val_p2) / (
                self.nbins_p2 - 1
            )
        else:
            self.bin_delta_p2 = 1

        # set values for the bin middles/edges
        self.bin_vals_p1 = (
            self.min_val_p1 + np.arange(self.nbins_p1) * self.bin_delta_p1
        )
        self.bin_edges_p1 = (
            self.min_val_p1 + (np.arange(self.nbins_p1 + 1) - 0.5) * self.bin_delta_p1
        )
        self.bin_vals_p2 = (
            self.min_val_p2 + np.arange(self.nbins_p2) * self.bin_delta_p2
        )
        self.bin_edges_p2 = (
            self.min_val_p2 + (np.arange(self.nbins_p2 + 1) - 0.5) * self.bin_delta_p2
        )

        # get PDF bin associated with each grid val
        pdf_bin_num_p1 = np.digitize(tgridvals_p1, self.bin_edges_p1)
        pdf_bin_num_p2 = np.digitize(tgridvals_p2, self.bin_edges_p2)

        # array to hold indices for each bin
        pdf_bin_indxs = [
            [0 for j in range(self.nbins_p2)] for i in range(self.nbins_p1)
        ]

        for i in range(self.nbins_p1):
            for j in range(self.nbins_p2):
                # find the indicies for the current bin
                (cur_bin_indxs,) = np.where(
                    (pdf_bin_num_p1 == (i + 1)) & (pdf_bin_num_p2 == (j + 1))
                )
                # save them
                pdf_bin_indxs[i][j] = cur_bin_indxs

        # transform the bin edges back to linear spacing if log spacing
        #  was asked for
        if logspacing_p1:
            self.bin_vals_p1 = np.power(10.0, self.bin_vals_p1)
            self.bin_edges_p1 = np.power(10.0, self.bin_edges_p1)
        if logspacing_p2:
            self.bin_vals_p2 = np.power(10.0, self.bin_vals_p2)
            self.bin_edges_p2 = np.power(10.0, self.bin_edges_p2)

        self.pdf_bin_indxs = pdf_bin_indxs

    def gen2d(self, gindxs, weights):
        """
        Compute the 2D posterior PDFs based on the nD probabilities

        Parameters
        ----------
        gindxs : ndarray
            1D `int` array with the indxs of the weights in the full model grid
        weights : ndarray
            1D `float` array with the fit probabilities (likelihood*prior)
            at each grid point

        Returns
        -------
        vals_2d : ndarray
            2D `float` array giving the bin pPDF values
        """

        _tgrid = np.zeros(self.n_gridvals)
        _tgrid[gindxs] = weights
        _vals_2d = np.zeros((self.nbins_p1, self.nbins_p2))
        for i in range(self.nbins_p1):
            for j in range(self.nbins_p2):
                if len(self.pdf_bin_indxs[i][j]) > 0:
                    _vals_2d[i, j] = np.sum(_tgrid[self.pdf_bin_indxs[i][j]])

        return _vals_2d
