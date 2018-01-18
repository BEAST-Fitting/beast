Core of probabilistic computations
==================================

If cython is available, the functions will automatically switch to their cython
compiled versions to be as fast as possible especially on large datasets.
Otherwise python/numpy version are used

Tests have shown that speeds of weave and cython codes are similar if
boundaries and negative indices checks are disabled.

likelihoods
-----------
This package implements many likelihoods based on the common chi2 statistics

	N_chi2            Computes a classic (non-reduced) chi2 with normal errors
	SN_chi2           Computes a chi2 with split errors pr non-symmetric errors
	N_logLikelihood   Computes a normal likelihood (default, symmetric errors)
	SN_logLikelihood  Computes a Split Normal likelihood (asymmetric errors)
	getNorm_lnP       Compute the norm of a log-likelihood (overflow robust)

common
------
This package includes all the optimized statistics independent from the
storage of the likelihoods

	percentile        Compute weighted percentiles of a quantity
	expectation       Compute the expectation value of a random variation

Helpers
-------
Mostly storage specific (HDF5) optimized functions

    arange                      Return evenly spaced values within a given interval.
    best_bins                   get non-uniform binning of unique values of Q.
    compute_uniform_prior       Compute the individual model weights to make a given flat prior.
    get_Q_from_node             returns a quantity from a HDF5 Node given its math expression.
    get_centers_from_bins       returns the bin centers from a list of edges.
    get_nclusters               returns the number of likelihoods stored in a file.
    nice_bins                   Define a prior bin width on a quantity q and
                                    iterate until none of the bins are purely empty.
    Q_expect                    Expectation values of any given grid property Q (incl. expression) but seds
    Q_percentile                Percentile values of any given grid property Q (incl. expression) but seds
    sed_expect                  Sed expectation values
    sed_percentile              SED percentile values

kde
---
Kernel density estimation tools

TODO:
    * [ ] add a prior class
        * [ ] defined by (un-)correlated quantities (from model grid)
        * [ ] compute the values
        * [ ] save/restore the values
        * [ ] need a uniform prior step (coded) in case of complex function

    * [x] add/rename percentile scripts

:author:        MF
:last update:   Fri Jun 14 15:33:09 PDT 2013
