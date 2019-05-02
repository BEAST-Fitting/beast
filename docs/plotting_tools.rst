##############
Plotting Tools
##############

There are `several scripts
<https://github.com/BEAST-Fitting/beast/tree/master/beast/plotting>`_ for making diagnostic plots.  Some are described here.

  * `plot_chi2_hist.py <https://github.com/BEAST-Fitting/beast/blob/master/beast/plotting/plot_chi2_hist.py>`_: Make a histogram of the best chi2 values (chi2=1 and the median chi2 are marked).  Note that there is no plot of reduced chi2, because it is mathematically difficult to define the number of degrees of freedom.  Inputs are the BEAST stats file and optionally the number of bins to use for the histogram.
  * `plot_cmd.py <https://github.com/BEAST-Fitting/beast/blob/master/beast/plotting/plot_cmd.py>`_: Make a color-magnitude diagram of the observations.  Inputs are the photometry file (which can be a `simulation <https://beast.readthedocs.io/en/latest/simulations.html#plotting>`_) and the three filters.
  * `plot_cmd_with_fits.py <https://github.com/BEAST-Fitting/beast/blob/master/beast/plotting/plot_cmd_with_fits.py>`_: Similar to above, but color-coding the data points using one of the parameters from the BEAST fitting.  Takes three additional inputs: a BEAST stats file, the parameter to use, and whether to apply color after taking the log10 of the parameter.
  * `plot_indiv_fit.py <https://github.com/BEAST-Fitting/beast/blob/master/beast/plotting/plot_indiv_fit.py>`_: For a given star, makes a multi-panel plot that shows the PDFs and best fits of each parameter, as well as an SED (similar to Figure 14 in `Gordon+16 <http://adsabs.harvard.edu/abs/2016ApJ...826..104G>`_).
  * `plot_mag_hist.py <https://github.com/BEAST-Fitting/beast/blob/master/beast/plotting/plot_mag_hist.py>`_: Make histograms of the magnitudes for each band in the photometry catalog.


