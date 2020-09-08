##############
Plotting Tools
##############

There are tools for making
diagnostic plots and visualizations.  If they are a command line script, then
the command needed to access them from anywhere with an installed version of
the beast is given along with a description.  If they are a module with
plotting functions, then the location in the repository is given.

Observation Models
------------------

- `beast plot_toothpick_details`:
  Shows the details of the creation of the toothpick observation model including
  the individual AST results.
  For more details see :mod:`~beast.plotting.plot_toothpick_details`.

- `beast plot_noisemodel.py`:
  Plot the bias and uncertainty as a function of flux (similar to Figure 12 in
  `Gordon+16 <https://ui.adsabs.harvard.edu/abs/2016ApJ...826..104G>`_).
  Multiple noise models can be overplotted, as long as they correspond to the
  same SED model grid.
  For more details see :mod:`~beast.plotting.plot_noisemodel`.

- `beast plot_ast_histogram`:
  Make a histogram of the input AST fluxes for each filter.  Optionally, also include
  a histogram of the SED grid for comparison.
  For more details see :mod:`~beast.plotting.plot_ast_histogram`.

- `beast plot_completeness`:
  Make a triangle plot with completeness averaged into 2D plots for each pair
  of parameters, and 1D plots along the diagonal for each individual parameter.
  For more details see :mod:`~beast.plotting.plot_completeness`.

Fitting
-------

- `beast plot_indiv_fit`:
  For a given star, makes a multi-panel plot that shows the PDFs and best fits
  of each parameter, as well as an SED (similar to Figure 14 in
  `Gordon+16 <https://ui.adsabs.harvard.edu/abs/2016ApJ...826..104G>`_).
  For more details see :mod:`~beast.plotting.plot_indiv_fit`.

- `beast plot_chi2_hist`:
  Make a histogram of the best chi2 values (chi2=1 and the median chi2 are
  marked).  Note that there is no plot of reduced chi2, because it is mathematically
  difficult to define the number of degrees of freedom.  Inputs are the BEAST stats
  file and optionally the number of bins to use for the histogram.
  For more details see :mod:`~beast.plotting.plot_chi2_hist`.

- `beast plot_indiv_pdfs`:
  For a given star, makes a triangle plot with all of the 2D PDFs.  Diagonals
  contain the 1D PDFs.
  For more details see :mod:`~beast.plotting.plot_indiv_pdfs`.

- `beast plot_param_err`:
  Reproduce the Figures 16-18 in Gordon et al. 2016. Make 2D histogram of 50th
  percentile values of each parameter against its uncertainty
  (=0.5x(percentile_84th-percentile_16th)) on the left columns. Make H-R Hess
  diagram colored coded by the uncertainty of a given parameter on the right
  column.
  For more details see :mod:`~beast.plotting.plot_param_err`.

- `beast plot_param_recovery`:
  Make a 2D histogram to compare simulated and recovered model parameters
  (similar to Figure 13 in `Gordon+16 <https://ui.adsabs.harvard.edu/abs/2016ApJ...826..104G>`_).
  If given multiple sets of files, can do additional panels to compare across
  noise models.
  For more details see :mod:`~beast.plotting.plot_param_recovery`.

- `beast plot_triangle`:
  Make a triangle/corner plot of all the parameters (p50) against each other.
  Diagonals contain histograms of each parameter.
  For more details see :mod:`~beast.plotting.plot_triangle`.

Data/Simulations
----------------

- 'python -m beast.plotting.plot_filters'
  Plot multiple filter response functions.
  For more details see :mod:`~beast.plotting.plot_filters`.

- `beast plot_cmd`:
  Make a color-magnitude diagram of the observations.  Inputs are the photometry
  file (which can be a `simulation <https://beast.readthedocs.io/en/latest/simulations.html#plotting>`_)
  and the three filters.
  For more details see :mod:`~beast.plotting.plot_cmd`.

- `beast plot_cmd_with_fits`:
  Similar to above, but color-coding the data points using one of the parameters
  from the BEAST fitting.  Takes three additional inputs: a BEAST stats file,
  the parameter to use, and whether to apply color after taking the log10 of the
  parameter.
  For more details see :mod:`~beast.plotting.plot_cmd_with_fits`.

- `beast plot_mag_hist`:
  Make histograms of the magnitudes for each band in the photometry catalog.
  For more details see :mod:`~beast.plotting.plot_mag_hist`.

- `beast.plotting.make_ds9_region_file`:
  Make a ds9 region file from an input fits catalog (use `region_file_fits`) or
  a list of artificial stars (use `region_file_txt`).  Can also choose a
  column+value as a cut to set two different region colors.
  For more details see :mod:`~beast.plotting.make_ds9_region_file`.
