###########
Simulations
###########

Simulations of observations are useful for testing the sensitivity
of a specific set of observations to BEAST model parameters.  This is
useful for quantifying the BEAST performance for proposed and actual
observations.

This is done using the
`beast.observationmodel.observations.gen_SimObs_from_sedgrid` function.
The script
`tools/simulate_obs.py` provides a commandline interface.  The module
uses already created BEAST physics and observation model grids
by sampling the full nD prior function that is part of the physics
model grid.  The observation model grid provides the information on
the flux uncertainties and biases as well as the completeness.

*********
Toothpick
*********

The files for the physicsgrid and obsgrid files are required inputs to
the script.  The output filename is also required.  Note that the extension
of this file will determine the type of file output (e.g. filebase.fits for
a FITS file).
The number of observations to simulate is given by the `--nsim` parameter.
The filter to use for the completeness function is given by the
`--compl_filter` parameter.

.. code-block:: console

   $ python simulate_obs.py physicsgrid obsgrid outfile \
                --nsim 200 --compl_filter f475w

The output file gives the simulated data in the observed data columns
identified in the physicsgrid file along with all the model parameters
from the physicsgrid file.  The names of the observed columns are
`band_rate` and the units are normalized Vega fluxes (to match how
the observed data are given).

*********
Truncheon
*********

The code does not handle the truncheon model at this point.  While this model
is doable in the BEAST, it has not been done yet due to several potentially
complex modeling questions for actually using it that might impact how the model
is implemented.

********
Plotting
********

To plot a color-magnitude diagram of the simulated observations, a
sample call from the command line may be:

.. code-block:: console

   $ python plot_cmd.py outfile.fits --mag1 F475W --mag2 F814W --magy F475W

where `outfile.fits` may be the output from `tools/simulate_obs.py`.
`mag1`-`mag2` is the color, and `magy` the magnitude.
By default the figure is saved as `outfile_plot.png` in the directory
of outfile.fits.

**************
Remove Filters
**************

One use case for simulations is to test the impact of specific filters
on the BEAST results.  One solution is to create multiple physics/observation
model grids, create simulations from each set of grids, and then fit the
simulations with the BEAST.  A quicker way to do this is to create the
physics/observation grid set with the full set of desired filters, create
the desired simulations, remove filters from the model and simulations as
needed, and then fit with the BEAST.  This has the benefit of the simulations
with different filter sets are exactly the same expect for the removed filters.

As an example, to remove the filters F275W and F336W from the simulated
observations contained in 'catfile' and the 'physgrid'/'obsgrid' set of models
use the following command.

.. code-block:: console

   $ python remove_filters.py catfile physgrid obsgrid outbase \
                --rm_filters F275W F336W

New physics/observation model grids and simulated observation files are
created as 'outbase_sed.grid.hd5', 'outbase_noisemodel.grid.hd5', and
'outbase_cat.fits'.
