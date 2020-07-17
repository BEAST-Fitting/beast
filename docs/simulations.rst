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
`tools/simulate_obs.py` can be run directly or using the `beast simulate_obs`
command once the beast has been installed.
Simulations require already created BEAST physics and observation model grids.
The physics model grid includes the ensemble parameters as these are the same
as the BEAST :ref:`beast_priors`.
If a different ensemble model is needed (e.g., with a different SFH), then a
new physics model (and possible observations model) will be needed.
The module
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
`--compl_filter` parameter (default=F475W).
Set `compl_filter=max` to use the max completeness value across all the filters.  
The SEDs are picked weighted by the product of the grid+prior weights
and the completeness from the noisemodel.  The grid+prior weights can be replaced
with either grid or prior weights by explicitly setting the `--weight_to_use`
parameter.

.. code-block:: console

   $ beast simulate_obs physicsgrid obsgrid outfile --nsim 200 --compl_filter=F475W

The output file gives the simulated data in the observed data columns
identified in the physicsgrid file along with all the model parameters
from the physicsgrid file.  The simulated observations in each band are given
as `band_flux` in physical units (ergs cm^-2 s^-1 A^-1),
`band_rate` as normalized Vega fluxes (`band_flux`/vega_flux to match how
the observed data are given), and `band_vega` as vega magnitudes with zero and
negative fluxes given as -99.999.
The physicsgrid values without noise/bias are given as `band_input_flux`,
`band_input_rate`, and `band_input_vega`.

When creating simulated observations, using the standard IMF mass prior will
skew your catalog to lower-mass stars.  If you wish to have similar weights for
stars of all masses, use a flat IMF and a log-flat age prior.  To do this,
set the mass prior to `{'name': 'flat'}` and the age prior to
`{'name': 'flat_log'}` in `beast_settings.txt` before creating the model grid.

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

   $ beast plot_cmd outfile.fits --mag1 F475W --mag2 F814W --mag3 F475W

where `outfile.fits` may be the output from `simulate_obs`.
`mag1`-`mag2` is the color, and `mag3` the magnitude.  If you would like to save
(rather than simply display) the figure, include ``--savefig png`` (or another
preferred file extension), and the figure will be saved as `outfile_plot.png` in
the directory of `outfile.fits`.

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
with different filter sets are exactly the same except for the removed filters.

As an example, to remove the filters F275W and F336W from the simulated
observations contained in 'catfile.fits' and the 'physgrid.hd5'/'obsgrid.hd5'
set of models use the following command.

.. code-block:: console

   $ python remove_filters.py catfile.fits --physgrid physgrid.hd5 \
        --obsgrid obsgrid.hd5 --outbase outbase --rm_filters F275W F336W

New physics/observation model grids and simulated observation files are
created as 'outbase_seds.grid.hd5', 'outbase_noisemodel.grid.hd5', and
'outbase_cat.fits'.
