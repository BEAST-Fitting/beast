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

.. code:: shell

   $ python simulate.py physicsgrid obsgrid outfile \
                --nsim 200 --compl_filter f475w

The output file gives the simulated data in the observed data columns
identified in the physicsgrid file along with all the model parameters
from the physicsgrid file.  The names of the observed columns are
`band_rate` and the units are normalized Vega fluxes (to match how
the observed data are given).

********
Trunchen
********

The code does not handle the trunchen model at this point.  It is
straightforward to extend the code for this model, but it has not
been done yet.
