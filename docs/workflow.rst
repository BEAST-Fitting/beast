
#################
Standard Workflow
#################

The workflow is setup to run the fitting on many sources efficiently by
splitting the full catalog into a number of smaller files.  This allows
distributing the fitting across cores.  There are manual steps to allow
for the refitting, fixing issues, etc without rerunning everything.  This
workflow has been tested on large (e.g., PHAT) and small (e.g. METAL)
datasets.

*****
Setup
*****

Working location
================

Setup a working location, usually a subdirectory. For reference, a
template is the 'metal_production' subdirectory in beast/examples.

In this location, at a minimum you will need the following files:

  * datamodel.py
  * run_beast_production.py: a "production" version of run_beast.py
        - Provides commandline options for sub region files
  * symbolic link to the beast directory in the beast repository

  .. code:: shell

     $ ln -s /location/beast/beast/ beast

Datamodel.py
============

Before running the BEAST, you will need to modify this file to specify
the required parameters for generating models and fitting data.
These parameters are described in the beast :ref:`setup documentation <beast_setup_datamodel>`.

****
Data
****

The data need to have source density information added as it is common
for the observation model (scatter and bias) to be strongly dependent
on source density due to crowding/confusion noise.

Adding source density to observations
=====================================

Create a new version of the observations that includes a column with the
source density.  The new observation file includes only sources that have
measurements in all bands (columns that match 'X_RATE').  In theory, sources
without measurements in all bands is the result of non-overlapping observations.
The BEAST is based on fitting sources with the same selection function,
in this case measurements in all bands.

A number of source density images are also created.  These include images
that map the source density of objects with zero fluxes in different bands
(or any band).

Command to create the observed catalog with source density column with
a pixel scale of 5 arcsec using the 'obscat.fits' catalog.

  .. code:: shell

     $ ./beast/tools/create_source_density_map.py --pixsize 5. obscat.fits

Split up observations by source density
---------------------------------------

The observed catalog should be split into separate files for each source
density.  In addition, each source density catalog is split into a set of
sub files to have at most 'n_per_file' sources.  The sources are sorted by
the 'sort_col' flux before splitting to put sources with similar brightness
together.  This splitting into sub files sorted by flux allows for trimming
the BEAST physics+observation model removing objects that are too bright
or too faint to fit any of the sources in the file.  In addition, this
allows for running the BEAST fitting in parallel with each sub file
on a different core.

Command to create the the source density split files

 .. code:: shell

    $ ./beast/tools/subdivide_obscat_by_source_density.py --n_per_file 6250 \
             --sort_col F475W_RATE obscat_with_sourceden.fits

*****
Model
*****

Physics model
=============

Generate the full physics model grid.  Needed for the fitting and generation of
the artificial star test (AST) inputs.  The '0 0' arguments are dummy values.

  .. code:: shell

     $ ./run_beast_production.py -p 0 0

Observation model
=================

The observation model is based on artificial star tests (ASTs).  ASTs are
artificial sources inserted into the observations and extracted with
the same software that was used for the observed photometry catalog.
This ensures that the observation model has the same selection
function as the data.

Create the AST input list
-------------------------

To be added.

Compute the ASTs
----------------

Done separately with the same code that was used to extract the source
photometry.


Split up the ASTs by source density
-----------------------------------

To be added.

Currently the workflow assumes a single AST file for all the source densities.

Create the observation models for each source density
-----------------------------------------------------

To be added.

Create a single observation model
---------------------------------

This assumes that the ASTs do not have a strong dependence on source
density.  This could be a good approximation if the source density does
not change much over the observation area or is low everywhere.
The '0 0' arguments are dummy values.

  .. code:: shell

     $ ./run_beast_production.py -o 0 0

******************
Trimming for speed
******************

Trim the full model grid for each source density split file
===========================================================

The physics+observation model can be trimmed of sources that are so bright or
so faint (compared to min/max flux in the observation file) that they will
by definition produce effectively zero likelihood fits.  Such trimming will
speed up the fitting.

The source density split sub files are organized such that the range of
fluxes is minimized in each sub file.  This allows for trimming and faster
fitting.

The trimming can take significant time to run.  In addition, reading in the
full physics+observation model can be slow and such reading can be minimized
by producing multiple trimmed models with a single read.  A specific tools is
provided to setup batch files for this trimming and to do the actual
trimming.

This code sets up batch files for submission to the 'at' queue on linux
(or similar) systems.  The projectname (e.g., 'PHAT') provides a portion
of the batch file names.  The datafile and astfile are the observed photometry
file (not sub files) and file with the ASTs in them.  A subdirection in the
project directory is created with a joblist file for submission to the batch
queue and smaller files used by the trimming code.

The joblist file can be split into smaller files if submission to multiple
cores is desired.  Use the 'split' commandline tool.

  .. code:: shell

     $ ./beast/tools/setup_batch_beast_trim.py projectname datafile astfile \
       --num_subtrim 5

Once the batch files are created, then the joblist can be submitted to the
queue.  The beast/tools/trim_many_via_obsdata.py code is called and trimmed
versions of the physics and observation models are created in the project
directory.

  .. code:: shell

     $ at -f project/trim_batch_jobs/XX_joblist now

*******
Fitting
*******

The fitting is done for each sub file separately.  Code in the tools directory
can be used to create the needed set of batch files for submission to a queue.
In addition, this code will check and see if the fitting has already been done
or was interrupted for the sub files.  Only sub files that have not been fit or
where the fitting was interrupted will be added to the batch files.  The number
of sub files to be run on each core is a command line argument (the runs will
are serial on the core).

  .. code:: shell

     $ ./beast/tools/setup_batch_beast_fit.py projectname datafile \
       --num_percore 2

The jobs can be submitted to the batch queue via:

  .. code:: shell

     $ at -f projectname/fit_batch_jobs/beast_batch_fit_X.joblist now

***************
Post-processing
***************

Create the merged stats file
============================

The stats (catalog of fit parameters) files can then be merged into a single
file for the region.  This only merges the stats output files, but not the
pdf1d or lnp files (see the next section).

  .. code:: shell

     $ beast/tools/merge_stats_file.py filebase

where the filebase where it is the first portion of the output stats filenames
(e.g., filebase_sdx-x_subx_stats.fits).

Reorganize the results into spatial region files
================================================

The output files from the BEAST with this workflow are organized by source
density and brightness.  This is not ideal for finding sources of interest
or performing ensemble processing.  A more useful organization is by spatial
region.  The large amount of BEAST output information makes it best to have
individual files for each spatial region.  Code to do this spatial reordering
is provided in two parts.  The 1st spatially reorders the results for each
source density/brightness BEAST run into files for each spatial region.  The
2nd condenses the multiple individual files for each spatial region into the
minimal set (stats, pdf1d, and lnp).

Divide each source density/brightness file into files of spatial regions
with 10"x10" pixels.

  .. code:: shell

     $ beast/tools/reorder_beast_results_spatial.py
        --stats_filename filebase_stats.fits
        --region_filebase filebase_
        --output_filebase spatial/filebase
        --reg_size 10.0

Condense the multiple files for each spatial region into the minimal set.
Each spatial region will have files containing the stats, pdf1d, and lnp
results for the stars in that region.

     $ beast/tools/condense_beast_results_spatial.py
        --filedir spatial
