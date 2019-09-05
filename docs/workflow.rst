
#################
Standard Workflow
#################

The workflow is setup to run the fitting on many sources efficiently by
splitting the full catalog into a number of smaller files.  This allows
distributing the fitting across cores.  There are manual steps to allow
for the refitting, fixing issues, etc without rerunning everything.  This
workflow has been tested on large (e.g., PHAT) and small (e.g. METAL)
datasets.

****************************
Production Conda Environment
****************************

Using a dedicated conda environment for production BEAST runs may be
desirable.  Such an environment provides a way to ensure that
production runs are reproducible by fixing the versions of all the
software used.  The instructions below assume that the `astroconda channel
<https://astroconda.readthedocs.io/>`_ is being used.

Create a conda environment.  Here we name it to include the BEAST version.

  .. code-block:: console

    $ conda create -n beast_v1.3.2 python=3.6

Activate the environment after all the packages are finished installing.

  .. code-block:: console

    $ source activate beast_v1.3.2

Install dependencies using conda (better for speed)

  .. code-block:: console

    $ conda install astropy scipy h5py matplotlib cython

Next, install the BEAST.  You have three options:

Option 1: Use pip to install the production version of the beast (currently v1.3.2)

  .. code-block:: console

    $ pip install beast==1.3.2

Option 2: Get the latest production branch, which can be ahead of pipy version

  .. code-block:: console

    $ pip install git+https://github.com/BEAST-Fitting/beast.git@v1.x

Option 3: If you'll be doing development, fork the beast (as described
`here <https://beast.readthedocs.io/en/latest/beast_development.html>`_\),
navigate into the first `beast` folder, and do this command.  Any changes
you make will be immediately reflected in your calls to the BEAST code. Note that
you can make separate environments for development and production modes.

  .. code-block:: console

    $ python setup.py develop

The BEAST production version is now ready for use.  Note, you need to
activate this conda environment every time you want to use this installed
version.

*****
Setup
*****

Working location
================

Setup a working location, usually a subdirectory. For reference, a
template is the 'metal_production' subdirectory in beast/examples.

In this location, at a minimum you will need the following files:

  * run_beast_production.py: a "production" version of run_beast.py
        - Provides commandline options for sub region files
  * beast_production_wrapper.py: assembles the commands below into a script
  * datamodel_template.py: a "production" version of datamodel.py that
    will have name/filter fields automatically filled in by beast_production_wrapper
  * symbolic link to the beast directory in the beast repository

  .. code-block:: console

     $ ln -s /location/beast/beast/ beast



Datamodel_template.py
=====================

Before running the BEAST, you will need to modify this file to specify
the required parameters for generating models and fitting data.
These parameters are described in the beast :ref:`setup documentation
<beast_setup_datamodel>`.  The fields for project, obsfile, astfile,
filters, and basefilters will be filled in by beast_production_wrapper.py.

***************************
Beast_production_wrapper.py
***************************

This is a wrapper for each of the commands described below.  You may
choose to run each of those commands individually, but this
conveniently packages them into one file.  If you use this wrapper, you
should edit several items in the file:

  * field_names: used to identify photometry files and create BEAST files
  * gst_filter_names: labels for the filters used in your photometry
    file (e.g., 'X_RATE')
  * beast_filter_names: the corresponding long names used by the BEAST
  * settings for the source density map: pixel size, filter, magnitude
    range
  * settings for the background map: pixel dimensions, reference image
  * settings for splitting the catalog by source density: filter,
    number of sources per file
  * settings for the trimming/fitting batch scripts: number of files, nice level

You can (and should!) read about the individual functions below before
running beast_production_wrapper:

  .. code-block:: console

     $ run beast_production_wrapper.py

The first thing it does is use datamodel_template.py to create a
datamodel.py file.  This will be imported as needed in the functions
called by the wrapper.  As noted above, five of the datamodel fields
will be updated, so ensure that the other fields in
datamodel_template.py have the desired values.

The wrapper will proceed through each of the functions below.  At
three points, you will need to manually run things independently of
the wrapper.  It will not continue running subsequent functions until
it finds that the necessary steps have been taken.

  * Creating ASTs (if a fake star catalog doesn't exist)
  * running the batch trimming scripts
  * running the batch fitting scripts

Once you have completed each of these, run the wrapper again.  It will
skip past the steps that it has already processed, and resume at the point
where you left off.  In the case of the batch scripts, if you only
partially completed them, it will re-generate new scripts for the
remaining trimming/fitting (and tell you which ones are new), and
pause again.

Note of warning: if you are using this wrapper for multiple fields,
check that the proper version of datamodel.py is in place before
running the batch trimming/fitting scripts.  For instance, if you have
recently used the wrapper to do part of the processing for field_A,
and you want to start the batch fitting script for field_B, re-run the
wrapper for field_B to make sure that datamodel.py refers to the
information for field_B.


****
Data
****

The data need to have source density information added as it is common
for the observation model (scatter and bias) to be strongly dependent
on source density due to crowding/confusion noise.  The background may
also be important in some situations, so there is code to calculate it as well.

Adding source density or background to observations
===================================================

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
a pixel scale of 5 arcsec using the 'datafile.fits' catalog.

  .. code-block:: console

     $ ./beast/tools/create_background_density_map.py sourceden -catfile datafile.fits --pixsize 5.

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

 .. code-block:: console

    $ ./beast/tools/subdivide_obscat_by_source_density.py --n_per_file 6250 \
             --sort_col F475W_RATE datafile_with_sourceden.fits



Adding background to observations
=================================

Create a new version of the observations that includes a column with the
background level.  This is done by calculating the median background for
stars that fall in each spatial bin.  The code will output a new catalog, an
hdf5 file with the background maps and grid information, and some
diagnostic plots.

Command to create the observed catalog with background column with a 15x15 pixel array using the 'datafile.fits' catalog and the 'image.fits' reference image.

  .. code-block:: console

     $ ./beast/tools/create_background_density_map.py background -catfile datafile.fits --npix 15 \
	     -reference image.fits

Plotting the background map onto a reference image
--------------------------------------------------

To check if the background (or source density) map makes sense, the 'tileplot' subcommand of the
same script can be used. If the output of one of the previous commands was 'map_name.hd5', then use

  .. code-block:: console

     $ ./beast/tools/create_background_density_map.py tileplot map_name.hd5 -image image.fits --colorbar 'background'

*****
Model
*****

Physics model
=============

Generate the full physics model grid.  Needed for the fitting and generation of
the artificial star test (AST) inputs.  The '0 0' arguments are dummy values.

  .. code-block:: console

     $ ./run_beast_production.py -p obscat.fits 0 0

Observation model
=================

The observation model is generally based on artificial star tests (ASTs).
ASTs are artificial sources inserted into the observations and extracted with
the same software that was used for the observed photometry catalog.
This ensures that the observation model has the same selection
function as the data.

There are 3 different flavors of observation models.

1. 'Splinter': A very simple (and likely not very good) model that assumes
   the noise is a fraction of the model SED flux and there is no bias.
   No ASTs are used.
2. 'Toothpick':  The AST results are assumed to be independent between
   different bands (even if they are not).  The ASTs results are binned
   in log(flux) bins and the average bias and standard deviation is tabulated
   and used to compute the bias and noise for each model in the physics grid.
3. 'Trunchen': The covariance between bands is measured using the AST results.
   The input AST SEDs are assumed to have been chosen from the BEAST
   physics model grid and are expected to sparsely sample the full model
   grid. The ASTs should be run simultaneously with all bands and it assumed that
   there are multiple ASTs run for the same model.  The covariance
   between the bands is approximated with a multi-variate Gaussian.
   The bias and a multi-variate Gaussian is computed for each model in the
   physic grid by interpolating between the sparse grid computed from the AST
   results.

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

  .. code-block:: console

     $ ./run_beast_production.py -o datafile.fits 0 0

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
by producing multiple trimmed models with a single read.  A specific tool is
provided to setup batch files for this trimming and to do the actual
trimming.

This code sets up batch files for submission to the 'at' queue on linux
(or similar) systems.  The projectname (e.g., 'PHAT') provides a portion
of the batch file names.  The datafile and astfile are the observed photometry
file (not sub files) and file with the ASTs in them.  The optional input
seds_fname can be used to specify the file with the physics model grid,
which overrides the default filename when you wish to use one model grid
for multiple fields. A subdirectory in the project directory is created with
a joblist file for submission to the batch queue and smaller files used by
the trimming code.

The joblist file can be split into smaller files if submission to multiple
cores is desired.  Use the 'split' commandline tool.  The optional 'nice'
input allows you to prepend a 'nice' option, expecially useful if
you're utilizing shared computing resources.

  .. code-block:: console

     $ ./beast/tools/setup_batch_beast_trim.py projectname datafile.fits \
          astfile.fits --num_subtrim 5 --nice 19

Once the batch files are created, then the joblist can be submitted to the
queue.  The beast/tools/trim_many_via_obsdata.py code is called and trimmed
versions of the physics and observation models are created in the project
directory.

  .. code-block:: console

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

  .. code-block:: console

     $ ./beast/tools/setup_batch_beast_fit.py projectname datafile.fits \
       --num_percore 2 --nice 19

The jobs can be submitted to the batch queue via:

  .. code-block:: console

     $ at -f projectname/fit_batch_jobs/beast_batch_fit_X.joblist now

***************
Post-processing
***************

Create the merged stats file
============================

The stats (catalog of fit parameters) files can then be merged into a single
file for the region.  This only merges the stats output files, but not the
pdf1d or lnp files (see the next section).

  .. code-block:: console

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

  .. code-block:: console

     $ beast/tools/reorder_beast_results_spatial.py
        --stats_filename filebase_stats.fits
        --region_filebase filebase_
        --output_filebase spatial/filebase
        --reg_size 10.0

Condense the multiple files for each spatial region into the minimal set.
Each spatial region will have files containing the stats, pdf1d, and lnp
results for the stars in that region.

  .. code-block:: console

     $ beast/tools/condense_beast_results_spatial.py
        --filedir spatial

You may wish to use these files as inputs for the `MegaBEAST <https://megabeast.readthedocs.io/en/latest/>`_.
