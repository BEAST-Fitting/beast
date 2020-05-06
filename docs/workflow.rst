.. _beast_standard_workflow:

#################
Standard Workflow
#################

The workflow is set up to run the fitting on many sources efficiently by
splitting the full catalog into a number of smaller files.  This allows
distributing the fitting across cores.  There are manual steps to allow
for the refitting, fixing issues, etc without rerunning everything.  This
workflow has been tested on large (e.g., PHAT) and small (e.g. METAL)
datasets.


*****
Setup
*****

Setup a working location. For reference, there are examples in `BEAST-Fitting/beast-examples <https://github.com/BEAST-Fitting/beast-examples>`.

In this location, you will need a datamodel.py parameter file.
These parameters are described in the BEAST :ref:`setup documentation
<beast_setup_datamodel>`.

The BEAST also has some tools for converting catalogs between file formats,
see :ref:`other tools<other_beast_tools>`.

**********************************
Source Density or Background Level
**********************************

The data generally need to have source density information added as it is common
for the observation model (scatter and bias) to be strongly dependent
on source density due to crowding/confusion noise.  The background may
also be important in some situations, so there is code to calculate it as well.

Adding source density to observations
=====================================

Create a new version of the observations that includes a column with the
source density.  The user chooses one band to use as the reference, and chooses
the magnitude range of sources to use for calculating the source density
(generally, this would be the range over which the catalog is complete).  The
user can also choose a band for which sources that have '[band]_FLAG == 99' are
ignored.

Three files are created.  The prefix is derived from the name of the input
photometry catalog.

* [prefix]_source_den_image.fits: an image of the source density map
* [prefix]_with_sourceden.fits: photometry catalog with an
  additional column that has the source density
* [prefix]_sourceden_map.hd5: contains information about the map grid

Command to create the observed catalog with source density column with
a pixel scale of 5 arcsec using the 'phot_catalog.fits' catalog.

  .. code-block:: console

     $ python -m beast.tools.create_background_density_map sourceden \
       -catfile phot_catalog.fits --pixsize 5.



Adding background to observations
=================================

Create a new version of the observations that includes a column with the
background level.  This is done by calculating the median background for
stars that fall in each spatial bin.  The code will output a new catalog, an
hdf5 file with the background maps and grid information, and some
diagnostic plots.

Command to create the observed catalog with background column with a 15x15 pixel
array using the 'phot_catalog.fits' catalog and the 'image.fits' reference image.

  .. code-block:: console

     $ python -m beast.tools.create_background_density_map background \
	     -catfile phot_catalog.fits --npix 15 -reference image.fits


To check if the background (or source density) map makes sense, the 'tileplot' subcommand of the
same script can be used. If the output of one of the previous commands was 'map_name.hd5', then use

  .. code-block:: console

     $ python -m beast.tools.create_background_density_map tileplot map_name.hd5 \
       -image image.fits --colorbar 'background'


*************
Physics model
*************

Generate the full physics model grid.  This is needed for both the fitting and
for generating the artificial star test (AST) inputs.  Note that you may want to
use a coarser model grid for the ASTs.

If you're creating a model grid that's so large it may not read into memory, you
can use subgrids, which splits the grid into more manageable pieces.

To create a physics model grid with 5 subgrids:

  .. code-block:: console

     $ python -m beast.tools.run.create_physicsmodel --nsubs=5

If you're running the BEAST on a survey in which different fields have different
filters, you may wish to save time by creating a master grid with all possible
filters and just copying out the subset of filters you need for each field.  To
do this, create a `datamodel.py` file with all relevant filters listed in
`filters` and `basefilters`, and run `create_physicsmodel` as above.  Then use
`remove_filters` to create each modified grid.  The list of filters to remove
will be determined by what's present in the input catalog file.  If you're using
subgrids, repeat the command for each subgrid.

  .. code-block:: console

     $ python -m beast.tools.remove_filters.py catfile.fits \
         --physgrid master_physgrid.hd5 --physgrid_outfile new_physgrid.hd5


If you would like to examine some or all of the grid values in a physics model,
you can use the `read_sed_data` function in `tools/read_beast_data.py`.  This
function can also be set to just extract the list of parameter names.


*********************
Artificial Star Tests
*********************

The observation model is based on artificial star tests (ASTs).  More details
about the BEAST AST code components can be found at :ref:`Artificial Star Input
Lists <beast_generating_asts>`.

The BEAST selects SEDs from the physics model grid with a technique that
minimizes the number of ASTs needed to allow the construction of a good
toothpick observation model.  For each band, the range of fluxes
in the model grid is split into bins (default=40, set by datamodel.ast_n_flux_bins),
and models are randomly selected.  The model is retained if there are fewer than
the set number of models (default=50, set by datamodel.ast_n_per_flux_bin) in
each of the relevant flux bins.

  .. code-block:: console

     $ python -m beast.tools.run.make_ast_inputs

While not recommended, it is possible to randomly select SEDs from the
physics model grid.

  .. code-block:: console

     $ python -m beast.tools.run.make_ast_inputs --random_seds

How the sources are placed in the image is determined by the ast_source_density_table
variable.

1. datamodel.ast_source_density_table is set to `filebase_sourceden_map.hd5`:
   For each source density or background bin, randomly place the SEDs
   within pixels of that bin.  Repeat for each of the bins.

2. datamodel.ast_source_density_table = None:
   Randomly choose a star from the photometry catalog, and place the
   artificial star nearby.  Repeat until all SEDs have been placed.

.. note::
   These ASTs should be processed with the same code that was used to extract the
   source photometry.


*******************
Edit/Split Catalogs
*******************

You may wish to remove artifacts from the photometry catalog.  If you do so, the
same criteria must be applied to the AST catalog.

The code to edit catalogs can do three different things:

* **Remove objects without full imaging coverage.** Note that the overlap is
  determined by eliminating sources with a flux of precisely 0 in any band.
  However, any sources with a flux of 0 in all bands are not removed, since
  that would indicate that an artificial star was not recovered (this
  criterion does not affect standard photometry catalogs, which do not have
  any sources with flux=0 in all bands).
* **Remove flagged sources.** This eliminates any source with `[filter]_FLAG=99`
  in the specified filter.  If that source has flux<0, it is not removed,
  because those sources are set by `dolphot` to have flag=99 regardless of
  quality.
* **Create ds9 region files.** If set, it will create a ds9 region file where
  good sources are green and removed sources are magenta.

Command to edit the files, both to remove flagged sources and eliminate sources
that don't have full imaging coverage, and to create ds9 region files:

  .. code-block:: console

    $ python -m beast.tools.cut_catalogs \
          phot_catalog_with_sourceden.fits phot_catalog_cut.fits \
          --input_ast_file ast_catalog.fits \
          --output_ast_file ast_catalog_cut.fits \
          --partial_overlap --region_file --flagged --flag_filter F475W


The observed catalog should be split into separate files for each source
density.  In addition, each source density catalog is split into a set of
sub files to have at most 'n_per_file' sources.  The sources are sorted by
the 'sort_col' flux before splitting to put sources with similar brightness
together.  This splitting into sub files sorted by flux allows for trimming
the BEAST physics+observation model, removing objects that are too bright
or too faint to fit any of the sources in the file.  In addition, this
allows for running the BEAST fitting in parallel with each sub file
on a different core.

Command to split both the catalog and AST files by source density:

  .. code-block:: console

    $ python -m beast.tools.split_catalog_using_map.py phot_catalog_cut.fits \
          ast_catalog_cut.fits phot_catalog_sourceden_map.hd5 --bin_width 1 \
          --n_per_file 6250 --sort_col F475W_RATE


*****************
Observation model
*****************

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
3. 'Truncheon': The covariance between bands is measured using the AST results.
   The input AST SEDs are assumed to have been chosen from the BEAST
   physics model grid and are expected to sparsely sample the full model
   grid. The ASTs should be run simultaneously with all bands and it assumed that
   there are multiple ASTs run for the same model.  The covariance
   between the bands is approximated with a multi-variate Gaussian.
   The bias and a multi-variate Gaussian is computed for each model in the
   physic grid by interpolating between the sparse grid computed from the AST
   results.

The code to compute the observation can be done with or without subgridding, and
with or without source density splitting.  Here are some examples:

  .. code-block:: console

     $ # with source density splitting and no subgridding
     $ python -m beast.tools.run.create_obsmodel --use_sd --nsubs 1
     $ # with source density splitting and 5 subgrids
     $ python -m beast.tools.run.create_obsmodel --use_sd --nsubs 5
     $ # no source density splitting or subgrids
     $ python -m beast.tools.run.create_obsmodel --nsubs 1

If you would like to examine some of all of the values in the observation model,
you can use the `read_noise_data` function in `tools/read_beast_data.py`.


******************
Trimming for speed
******************

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

This code sets up batch files for submission to the 'at' queue on linux or
similar systems (such as slurm).  The projectname (e.g., 'PHAT') provides a portion
of the batch file names.  The datafile and astfile are the observed photometry
file (not sub files) and file with the ASTs in them.  The optional input
seds_fname can be used to specify the file with the physics model grid,
which overrides the default filename when you wish to use one model grid
for multiple fields. A subdirectory in the project directory is created with
a joblist file for submission to the batch queue and smaller files used by
the trimming code.

The joblist file can be split into smaller files if submission to multiple
cores is desired.  Use the 'num_subtrim' commandline tool.  The optional 'nice'
input allows you to prepend a 'nice' option, especially useful if
you're utilizing shared computing resources.

  .. code-block:: console

     $ python -m beast.tools.setup_batch_beast_trim projectname phot_catalog_cut.fits \
          ast_catalog_cut.fits --num_subtrim 5 --nice 19

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

     $ python -m beast.tools.setup_batch_beast_fit.py --num_percore 2 --nice 19 \
           --use_sd 1 --nsubs 5 --pdf2d_param_list Av M_ini logT

The jobs can be submitted to the batch queue via:

  .. code-block:: console

     $ at -f projectname/fit_batch_jobs/beast_batch_fit_X.joblist now

The fitting yields several output files (which are described in detail
:doc:`here <outputs>`):

* `*_stats.fits`: Statistics for each of the fitted and derived parameters,
  including the 16th/50th/84th percentiles, mean, and expectation value
* `*_pdf1d.fits`: Marginalized 1D PDFs for each of the fitted and derived
  parameters
* `*_pdf2d.fits`: Marginalized 2D PDFs for pairs of parameters.  If
  `pdf2d_param_list` is set to `None`, 2D PDFs will not be generated.  The
  default set is the 7 main BEAST parameters, but any parameters in the grid can
  be chosen.
* `*_lnp.hd5`: Sparsely sampled log likelihoods

The contents of the `lnp` file can be easily accessed with the `read_lnp_data`
function in `tools/read_beast_data.py`, which converts the hdf5 file structure
into a dictionary.  If you need the SED grid values associated with the saved
lnP points, use the `get_lnp_grid_vals` function in the same file.


***************
Post-processing
***************

Create the merged stats file
============================

The stats files (catalog of fit parameters) can then be merged into a single
file for the field.  The 1D PDF and lnP files are merged across subgrids, but
not yet across source density or background bins.  Merging 2D PDFs has not yet
been implemented.

  .. code-block:: console

     $ python -m beast.tools.run.merge_files --use_sd 1


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

     $ python -m beast.tools.reorder_beast_results_spatial
        --stats_filename filebase_stats.fits
        --region_filebase filebase_
        --output_filebase spatial/filebase
        --reg_size 10.0

Condense the multiple files for each spatial region into the minimal set.
Each spatial region will have files containing the stats, pdf1d, and lnp
results for the stars in that region.

  .. code-block:: console

     $ python -m beast.tools.condense_beast_results_spatial
        --filedir spatial

You may wish to use these files as inputs for the `MegaBEAST <https://megabeast.readthedocs.io/en/latest/>`_.


**************
Python wrapper
**************

This is a wrapper for each of the commands described above:
`beast/examples/production_runs_2019/beast_production_wrapper.py`

You may choose to run each of the above commands individually, but this
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

You can (and should!) read about the individual functions above before
running beast_production_wrapper:

  .. code-block:: console

     $ python beast_production_wrapper

The first thing it does is use datamodel_template.py to create a
datamodel.py file.  You will need to modify datamodel_template.py file to
specify the required parameters for generating models and fitting data.
datamodel.py will be imported as needed in the functions
called by the wrapper.  Four of the datamodel fields (project, obsfile,
filters, and basefilters) will be filled in by beast_production_wrapper.py,
so ensure that the other fields in datamodel_template.py have the desired values.

The wrapper will proceed through each of the functions above.  At
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

*************
Using `slurm`
*************

Many of the steps described above require considerable computational resources,
especially if your grid is large.  If you're running on `XSEDE <https://www.xsede.org/>`_
or another system that uses the slurm queue, you may wish to use
`write_sbatch_file.py`.  This will create a job file that can be submitted with ``sbatch``.
More information about how this file is constructed can be found in the TACC user guide
`here <https://portal.tacc.utexas.edu/archives/stampede#slurm-job-control>`_.

Here is an example call to `write_sbatch_file.py` that shows some of its
functionality.

 .. code-block:: console

    $ # create submission script
    $ python -m beast.tools.write_sbatch_file \
      'sbatch_file.script' './path/to/job/beast_batch_fit_X.joblist' \
      '/path/to/files/projectname/' \
      --modules 'module load anaconda3' 'source activate beast_v1.3' \
      --queue LM --run_time 2:30:00 --mem 250GB


This creates a file ``sbatch_file.script`` with these contents:

 .. code-block:: console

    #!/bin/bash

    #SBATCH -J beast      # Job name
    #SBATCH -p LM            # Queue name
    #SBATCH -t 2:30:00      # Run time (hh:mm:ss)
    #SBATCH --mem 250GB      # Requested memory

    # move to appropriate directory
    cd /path/to/files/projectname/

    # Load any necessary modules
    # Loading modules in the script ensures a consistent environment.
    module load anaconda3
    source activate beast_v1.3

    # Launch a job
    ./path/to/job/beast_batch_fit_X.joblist


Then the file can be submitted:

 .. code-block:: console

    $ sbatch sbatch_file.script
