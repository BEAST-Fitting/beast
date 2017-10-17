
*****************
Standard Workflow
*****************

Describe the assumptions about the data, etc.

inc. location of datamodel.py and link to beast source code

Mention example code location

Workflow is setup to run the fitting on many sources efficiently by
splitting the full catalog into a number of smaller files.  This allows
distributing the fitting across cores.  There are manual steps to allow
for the refitting, fixing issues, etc without rerunning everything.  This
workflow has been tested on large (e.g., PHAT) and small (e.g. METAL)
datasets.

Setup working location
----------------------

Setup a working location, usually a subdirectory

Use the example in the beast repository as a template.

In this location, at a minimum you will need the following files:

  * datamodel.py
  * run_beast.py
  * symbolic link to the beast directory in the beast repository

  .. code:: shell

     $ ln -s /location/beast/beast/ beast

Setup datamodel.py
------------------

  * set the name for the project
  * set the survey name
  * update the 3 lists of filters
  * set the physics model grid parameters

Full physics model grid
-----------------------

Generate the full model grid.  Needed for the fitting and generation of
the artifical star test (AST) inputs.

  .. code:: shell

     $ ./run_beast.py -p

Create the AST input list
-------------------------

TBA

Compute the ASTs
----------------

Done separately with whatever code was used to extract the source photometry.
     
Source density map
------------------
  
Create a new version of the observations that includes a column with the
source density.  This also creates a source density image.

Code to create this source density map with a pixel scale of 5 arcsec using
the 'obscat.fits' file of observations.

  .. code:: shell

     $ ./beast/tools/create_source_density_map.py --pixsize 5 obscat.fits
    
Split up observations by source density
---------------------------------------



  * use tools/subdivide_obscat_by_source_density.py
  * creates smaller files allowing for the fitting grid to be smaller (trimmed)

 .. code:: shell

    $ ./beast/tools/subdivide_obscat_by_source_density.py --n_per_file 6250 \
             -sort_col 
    
Split up the ASTs by source density
-----------------------------------
  
  * TBD

Create the observation models for each source density
-----------------------------------------------------
  
  * TBD

  .. code:: shell

     ./run_beast.py -
    
Trim the full model grid for each source density split file
-----------------------------------------------------------

  * use tools/setup_batch_beast_trim.py
  * creates a set of batch files for submission (use 'at -f filename' to submit)
  
  * trimming done such that models are are much too bright or faint are removed
    as they will always give "zero" likelihood
  * files are sorted by brightness and this allows for more trimming of grid
  * smaller grids mean faster fits

Do the fitting
--------------
  
  * each source density split file run with specific trimmed physics and 
    observation model files

Create the merged stats file
----------------------------

  * use tools/merge_stats_file.py
    
Reorganize the results into spatial region files
------------------------------------------------
  
  * TBD (files need to move from megabeast to beast repository)
  * needed for megabeast as well as most other BEAST work
