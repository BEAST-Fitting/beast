
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

Setup datamodel.py
------------------

  * set the name for the project
  * set the survey name
  * update the 3 lists of filters
  * set the physics model grid

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
  
  * use tools/create_source_density_map.py
  * also creates an updated observation file that includes source density as a
    seperate entry

Split up observations by source density
---------------------------------------
  
  * use tools/subdivide_obscat_by_source_density.py
  * creates smaller files allowing for the fitting grid to be smaller (trimmed)

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
