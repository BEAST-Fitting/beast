
Standard Workflow
=================

(assumes that the ASTs already run)

- Setup datamodel.py

  * set the name for the project
  * set the survey name
  * update the 3 lists of filters
  * set the physics model grid

- Create the full physics model grid

  * use 'run_beast.py -p'

- Generate the source density map
  
  * use tools/create_source_density_map.py
  * also creates an updated observation file that includes source density as a
    seperate entry

- Split up observations by source density
  
  * use tools/subdivide_obscat_by_source_density.py
  * creates smaller files allowing for the fitting grid to be smaller (trimmed)

- Split up the ASTs by source density
  
  * TBD

- Create the observation models for each source density
  
  * TBD

- Trim the full model grid for each source density split file

  * use tools/setup_batch_beast_trim.py
  * creates a set of batch files for submission (use 'at -f filename' to submit)
  
  * trimming done such that models are are much too bright or faint are removed
    as they will always give "zero" likelihood
  * files are sorted by brightness and this allows for more trimming of grid
  * smaller grids mean faster fits

- Do the fitting
  
  * each source density split file run with specific trimmed physics and 
    observation model files

- Merge small stats files into one large one

  * use tools/merge_beast_stats.py

- Create the merged stats file

  * use tools/merge_stats_file.py
    
- Reorganize the results into spatial region files
  
  * TBD (files need to move from megabeast to beast repository)
  * needed for megabeast as well as most other BEAST work
