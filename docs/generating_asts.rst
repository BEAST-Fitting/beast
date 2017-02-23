
Make Artificial Star Input Lists
================================

Generating BEAST-friendly lists of artificial star tests

1) Run "run_small.py -p" to produce the physics model grid file "project_name_seds.grid.hd5".
2) Run "run_small.py -a"   This will use the datamodel to find everything it needs to make ASTs (filters, limits, SED grid, etc.).  It will produce a list of fake stars in all bands using the datamodel photometry catalog to trim the inputs at the proper magnitudes.  Currently, this script generates fake stars uniformly sampling log(age) space and randomly drawing from the metallicities in the model grid.

Functions
=========

mag_limits: Determines the magnitude limits for the models in each filter in the photometry file.

pick_models:  Samples the model grid and outputs models that fit within the mag limits.

pick_positions: Uses the observed stellar catalog to distribution the artificial stars in a similar spatial pattern to the observed catalog

Parameters
==========
   Input parameters in datamodel needed for creating AST lists:
   The magnitude cuts are set to be 1.0 mag fainter than the 90th percentile magnitude in the photometry catalog in each band
   limits # User can set limits = [number] to change the number of mags fainter than the 90th percentile faintest star in the photometry catalog to be used for the mag cut. (Default = 1)
   limits # # # ... User can set limits [space-separated list of numbers] to set custom faint end limits (one value for each band).
   N_models_per_age # Number of models to pick per age bin (Default = 70). Because this is evenly distributed in log age, it is important that the user understand how many models are necessary for reasonable sampling of reddening and metallicity space.  This typically is at least 70 models per age bin.
   Nfilters # Number of filters that must be above the magnitude limit for an AST to be included in the list (Default = 3)
   Nrealize # Number of Realizations of each included AST to be put into the list. (Default = 20)

Returns
=======

Table of fake star magnitudes for all bands in the datamodel photometry file.
The table will have N_models_per_age * Nrealize lines.
