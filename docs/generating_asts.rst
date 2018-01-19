###########################
Artificial Star Input Lists
###########################

The BEAST requires artificial star tests (ASTs) to produce a noise model.  The AST input list software generates lists of magnitudes and (if desired) positions for ASTs that can be injected into the observed imaging and then re-photometered to assess the photometric bias, uncertainty, and completeness as a function of the model grid.  The output from this software must be run through the same photometry routine (typically DOLPHOT) as used for the photometry measurements themselves.

Once the input lists have been run through the user's photometry program and each input magnitude has an associated output magnitude (or non-detection value), those results can be used as input ASTs for the building the BEAST noise model.

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
ast_models_selected_per_age : integer
Number of models to pick per age (Default = 70).

ast_bands_above_maglimit : integer
Number of filters that must be above the magnitude limit
for an AST to be included in the list (Default = 3)

ast_realization_per_model : integer
Number of Realizations of each included AST model
to be put into the list. (Default = 20)

ast_maglimit : float (single value or array with one value per filter)

1. option 1: [number] to change the number of mags fainter than the 90th percentile
   faintest star in the photometry catalog to be used for the mag cut.
   (Default = 1)

2. option 2: [space-separated list of numbers] to set custom faint end limits
   (one value for each band).

ast_with_positions :  (bool,optional)
If True, the ast list is produced with X,Y positions.
If False, the ast list is produced with only magnitudes.

ast_pixel_distribution : float (optional)
(Used if ast_with_positions is True), minimum pixel separation between AST
position and catalog star used to determine the AST spatial distribution.

ast_reference_image : string (optional, but required if ast_with_positions
is True and no X and Y information  is present in the photometry catalog)
Name of the reference image used by DOLPHOT when running the measured
photometry.

Returns
=======

Table of fake star magnitudes for all bands in the datamodel photometry file.
The file will be in ascii format in the project directory, and it will have the
name: [project]/[project]_inputAST.txt

The table will have ast_models_selected_per_age * ast_realization_per_model lines.
If ast_with_positions is True then each line will start with 0 1 X Y, which are the first
four columns required by DOLPHOT to define the input star position.
