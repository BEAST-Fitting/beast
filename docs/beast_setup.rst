####################
Setting Up the BEAST
####################

Basics
======

1) Define project and grid input parameters in `beast_settings.txt`.  (This is the
   default file name used throughout the documentation and examples, but it can
   be named whatever you'd like.)

2) Execute BEAST Run using ``python run_beast.py`` with appropriate task flags

   * Default Full Stack Run: ``python run_beast.py -p -o -t -f``

.. _beast_setup_settings:

BEAST Settings
==============

Before running the BEAST, you will need to modify `beast_settings.txt` to specify
the required parameters for generating models and fitting data. These parameters
(and example values) are described below. An example `beast_settings.txt` can be
found `here <https://github.com/BEAST-Fitting/beast-examples/blob/master/phat_small/beast_settings.txt>`_.

Project Details
---------------

* ``project``: pathname of of working subdirectory.
* ``surveyname``: name of the survey.
* ``filters``: names of photometric filter passbands (matching library names).
* ``basefilters``: short versions of passband names.
* ``obs_colnames``: column names in ``obsfile`` for observed fluxes. The input data MUST be in fluxes, NOT in magnitudes and the fluxes MUST be in normalized Vega units.
* ``obsfile``: filename for input flux data.

Physics Model Grid Definition
-----------------------------

The BEAST generates a grid of dust extinguished stellar models based on input parameters
from `beast_settings.txt`.
See :ref:`BEAST grid inputs <beast_grid_inputs>` for details on model libraries. 
For more on setting up priors, see :ref:`BEAST priors <beast_priors>`.

General Parameters
^^^^^^^^^^^^^^^^^^

* ``n_subgrid``: number of sub-grids to use (1 means no subgrids), useful for when
  the physics model grid is too large to read into memory.
* ``velocity`` : heliocentric velocity of a galaxy (e.g., -300 km/s for M31).
* ``distances``: distance grid range parameters. ``[min, max, step]``, or ``[fixed number]``.
* ``distance_unit``: specify magnitude (``units.mag``) or a length unit.
* ``distance_prior_model``: specify a prior for distance parameter.

Stellar parameters
^^^^^^^^^^^^^^^^^^

There are a set of parameters parameters always required and then other parameters that 
are needed depending on if isochrones or evolutionary tracks are used.  Evolutionary tracks
have a user specified mass range with each mass having an age spacing specific that resolved
the evolution for that mass.  Isochrones have a user specified age range with each age having 
a mass spacing set by the website used.

* Always

  - ``age_prior_model``: specify a prior for age parameter.
  - ``mass_prior_model``: specify a stellar IMF.
  - ``z``: metallicity grid points.
  - ``met_prior_model``: specify a prior for metallicity parameter.
  - ``oiso``: isochrone or evolutionary model grid. See :ref:`BEAST grid inputs <beast_grid_inputs>`
  - ``osl``: stellar library definition. See :ref:`BEAST grid inputs <beast_grid_inputs>` for choices.

* Evolutionary Tracks

  - ``logmass``: log mass grid range parameters (min, max, step).
  - ``condense``: boolean if the age spacing should be condensed for each initial mass based on the following parameters.
  - ``condense_logT_delta``: requested spacing in logT along individual evolutionary track
  - ``condense_logL_delta``: requested spacing in logL along individual evolutionary track

Setting ``condense == True`` means the spacing for each mass is checked and points removed that have a delta
in logT and logL that are less than the specified values.

* Isochrones

  - ``logt``: log age grid range parameters (min, max, step).

The mass spacing at each age is not user controllable.  It is set by the website from which the 
isochrone is downloaded.

Dust parameters
^^^^^^^^^^^^^^^

The dust extinction model can be set to a single model or a combination of two models.
[TBA], picking between a linear and log A(V) spacing.

Example of a single model: 

.. code-block:: python
  
     extLaw = extinction.Generalized_DustExt(curve='G23')

Example of a mixture model: 

.. code-block:: python

     extA = extinction.Generalized_DustExt(curve='G23')
     extB = extinction.Generalized_DustExt(curve='G03_SMCBar')
     extLaw = extinction.Generalized_RvFALaw(ALaw=extA, BLaw=extB)

* Always

  - ``extLaw``: extinction law definition.
  - ``avs``: dust column in magnitudes (A_V) grid range parameters (min, max, step).
  - ``av_prior_model``: prior for A_V parameter.
  - ``rvs``: average dust grain size grid (R_V) range parameters (min, max, step).
  - ``rv_prior_model``: prior for R_V parameter.

* Mixture

  - ``fAs``: mixture factor between "MW" and "SMCBar" extinction curves (f_A) grid range parameters (min, max, step).
  - ``fA_prior_model``: prior for f_A parameter.

Artificial Star Test (AST) File Parameters
------------------------------------------

The BEAST generates artificial star test (AST) input files based on additional
input parameters from beast_settings.txt.  The ASTs are generated from the physics model
grid

* ``ast_models_selected_per_age``: number of models to pick per age (default = 70).
* ``ast_bands_above_maglimit``: number of filters that must be above the magnitude limit for an AST to be included in the list (default = 3).
* ``ast_realization_per_model``: number of realizations of each included AST model to be put into the list (default = 20).
* ``ast_maglimit``: two options: (1) number of magnitudes fainter than the 90th percentile faintest star in the photometry catalog to be used for the mag cut (default = 1); (2) custom faint end limits (space-separated list of numbers, one for each band).
* ``ast_with_positions``:  (optional; bool) if True, the AST list is produced with X,Y positions. If False, the AST list is produced with only magnitudes.
* ``ast_density_table``: (optional; string) name of density table, containing either the source density map or the background density map. If supplied, the ASTs will be repeated for each density bin in the table (default = None).
* ``ast_N_bins``: (optional; int) number of source or background bins that you want ASTs repeated over.
* ``ast_pixel_distribution``: (optional; float) minimum pixel separation between AST position and catalog star used to determine the AST spatial distribution. Used if ast_with_positions is True.
* ``ast_reference_image``: (optional; string) name of the reference image used by DOLPHOT when running the measured photometry. Required if ast_with_positions is True and no X,Y information is present in the photometry catalog.
* ``ast_reference_image_hdu_extension``: (optional; int) extension number of the reference image file where the WCS information is stored. Required if ast_with_positions is True and no X,Y information is present in the photometry catalog.
* ``ast_coord_boundary``: (optional; list of two arrays) if supplied, these RA/Dec coordinates will be used to limit the region over which ASTs are generated (default = None).
* ``ast_erode_selection_region``: (optional; float) To avoid placing ASTs near the edge of the image, set this to the number of arcseconds (default=0.5, which is ~10 pixels for WFC3/UVIS) to shrink the allowed AST placement region.  This is applied by doing an erosion to both ast_coord_boundary (if set) and a convex hull around the photometry catalog.
* ``astfile``:  pathname to the AST files (single camera ASTs).
* ``ast_colnames``:  names of columns for filters in the AST catalog (default is the basefilter list).
* ``noisefile`` : pathname to the output noise model file.
* ``absflux_a_matrix`` : absolute flux calibration covariance matrix for HST specfic filters.

Optional Features
-----------------

Add additional filters to grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Define list of filternames as ``additional_filters`` and alter ``add_spectral_properties`` call:

``add_spectral_properties_kwargs = dict(filternames=filters + additional_filters)``

Allow non-interrupting warnings in verify_params
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Set ``allow_verify_warnings`` boolean variable in beast_settings.txt to allow non-interrupting warnings. Default: raise UserWarning exception.

``allow_verify_warnings = True``

Remove constant SFH prior
^^^^^^^^^^^^^^^^^^^^^^^^^
Add ``prior_kwargs`` to beast_settings.txt:

``prior_kwargs = dict(constantSFR=False)``

Add kwargs defining code block before ``add_stellar_priors()`` call in run_beast.py:

.. code-block:: python

  if hasattr(settings, 'prior_kwargs'):
    prior_kwargs = settings.prior_kwargs
  else:
    prior_kwargs = {}

Enable Exponential Av Prior
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set ``av_prior_model`` in beast_settings.txt:

``av_prior_model = {'name': 'exponential', 'a': 2.0, 'N': 4.0}``


BEAST Filters
=============

The filters are defined in ``beast/libs/filters.hd5``. The file
has been updated in Feb 2024 using stsynphot (HST/GALEX) and 
pandeia (JWST) to have correct, 
total throughput for HST filters and to remove unused filters. 
The file contains two groups:

* ``content``: fields are ``TABLENAME`` (string), ``OBSERVATORY``
  (string), ``INSTRUMENT`` (string), ``NORM`` (float), ``CWAVE`` (float),
  ``PWAVE`` (float), ``COMMENT`` (string)

* ``filters`` has a group for each filter, with the same names as
  ``TABLENAME``.  The groups contain a dataset with the fields
  ``WAVELENGTH`` (float array, in Angstroms) and ``THROUGHPUT``
  (float array).

The filters currently included in the BEAST filter library are as follows.

Please do not forget updating ``beast/libs/vega.hd5`` as well when making 
any updates in ``beast/libs/filters.hd5``. Vega fluxes and magnitudes in 
updated filters need to be correspondingly recomputed and saved in vega.hd5.
See :doc:`beast_filters` for the full current set of included filters.
