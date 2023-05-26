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

Artificial Star Test (AST) File Parameters
------------------------------------------

The BEAST generates artificial star test (AST) input files based on additional
input parameters from beast_settings.txt.

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

Grid Definition Parameters
--------------------------

The BEAST generates a grid of stellar models based on aditional input parameters
from `beast_settings.txt`. See <beast_grid_inputs.rst> for details on model libraries.
For more on setting up priors, see :ref:`BEAST priors <beast_priors>`.

* ``n_subgrid``: number of sub-grids to use (1 means no subgrids), useful for when
  the physics model grid is too large to read into memory.
* ``velocity`` : heliocentric velocity of a galaxy (e.g., -300 km/s for M31).
* ``distances``: distance grid range parameters. ``[min, max, step]``, or ``[fixed number]``.
* ``distance_unit``: specify magnitude (``units.mag``) or a length unit.
* ``distance_prior_model``: specify a prior for distance parameter.
* Stellar parameters

  - ``logt``: age grid range parameters (min, max, step).
  - ``age_prior_model``: specify a prior for age parameter.
  - ``mass_prior_model``: specify a stellar IMF.
  - ``z``: metallicity grid points. For PARSEC, 1e-4<=z<=0.02; For MIST, -4.0<=[Z/H]<=0.5
  - ``met_prior_model``: specify a prior for metallicity parameter.
  - ``oiso``: isochrone model grid. Current choices: Padova or MIST. Default: PARSEC+CALIBRI: ``oiso = isochrone.PadovaWeb()``
  - ``osl``: stellar library definition. Options include Kurucz, Tlusty, BTSettl, Munari, Elodie and BaSel. You can also generate an object from the union of multiple individual libraries: ``osl = stellib.Tlusty() + stellib.Kurucz()``

* Dust parameters

  - ``extLaw``: extinction law definition.
  - ``avs``: dust column in magnitudes (A_V) grid range parameters (min, max, step).
  - ``av_prior_model``: specify a prior for A_V parameter.
  - ``rvs``: average dust grain size grid (R_V) range parameters (min, max, step).
  - ``rv_prior_model``: specify a prior for R_V parameter.
  - ``fAs``: mixture factor between "MW" and "SMCBar" extinction curves (f_A) grid range parameters (min, max, step).
  - ``avs``: dust column in magnitudes (A_V) grid range parameters (min, max, step).
  - ``av_prior_model``: specify a prior for A_V parameter.
  - ``rvs``: average dust grain size grid (R_V) range parameters (min, max, step).
  - ``rv_prior_model``: specify a prior for R_V parameter.
  - ``fAs``: mixture factor between "MW" and "SMCBar" extinction curves (f_A) grid range parameters (min, max, step).
  - ``fA_prior_model``: specify a prior for f_A parameter.


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
has been updated in March, 2022 using stsynphot to have correct, 
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
udpated filters need to be correspondingly recomputed and saved in vega.hd5.

+--------------------------+
| HST_WFC3_F218W           |
+--------------------------+
| HST_WFC3_F225W           |
+--------------------------+
| HST_WFC3_F275W           |
+--------------------------+
| HST_WFC3_F336W           |
+--------------------------+
| HST_WFC3_F390M           |
+--------------------------+
| HST_WFC3_F390W           |
+--------------------------+
| HST_WFC3_F410M           |
+--------------------------+
| HST_WFC3_F438W           |
+--------------------------+
| HST_WFC3_F467M           |
+--------------------------+
| HST_WFC3_F475W           |
+--------------------------+
| HST_WFC3_F547M           |
+--------------------------+
| HST_WFC3_F555W           |
+--------------------------+
| HST_WFC3_F606W           |
+--------------------------+
| HST_WFC3_F621M           |
+--------------------------+
| HST_WFC3_F625W           |
+--------------------------+
| HST_WFC3_F689M           |
+--------------------------+
| HST_WFC3_F763M           |
+--------------------------+
| HST_WFC3_F775W           |
+--------------------------+
| HST_WFC3_F814W           |
+--------------------------+
| HST_WFC3_F845M           |
+--------------------------+
| HST_WFC3_F098M           |
+--------------------------+
| HST_WFC3_F105W           |
+--------------------------+
| HST_WFC3_F110W           |
+--------------------------+
| HST_WFC3_F125W           |
+--------------------------+
| HST_WFC3_F127M           |
+--------------------------+
| HST_WFC3_F139M           |
+--------------------------+
| HST_WFC3_F140W           |
+--------------------------+
| HST_WFC3_F153M           |
+--------------------------+
| HST_WFC3_F160W           |
+--------------------------+
| HST_WFPC2_F122M          |
+--------------------------+
| HST_WFPC2_F157W          |
+--------------------------+
| HST_WFPC2_F336W          |
+--------------------------+
| HST_WFPC2_F410M          |
+--------------------------+
| HST_WFPC2_F467M          |
+--------------------------+
| HST_WFPC2_F547M          |
+--------------------------+
| HST_WFPC2_F439W          |
+--------------------------+
| HST_WFPC2_F569W          |
+--------------------------+
| HST_WFPC2_F675W          |
+--------------------------+
| HST_WFPC2_F791W          |
+--------------------------+
| HST_WFPC2_F170W          |
+--------------------------+
| HST_WFPC2_F185W          |
+--------------------------+
| HST_WFPC2_F218W          |
+--------------------------+
| HST_WFPC2_F255W          |
+--------------------------+
| HST_WFPC2_F300W          |
+--------------------------+
| HST_WFPC2_F380W          |
+--------------------------+
| HST_WFPC2_F555W          |
+--------------------------+
| HST_WFPC2_F622W          |
+--------------------------+
| HST_WFPC2_F450W          |
+--------------------------+
| HST_WFPC2_F606W          |
+--------------------------+
| HST_WFPC2_F702W          |
+--------------------------+
| HST_WFPC2_F814W          |
+--------------------------+
| HST_ACS_WFC_F435W        |
+--------------------------+
| HST_ACS_WFC_F475W        |
+--------------------------+
| HST_ACS_WFC_F550M        |
+--------------------------+
| HST_ACS_WFC_F555W        |
+--------------------------+
| HST_ACS_WFC_F606W        |
+--------------------------+
| HST_ACS_WFC_F625W        |
+--------------------------+
| HST_ACS_WFC_F775W        |
+--------------------------+
| HST_ACS_WFC_F814W        |
+--------------------------+
| GALEX_FUV                |
+--------------------------+
| GALEX_NUV                |
+--------------------------+
