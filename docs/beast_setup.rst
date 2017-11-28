
..  _beast_setup:

Setting Up the BEAST
====================

1) Define project and grid input parameters in datamodel.py

2) Execute BEAST Run using ``python run_beast.py`` with appropriate task flags

   * Default Full Stack Run: ``python run_beast.py -p -o -t -f``

BEAST Data Model
================

Before running the BEAST, you will need to modify datamodel.py to specify the required parameters for generating models and fitting data. These parameters (and example values) are described below.

Project Details
---------------

* ``project``: pathname of of working subdirectory.
* ``filters``: names of photometric filter passbands (matching library names).
* ``basefilters``: short versions of passband names.
* ``obsfile``: filename for input flux data.
* ``obs_colnames``: column names in ``obsfile`` for observed fluxes. The input data MUST be in fluxes, NOT in magnitudes and the fluxes MUST be in normalized Vega units.
* ``distanceModulus``: distance modulus to target in magnitudes.


Artificial Star Test (AST) File Parameters
------------------------------------------
          
The BEAST generates artificial star test (AST) input files based on additional
input parameters from datamodel.py. 

* ``ast_models_selected_per_age``: number of models to pick per age (default = 70).
* ``ast_bands_above_maglimit``: number of filters that must be above the magnitude limit for an AST to be included in the list (default = 3).
* ``ast_realization_per_model``: number of realizations of each included AST model to be put into the list (default = 20).
* ``ast_maglimit``: two options: (1) number of magnitudes fainter than the 90th percentile faintest star in the photometry catalog to be used for the mag cut (default = 1); (2) custom faint end limits (space-separated list of numbers, one for each band). 
* ``ast_with_positions``:  (optional; bool) if True, the AST list is produced with X,Y positions. If False, the AST list is produced with only magnitudes.
* ``ast_pixel_distribution``: (optional; float) minimum pixel separation between AST position and catalog star used to determine the AST spatial distribution. Used if ast_with_positions is True.
* ``ast_reference_image``: (optional; string)	name of the reference image used by DOLPHOT when running the measured photometry. Required if ast_with_positions is True and no X,Y information is present in the photometry catalog.        
* ``astfile``:  pathname to the AST files (single camera ASTs).
* ``noisefile`` : pathname to the output noise model file.
* ``distanceModulus``: distance modulus to the galaxy, in magnitudes.

Grid Definition Parameters
--------------------------

The BEAST generates a grid of stellar models based on aditional input parameters
from datamodel.py.

* ``logt``: age grid range parameters (min, max, step).
* ``z``: metallicity grid points.
* ``oiso``: isochrone model grid. Current choices: Padova or MIST. Default: PARSEC+CALIBRI: ``oiso = isochrone.PadovaWeb(modeltype='parsec12s', filterPMS=True)``
* ``osl``: stellar library definition. Options include Kurucz, `Tlusty`_, `BTSettl`_, Munari, Elodie and BaSel. You can also generate an object from the union of multiple individual libraries: ``osl = stellib.Tlusty() + stellib.Kurucz()``

* ``extLaw``: extinction law definition.

* ``avs``: dust column in magnitudes (A_V) grid range parameters (min, max, step).
* ``rvs``: average dust grain size grid (R_V) range parameters (min, max, step).
* ``fAs``: mixture factor between "MW" and "SMCBar" extinction curves (f_A) grid range parameters (min, max, step).
* ``*_prior_model``: prior model definitions for dust parameters (A_V, R_V, f_A). Default: flat prior.

Optional Features
-----------------

Add additional filters to grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Define list of filternames as ``additional_filters`` and alter ``add_spectral_properties`` call:

``add_spectral_properties_kwargs = dict(filternames=filters + additional_filters)``

Skip verify_params exit
^^^^^^^^^^^^^^^^^^^^^^^
Add ``noexit=True`` keyword to ``verify_input_format()`` call in run_beast.py:

``verify_params.verify_input_format(datamodel, noexit=True)``

Remove constant SFH prior
^^^^^^^^^^^^^^^^^^^^^^^^^
Add ``prior_kwargs`` to datamodel.py:

``prior_kwargs = dict(constantSFR=False)``

Add kwargs defining code block before ``add_stellar_priors()`` call in run_beast.py:

.. code-block:: python

  if hasattr(datamodel, 'prior_kwargs'):
    prior_kwargs = datamodel.prior_kwargs
  else:
    prior_kwargs = {}

Enable Exponential Av Prior
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set ``av_prior_model`` in datamodel.py:

``av_prior_model = {'name': 'exponential', 'a': 2.0, 'N': 4.0}``

          
.. _BTSettl:  https://phoenix.ens-lyon.fr/Grids/BT-Settl/
.. _TLusty:  http://nova.astro.umd.edu/Tlusty2002/database/
.. _Munari:  http://archives.pd.astro.it/2500-10500/
.. _BaSel:  http://www.astro.unibas.ch/BaSeL_files/BaSeL2_2.tar.gz
