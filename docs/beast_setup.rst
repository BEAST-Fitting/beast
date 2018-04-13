####################
Setting Up the BEAST
####################

Basics
======

1) Define project and grid input parameters in datamodel.py

2) Execute BEAST Run using ``python run_beast.py`` with appropriate task flags

   * Default Full Stack Run: ``python run_beast.py -p -o -t -f``

.. _beast_setup_datamodel:

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

Grid Definition Parameters
--------------------------

The BEAST generates a grid of stellar models based on aditional input parameters
from datamodel.py.

* ``distances``: distance grid range parameters. ``[min, max, step]``, or ``[fixed number]``.
* ``distance_unit``: specify magnitude (``units.mag``) or a length unit
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


BEAST Filters
=============

The filters are defined in ``beast/libs/filters.hd5``.  The file
contains two groups:

* ``content``: fields are ``TABLENAME`` (string), ``OBSERVATORY``
  (string), ``INSTRUMENT`` (string), ``NORM`` (float), ``CWAVE`` (float),
  ``PWAVE`` (float), ``COMMENT`` (string)

* ``filters`` has a group for each filter, with the same names as
  ``TABLENAME``.  The groups contain a dataset with the fields
  ``WAVELENGTH`` (float array, in Angstroms) and ``THROUGHPUT``
  (float array).

The filters currently included in the BEAST filter library are as follows.

+--------------------------+
| GROUND_JOHNSON_U         |
+--------------------------+
| GROUND_JOHNSON_B         |
+--------------------------+
| GROUND_JOHNSON_V         |
+--------------------------+
| GROUND_COUSINS_R         |
+--------------------------+
| GROUND_COUSINS_I         |
+--------------------------+
| GROUND_BESSELL_J         |
+--------------------------+
| GROUND_BESSELL_H         |
+--------------------------+
| GROUND_BESSELL_K         |
+--------------------------+
| HST_NIC2_F110W           |
+--------------------------+
| HST_NIC2_F160W           |
+--------------------------+
| HST_NIC2_F205W           |
+--------------------------+
| HST_WFPC2_F218W          |
+--------------------------+
| HST_ACS_HRC_F220W        |
+--------------------------+
| HST_ACS_HRC_F250W        |
+--------------------------+
| HST_WFPC2_F255W          |
+--------------------------+
| HST_WFPC2_F300W          |
+--------------------------+
| HST_ACS_HRC_F330W        |
+--------------------------+
| HST_WFPC2_F336W          |
+--------------------------+
| HST_ACS_HRC_F344N        |
+--------------------------+
| HST_ACS_HRC_F435W        |
+--------------------------+
| HST_ACS_WFC_F435W        |
+--------------------------+
| HST_WFPC2_F439W          |
+--------------------------+
| HST_WFPC2_F450W          |
+--------------------------+
| HST_ACS_HRC_F475W        |
+--------------------------+
| HST_ACS_WFC_F475W        |
+--------------------------+
| HST_ACS_HRC_F502N        |
+--------------------------+
| HST_ACS_WFC_F502N        |
+--------------------------+
| HST_ACS_HRC_F550M        |
+--------------------------+
| HST_ACS_WFC_F550M        |
+--------------------------+
| HST_ACS_HRC_F555W        |
+--------------------------+
| HST_ACS_WFC_F555W        |
+--------------------------+
| HST_WFPC2_F555W          |
+--------------------------+
| HST_ACS_HRC_F606W        |
+--------------------------+
| HST_ACS_WFC_F606W        |
+--------------------------+
| HST_WFPC2_F606W          |
+--------------------------+
| HST_WFPC2_F622W          |
+--------------------------+
| HST_ACS_HRC_F625W        |
+--------------------------+
| HST_ACS_WFC_F625W        |
+--------------------------+
| HST_ACS_HRC_F658N        |
+--------------------------+
| HST_ACS_WFC_F658N        |
+--------------------------+
| HST_ACS_HRC_F660N        |
+--------------------------+
| HST_ACS_WFC_F660N        |
+--------------------------+
| HST_WFPC2_F675W          |
+--------------------------+
| HST_ACS_HRC_F775W        |
+--------------------------+
| HST_ACS_WFC_F775W        |
+--------------------------+
| HST_WFPC2_F791W          |
+--------------------------+
| HST_ACS_HRC_F814W        |
+--------------------------+
| HST_ACS_WFC_F814W        |
+--------------------------+
| HST_WFPC2_F814W          |
+--------------------------+
| HST_ACS_HRC_F850LP       |
+--------------------------+
| HST_ACS_WFC_F850LP       |
+--------------------------+
| HST_WFPC2_F850LP         |
+--------------------------+
| HST_ACS_HRC_F892N        |
+--------------------------+
| HST_ACS_WFC_F892N        |
+--------------------------+
| CFHT_CFH12K_CFH7406      |
+--------------------------+
| CFHT_CFH12K_CFH7504      |
+--------------------------+
| CFHT_MEGAPRIME_CFH7605   |
+--------------------------+
| CFHT_MEGAPRIME_CFH7701   |
+--------------------------+
| CFHT_MEGAPRIME_CFH7803   |
+--------------------------+
| CFHT_WIRCAM_CFH8002      |
+--------------------------+
| CFHT_WIRCAM_CFH8101      |
+--------------------------+
| CFHT_WIRCAM_CFH8102      |
+--------------------------+
| CFHT_WIRCAM_CFH8103      |
+--------------------------+
| CFHT_WIRCAM_CFH8104      |
+--------------------------+
| CFHT_WIRCAM_CFH8201      |
+--------------------------+
| CFHT_WIRCAM_CFH8202      |
+--------------------------+
| CFHT_WIRCAM_CFH8203      |
+--------------------------+
| CFHT_WIRCAM_CFH8204      |
+--------------------------+
| CFHT_WIRCAM_CFH8301      |
+--------------------------+
| CFHT_WIRCAM_CFH8302      |
+--------------------------+
| CFHT_WIRCAM_CFH8303      |
+--------------------------+
| CFHT_WIRCAM_CFH8304      |
+--------------------------+
| CFHT_WIRCAM_CFH8305      |
+--------------------------+
| CFHT_MEGAPRIME_CFH9301   |
+--------------------------+
| CFHT_MEGAPRIME_CFH9401   |
+--------------------------+
| CFHT_MEGAPRIME_CFH9601   |
+--------------------------+
| CFHT_MEGAPRIME_CFH9701   |
+--------------------------+
| CFHT_MEGAPRIME_CFH9801   |
+--------------------------+
| HST_WFC3_F098M           |
+--------------------------+
| HST_WFC3_F105W           |
+--------------------------+
| HST_WFC3_F110W           |
+--------------------------+
| HST_WFC3_F125W           |
+--------------------------+
| HST_WFC3_F126N           |
+--------------------------+
| HST_WFC3_F127M           |
+--------------------------+
| HST_WFC3_F128N           |
+--------------------------+
| HST_WFC3_F130N           |
+--------------------------+
| HST_WFC3_F132N           |
+--------------------------+
| HST_WFC3_F139M           |
+--------------------------+
| HST_WFC3_F140W           |
+--------------------------+
| HST_WFC3_F153M           |
+--------------------------+
| HST_WFC3_F160W           |
+--------------------------+
| HST_WFC3_F164N           |
+--------------------------+
| HST_WFC3_F167N           |
+--------------------------+
| HST_WFC3_F200LP          |
+--------------------------+
| HST_WFC3_F218W           |
+--------------------------+
| HST_WFC3_F225W           |
+--------------------------+
| HST_WFC3_F275W           |
+--------------------------+
| HST_WFC3_F280N           |
+--------------------------+
| HST_WFC3_F300X           |
+--------------------------+
| HST_WFC3_F336W           |
+--------------------------+
| HST_WFC3_F343N           |
+--------------------------+
| HST_WFC3_F350LP          |
+--------------------------+
| HST_WFC3_F373N           |
+--------------------------+
| HST_WFC3_F390M           |
+--------------------------+
| HST_WFC3_F390W           |
+--------------------------+
| HST_WFC3_F395N           |
+--------------------------+
| HST_WFC3_F410M           |
+--------------------------+
| HST_WFC3_F438W           |
+--------------------------+
| HST_WFC3_F467M           |
+--------------------------+
| HST_WFC3_F469N           |
+--------------------------+
| HST_WFC3_F475W           |
+--------------------------+
| HST_WFC3_F475X           |
+--------------------------+
| HST_WFC3_F487N           |
+--------------------------+
| HST_WFC3_F502N           |
+--------------------------+
| HST_WFC3_F547M           |
+--------------------------+
| HST_WFC3_F555W           |
+--------------------------+
| HST_WFC3_F600LP          |
+--------------------------+
| HST_WFC3_F606W           |
+--------------------------+
| HST_WFC3_F621M           |
+--------------------------+
| HST_WFC3_F625W           |
+--------------------------+
| HST_WFC3_F631N           |
+--------------------------+
| HST_WFC3_F645N           |
+--------------------------+
| HST_WFC3_F656N           |
+--------------------------+
| HST_WFC3_F657N           |
+--------------------------+
| HST_WFC3_F658N           |
+--------------------------+
| HST_WFC3_F665N           |
+--------------------------+
| HST_WFC3_F673N           |
+--------------------------+
| HST_WFC3_F680N           |
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
| HST_WFC3_F850LP          |
+--------------------------+
| HST_WFC3_F953N           |
+--------------------------+
| HST_WFC3_FQ232N          |
+--------------------------+
| HST_WFC3_FQ243N          |
+--------------------------+
| HST_WFC3_FQ378N          |
+--------------------------+
| HST_WFC3_FQ387N          |
+--------------------------+
| HST_WFC3_FQ422M          |
+--------------------------+
| HST_WFC3_FQ436N          |
+--------------------------+
| HST_WFC3_FQ437N          |
+--------------------------+
| HST_WFC3_FQ492N          |
+--------------------------+
| HST_WFC3_FQ508N          |
+--------------------------+
| HST_WFC3_FQ575N          |
+--------------------------+
| HST_WFC3_FQ619N          |
+--------------------------+
| HST_WFC3_FQ634N          |
+--------------------------+
| HST_WFC3_FQ672N          |
+--------------------------+
| HST_WFC3_FQ674N          |
+--------------------------+
| HST_WFC3_FQ727N          |
+--------------------------+
| HST_WFC3_FQ750N          |
+--------------------------+
| HST_WFC3_FQ889N          |
+--------------------------+
| HST_WFC3_FQ906N          |
+--------------------------+
| HST_WFC3_FQ924N          |
+--------------------------+
| HST_WFC3_FQ937N          |
+--------------------------+
| HST_NIC3_F108N           |
+--------------------------+
| HST_NIC3_F110W           |
+--------------------------+
| HST_NIC3_F113N           |
+--------------------------+
| HST_NIC3_F150W           |
+--------------------------+
| HST_NIC3_F160W           |
+--------------------------+
| HST_NIC3_F164N           |
+--------------------------+
| HST_NIC3_F166N           |
+--------------------------+
| HST_NIC3_F175W           |
+--------------------------+
| HST_NIC3_F187N           |
+--------------------------+
| HST_NIC3_F190N           |
+--------------------------+
| HST_NIC3_F196N           |
+--------------------------+
| HST_NIC3_F200N           |
+--------------------------+
| HST_NIC3_F205M           |
+--------------------------+
| HST_NIC3_F212N           |
+--------------------------+
| HST_NIC3_F215N           |
+--------------------------+
| HST_NIC3_F222M           |
+--------------------------+
| HST_NIC3_F240M           |
+--------------------------+
| CFHT_MEGAPRIME_CFH9702   |
+--------------------------+
| HST_WFPC2_F170W          |
+--------------------------+
| GALEX_FUV                |
+--------------------------+
| GALEX_NUV                |
+--------------------------+
| GROUND_2MASS_J           |
+--------------------------+
| GROUND_2MASS_H           |
+--------------------------+
| GROUND_2MASS_Ks          |
+--------------------------+
| SPITZER_IRAC_36          |
+--------------------------+
| SPITZER_IRAC_45          |
+--------------------------+
| SPITZER_IRAC_58          |
+--------------------------+
| SPITZER_IRAC_80          |
+--------------------------+
| WISE_RSR_W1              |
+--------------------------+
| WISE_RSR_W2              |
+--------------------------+
| WISE_RSR_W3              |
+--------------------------+
| WISE_RSR_W4              |
+--------------------------+
| GROUND_SDSS_U            |
+--------------------------+
| GROUND_SDSS_G            |
+--------------------------+
| GROUND_SDSS_R            |
+--------------------------+
| GROUND_SDSS_I            |
+--------------------------+
| GROUND_SDSS_Z            |
+--------------------------+


.. _BTSettl:  https://phoenix.ens-lyon.fr/Grids/BT-Settl/
.. _TLusty:  http://nova.astro.umd.edu/Tlusty2002/database/
.. _Munari:  http://archives.pd.astro.it/2500-10500/
.. _BaSel:  http://www.astro.unibas.ch/BaSeL_files/BaSeL2_2.tar.gz
