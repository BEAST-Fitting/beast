Running the BEAST
=================

1) Define project and grid input parameters in datamodel.py

2) Execute BEAST Run using ``python run_beast.py`` with appropriate task flags

   * Default Full Stack Run: ``python run_beast.py -p -o -t -f``

BEAST Data Model
================

Primary Parameters
------------------

Project Details

* ``project``: sets name of working subdirectory
* ``filters``: names of photometric filter passbands (matching library names)
* ``basefilters``: short versions of passband names
* ``obsfile``: filename for input flux data
* ``obs_colnames``: column names in ``obsfile`` for observed fluxes
* ``astfile``: filename for AST results file
* ``distanceModulus``: distance modulus to target in mags

Grid Definition Parameters

* ``logt``: age grid range parameters (min, max, step)
* ``z``: metallicity grid points
* ``oiso``: stellar model definition
* ``osl``: stellar library definition
* ``extLaw``: extinction law definition
* ``avs``: A_V grid range parameters (min, max, step)
* ``rvs``: R_V grid range parameters (min, max, step)
* ``fAs``: f_A grid range parameters (min, max, step)
* ``*_prior_model``: prior model definitions for dust parameters (default: flat prior)

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

.. codeblock:: python

  if hasattr(datamodel, 'prior_kwargs'):
    prior_kwargs = datamodel.prior_kwargs
  else:
    prior_kwargs = {}

Enable Exponential Av Prior
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set ``av_prior_model`` in datamodel.py:

``av_prior_model = {'name': 'exponential', 'a': 2.0, 'N': 4.0}``
