#####
BEAST
#####

The Bayesian Extinction and Stellar Tool (BEAST) fits the ultraviolet to
near-infrared photometric SEDs of stars to extract stellar and
dust extinction parameters.
The stellar parameters are age (t), mass (M), metallicity (M), and distance (d).
The dust extinction parameters are dust column (Av), average grain size (Rv),
and mixing between type A and B extinction curves (fA).

The full details of the BEAST are provide by
Gordon et al. (2016, ApJ, 826, 104).
<http://adsabs.harvard.edu/abs/2016ApJ...826..104G>

User Documentation
==================

.. toctree::
   :maxdepth: 2

   BEAST run setup details <beast_setup.rst>
   Example production run workflow <workflow.rst>
   Running in parallel by using subgrids <subgrid_parallelism.rst>
   Generating AST inputs <generating_asts.rst>
   Format of BEAST grid files <beast_grid_format.rst>

Installation
============

.. toctree::
   :maxdepth: 2

   How to install <install.rst>

Developer Documentation
=======================

.. toctree::
   :maxdepth: 2

   How to contribute <beast_development.rst>

Reporting Issues
================

If you have found a bug in ``beast`` please report it by creating a
new issue on the ``beast`` `GitHub issue tracker
<https://github.com/BEAST-Fitting/beast/issues>`_.

Please include an example that demonstrates the issue sufficiently so that
the developers can reproduce and fix the problem. You may also be asked to
provide information about your operating system and a full Python
stack trace.  The developers will walk you through obtaining a stack
trace if it is necessary.

Contributing
============

Like the `Astropy`_ project, ``beast`` is made both by and for its
users.  We accept contributions at all levels, spanning the gamut from
fixing a typo in the documentation to developing a major new feature.
We welcome contributors who will abide by the `Python Software
Foundation Code of Conduct
<https://www.python.org/psf/codeofconduct/>`_.

``beast`` follows the same workflow and coding guidelines as
`Astropy`_.  The following pages will help you get started with
contributing fixes, code, or documentation (no git or GitHub
experience necessary):

* `How to make a code contribution <http://astropy.readthedocs.io/en/stable/development/workflow/development_workflow.html>`_

* `Coding Guidelines <http://docs.astropy.io/en/latest/development/codeguide.html>`_

* `Try the development version <http://astropy.readthedocs.io/en/stable/development/workflow/get_devel_version.html>`_

* `Developer Documentation <http://docs.astropy.org/en/latest/#developer-documentation>`_


For the complete list of contributors please see the `beast
contributors page on Github
<https://github.com/BEAST-Fitting/beast/graphs/contributors>`_.

Reference API
=============
.. toctree::
   :maxdepth: 1

   physicsmodel_api.rst
   observationmodel_api.rst
   fitting_api.rst
