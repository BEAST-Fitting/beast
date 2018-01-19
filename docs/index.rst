#####
BEAST
#####

The Bayesian Extinction and Stellar Tool (BEAST) fits the ultraviolet to
near-infrared photometric SEDs of stars to extract stellar and
dust extinction parameters.
The stellar parameters are age (t), mass (M), and metallicity (M).
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
   Generating AST inputs <generating_asts.rst>
   Format of BEAST grid files <beast_grid_format.rst>

Installation
============

.. toctree::
   :maxdepth: 2

   How to intall <install.rst>

Developer Documentation
=======================

.. toctree::
   :maxdepth: 2

   How to contribute <beast_development.rst>

Repository
==========

Github: <https://github.com/BEAST-Fitting/beast>


Reference API
=============
.. toctree::
   :maxdepth: 1

   physicsmodel_api.rst
   observationmodel_api.rst
   fitting_api.rst
