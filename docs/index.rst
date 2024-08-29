#####
BEAST
#####

The Bayesian Extinction and Stellar Tool (BEAST) is a Python package that fits
the ultraviolet to near-infrared photometric SEDs of stars to extract stellar
and dust extinction parameters.
The stellar parameters are age (t), mass (M), metallicity (M), and distance (d).
The dust extinction parameters are dust column (Av), average grain size (Rv),
and mixing between type A and B extinction curves (fA).

The full details of the BEAST are provided by
Gordon et al. (2016, ApJ, 826, 104): http://adsabs.harvard.edu/abs/2016ApJ...826..104G

Getting started
===============

.. toctree::
   :maxdepth: 1

   Installation instructions <install.rst>
   Running an example <example.rst>

User Documentation
==================

Basics:

.. toctree::
   :maxdepth: 1

   Run setup <beast_setup.rst>
   Filters supported <beast_filters.rst>
   Photometry files <photometry_files.rst>
   Output files <outputs.rst>

Detals:

.. toctree::
   :maxdepth: 1

   Graphical Models <beast_graphical_model.rst>
   Stellar/Extinction Priors <beast_priors.rst>
   Generating AST inputs <generating_asts.rst>
   Example production run workflow <workflow.rst>
   Running in parallel by using subgrids <subgrid_parallelism.rst>
   Plotting Tools <plotting_tools.rst>
   Analysis Tools <analysis_tools.rst>
   Other Tools <other_tools.rst>
   Format of BEAST grid files <beast_grid_format.rst>
   Details on BEAST libraries for grid <beast_grid_inputs.rst>
   Known issues <beast_issues.rst>


Publications
============
A list of publications using the BEAST, can be found :ref:`here<publications>`.


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

Like the `Astropy`_ project, the ``beast`` is made both by and for its
users.  We accept contributions at all levels, spanning the gamut from fixing a
typo in the documentation to developing a major new feature. We welcome
contributors who will abide by the `Python Software Foundation Code of Conduct
<https://www.python.org/psf/conduct/>`_.

More details on how to contribute to the ``beast`` can be found in the Developer Documentation:

.. toctree::
   :maxdepth: 1

   BEAST Contribution Workflow <beast_development.rst>
   Internal Classes <development_classes.rst>

For the complete list of contributors, please see the `beast
contributors page on Git/hub
<https://github.com/BEAST-Fitting/beast/graphs/contributors>`_.


Reference API
=============
.. toctree::
   :maxdepth: 1

   physicsmodel_api.rst
   observationmodel_api.rst
   fitting_api.rst
   plotting_api.rst
   tools_api.rst
