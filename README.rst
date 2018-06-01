BEAST
=====

.. image:: http://readthedocs.org/projects/beast/badge/?version=latest
   :target: http://beast.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://travis-ci.org/BEAST-Fitting/beast.svg?branch=master
    :target: https://travis-ci.org/BEAST-Fitting/beast
    :alt: Build Status

.. image:: https://coveralls.io/repos/github/BEAST-Fitting/beast/badge.svg?branch=master
    :target: https://coveralls.io/github/BEAST-Fitting/beast?branch=master
    :alt: Coverage Status

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy

.. image:: http://img.shields.io/badge/arXiv-1606.06182-orange.svg?style=flat
    :target: https://arxiv.org/abs/1606.06182
    :alt: arXiv paper

The Bayesian Extinction and Stellar Tool (BEAST) fits the ultraviolet to
near-infrared photometric SEDs of stars to extract stellar and
dust extinction parameters.
The stellar parameters are age (t), mass (M), metallicity (M), and distance (d).
The dust extinction parameters are dust column (Av), average grain size (Rv),
and mixing between type A and B extinction curves (fA).

The full details of the BEAST are provided by
`Gordon et al. (2016, ApJ, 826, 104)
<http://adsabs.harvard.edu/abs/2016ApJ...826..104G>`_.

Documentation
-------------

Details of installing, running, and contributing to the BEAST are at
<http://beast.readthedocs.io>.

Contributors
------------

BEAST contributors 2016 and before (BEAST paper authorship):
Karl D. Gordon, Morgan Fouesneau, Heddy Arab, Kirill Tchernyshyov,
Daniel R. Weisz, Julianne J. Dalcanton, Benjamin F. Williams,
Eric F. Bell, Lucianna Bianchi, Martha Boyer, Yumi Choi, Andrew Dolphin,
Leo Girardi, David W. Hogg, Jason S. Kalirai, Maria Kapala,
Alexia R. Lewis, Hans-Walter Rix, Karin Sandstrom, and Evan D. Skillman

Direct code contributors (including new contributors since 2016):
<https://github.com/BEAST-fitting/beast/graphs/contributors>

Attribution
-----------

Please cite `Gordon et al. (2016, ApJ, 826, 104)
<http://adsabs.harvard.edu/abs/2016ApJ...826..104G>`_
if you find this code useful in your research.
The BibTeX entry for the paper is::

    @ARTICLE{2016ApJ...826..104G,
      author = {{Gordon}, K.~D. and {Fouesneau}, M. and {Arab}, H. and {Tchernyshyov}, K. and
          {Weisz}, D.~R. and {Dalcanton}, J.~J. and {Williams}, B.~F. and
          {Bell}, E.~F. and {Bianchi}, L. and {Boyer}, M. and {Choi}, Y. and
          {Dolphin}, A. and {Girardi}, L. and {Hogg}, D.~W. and {Kalirai}, J.~S. and
          {Kapala}, M. and {Lewis}, A.~R. and {Rix}, H.-W. and {Sandstrom}, K. and
          {Skillman}, E.~D.},
      title = "{The Panchromatic Hubble Andromeda Treasury. XV. The BEAST: Bayesian Extinction and Stellar Tool}",
      journal = {\apj},
      archivePrefix = "arXiv",
      eprint = {1606.06182},
      keywords = {dust, extinction, galaxies: individual: M31, methods: data analysis, methods: statistical, stars: fundamental parameters},
      year = 2016,
      month = aug,
      volume = 826,
      eid = {104},
      pages = {104},
      doi = {10.3847/0004-637X/826/2/104},
      adsurl = {http://adsabs.harvard.edu/abs/2016ApJ...826..104G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

In Development!
---------------

This code is currently in active development.

Contributing
------------

Please open a new issue or new pull request for bugs, feedback, or new features
you would like to see.   If there is an issue you would like to work on, please
leave a comment and we will be happy to assist.   New contributions and
contributors are very welcome!

New to github or open source projects?  If you are unsure about where to start
or haven't used github before, please feel free to contact `@karllark`.
Want more information about how to make a contribution?  Take a look at
the astropy `contributing`_ and `developer`_ documentation.

Feedback and feature requests?   Is there something missing you would like
to see?  Please open an issue or send an email to  `@karllark`.
BEAST follows the `Astropy Code of Conduct`_ and strives to provide a
welcoming community to all of our users and contributors.


License
-------

The BEAST is licensed under a 3-clause BSD style license (see the
``licenses/LICENSE.rst`` file).


.. _AstroPy: http://www.astropy.org/
.. _contributing: http://docs.astropy.org/en/stable/index.html#contributing
.. _developer: http://docs.astropy.org/en/stable/index.html#developer-documentation
.. _Astropy Code of Conduct:  http://www.astropy.org/about.html#codeofconduct
