BEAST
=====

The Bayesian Extinction and Stellar Tool (BEAST) fits the ultraviolet to
near-infrared photometric SEDs of stars to extract stellar and
dust extinction parameters.
The stellar parameters are age (t), mass (M), metallicity (M), and distance (d).
The dust extinction parameters are dust column (Av), average grain size (Rv),
and mixing between type A and B extinction curves (fA).

The full details of the BEAST are provided by
`Gordon et al. (2016, ApJ, 826, 104)
<https://ui.adsabs.harvard.edu/abs/2016ApJ...826..104G/abstract>`_.

Build Status/Checks
-------------------

.. image:: http://readthedocs.org/projects/beast/badge/?version=latest
   :target: http://beast.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://github.com/BEAST-Fitting/beast/workflows/Python%20Tests/badge.svg
   :target: https://github.com/BEAST-Fitting/beast/actions/
   :alt: Test Status

.. image:: https://codecov.io/gh/BEAST-Fitting/beast/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/BEAST-Fitting/beast
   :alt: Test Coverage Status

.. image:: https://img.shields.io/lgtm/grade/python/g/BEAST-Fitting/beast.svg?logo=lgtm&logoWidth=18
    :target: https://lgtm.com/projects/g/BEAST-Fitting/beast/context:python
    :alt: Language grade: Python

.. image:: https://app.codacy.com/project/badge/Grade/3ebbffbfb2634e6fad620b4931be3cbc
    :target: https://www.codacy.com/gh/BEAST-Fitting/beast/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=BEAST-Fitting/beast&amp;utm_campaign=Badge_Grade
    :alt: Codacy grade

Packaging
---------

.. image:: https://badge.fury.io/py/beast.svg
   :target: https://badge.fury.io/py/beast

.. image:: https://img.shields.io/badge/ascl-1908.013-blue.svg?colorB=262255
   :target: http://ascl.net/1908.013
   :alt: ascl:1908.013

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy

Documentation
-------------

Details of installing, running, and contributing to the BEAST are at
<http://beast.readthedocs.io>.

.. image:: https://img.shields.io/badge/ApJ-Gordon%20et%20al.%202016,%20ApJ,%20826,%20104-brightgreen
    :target: https://ui.adsabs.harvard.edu/abs/2016ApJ...826..104G/abstract
    :alt: ApJ paper

.. image:: http://img.shields.io/badge/arXiv-1606.06182-orange.svg?style=flat
    :target: https://arxiv.org/abs/1606.06182
    :alt: arXiv paper

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
<https://ui.adsabs.harvard.edu/abs/2016ApJ...826..104G/abstract>`_
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

We love contributions! beast is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

*This disclaimer was originally written by*
`Adrienne Lowe <https://github.com/adriennefriend>`_ *for a*
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, *and was adapted by
the BEAST based on its use in the README file for the*
`MetPy project <https://github.com/Unidata/MetPy>`_.

License
-------

This project is Copyright (c) Karl Gordon and BEAST Team and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.

.. _AstroPy: http://www.astropy.org/
.. _contributing: http://docs.astropy.org/en/stable/index.html#contributing
.. _developer: http://docs.astropy.org/en/stable/index.html#developer-documentation
.. _Astropy Code of Conduct:  http://www.astropy.org/about.html#codeofconduct
